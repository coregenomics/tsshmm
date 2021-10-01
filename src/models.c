/** @file

    @brief Generate pre-designed hidden Markov models.
 */

#include <assert.h>

#include "models.h"

/** Transcription start site discovery using 7-state HMM.

    This model was published in Core et al 2014 (doi:10.1038/ng.3142),
    providing the means to use signal and background nascent RNA
    sequencing data and use regions of enrichment, depletion and no
    signal to predict promoter regions.  The promoter regions are then
    used to find transcription start sites, although that processing
    step downstream from this HMM; the goal of the HMM is to find
    putative gene promoters and enhancers containing TSSs.

    @param model Output pointer to initialized 7-state TSS HMM.
    @param proseq Whether to use PRO-seq as background instead of GRO-cap.
 */
void
model_tsshmm(ghmm_dmodel** model, int proseq)
{
  /* TSS HMM from doi:10.1038/ng.3142 supp. fig. 2 */
  enum { NO_SIGNAL, ENRICHED, DEPLETED, N_EMIS };
  enum { B, N1, N2, N3, P1, P2, P3 };
  /* Extend the model to support a gene body (GB) state. */
  const int GB = P3 + 1;
  const int N_STATES = proseq ? P3 + 2 : P3 + 1;

  /* Initialize model. */
  if (! proseq) {
    /*                  B N1 N2 N3 P1 P2 P3 */
    int degree_out[] = {3, 1, 2, 1, 2, 3, 2};
    int degree_in[] =  {3, 1, 2, 1, 3, 2, 2};
    *model = ghmm_dmodel_calloc(N_EMIS, N_STATES,
				GHMM_kDiscreteHMM + GHMM_kTiedEmissions,
				degree_in, degree_out);
  } else {
    /*                  B N1 N2 N3 P1 P2 P3 GB */
    int degree_out[] = {3, 1, 2, 2, 2, 3, 3, 2};
    int degree_in[] =  {4, 1, 2, 1, 3, 2, 2, 3};
    *model = ghmm_dmodel_calloc(N_EMIS, N_STATES,
				GHMM_kDiscreteHMM + GHMM_kTiedEmissions,
				degree_in, degree_out);
  }

  /* Shorthand model variables with more comprehensible symbols. */
  ghmm_dstate* states = (*model)->s;
  int* tied_to = (*model)->tied_to;

  /* Emissions. */
  states[B].b[0] = 0.90;	/* NO_SIGNAL */
  states[B].b[1] = 0.05;	/* ENRICHED */
  states[B].b[2] = 0.05;	/* DEPLETED */
  tied_to[0] = GHMM_kUntied;
  for (int N = N1; N <= N3; ++N) {
    states[N].b[0] = 0.09;
    states[N].b[1] = 0.90;
    states[N].b[2] = 0.01;
    tied_to[N] = N1;
  }
  for (int P = P1; P <= P3; ++P) {
    states[P].b[0] = 0.10;
    states[P].b[1] = 0.45;
    states[P].b[2] = 0.45;
    tied_to[P] = P1;
  }
  if (proseq) {
    states[GB].b[0] = 0.05;
    states[GB].b[1] = 0.05;
    states[GB].b[2] = 0.90;
  }

  /* Transitions: B. */
  states[B].out_states = 3;
  states[B].out_id[0] = B;  states[B].out_a[0] = 0.99;
  states[B].out_id[1] = N1; states[B].out_a[1] = 0.005;
  states[B].out_id[2] = P1; states[B].out_a[2] = 0.005;
  states[B].in_states = 3;
  states[B].in_id[0]  = B;  states[B].in_a[0]  = 0.99;
  states[B].in_id[1]  = N3; states[B].in_a[1]  = 1.0;
  states[B].in_id[2]  = P3; states[B].in_a[2]  = 0.5;
  if (proseq) {
    states[B].in_states = 4;
    states[B].in_a[1]  = 0.25;
    states[B].in_a[2]  = 0.125;
    states[B].in_a[3]  = 1.0;
    states[B].in_id[3]  = GB;
  }

  /* Transitions: N1. */
  states[N1].out_states = 1;
  states[N1].out_id[0] = N2; states[N1].out_a[0] = 1.0;
  states[N1].in_states = 1;
  states[N1].in_id[0] = B;  states[N1].in_a[0] = 0.99;

  /* Transitions: N2. */
  states[N2].out_states = 2;
  states[N2].out_id[0] = N2; states[N2].out_a[0] = 0.5;
  states[N2].out_id[1] = N3; states[N2].out_a[1] = 0.5;
  states[N2].in_states = 2;
  states[N2].in_id[0] = N1; states[N2].in_a[0] = 1.0;
  states[N2].in_id[1] = N2; states[N2].in_a[1] = 0.5;

  /* Transitions: N3. */
  states[N3].out_states = 1;
  states[N3].out_id[0] = B; states[N3].out_a[0] = 1.0;
  if (proseq) {
    states[N3].out_states = 2;
    states[N3].out_a[0] = 0.25;
    states[N3].out_a[1] = 0.75;
    states[N3].out_id[1] = GB;
  }
  states[N3].in_states = 1;
  states[N3].in_id[0] = N2; states[N3].in_a[0] = 0.5;

  /* Transitions: P1. */
  states[P1].out_states = 2;
  states[P1].out_id[0] = P1; states[P1].out_a[0] = 0.5;
  states[P1].out_id[1] = P2; states[P1].out_a[1] = 0.5;
  states[P1].in_states = 3;
  states[P1].in_id[0] = B;  states[P1].in_a[0] = 0.005;
  states[P1].in_id[1] = P1; states[P1].in_a[1] = 0.5;
  states[P1].in_id[2] = P2; states[P1].in_a[2] = 0.45;

  /* Transitions: P2. */
  states[P2].out_states = 3;
  states[P2].out_id[0] = P1; states[P2].out_a[0] = 0.45;
  states[P2].out_id[1] = P2; states[P2].out_a[1] = 0.1;
  states[P2].out_id[2] = P3; states[P2].out_a[2] = 0.45;
  states[P2].in_states = 2;
  states[P2].in_id[0] = P1; states[P2].in_a[0] = 0.5;
  states[P2].in_id[1] = P2; states[P2].in_a[1] = 0.1;

  /* Transitions: P3. */
  states[P3].out_states = 2;
  states[P3].out_id[0] = B;  states[P3].out_a[0] = 0.5;
  states[P3].out_id[1] = P3; states[P3].out_a[1] = 0.5;
  if (proseq) {
    states[P3].out_states = 3;
    states[P3].out_a[0] = 0.125;
    states[P3].out_a[2] = 0.375;
    states[P3].out_id[2] = GB;
  }
  states[P3].in_states = 2;
  states[P3].in_id[0] = P2; states[P3].in_a[0] = 0.45;
  states[P3].in_id[1] = P3; states[P3].in_a[1] = 0.5;

  if (proseq) {
    /* Transitions: GB */
    states[GB].out_states = 2;
    states[GB].out_id[0] = B;  states[GB].out_a[0] = 0.01;
    states[GB].out_id[1] = GB; states[GB].out_a[1] = 0.99;
    states[GB].in_states = 3;
    states[GB].in_id[0] = N3; states[GB].in_a[0] = 0.75;
    states[GB].in_id[1] = P3; states[GB].in_a[0] = 0.375;
    states[GB].in_id[2] = GB; states[GB].in_a[0] = 0.99;
  }

  /* Starting probabilities and fixed parameters. */
  for (int i = 0; i < N_STATES; ++i) {
    states[i].pi = 0.0;
    states[i].fix = 0;
  }
  states[B].pi = 1.0;

  /* Other model settings. */
  (*model)->prior = -1;
  (*model)->name = NULL;
  (*model)->silent = NULL;
  (*model)->maxorder = 0;
  (*model)->emission_history = 0;
  (*model)->order = NULL;
  (*model)->bp = NULL;
  (*model)->background_id = NULL;
  (*model)->topo_order =NULL;
  (*model)->topo_order_length = 0;
  (*model)->pow_lookup = NULL;
  (*model)->label = NULL;
  (*model)->label_alphabet = NULL;
  (*model)->alphabet = NULL;

  /* Validate the model. */
  assert(ghmm_dmodel_check(*model) == 0);
}
