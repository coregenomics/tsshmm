/** @file

    @brief Generate pre-designed hidden Markov models.
 */

#include <R_ext/Arith.h>

#include "models.h"


/** Initialize and validate an HMM.

    @param model Output pointer to initialized HMM.
    @param is_valid Output whether the HMM is valid.
    @param dim Size-2 array of the number of states and emissions.
    @param trans Transitions flat matrix.
    @param emis Emissions flat matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
 */
void
model_init(ghmm_dmodel** model, int* is_valid, int* dim, double* trans,
	   double* emis, int* emis_tied, double* start)
{
  const int n_states = dim[0];
  const int n_emis = dim[1];

  /* Initialize model. */
  int degree_out[n_states];
  int degree_in[n_states];
  for (int i = 0; i < n_states; ++i) {
    degree_out[i] = 0;
    degree_in[i] = 0;
    for (int j = 0; j < n_states; ++j) {
      if (! isnan(trans[i * n_states + j])) {
	++degree_out[i];
      }
      if (! isnan(trans[j * n_states + i])) {
	++degree_in[i];
      }
    }
  }
  *model = ghmm_dmodel_calloc(n_emis, n_states,
			      GHMM_kDiscreteHMM + GHMM_kTiedEmissions,
			      degree_in, degree_out);

  /* Shorthand model variables with more comprehensible symbols. */
  ghmm_dstate* states = (*model)->s;
  int* tied_to = (*model)->tied_to;

  for (int i = 0; i < n_states; ++i) {
    /* Transitions. */
    states[i].out_states = degree_out[i];
    states[i].in_states = degree_in[i];
    int out = 0;
    int in = 0;
    for (int j = 0; j < n_states; ++j) {
      if (! isnan(trans[i * n_states + j])) {
	states[i].out_id[out] = j;
	states[i].out_a[out] = trans[i * n_states + j];
	++out;
      }
      if (! isnan(trans[j * n_states + i])) {
	states[i].in_id[in] = j;
	states[i].in_a[in] = trans[j * n_states + i];
	++in;
      }
    }

    /* Emissions. */
    tied_to[i] = emis_tied[i] - 1;
    for (int j = 0; j < n_emis; ++j)
      states[i].b[j] = emis[i * n_emis + j];

    /* Starting probabilities. */
    states[i].pi = start[i];

    /* No fixed parameters. */
    states[i].fix = 0;
  }

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

  *is_valid = ghmm_dmodel_check(*model) == 0;
}
