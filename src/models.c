/** @file

    @brief Generate pre-designed hidden Markov models.
 */

#include <R_ext/Arith.h>

#include "models.h"


/** Initialize and validate an HMM.

    @param model Output pointer to initialized HMM.
    @param is_valid Output whether the HMM is valid.
    @param n_states Number of states.
    @param n_emis Number of emissions.
    @param trans Transitions matrix.
    @param emis Emissions matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
 */
void
model_init(ghmm_dmodel** model, int* is_valid, int n_states, int n_emis,
	   double* trans, double* emis, int* emis_tied, double* start)
{
  /* Initialize model. */
  int degree_out[n_states];
  int degree_in[n_states];
  for (int j = 0; j < n_states; ++j) {
    degree_out[j] = 0;
    degree_in[j] = 0;
    for (int i = 0; i < n_states; ++i) {
      if (! isnan(trans[i * n_states + j])) {
	++degree_out[j];
      }
      if (! isnan(trans[j * n_states + i])) {
	++degree_in[j];
      }
    }
  }
  *model = ghmm_dmodel_calloc(n_emis, n_states,
			      GHMM_kDiscreteHMM + GHMM_kTiedEmissions,
			      degree_in, degree_out);

  /* Shorthand model variables with more comprehensible symbols. */
  ghmm_dstate* states = (*model)->s;
  int* tied_to = (*model)->tied_to;

  for (int j = 0; j < n_states; ++j) {
    /* Transitions. */
    states[j].out_states = degree_out[j];
    states[j].in_states = degree_in[j];
    int out = 0;
    int in = 0;
    for (int i = 0; i < n_states; ++i) {
      if (! isnan(trans[i * n_states + j])) {
	states[j].out_id[out] = i;
	states[j].out_a[out] = trans[i * n_states + j];
	++out;
      }
      if (! isnan(trans[j * n_states + i])) {
	states[j].in_id[in] = i;
	states[j].in_a[in] = trans[j * n_states + i];
	++in;
      }
    }

    /* Emissions. */
    tied_to[j] = emis_tied[j] - 1;
    for (int i = 0; i < n_emis; ++i)
      states[j].b[i] = emis[i * n_states + j];

    /* Starting probabilities. */
    states[j].pi = start[j];

    /* No fixed parameters. */
    states[j].fix = 0;
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

/** Read model transition matrix.

    Not all transitions are valid; invalid transitions are filled with NA
    values.

    @param trans Output transition flat matrix.
    @param model Pointer to initialized HMM.
 */
void
model_trans(double* trans, ghmm_dmodel* model) {
  const int n_states = model->N;
  for (int j = 0; j < n_states; ++j)
    for (int i = 0; i < n_states; ++i)
      trans[i * n_states + j] =
	ghmm_dmodel_check_transition(model, j, i) ?
	ghmm_dmodel_get_transition(model, j, i) : NA_REAL;
}

/** Read model emission matrix.

    @param emis Output emission flat matrix.
    @param model Pointer to initialized HMM.
 */
void
model_emis(double* emis, ghmm_dmodel* model) {
  const int n_states = model->N;
  const int n_emis = model->M;
  for (int j = 0; j < n_states; ++j)
    for (int i = 0; i < n_emis; ++i)
      emis[i * n_states + j] = model->s[j].b[i];
}
