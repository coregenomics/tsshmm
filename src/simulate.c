/** @file

    @brief Simulate observations from a hidden Markov model.
 */

#include <R_ext/Random.h>

#include "simulate.h"

/** Generate a random deviate and return index of probability that exceeds it.

    @param model HMM to lookup probabilities.
    @param n Number of indices to loop and check.
    @param from Previous value required to lookup model probability.
    @param f Lookup function for model emission, transition, or start state.
    @return Index of probable emission, sparse transition, or starting state.
 */
int
sample(ghmm_dmodel* model, int n, int from, double (*f)(ghmm_dmodel*, int, int))
{
  GetRNGstate();
  double deviate = unif_rand();
  PutRNGstate();

  double sum = 0.0;
  int to = 0;
  for (to = 0; to < n; ++to) {
    sum += (*f)(model, from, to);
    if (sum >= deviate)
      break;
  }
  return to;
}

/** Lookup function for start state.

    @param model HMM to lookup probabilities.
    @param unused Integer value for consistent lookup function signature.
    @param state Initial state index to lookup.
    @return Probability of starting state.
 */
double
prob_start(ghmm_dmodel* model, int unused, int state)
{
  (void)unused; /* Canonical suppression of unused variable compiler warning. */
  return model->s[state].pi;
}

/** Lookup function for model emission.

    @param model HMM to lookup probabilities.
    @param state Transition state index to lookup.
    @param ob Emission state index to lookup.
    @return Probability of emission.
 */
double
prob_emis(ghmm_dmodel* model, int state, int ob)
{
  return model->s[state].b[ob];
}

/** Lookup function for model transition.

    @param model HMM to lookup probabilities.
    @param from Transition state index to lookup.
    @param to Transition state index to lookup.
    @return Probability of transition.
 */
double
prob_trans(ghmm_dmodel* model, int from, int to)
{
  return model->s[from].out_a[to];
}

/** Generate and point C GHMM generated sequences to R integer vectors.

    Our version of ghmm_dmodel_generate_sequences() that reuses the R vector.

    @param model HMM to simulate from.
    @param obs Output encoded integer observations.
    @param nrow Number of obs matrix rows.
    @param ncol Number of obs matrix columns.
*/
void
simulate(ghmm_dmodel* model, int* obs, int ncol, int nrow)
{
  /* Our function definition swaps the ncol and nrow values so that we can
     write normal looking row-major code for a column-major S-expression. */

  for (int i = 0; i < nrow; ++i) {
    /* Sample an initial state and observation. */
    int state = sample(model, model->N, -1, prob_start);
    obs[ncol * i] = sample(model, model->M, state, prob_emis);

    for (int j = 1; j < ncol; ++j) {
      /* Get next state and observation. */
      int idx = sample(model, model->s[state].out_states, state, prob_trans);
      /* Note that the state out_* arrays are sparse, which is why we have to
	 dereference the output array "idx" to get the actual state. */
      state = model->s[state].out_id[idx];
      obs[ncol * i + j] = sample(model, model->M, state, prob_emis);
    }
  }
}
