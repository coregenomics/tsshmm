/** @file

    @brief Wrap C functions to accept R S-expression object arguments.
 */

#include "R_wrap_tsshmm.h"
#include "models.h"
#include "train.h"
#include "tss.h"
#include "viterbi.h"

/** Frees the memory of a TSS hidden Markov model.

    @param external_pointer Pointer to the TSS hidden Markov model C object.
 */
void
C_model_destroy(SEXP external_pointer)
{
  ghmm_dmodel* model = R_ExternalPtrAddr(external_pointer);
  R_ClearExternalPtr(external_pointer);
  ghmm_dmodel_free(&model);
}

/** Allocates a TSS hidden Markov model for training and Viterbi decoding.

    @param external_pointer Pointer to the TSS hidden Markov model C object.
    @return The nil object
*/
SEXP
C_model_tsshmm(SEXP external_pointer)
{
  ghmm_dmodel* model = NULL;
  model_tsshmm(&model);
  R_SetExternalPtrAddr(external_pointer, model);
  /* Destroy the model when the R object is garbage collected, including on
     exit of R. */
  R_RegisterCFinalizerEx(external_pointer, C_model_destroy, TRUE);
  return R_NilValue;
}

/** Run Baum-Welch training using streaming logic for large input data.

    Run a single step of Buam-Welch training, because we cannot fit the entire
    training sequence from realistic data in RAM.

    @param converged Output of 0 if converged and -1 otherwise.
    @param model HMM to train.
    @param obs Encoded integer observations.
    @param lengths Segmentation of observations to allow discontiguous training.
    @return The nil object
*/
SEXP
C_train(SEXP converged, SEXP model, SEXP obs, SEXP lengths)
{
  train(INTEGER(converged),
	R_ExternalPtrAddr(model),
	INTEGER(obs),
	INTEGER(lengths),
	LENGTH(lengths));
  return R_NilValue;
}

/** Find peaks using 3 basepair tie-breaking.

    @param indices_peak Output integer vector to store peak locations.
    @param groups Sequential number indicating signal count group membership.
    @param indices_signal Absolute locations of signal counts.
    @param starts_signal Relative locations of signal counts relative to region.
    @param scores_signal Values of signal counts.
    @param prefer_last Boolean whether to effectively scan from the last count.
    @return The nil object
 */
SEXP
C_tss(SEXP indices_peak, SEXP groups, SEXP indices_signal, SEXP starts_signal,
      SEXP scores_signal, SEXP prefer_last)
{
  tss(INTEGER(indices_peak),
      INTEGER(groups),
      INTEGER(indices_signal),
      INTEGER(starts_signal),
      INTEGER(scores_signal),
      LENGTH(scores_signal),
      INTEGER(prefer_last));
  return R_NilValue;
}

/** Return the most probably Viterbi path.

    @param hidden_states Output of most probably hidden state path.
    @param observations Encoded integer observations.
    @param lengths Segmentation of observations to allow parallel calculation.
    @return The nil object
 */
SEXP
C_viterbi(SEXP hidden_states, SEXP observations, SEXP lengths)
{
  /* Reuse the trellis by allocating it once using the longest observation. */
  int len_max = INTEGER(lengths)[0];
  int* cumsum = Calloc(LENGTH(lengths), int);
  cumsum[0] = 0;
  for (int i = 1; i < LENGTH(lengths); ++i) {
    cumsum[i] = cumsum[i-1] + INTEGER(lengths)[i-1];
    if (INTEGER(lengths)[i] > len_max) {
      len_max = INTEGER(lengths)[i];
    }
  }
#pragma omp parallel
  {
    trellis_t* trellis = NULL;
    trellis_init(&trellis, len_max);

#pragma omp for nowait
    for (int i = 0; i < LENGTH(lengths); ++i) {
      int* obs = INTEGER(observations) + cumsum[i];
      int len = INTEGER(lengths)[i];
      int* ret = INTEGER(hidden_states) + cumsum[i];

      /* Run Viterbi. */
      viterbi_fill_trellis(trellis, obs, len);
      viterbi_choose_path(ret, trellis, len);
    }

    /* Cleanup. */
    trellis_destroy(&trellis, len_max);
  }
  Free(cumsum);

  return R_NilValue;
}
