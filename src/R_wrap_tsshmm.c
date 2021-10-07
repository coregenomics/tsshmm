/** @file

    @brief Wrap C functions to accept R S-expression object arguments.
 */

#include "R_wrap_tsshmm.h"
#include "models.h"
#include "parameters.h"
#include "simulate.h"
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
    @param proseq Whether to use PRO-seq as background instead of GRO-cap.
    @return The nil object
*/
SEXP
C_model_tsshmm(SEXP external_pointer, SEXP proseq)
{
  ghmm_dmodel* model = NULL;
  model_tsshmm(&model, *INTEGER(proseq));
  R_SetExternalPtrAddr(external_pointer, model);
  /* Destroy the model when the R object is garbage collected, including on
     exit of R. */
  R_RegisterCFinalizerEx(external_pointer, C_model_destroy, TRUE);
  return R_NilValue;
}

/** Read model sizes to allocate storage for R matrices.

    @param n_states Output number of model states.
    @param n_emis Output number of model emissions.
    @param model Pointer to initialized HMM.
    @return The nil object
 */
SEXP
C_model_sizes(SEXP n_states, SEXP n_emis, SEXP model)
{
  model_sizes(INTEGER(n_states), INTEGER(n_emis), R_ExternalPtrAddr(model));
  return R_NilValue;
}

/** Read model transition and emission matrices.

    Not all transitions are valid; invalid transitions are filled with NA
    values.

    @param trans Output transition flat matrix.
    @param emis Output emission flat matrix.
    @param model Pointer to initialized HMM.
    @return The nil object
 */
SEXP
C_model_matrices(SEXP trans, SEXP emis, SEXP model)
{
  model_matrices(REAL(trans), REAL(emis), R_ExternalPtrAddr(model));
  return R_NilValue;
}

/** Write model transition and emission matrices.

    @param trans Input transition flat matrix.
    @param emis Input emission flat matrix.
    @param model Pointer to initialized HMM to update.
    @return The nil object
 */
SEXP
C_model_set_matrices(SEXP trans, SEXP emis, SEXP model)
{
  model_set_matrices(REAL(trans), REAL(emis), R_ExternalPtrAddr(model));
  return R_NilValue;
}


/** Read tied emissions from the model.

    @param tied_emis Output 1-based indices of tied emissions.
    @param model Pointer to initialized HMM.
    @return The nil object
 */
SEXP
C_model_tied_emis(SEXP tied_emis, SEXP model)
{
  model_tied_emis(INTEGER(tied_emis), R_ExternalPtrAddr(model));
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

/** Simulate data from a hidden Markov model.

    @param obs Output encoded integer observations.
    @param model HMM to simulate from.
    @return The nil object
*/
SEXP
C_simulate(SEXP obs, SEXP model)
{
  simulate(R_ExternalPtrAddr(model), INTEGER(obs), nrows(obs), ncols(obs));
  return R_NilValue;
}

/** Return the most probably Viterbi path.

    @param hidden_states Output of most probably hidden state path.
    @param model Pointer to the trained HMM.
    @param observations Encoded integer observations.
    @param lengths Segmentation of observations to allow parallel calculation.
    @return The nil object
 */
SEXP
C_viterbi(SEXP hidden_states, SEXP model, SEXP observations, SEXP lengths)
{
  viterbi(INTEGER(hidden_states),
	  R_ExternalPtrAddr(model),
	  INTEGER(observations),
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
