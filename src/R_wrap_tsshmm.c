/** @file

    @brief Wrap C functions to accept R S-expression object arguments.
 */

#include "R_wrap_tsshmm.h"
#include "models.h"
#include "simulate.h"
#include "train.h"
#include "tss.h"
#include "viterbi.h"


/** Decorate R's error() with file, function, and line. */
#define lerror(fmt) error("%s %s:%d: " fmt,__FILE__,__func__,__LINE__)


/** Allocates a TSS hidden Markov model for training and Viterbi decoding.

    @param trans Transitions matrix.
    @param emis Emissions matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
    @return Logical size-1 vector whether the HMM is valid.
*/
SEXP
C_is_model_valid(SEXP trans, SEXP emis, SEXP emis_tied, SEXP start)
{
  SEXP is_valid = PROTECT(allocVector(LGLSXP, 1));
  ghmm_dmodel* model = NULL;
  model_init(&model, INTEGER(is_valid), nrows(emis), ncols(emis), REAL(trans),
	     REAL(emis), INTEGER(emis_tied), REAL(start));
  ghmm_dmodel_free(&model);
  UNPROTECT(1);
  return is_valid;
}

/** Run Baum-Welch training for input data.

    Run Buam-Welch training until convergence, for a training sequence that can
    be entirely held in RAM.

    @param obs Encoded integer observations.
    @param lengths Segmentation of observations to allow discontiguous training.
    @param trans Transitions matrix.
    @param emis Emissions matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
    @return The new model and whether training converged.
*/
SEXP
C_train(SEXP obs, SEXP lengths, SEXP trans, SEXP emis, SEXP emis_tied,
	SEXP start)
{
  ghmm_dmodel* model = NULL;
  int is_valid = 0;
  model_init(&model, &is_valid, nrows(emis), ncols(emis), REAL(trans),
	     REAL(emis), INTEGER(emis_tied), REAL(start));
  if (! is_valid) {
    ghmm_dmodel_free(&model);
    lerror("model provided is invalid!  See reasons above.");
  }
  SEXP converged = PROTECT(allocVector(LGLSXP, 1));
  train(INTEGER(converged),
	model,
	INTEGER(obs),
	INTEGER(lengths),
	LENGTH(lengths));
  /* Copy the trained parameters from the model. */
  SEXP trans_trained =
    PROTECT(allocMatrix(REALSXP, nrows(trans), ncols(trans)));
  model_trans(REAL(trans_trained), model);
  SEXP emis_trained =
    PROTECT(allocMatrix(REALSXP, nrows(emis), ncols(emis)));
  model_emis(REAL(emis_trained), model);
  ghmm_dmodel_free(&model);
  /* Return convergence and trained parameters as a list object. */
  SEXP list = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(list, 0, converged);
  SET_VECTOR_ELT(list, 1, trans_trained);
  SET_VECTOR_ELT(list, 2, emis_trained);
  UNPROTECT(4);
  return list;
}

/** Simulate data from a hidden Markov model.

    @param dim Size-2 array of the desired output rows and columns.
    @param trans Transitions matrix.
    @param emis Emissions matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
    @return Encoded integer observations.
*/
SEXP
C_simulate(SEXP dim, SEXP trans, SEXP emis, SEXP emis_tied, SEXP start)
{
  ghmm_dmodel* model = NULL;
  int is_valid = 0;
  model_init(&model, &is_valid, nrows(emis), ncols(emis), REAL(trans),
	     REAL(emis), INTEGER(emis_tied), REAL(start));
  if (! is_valid) {
    ghmm_dmodel_free(&model);
    lerror("model provided is invalid!  See reasons above.");
  }
  SEXP obs =
    PROTECT(allocMatrix(INTSXP, INTEGER(dim)[0], INTEGER(dim)[1]));
  simulate(model, INTEGER(obs), nrows(obs), ncols(obs));
  ghmm_dmodel_free(&model);
  UNPROTECT(1);
  return obs;
}

/** Return the most probably Viterbi path.

    @param obs Encoded integer observations.
    @param lengths Segmentation of observations to allow parallel calculation.
    @param trans Transitions matrix.
    @param emis Emissions matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
    @return The most probable hidden state path.
 */
SEXP
C_viterbi(SEXP obs, SEXP lengths, SEXP trans, SEXP emis, SEXP emis_tied,
	  SEXP start)
{
  ghmm_dmodel* model = NULL;
  int is_valid = 0;
  model_init(&model, &is_valid, nrows(emis), ncols(emis), REAL(trans),
	     REAL(emis), INTEGER(emis_tied), REAL(start));
  if (! is_valid) {
    ghmm_dmodel_free(&model);
    lerror("model provided is invalid!  See reasons above.");
  }
  int sum = 0;
  for (int i = 0; i < length(lengths); ++i) {
    sum += INTEGER(lengths)[i];
  }
  SEXP hidden_states = PROTECT(allocVector(INTSXP, sum));
  viterbi(INTEGER(hidden_states),
	  model,
	  INTEGER(obs),
	  INTEGER(lengths),
	  LENGTH(lengths));
  ghmm_dmodel_free(&model);
  UNPROTECT(1);
  return hidden_states;
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
