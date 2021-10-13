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

    @param is_valid Output whether the HMM is valid.
    @param dim Size-2 array of the number of states and emissions.
    @param trans Transitions flat matrix.
    @param emis Emissions flat matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
    @return The nil object
*/
SEXP
C_is_model_valid(SEXP is_valid, SEXP dim, SEXP trans, SEXP emis,
		 SEXP emis_tied, SEXP start)
{
  ghmm_dmodel* model = NULL;
  model_init(&model, INTEGER(is_valid), INTEGER(dim), REAL(trans), REAL(emis),
	     INTEGER(emis_tied), REAL(start));
  ghmm_dmodel_free(&model);
  return R_NilValue;
}

/** Run Baum-Welch training for input data.

    Run Buam-Welch training until convergence, for a training sequence that can
    be entirely held in RAM.

    @param converged Output of 0 if converged and -1 otherwise.
    @param obs Encoded integer observations.
    @param lengths Segmentation of observations to allow discontiguous training.
    @param dim Size-2 array of the number of states and emissions.
    @param trans Transitions flat matrix.
    @param emis Emissions flat matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
    @return The nil object
*/
SEXP
C_train(SEXP converged, SEXP obs, SEXP lengths, SEXP dim, SEXP trans,
	SEXP emis, SEXP emis_tied, SEXP start)
{
  ghmm_dmodel* model = NULL;
  int is_valid = 0;
  INTEGER(converged)[0] = 0;
  model_init(&model, &is_valid, INTEGER(dim), REAL(trans), REAL(emis),
	     INTEGER(emis_tied), REAL(start));
  if (! is_valid) {
    ghmm_dmodel_free(&model);
    lerror("model provided is invalid!  See reasons above.");
  }
  train(INTEGER(converged),
	model,
	INTEGER(obs),
	INTEGER(lengths),
	LENGTH(lengths));
  ghmm_dmodel_free(&model);
  return R_NilValue;
}

/** Simulate data from a hidden Markov model.

    @param obs Output encoded integer observations.
    @param dim Size-2 array of the number of states and emissions.
    @param trans Transitions flat matrix.
    @param emis Emissions flat matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
    @return The nil object
*/
SEXP
C_simulate(SEXP obs, SEXP dim, SEXP trans, SEXP emis, SEXP emis_tied,
	   SEXP start)
{
  ghmm_dmodel* model = NULL;
  int is_valid = 0;
  model_init(&model, &is_valid, INTEGER(dim), REAL(trans), REAL(emis),
	     INTEGER(emis_tied), REAL(start));
  if (! is_valid) {
    ghmm_dmodel_free(&model);
    lerror("model provided is invalid!  See reasons above.");
  }
  simulate(model, INTEGER(obs), nrows(obs), ncols(obs));
  ghmm_dmodel_free(&model);
  return R_NilValue;
}

/** Return the most probably Viterbi path.

    @param hidden_states Output of most probably hidden state path.
    @param observations Encoded integer observations.
    @param lengths Segmentation of observations to allow parallel calculation.
    @param dim Size-2 array of the number of states and emissions.
    @param trans Transitions flat matrix.
    @param emis Emissions flat matrix.
    @param emis_tied Tied emissions indices.
    @param start Starting probability of states.
    @return The nil object
 */
SEXP
C_viterbi(SEXP hidden_states, SEXP observations, SEXP lengths, SEXP dim,
	  SEXP trans, SEXP emis, SEXP emis_tied, SEXP start)
{
  ghmm_dmodel* model = NULL;
  int is_valid = 0;
  model_init(&model, &is_valid, INTEGER(dim), REAL(trans), REAL(emis),
	     INTEGER(emis_tied), REAL(start));
  if (! is_valid) {
    ghmm_dmodel_free(&model);
    lerror("model provided is invalid!  See reasons above.");
  }
  viterbi(INTEGER(hidden_states),
	  model,
	  INTEGER(observations),
	  INTEGER(lengths),
	  LENGTH(lengths));
  ghmm_dmodel_free(&model);
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
