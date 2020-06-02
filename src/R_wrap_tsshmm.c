#include "R_wrap_tsshmm.h"
#include "tss.h"
#include "viterbi.h"

/* Wrap C functions to accept R S-expression object arguments. */

SEXP
C_tss(SEXP indices_peak, SEXP groups, SEXP indices_signal, SEXP starts_signal,
      SEXP scores_signal)
{
  tss(INTEGER(indices_peak),
      INTEGER(groups),
      INTEGER(indices_signal),
      INTEGER(starts_signal),
      INTEGER(scores_signal),
      LENGTH(scores_signal));
  return R_NilValue;
}


SEXP
C_viterbi(SEXP hidden_states,
	  SEXP observations)
{
  viterbi(INTEGER(hidden_states),
	  INTEGER(observations),
	  LENGTH(observations));
  return R_NilValue;
}
