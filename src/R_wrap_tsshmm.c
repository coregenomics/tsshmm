#include "R_wrap_tsshmm.h"
#include "viterbi.h"

/* Wrap C functions to accept R S-expression object arguments. */

SEXP
C_viterbi(SEXP hidden_states,
	  SEXP observations)
{
  viterbi(INTEGER(hidden_states),
	  INTEGER(observations),
	  LENGTH(observations));
  return R_NilValue;
}
