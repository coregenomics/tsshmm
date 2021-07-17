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


SEXP
C_viterbi_vectorized(SEXP hidden_states,
		     SEXP observations)
{
  /* Only allocate the trellis once using the longest observation. */
  int len_max = -1;
  for (int i = 0; i < length(observations); ++i) {
    if (LENGTH(VECTOR_ELT(observations, i)) > len_max) {
      len_max = LENGTH(VECTOR_ELT(observations, i));
    }
  }
  trellis_t* trellis = NULL;
  trellis_init(&trellis, len_max);

  /* Run Viterbi. */
  for (int i = 0; i < length(observations); ++i) {
    int* obs = INTEGER(VECTOR_ELT(observations, i));
    int len = LENGTH(VECTOR_ELT(observations, i));
    int* ret = INTEGER(VECTOR_ELT(hidden_states, i));

    viterbi_fill_trellis(trellis, obs, len);
    viterbi_choose_path(ret, trellis, len);
  }

  /* Cleanup. */
  trellis_destroy(&trellis, len_max);

  return R_NilValue;
}
