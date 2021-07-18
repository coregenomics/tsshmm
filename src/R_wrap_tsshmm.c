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
		     SEXP observations,
		     SEXP lengths)
{
  /* Only allocate the trellis once using the longest observation. */
  int len_max = -1;
  for (int i = 0; i < LENGTH(lengths); ++i) {
    if (INTEGER(lengths)[i] > len_max) {
      len_max = INTEGER(lengths)[i];
    }
  }
  trellis_t* trellis = NULL;
  trellis_init(&trellis, len_max);

  /* Run Viterbi. */
  int cumsum = 0;
  for (int i = 0; i < LENGTH(lengths); ++i) {
    int* obs = INTEGER(observations) + cumsum;
    int len = INTEGER(lengths)[i];
    int* ret = INTEGER(hidden_states) + cumsum;

#ifdef DEBUG_VITERBI_VECTORIZED
#include <R_ext/RS.h>
    REprintf("\n     i: %d\n", i);
    REprintf("cumsum: %d\n", cumsum);
    REprintf("   len: %d\n", len);
    REprintf("  *obs: [ 0x%X ] ", obs);
    for (int j = 0; j < len; ++j) {
      REprintf("%d ", obs[j]);
    }
    REprintf("\n");
#endif
    viterbi_fill_trellis(trellis, obs, len);
    viterbi_choose_path(ret, trellis, len);
#ifdef DEBUG_VITERBI_VECTORIZED
    REprintf("  *ret: [ 0x%X ] ", ret);
    for (int j = 0; j < len; ++j) {
      REprintf("%d ", ret[j]);
    }
    REprintf("\n");
#endif

    cumsum += len;
  }

  /* Cleanup. */
  trellis_destroy(&trellis, len_max);

  return R_NilValue;
}

SEXP
C_viterbi_vectorized_mt(SEXP hidden_states,
			SEXP observations,
			SEXP lengths)
{
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

      viterbi_fill_trellis(trellis, obs, len);
      viterbi_choose_path(ret, trellis, len);
    }

    trellis_destroy(&trellis, len_max);
  }
  Free(cumsum);

  return R_NilValue;
}
