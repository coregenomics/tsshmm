#include <Rdefines.h>

/* Wrap C functions to accept R S-expression object arguments. */

SEXP C_tss(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_viterbi(SEXP, SEXP);
SEXP C_viterbi_vectorized(SEXP, SEXP);
