#include <Rdefines.h>

/* Wrap C functions to accept R S-expression object arguments. */

SEXP C_model_tsshmm(SEXP);
SEXP C_tss(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_viterbi(SEXP, SEXP, SEXP);
