/** @file

    @brief Wrap C functions to accept R S-expression object arguments.
 */

#ifndef R_WRAP_TSSHMM
#define R_WRAP_TSSHMM

#include <Rdefines.h>

SEXP C_model_tsshmm(SEXP);
SEXP C_model_sizes(SEXP, SEXP, SEXP);
SEXP C_model_matrices(SEXP, SEXP, SEXP);
SEXP C_model_set_matrices(SEXP, SEXP, SEXP);
SEXP C_train(SEXP, SEXP, SEXP, SEXP);
SEXP C_viterbi(SEXP, SEXP, SEXP, SEXP);
SEXP C_tss(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif  /* R_WRAP_TSSHMM */
