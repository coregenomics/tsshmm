/** @file

    @brief Explicitly register C routines for R.
 */

#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>   /* attribute_visible */

#include "R_wrap_tsshmm.h"

/** Name and address of C routine with number of S-expression arguments. */
#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

/** Table of C routine entry points addresses in the dynamic library. */
static const R_CallMethodDef call_methods[] =
  {
   CALLDEF(C_model_tsshmm, 1),
   CALLDEF(C_train, 4),
   CALLDEF(C_tss, 6),
   CALLDEF(C_viterbi, 3),
   { NULL, NULL, 0 }            /* Terminating NULL entry. */
  };

/** Register C routines for R and disable slow character vector .C calls. */
void
attribute_visible
R_init_tsshmm(DllInfo *info)
{
  R_registerRoutines(info, NULL, call_methods, NULL, NULL);
  /* These two settings prevent .C calls from R to use character string to
     access to speed up searches of routines and require explicit namespace
     access from other R packages. */
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
