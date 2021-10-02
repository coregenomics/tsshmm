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
   CALLDEF(C_model_tsshmm, 2),
   CALLDEF(C_model_sizes, 3),
   CALLDEF(C_model_matrices, 3),
   CALLDEF(C_model_set_matrices, 3),
   CALLDEF(C_model_tied_emis, 2),
   CALLDEF(C_train, 4),
   CALLDEF(C_viterbi, 4),
   CALLDEF(C_tss, 6),
   { NULL, NULL, 0 }            /* Terminating NULL entry. */
  };

/** Register C routines for R and disable slow character vector .C calls.

    @param info Handle to this dynamic library so that this library's
    entry points can be registered.
*/
void
attribute_visible
R_init_tsshmm(DllInfo *info)
{
  R_registerRoutines(info, NULL, call_methods, NULL, NULL);
  /* These next two settings forbid .C calls in R from using character strings,
     which speeds up searches of routines.  The symbols additionally require
     explicit namespace access from other R packages. */
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
