#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>   /* attribute_visible */

#include "R_wrap_tsshmm.h"

/* Register C routines for R. */

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef call_methods[] =
  {
   CALLDEF(C_model_tsshmm, 1),
   CALLDEF(C_tss, 6),
   CALLDEF(C_viterbi, 3),
   { NULL, NULL, 0 }            /* Terminating NULL entry. */
  };

void
attribute_visible
R_init_tsshmm(DllInfo *info)
{
  R_registerRoutines(info, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
