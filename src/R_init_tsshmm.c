#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>	/* attribute_visible */

#include "R_wrap_tsshmm.h"

/* Register C routines for R. */

static const R_CallMethodDef call_methods[] =
  {
   { "C_viterbi", (DL_FUNC) &C_viterbi, 2 },
   { NULL, NULL, 0 }		/* Terminating NULL entry. */
  };

void
attribute_visible
R_init_tsshmm(DllInfo *info)
{
  R_registerRoutines(info, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
