/** @file

    @brief Decode hidden states from observations using trained hidden
    Markov model.
 */

#include <stdlib.h>		/* free */

#include <R_ext/RS.h>
#include <ghmm/viterbi.h>

#include "viterbi.h"

/** Return the most probable Viterbi path.

    @param hidden_states Optimal Viterbi path of hidden states for each
    observation.
    @param model Pointer to the trained HMM.
    @param obs Encoded integer observations.
    @param lengths Segmentation of observations to allow parallel calculation.
    @param n Number of observations.
 */
void
viterbi(int* hidden_states, ghmm_dmodel* model, int* obs, int* lengths, int n)
{
  int* cumsum = Calloc(n, int);
  cumsum[0] = 0;
  for (int i = 1; i < n; ++i) {
    cumsum[i] = cumsum[i-1] + lengths[i-1];
  }

#pragma omp parallel
  {
#pragma omp for nowait
    for (int i = 0; i < n; ++i) {
      int* ret = hidden_states + cumsum[i];
      int path_len = -1;
      double log_p = 0;

      /* Pre-allocating the temporary memory used to calculate the Viterbi path
	 would be greatly reduce the number of allocations to the number of
	 threads instead of allocating n times.  But this would require
	 rewriting ghmm_dmodel_viterbi() */
      int* path =
	ghmm_dmodel_viterbi(model, obs + cumsum[i], lengths[i], &path_len,
			    &log_p);

      /* Copy results and cleanup. */
      for (int j = 0; j < lengths[i]; ++j) {
	ret[j] = path[j];
      }
      /* It's not clear whether we can use R's Free() on memory not allocated
	 by Calloc(). */
      free(path);
      path = NULL;
    }
  }

  Free(cumsum);
}
