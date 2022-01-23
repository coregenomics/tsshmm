/** @file

    @brief Patch GHMM reestimate.c to renormalize after each Baum-Welch loop.
 */

#ifndef TRAIN_GHMM_NORMALIZE_PATCH_H
#define TRAIN_GHMM_NORMALIZE_PATCH_H

#include <ghmm/model.h>

/* The "real" one-line patch happens in the reestimate_one_step() static
   function which is private to GHMM's reestimate.c, but these 2 functions
   below are for consistency keeping the GHMM's higher level API. */

int ghmm_dmodel_baum_welch_normalize(ghmm_dmodel* model, ghmm_dseq* seq);
int ghmm_dmodel_baum_welch_nstep_normalize(ghmm_dmodel* model, ghmm_dseq* seq,
					   int max_iter,
					   double likelihood_delta);

#endif /* TRAIN_GHMM_NORMALIZE_PATCH_H */
