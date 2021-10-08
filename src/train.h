/** @file

    @brief Train pre-defined hidden Markov models with large genomic datasets.
 */

#ifndef TRAIN_H
#define TRAIN_H

#include <ghmm/model.h>

void train(int* converged, ghmm_dmodel* model, int* obs, int* lengths, int n);
void train_loop(int* converged, ghmm_dmodel* model, int* obs, int* lengths,
		int n);

#endif  /* TRAIN_H */
