/** @file

    @brief Generate pre-designed hidden Markov models.
 */

#ifndef MODELS_H
#define MODELS_H

#include <ghmm/model.h>

void model_init(ghmm_dmodel** model, int* is_valid, int* dim, double* trans,
		double* emis, int* emis_tied, double* start);

#endif	/* MODELS_H */
