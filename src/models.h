/** @file

    @brief Generate pre-designed hidden Markov models.
 */

#ifndef MODELS_H
#define MODELS_H

#include <ghmm/model.h>

void model_init(ghmm_dmodel** model, int* is_valid, int n_states, int n_emis,
		double* trans, double* emis, int* emis_tied, double* start);
void model_trans(double* trans, ghmm_dmodel* model);
void model_emis(double* emis, ghmm_dmodel* model);

#endif	/* MODELS_H */
