/** @file

    @brief Read current parameters of hidden Markov models.
 */

#ifndef MODEL_READ_H
#define MODEL_READ_H

#include <ghmm/model.h>

void model_sizes(int* n_states, int* n_emis, ghmm_dmodel* model);
void model_matrices(double* trans, double* emis, ghmm_dmodel* model);
void model_set_matrices(double* trans, double* emis, ghmm_dmodel* model);

#endif	/* MODEL_READ_H */
