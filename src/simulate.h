/** @file

    @brief Simulate observations from a hidden Markov model.
 */

#ifndef SIMULATE_H
#define SIMULATE_H

#include <ghmm/model.h>

void simulate(ghmm_dmodel* model, int* obs, int ncol, int nrow);

#endif	/* SIMULATE_H */
