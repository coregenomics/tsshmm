/** @file

    @brief Decode hidden states from observations using trained hidden
    Markov model.
 */

#ifndef VITERBI_H
#define VITERBI_H

#include <ghmm/model.h>

void viterbi(int* hidden_states, ghmm_dmodel* model, int* obs, int* lengths,
	     int n);

#endif  /* VITERBI_H */
