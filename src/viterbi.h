#ifndef VITERBI_H
#define VITERBI_H

#include <stddef.h>		/* size_t */

/* Viterbi algorithm. */

/* The LENGTH macro ultimately returns an int. */
void viterbi(int* ret, int* obs, int len);

#endif	/* VITERBI_H */
