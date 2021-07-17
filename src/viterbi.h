#ifndef VITERBI_H
#define VITERBI_H

/* Viterbi algorithm. */

/* The LENGTH macro ultimately returns an int. */
void viterbi(int* ret, int* obs, int len);

typedef struct trellis_st trellis_t;
void trellis_init(trellis_t** trellis, int len);
void trellis_destroy(trellis_t** trellis, int len);
void viterbi_fill_trellis(trellis_t* trellis, int* obs, int len);
void viterbi_choose_path(int* ret, trellis_t* trellis, int len);

#endif  /* VITERBI_H */
