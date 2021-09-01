/** @file

    @brief Decode hidden states using TSS HMM model.
 */

#ifndef VITERBI_H
#define VITERBI_H

/** Viterbi trellis type. */
typedef struct trellis_st trellis_t;
void trellis_init(trellis_t** trellis, int len);
void trellis_destroy(trellis_t** trellis, int len);
void viterbi_fill_trellis(trellis_t* trellis, int* obs, int len);
void viterbi_choose_path(int* ret, trellis_t* trellis, int len);

#endif  /* VITERBI_H */
