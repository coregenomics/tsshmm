#include <R_ext/Arith.h>

#include "model_read.h"

/** Read model sizes to allocate storage for R matrices.

    @param n_states Output number of model states.
    @param n_emis Output number of model emissions.
    @param model Pointer to initialized HMM.
 */
void
model_sizes(int* n_states, int* n_emis, ghmm_dmodel* model)
{
  *n_states = model->N;
  *n_emis = model->M;
}

/** Read model transition and emission matrices.

    Not all transitions are valid; invalid transitions are filled with NA
    values.

    @param trans Output transition flat matrix.
    @param emis Output emission flat matrix.
    @param model Pointer to initialized HMM.
 */
void
model_matrices(double* trans, double* emis, ghmm_dmodel* model)
{
  for (int i = 0; i < model->N; ++i) {
    for (int j = 0; j < model->N; ++j)
      trans[i * model->N + j] =
	ghmm_dmodel_check_transition(model, i, j) ?
	ghmm_dmodel_get_transition(model, i, j) : NA_REAL;

    for (int j = 0; j < model->M; ++j)
      emis[i * model->M + j] = model->s[i].b[j];
  }
}
