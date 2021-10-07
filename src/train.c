/** @file

    @brief Train pre-defined hidden Markov models with large genomic datasets.
 */

#include <R_ext/RS.h>

#include <ghmm/reestimate.h>

/** Allocate and point C GHMM discrete sequence type to R integer vectors.

    Our version of ghmm_dseq_calloc() that reuses the R vector.

    @param seq Output GHMM discrete sequence object of observations.
    @param obs Encoded integer observations.
    @param lengths Segmentation of observations to allow discontiguous training.
    @param n Number of observation sequences.
*/
void
sequence_alloc(ghmm_dseq **seq, int* obs, int* lengths, int n)
{
    *seq = Calloc(1, ghmm_dseq);
    (*seq)->seq_number = n;
    (*seq)->seq_len = lengths;
    (*seq)->seq = Calloc(n, int*);
    (*seq)->seq_w = Calloc(n, double);
    int cumsum = 0;
    for (int i = 0; i < n; ++i) {
      (*seq)->seq[i] = obs + cumsum;
      cumsum += lengths[i];
      (*seq)->seq_w[i] = 1.0;
    }
}

/** Free C GHMM discrete sequence type from pointing to R vectors.

    Our version of ghmm_dseq_free() that preserves the R vector.

    @param seq Output GHMM discrete sequence object of observations.
*/
void
sequence_free(ghmm_dseq **seq)
{
  Free((*seq)->seq_w);
  Free((*seq)->seq);
  Free(*seq);
}

/** Run Baum-Welch training using streaming logic for large input data.

    Run a single step of Buam-Welch training, because we cannot fit the entire
    training sequence from realistic data in RAM.

    @param converged Output of 0 if converged and -1 otherwise.
    @param model HMM to train.
    @param obs Encoded integer observations.
    @param lengths Segmentation of observations to allow discontiguous training.
    @param n Number of observation sequences.
*/
void
train(int* converged, ghmm_dmodel* model, int* obs, int* lengths, int n)
{
  /* Pad our */
  ghmm_dseq *seq = NULL;
  sequence_alloc(&seq, obs, lengths, n);
  /* Set the likelihood-delta to zero to ignore the likelihood-delta condition
     to complete a single pass of all the data instead of terminating
     prematurely. */
  *converged = ghmm_dmodel_baum_welch_nstep(model, seq, 1, 0);
  /* Update initial states to normalized B->N1 and B-P1 transitions. */
  enum { B, N1, P1 = 4 };
  double sum = model->s[B].out_a[1] + model->s[B].out_a[2];
  model->s[N1].pi = model->s[B].out_a[1] / sum;
  model->s[P1].pi = model->s[B].out_a[2] / sum;
  sequence_free(&seq);
}
