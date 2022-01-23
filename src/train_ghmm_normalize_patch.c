/** @file

    @brief Patch GHMM reestimate.c to renormalize after each Baum-Welch loop.
 */

#include <float.h>
#include <math.h>

#include <ghmm/reestimate.h>
#include <ghmm/foba.h>

#include "ghmm_internals.h"
#include "train_ghmm_normalize_patch.h"

/* ----------------------------------------------------------------------------
   Pariksheet: The patch we needed to apply is a one-line change to function
   reestimate_one_step(), namely, the addition of ghmm_model_normalize().
   However, a lot of code in GHMM's reestimate.c is declared as static and
   therefore private to that file.  This required recopying static code from
   reestimate.c here as-is.  There has been some light editing to the
   reestimate.c code:

   1. Added doxygen comments for modified function reestimate_one_step(), and
      new functions ghmm_dmodel_baum_welch_normalize() and
      ghmm_dmodel_baum_welch_nstep_normalize().  Other static functions that
      are required to be included have been doxygen commented as "unchanged".

   2. Rewrote unclear comments with consistent tense and better grammar.

   3. Removed dead code that had been commented out, and from conditionals such
      as if (1) {...} else {...}.

   5. Wrapped long lines to be within 79 columns

   6. Fixed compiler warnings.

   Additionally the private header, ghmm_internals.h, has been symlinked to the
   bundled headers, and the dependencies of ghmm_internals.h, model.h and
   mes.h, have been included because of their relative include paths.
   ------------------------------------------------------------------------- */

/** Intermediate model update calculated values unchanged from reestimate.c. */
typedef struct local_store_t {
  /** Numerators of initial probabilities. */
  double* pi_num;
  /** Denominator of initial probabilities. */
  double pi_denom;
  /** Row numerators of transition probabilities. */
  double** a_num;
  /** Row denomenators of transition probabilities. */
  double* a_denom;
  /** Row numerators of emission probabilities. */
  double** b_num;
  /** Row denomenators of emission probabilities. */
  double** b_denom;
} local_store_t;


/** @def CUR_PROC 

    @brief Function name to be logged by GHMM's LOC macro.
*/


/** Unchanged copy of static function from restimate.c */
static int reestimate_free (local_store_t ** r, int N)
{
# define CUR_PROC "reestimate_free"
  int i;
  mes_check_ptr (r, return (-1));
  if (!*r)
    return (0);
  m_free ((*r)->pi_num);

  if ((*r)->a_num){
    for (i = 0; i < N; i++){
      m_free ((*r)->a_num[i]);
    }  
  }  
  m_free ((*r)->a_num);
  m_free ((*r)->a_denom);

  if ((*r)->b_num){
    for (i = 0; i < N; i++){
      m_free ((*r)->b_num[i]);
    }
  }
  m_free ((*r)->b_num);
  if ((*r)->b_denom){
    for (i = 0; i < N; i++){
      m_free ((*r)->b_denom[i]);
    }  
  }
  m_free ((*r)->b_denom);

  m_free (*r);
  return (0);
# undef CUR_PROC
} /* reestimate_free */

/** Unchanged copy of static function from restimate.c */
static local_store_t *reestimate_alloc (const ghmm_dmodel * mo)
{
# define CUR_PROC "reestimate_alloc"
  int i;
  local_store_t *r = NULL;
  ARRAY_CALLOC(r, 1);

  ARRAY_CALLOC(r->pi_num, mo->N);
  ARRAY_CALLOC(r->a_num, mo->N);
  for (i = 0; i < mo->N; i++)
    ARRAY_CALLOC(r->a_num[i], mo->s[i].out_states);
  ARRAY_CALLOC(r->a_denom, mo->N);

  ARRAY_MALLOC(r->b_num, mo->N);
  ARRAY_MALLOC(r->b_denom, mo->N);

  if (mo->model_type & GHMM_kHigherOrderEmissions)
    for (i = 0; i < mo->N; i++) {
      ARRAY_CALLOC(r->b_num[i], ghmm_ipow(mo, mo->M, mo->order[i] + 1));
      ARRAY_CALLOC(r->b_denom[i], ghmm_ipow (mo, mo->M, mo->order[i]));
    }
  else
    for (i = 0; i < mo->N; i++) {
      ARRAY_CALLOC(r->b_num[i], mo->M);
      ARRAY_CALLOC(r->b_denom[i], 1);
    }

  return (r);
STOP:				/* Label STOP from ARRAY_[CM]ALLOC */
  reestimate_free(&r, mo->N);
  return (NULL);
# undef CUR_PROC
} /* reestimate_alloc */

/** Unchanged copy of static function from restimate.c */
static int reestimate_init(local_store_t * r, const ghmm_dmodel * mo) {
# define CUR_PROC "reestimate_init"
  int i, j, m, size;
  r->pi_denom = 0.0;

  for (i=0; i<mo->N; i++) {
    r->pi_num[i] = 0.0;
    r->a_denom[i] = 0.0;
    for (j=0; j<mo->s[i].out_states; j++)
      r->a_num[i][j] = 0.0;

    if (mo->model_type & GHMM_kHigherOrderEmissions) {
      size = ghmm_ipow(mo, mo->M, mo->order[i]);
      for (m=0; m<size; m++)
	r->b_denom[i][m] = 0.0;
      size *= mo->M;
      for (m=0; m<size; m++)
	r->b_num[i][m] = 0.0;
    } else {
      r->b_denom[i][0] = 0.0;
      for (m=0; m<mo->M; m++)
	r->b_num[i][m] = 0.0;
    }
  }
  return (0);
# undef CUR_PROC
} /* reestimate_init */

/** Unchanged copy of static function from restimate.c */
static int reestimate_setlambda (local_store_t * r, ghmm_dmodel * mo)
{
# define CUR_PROC "reestimate_setlambda"
  int res = -1;
  int h, i, j, m, l, j_id, positive;
  double factor, p_i;
  int hist, col, size, reachable;

  mes_check_0 (r->pi_denom, goto STOP) {};

  for (i = 0; i < mo->N; i++) {
    reachable = 1;
    positive = 0;

    /* Pi */
    mo->s[i].pi = r->pi_num[i] / r->pi_denom;

    /* A */
    /* Note: The denominator might be 0; never reached that state? */
    p_i = 0.0;
    if (r->a_denom[i] < GHMM_EPS_PREC) {
      for (h = 0; h < mo->s[i].in_states; h++) {
        p_i += mo->s[i].in_a[h];
      }
      if (p_i == 0.0) {
        if (h == 0)
          GHMM_LOG_PRINTF(LINFO, LOC,
			  "State %d can't be reached (no in_states)", i);
        else
          GHMM_LOG_PRINTF(LINFO, LOC,
			  "State %d can't be reached (prob = 0.0)", i);
        reachable = 0;
      }
      factor = 0.0;
    }
    else
      factor = (1 / r->a_denom[i]);


    for (j = 0; j < mo->s[i].out_states; j++) {
      /* Test that the denominator < numerator. */
      if ((r->a_denom[i] - r->a_num[i][j]) <= -GHMM_EPS_PREC) {
        GHMM_LOG(LCONVERTED, " numerator > denom!\n");
      }
      mo->s[i].out_a[j] = r->a_num[i][j] * factor;
      if (r->a_num[i][j] >= GHMM_EPS_PREC)
        positive = 1;
      /* Update in_a. */
      l = 0;
      j_id = mo->s[i].out_id[j];
      while (l < mo->s[j_id].in_states)
        if (mo->s[j_id].in_id[l] == i)
          break;
        else
          l++;
      if (l == mo->s[j_id].in_states) {
        GHMM_LOG(LCONVERTED, "no corresponding in_a to out_a found!\n");
        goto STOP;
      }
      mo->s[j_id].in_a[l] = mo->s[i].out_a[j];
    }

    /* If the parameter is fixed, continue to next state. */
    if (mo->s[i].fix)
      continue;

    /* B */
    if (mo->model_type & GHMM_kHigherOrderEmissions)
      size = ghmm_ipow (mo, mo->M, mo->order[i]);
    else
      size = 1;
    /* If all in_a's are zero, the state can't be reached.
       Set all b's to -1.0 */
    if (!reachable) {
      for (hist = 0; hist < size; hist++) {
        col = hist * mo->M;
        for (m = col; m < col + mo->M; m++) {
          mo->s[i].b[m] = -1.0;
        }
      }
    }
    else {
      for (hist = 0; hist < size; hist++) {
        /* If the denominator is very small, we have not seen many emissions
           in this state with this history.
           We are conservative and just skip them. */
        if (r->b_denom[i][hist] < GHMM_EPS_PREC)
          continue;
        else
          factor = (1.0 / r->b_denom[i][hist]);

        positive = 0;
	/* Test that the denominator < numerator. */
        col = hist * mo->M;
        for (m = col; m < col + mo->M; m++) {
          if ((r->b_denom[i][hist] - r->b_num[i][m]) <= -GHMM_EPS_PREC) {
            GHMM_LOG_PRINTF(LCONVERTED, LOC,
			    "numerator b (%.4f) > denom (%.4f)!\n",
                            r->b_num[i][m], r->b_denom[i][hist]);
          }
          
          mo->s[i].b[m] = r->b_num[i][m] * factor;
          if (mo->s[i].b[m] >= GHMM_EPS_PREC)
            positive = 1;
        }

        if (!positive) {
          GHMM_LOG_PRINTF(LCONVERTED, LOC,
			  "All numerator b[%d][%d-%d] == 0 (denom = %g)!\n",
                          i, col, col + mo->M, r->b_denom[i][hist]);
        }
      }	/* for each history */
    }
  } /* for (i = 0 .. < mo->N)  */

  res = 0;
  if (mo->model_type & GHMM_kTiedEmissions)
    ghmm_dmodel_update_tie_groups (mo);
STOP:				/* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
# undef CUR_PROC
} /* reestimate_setlambda */

/** Single iteration of Baum-Welch with model renormalization.

    Modified reestimate_one_step() that normalizes after each loop of parameter
    updates.  Training over hundreds of millions of observations eventually
    causes the row sums of transition and emission probability parameters to
    barely fail the GHMM_EPS_PREC row sums tests of ghmm_dmodel_check(), and
    therefore those probability row sums need to be renormalied periodically.
    Renormalizing is quick, and therefore applied at each iteration.

    @param mo HMM to train.
    @param r Intermediate model update calculated values.
    @param seq_number Number of encoded integer observations.
    @param seq_length Lengths of encoded integer observations.
    @param O Encoded integer observations.
    @param log_p Log-likelihood of the data given the paramters.
    @param seq_w Weights of encoded integer observations.
    @return 0 on success, or -1 on error.
 */
static int reestimate_one_step (ghmm_dmodel * mo, local_store_t * r,
				int seq_number, int *seq_length, int **O,
				double *log_p, double *seq_w)
{
# define CUR_PROC "reestimate_one_step"
  int res = -1;
  int k, i, j, t, j_id, valid=0;
  int e_index;
  int errors;
  double **alpha = NULL;
  double **beta = NULL;
  double *scale = NULL;
  int T_k=0;
  double gamma;
  double log_p_k;

  /* first set maxorder to zero if model_type & kHigherOrderEmissions is FALSE 

     TODO XXX use model->maxorder only
     if model_type & kHigherOrderEmissions is TRUE */
  if (!(mo->model_type & GHMM_kHigherOrderEmissions))
    mo->maxorder = 0;

  *log_p = 0.0;
  /* Loop over all sequences. */
  for (k = 0; k < seq_number; k++) {
    mo->emission_history = 0;
    T_k = seq_length[k];        /* Current seq length. */

    /* Initialize matrices and vector dependent on T_k. */
    if (ighmm_reestimate_alloc_matvek (&alpha, &beta, &scale, T_k, mo->N) == -1) {
      GHMM_LOG_QUEUED(LCONVERTED);
      goto FREE;
    }
    if (ghmm_dmodel_forward (mo, O[k], T_k, alpha, scale, &log_p_k) == -1) {
      GHMM_LOG_QUEUED(LCONVERTED);
      goto FREE;
    }

    if (log_p_k != +1) {        /* O[k] can be generated. */
      *log_p += log_p_k;
      valid = 1;
      
      if (ghmm_dmodel_backward (mo, O[k], T_k, beta, scale) == -1) {
        GHMM_LOG_QUEUED(LCONVERTED);
        goto FREE;
      }

      /* Loop over all states. */
      for (i = 0; i < mo->N; i++) {
        /* Pi */
        r->pi_num[i] += seq_w[k] * alpha[0][i] * beta[0][i];
        r->pi_denom  += seq_w[k] * alpha[0][i] * beta[0][i];

        for (t=0; t < T_k-1; t++) {
          /* B */
          if (!mo->s[i].fix) {
            e_index = get_emission_index(mo, i, O[k][t], t);
            if (e_index != -1) {
              gamma = seq_w[k] * alpha[t][i] * beta[t][i];
              r->b_num[i][e_index] += gamma;
              r->b_denom[i][e_index / (mo->M)] += gamma;
            }
          }
          update_emission_history(mo, O[k][t]);

          /* A */
          r->a_denom[i] += (seq_w[k] * alpha[t][i] * beta[t][i]);
          for (j=0; j < mo->s[i].out_states; j++) {
            j_id = mo->s[i].out_id[j];
            e_index = get_emission_index(mo, j_id, O[k][t+1], t+1);
            if (e_index != -1)
              r->a_num[i][j] += (seq_w[k] * alpha[t][i] * mo->s[i].out_a[j]
                                 * mo->s[j_id].b[e_index] * beta[t+1][j_id]
                                 * (1.0 / scale[t + 1])); /* c[t]=1/scale[t] */
          }
        }
        /* B: last iteration for t == T_k-1 */
        t = T_k - 1;
        if (!mo->s[i].fix) {
          e_index = get_emission_index (mo, i, O[k][t], t);
          if (e_index != -1) {
            gamma = seq_w[k] * alpha[t][i] * beta[t][i];
            r->b_num[i][e_index] += gamma;
            r->b_denom[i][e_index / (mo->M)] += gamma;
          }
        }
      }
    } /* if (log_p_k != +1) */
    else {
      GHMM_LOG_PRINTF(LWARN, LOC, "O(%d) can't be built from model mo!\n", k);
    }

    ighmm_reestimate_free_matvek(alpha, beta, scale, T_k);
  } /* for (k = 0; k < seq_number; k++) */

  if (valid) {
    /* New parameter lambda: set directly in the model. */
    if (reestimate_setlambda(r, mo) == -1) {
      GHMM_LOG_QUEUED(LCONVERTED);
      goto STOP;
    }
    
    /* ------------------------------------------------------------------------
       Pariksheet: Adding just this one line below required adding everything
       else in this file!  But this patching approach was necessary so that any
       system installed GHMM can be used instead of compiling the bundled copy.
       --------------------------------------------------------------------- */
    ghmm_dmodel_normalize(mo);

    if ((errors = ghmm_dmodel_check(mo))) {
      GHMM_LOG_PRINTF(LERROR, LOC, "Reestimated model is invalid, "
                      "model_check found %d errors", -errors);
      goto STOP;
    }
  }
  else {                        /* NO sequence can be built from model mo! */
    *log_p = +1;
  }
  res = 0;

STOP:
   return (res);
FREE:
   ighmm_reestimate_free_matvek(alpha, beta, scale, T_k);
   return (res);
# undef CUR_PROC
} /* reestimate_one_step */

/** Run Baum-Welch with model renormalization.

    Clone of ghmm_dmodel_baum_welch() that replaces the call to
    ghmm_dmodel_baum_welch_nstep() with our equivalent that normalizes after
    each loop of parameter updates.

    @param model HMM to train.
    @param seq Encoded integer observations.
    @return 0 on success, or -1 on error.
 */
int
ghmm_dmodel_baum_welch_normalize(ghmm_dmodel* model, ghmm_dseq* seq)
{
# define CUR_PROC "ghmm_dmodel_baum_welch_normalize"
  return ghmm_dmodel_baum_welch_nstep_normalize(model,
						seq,
						GHMM_MAX_ITER_BW,
						GHMM_EPS_ITER_BW);
# undef CUR_PROC
} /* ghmm_dmodel_baum_welch_normalize */

/** Loop Baum-Welch with model renormalization.

    Clone of ghmm_dmodel_baum_welch_nstep() that replaces the call to the
    static function reestimate_one_step() with our equivalent that normalizes
    after each loop of parameter updates.

    @param mo HMM to train.
    @param sq Encoded integer observations.
    @param max_step Baum-Welch iteration limit.
    @param likelihood_delta Minimum relative likelihood improvement to continue.
    @return 0 on success, or -1 on error.
 */
int
ghmm_dmodel_baum_welch_nstep_normalize(ghmm_dmodel* mo, ghmm_dseq* sq,
				       int max_step, double likelihood_delta)
{
# define CUR_PROC "ghmm_dmodel_baum_welch_normalize"
  int n, k, valid;
  double log_p, log_p_old, log_p_k, diff;
  local_store_t *r = NULL;
  int res = -1;
  log_p_old = -DBL_MAX;
  n = 1;

  /* Allocate the local store of model parameter intermediary calculations. */
  r = reestimate_alloc (mo);
  if (!r) {
    goto STOP;
  };

  /* Main loop of the Baum-Welch algorithm. */
  while (n <= max_step) {
    
    if (reestimate_one_step(mo, r, sq->seq_number, sq->seq_len, sq->seq,
			    &log_p, sq->seq_w) == -1) {
      GHMM_LOG_PRINTF(LCONVERTED, LOC,
		      "reestimate_one_step false (%d.step)\n", n);
      goto STOP;
    }

    if (log_p == +1) {
      GHMM_LOG(LERROR,
	       "Reestimate stopped: No sequence can be built from model mo!");
      break;
    }

    diff = log_p - log_p_old;
    /* Error in convergence? */
    if (diff < -GHMM_EPS_PREC) {
      GHMM_LOG_PRINTF(LCONVERTED, LOC,
		      "No convergence: log P < log P-old! (n=%d)\n", n);
      goto STOP;
    }
    else if (log_p > GHMM_EPS_PREC) {
      GHMM_LOG_PRINTF(LCONVERTED, LOC,
		      "No convergence: log P > 0! (n=%d)\n", n);
      goto STOP;
    }

    /* Stop iterations? */
    if (diff < fabs ((double) likelihood_delta * log_p)) {
      GHMM_LOG_PRINTF(LINFO, LOC, "Convergence after %d steps", n);
      break;
    }
    else {
      /* For the next iteration. */
      log_p_old = log_p;
      reestimate_init (r, mo);  /* Sets all fields to zero. */
      n++;
    }
  } /* while (n <= MAX_ITER) */

  /* log_p of the reestimated model. */
  log_p = 0.0;
  valid = 0;
  for (k = 0; k < sq->seq_number; k++) {
    if (ghmm_dmodel_logp (mo, sq->seq[k], sq->seq_len[k], &log_p_k) == -1) {
      GHMM_LOG_QUEUED(LCONVERTED);
      goto STOP;
    }
    if (log_p_k != +1) {
      log_p += log_p_k;
      valid = 1;
    }
  }
  if (!valid)
    log_p = +1;

  res = 0;

STOP:				/* Label STOP from ARRAY_[CM]ALLOC */
  reestimate_free (&r, mo->N);
  return res;
# undef CUR_PROC  
} /* ghmm_dmodel_baum_welch_nstep_normalize */
