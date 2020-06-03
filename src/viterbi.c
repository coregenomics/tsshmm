#include <R_ext/Arith.h>
#include <R_ext/RS.h>

#include "viterbi.h"

/* TSS hidden Markov model from Core et al 2014 doi:10.1038/ng.3142

   The model is intentionally hard coded to allow multiple calls to the Viterbi
   algorithm without a model function parameter, which would anyway be a nasty
   singleton a.k.a. global object.

   We can get away with hard coding the model because of the sparse transition
   states; it makes no sense to train the model because that would create
   non-sensical transitions between the peaked and non-peaked states, etc.

   The model was previously implemented using RcppHMM at
   https://gitlab.com/coregenomics/evolength/-/blob/master/R/models.R and is
   rewritten in C here to more tightly integrate several of the slower
   computation tasks.

   Also yes, the C algorithm below has been tested against several HMMs; you
   can see the unit tests at https://gitlab.com/omsai/viterbi/
 */

typedef enum { NO_SIGNAL, ENRICHED, DEPLETED, N_OBS } OBS;
typedef enum { B, N1, N2, N3, P1, P2, P3, N_STATES } STATES;

/* Viterbi algorithm. */

typedef struct trellis_node_st {
  int prev;
  double log_prob;
} trellis_node_t;

typedef struct trellis_st {
  trellis_node_t** nodes;
} trellis_t;

void
trellis_init(trellis_t** trellis, int len)
{
  *trellis = Calloc(1, trellis_t);
  (*trellis)->nodes = Calloc(len, trellis_node_t*);
  for (int i = 0; i < len; ++i)
    (*trellis)->nodes[i] = Calloc(N_STATES, trellis_node_t);
}

void
trellis_destroy(trellis_t** trellis, int len)
{
  for (int i = 0; i < len; ++i)
    Free((*trellis)->nodes[i]);
  Free((*trellis)->nodes);
  Free((*trellis));
}

void
viterbi_fill_trellis(trellis_t* trellis, int* obs, int len)
{
  /* Log transform probabilities to maintain numeric precision during long
     series of multiplications.  And, yes, it turns out that taking log(0.0) in
     our sparse state transition matrix does not adversely affect the
     calculations from inspecting the trellis below; running viterbi using
     log1p() instead produces invalid results.

     Taking log(0.0) and ignoring the -inf is done in several other
     implementations of Viterbi decoding, such as the Octave Forge statistics
     package.  Unlike higher languages, C allows you to do math with infinities
     because one can ignore the FE_DIVBYZERO C exception and adding to
     inifinite values is possible in C.

     These const log value aliases are to make the matrices more legible.  We
     also want to suppress the log(0) wrongmathcall warning from GCC.

     NB: the probability rows must all sum to 1.
  */
  const double L0    = 0.0 - INFINITY;
  const double L0005 = log(0.005);
  const double L001  = log(0.01);
  const double L005  = log(0.05);
  const double L009  = log(0.09);
  const double L01   = log(0.1);
  const double L045  = log(0.45);
  const double L05   = log(0.5);
  const double L09   = log(0.9);
  const double L099  = log(0.99);
  const double L1    = log(1.0);

  /* Sparse transition matrix:

           B    N1  N2  N3    P1  P2   P3
     B  0.99 0.005 0.0 0.0 0.005 0.0 0.00
     N1 0.00 0.000 1.0 0.0 0.000 0.0 0.00
     N2 0.00 0.000 0.5 0.5 0.000 0.0 0.00
     N3 1.00 0.000 0.0 0.0 0.000 0.0 0.00
     P1 0.00 0.000 0.0 0.0 0.500 0.5 0.00
     P2 0.00 0.000 0.0 0.0 0.450 0.1 0.45
     P3 0.50 0.000 0.0 0.0 0.000 0.0 0.50
  */
  const double trans[N_STATES][N_STATES] =
    {/*   B     N1   N2   N3     P1   P2    P3 */
     { L099, L0005,  L0,  L0, L0005,  L0,   L0 },
     /* Non-peaked TSS transitions. */
     {   L0,    L0,  L0,  L1,    L0,  L0,   L0 },
     {   L0,    L0, L05, L05,    L0,  L0,   L0 },
     {   L1,    L0,  L0,  L0,    L0,  L0,   L0 },
     /* Peaked TSS transitions. */
     {   L0,    L0,  L0,  L0,   L05, L05,   L0 },
     {   L0,    L0,  L0,  L0,  L045, L01, L045 },
     {  L05,    L0,  L0,  L0,    L0,  L0,  L05 },
    };
#ifdef PARANOID
  #include <assert.h>
  REprintf("Validating transition matrix:\n");
  for (int i = 0; i < N_STATES; ++i) {
    double sum = 0.0;
    for (int j = 0; j < N_STATES; ++j) {
      sum += trans[i][j];
      REprintf("%6.4f ", trans[i][j]);
    }
    REprintf("\n");
    assert(sum == 1.0);
  }
#endif

  /* Dense emission matrix.

      no signal enriched depleted
   B       0.90     0.05     0.05
   N1      0.09     0.90     0.01
   N2      0.09     0.90     0.01
   N3      0.09     0.90     0.01
   P1      0.10     0.45     0.45
   P2      0.10     0.45     0.45
   P3      0.10     0.45     0.45
  */
  const double emis[N_STATES][N_OBS] =
    {/*     NO_SIGNAL, ENRICHED, DEPLETED */
     /* B  */ {  L09,      L005,     L005 },
     /*         Non-peaked TSS */
     /* N1 */ { L009,       L09,     L001 },
     /* N2 */ { L009,       L09,     L001 },
     /* N3 */ { L009,       L09,     L001 },
     /*         Peaked TSS */
     /* P1 */ {  L01,      L045,     L045 },
     /* P2 */ {  L01,      L045,     L045 },
     /* P3 */ {  L01,      L045,     L045 },
    };

  /* Start states */
  const double start[N_STATES] =
    {
     L1,                        /* Background */
     L0, L0, L0,                /* Non-peaked TSS */
     L0, L0, L0                 /* Peaked TSS */
    };

  /* Initialize starting values. */
  for (int state = 0; state < N_STATES; ++state) {
     /* We don't care about the -1 sentinel value, because we never use it when
        choosing the Viterbi path.  It's there for flavor when inspecting the
        trellis in a debugger. */
    trellis->nodes[0][state].prev = -1;
    trellis->nodes[0][state].log_prob = start[state] + emis[state][obs[0]];
  }

  /* Run Viterbi. */
  for (int i = 1; i < len; ++i) {
    for (int state_cur = 0; state_cur < N_STATES; ++state_cur) {
      /* For a given state, choose the highest transition probability from the
         previous states.*/
      double prob_trans_max = (trellis->nodes[i-1][0].log_prob +
                               trans[0][state_cur]);
      int prev = 0;
      for (int state_prev = 1; state_prev < N_STATES; ++state_prev) {
        double prob_trans = (trellis->nodes[i-1][state_prev].log_prob +
                             trans[state_prev][state_cur]);
        if (prob_trans > prob_trans_max) {
          prob_trans_max = prob_trans;
          prev = state_prev;
        }
      }
      /* Add the emission probability. */
      trellis->nodes[i][state_cur].prev = prev;
      trellis->nodes[i][state_cur].log_prob = (prob_trans_max +
                                               emis[state_cur][obs[i]]);
    }
  }
}

void
viterbi_choose_path(int* ret, trellis_t* trellis, int len)
{
  /* Get the most probable state from the end. */
  double prob_max = 0.0 - INFINITY;
  int prev = -1;
  for (int state = 0; state < N_STATES; ++state)
    if (trellis->nodes[len-1][state].log_prob > prob_max) {
      prob_max = trellis->nodes[len-1][state].log_prob;
      prev = state;
    }
  ret[len-1] = prev;
  /* Backtrack for most probable states. */
  for (int i = len - 2; i >= 0; --i) {
    ret[i] = trellis->nodes[i+1][prev].prev;
    prev = trellis->nodes[i+1][prev].prev;
  }
}

void
viterbi(int* ret, int* obs, int len)
{
  trellis_t* trellis = NULL;
  trellis_init(&trellis, len);

  viterbi_fill_trellis(trellis, obs, len);
  viterbi_choose_path(ret, trellis, len);

  trellis_destroy(&trellis, len);
}
