#include "tss.h"

void
tss(int* indices_peak, int* groups, int* indices_signal, int* starts_signal,
    int* scores_signal, int len)
{
  int score_prev = -1;
  int group_prev = -1;
  for (int g = -1, i = 0; i < len; ++i) {
    int group = groups[i];
    /* Reset counters on reaching the next group.  Assume that the "groups"
       array from IRanges::findOverlaps() gives us contiguous group indices,
       but don't use the indices directly in case a particular entry is missing
       hits.
    */
    if (group_prev != group) {
      group_prev = group;
      ++g;
      indices_peak[g] = indices_signal[i];
      score_prev = scores_signal[i];
    }
    /* Calculate 2bp adjacency bonus. */
    int is_adjacent = 0;
    int bonus = 0;
    /* Look ahead. */
    if (i < len - 1 &&		    /* Boundary limit for the next check. */
	groups[i] == groups[i+1] && /* Compare only within the group. */
	starts_signal[i+1] - starts_signal[i] == 1) { /* Adjacency check. */
      is_adjacent = 1;
      bonus = scores_signal[i+1];
    }
    /* Look behind. */
    if (i > 0 &&		    /* Boundary limit for the next check. */
	groups[i-1] == groups[i] && /* Compare only within the group. */
	starts_signal[i] - starts_signal[i-1] == 1) { /* Adjacency check. */
      is_adjacent = 1;
      bonus += scores_signal[i-1];
    }
    if (! is_adjacent) {
      bonus = 0;
    }
    /* Find max value. */
    if (scores_signal[i] + bonus > score_prev) {
      indices_peak[g] = indices_signal[i];
      score_prev = scores_signal[i] + bonus;
    }
  }
}
