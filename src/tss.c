/** @file

    @brief Find peaks using 3 basepair tie-breaking.
 */

#include "tss.h"

/** Find peaks using 3 basepair tie-breaking.

    The original method for extracting the TSS was using 2 basepair windows,
    but note that this detail was not documented in Core 2014.  Preserving the
    intention of that method to use neighboring counts to break ties, this
    method uses 3 basepair overlapping windows, where the preceding and
    succeeding basepair counts act as a "bonus" value to help break ties.

    @param indices_peak Output integer vector to store peak locations.
    @param groups Sequential number indicating signal count group membership.
    @param indices_signal Absolute locations of signal counts.
    @param starts_signal Relative locations of signal counts relative to region.
    @param scores_signal Values of signal counts.
    @param len Number of members in signal group.
    @param prefer_last Boolean whether to effectively scan from the last count.
 */
void
tss(int* indices_peak, int* groups, int* indices_signal, int* starts_signal,
    int* scores_signal, int len, int* prefer_last)
{
  int score_prev = -1;
  int bonus_prev = 0;
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
      bonus_prev = 0;
    }
    /* Calculate 2bp adjacency bonus. */
    int bonus = 0;
    /* Look ahead. */
    if (i < len - 1 &&              /* Boundary limit for the next check. */
        groups[i] == groups[i+1] && /* Compare only within the group. */
        starts_signal[i+1] - starts_signal[i] == 1) { /* Adjacency check. */
      bonus = scores_signal[i+1];
    }
    /* Look behind. */
    if (i > 0 &&                    /* Boundary limit for the next check. */
        groups[i-1] == groups[i] && /* Compare only within the group. */
        starts_signal[i] - starts_signal[i-1] == 1) { /* Adjacency check. */
      bonus += scores_signal[i-1];
    }
    /* Find max value. */
    if (scores_signal[i] > score_prev ||
        ((scores_signal[i] == score_prev &&
	  (prefer_last[g] ?
	   bonus >= bonus_prev : /* Last peak on negative strand, or */
	   bonus > bonus_prev)	 /* first peak on positive strand. */
	  )))  {
      indices_peak[g] = indices_signal[i];
      score_prev = scores_signal[i];
      bonus_prev = bonus;
    }
  }
}
