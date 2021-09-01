/** @file

    @brief Find peaks using 3 basepair tie-breaking.
 */

#ifndef TSS_H
#define TSS_H

void tss(int* indices_peak, int* groups, int* indices_signal,
         int* starts_signal, int* scores_signal, int len, int* prefer_last);

#endif	/* TSS_H */
