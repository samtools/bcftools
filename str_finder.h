#ifndef _STR_FINDER_H_
#define _STR_FINDER_H_

#include "utlist.h"

typedef struct rep_ele {
    int start, end, rep_len;
    struct rep_ele *prev;
    struct rep_ele *next;
} rep_ele;

/*
 * Finds repeated homopolymers up to 8-mers.
 *
 * If lower_only is true then it only adds STRs for regions that
 * contain at least one lower-case base. This can be used as a marker
 * for looking for specific types of repeats.
 * (One use for this is to only mark STRs that overlap a heterozygous
 * indel region.)
 *
 * Returns a list of rep_ele structs holding the start,end tuples of repeats;
 *         NULL on failure.
 */
rep_ele *find_STR(char *cons, int len, int lower_only);

/*
 * Returns an array of STR vs no-STR values.
 *         0  => non repetitive.
 *         1+ => repeat with consecutive bit-number for repeat size.
 *
 * Eg:  AGGGGAGGAGAAGAC
 *       1111  1111
 *         2222222
 *              444444
 * =>   011331137754440
 */
char *cons_mark_STR(char *cons, int len, int lower_only);

#endif /* _STR_FINDER_H_ */
