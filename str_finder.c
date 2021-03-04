#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>

#include "str_finder.h"
#include "utlist.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

typedef unsigned char uc;

//static int L[256] = {
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*   0-15 */
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  16 */
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  32 */
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  48 */
//    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, /*  64 */
//    0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /*  80 */
//    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, /*  96 */
//    0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 112-127 */
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 128 */
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  /* 255 */
//};

static void add_rep(rep_ele **list, char *cons, int clen, int pos, int rlen,
		    int lower_only, int w) {
    rep_ele *el, *tmp, *prev;
    char *cp1, *cp2, *cp_end;
    int i;

    // Already handled this in previous overlap?
    if (*list) {
	tmp = DL_TAIL(*list);
	if (tmp->start <= pos-rlen*2+1 && tmp->end >= pos)
	    return;
    }

    // Find current and last occurence of repeated word.

    cp2 = &cons[pos+1];
    // If unpadded, this is quicker: cp1 = &cons[pos+1-rlen];

    for (cp1 = &cons[pos], i = 1; i < rlen; cp1--) // compensate for pads
	if (*cp1 == '*')
	    continue;
	else
	    i++;
    while (*cp1 == '*')
	cp1--;


    // Scan ahead to see how much further it goes.
    cp_end = &cons[clen];
    while (cp2 < cp_end) {
	while (*cp1 == '*') cp1++;
	while (*cp2 == '*') cp2++;

	//if (L[(uc)*cp1] != L[(uc)*cp2])
	if (*cp1 != *cp2)
	    break;

	w<<=2;
	//w|=L[(uc)*cp2];
	w|=*cp2;
	cp1++;
	cp2++;
    }

    if (!(el = malloc(sizeof(*el))))
	return;

    el->end   = pos + cp2-&cons[pos+1];
    el->rep_len = rlen;
    pos++;
    while (rlen--) {
	while (cons[--pos] == '*');
	while (cons[--pos] == '*');
    }
    //pos++;
    while (pos > 1 && cons[pos-1] == '*') pos--;
    el->start = pos;

    // Check it meets the lower-case only criteria
    if (lower_only) {
	int lc = 0;
	for (i = el->start; i <= el->end; i++) {
	    if (islower(cons[i])) {
		lc = 1;
		break;
	    }
	}

	if (!lc) {
	    free(el);
	    return;
	}
    }

    // Remove any older items on the list that are entirely contained within el
    if (*list) {
	tmp = DL_TAIL(*list);
	do {
	    prev = tmp->prev;
	    if (tmp->end < el->start)
		break;

	    if (tmp->start >= el->start) {
		DL_DELETE(*list, tmp);
		free(tmp);
	    }

	    if (tmp == DL_HEAD(*list))
		break;
	    tmp = prev;
	} while (*list);
    }

    DL_APPEND(*list, el);

    return;
}

/*
 * Finds repeated homopolymers up to 8-mers.
 *
 * Returns a list of rep_ele structs holding the start,end tuples of repeats;
 *         NULL on failure.
 */
rep_ele *find_STR(char *cons, int len, int lower_only) {
    int i, j;
    uint32_t w = 0;
    rep_ele *reps = NULL;

    for (i = j = 0; i < len && j < 15; i++) {
	if (cons[i] == '*') continue;

	w <<= 2;
	//w |= L[(uc)cons[i]];;
	w |= cons[i];
	//printf("%3d %c w=%08x\n", i, cons[i], w);
	if (j>= 1 && (w&0x0003) == ((w>> 2)&0x0003))
	    add_rep(&reps, cons, len, i, 1, lower_only, w);
	if (j>= 3 && (w&0x000f) == ((w>> 4)&0x000f))
	    add_rep(&reps, cons, len, i, 2, lower_only, w);
	if (j>= 5 && (w&0x003f) == ((w>> 6)&0x003f))
	    add_rep(&reps, cons, len, i, 3, lower_only, w);
	if (j>= 7 && (w&0x00ff) == ((w>> 8)&0x00ff))
	    add_rep(&reps, cons, len, i, 4, lower_only, w);
	if (j>= 9 && (w&0x03ff) == ((w>>10)&0x03ff))
	    add_rep(&reps, cons, len, i, 5, lower_only, w);
	if (j>=11 && (w&0x0fff) == ((w>>12)&0x0fff))
	    add_rep(&reps, cons, len, i, 6, lower_only, w);
	if (j>=13 && (w&0x3fff) == ((w>>14)&0x3fff))
	    add_rep(&reps, cons, len, i, 7, lower_only, w);

	j++;
    }

    for (; i < len; i++) {	
	if (cons[i] == '*') continue;

	w <<= 2;
	//w |= L[(uc)cons[i]];
	w |= cons[i];
	//printf("%3d %c w=%08x\n", i, cons[i], w);
	if ((w&0xffff) == ((w>>16)&0xffff)) 
	    add_rep(&reps, cons, len, i, 8, lower_only, w);
	else if ((w&0x3fff) == ((w>>14)&0x3fff)) 
	    add_rep(&reps, cons, len, i, 7, lower_only, w);
	else if ((w&0x0fff) == ((w>>12)&0x0fff)) 
	    add_rep(&reps, cons, len, i, 6, lower_only, w);
	else if ((w&0x03ff) == ((w>>10)&0x03ff)) 
	    add_rep(&reps, cons, len, i, 5, lower_only, w);
	else if ((w&0x00ff) == ((w>> 8)&0x00ff)) 
	    add_rep(&reps, cons, len, i, 4, lower_only, w);
	else if ((w&0x003f) == ((w>> 6)&0x003f)) 
	    add_rep(&reps, cons, len, i, 3, lower_only, w);
	else if ((w&0x000f) == ((w>> 4)&0x000f)) 
	    add_rep(&reps, cons, len, i, 2, lower_only, w);
	else if ((w&0x0003) == ((w>> 2)&0x0003)) 
	    add_rep(&reps, cons, len, i, 1, lower_only, w);
    }

    return reps;
}

/* -----------------------------------------------------------------------------
 * Computes repeat regions in the consensus and then provides a bit mask
 * indicating the extend of the STRs.
 *
 * The purpose of this is to identify where a read needs to span the entire
 * region in order to validate how many copies of a repeat word are present.
 * This only really has a major impact when indels are involved.
 *
 * For example, given this multiple alignment:
 *
 * S1 GATCGGACGAGAG
 * S2 GATCGGACGAGAGAGAGAGAGT
 * S3 GATCGGACGAGAGAGAGAG**TCGGAC
 * S4     GGACGAGAGAGAGAGAGTCGGAC
 * S5        CGAGAGAGAGAG**TCGGAC
 * S6              AGAGAGAGTCGGAC
 *
 * We have subseq of GAGAGAGAGAG** vs GAGAGAGAGAGAG. The first and last
 * (S1 and S6) sequences do not span and so we do not know which allele they
 * match. Specifically as the pad is at the right hand end, the alignment of
 * S6 gives incorrect weight to the consensus as it is stating AG when it
 * may actually be ** at that point.
 *
 * By identifying the repeats we can soft clip as follows:
 *
 * S1 GATCGGACgagag
 * S2 GATCGGACGAGAGAGAGAGAGT
 * S3 GATCGGACGAGAGAGAGAG**TCGGAC
 * S4     GGACGAGAGAGAGAGAGTCGGAC
 * S5        CGAGAGAGAGAG**TCGGAC
 * S6              agagagagTCGGAC
 *
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
char *cons_mark_STR(char *cons, int len, int lower_only) {
    rep_ele *reps, *elt, *tmp;
    char *str;

    str = calloc(1, len);
    reps = find_STR(cons, len, lower_only);

    DL_FOREACH_SAFE(reps, elt, tmp) {
	int i, v = 0;
	
	//printf("%2d .. %2d %.*s\n", elt->start, elt->end,
	//       elt->end - elt->start+1, &cons[elt->start]);

	// What is there?
	for (i = MAX(elt->start-1,0); i <= MIN(elt->end+1,len-1); i++)
	    v |= str[i];

	for (i = 0; i < 8; i++) {
	    if (!(v&(1<<i)))
		break;
	}
	v = (i == 8) ? 1 : (1<<i);

	// Add new if available, or just overload 1 if not
	for (i = elt->start; i <= elt->end; i++)
	    str[i] |= v;

	DL_DELETE(reps, elt);
	free(elt);
    }

    return str;
}

#ifdef TEST_MAIN
int main(int argc, char **argv) {
    rep_ele *reps, *elt, *tmp;
    int len = strlen(argv[1]), pos = -1;
    if (argc > 2) pos = atoi(argv[2]);

    reps = find_STR(argv[1], len, 0);

    DL_FOREACH_SAFE(reps, elt, tmp) {
	if (pos != -1 &&
	    !(pos >= elt->start &&
	      pos <= elt->end)) {
	    DL_DELETE(reps, elt);
	    free(elt);
	    continue;
	}

	printf("%2d %2d .. %2d %.*s\n", elt->rep_len, elt->start, elt->end,
	       elt->end - elt->start+1, &argv[1][elt->start]);
	DL_DELETE(reps, elt);
	free(elt);
    }

    //char *str;
    //int i;
    //str = cons_mark_STR(argv[1], len, 1);
    //for (i = 0; i < len; i++) {
    //	printf("%3d %c %d\n", i, argv[1][i], str[i]);
    //}

    return 0;
}
#endif

#ifdef TEST_MAIN2
#define BS 1024*1024
static unsigned char *load(uint64_t *lenp) {
    unsigned char *data = NULL;
    uint64_t dsize = 0;
    uint64_t dcurr = 0;
    signed int len;

    do {
	if (dsize - dcurr < BS) {
	    unsigned char *tmp;
	    dsize = dsize ? dsize * 2 : BS;
	    if (!(tmp = realloc(data, dsize))) {
		free(data);
		return NULL;
	    }
	    data = tmp;
	}

	len = read(0, data + dcurr, BS);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
    }

    *lenp = dcurr;
    return data;
}

int main(int argc, char **argv) {
    rep_ele *reps, *elt, *tmp;
    char *str;
    uint64_t in_len;
    int count = 0;

    str = load(&in_len);

    reps = find_STR(str, in_len, 0);

    DL_FOREACH_SAFE(reps, elt, tmp) {
	printf("%2d .. %2d %.*s\n", elt->start, elt->end,
	       elt->end - elt->start+1, &str[elt->start]);
	count++;
	DL_DELETE(reps, elt);
	free(elt);
    }

    printf("Found %d reps\n", count);

    //uint64_t i;
    //str = cons_mark_STR(argv[1], in_len, 1);
    //for (i = 0; i < in_len; i++) {
    //	printf("%3d %c %d\n", i, argv[1][i], str[i]);
    //}

    return 0;
}
#endif
