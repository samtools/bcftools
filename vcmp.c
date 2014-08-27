/*  vcmp.c -- reference allele utility functions.

    Copyright (C) 2013 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include "vcmp.h"

struct _vcmp_t
{
    char *dref;
    int ndref, mdref;   // ndref: positive when ref1 longer, negative when ref2 is longer
    int nmatch;
};

vcmp_t *vcmp_init()
{
    return (vcmp_t*)calloc(1,sizeof(vcmp_t));
}

void vcmp_destroy(vcmp_t *vcmp)
{
    free(vcmp->dref);
    free(vcmp);
}

int vcmp_set_ref(vcmp_t *vcmp, char *ref1, char *ref2)
{
    vcmp->ndref = 0;

    char *a = ref1, *b = ref2;
    while ( *a && *b && *a==*b ) { a++; b++; }
    if ( !*a && !*b ) return 0;
    if ( *a && *b ) return -1;  // refs not compatible

    if ( *a )   // ref1 is longer
    {
        vcmp->nmatch = b-ref2;
        while ( *a ) a++;
        vcmp->ndref = (a-ref1) - vcmp->nmatch;
        hts_expand(char,vcmp->ndref+1,vcmp->mdref,vcmp->dref);
        memcpy(vcmp->dref,ref1+vcmp->nmatch,vcmp->ndref);
        vcmp->dref[vcmp->ndref] = 0;
        return 0;
    }

    // ref2 is longer
    vcmp->nmatch = a-ref1;
    while ( *b ) b++;
    vcmp->ndref = (b-ref2) - vcmp->nmatch;
    hts_expand(char,vcmp->ndref+1,vcmp->mdref,vcmp->dref);
    memcpy(vcmp->dref,ref2+vcmp->nmatch,vcmp->ndref);
    vcmp->dref[vcmp->ndref] = 0;
    vcmp->ndref *= -1;
    return 0;
}

int vcmp_find_allele(vcmp_t *vcmp, char **als1, int nals1, char *al2)
{
    int i, j;
    for (i=0; i<nals1; i++)
    {
        char *a = als1[i], *b = al2;
        while ( *a && *b && *a==*b ) { a++; b++; }
        if ( *a && *b ) continue;   // mismatch
        if ( !vcmp->ndref )
        {
            if ( !*a && !*b ) break;    // found
            continue;
        }

        // the prefixes match
        if ( *a )
        {
            if ( vcmp->ndref<0 ) continue;
            for (j=0; j<vcmp->ndref; j++)
                if ( !a[j] || a[j]!=vcmp->dref[j] ) break;
            if ( j!=vcmp->ndref || a[j] ) continue;
            break;  // found
        }

        if ( vcmp->ndref>0 ) continue;
        for (j=0; j<-vcmp->ndref; j++)
            if ( !b[j] || b[j]!=vcmp->dref[j] ) break;
        if ( j!=-vcmp->ndref || b[j] ) continue;
        break;  // found
    }
    if (i==nals1) return -1;
    return i;
}


