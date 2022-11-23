/*  mileup2.c -- mpileup v2.0 API

   Copyright (C) 2022 Genome Research Ltd.

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
   THE SOFTWARE.

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <htslib/hts.h>
#include <sys/stat.h>
#include "mpileup.h"

static int mplp_add_bam(mpileup_t *mplp, char *bam, int is_file);

struct _mpileup_t
{
    char *fasta_ref, *regions, *regions_fname, *targets, *targets_fname, *samples, *samples_fname, *read_groups_fname;
    int smart_overlaps, smpl_ignore_rg, max_dp_per_file, adjust_mq, min_mq, min_bq, max_bq, delta_bq;
    int min_realn_dp, max_realn_dp, max_realn_len;
    float min_realn_frac;
    char **bam_names;
    int nbam_names;
};

mpileup_t *mpileup_init(void)
{
    mpileup_t *mplp = (mpileup_t*) calloc(1,sizeof(mpileup_t));
    mpileup_set_opt(mplp,int,MAX_DP_PER_FILE,250);
    mpileup_set_opt(mplp,int,MIN_BQ,1);
    mpileup_set_opt(mplp,int,MAX_BQ,60);
    mpileup_set_opt(mplp,int,DELTA_BQ,30);
    mpileup_set_opt(mplp,float,MIN_REALN_FRAC,0.05);
    mpileup_set_opt(mplp,int,MIN_REALN_DP,2);
    mpileup_set_opt(mplp,int,MAX_REALN_DP,250);
    mpileup_set_opt(mplp,int,MAX_REALN_LEN,500);
    return mplp;
}

int mpileup_set(mpileup_t *mplp, mpileup_opt_t key, void *value)
{
    switch (key)
    {
        case FASTA_REF: free(mplp->fasta_ref); mplp->fasta_ref = strdup(*((char**)value)); break;
        case BAM: return mplp_add_bam(mplp,*((char**)value),0); break;
        case BAM_FNAME: return mplp_add_bam(mplp,*((char**)value),1); break;
        case REGIONS: free(mplp->regions); mplp->regions = strdup(*((char**)value)); break;
        case REGIONS_FNAME: free(mplp->regions_fname); mplp->regions_fname = strdup(*((char**)value)); break;
        case TARGETS: free(mplp->targets); mplp->targets = strdup(*((char**)value)); break;
        case TARGETS_FNAME: free(mplp->targets_fname); mplp->targets_fname = strdup(*((char**)value)); break;
        case SAMPLES: free(mplp->samples); mplp->samples = strdup(*((char**)value)); break;
        case SAMPLES_FNAME: free(mplp->samples_fname); mplp->samples_fname = strdup(*((char**)value)); break;
        case READ_GROUPS_FNAME: free(mplp->read_groups_fname); mplp->read_groups_fname = strdup(*((char**)value)); break;
        case SMART_OVERLAPS: mplp->smart_overlaps = *((int*)value); break;
        case SMPL_IGNORE_RG: mplp->smpl_ignore_rg = *((int*)value); break;
        case MAX_DP_PER_FILE: mplp->max_dp_per_file = *((int*)value); break;
        case ADJUST_MQ: mplp->adjust_mq = *((int*)value); break;
        case MIN_MQ: mplp->min_mq = *((int*)value); break;
        case MIN_BQ: mplp->min_bq = *((int*)value); break;
        case MAX_BQ: mplp->max_bq = *((int*)value); break;
        case DELTA_BQ: mplp->delta_bq = *((int*)value); break;
        case MIN_REALN_FRAC: mplp->min_realn_frac = *((float*)value); break;
        case MIN_REALN_DP: mplp->min_realn_dp = *((int*)value); break;
        case MAX_REALN_DP: mplp->max_realn_dp = *((int*)value); break;
        case MAX_REALN_LEN: mplp->max_realn_len = *((int*)value); break;
        default: hts_log_error("Todo: mpileup_set key=%d",(int)key); return 1; break;
    }
    return 0;
}

void *mpileup_get(mpileup_t *mplp, mpileup_opt_t key)
{
    switch (key)
    {
        case SMART_OVERLAPS: return &mplp->smart_overlaps; break;
        case SMPL_IGNORE_RG: return &mplp->smpl_ignore_rg; break;
        case MAX_DP_PER_FILE: return &mplp->max_dp_per_file; break;
        case ADJUST_MQ: return &mplp->adjust_mq; break;
        case MIN_MQ: return &mplp->min_mq; break;
        case MIN_BQ: return &mplp->min_bq; break;
        case MAX_BQ: return &mplp->max_bq; break;
        case DELTA_BQ: return &mplp->delta_bq; break;
        case MIN_REALN_FRAC: return &mplp->min_realn_frac; break;
        case MIN_REALN_DP: return &mplp->min_realn_dp; break;
        case MAX_REALN_DP: return &mplp->max_realn_dp; break;
        case MAX_REALN_LEN: return &mplp->max_realn_len; break;
        default: hts_log_error("Todo: mpileup_get key=%d",(int)key); return NULL; break;
    }
    return NULL;
}

void mpileup_destroy(mpileup_t *mplp)
{
    int i;
    for (i=0; i<mplp->nbam_names; i++) free(mplp->bam_names[i]);
    free(mplp->bam_names);
    free(mplp->fasta_ref);
    free(mplp->regions);
    free(mplp->regions_fname);
    free(mplp->targets);
    free(mplp->targets_fname);
    free(mplp->samples);
    free(mplp->samples_fname);
    free(mplp->read_groups_fname);
    free(mplp);
}

static int is_url(const char *str)
{
    static const char uri_scheme_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+.-";
    return str[strspn(str, uri_scheme_chars)] == ':';
}

static int mplp_add_bam(mpileup_t *mplp, char *bam, int is_file)
{
    struct stat sbuf;
    if ( is_file )
    {
        int i, nnames = 0;
        char **names = hts_readlist(bam, 1, &nnames);
        for (i=0; i<nnames; i++)
        {
            if ( !(is_url(names[i]) || stat(names[i],&sbuf)==0) )
            {
                // no such file, check if it is safe to print its name
                char *tmp = names[i];
                while ( *tmp )
                {
                    if ( !isprint(*tmp) ) break;
                    tmp++;
                }
                if ( !*tmp )
                    hts_log_error("The file list \"%s\" appears broken, could not locate file: \"%s\"", bam,names[i]);
                else
                    hts_log_error("Does the file \"%s\" really contain a list of files and do all exist?", bam);
                break;
            }
            mplp->nbam_names++;
            char **tmp = (char**)realloc(mplp->bam_names,mplp->nbam_names*sizeof(*mplp->bam_names));
            if ( !tmp )
            {
                hts_log_error("Could not allocate %zu bytes",mplp->nbam_names*sizeof(*mplp->bam_names));
                break;
            }
            mplp->bam_names = tmp;
            mplp->bam_names[mplp->nbam_names-1] = strdup(names[i]);
        }
        for (i=0; i<nnames; i++) free(names[i]);
        free(names);
        return i==nnames ? 0 : 1;
    }
    mplp->nbam_names++;
    char **tmp = (char**)realloc(mplp->bam_names,mplp->nbam_names*sizeof(*mplp->bam_names));
    if ( !tmp )
    {
        hts_log_error("Could not allocate %zu bytes",mplp->nbam_names*sizeof(*mplp->bam_names));
        return 1;
    }
    mplp->bam_names = tmp;
    mplp->bam_names[mplp->nbam_names-1] = strdup(bam);
    return 0;
}

