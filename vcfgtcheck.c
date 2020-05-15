/*  vcfgtcheck.c -- Check sample identity.

    Copyright (C) 2013-2020 Genome Research Ltd.

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
#include <stdarg.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include <sys/time.h>
#include "bcftools.h"
#include "hclust.h"

typedef struct
{
    bcf_srs_t *files;           // first reader is the query VCF - single sample normally or multi-sample for cross-check
    bcf_hdr_t *gt_hdr, *qry_hdr; // VCF with genotypes to compare against and the query VCF
    char *cwd, **argv, *gt_samples, *qry_samples, *regions, *targets, *qry_fname, *gt_fname;
    int argc, gt_samples_is_file, qry_samples_is_file, regions_is_file, targets_is_file;
    int qry_use_GT,gt_use_GT, nqry_smpl,ngt_smpl, *qry_smpl,*gt_smpl;
    uint32_t *ndiff,*ncnt,ncmp, npairs;
    int32_t *qry_arr,*gt_arr, nqry_arr,ngt_arr;
    double *hwe_prob;
    double min_inter_err, max_intra_err;
    int all_sites, hom_only, ntop, cross_check, calc_hwe_prob, sort_by_hwe;
    FILE *fp;
}
args_t;

static void set_cwd(args_t *args)
{
    int i;
    char *buf;
    size_t nbuf = 500;
    args->cwd = (char*) malloc(sizeof(char)*nbuf);
    for (i=0; i<5; i++)
    {
        if ( (buf = getcwd(args->cwd, nbuf)) ) break;
        nbuf *= 2;
        args->cwd = (char*) realloc(args->cwd, sizeof(char)*nbuf);
    }
    assert(buf);
}
static void print_header(args_t *args, FILE *fp)
{
    fprintf(fp, "# This file was produced by bcftools (%s+htslib-%s), the command line was:\n", bcftools_version(), hts_version());
    fprintf(fp, "# \t bcftools %s ", args->argv[0]);
    int i;
    for (i=1; i<args->argc; i++)
        fprintf(fp, " %s",args->argv[i]);
    fprintf(fp, "\n# and the working directory was:\n");
    fprintf(fp, "# \t %s\n#\n", args->cwd);
}

static int cmp_int(const void *_a, const void *_b)
{
    int a = *((int*)_a);
    int b = *((int*)_b);
    if ( a < b ) return -1;
    if ( a > b ) return 1;
    return 0;
}

static void init_data(args_t *args)
{
    args->files = bcf_sr_init();
    if ( args->regions && bcf_sr_set_regions(args->files, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n", args->regions);
    if ( args->targets && bcf_sr_set_targets(args->files, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n", args->targets);

    if ( args->gt_fname ) bcf_sr_set_opt(args->files, BCF_SR_REQUIRE_IDX);
    if ( !bcf_sr_add_reader(args->files,args->qry_fname) ) error("Failed to open %s: %s\n", args->qry_fname,bcf_sr_strerror(args->files->errnum));
    if ( args->gt_fname && !bcf_sr_add_reader(args->files, args->gt_fname) )
        error("Failed to read from %s: %s\n", !strcmp("-",args->gt_fname)?"standard input":args->gt_fname,bcf_sr_strerror(args->files->errnum));

    args->qry_hdr = bcf_sr_get_header(args->files,0);
    if ( !bcf_hdr_nsamples(args->qry_hdr) ) error("No samples in %s?\n", args->qry_fname);
    if ( args->gt_fname )
    {
        args->gt_hdr = bcf_sr_get_header(args->files,1);
        if ( !bcf_hdr_nsamples(args->gt_hdr) ) error("No samples in %s?\n", args->gt_fname);
    }

    // Determine whether GT or PL will be used
    if ( args->qry_use_GT==-1 ) // not set by -u, qry uses PL by default
    {
        if ( bcf_hdr_id2int(args->qry_hdr,BCF_DT_ID,"PL")>=0 )
            args->qry_use_GT = 0;
        else if ( bcf_hdr_id2int(args->qry_hdr,BCF_DT_ID,"GT")>=0 )
            args->qry_use_GT = 1;
        else
            error("[E::%s] Neither PL nor GT tag is present in the header of %s\n", __func__, args->qry_fname);
    }
    else if ( args->qry_use_GT==0 && bcf_hdr_id2int(args->qry_hdr,BCF_DT_ID,"GT")<0 )
        error("[E::%s] The GT tag is not present in the header of %s\n", __func__, args->qry_fname);
    else if ( bcf_hdr_id2int(args->qry_hdr,BCF_DT_ID,"PL")<0 )
        error("[E::%s] The PL tag is not present in the header of %s\n", __func__, args->qry_fname);

    if ( args->gt_hdr )
    {
        if ( args->gt_use_GT==-1 ) // not set by -u, gt uses GT by default
        {
            if ( bcf_hdr_id2int(args->gt_hdr,BCF_DT_ID,"GT")>=0 )
                args->gt_use_GT = 1;
            else if ( bcf_hdr_id2int(args->gt_hdr,BCF_DT_ID,"PL")>=0 )
                args->gt_use_GT = 0;
            else
                error("[E::%s] Neither PL nor GT tag is present in the header of %s\n", __func__, args->gt_fname);
        }
        else if ( args->gt_use_GT==0 && bcf_hdr_id2int(args->gt_hdr,BCF_DT_ID,"GT")<0 )
            error("[E::%s] The GT tag is not present in the header of %s\n", __func__, args->gt_fname);
        else if ( bcf_hdr_id2int(args->gt_hdr,BCF_DT_ID,"PL")<0 )
            error("[E::%s] The PL tag is not present in the header of %s\n", __func__, args->gt_fname);
    }
    else
        args->gt_use_GT = args->qry_use_GT;

    // Prepare samples
    int i;
    args->nqry_smpl = bcf_hdr_nsamples(args->qry_hdr);
    if ( args->qry_samples )
    {
        char **tmp = hts_readlist(args->qry_samples, args->qry_samples_is_file, &args->nqry_smpl);
        if ( !tmp || !args->nqry_smpl ) error("Failed to parse %s\n", args->qry_samples);
        args->qry_smpl = (int*) malloc(sizeof(*args->qry_smpl)*args->nqry_smpl);
        for (i=0; i<args->nqry_smpl; i++)
        {
            int idx = bcf_hdr_id2int(args->qry_hdr, BCF_DT_SAMPLE, tmp[i]);
            if ( idx<0 ) error("No such sample in %s: [%s]\n",args->qry_fname,tmp[i]);
            // todo: add a check to prevent duplicates
            args->qry_smpl[i] = idx;
            free(tmp[i]);
        }
        free(tmp);
        qsort(args->qry_smpl,args->nqry_smpl,sizeof(*args->qry_smpl),cmp_int);
    }
    if ( args->gt_hdr )
    {
        args->ngt_smpl = bcf_hdr_nsamples(args->gt_hdr);
        if ( args->gt_samples )
        {
            char **tmp = hts_readlist(args->gt_samples, args->gt_samples_is_file, &args->ngt_smpl);
            if ( !tmp || !args->ngt_smpl ) error("Failed to parse %s\n", args->gt_samples);
            args->gt_smpl = (int*) malloc(sizeof(*args->gt_smpl)*args->ngt_smpl);
            for (i=0; i<args->ngt_smpl; i++)
            {
                int idx = bcf_hdr_id2int(args->gt_hdr, BCF_DT_SAMPLE, tmp[i]);
                if ( idx<0 ) error("No such sample in %s: [%s]\n",args->gt_fname,tmp[i]);
                // todo: add a check to prevent duplicates
                args->gt_smpl[i] = idx;
                free(tmp[i]);
            }
            free(tmp);
            qsort(args->gt_smpl,args->ngt_smpl,sizeof(*args->gt_smpl),cmp_int);
        }
    }
    else if ( args->gt_samples )
    {
        char **tmp = hts_readlist(args->gt_samples, args->gt_samples_is_file, &args->ngt_smpl);
        if ( !tmp || !args->ngt_smpl ) error("Failed to parse %s\n", args->gt_samples);
        args->gt_smpl = (int*) malloc(sizeof(*args->gt_smpl)*args->ngt_smpl);
        for (i=0; i<args->ngt_smpl; i++)
        {
            int idx = bcf_hdr_id2int(args->qry_hdr, BCF_DT_SAMPLE, tmp[i]);
            if ( idx<0 ) error("No such sample in %s: [%s]\n",args->gt_fname,tmp[i]);
            // todo: add a check to prevent duplicates
            args->gt_smpl[i] = idx;
            free(tmp[i]);
        }
        free(tmp);
        qsort(args->gt_smpl,args->ngt_smpl,sizeof(*args->gt_smpl),cmp_int);
    }
    else
    {
        args->ngt_smpl = args->nqry_smpl;
        args->gt_smpl  = args->qry_smpl;
        args->cross_check = 1;
    }

    // The data arrays
    args->npairs = args->cross_check ? args->nqry_smpl*(args->nqry_smpl+1)/2 : args->ngt_smpl*args->nqry_smpl;
    args->ndiff = (uint32_t*) calloc(args->npairs,sizeof(*args->ndiff));    // number of differing genotypes for each pair of samples
    args->ncnt  = (uint32_t*) calloc(args->npairs,sizeof(*args->ncnt));     // number of comparisons performed (non-missing data)
    if ( !args->ncnt ) error("Error: failed to allocate %.1f Mb\n", args->npairs*sizeof(*args->ncnt)/1e6);
    if ( args->calc_hwe_prob )
    {
        // prob of the observed sequence of matches given site AFs and HWE
        args->hwe_prob = (double*) calloc(args->npairs,sizeof(*args->hwe_prob));
        if ( !args->hwe_prob ) error("Error: failed to allocate %.1f Mb. Run with --no-HWE-prob to save some memory.\n", args->npairs*sizeof(*args->hwe_prob)/1e6);
    }

    args->fp = stdout;
    print_header(args, args->fp);
}

static void destroy_data(args_t *args)
{
    free(args->hwe_prob);
    free(args->cwd);
    free(args->qry_arr);
    if ( args->gt_hdr ) free(args->gt_arr);
    free(args->ndiff);
    free(args->ncnt);
    free(args->qry_smpl);
    if ( args->gt_smpl!=args->qry_smpl ) free(args->gt_smpl);
    bcf_sr_destroy(args->files);
}

/*
   Return -1 on missing data, 0 on mismatch, 1 on match.
   Note that:
        - currently only diploid, non-missing values are considered
        - with PLs we only compare whether the most likely PL=P(D|G) values match
*/
#define _SLOWER_BRANCH 0
#if _SLOWER_BRANCH
static inline int match_GT_GT(int32_t *aptr, int32_t *bptr, int *adsg, int *bdsg)
{
    if ( bcf_gt_is_missing(aptr[0]) || bcf_gt_is_missing(aptr[1]) || aptr[1]==bcf_int32_vector_end ) return -1;
    if ( bcf_gt_is_missing(bptr[0]) || bcf_gt_is_missing(bptr[1]) || bptr[1]==bcf_int32_vector_end ) return -1;
    *adsg = (bcf_gt_allele(aptr[0])?1:0) + (bcf_gt_allele(aptr[1])?1:0);
    *bdsg = (bcf_gt_allele(bptr[0])?1:0) + (bcf_gt_allele(bptr[1])?1:0);
    return *adsg==*bdsg ? 1 : 0;
}
static inline int match_GT_PL(int32_t *aptr, int32_t *bptr, int *adsg, int *bdsg)
{
    if ( bcf_gt_is_missing(aptr[0]) || bcf_gt_is_missing(aptr[1]) || aptr[1]==bcf_int32_vector_end ) return -1;
    if ( bptr[0]==bcf_int32_missing || bptr[1]==bcf_int32_missing || bptr[2]==bcf_int32_missing || bptr[1]==bcf_int32_vector_end || bptr[2]==bcf_int32_vector_end ) return -1;
    *adsg = (bcf_gt_allele(aptr[0])?1:0) + (bcf_gt_allele(aptr[1])?1:0);
    *bdsg = 0;
    int i, min = bptr[0];
    for (i=1; i<3; i++)
        if ( min > bptr[i] ) { min = bptr[i]; *bdsg = i; }
    return min==bptr[*adsg] ? 1 : 0;
}
static inline int match_PL_PL(int32_t *aptr, int32_t *bptr, int *adsg, int *bdsg)
{
    if ( aptr[0]==bcf_int32_missing || aptr[1]==bcf_int32_missing || aptr[2]==bcf_int32_missing || aptr[1]==bcf_int32_vector_end || aptr[2]==bcf_int32_vector_end ) return -1;
    if ( bptr[0]==bcf_int32_missing || bptr[1]==bcf_int32_missing || bptr[2]==bcf_int32_missing || bptr[1]==bcf_int32_vector_end || bptr[2]==bcf_int32_vector_end ) return -1;
    int i, amin = aptr[0], bmin = bptr[0];
    *adsg = 0; *bdsg = 0;
    for (i=1; i<3; i++)
        if ( amin > aptr[i] ) { amin = aptr[i]; *adsg = i; }
    for (i=1; i<3; i++)
        if ( bmin > bptr[i] ) { bmin = bptr[i]; *bdsg = i; }
    for (i=0; i<3; i++)
        if ( aptr[i]==amin && bptr[i]==bmin ) { *adsg = *bdsg = i; return 1; }
    return 0;
}
#else   /* faster branch for missing data */
#define HAS_GT(ptr) (!bcf_gt_is_missing(ptr[0]) && !bcf_gt_is_missing(ptr[1]) && ptr[1]!=bcf_int32_vector_end)
#define HAS_PL(ptr) (ptr[0]!=bcf_int32_missing && ptr[1]!=bcf_int32_missing && ptr[2]!=bcf_int32_missing && ptr[1]!=bcf_int32_vector_end && ptr[2]!=bcf_int32_vector_end)
#define MIN_PL(ptr) ptr[0]<ptr[1]?(ptr[0]<ptr[2]?ptr[0]:ptr[2]):(ptr[1]<ptr[2]?ptr[1]:ptr[2])
#define DSG_PL(ptr) ptr[0]<ptr[1]?(ptr[0]<ptr[2]?0:2):(ptr[1]<ptr[2]?1:2)
#define DSG_GT(ptr) (bcf_gt_allele(ptr[0])?1:0) + (bcf_gt_allele(ptr[1])?1:0)
#endif

static void process_line(args_t *args)
{
    int nqry1, ngt1;

    bcf1_t *gt_rec, *qry_rec = bcf_sr_get_line(args->files,0);   // the query file
    if ( args->qry_use_GT )
    {
        if ( (nqry1=bcf_get_genotypes(args->qry_hdr,qry_rec,&args->qry_arr,&args->nqry_arr)) <= 0 ) return;
        if ( nqry1 != 2*bcf_hdr_nsamples(args->qry_hdr) ) return;    // only diploid data for now
        nqry1 = 2;
    }
    else
    {
        if ( (nqry1=bcf_get_format_int32(args->qry_hdr,qry_rec,"PL",&args->qry_arr,&args->nqry_arr)) <= 0 ) return;
        if ( nqry1 != 3*bcf_hdr_nsamples(args->qry_hdr) ) return;    // not diploid
        nqry1 = 3;
    }

    if ( args->gt_hdr )
    {
        gt_rec = bcf_sr_get_line(args->files,1);
        if ( args->gt_use_GT )
        {
            if ( (ngt1=bcf_get_genotypes(args->gt_hdr,gt_rec,&args->gt_arr,&args->ngt_arr)) <= 0 ) return;
            if ( ngt1 != 2*bcf_hdr_nsamples(args->gt_hdr) ) return;    // not diploid
            ngt1 = 2;
        }
        else
        {
            if ( (ngt1=bcf_get_format_int32(args->gt_hdr,gt_rec,"PL",&args->gt_arr,&args->ngt_arr)) <= 0 ) return;
            if ( ngt1 != 3*bcf_hdr_nsamples(args->gt_hdr) ) return;    // not diploid
            ngt1 = 3;
        }
    }
    else
    {
        ngt1 = nqry1;
        args->gt_arr = args->qry_arr;
    }

    int ac[2];
    if ( args->gt_hdr )
    {
        if ( bcf_calc_ac(args->gt_hdr, gt_rec, ac, BCF_UN_INFO|BCF_UN_FMT)!=1 ) error("todo: bcf_calc_ac() failed\n");
    }
    else if ( bcf_calc_ac(args->qry_hdr, qry_rec, ac, BCF_UN_INFO|BCF_UN_FMT)!=1 ) error("todo: bcf_calc_ac() failed\n");

    args->ncmp++;

    double af,hwe[3];
    if ( args->calc_hwe_prob )
    {
        const double min_af = 1e-3;
        af = (double)ac[1]/(ac[0]+ac[1]);
        hwe[0] = af>min_af ? -log(af*af) : -log(min_af*min_af);
        hwe[1] = af>min_af && af<1-min_af ? -log(2*af*(1-af)) : -log(2*min_af*(1-min_af));
        hwe[2] = af<1-min_af ? -log((1-af)*(1-af)) : -log(min_af*min_af);
    }

#if _SLOWER_BRANCH
    int i,j,idx=0;
    for (i=0; i<args->nqry_smpl; i++)
    {
        int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
        int ngt  = args->cross_check ? i : args->ngt_smpl;     // two files or a sub-diagnoal cross-check mode?
        for (j=0; j<ngt; j++)
        {
            int igt = args->gt_smpl ? args->gt_smpl[j] : j;
            int32_t *aptr = args->qry_arr + iqry*nqry1;
            int32_t *bptr = args->gt_arr + igt*ngt1;
            int match, qry_dsg, ign;
            if ( args->qry_use_GT && args->gt_use_GT )
                match = match_GT_GT(aptr,bptr,&qry_dsg,&ign);
            else if ( !args->qry_use_GT && !args->gt_use_GT )
                match = match_PL_PL(aptr,bptr,&qry_dsg,&ign);
            else if ( args->qry_use_GT )
                match = match_GT_PL(aptr,bptr,&qry_dsg,&ign);
            else
                match = match_GT_PL(bptr,aptr,&ign,&qry_dsg);
            if ( match>=0 ) 
            {
                if ( !match ) args->ndiff[idx]++;
                else if ( args->calc_hwe_prob ) args->hwe_prob[idx] += hwe[qry_dsg];
                args->ncnt[idx]++;
            }
            idx++;
        }
    }
#else
    int i,j, idx = 0;
    for (i=0; i<args->nqry_smpl; i++)
    {
        int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
        int ngt  = args->cross_check ? i : args->ngt_smpl;     // two files or a sub-diagnoal cross-check mode?
        int32_t *aptr = args->qry_arr + iqry*nqry1;
        int32_t aval, qry_dsg;
        if ( args->qry_use_GT )
        {
            if ( !HAS_GT(aptr) ) { idx += ngt; continue; }
            aval = qry_dsg = DSG_GT(aptr);
        }
        else
        {
            if ( !HAS_PL(aptr) ) { idx += ngt; continue; }
            aval = MIN_PL(aptr);
            qry_dsg = DSG_PL(aptr);
        }

        for (j=0; j<ngt; j++)
        {
            int igt = args->gt_smpl ? args->gt_smpl[j] : j;
            int32_t *bptr = args->gt_arr + igt*ngt1;
            int32_t bval;
            if ( args->gt_use_GT )
            {
                if ( !HAS_GT(bptr) ) { idx++; bptr[0] = bcf_gt_missing; continue; }
                bval = DSG_GT(bptr);
            }
            else
            {
                if ( !HAS_PL(bptr) ) { idx++; bptr[0] = bcf_int32_missing; continue; }
                bval = MIN_PL(bptr);
            }

            int match;
            if ( args->qry_use_GT )
            {
                if ( args->gt_use_GT )
                    match = aval==bval ? 1 : 0;
                else
                    match = bptr[aval]==bval ? 1 : 0;
            }
            else
            {
                if ( args->gt_use_GT )
                    match = aptr[bval]==aval ? 1 : 0;
                else
                {
                    int k;
                    match = 0;
                    for (k=0; k<3; k++)
                        if ( aptr[k]==aval && bptr[k]==bval ) { match = 1; break; }
                }
            }
            if ( !match ) args->ndiff[idx]++;
            else if ( args->calc_hwe_prob ) args->hwe_prob[idx] += hwe[qry_dsg];
            args->ncnt[idx]++;
            idx++;
        }
    }
#endif
}


typedef struct
{
    int ism, idx;
    double val;
}
idbl_t;
static int cmp_idbl(const void *_a, const void *_b)
{
    idbl_t *a = (idbl_t*)_a;
    idbl_t *b = (idbl_t*)_b;
    if ( a->val < b->val ) return -1;
    if ( a->val > b->val ) return 1;
    return 0;
}
static void report(args_t *args)
{
    fprintf(args->fp,"# DC, discordance:\n");
    fprintf(args->fp,"#     - query sample\n");
    fprintf(args->fp,"#     - genotyped sample\n");
    fprintf(args->fp,"#     - discordance (number of mismatches; smaller is better)\n");
    fprintf(args->fp,"#     - negative log of HWE probability at matching sites (bigger is better)\n");
    fprintf(args->fp,"#     - number of sites compared (bigger is better)\n");
    fprintf(args->fp,"#DC\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]-log P(HWE)\t[6]Number of sites compared\n");

    int trim = args->ntop;
    if ( !args->ngt_smpl && args->nqry_smpl <= args->ntop ) trim = 0;
    if ( args->ngt_smpl && args->ngt_smpl <= args->ntop  ) trim = 0;

    if ( !trim )
    {
        int i,j,idx=0;
        for (i=0; i<args->nqry_smpl; i++)
        {
            int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
            int ngt  = args->cross_check ? i : args->ngt_smpl;
            for (j=0; j<ngt; j++)
            {
                int igt = args->gt_smpl ? args->gt_smpl[j] : j;
                fprintf(args->fp,"DC\t%s\t%s\t%u\t%e\t%u\n",
                        args->qry_hdr->samples[iqry],
                        args->gt_hdr?args->gt_hdr->samples[igt]:args->qry_hdr->samples[igt],
                        args->ndiff[idx],
                        args->calc_hwe_prob ? args->hwe_prob[idx] : 0,
                        args->ncnt[idx]);
                idx++;
            }
        }
    }
    else if ( !args->cross_check )
    {
        idbl_t *arr = (idbl_t*)malloc(sizeof(*arr)*args->ngt_smpl);
        int i,j;
        for (i=0; i<args->nqry_smpl; i++)
        {
            int idx  = i*args->ngt_smpl;
            for (j=0; j<args->ngt_smpl; j++)
            {
                if ( args->sort_by_hwe )
                    arr[j].val = -args->hwe_prob[idx];
                else
                    arr[j].val = args->ncnt[idx] ? (double)args->ndiff[idx]/args->ncnt[idx] : 0;
                arr[j].ism = j;
                arr[j].idx = idx;
                idx++;
            }
            qsort(arr, args->ngt_smpl, sizeof(*arr), cmp_idbl);
            int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
            for (j=0; j<args->ntop; j++)
            {
                int idx = arr[j].idx;
                int igt = args->gt_smpl ? args->gt_smpl[arr[j].ism] : arr[j].ism;
                fprintf(args->fp,"DC\t%s\t%s\t%u\t%e\t%u\n",
                        args->qry_hdr->samples[iqry],
                        args->gt_hdr->samples[igt],
                        args->ndiff[idx],
                        args->calc_hwe_prob ? args->hwe_prob[idx] : 0,
                        args->ncnt[idx]);
            }
        }
        free(arr);
    }
    else
    {
        int narr = args->nqry_smpl-1;
        idbl_t *arr = (idbl_t*)malloc(sizeof(*arr)*narr);
        int i,j,k,idx;
        for (i=0; i<args->nqry_smpl; i++)
        {
            k = 0, idx = i*(i-1)/2;
            for (j=0; j<i; j++)
            {
                if ( args->sort_by_hwe )
                    arr[k].val = -args->hwe_prob[idx];
                else
                    arr[k].val = args->ncnt[idx] ? (double)args->ndiff[idx]/args->ncnt[idx] : 0;
                arr[k].ism = j;
                arr[k].idx = idx;
                idx++;
                k++;
            }
            for (; j<narr; j++)
            {
                idx = j*(j+1)/2 + i;
                if ( args->sort_by_hwe )
                    arr[k].val = -args->hwe_prob[idx];
                else
                    arr[k].val = args->ncnt[idx] ? (double)args->ndiff[idx]/args->ncnt[idx] : 0;
                arr[k].ism = j + 1;
                arr[k].idx = idx;
                k++;
            }
            qsort(arr, narr, sizeof(*arr), cmp_idbl);
            int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
            for (j=0; j<args->ntop; j++)
            {
                if ( i <= arr[j].ism ) continue;
                int idx = arr[j].idx;
                int igt = args->qry_smpl ? args->qry_smpl[arr[j].ism] : arr[j].ism;
                fprintf(args->fp,"DC\t%s\t%s\t%u\t%e\t%u\n",
                        args->qry_hdr->samples[iqry],
                        args->qry_hdr->samples[igt],
                        args->ndiff[idx],
                        args->calc_hwe_prob ? args->hwe_prob[idx] : 0,
                        args->ncnt[idx]);
            }
        }
        free(arr);
    }

    fclose(args->fp);
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Check sample identity. With no -g BCF given, multi-sample cross-check is performed.\n");
    fprintf(stderr, "Usage:   bcftools gtcheck [options] [-g <genotypes.vcf.gz>] <query.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -a, --all-sites                    output comparison for all sites\n");
    fprintf(stderr, "    -c, --cluster MIN,MAX              min inter- and max intra-sample error [0.23,-0.3]\n");
    fprintf(stderr, "    -g, --genotypes FILE               genotypes to compare against\n");
    fprintf(stderr, "    -H, --homs-only                    homozygous genotypes only (useful for low coverage data)\n");
    fprintf(stderr, "        --n-matches INT                print only top INT matches for each sample, 0 for unlimited. Use negative value\n");
    fprintf(stderr, "                                            to sort by HWE probability rather than the number of discordant sites [0]\n");
    fprintf(stderr, "        --no-HWE-prob                  disable calculation of HWE probability\n");
    fprintf(stderr, "    -r, --regions REGION               restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file FILE            restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --samples [qry|gt]:LIST        list of query or -g samples (by default all samples are compared)\n");
    fprintf(stderr, "    -S, --samples-file [qry|gt]:FILE   file with the query or -g samples to compare\n");
    fprintf(stderr, "    -t, --targets REGION               similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file FILE            similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "    -u, --use TAG1[,TAG2]              which tag to use in the query file (TAG1) and the -g (TAG2) files [PL,GT]\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "   # Are there any matching samples in file A and B?\n");
    fprintf(stderr, "   bcftools gtcheck -g A.bcf B.bcf > out.txt\n");
    fprintf(stderr, "   cat out.txt | ... \n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfgtcheck(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv; set_cwd(args);
    args->qry_use_GT = -1;
    args->gt_use_GT  = -1;
    args->calc_hwe_prob = 1;

    // In simulated sample swaps the minimum error was 0.3 and maximum intra-sample error was 0.23
    //    - min_inter: pairs with smaller err value will be considered identical 
    //    - max_intra: pairs with err value bigger than abs(max_intra_err) will be considered
    //                  different. If negative, the cutoff may be heuristically lowered
    args->min_inter_err =  0.23;
    args->max_intra_err = -0.3;

    static struct option loptions[] =
    {
        {"use",1,0,'u'},
        {"cluster",1,0,'c'},
        {"GTs-only",1,0,'G'},
        {"all-sites",0,0,'a'},
        {"homs-only",0,0,'H'},
        {"help",0,0,'h'},
        {"genotypes",1,0,'g'},
        {"plot",1,0,'p'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"n-matches",1,0,2},
        {"no-HWE-prob",0,0,3},
        {"target-sample",1,0,4},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {0,0,0,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "hg:p:s:S:Hr:R:at:T:G:c:u:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'u':
                {
                    int i,nlist;
                    char **list = hts_readlist(optarg, 0, &nlist);
                    if ( !list || nlist<=0 || nlist>2 ) error("Failed to parse --use %s\n", optarg);
                    if ( !strcasecmp("GT",list[0]) ) args->qry_use_GT = 1;
                    else if ( !strcasecmp("PL",list[0]) ) args->qry_use_GT = 0;
                    else error("Failed to parse --use %s; only GT and PL are supported\n", optarg);
                    if ( nlist==2 )
                    {
                        if ( !strcasecmp("GT",list[1]) ) args->gt_use_GT = 1;
                        else if ( !strcasecmp("PL",list[1]) ) args->gt_use_GT = 0;
                        else error("Failed to parse --use %s; only GT and PL are supported\n", optarg);
                    }
                    else args->gt_use_GT = args->qry_use_GT;
                    for (i=0; i<nlist; i++) free(list[i]);
                    free(list);
                }
                break;
            case 2 :
                args->ntop = strtol(optarg,&tmp,10);
                if ( !tmp || *tmp ) error("Could not parse: --n-matches %s\n", optarg);
                if ( args->ntop < 0 )
                {
                    args->sort_by_hwe = 1;
                    args->ntop *= -1;
                }
                break;
            case 3 : args->calc_hwe_prob = 0; break;
            case 4 : error("The option -S, --target-sample has been deprecated\n"); break;
            case 'c':
                args->min_inter_err = strtod(optarg,&tmp);
                if ( *tmp )
                {
                    if ( *tmp!=',') error("Could not parse: -c %s\n", optarg);
                    args->max_intra_err = strtod(tmp+1,&tmp);
                    if ( *tmp ) error("Could not parse: -c %s\n", optarg);
                }
                break;
            case 'G': error("The option -G, --GTs-only has been deprecated\n"); break;
            case 'a': args->all_sites = 1; break;
            case 'H': args->hom_only = 1; error("todo: -H\n"); break;
            case 'g': args->gt_fname = optarg; break;
//            case 'p': args->plot = optarg; break;
            case 's':
                if ( !strncasecmp("gt:",optarg,3) ) args->gt_samples = optarg+3;
                else if ( !strncasecmp("qry:",optarg,4) ) args->qry_samples = optarg+4;
                else error("Which one? Query samples (qry:%s) or genotype samples (gt:%s)?\n",optarg,optarg);
                break;
            case 'S': 
                if ( !strncasecmp("gt:",optarg,3) ) args->gt_samples = optarg+3, args->gt_samples_is_file = 1;
                else if ( !strncasecmp("qry:",optarg,4) ) args->qry_samples = optarg+4, args->qry_samples_is_file = 1;
                else error("Which one? Query samples (qry:%s) or genotype samples (gt:%s)?\n",optarg,optarg);
                break;
            case 'r': args->regions = optarg; break;
            case 'R': args->regions = optarg; args->regions_is_file = 1; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; args->targets_is_file = 1; break;
            case 'h':
            case '?': usage(); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->qry_fname = "-";  // reading from stdin
        else usage();   // no files given
    }
    else args->qry_fname = argv[optind];
    if ( argc>optind+1 )  usage();  // too many files given

    init_data(args);

    int ret;
    while ( (ret=bcf_sr_next_line(args->files)) )
    {
        if ( args->gt_hdr && ret!=2 ) continue;     // not a cross-check mode and lines don't match

        // time one record to give the user an estimate with very big files
        struct timeval t0, t1;
        if ( !args->ncmp )  gettimeofday(&t0, NULL);

        process_line(args);

        if ( args->ncmp==1 )
        {
            gettimeofday(&t1, NULL);
            double delta = (t1.tv_sec - t0.tv_sec) * 1e6 + (t1.tv_usec - t0.tv_usec);
            fprintf(stderr,"INFO:\tTime required to process one record .. %f seconds\n",delta/1e6);
            fprintf(args->fp,"INFO\tTime required to process one record .. %f seconds\n",delta/1e6);
        }
    }
    report(args);

    destroy_data(args);
    free(args);
    return 0;
}

