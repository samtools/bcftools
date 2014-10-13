/* 
    Copyright (C) 2014 Genome Research Ltd.

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
#include <getopt.h>
#include <stdarg.h>
#include <stdint.h>
#include <htslib/vcf.h>
#include <htslib/regidx.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "ploidy.h"

typedef struct
{
    int nsites, nsex, *sex2ploidy, dflt_ploidy, max_ploidy, guess;
    int ncounts, *counts, nsample, verbose;
    float *sex2prob, min_hets;
    int32_t *gts, ngts;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
    ploidy_t *ploidy;
}
args_t;

const char *about(void)
{
    return "Determine sample sex by checking genotypes in haploid regions.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Determine sample sex by either checking the presence of haploid/diploid\n"
        "       genotypes (requires correct ploidy in the VCF) or by counting homs/hets\n"
        "       in haploid regions.\n"
        "Usage: bcftools +vcf2sex <file.vcf.gz> -- [Plugin Options]\n"
        "Plugin options:\n"
        "   -g, --guess             do not trust genotypes ploidy, count hom/hets\n"
        "   -m, --min-hets <float>  minimum fraction of hets in diploid regions [0.05]\n"
        "   -n, --nsites <int>      number of sites to check per region (ignored with -g) [10]\n"
        "   -p, --ploidy <file>     space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY\n"
        "\n"
        "Example:\n"
        "   # Default ploidy, if -p not given. Unlisted regions have ploidy 2\n"
        "   X 1 60000 M 1\n"
        "   X 2699521 154931043 M 1\n"
        "   Y 1 59373566 M 1\n"
        "   Y 1 59373566 F 0\n"
        "   \n"
        "   bcftools +vcf2sex in.vcf.gz\n"
        "   bcftools +vcf2sex in.vcf.gz -- -n 10\n"
        "\n";
}

int process_region_precise(args_t *args, char *seq, regitr_t *itr)
{
    int k = 1;
    uint32_t start = itr->reg[itr->i].start, end = itr->reg[itr->i].end;
    while ( itr->i+k<itr->n && start==itr->reg[itr->i+k].start && end==itr->reg[itr->i+k].end ) k++;
    
    int ret = ploidy_query(args->ploidy, seq, start, args->sex2ploidy, NULL, NULL);
    assert(ret);

    memset(args->counts,0,args->ncounts*sizeof(int));

    // Select 'nsites' sites spaced so that they evenly cover the whole region 
    // to get a representative sample. We index-jump as we should be checking
    // a few sites only.
    int i, rid = -1, pos, prev_pos = -1, ismpl;
    for (i=0; i<args->nsites; i++)
    {
        rid = -1;
        pos = ((i+1.0)/(args->nsites+1))*(end - start) + start;
        if ( i>0 && pos <= prev_pos ) continue;     // the vcf is too sparse
        if ( bcf_sr_seek(args->sr,seq,pos)!=0 ) return k;   // sequence not present
        if ( !bcf_sr_next_line(args->sr) ) return k;        // no sites found
        bcf1_t *rec = bcf_sr_get_line(args->sr,0);
        if ( rid==-1 ) rid = rec->rid;
        if ( rid!=rec->rid || rec->pos > end ) break;
        prev_pos = rec->pos;

        int ngts = bcf_get_genotypes(args->hdr,rec,&args->gts,&args->ngts);
        ngts /= args->nsample;
        for (ismpl=0; ismpl<args->nsample; ismpl++)
        {
            int32_t *gts = args->gts + ngts*ismpl;
            int igt, ploidy = 0;
            for (igt=0; igt<ngts; igt++)
                if ( gts[igt]==bcf_gt_missing || gts[igt]==bcf_int32_missing || gts[igt]==bcf_int32_vector_end ) break;
                else ploidy++;
            args->counts[ismpl*(args->max_ploidy+1) + ploidy]++;
            if ( args->verbose )
                fprintf(stderr,"%s:%d\t%s\tploidy=%d\n", seq,rec->pos+1,args->hdr->samples[ismpl],ploidy);
        }
    }

    for (ismpl=0; ismpl<args->nsample; ismpl++)
    {
        float sum = 0, *probs = args->sex2prob + ismpl*args->nsex;
        int *counts = args->counts + ismpl*(args->max_ploidy+1);
        for (i=0; i<args->max_ploidy+1; i++) sum += counts[i];
        if ( !sum ) continue;
        for (i=0; i<args->nsex; i++)
        {
            int ploidy = args->sex2ploidy[i];
            probs[i] *= counts[ploidy]/sum;
        }
    }

    return k;
}

int process_region_guess(args_t *args, char *seq, regitr_t *itr)
{
    int ismpl, k = 1;
    uint32_t start = itr->reg[itr->i].start, end = itr->reg[itr->i].end;
    while ( itr->i+k<itr->n && start==itr->reg[itr->i+k].start && end==itr->reg[itr->i+k].end ) k++;
    
    int ret = ploidy_query(args->ploidy, seq, start, args->sex2ploidy, NULL, NULL);
    assert(ret);

    memset(args->counts,0,args->ncounts*sizeof(int));

    if ( bcf_sr_seek(args->sr,seq,start)!=0 ) return k;   // sequence not present
    int rid = bcf_hdr_name2id(args->hdr,seq);
    while ( bcf_sr_next_line(args->sr) )
    {
        bcf1_t *rec = bcf_sr_get_line(args->sr,0);
        if ( rec->rid!=rid || rec->pos > end ) break;

        bcf_fmt_t *fmt = bcf_get_fmt(args->hdr, rec, "GT");
        for (ismpl=0; ismpl<args->nsample; ismpl++)
        {
            int gt = bcf_gt_type(fmt, ismpl, NULL,NULL);
            if ( gt==GT_UNKN ) args->counts[ismpl*3+0]++;       // missing
            else if ( gt==GT_HET_RA || gt==GT_HET_AA ) args->counts[ismpl*3+1]++;   // het
            else args->counts[ismpl*3+2]++; // hom
        }
    }

    for (ismpl=0; ismpl<args->nsample; ismpl++)
    {
        float sum = 0, *probs = args->sex2prob + ismpl*args->nsex;
        int i, *counts = args->counts + ismpl*(args->max_ploidy+1);
        float fhet = (counts[1]+counts[2]) ? (float)counts[1]/(counts[1]+counts[2]) : 0;
        for (i=0; i<args->max_ploidy+1; i++) sum += counts[i];
        for (i=0; i<args->nsex; i++)
        {
            // a very simple heuristics to determine sex by counting hets/homs/missing sites,
            // in human nhet/nhom ~ 0.2
            int ploidy = args->sex2ploidy[i];
            float prob = 1;
            if ( ploidy==0 )
                prob = sum ? counts[0] / sum : 1;   // fraction of missing sites
            else if ( ploidy==1 )
            {
                if ( counts[1]+counts[2] ) prob = fhet > args->min_hets ? 0.1 : 0.9;
                prob *= sum ? 1 - counts[0] / sum : 1./args->nsex;
            }
            else 
            {
                if ( counts[1]+counts[2] ) prob = fhet > args->min_hets ? 0.9 : 0.1;
                prob *= sum ? 1 - counts[0] / sum : 1./args->nsex;
            }
            probs[i] *= prob;
        }
        if ( args->verbose )
            printf("DBG\t%s:%d-%d\t%s\t%f\t%d\t%d\t%d\n", seq,start+1,end+1,args->hdr->samples[ismpl], fhet,counts[0],counts[1],counts[2]);
    }

    return k;
}

int run(int argc, char **argv)
{
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->nsites = 10;
    args->min_hets = 0.05;
    static struct option loptions[] =
    {
        {"verbose",1,0,'v'},
        {"ploidy",1,0,'p'},
        {"nsites",1,0,'n'},
        {"guess",0,0,'g'},
        {"min-hets",1,0,'m'},
        {0,0,0,0}
    };
    char c, *tmp, *ploidy_fname = NULL;
    while ((c = getopt_long(argc, argv, "p:n:gm:v",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'v': args->verbose = 1; break; 
            case 'g': args->guess = 1; break;
            case 'm': 
                args->min_hets = strtod(optarg,&tmp); 
                if ( *tmp ) error("Unexpected argument to --min-hets: %s\n", optarg);
                break; 
            case 'p': ploidy_fname = optarg; break; 
            case 'n': 
                args->nsites = strtol(optarg,&tmp,10); 
                if (*tmp) error("Unexpected argument to --nsites: %s\n", optarg); break; 
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }

    args->sr = bcf_sr_init();
    args->sr->require_index = 1;
    if ( !argv[0] || !bcf_sr_add_reader(args->sr,argv[0]) ) error("%s", usage());
    args->hdr = args->sr->readers[0].header;
    args->nsample = bcf_hdr_nsamples(args->hdr);
 
    args->dflt_ploidy = 2;
    if ( ploidy_fname )
    {
        args->ploidy = ploidy_init(ploidy_fname, args->dflt_ploidy);
        if ( !args->ploidy ) error("Could not read %s\n", ploidy_fname);
    }
    else
    {
        args->ploidy = ploidy_init_string(
                "X 1 60000 M 1\n"
                "X 2699521 154931043 M 1\n"
                "Y 1 59373566 M 1\n"
                "Y 1 59373566 F 0\n", args->dflt_ploidy);
    }
    args->nsex = ploidy_nsex(args->ploidy);
    args->sex2ploidy = (int*) malloc(sizeof(int)*args->nsex);
    args->max_ploidy = ploidy_max(args->ploidy);
    if ( args->guess && args->max_ploidy > 2 ) error("Sorry, ploidy %d not supported with -g\n", args->max_ploidy);
    args->ncounts = args->nsample * ((args->max_ploidy>2 ? args->max_ploidy : 2)+1);
    args->counts = (int*) malloc(sizeof(int)*args->ncounts);
    args->sex2prob = (float*) calloc(args->nsample*args->nsex,sizeof(float));

    int i, nseq;
    for (i=0; i<args->nsample*args->nsex; i++) args->sex2prob[i] = 1;

    if ( args->verbose && args->guess )
        printf("# [1]DBG\t[2]Region\t[3]Sample\t[4]HET fraction\t[5]nHet\t[6]nHom\t[7]nMissing\n");

    regidx_t *idx = ploidy_regions(args->ploidy);
    char **seqs = regidx_seq_names(idx, &nseq);
    for (i=0; i<nseq; i++)
    {
        regitr_t itr;
        regidx_overlap(idx, seqs[i], 0, UINT32_MAX, &itr);
        while ( itr.i < itr.n )
            if ( args->guess )
                itr.i += process_region_guess(args, seqs[i], &itr);
            else
                itr.i += process_region_precise(args, seqs[i], &itr);
    }

    for (i=0; i<args->nsample; i++)
    {
        int j, jmax = 0;
        float max = 0, sum = 0;
        for (j=0; j<args->nsex; j++)
        {
            sum += args->sex2prob[i*args->nsex+j];
            if ( max < args->sex2prob[i*args->nsex+j] )
            {
                jmax = j;
                max = args->sex2prob[i*args->nsex+j];
            }
        }
        if ( args->verbose )
            printf("%s\t%s\t%f\n", args->hdr->samples[i],ploidy_id2sex(args->ploidy,jmax),args->sex2prob[i*args->nsex+jmax]/sum);
        else
            printf("%s\t%s\n", args->hdr->samples[i],ploidy_id2sex(args->ploidy,jmax));
    }
   
    bcf_sr_destroy(args->sr);
    ploidy_destroy(args->ploidy);
    free(args->sex2ploidy);
    free(args->counts);
    free(args->gts);
    free(args->sex2prob);
    free(args);
    return 0;
}




