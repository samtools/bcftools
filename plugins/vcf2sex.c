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
#include <inttypes.h>
#include "bcftools.h"
#include "ploidy.h"

#define GUESS_GT 1
#define GUESS_PL 2
#define GUESS_GL 3

typedef struct
{
    uint64_t nhet, nhom, nmiss;
}
count_t;

typedef struct
{
    char *chr;
    uint32_t start, end;
    int *sex2ploidy;    // sex ploidies
    count_t *counts;    // per-sample counts: counts[isample]
}
reg_stats_t;

typedef struct
{
    int nsites, nsex, *sex2ploidy, dflt_ploidy, max_ploidy, guess;
    count_t *bg_counts;         // background het/hom counts for regions with the same ploidy in all sexes
    reg_stats_t *reg_stats;     // counts for all regions, used with -g
    int nreg_stats, nsample, verbose;
    int *counts, ncounts;       // number of observed GTs with given ploidy, used when -g is not given
    float *sex2prob, min_hets;
    int32_t *gts, ngts, *pls, npls;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
    ploidy_t *ploidy;
    char *background;
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
        "About: Determine sample sex by checking the presence of haploid/diploid genotypes\n"
        "       (requires correct ploidy in the VCF). Alternatively, determine the sex by\n"
        "       counting homs/hets (-g GT) or most likely hom/het genotypes (-g PL) in haploid\n"
        "       regions. With -g, a region is deemed haploid when the fraction of hets in the\n"
        "       region is significantly smaller than in the background region. For example,\n"
        "       the default \"-m 0.3\" requires the het rate to be at least 30% of the background\n"
        "       rate. If the background region is disabled (\"-b -\"), the threshold -m is interpreted\n"
        "       directly as the minimum fraction of hets in the region.\n"
        "\n"
        "Usage: bcftools +vcf2sex <file.vcf.gz> -- [Plugin Options]\n"
        "Plugin options:\n"
        "   -b, --background <region>   diploid region to determine normal hom/hets counts [X:60001-2699520]\n"
        "   -g, --guess <tag>           determine ploidy by counting hom/hets (GT) or most likely genotypes (PL or GL)\n"
        "   -m, --min-hets <float>      minimum fraction of hets in diploid regions [0.3]\n"
        "   -n, --nsites <int>          number of sites to check per region (ignored with -g) [10]\n"
        "   -p, --ploidy <file>         space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY\n"
        "   -v, --verbose               print debugging information\n"
        "\n"
        "Example:\n"
        "   # Default ploidy, if -p not given. Unlisted regions have ploidy 2\n"
        "   X 1 60000 M 1\n"
        "   X 2699521 154931043 M 1\n"
        "   Y 1 59373566 M 1\n"
        "   Y 1 59373566 F 0\n"
        "   \n"
        "   bcftools +vcf2sex in.vcf.gz\n"
        "   bcftools +vcf2sex -n 10 in.vcf.gz\n"
        "   bcftools +vcf2sex -g GT in.vcf.gz\n"
        "\n";
}

reg_stats_t *expand_regs(args_t *args, char *seq, uint32_t start, uint32_t end)
{
    args->nreg_stats++;
    args->reg_stats = (reg_stats_t*) realloc(args->reg_stats,sizeof(reg_stats_t)*args->nreg_stats);
    reg_stats_t *stats = args->reg_stats + args->nreg_stats-1;
    stats->chr = strdup(seq);
    stats->start = start;
    stats->end = end;
    stats->sex2ploidy = (int*) malloc(sizeof(int)*args->nsex);
    memcpy(stats->sex2ploidy,args->sex2ploidy,sizeof(int)*args->nsex);
    stats->counts = (count_t*) calloc(args->nsample,sizeof(count_t));
    return stats;
}
void destroy_regs(args_t *args)
{   
    int i;
    for (i=0; i<args->nreg_stats; i++)
    {
        free(args->reg_stats[i].chr);
        free(args->reg_stats[i].counts);
        free(args->reg_stats[i].sex2ploidy);
    }
    free(args->reg_stats);
    args->nreg_stats = 0;
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
            {
                if ( gts[igt]==bcf_int32_vector_end || bcf_gt_is_missing(gts[igt]) ) break;
                else ploidy++;
            }
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
    int kitr = 1;
    uint32_t start = 0, end = INT_MAX;
    reg_stats_t *stats = NULL;

    // set the start and the end position
    if ( itr )
    {
        start = itr->reg[itr->i].start;
        end   = itr->reg[itr->i].end;

        // flush all records with the same coordinates
        while ( itr->i+kitr<itr->n && start==itr->reg[itr->i+kitr].start && end==itr->reg[itr->i+kitr].end ) kitr++;

        int min,max,ret = ploidy_query(args->ploidy, seq, start, args->sex2ploidy, &min, &max);
        assert(ret);
        stats = expand_regs(args, seq,start,end);
    }
    else
    {
        // background region
        int spos, epos;
        const char *ptr = hts_parse_reg(args->background, &spos, &epos);
        if ( !ptr )
            error("Could not parse the region: %s\n", args->background);
        seq = (char*) malloc(ptr - args->background + 1);
        memcpy(seq,args->background,ptr-args->background);
        seq[ptr-args->background] = 0;
        start = spos;
        end   = epos;
    }

    if ( bcf_sr_seek(args->sr,seq,start)!=0 ) 
    {
        // sequence not present
        if ( !itr ) free(seq);
        return kitr;
    }

    int ismpl, rid = bcf_hdr_name2id(args->hdr,seq);
    if ( !itr ) free(seq);

    while ( bcf_sr_next_line(args->sr) )
    {
        bcf1_t *rec = bcf_sr_get_line(args->sr,0);
        if ( rec->rid!=rid || rec->pos > end ) break;

        if ( args->guess & GUESS_GT )   // use GTs to guess the ploidy
        {
            bcf_fmt_t *fmt = bcf_get_fmt(args->hdr, rec, "GT");
            if ( !fmt ) continue;
            for (ismpl=0; ismpl<args->nsample; ismpl++)
            {
                count_t *counts = stats ? &stats->counts[ismpl] : &args->bg_counts[ismpl];
                int gt = bcf_gt_type(fmt, ismpl, NULL,NULL);
                if ( gt==GT_UNKN ) counts->nmiss++;
                else if ( gt==GT_HET_RA || gt==GT_HET_AA ) counts->nhet++;
                else counts->nhom++;
            }
        }
        else    // use PLs to guess the ploidy
        {
            int gl2pl = args->guess & GUESS_PL ? 1 : -1;
            int npl = bcf_get_format_int32(args->hdr,rec,args->guess&GUESS_PL?"PL":"GL",&args->pls,&args->npls);
            if ( npl<=0 ) continue;
            npl /= args->nsample;
            for (ismpl=0; ismpl<args->nsample; ismpl++)
            {
                int32_t *ptr = args->pls + ismpl*npl;
                int phom = INT_MAX, phet = INT_MAX, ial, jal, k = 0;
                for (ial=0; ial<rec->n_allele; ial++)
                {
                    for (jal=0; jal<ial; jal++)
                    {
                        if ( ptr[k] == bcf_int32_missing || ptr[k] == bcf_int32_vector_end )  break;
                        ptr[k] *= gl2pl;
                        if ( phet > ptr[k] ) phet = ptr[k];
                        k++;
                    }
                    if ( ptr[k] == bcf_int32_missing || ptr[k] == bcf_int32_vector_end )  break;
                    ptr[k] *= gl2pl;
                    if ( phom > ptr[k] ) phom = ptr[k];
                    k++;
                }
                count_t *counts = stats ? &stats->counts[ismpl] : &args->bg_counts[ismpl];
                if ( k == rec->n_allele ) counts->nhom++;   // haploid
                else if ( phet == phom || k != rec->n_allele*(rec->n_allele+1)/2 ) counts->nmiss++;
                else if ( phet < phom ) counts->nhet++;
                else counts->nhom++;
            }
        }
    }
    return kitr;
}

void sex2prob_guess(args_t *args)
{
    int ismpl, ireg; 

    // get numbers from the background region
    if ( args->background )
    {
        process_region_guess(args, NULL, NULL);

        if ( args->verbose )
            printf("# [1]BGR\t[2]Region\t[3]Sample\t[4]Het fraction\t[5]nHet\t[6]nHom\t[7]nMissing\n");

        for (ismpl=0; ismpl<args->nsample; ismpl++)
        {
            uint64_t nhom  = args->bg_counts[ismpl].nhom;
            uint64_t nhet  = args->bg_counts[ismpl].nhet;
            uint64_t nmiss = args->bg_counts[ismpl].nmiss;
            if ( nhom+nhet==0 )
                fprintf(stderr,"Warning: The sample %s has no variants in the background region %s\n", args->hdr->samples[ismpl],args->background);

            float bg_het = (float)nhet/(nhom+nhet);
            if ( args->verbose )
                printf("BGR\t%s\t%s\t%f\t%"PRId64"\t%"PRId64"\t%"PRId64"\n", args->background,args->hdr->samples[ismpl],bg_het,nhet,nhom,nmiss);
        }
    }

    // a very simple heuristics to determine sex by counting hets/homs/missing sites
    for (ireg=0; ireg<args->nreg_stats; ireg++)
    {
        reg_stats_t *stats = &args->reg_stats[ireg];
        for (ismpl=0; ismpl<args->nsample; ismpl++)
        {
            uint64_t nhom  = stats->counts[ismpl].nhom;
            uint64_t nhet  = stats->counts[ismpl].nhet;
            uint64_t nmiss = stats->counts[ismpl].nmiss;
            float fhet   = nhom+nhet ? (float)nhet/(nhom+nhet) : 0;

            float bg_fhet = -1;
            if ( args->background )
            {
                uint64_t bg_nhom  = args->bg_counts[ismpl].nhom;
                uint64_t bg_nhet  = args->bg_counts[ismpl].nhet;
                bg_fhet = bg_nhom+bg_nhet ? (float)bg_nhet/(bg_nhom+bg_nhet) : 0;
            }

            if ( args->verbose )
                printf("REG\t%s:%d-%d\t%s\t%f\t%"PRId64"\t%"PRId64"\t%"PRId64"\n", stats->chr,stats->start+1,stats->end+1,args->hdr->samples[ismpl],fhet,nhet,nhom,nmiss);

            int i, ntot = nhom + nhet + nmiss;
            if ( !ntot ) continue;

            float *probs = args->sex2prob + ismpl*args->nsex;
            for (i=0; i<args->nsex; i++)
            {
                int ploidy = stats->sex2ploidy[i];
                float prob;
                if ( ploidy==0 )
                    prob = (float) nmiss / ntot;   // fraction of missing sites
                else if ( ploidy==1 )
                {
                    // NB: these numbers (0.1,0.9) are made up, to be improved
                    if ( bg_fhet<0 )
                        prob = fhet > args->min_hets ? 0.1 : 0.9;
                    else
                        prob = fhet > args->min_hets*bg_fhet ? 0.1 : 0.9;
                    prob *= 1 - (float)nmiss / ntot;
                }
                else 
                {
                    if ( bg_fhet<0 )
                        prob = fhet > args->min_hets ? 0.9 : 0.1;
                    else
                        prob = fhet > args->min_hets*bg_fhet ? 0.9 : 0.1;
                    prob *= 1 - (float)nmiss / ntot;
                }
                probs[i] *= prob;
            }
        }
    }
}

int run(int argc, char **argv)
{
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->nsites = 10;
    args->min_hets = 0.3;
    args->background = "X:60001-2699520";
    static struct option loptions[] =
    {
        {"verbose",1,0,'v'},
        {"ploidy",1,0,'p'},
        {"nsites",1,0,'n'},
        {"guess",1,0,'g'},
        {"min-hets",1,0,'m'},
        {"background",1,0,'b'},
        {0,0,0,0}
    };
    char *tmp, *ploidy_fname = NULL;
    int c;
    while ((c = getopt_long(argc, argv, "p:n:g:m:vb:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'b': 
                if ( !strcmp("-",optarg) ) args->background = NULL;
                else args->background = optarg; 
                break; 
            case 'v': args->verbose = 1; break; 
            case 'g':
                if ( !strcasecmp(optarg,"GT") ) args->guess = GUESS_GT;
                else if ( !strcasecmp(optarg,"PL") ) args->guess = GUESS_PL;
                else if ( !strcasecmp(optarg,"GL") ) args->guess = GUESS_GL;
                else error("The argument not recognised, expected --guess GT, --guess PL or --guess GL: %s\n", optarg);
                break;
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
    if ( optind==argc ) error("%s", usage());
    if ( !bcf_sr_add_reader(args->sr,argv[optind]) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
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
    args->bg_counts = (count_t*) calloc(args->nsample,sizeof(count_t));
    args->sex2prob = (float*) calloc(args->nsample*args->nsex,sizeof(float));

    int i, nseq;
    for (i=0; i<args->nsample*args->nsex; i++) args->sex2prob[i] = 1;

    if ( args->verbose && args->guess )
        printf("# [1]REG\t[2]Region\t[3]Sample\t[4]Het fraction\t[5]nHet\t[6]nHom\t[7]nMissing\n");

    // First get the counts from expected haploid regions
    regidx_t *idx = ploidy_regions(args->ploidy);
    char **seqs = regidx_seq_names(idx, &nseq);
    for (i=0; i<nseq; i++)
    {
        regitr_t itr;
        regidx_overlap(idx, seqs[i], 0, UINT32_MAX, &itr);
        while ( itr.i < itr.n )
        {
            if ( args->guess )
                itr.i += process_region_guess(args, seqs[i], &itr);
            else
                itr.i += process_region_precise(args, seqs[i], &itr);
        }
    }
    // Get the counts from a PAR (the background diploid region) and see if the fraction
    // of hets is different
    if ( args->guess ) sex2prob_guess(args);

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
    destroy_regs(args);
    free(args->sex2ploidy);
    free(args->counts);
    free(args->bg_counts);
    free(args->gts);
    free(args->pls);
    free(args->sex2prob);
    free(args);
    return 0;
}




