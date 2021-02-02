/* The MIT License

   Copyright (c) 2018-2021 Genome Research Ltd.

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
#include <unistd.h>     // for isatty
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

static int NVAF = 20;
static int MAX_LEN = 20;

static inline int len2bin(int len)
{
    if ( len < -MAX_LEN ) return 0;
    if ( len > MAX_LEN ) return 2*MAX_LEN;
    return MAX_LEN + len;
}
HTS_UNUSED static inline int bin2len(int bin)
{
    return bin - MAX_LEN;
}
static inline int vaf2bin(float vaf)
{
    return vaf*(NVAF-1);
}
HTS_UNUSED static inline float bin2vaf(int bin)
{
    return (float)bin/(NVAF-1);
}

typedef struct
{
    uint32_t
        *nvaf,                  // number of indels genotypes with low VAF (<=0.2) and high VAF (>0.2); use vaf2bin and bin2vaf
        *nlen,                  // length distribution (-MAX_LEN,MAX_LEN); use len2bin and bin2len; site-wise unless samples are present
        npass_gt,               // number of indel genotypes passing the filter
        npass,                  // number of sites passing the filter
        nsites,                 // number of sites total
        nins, ndel,             // number of insertions and deletions, site-wise, not genotype-wise
        nframeshift, ninframe,  // site-wise
        *nfrac;                 // number of het indels contributing to dfrac
    double
        *dfrac;                 // minor allele fraction at HET indel genotypes, determined from FORMAT/AD
}
stats_t;

typedef struct
{
    stats_t stats;
    filter_t *filter;
    char *expr;
}
flt_stats_t;

#define iCHILD  0
#define iFATHER 1
#define iMOTHER 2

typedef struct
{
    int idx[3];     // VCF sample index for father, mother and child
    int pass;       // do all three pass the filters?
}
trio_t;

typedef struct
{
    int argc, filter_logic, regions_is_file, targets_is_file;
    int nflt_str;
    char *filter_str, **flt_str;
    char **argv, *output_fname, *fname, *regions, *targets, *csq_tag, *ped_fname;
    trio_t *trio;
    int ntrio, mtrio;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
    flt_stats_t *filters;
    int nfilters, nsmpl;
    char *csq_str;
    int32_t *gt_arr, *ad_arr, *ac;
    int mgt_arr, mad_arr, mac, mcsq_str;
    int ngt, ngt1, nad, nad1;
    int allow_alt2ref_DNMs;     // is "0/0 0/1 1/1" (child,father,mother) a valid DNM?
}
args_t;

args_t args;

const char *about(void)
{
    return "Calculate indel stats scanning over a range of thresholds simultaneously.\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Calculates indel stats. Use curly brackets to scan a range of values simultaneously\n"
        "Usage: bcftools +indel-stats [Plugin Options]\n"
        "Plugin options:\n"
        "       --alt2ref-DNM           consider GT errors such as 0/1 + 1/1 -> 0/0 a valid DNM\n"
        "   -c, --csq-tag STR           VEP or BCSQ tag to determine inframe and frameshift variants [CSQ]\n"
        "   -e, --exclude EXPR          exclude sites and samples for which the expression is true\n"
        "   -i, --include EXPR          include sites and samples for which the expression is true\n"
        "       --max-len INT           maximum indel length to consider [20]\n"
        "       --nvaf INT              number of variant allele frequency bins [20]\n"
        "   -o, --output FILE           output file name [stdout]\n"
        "   -p, --ped FILE              limit the stats to de novo indels\n"
        "   -r, --regions REG           restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE     restrict to regions listed in a file\n"
        "   -t, --targets REG           similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE     similar to -R but streams rather than index-jumps\n"
        "\n"
        "Example:\n"
        "   bcftools +indel-stats -i 'GQ>{10,20,30,40,50}' file.bcf\n"
        "\n";
}

static void parse_filters(args_t *args)
{
    if ( !args->filter_str ) return;
    int mflt = 1;
    args->nflt_str = 1;
    args->flt_str  = (char**) malloc(sizeof(char*));
    args->flt_str[0] = strdup(args->filter_str);
    while (1)
    {
        int i, expanded = 0;
        for (i=args->nflt_str-1; i>=0; i--)
        {
            char *exp_beg = strchr(args->flt_str[i], '{');
            if ( !exp_beg ) continue;
            char *exp_end = strchr(exp_beg+1, '}');
            if ( !exp_end ) error("Could not parse the expression: %s\n", args->filter_str);
            char *beg = exp_beg+1, *mid = beg;
            while ( mid<exp_end )
            {
                while ( mid<exp_end && *mid!=',' ) mid++;
                kstring_t tmp = {0,0,0};
                kputsn(args->flt_str[i], exp_beg - args->flt_str[i], &tmp);
                kputsn(beg, mid - beg, &tmp);
                kputs(exp_end+1, &tmp);
                args->nflt_str++;
                hts_expand(char*, args->nflt_str, mflt, args->flt_str);
                args->flt_str[args->nflt_str-1] = tmp.s;
                beg = ++mid;
            }
            expanded = 1;
            free(args->flt_str[i]);
            memmove(&args->flt_str[i], &args->flt_str[i+1], (args->nflt_str-i-1)*sizeof(*args->flt_str));
            args->nflt_str--;
            args->flt_str[args->nflt_str] = NULL;
        }
        if ( !expanded ) break;
    }
    
    fprintf(stderr,"Collecting data for %d filtering expressions\n", args->nflt_str);
}

static int cmp_trios(const void *_a, const void *_b)
{
    trio_t *a = (trio_t *) _a;
    trio_t *b = (trio_t *) _b;
    int i;
    int amin = a->idx[0];
    for (i=1; i<3; i++)
        if ( amin > a->idx[i] ) amin = a->idx[i];
    int bmin = b->idx[0];
    for (i=1; i<3; i++)
        if ( bmin > b->idx[i] ) bmin = b->idx[i];
    if ( amin < bmin ) return -1;
    if ( amin > bmin ) return 1;
    return 0;
}
static void parse_ped(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL;
    do
    {
        // familyID    sampleID paternalID maternalID sex   phenotype   population relationship   siblings   secondOrder   thirdOrder   children    comment
        // BB03    HG01884 HG01885 HG01956 2   0   ACB child   0   0   0   0
        int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<4 ) error("Could not parse the ped file: %s\n", str.s);

        int father = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[2]]);
        if ( father<0 ) continue;
        int mother = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[3]]);
        if ( mother<0 ) continue;
        int child = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( child<0 ) continue;

        args->ntrio++;
        hts_expand0(trio_t,args->ntrio,args->mtrio,args->trio);
        trio_t *trio = &args->trio[args->ntrio-1];
        trio->idx[iFATHER] = father;
        trio->idx[iMOTHER] = mother;
        trio->idx[iCHILD]  = child;
    }
    while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    fprintf(stderr,"Identified %d complete trios in the VCF file\n", args->ntrio);
    if ( !args->ntrio ) error("No complete trio identified\n");

    // sort the sample by index so that they are accessed more or less sequentially
    qsort(args->trio,args->ntrio,sizeof(trio_t),cmp_trios);
    
    free(str.s);
    free(off);
    if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,fname);
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->regions )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n",args->regions);
    }
    if ( args->targets && bcf_sr_set_targets(args->sr, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->targets);
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    if ( args->ped_fname )
        parse_ped(args, args->ped_fname);

    parse_filters(args);

    int i;
    if ( !args->nflt_str )
    {
        args->filters = (flt_stats_t*) calloc(1, sizeof(flt_stats_t));
        args->nfilters = 1;
        args->filters[0].expr = strdup("all");
        args->filters[0].stats.nvaf  = (uint32_t*) calloc(NVAF, sizeof(uint32_t));
        args->filters[0].stats.nlen  = (uint32_t*) calloc(MAX_LEN*2+1, sizeof(uint32_t));
        args->filters[0].stats.nfrac = (uint32_t*) calloc(MAX_LEN*2+1, sizeof(uint32_t));
        args->filters[0].stats.dfrac = (double*) calloc(MAX_LEN*2+1, sizeof(double));
    }
    else
    {
        args->nfilters = args->nflt_str;
        args->filters = (flt_stats_t*) calloc(args->nfilters, sizeof(flt_stats_t));
        for (i=0; i<args->nfilters; i++)
        {
            args->filters[i].filter = filter_init(args->hdr, args->flt_str[i]);
            args->filters[i].expr   = strdup(args->flt_str[i]);
            args->filters[i].stats.nvaf  = (uint32_t*) calloc(NVAF, sizeof(uint32_t));
            args->filters[i].stats.nlen  = (uint32_t*) calloc(MAX_LEN*2+1, sizeof(uint32_t));
            args->filters[i].stats.nfrac = (uint32_t*) calloc(MAX_LEN*2+1, sizeof(uint32_t));
            args->filters[i].stats.dfrac = (double*) calloc(MAX_LEN*2+1, sizeof(double));

            // replace tab's with spaces so that the output stays parsable
            char *tmp = args->filters[i].expr;
            while ( *tmp )
            { 
                if ( *tmp=='\t' ) *tmp = ' '; 
                tmp++; 
            }
        }
    }
    args->nsmpl = bcf_hdr_nsamples(args->hdr);
}
static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nfilters; i++)
    {
        if ( args->filters[i].filter ) filter_destroy(args->filters[i].filter);
        free(args->filters[i].stats.nvaf);
        free(args->filters[i].stats.nlen);
        free(args->filters[i].stats.nfrac);
        free(args->filters[i].stats.dfrac);
        free(args->filters[i].expr);
    }
    free(args->filters);
    for (i=0; i<args->nflt_str; i++) free(args->flt_str[i]);
    free(args->flt_str);
    bcf_sr_destroy(args->sr);
    free(args->ac);
    free(args->trio);
    free(args->csq_str);
    free(args->gt_arr);
    free(args->ad_arr);
    free(args);
}
static void report_stats(args_t *args)
{
    int i = 0,j;
    FILE *fh = !args->output_fname || !strcmp("-",args->output_fname) ? stdout : fopen(args->output_fname,"w");
    if ( !fh ) error("Could not open the file for writing: %s\n", args->output_fname);
    fprintf(fh,"# CMD line shows the command line used to generate this output\n");
    fprintf(fh,"# DEF lines define expressions for all tested thresholds\n");
    fprintf(fh,"# SN* summary number for every threshold:\n");
    fprintf(fh,"#   %d) SN*, filter id\n", ++i);
    fprintf(fh,"#   %d) number of samples (or trios with -p)\n", ++i);
    fprintf(fh,"#   %d) number of indel sites total\n", ++i);
    fprintf(fh,"#   %d) number of indel sites that pass the filter (and, with -p, have a de novo indel)\n", ++i);
    fprintf(fh,"#   %d) number of indel genotypes that pass the filter (and, with -p, are de novo)\n", ++i);
    fprintf(fh,"#   %d) number of insertions (site-wise, not genotype-wise)\n", ++i);
    fprintf(fh,"#   %d) number of deletions (site-wise, not genotype-wise)\n", ++i);
    fprintf(fh,"#   %d) number of frameshifts (site-wise, not genotype-wise)\n", ++i);
    fprintf(fh,"#   %d) number of inframe indels (site-wise, not genotype-wise)\n", ++i);
    fprintf(fh,"#\n");
    i = 0;
    fprintf(fh,"# DVAF* lines report indel variant allele frequency (VAF) distribution for every threshold,\n");
    fprintf(fh,"#   k-th bin corresponds to the frequency k/(nVAF-1):\n");
    fprintf(fh,"#   %d) DVAF*, filter id\n", ++i);
    fprintf(fh,"#   %d) nVAF, number of bins which split the [0,1] VAF interval.\n", ++i);
    fprintf(fh,"#   %d-%d) counts of indel genotypes in the VAF bin. For non-reference hets, the VAF of the less supported allele is recorded\n", i+1, i+NVAF);
    fprintf(fh,"#\n");
    i = 0;
    fprintf(fh,"# DLEN* lines report indel length distribution for every threshold. When genotype fields are available,\n");
    fprintf(fh,"#   the counts correspond to the number of genotypes, otherwise the number of sites are given.\n");
    fprintf(fh,"#   The k-th bin corresponds to the indel size k-MAX_LEN, negative for deletions, positive for insertions.\n");
    fprintf(fh,"#   The first/last bin contains also all deletions/insertions larger than MAX_LEN:\n");
    fprintf(fh,"#   %d) DLEN*, filter id\n", ++i);
    fprintf(fh,"#   %d) maximum indel length\n", ++i);
    fprintf(fh,"#   %d-%d) counts of indel lengths (-max,..,0,..,max), all unique alleles in a genotype are recorded (alt hets increase the counters 2x, alt homs 1x)\n", i+1, i+MAX_LEN*2+1);
    fprintf(fh,"#\n");
    i = 0;
    fprintf(fh,"# DFRAC* lines report the mean minor allele fraction at HET indel genotypes as a function of indel size.\n");
    fprintf(fh,"#   The format is the same as for DLEN:\n");
    fprintf(fh,"#   %d) DFRAC*, filter id\n", ++i);
    fprintf(fh,"#   %d) maximum indel length\n", ++i);
    fprintf(fh,"#   %d-%d) mean fraction at indel lengths (-max,..,0,..,max)\n", i+1, i+MAX_LEN*2+1);
    fprintf(fh,"#\n");
    i = 0;
    fprintf(fh,"# NFRAC* lines report the number of indels informing the DFRAC distribution.\n");
    fprintf(fh,"#   %d) NFRAC*, filter id\n", ++i);
    fprintf(fh,"#   %d) maximum indel length\n", ++i);
    fprintf(fh,"#   %d-%d) counts at indel lengths (-max,..,0,..,max)\n", i+1, i+MAX_LEN*2+1);
    fprintf(fh,"#\n");
    fprintf(fh, "CMD\t%s", args->argv[0]);
    for (i=1; i<args->argc; i++) fprintf(fh, " %s",args->argv[i]);
    fprintf(fh, "\n");
    for (i=0; i<args->nfilters; i++)
    {
        flt_stats_t *flt = &args->filters[i];
        fprintf(fh,"DEF\tFLT%d\t%s\n", i, flt->expr);
    }
    for (i=0; i<args->nfilters; i++)
    {
        stats_t *stats = &args->filters[i].stats;

        fprintf(fh,"SN%d", i);
        fprintf(fh,"\t%u", args->ntrio ? args->ntrio : args->nsmpl);
        fprintf(fh,"\t%u", stats->nsites);
        fprintf(fh,"\t%u", stats->npass);
        fprintf(fh,"\t%u", stats->npass_gt);
        fprintf(fh,"\t%u", stats->nins);
        fprintf(fh,"\t%u", stats->ndel);
        fprintf(fh,"\t%u", stats->nframeshift);
        fprintf(fh,"\t%u", stats->ninframe);
        fprintf(fh,"\n");

        fprintf(fh,"DVAF%d", i);
        fprintf(fh,"\t%d", NVAF);
        for (j=0; j<NVAF; j++) fprintf(fh,"\t%u",stats->nvaf[j]);
        fprintf(fh,"\n");

        fprintf(fh,"DLEN%d", i);
        fprintf(fh,"\t%d", MAX_LEN);
        for (j=0; j<MAX_LEN*2+1; j++) fprintf(fh,"\t%u",stats->nlen[j]);
        fprintf(fh,"\n");

        fprintf(fh,"DFRAC%d", i);
        fprintf(fh,"\t%d", MAX_LEN);
        for (j=0; j<MAX_LEN*2+1; j++)
            if ( stats->nfrac[j] ) fprintf(fh,"\t%.2f",stats->dfrac[j]/stats->nfrac[j]);
            else fprintf(fh,"\t.");
        fprintf(fh,"\n");

        fprintf(fh,"NFRAC%d", i);
        fprintf(fh,"\t%d", MAX_LEN);
        for (j=0; j<MAX_LEN*2+1; j++) fprintf(fh,"\t%u",stats->nfrac[j]);
        fprintf(fh,"\n");
    }
    if ( fclose(fh)!=0 ) error("Close failed: %s\n", (!args->output_fname || !strcmp("-",args->output_fname)) ? "stdout" : args->output_fname);
}

static inline int parse_genotype(int32_t *arr, int ngt1, int idx, int als[2])
{
    int32_t *ptr = arr + ngt1 * idx;
    if ( bcf_gt_is_missing(ptr[0]) ) return -1;
    als[0] = bcf_gt_allele(ptr[0]);

    if ( ngt1==1 || ptr[1]==bcf_int32_vector_end ) { ptr[1] = ptr[0]; return -2; }

    if ( bcf_gt_is_missing(ptr[1]) ) return -1;
    als[1] = bcf_gt_allele(ptr[1]);

    return 0;
}

static inline void update_indel_stats(args_t *args, bcf1_t *rec, stats_t *stats, int ismpl, int *als)
{
    int j;
    if ( als[0] >= args->nad1 || als[1] >= args->nad1 ) error("Incorrect GT allele at %s:%"PRId64" .. %d/%d\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,als[0],als[1]);
    int32_t *ad_ptr = args->ad_arr + ismpl*args->nad1;

    // find the allele with most support
    uint32_t ntot = 0;
    for (j=0; j<args->nad1; j++)
    {
        if ( ad_ptr[j]==bcf_int32_missing ) continue;
        if ( ad_ptr[j]==bcf_int32_vector_end ) break;
        ntot += ad_ptr[j];
    }
    if ( !ntot ) return;

    // Find the alternate allele fraction, total and relative. Set al0 to be the more frequent indel allele.
    // The genotypes have been already sanitized in parse_genotype().
    int al0 = als[0], al1 = als[1];
    if ( !(bcf_get_variant_type(rec,al0) & VCF_INDEL) )
    {
        if ( !(bcf_get_variant_type(rec,al1) & VCF_INDEL) ) error("FIXME: this should not happen .. %s:%"PRId64" .. %d/%d\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,al0,al1);
        al0 = als[1]; al1 = als[0];
    }
    else if ( (bcf_get_variant_type(rec,al1) & VCF_INDEL) && al0!=al1 )
    {
        // Select the more frequent indel allele.
        if ( ad_ptr[al0] < ad_ptr[al1] ) al0 = als[1], al1 = als[0];

        // Record length of both indel alleles
        int bin = len2bin(rec->d.var[al1].n);
        if ( bin >= 0 ) stats->nlen[bin]++;
    }

    float vaf = (float)ad_ptr[al0] / ntot;
    int bin = vaf2bin(vaf);
    stats->nvaf[bin]++;

    // al0 is now the major indel allele
    int len_bin = len2bin(rec->d.var[al0].n);
    if ( len_bin < 0 ) return;
    stats->nlen[len_bin]++;

    if ( al0!=al1 )
    {
        ntot = ad_ptr[al0] + ad_ptr[al1];
        if ( ntot )
        {
            stats->nfrac[len_bin]++;
            stats->dfrac[len_bin]+= (double)ad_ptr[al0] / ntot;
        }
    }
}

static void process_record(args_t *args, bcf1_t *rec, flt_stats_t *flt)
{
    int i,j;
    uint8_t *smpl_pass = NULL;

    stats_t *stats = &flt->stats;
    stats->nsites++;

    // Find out which samples/trios pass and if the site passes
    if ( flt->filter )
    {
        int pass_site = filter_test(flt->filter, rec, (const uint8_t**) &smpl_pass);
        if ( args->ntrio )
        {
            if ( args->filter_logic & FLT_EXCLUDE )
            {
                if ( pass_site )
                {
                    if ( !smpl_pass ) return;
                    pass_site = 0;
                    for (i=0; i<args->ntrio; i++)
                    {
                        int pass_trio = 1;
                        for (j=0; j<3; j++)
                        {
                            int idx = args->trio[i].idx[j];
                            if ( smpl_pass[idx] ) { pass_trio = 0; break; }
                        }
                        args->trio[i].pass = pass_trio;
                        if ( pass_trio ) pass_site = 1;
                    }
                    if ( !pass_site ) return;
                }
                else
                    for (i=0; i<args->ntrio; i++) args->trio[i].pass = 1;
            }
            else if ( !pass_site ) return;
            else if ( smpl_pass )
            {
                pass_site = 0;
                for (i=0; i<args->ntrio; i++)
                {
                    int pass_trio = 1;
                    for (j=0; j<3; j++)
                    {
                        int idx = args->trio[i].idx[j];
                        if ( !smpl_pass[idx] ) { pass_trio = 0; break; }
                    }
                    args->trio[i].pass = pass_trio;
                    if ( pass_trio ) pass_site = 1;
                }
                if ( !pass_site ) return;
            }
            else
                for (i=0; i<args->ntrio; i++) args->trio[i].pass = 1;
        }
        else
        {
            if ( args->filter_logic & FLT_EXCLUDE )
            {
                if ( pass_site )
                {
                    if ( !smpl_pass ) return;
                    pass_site = 0;
                    for (i=0; i<args->nsmpl; i++)
                    {
                        if ( smpl_pass[i] ) smpl_pass[i] = 0;
                        else { smpl_pass[i] = 1; pass_site = 1; }
                    }
                    if ( !pass_site ) return;
                }
                else
                    for (i=0; i<args->nsmpl; i++) smpl_pass[i] = 1;
            }
            else if ( !pass_site ) return;
        }
    }

    args->ngt = 0;
    if ( args->nsmpl )
    {
        // Get the genotypes
        args->ngt = bcf_get_genotypes(args->hdr, rec, &args->gt_arr, &args->mgt_arr);
        args->ngt1 = args->ngt / rec->n_sample;

        if ( args->ngt>0 )
        {
            // Get the AD counts
            args->nad = bcf_get_format_int32(args->hdr, rec, "AD", &args->ad_arr, &args->mad_arr);
            args->nad1 = args->nad / rec->n_sample;
            if ( args->nad>0 && args->nad1 != rec->n_allele ) error("Incorrect number of FORMAT/AD values at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        }
    }

    // Is there a star allele? Don't count overlapping deletions twice
    int star_allele = -1;
    for (i=1; i<rec->n_allele; i++)
        if ( !rec->d.allele[i][1] && rec->d.allele[i][0]=='*' ) { star_allele = i; break; }


    if ( args->ngt>0 && args->ntrio )
    {
        int is_dnm = 0;
        for (i=0; i<args->ntrio; i++)
        {
            if ( flt->filter && !args->trio[i].pass ) continue;

            // Determine the alternate allele and the genotypes, skip if any of the alleles is missing.
            // the order is: child, father, mother
            int als[6], *als_child = als, *als_father = als+2, *als_mother = als+4; 
            if ( parse_genotype(args->gt_arr, args->ngt1, args->trio[i].idx[iCHILD], als_child) < 0 ) continue;
            if ( parse_genotype(args->gt_arr, args->ngt1, args->trio[i].idx[iFATHER], als_father) < 0 ) continue;
            if ( parse_genotype(args->gt_arr, args->ngt1, args->trio[i].idx[iMOTHER], als_mother) < 0 ) continue;

            // Is it a DNM?
            if ( !args->allow_alt2ref_DNMs && als_child[0]==0 && als_child[1]==0 ) continue;
            if ( (als_child[0]==als_father[0] || als_child[0]==als_father[1]) && (als_child[1]==als_mother[0] || als_child[1]==als_mother[1]) ) continue;
            if ( (als_child[1]==als_father[0] || als_child[1]==als_father[1]) && (als_child[0]==als_mother[0] || als_child[0]==als_mother[1]) ) continue;
            if ( als_child[0]==star_allele || als_child[1]==star_allele ) continue;     // don't count the same event multiple times
            if ( als_father[0]==star_allele || als_father[1]==star_allele ) continue;
            if ( als_mother[0]==star_allele || als_mother[1]==star_allele ) continue;

            int child_is_indel = (bcf_get_variant_type(rec,als_child[0]) & VCF_INDEL) || (bcf_get_variant_type(rec,als_child[1]) & VCF_INDEL) ? 1 : 0;

            if ( !args->allow_alt2ref_DNMs )
            {
                if ( !child_is_indel ) continue;
            }
            else 
            {
                if ( !child_is_indel &&
                     !(bcf_get_variant_type(rec,als_father[0]) & VCF_INDEL) &&
                     !(bcf_get_variant_type(rec,als_father[1]) & VCF_INDEL) &&
                     !(bcf_get_variant_type(rec,als_mother[0]) & VCF_INDEL) &&
                     !(bcf_get_variant_type(rec,als_mother[1]) & VCF_INDEL) ) continue; // not an indel, in any sample
            }

            if ( child_is_indel )
                update_indel_stats(args, rec, stats, args->trio[i].idx[iCHILD], als_child);

            //printf("MERR\t%s\t%d\t%s\n", bcf_seqname(args->hdr,rec),rec->pos+1,args->hdr->samples[args->trio[i].idx[iCHILD]]);

            stats->npass_gt++;

            is_dnm = 1;
        }
        if ( !is_dnm ) return;
    }
    else if ( args->ngt>0 )
    {
        for (i=0; i<args->nsmpl; i++)
        {
            if ( smpl_pass && !smpl_pass[i] ) continue;

            // Determine the alternate allele and the genotypes, skip if any of the alleles is missing.
            int als[2] = {0,0};
            int ret = parse_genotype(args->gt_arr, args->ngt1, i, als);
            if ( ret==-1 ) continue;    // missing genotype
            if ( !(bcf_get_variant_type(rec,als[0]) & VCF_INDEL) && !(bcf_get_variant_type(rec,als[1]) & VCF_INDEL) ) continue;     // not an indel

            update_indel_stats(args, rec, stats, i, als);

            stats->npass_gt++;
        }
    }

    if ( bcf_get_info_string(args->hdr,rec,args->csq_tag,&args->csq_str,&args->mcsq_str) > 0 )
    {
        if ( strstr(args->csq_str,"inframe") ) stats->ninframe++;
        if ( strstr(args->csq_str,"frameshift") ) stats->nframeshift++;
    }

    for (i=1; i<rec->n_allele; i++)
    {
        if ( !(bcf_get_variant_type(rec,i) & VCF_INDEL) ) continue;
        if ( rec->d.var[i].n < 0 ) stats->ndel++;
        else if ( rec->d.var[i].n > 0 ) stats->nins++;
        if ( args->ngt <= 0 )
        {
            int bin = len2bin(rec->d.var[i].n);
            if ( bin >= 0 ) stats->nlen[bin]++;
        }
    }
    stats->npass++;
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->csq_tag = "CSQ";
    static struct option loptions[] =
    {
        {"max-len",required_argument,0,1},
        {"nvaf",required_argument,0,2},
        {"alt2ref-DNM",no_argument,0,3},
        {"ped",required_argument,0,'p'},
        {"csq-tag",required_argument,0,'c'},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"output",required_argument,NULL,'o'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {NULL,0,NULL,0}
    };
    char *tmp;
    int c, i;
    while ((c = getopt_long(argc, argv, "o:s:i:e:r:R:t:T:c:p:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case  1 :
                MAX_LEN = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: --max-len %s\n", optarg);
                if ( MAX_LEN<=0 ) error("Expected value bigger than 0 --max-len\n");
                break;
            case  2 :
                NVAF = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: --max-len %s\n", optarg);
                if ( NVAF<0 || NVAF>1 ) error("Expected value from the interval [0,1] with --nvaf\n");
                break;
            case  3 : args->allow_alt2ref_DNMs = 1; break;
            case 'p': args->ped_fname = optarg; break;
            case 'c': args->csq_tag = optarg; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; args->targets_is_file = 1; break;
            case 'r': args->regions = optarg; break;
            case 'R': args->regions = optarg; args->regions_is_file = 1; break;
            case 'o': args->output_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s",usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s",usage_text());
    else args->fname = argv[optind];

    init_data(args);

    while ( bcf_sr_next_line(args->sr) )
    {
        bcf1_t *rec = bcf_sr_get_line(args->sr,0);
        if ( !(bcf_get_variant_types(rec) & VCF_INDEL) ) continue;
        for (i=0; i<args->nfilters; i++)
            process_record(args, rec, &args->filters[i]);
    }

    report_stats(args);
    destroy_data(args);

    return 0;
}
