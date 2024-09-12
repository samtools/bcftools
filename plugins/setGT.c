/*  plugins/setGT.c -- set gentoypes to given values

    Copyright (C) 2015-2024 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kfunc.h>
#include <htslib/hts_os.h>
#include <inttypes.h>
#include <getopt.h>
#include <ctype.h>
#include "bcftools.h"
#include "filter.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef int (*cmp_f)(double a, double b);

static int cmp_eq(double a, double b) { return a==b ? 1 : 0; }
static int cmp_le(double a, double b) { return a<=b ? 1 : 0; }
static int cmp_ge(double a, double b) { return a>=b ? 1 : 0; }
static int cmp_lt(double a, double b) { return a<b ? 1 : 0; }
static int cmp_gt(double a, double b) { return a>b ? 1 : 0; }

typedef struct
{
    bcf_hdr_t *in_hdr, *out_hdr;
    int32_t *gts, mgts, *iarr, miarr, *xarr, mxarr;
    int *arr, marr;
    uint64_t nchanged;
    int tgt_mask, new_mask, new_gt;
    filter_t *filter;
    char *filter_str;
    struct {
        int m_allele, M_allele, x_vaf_allele, *gt, *phased, ploidy;
        char *gt_str;
    } custom;
    int filter_logic, rand_seed;
    uint8_t *smpl_pass;
    double binom_val, rand_frac;
    char *binom_tag;
    cmp_f binom_cmp;
}
args_t;

args_t *args = NULL;

#define GT_MISSING   1
#define GT_PARTIAL  (1<<1)
#define GT_REF      (1<<2)
#define GT_MAJOR    (1<<3)
#define GT_PHASED   (1<<4)
#define GT_UNPHASED (1<<5)
#define GT_ALL      (1<<6)
#define GT_QUERY    (1<<7)
#define GT_BINOM    (1<<8)
#define GT_MINOR    (1<<9)
#define GT_CUSTOM   (1<<10)
#define GT_X_VAF    (1<<11)
#define GT_RAND     (1<<12)

#define MINOR_ALLELE -1
#define MAJOR_ALLELE -2
#define X_VAF_ALLELE -3

const char *about(void)
{
    return "Set genotypes: partially missing to missing, missing to ref/major allele, etc.\n";
}

const char *usage(void)
{
    return
        "About: Sets genotypes. The target genotypes can be specified as:\n"
        "           ./.     .. completely missing (\".\" or \"./.\", depending on ploidy)\n"
        "           ./x     .. partially missing (e.g., \"./0\" or \".|1\" but not \"./.\")\n"
        "           .       .. partially or completely missing\n"
        "           a       .. all genotypes\n"
        "           b       .. heterozygous genotypes failing two-tailed binomial test (example below)\n"
        "           q       .. select genotypes using -i/-e options\n"
        "           r:FLOAT .. select randomly a proportion of FLOAT genotypes (can be combined with other modes)\n"
        "       and the new genotype can be one of:\n"
        "           .       .. missing (\".\" or \"./.\", keeps ploidy)\n"
        "           0       .. reference allele (e.g. 0/0 or 0, keeps ploidy)\n"
        "           c:GT    .. custom genotype (e.g. 0/0, 0, 0/1, m/M, 0/X overrides ploidy)\n"
        "           m       .. minor (the second most common) allele as determined from INFO/AC or FMT/GT (e.g. 1/1 or 1, keeps ploidy)\n"
        "           M       .. major allele as determined from INFO/AC or FMT/GT (e.g. 1/1 or 1, keeps ploidy)\n"
        "           X       .. allele with bigger read depth as determined from FMT/AD\n"
        "           p       .. phase genotype (0/1 becomes 0|1)\n"
        "           u       .. unphase genotype and sort by allele (1|0 becomes 0/1)\n"
        "Usage: bcftools +setGT [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -e, --exclude EXPR        Exclude a genotype if true (requires -t q)\n"
        "   -i, --include EXPR        include a genotype if true (requires -t q)\n"
        "   -n, --new-gt TYPE         Genotypes to set, see above\n"
        "   -s, --seed INT            Random seed to use with -t r [0]\n"
        "   -t, --target-gt TYPE      Genotypes to change, see above\n"
        "\n"
        "Example:\n"
        "   # set missing genotypes (\"./.\") to phased ref genotypes (\"0|0\")\n"
        "   bcftools +setGT in.vcf -- -t . -n 0p\n"
        "\n"
        "   # set missing genotypes with DP>0 and GQ>20 to ref genotypes (\"0/0\")\n"
        "   bcftools +setGT in.vcf -- -t q -n 0 -i 'GT=\".\" && FMT/DP>0 && GQ>20'\n"
        "\n"
        "   # set partially missing genotypes to completely missing\n"
        "   bcftools +setGT in.vcf -- -t ./x -n .\n"
        "\n"
        "   # set heterozygous genotypes to 0/0 if binom.test(nAlt,nRef+nAlt,0.5)<1e-3\n"
        "   bcftools +setGT in.vcf -- -t \"b:AD<1e-3\" -n 0\n"
        "\n"
        "   # force unphased heterozygous genotype if binom.test(nAlt,nRef+nAlt,0.5)>0.1\n"
        "   bcftools +setGT in.vcf -- -t ./x -n c:'m/M'\n"
        "\n";
}

static void _parse_binom_expr_error(char *str)
{
    error(
            "Error parsing the expression: %s\n"
            "Expected TAG CMP VAL, where\n"
            "   TAG .. one of the format tags\n"
            "   CMP .. operator, one of <, <=, >, >=\n"
            "   VAL .. value\n"
            "For example:\n"
            "   bcftools +setGT in.vcf -- -t \"b:AD>1e-3\" -n 0\n"
            "\n", str
         );
}
void parse_binom_expr(args_t *args, char *str)
{
    if ( str[1]!=':' ) _parse_binom_expr_error(str);

    char *beg = str+2;
    while ( *beg && isspace(*beg) ) beg++;
    if ( !*beg ) _parse_binom_expr_error(str);
    char *end = beg;
    while ( *end )
    {
        if ( isspace(*end) || *end=='<' || *end=='=' || *end=='>' ) break;
        end++;
    }
    if ( !*end ) _parse_binom_expr_error(str);
    args->binom_tag = (char*) calloc(1,end-beg+1);
    memcpy(args->binom_tag,beg,end-beg);
    int tag_id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,args->binom_tag);
    if ( !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The FORMAT tag \"%s\" is not present in the VCF\n", args->binom_tag);

    while ( *end && isspace(*end) ) end++;
    if ( !*end ) _parse_binom_expr_error(str);

    if ( !strncmp(end,"<=",2) ) { args->binom_cmp = cmp_le; beg = end+2; }
    else if ( !strncmp(end,">=",2) ) { args->binom_cmp = cmp_ge; beg = end+2; }
    else if ( !strncmp(end,"==",2) ) { args->binom_cmp = cmp_eq; beg = end+2; }
    else if ( !strncmp(end,"<",1) ) { args->binom_cmp = cmp_lt; beg = end+1; }
    else if ( !strncmp(end,">",1) ) { args->binom_cmp = cmp_gt; beg = end+1; }
    else if ( !strncmp(end,"=",1) ) { args->binom_cmp = cmp_eq; beg = end+1; }
    else _parse_binom_expr_error(str);

    while ( *beg && isspace(*beg) ) beg++;
    if ( !*beg ) _parse_binom_expr_error(str);

    args->binom_val = strtod(beg, &end);
    while ( *end && isspace(*end) ) end++;
    if ( *end ) _parse_binom_expr_error(str);

    args->tgt_mask |= GT_BINOM;
    return;
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t*) calloc(1,sizeof(args_t));
    args->in_hdr  = in;
    args->out_hdr = out;

    char *tmp;
    int c;
    static struct option loptions[] =
    {
        {"include",required_argument,NULL,'i'},
        {"exclude",required_argument,NULL,'e'},
        {"new-gt",required_argument,NULL,'n'},
        {"target-gt",required_argument,NULL,'t'},
        {"seed",required_argument,NULL,'s'},
        {NULL,0,NULL,0}
    };
    while ((c = getopt_long(argc, argv, "?hn:t:i:e:s:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 's':
                args->rand_seed = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse: -s %s\n",optarg);
                break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'n': args->new_mask = 0;
                if ( strchr(optarg,'.') ) args->new_mask |= GT_MISSING;
                if ( strchr(optarg,'0') ) args->new_mask |= GT_REF;
                if ( strchr(optarg,'m') ) args->new_mask |= GT_MINOR;
                if ( strchr(optarg,'M') ) args->new_mask |= GT_MAJOR;
                if ( strchr(optarg,'X') ) args->new_mask |= GT_X_VAF;
                if ( strchr(optarg,'p') ) args->new_mask |= GT_PHASED;
                if ( strchr(optarg,'u') ) args->new_mask |= GT_UNPHASED;
                if ( !strncmp(optarg,"c:",2) ) { args->new_mask |= GT_CUSTOM; args->custom.gt_str = optarg; }
                if ( args->new_mask==0 ) error("Unknown parameter to --new-gt: %s\n", optarg);
                break;
            case 't':
                if ( !strcmp(optarg,".") ) args->tgt_mask |= GT_MISSING|GT_PARTIAL;
                if ( !strcmp(optarg,"./x") ) args->tgt_mask |= GT_PARTIAL;
                if ( !strcmp(optarg,"./.") ) args->tgt_mask |= GT_MISSING;
                if ( !strcmp(optarg,"a") ) args->tgt_mask |= GT_ALL;
                if ( !strcmp(optarg,"q") ) args->tgt_mask |= GT_QUERY;
                if ( !strcmp(optarg,"?") ) args->tgt_mask |= GT_QUERY;        // for backward compatibility
                if ( !strncmp(optarg,"r:",2) )
                {
                    args->rand_frac = strtod(optarg+2,&tmp);
                    if ( *tmp ) error("Could not parse: -t %s\n", optarg);
                    if ( args->rand_frac<=0 || args->rand_frac>=1 ) error("Expected value between 0 and 1 with -t\n");
                    args->tgt_mask |= GT_RAND;
                }
                if ( strchr(optarg,'b') ) parse_binom_expr(args, strchr(optarg,'b'));
                if ( args->tgt_mask==0 ) error("Unknown parameter to --target-gt: %s\n", optarg);
                break;
            case 'h':
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }

    if ( !args->new_mask ) error("Expected -n option\n");
    if ( !args->tgt_mask ) error("Expected -t option\n");
    if ( args->tgt_mask & GT_RAND )
    {
        if ( args->tgt_mask==GT_RAND ) args->tgt_mask |= GT_ALL;
        hts_srand48(args->rand_seed);
    }

    if ( args->new_mask & GT_MISSING ) args->new_gt = bcf_gt_missing;
    if ( args->new_mask & GT_REF ) args->new_gt = args->new_mask&GT_PHASED ? bcf_gt_phased(0) : bcf_gt_unphased(0);
    if ( args->new_mask & GT_CUSTOM )
    {
        char *ptr = args->custom.gt_str + 2;
        while ( *ptr )
        {
            int allele;
            if ( *ptr=='m' ) { allele = MINOR_ALLELE; args->new_mask |= GT_MINOR; }
            else if ( *ptr=='M' ) { allele = MAJOR_ALLELE; args->new_mask |= GT_MAJOR; }
            else if ( *ptr=='X' ) { allele = X_VAF_ALLELE; args->new_mask |= GT_X_VAF; }
            else if ( *ptr=='/' || *ptr=='|' )
            {
                if ( !args->custom.ploidy ) error("Could not parse the genotype: %s\n",args->custom.gt_str);
                args->custom.phased[args->custom.ploidy-1] = *ptr=='|' ? 1 : 0;
                ptr++;
                continue;
            }
            else
            {
                char *end;
                allele = strtol(ptr,&end,10);
                if ( end==ptr || (*end && *end!='/' && *end!='|') ) error("Could not parse the genotype: %s\n",args->custom.gt_str);
                ptr = end - 1;
            }
            args->custom.ploidy++;
            args->custom.gt = realloc(args->custom.gt,sizeof(*args->custom.gt)*args->custom.ploidy);
            args->custom.phased = realloc(args->custom.phased,sizeof(*args->custom.phased)*args->custom.ploidy);
            args->custom.gt[args->custom.ploidy-1] = allele;
            args->custom.phased[args->custom.ploidy-1] = 0;
            ptr++;
        }
        // the phasing sign comes before the allele, move to the correct one
        int i;
        for (i=args->custom.ploidy-1; i>0; i--)
            args->custom.phased[i] = args->custom.phased[i-1];
    }

    if ( args->filter_str  && !(args->tgt_mask&GT_QUERY) ) error("Expected -tq with -i/-e\n");
    if ( !args->filter_str && args->tgt_mask&GT_QUERY ) error("Expected -i/-e with -tq\n");
    if ( args->filter_str ) args->filter = filter_init(in,args->filter_str);

    // Check the existence of FORMAT/GT tag, add it if not present
    int id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"GT");
    if ( !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,id) )
        bcf_hdr_printf(args->out_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if ( (args->new_mask & GT_X_VAF) && !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"AD")) )
        error("Error: the FORMAT/AD annotation does exist, cannot run with --new-gt %s\n",args->custom.gt_str);

    return 0;
}

// sets GT phase for a single sample, ngts is the ploidy
static inline int phase_gt(int32_t *ptr, int ngts)
{
    int j, changed = 0;
    for (j=0; j<ngts; j++)
    {
        if ( ptr[j]==bcf_int32_vector_end ) break;
        if ( bcf_gt_is_phased(ptr[j]) ) continue;
        ptr[j] = bcf_gt_phased(bcf_gt_allele(ptr[j]));    // add phasing; this may need a fix, I think the flag should be set for one allele only?
        changed++;
    }
    return changed;
}

// unphase GT for a single sample, ngts is the ploidy
static inline int unphase_gt(int32_t *ptr, int ngts)
{
    int j, changed = 0;
    for (j=0; j<ngts; j++)
    {
        if ( ptr[j]==bcf_int32_vector_end ) break;
        if ( !bcf_gt_is_phased(ptr[j]) ) continue;
        ptr[j] = bcf_gt_unphased(bcf_gt_allele(ptr[j]));    // remove phasing
        changed++;
    }

    // insertion sort
    int k, l;
    for (k=1; k<j; k++)
    {
        int32_t x = ptr[k];
        l = k;
        while ( l>0 && ptr[l-1]>x )
        {
            ptr[l] = ptr[l-1];
            l--;
        }
        ptr[l] = x;
    }
    return changed;
}

// sets GT for a single sample, ngts is the ploidy, allele
static inline int set_gt(int32_t *ptr, int ngts, int allele)
{
    int j, changed = 0;
    for (j=0; j<ngts; j++)
    {
        if ( ptr[j]==bcf_int32_vector_end ) break;
        if ( ptr[j] != allele ) changed++;
        ptr[j] = allele;
    }
    return changed;
}

// sets GT for a single sample using custom.[mMX]_allele; ngts is the ploidy, nals is the number of REF+ALT alleles
static inline int set_gt_custom(args_t *args, int32_t *ptr, int ngts, int nals)
{
    int i, changed = 0, new_allele = 0;
    for (i=0; i<args->custom.ploidy; i++)
    {
        if ( args->custom.gt[i]==MINOR_ALLELE ) new_allele = args->custom.m_allele;
        else if ( args->custom.gt[i]==MAJOR_ALLELE ) new_allele = args->custom.M_allele;
        else if ( args->custom.gt[i]==X_VAF_ALLELE ) new_allele = args->custom.x_vaf_allele==bcf_gt_missing ? nals : bcf_gt_allele(args->custom.x_vaf_allele);
        else new_allele = args->custom.gt[i];
        if ( new_allele >= nals ) // cannot set, the requested index is bigger than there are alleles in ALT
            new_allele = bcf_gt_missing;
        else
            new_allele = args->custom.phased[i] ? bcf_gt_phased(new_allele) : bcf_gt_unphased(new_allele);
        if ( ptr[i]==new_allele ) continue;
        ptr[i] = new_allele;
        changed = 1;
    }
    while ( i<ngts )
    {
        if ( !changed && ptr[i]!=bcf_int32_vector_end ) changed = 1;
        ptr[i++] = bcf_int32_vector_end;
    }
    return changed;
}

static inline int random_draw(args_t *args)
{
    return hts_drand48() > args->rand_frac ? 1 : 0; // reversed random draw
}

bcf1_t *process(bcf1_t *rec)
{
    if ( !rec->n_sample ) return rec;

    int ngts = bcf_get_genotypes(args->in_hdr, rec, &args->gts, &args->mgts);
    ngts /= rec->n_sample;
    int i, j, changed = 0;

    if ( args->new_mask & GT_CUSTOM )
    {
        if ( args->custom.ploidy > ngts ) // increased ploidy, expand the array
        {
            if ( args->mgts < args->custom.ploidy * rec->n_sample )
            {
                args->mgts = args->custom.ploidy * rec->n_sample;
                args->gts  = (int32_t*)realloc(args->gts,args->mgts*sizeof(*args->gts));
            }
            for (i=0; i<rec->n_sample; i++)
            {
                int32_t *src = args->gts + (rec->n_sample-i-1)*ngts;
                int32_t *dst = args->gts + (rec->n_sample-i-1)*args->custom.ploidy;
                for (j=0; j<ngts; j++) dst[j] = src[j];
                for (; j<args->custom.ploidy; j++) dst[j] = bcf_int32_vector_end;
            }
            ngts = args->custom.ploidy;
        }
    }

    // Calculating allele frequency for each allele and determining major allele
    // only do this if use_major is true
    if ( args->new_mask & GT_MAJOR )
    {
        int maxAC = -1, majorAllele = -1;
        hts_expand(int,rec->n_allele,args->marr,args->arr);
        int ret = bcf_calc_ac(args->in_hdr,rec,args->arr,BCF_UN_FMT);
        if ( ret<= 0 )
            error("Could not calculate allele count at %s:%"PRId64"\n", bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);

        for (i=0; i < rec->n_allele; ++i)
        {
            if (args->arr[i] > maxAC)
            {
                maxAC = args->arr[i];
                majorAllele = i;
            }
        }

        // replacing new_gt by major allele
        args->new_gt = args->new_mask & GT_PHASED ?  bcf_gt_phased(majorAllele) : bcf_gt_unphased(majorAllele);
        if ( args->new_mask & GT_CUSTOM ) args->custom.M_allele = majorAllele;
    }
    if ( args->new_mask & GT_MINOR )
    {
        hts_expand(int,rec->n_allele,args->marr,args->arr);
        int ret = bcf_calc_ac(args->in_hdr,rec,args->arr,BCF_UN_FMT);
        if ( ret<= 0 )
            error("Could not calculate allele count at %s:%"PRId64"\n", bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);

        int imax = 0;
        for (i=1; i<rec->n_allele; i++)
            if ( args->arr[imax] < args->arr[i] ) imax = i;
        int imax2 = imax>0 ? 0 : (rec->n_allele>1 ? 1 : 0);
        for (i=0; i<rec->n_allele; i++)
            if ( i!=imax && args->arr[imax2] < args->arr[i] ) imax2 = i;

        // replacing new_gt by major allele
        args->new_gt = args->new_mask & GT_PHASED ?  bcf_gt_phased(imax2) : bcf_gt_unphased(imax2);
        if ( args->new_mask & GT_CUSTOM ) args->custom.m_allele = imax2;
    }
    if ( args->new_mask & GT_X_VAF )
    {
        if ( bcf_get_format_int32(args->in_hdr,rec,"AD",&args->xarr,&args->mxarr)==rec->n_allele*rec->n_sample )
        {
            int32_t *src = args->xarr, *dst = args->xarr;
            for (i=0; i<rec->n_sample; i++)
            {
                int jmax = -1;
                for (j=0; j<rec->n_allele; j++)
                {
                    if ( src[j]==bcf_int32_vector_end ) break;
                    if ( src[j]==bcf_int32_missing  ) continue;
                    if ( j==0 || src[jmax] < src[j] ) jmax = j;
                }
                *dst = jmax==-1 ? bcf_gt_missing : bcf_gt_unphased(jmax);
                src += rec->n_allele;
                dst++;
            }
        }
        else return rec;
    }

    int nbinom = 0;
    if ( args->tgt_mask & GT_BINOM )
    {
        nbinom = bcf_get_format_int32(args->in_hdr, rec, args->binom_tag, &args->iarr, &args->miarr);
        if ( nbinom<0 ) nbinom = 0;
        nbinom /= rec->n_sample;
    }

    // replace gts
    if ( nbinom && ngts>=2 )    // only diploid genotypes are considered: higher ploidy ignored further, haploid here
    {
        if ( args->filter ) filter_test(args->filter,rec,(const uint8_t **)&args->smpl_pass);
        for (i=0; i<rec->n_sample; i++)
        {
            if ( args->smpl_pass )
            {
                if ( !args->smpl_pass[i] && args->filter_logic==FLT_INCLUDE ) continue;
                if (  args->smpl_pass[i] && args->filter_logic==FLT_EXCLUDE ) continue;
            }
            int32_t *ptr = args->gts + i*ngts;
            if ( bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1]) || ptr[1]==bcf_int32_vector_end ) continue;
            if ( ptr[0]==ptr[1] ) continue; // a hom
            int ia = bcf_gt_allele(ptr[0]);
            int ib = bcf_gt_allele(ptr[1]);
            if ( ia>=nbinom || ib>=nbinom )
                error("The sample %s has incorrect number of %s fields at %s:%"PRId64"\n",
                        args->in_hdr->samples[i],args->binom_tag,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);

            double prob = calc_binom_two_sided(args->iarr[i*nbinom+ia],args->iarr[i*nbinom+ib],0.5);
            if ( !args->binom_cmp(prob,args->binom_val) ) continue;
            if ( args->tgt_mask&GT_RAND && random_draw(args) ) continue;

            if ( args->new_mask&GT_UNPHASED )
                changed += unphase_gt(ptr, ngts);
            else if ( args->new_mask==GT_PHASED )
                changed += phase_gt(ptr, ngts);
            else if ( args->new_mask & GT_CUSTOM )
            {
                if ( args->new_mask & GT_X_VAF ) args->custom.x_vaf_allele = args->xarr[i];
                changed += set_gt_custom(args, ptr, ngts, rec->n_allele);
            }
            else if ( args->new_mask & GT_X_VAF )
                changed += set_gt(ptr, ngts, args->xarr[i]);
            else
                changed += set_gt(ptr, ngts, args->new_gt);
        }
    }
    else if ( args->tgt_mask&GT_QUERY )
    {
        int pass_site = filter_test(args->filter,rec,(const uint8_t **)&args->smpl_pass);
        if ( args->filter_logic==FLT_EXCLUDE )
        {
            // -i can include a site but exclude a sample, -e exclude a site but include a sample
            if ( pass_site )
            {
                if ( !args->smpl_pass ) return rec;
                pass_site = 0;
                for (i=0; i<rec->n_sample; i++)
                {
                    if ( args->smpl_pass[i] ) args->smpl_pass[i] = 0;
                    else { args->smpl_pass[i] = 1; pass_site = 1; }
                }
                if ( !pass_site ) return rec;
            }
            else if ( args->smpl_pass )
                for (i=0; i<rec->n_sample; i++) args->smpl_pass[i] = 1;
        }
        else if ( !pass_site ) return rec;

        for (i=0; i<rec->n_sample; i++)
        {
            if ( args->smpl_pass && !args->smpl_pass[i] ) continue;
            if ( args->tgt_mask&GT_RAND && random_draw(args) ) continue;
            if ( args->new_mask&GT_UNPHASED )
                changed += unphase_gt(args->gts + i*ngts, ngts);
            else if ( args->new_mask==GT_PHASED )
                changed += phase_gt(args->gts + i*ngts, ngts);
            else if ( args->new_mask & GT_CUSTOM )
            {
                if ( args->new_mask & GT_X_VAF ) args->custom.x_vaf_allele = args->xarr[i];
                changed += set_gt_custom(args, args->gts + i*ngts, ngts, rec->n_allele);
            }
            else if ( args->new_mask & GT_X_VAF )
                changed += set_gt(args->gts + i*ngts, ngts, args->xarr[i]);
            else
                changed += set_gt(args->gts + i*ngts, ngts, args->new_gt);
        }
    }
    else
    {
        for (i=0; i<rec->n_sample; i++)
        {
            int ploidy = 0, nmiss = 0;
            int32_t *ptr = args->gts + i*ngts;
            for (j=0; j<ngts; j++)
            {
                if ( ptr[j]==bcf_int32_vector_end ) break;
                ploidy++;
                if ( bcf_gt_is_missing(ptr[j]) ) nmiss++;
            }

            int do_set = 0;
            if ( args->tgt_mask&GT_ALL ) do_set = 1;
            else if ( args->tgt_mask&GT_PARTIAL && nmiss ) do_set = 1;
            else if ( args->tgt_mask&GT_MISSING && ploidy==nmiss ) do_set = 1;

            if ( !do_set ) continue;
            if ( args->tgt_mask&GT_RAND && random_draw(args) ) continue;

            if ( args->new_mask&GT_UNPHASED )
                changed += unphase_gt(ptr, ngts);
            else if ( args->new_mask==GT_PHASED )
                changed += phase_gt(ptr, ngts);
            else if ( args->new_mask & GT_CUSTOM )
            {
                if ( args->new_mask & GT_X_VAF ) args->custom.x_vaf_allele = args->xarr[i];
                changed += set_gt_custom(args, args->gts + i*ngts, ngts, rec->n_allele);
            }
            else if ( args->new_mask & GT_X_VAF )
                changed += set_gt(ptr, ngts, args->xarr[i]);
            else
                changed += set_gt(ptr, ngts, args->new_gt);
        }
    }
    args->nchanged += changed;
    if ( changed )
    {
        int ret = bcf_update_genotypes(args->out_hdr, rec, args->gts, ngts*rec->n_sample);
        if ( ret!=0 ) error("Error: failed to update genotypes at %s:%"PRIhts_pos"\n",bcf_seqname(args->in_hdr,rec),rec->pos+1);
    }
    return rec;
}

void destroy(void)
{
    fprintf(stderr,"Filled %"PRId64" alleles\n", args->nchanged);
    free(args->custom.gt);
    free(args->custom.phased);
    free(args->binom_tag);
    if ( args->filter ) filter_destroy(args->filter);
    free(args->arr);
    free(args->iarr);
    free(args->xarr);
    free(args->gts);
    free(args);
}


