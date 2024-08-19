/* The MIT License

   Copyright (c) 2015-2023 Genome Research Ltd.

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
#include <strings.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/vcf.h>
#include <htslib/khash_str2int.h>
#include <htslib/kbitset.h>
#include "bcftools.h"
#include "filter.h"

#define SET_AN      (1<<0)
#define SET_AC      (1<<1)
#define SET_AC_Hom  (1<<2)
#define SET_AC_Het  (1<<3)
#define SET_AC_Hemi (1<<4)
#define SET_AF      (1<<5)
#define SET_NS      (1<<6)
#define SET_MAF     (1<<7)
#define SET_HWE     (1<<8)
#define SET_ExcHet  (1<<9)
#define SET_FUNC    (1<<10)
#define SET_END     (1<<11)
#define SET_TYPE    (1<<12)
#define SET_VAF     (1<<13)
#define SET_VAF1    (1<<14)

typedef struct _args_t args_t;
typedef struct _ftf_t ftf_t;
typedef struct _pop_t pop_t;
typedef int (*fill_tag_f)(args_t *, bcf1_t *, pop_t *, ftf_t *);
struct _ftf_t
{
    char *src_tag, *dst_tag;
    fill_tag_f func;
    float *fval;
    int32_t *ival;
    int nfval, nival;
    int type, len, cnt;     // VCF type (e.g. BCF_HT_REAL); length (BCF_VL_FIXED); count (e.g. 1)
    int info;               // filling INFO or FORMAT?
    filter_t *filter;
};

typedef struct
{
    int nhom, nhet, nhemi, nac;
}
counts_t;

struct _pop_t
{
    int ns;
    int ncounts, mcounts;
    counts_t *counts;
    char *name, *suffix;
    int nsmpl, *smpl;
    ftf_t *ftf;
    int nftf;
};

struct _args_t
{
    bcf_hdr_t *in_hdr, *out_hdr;
    uint32_t tags, warned;
    int npop, drop_missing, gt_id;
    pop_t *pop, **smpl2pop;
    float *farr;
    int32_t *iarr, niarr, miarr, nfarr, mfarr, unpack;
    double *hwe_probs;
    int mhwe_probs;
    kstring_t str;
    kbitset_t *bset;
};

static args_t *args;

const char *about(void)
{
    return "Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS; FORMAT/VAF and more.\n";
}

const char *usage(void)
{
    return
        "\n"
        "About: Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS\n"
        "       FORMAT tag VAF, custom INFO/TAG=func(FMT/TAG).\n"
        "       See examples below, run with -l for detailed description.\n"
        "Usage: bcftools +fill-tags [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -d, --drop-missing          do not count half-missing genotypes \"./1\" as hemizygous\n"
        "   -l, --list-tags             list available tags with description\n"
        "   -t, --tags LIST             list of output tags, \"all\" for all tags\n"
        "   -S, --samples-file FILE     list of samples (first column) and comma-separated list of populations (second column)\n"
        "\n"
        "Example:\n"
        "   # Print a detailed list of available tags\n"
        "   bcftools +fill-tags -- -l\n"
        "\n"
        "   # Fill INFO/AN and INFO/AC\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t AN,AC\n"
        "\n"
        "   # Fill (almost) all available tags\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t all\n"
        "\n"
        "   # Calculate HWE for sample groups (possibly multiple) read from a file\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -S sample-group.txt -t HWE\n"
        "\n"
        "   # Calculate total read depth (INFO/DP) from per-sample depths (FORMAT/DP)\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t 'DP:1=int(sum(FORMAT/DP))'\n"
        "\n"
        "   # Calculate per-sample read depth (FORMAT/DP) from per-sample allelic depths (FORMAT/AD)\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t 'FORMAT/DP:1=int(smpl_sum(FORMAT/AD))'\n"
        "\n"
        "   # Annotate with allelic fraction\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t FORMAT/VAF\n"
        "\n";
}

void parse_samples(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    void *pop2i = khash_str2int_init();
    void *smpli = khash_str2int_init();
    kstring_t str = {0,0,0};

    int moff = 0, *off = NULL, nsmpl = 0;
    while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 )
    {
        // NA12400 GRP1
        // NA18507 GRP1,GRP2
        char *pop_names = str.s + str.l - 1;
        while ( pop_names >= str.s && isspace(*pop_names) ) pop_names--;
        if ( pop_names <= str.s ) error("Could not parse the file: %s\n", str.s);
        pop_names[1] = 0;   // trailing spaces
        while ( pop_names >= str.s && !isspace(*pop_names) ) pop_names--;
        if ( pop_names <= str.s ) error("Could not parse the file: %s\n", str.s);

        char *smpl = pop_names++;
        while ( smpl >= str.s && isspace(*smpl) ) smpl--;
        if ( smpl <= str.s+1 ) error("Could not parse the file: %s\n", str.s);
        smpl[1] = 0;
        smpl = str.s;

        int ismpl = bcf_hdr_id2int(args->in_hdr,BCF_DT_SAMPLE,smpl);
        if ( ismpl<0 )
        {
            fprintf(stderr,"Warning: The sample not present in the VCF: %s\n",smpl);
            continue;
        }
        if ( khash_str2int_has_key(smpli,smpl) )
        {
            fprintf(stderr,"Warning: The sample is listed twice in %s: %s\n",fname,smpl);
            continue;
        }
        khash_str2int_inc(smpli,strdup(smpl));

        int i,npops = ksplit_core(pop_names,',',&moff,&off);
        for (i=0; i<npops; i++)
        {
            char *pop_name = &pop_names[off[i]];
            if ( !khash_str2int_has_key(pop2i,pop_name) )
            {
                pop_name = strdup(pop_name);
                khash_str2int_set(pop2i,pop_name,args->npop);
                args->npop++;
                args->pop = (pop_t*) realloc(args->pop,args->npop*sizeof(*args->pop));
                memset(args->pop+args->npop-1,0,sizeof(*args->pop));
                args->pop[args->npop-1].name = pop_name;
                args->pop[args->npop-1].suffix = (char*)malloc(strlen(pop_name)+2);
                memcpy(args->pop[args->npop-1].suffix+1,pop_name,strlen(pop_name)+1);
                args->pop[args->npop-1].suffix[0] = '_';
            }
            int ipop = 0;
            khash_str2int_get(pop2i,pop_name,&ipop);
            pop_t *pop = &args->pop[ipop];
            pop->nsmpl++;
            pop->smpl = (int*) realloc(pop->smpl,pop->nsmpl*sizeof(*pop->smpl));
            pop->smpl[pop->nsmpl-1] = ismpl;
        }
        nsmpl++;
    }

    if ( nsmpl != bcf_hdr_nsamples(args->in_hdr) )
        fprintf(stderr,"Warning: %d samples in the list, %d samples in the VCF.\n", nsmpl,bcf_hdr_nsamples(args->in_hdr));

    if ( !args->npop ) error("No populations given?\n");

    khash_str2int_destroy(pop2i);
    khash_str2int_destroy_free(smpli);
    free(str.s);
    free(off);
    if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,fname);
}

void init_pops(args_t *args)
{
    int i,j, nsmpl;

    // add the population "ALL", which is a summary population for all samples
    args->npop++;
    args->pop = (pop_t*) realloc(args->pop,args->npop*sizeof(*args->pop));
    memset(args->pop+args->npop-1,0,sizeof(*args->pop));
    args->pop[args->npop-1].name   = strdup("");
    args->pop[args->npop-1].suffix = strdup("");

    nsmpl = bcf_hdr_nsamples(args->in_hdr);
    args->smpl2pop = (pop_t**) calloc(nsmpl*(args->npop+1),sizeof(pop_t*));
    for (i=0; i<nsmpl; i++)
        args->smpl2pop[i*(args->npop+1)] = &args->pop[args->npop-1];

    for (i=0; i<args->npop; i++)
    {
        for (j=0; j<args->pop[i].nsmpl; j++)
        {
            int ismpl = args->pop[i].smpl[j];
            pop_t **smpl2pop = &args->smpl2pop[ismpl*(args->npop+1)];
            while (*smpl2pop) smpl2pop++;
            *smpl2pop = &args->pop[i];
        }
    }
}

void ftf_destroy(pop_t *pop)
{
    int i;
    for (i=0; i<pop->nftf; i++)
    {
        ftf_t *ftf = &pop->ftf[i];
        free(ftf->src_tag);
        free(ftf->dst_tag);
        free(ftf->fval);
        free(ftf->ival);
        if ( ftf->filter ) filter_destroy(ftf->filter);
    }
    free(pop->ftf);
}

#define float_set_from_double(dst, src) \
    if (bcf_double_is_missing_or_vector_end(src)) bcf_float_set_missing(dst); else (dst) = (src)

static inline int32_t int32_from_double(double src)
{
    if (bcf_double_is_missing_or_vector_end(src)) return bcf_int32_missing;
    else return (int)src;
}

int ftf_filter_expr(args_t *args, bcf1_t *rec, pop_t *pop, ftf_t *ftf)
{
    int j,k, nval, nval1;

    args->str.l = 0;
    ksprintf(&args->str, "%s%s", ftf->dst_tag,pop->suffix);

    filter_test(ftf->filter, rec, NULL);
    const double *val = filter_get_doubles(ftf->filter,&nval,&nval1);

    int ret;
    if ( ftf->info )
    {
        int nfill = ftf->len==BCF_VL_FIXED ? ftf->cnt : nval;
        int ncopy = nval < nfill ? nval : nfill;    // number of values available, the rest is filled with missing values
        if ( ftf->type==BCF_HT_REAL )
        {
            hts_expand(float,nfill,ftf->nfval,ftf->fval);
            for (j=0; j<ncopy; j++) float_set_from_double(ftf->fval[j], val[j]);
            for (; j<nfill; j++) bcf_float_set_missing(ftf->fval[j]);
            ret = bcf_update_info_float(args->out_hdr,rec,args->str.s,ftf->fval,nfill);
        }
        else
        {
            hts_expand(int32_t,nfill,ftf->nival,ftf->ival);
            for (j=0; j<ncopy; j++) ftf->ival[j] = int32_from_double(val[j]);
            for (; j<nfill; j++) ftf->ival[j] = bcf_int32_missing;
            ret = bcf_update_info_int32(args->out_hdr,rec,args->str.s,ftf->ival,nfill);
        }
    }
    else
    {
        int nfill = ftf->len==BCF_VL_FIXED ? ftf->cnt : nval1;
        int ncopy = nval1 < nfill ? nval1 : nfill;
        if ( ftf->type==BCF_HT_REAL )
        {
            hts_expand(float,nfill*rec->n_sample,ftf->nfval,ftf->fval);
            for (k=0; k<rec->n_sample; k++)
            {
                float *dst  = ftf->fval + k*nfill;
                const double *src = val + k*nval1;
                for (j=0; j<ncopy; j++) float_set_from_double(dst[j], src[j]);
                for (; j<nfill; j++) bcf_float_set_missing(dst[j]);
            }
            ret = bcf_update_format_float(args->out_hdr,rec,args->str.s,ftf->fval,nfill*rec->n_sample);
        }
        else
        {
            hts_expand(int32_t,nfill*rec->n_sample,ftf->nival,ftf->ival);
            for (k=0; k<rec->n_sample; k++)
            {
                int32_t *dst = ftf->ival + k*nfill;
                const double *src  = val + k*nval1;
                for (j=0; j<ncopy; j++) dst[j] = int32_from_double(src[j]);
                for (; j<nfill; j++) dst[j] = bcf_int32_missing;
            }
            ret = bcf_update_format_int32(args->out_hdr,rec,args->str.s,ftf->ival,nfill*rec->n_sample);
        }
    }
    if ( ret )
        error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
    return 0;
}

// This is implemented as a macro so the compiler can properly validate the
// printf format string.
#define hdr_append(args, fmt) \
do { \
    int i; \
    for (i=0; i<args->npop; i++) \
        bcf_hdr_printf(args->out_hdr, fmt, args->pop[i].suffix,*args->pop[i].name ? " in " : "",args->pop[i].name); \
} while (0)


int parse_func_pop(args_t *args, pop_t *pop, char *tag_expr, char *expr)
{
    pop->nftf++;
    pop->ftf = (ftf_t *)realloc(pop->ftf,sizeof(*pop->ftf)*pop->nftf);
    ftf_t *ftf = &pop->ftf[ pop->nftf - 1 ];
    memset(ftf,0,sizeof(ftf_t));
    ftf->type = BCF_HT_REAL;
    ftf->len  = BCF_VL_VAR;
    ftf->cnt  = 0;
    ftf->info = 1;

    char *end = strchr(tag_expr,'=');
    char *beg = tag_expr;
    if ( !strncasecmp("info/",tag_expr,5) ) beg += 5;
    else if ( !strncasecmp("format/",tag_expr,7) ) { beg += 7; ftf->info = 0; }
    else if ( !strncasecmp("fmt/",tag_expr,4) ) { beg += 4; ftf->info = 0; }
    ftf->dst_tag = (char*)calloc(end-beg+1,1);
    memcpy(ftf->dst_tag, beg, end-beg);
    char *tmp, *cnt = strchr(ftf->dst_tag,':');
    if ( cnt )
    {
        *cnt = 0;
        ftf->cnt = strtol(cnt+1, &tmp, 10);
        if ( *tmp ) error("Could not parse the expression: %s\n",tag_expr);
        ftf->len = BCF_VL_FIXED;
    }

    char *filter = NULL;
    if ( expr[strlen(expr)-1]==')' )
    {
        if ( !strncasecmp(expr,"int(",4) ) { filter = strdup(expr+4); ftf->type = BCF_HT_INT; }
        else if ( !strncasecmp(expr,"integer(",8) ) { filter = strdup(expr+8); ftf->type = BCF_HT_INT; }
        else if ( !strncasecmp(expr,"float(",6) ) { filter = strdup(expr+6); ftf->type = BCF_HT_REAL; }
    }
    if ( filter )
        filter[strlen(filter)-1] = 0;
    else
        filter = strdup(expr);

    ftf->func = ftf_filter_expr;
    uint8_t *samples = (uint8_t*)calloc(bcf_hdr_nsamples(args->in_hdr),1);

    args->str.l = 0;
    ksprintf(&args->str, "%s%s", ftf->dst_tag,pop->suffix);
    int id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,args->str.s);
    int update_hdr = 1;
    if ( ftf->info && bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_INFO,id) )
    {
        if ( bcf_hdr_id2length(args->in_hdr,BCF_HL_INFO,id)!=ftf->len )
            error("Error: the field INFO/%s already exists and its Number definition is different\n",args->str.s);
        if ( ftf->len==BCF_VL_FIXED && bcf_hdr_id2number(args->in_hdr,BCF_HL_INFO,id)!=ftf->cnt )
            error("Error: the field INFO/%s already exists and its Number definition is different from %d\n",args->str.s,ftf->cnt);
        if ( bcf_hdr_id2type(args->in_hdr,BCF_HT_INT,id)!=ftf->type )
            error("Error: the field INFO/%s already exists and its Type definition is different\n",args->str.s);
        update_hdr = 0;
    }
    else if ( !ftf->info && bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,id) )
    {
        if ( bcf_hdr_id2length(args->in_hdr,BCF_HL_FMT,id)!=ftf->len )
            error("Error: the field FORMAT/%s already exists and its Number definition is different\n",args->str.s);
        if ( ftf->len==BCF_VL_FIXED && bcf_hdr_id2number(args->in_hdr,BCF_HL_FMT,id)!=ftf->cnt )
            error("Error: the field FORMAT/%s already exists and its Number definition is different from %d\n",args->str.s,ftf->cnt);
        if ( bcf_hdr_id2type(args->in_hdr,BCF_HT_INT,id)!=ftf->type )
            error("Error: the field FORMAT/%s already exists and its Type definition is different\n",args->str.s);
        update_hdr = 0;
    }
    if ( update_hdr )
    {
        args->str.l = 0;
        ksprintf(&args->str, "##%s=<ID=%s%s,Number=", ftf->info?"INFO":"FORMAT",ftf->dst_tag,pop->suffix);
        if ( ftf->len==BCF_VL_FIXED ) kputw(ftf->cnt, &args->str);
        else kputc('.', &args->str);
        kputs(",Type=", &args->str);
        if ( ftf->type==BCF_HT_INT ) kputs("Integer", &args->str);
        else kputs("Float", &args->str);
        kputs(",Description=\"Added by +fill-tags expression ", &args->str);
        // escape quotes
        char *tmp = tag_expr;
        while ( *tmp )
        {
            if ( *tmp=='"' ) kputc('\\', &args->str);
            kputc(*tmp, &args->str);
            tmp++;
        }
        if ( *pop->name )
        {
            kputs(" in ", &args->str);
            kputs(pop->name, &args->str);
        }
        kputs("\">", &args->str);
        bcf_hdr_append(args->out_hdr, args->str.s);
    }
    ftf->filter = filter_init(args->in_hdr, filter);
    args->unpack |= filter_max_unpack(ftf->filter);
    memset(samples,0,sizeof(*samples)*bcf_hdr_nsamples(args->in_hdr));
    if ( *pop->name )
    {
        int j;
        for (j=0; j<pop->nsmpl; j++) samples[ pop->smpl[j] ] = 1;
        filter_set_samples(ftf->filter, samples);
    }
    free(samples);
    free(filter);
    return SET_FUNC;
}
int parse_func(args_t *args, char *tag_expr, char *expr)
{
    int i, ret = 0;
    for (i=0; i<args->npop; i++)
        ret |= parse_func_pop(args,&args->pop[i],tag_expr,expr);
    return ret;
}
uint32_t parse_tags(args_t *args, const char *str)
{
    if ( !args->in_hdr ) error("%s", usage());

    args->warned = 0;
    uint32_t flag = 0;
    int i, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags), *ptr;
    for(i=0; i<n_tags; i++)
    {
        if ( !strcasecmp(tags[i],"all") )
        {
            flag |= ~(SET_END|SET_TYPE);
            // include F_MISSING as part of 'all', which requires explicitly
            // initialising it as a filter expression not just setting a
            // bitfield flag.
            flag |= parse_func(args,"F_MISSING=F_MISSING","F_MISSING");
            args->warned = ~(SET_END|SET_TYPE);
            args->unpack |= BCF_UN_FMT;
        }
        else if ( !strcasecmp(tags[i],"AN") || !strcasecmp(tags[i],"INFO/AN") ) { flag |= SET_AN; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"AC") || !strcasecmp(tags[i],"INFO/AC")  ) { flag |= SET_AC; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"NS") || !strcasecmp(tags[i],"INFO/NS")  ) { flag |= SET_NS; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"AC_Hom") || !strcasecmp(tags[i],"INFO/AC_Hom")  ) { flag |= SET_AC_Hom; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"AC_Het") || !strcasecmp(tags[i],"INFO/AC_Het")  ) { flag |= SET_AC_Het; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"AC_Hemi") || !strcasecmp(tags[i],"INFO_Hemi")  ) { flag |= SET_AC_Hemi; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"AF") || !strcasecmp(tags[i],"INFO/AF")  ) { flag |= SET_AF; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"MAF") || !strcasecmp(tags[i],"INFO/MAF")  ) { flag |= SET_MAF; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"HWE") || !strcasecmp(tags[i],"INFO/HWE")  ) { flag |= SET_HWE; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"ExcHet") || !strcasecmp(tags[i],"INFO/ExcHet")  ) { flag |= SET_ExcHet; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"VAF") || !strcasecmp(tags[i],"FORMAT/VAF") ) { flag |= SET_VAF; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"VAF1") || !strcasecmp(tags[i],"FORMAT/VAF1") ) { flag |= SET_VAF1; args->unpack |= BCF_UN_FMT; }
        else if ( !strcasecmp(tags[i],"END") || !strcasecmp(tags[i],"INFO/END")  ) flag |= SET_END;
        else if ( !strcasecmp(tags[i],"TYPE") || !strcasecmp(tags[i],"INFO/TYPE")  ) flag |= SET_TYPE;
        else if ( !strcasecmp(tags[i],"F_MISSING") || !strcasecmp(tags[i],"INFO/F_MISSING")  ) { flag |= parse_func(args,"F_MISSING=F_MISSING","F_MISSING"); args->unpack |= BCF_UN_FMT; }
        else if ( (ptr=strchr(tags[i],'=')) ) { flag |= parse_func(args,tags[i],ptr+1);  args->unpack |= BCF_UN_FMT; }
        else
        {
            fprintf(stderr,"Error parsing \"--tags %s\": the tag \"%s\" is not supported\n", str,tags[i]);
            exit(1);
        }
        free(tags[i]);
    }
    if (n_tags) free(tags);
    return flag;
}

void list_tags(void)
{
    error(
        "INFO/AC        Number:A  Type:Integer  ..  Allele count in genotypes\n"
        "INFO/AC_Hom    Number:A  Type:Integer  ..  Allele counts in homozygous genotypes\n"
        "INFO/AC_Het    Number:A  Type:Integer  ..  Allele counts in heterozygous genotypes\n"
        "INFO/AC_Hemi   Number:A  Type:Integer  ..  Allele counts in hemizygous genotypes\n"
        "INFO/AF        Number:A  Type:Float    ..  Allele frequency from FMT/GT or AC,AN if FMT/GT is not present\n"
        "INFO/AN        Number:1  Type:Integer  ..  Total number of alleles in called genotypes\n"
        "INFO/ExcHet    Number:A  Type:Float    ..  Test excess heterozygosity; 1=good, 0=bad\n"
        "INFO/END       Number:1  Type:Integer  ..  End position of the variant\n"
        "INFO/F_MISSING Number:1  Type:Float    ..  Fraction of missing genotypes (all samples, experimental)\n"
        "INFO/HWE       Number:A  Type:Float    ..  HWE test (PMID:15789306); 1=good, 0=bad\n"
        "INFO/MAF       Number:1  Type:Float    ..  Frequency of the second most common allele\n"
        "INFO/NS        Number:1  Type:Integer  ..  Number of samples with data\n"
        "INFO/TYPE      Number:.  Type:String   ..  The record type (REF,SNP,MNP,INDEL,etc)\n"
        "FORMAT/VAF     Number:A  Type:Float    ..  The fraction of reads with the alternate allele, requires FORMAT/AD or ADF+ADR\n"
        "FORMAT/VAF1    Number:1  Type:Float    ..  The same as FORMAT/VAF but for all alternate alleles cumulatively\n"
        "TAG:Number=Type(EXPR)                  ..  Experimental support for user expressions such as DP:1=int(sum(DP))\n"
        "               If Number and Type are not given (e.g. DP=sum(DP)), variable number (Number=.) of floating point\n"
        "               values (Type=Float) will be used.\n"
        );
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t*) calloc(1,sizeof(args_t));
    args->in_hdr  = in;
    args->out_hdr = out;
    char *samples_fname = NULL, *tags_str = "all";
    static struct option loptions[] =
    {
        {"list-tags",0,0,'l'},
        {"drop-missing",0,0,'d'},
        {"tags",1,0,'t'},
        {"samples-file",1,0,'S'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?ht:dS:l",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'l': list_tags(); break;
            case 'd': args->drop_missing = 1; break;
            case 't': tags_str = optarg; break;
            case 'S': samples_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }

    if ( optind != argc ) error("%s",usage());
    if ( !in || !out ) error("%s",usage());

    args->gt_id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"GT");
    if ( args->unpack&BCF_UN_FMT && args->gt_id<0 ) error("Error: GT field is not present\n");

    if ( samples_fname ) parse_samples(args, samples_fname);
    init_pops(args);

    args->tags |= parse_tags(args,tags_str);

    if ( args->tags & SET_AN ) hdr_append(args, "##INFO=<ID=AN%s,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes%s%s\">");
    if ( args->tags & SET_AC ) hdr_append(args, "##INFO=<ID=AC%s,Number=A,Type=Integer,Description=\"Allele count in genotypes%s%s\">");
    if ( args->tags & SET_NS ) hdr_append(args, "##INFO=<ID=NS%s,Number=1,Type=Integer,Description=\"Number of samples with data%s%s\">");
    if ( args->tags & SET_AC_Hom ) hdr_append(args, "##INFO=<ID=AC_Hom%s,Number=A,Type=Integer,Description=\"Allele counts in homozygous genotypes%s%s\">");
    if ( args->tags & SET_AC_Het ) hdr_append(args, "##INFO=<ID=AC_Het%s,Number=A,Type=Integer,Description=\"Allele counts in heterozygous genotypes%s%s\">");
    if ( args->tags & SET_AC_Hemi ) hdr_append(args, "##INFO=<ID=AC_Hemi%s,Number=A,Type=Integer,Description=\"Allele counts in hemizygous genotypes%s%s\">");
    if ( args->tags & SET_AF ) hdr_append(args, "##INFO=<ID=AF%s,Number=A,Type=Float,Description=\"Allele frequency%s%s\">");
    if ( args->tags & SET_MAF ) hdr_append(args, "##INFO=<ID=MAF%s,Number=1,Type=Float,Description=\"Frequency of the second most common allele%s%s\">");
    if ( args->tags & SET_HWE ) hdr_append(args, "##INFO=<ID=HWE%s,Number=A,Type=Float,Description=\"HWE test%s%s (PMID:15789306); 1=good, 0=bad\">");
    if ( args->tags & SET_END ) bcf_hdr_printf(args->out_hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">");
    if ( args->tags & SET_TYPE ) bcf_hdr_printf(args->out_hdr, "##INFO=<ID=TYPE,Number=.,Type=String,Description=\"Variant type\">");
    if ( args->tags & SET_ExcHet ) hdr_append(args, "##INFO=<ID=ExcHet%s,Number=A,Type=Float,Description=\"Test excess heterozygosity%s%s; 1=good, 0=bad\">");
    if ( args->tags & SET_VAF ) bcf_hdr_append(args->out_hdr, "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"The fraction of reads with alternate allele (nALT/nSumAll)\">");
    if ( args->tags & SET_VAF1 ) bcf_hdr_append(args->out_hdr, "##FORMAT=<ID=VAF1,Number=1,Type=Float,Description=\"The fraction of reads with alternate alleles (nSumALT/nSumAll)\">");

    int i, max_unpack_bit = 0;
    for (i=0; i<=3; i++)
        if ( args->unpack & (1<<i) ) max_unpack_bit = i;
    for (i=0; i<=max_unpack_bit; i++)
        args->unpack |= 1<<i;

    return 0;
}

/*
    Wigginton 2005, PMID: 15789306

    nref .. number of reference alleles
    nalt .. number of alt alleles
    nhet .. number of het genotypes, assuming number of genotypes = (nref+nalt)*2

*/

void calc_hwe(args_t *args, int nref, int nalt, int nhet, float *p_hwe, float *p_exc_het)
{
    int ngt   = (nref+nalt) / 2;
    int nrare = nref < nalt ? nref : nalt;

    // sanity check: there is odd/even number of rare alleles iff there is odd/even number of hets
    if ( (nrare & 1) ^ (nhet & 1) ) error("nrare/nhet should be both odd or even: nrare=%d nref=%d nalt=%d nhet=%d\n",nrare,nref,nalt,nhet);
    if ( nrare < nhet ) error("Fewer rare alleles than hets? nrare=%d nref=%d nalt=%d nhet=%d\n",nrare,nref,nalt,nhet);
    if ( (nref+nalt) & 1 ) error("Expected diploid genotypes: nref=%d nalt=%d\n",nref,nalt);

    // initialize het probs
    hts_expand(double,nrare+1,args->mhwe_probs,args->hwe_probs);
    memset(args->hwe_probs, 0, sizeof(*args->hwe_probs)*(nrare+1));
    double *probs = args->hwe_probs;

    // start at midpoint
    int mid = (double)nrare * (nref + nalt - nrare) / (nref + nalt);

    // check to ensure that midpoint and rare alleles have same parity
    if ( (nrare & 1) ^ (mid & 1) ) mid++;

    int het = mid;
    int hom_r  = (nrare - mid) / 2;
    int hom_c  = ngt - het - hom_r;
    double sum = probs[mid] = 1.0;

    for (het = mid; het > 1; het -= 2)
    {
        probs[het - 2] = probs[het] * het * (het - 1.0) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0));
        sum += probs[het - 2];

        // 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
        hom_r++;
        hom_c++;
    }

    het = mid;
    hom_r = (nrare - mid) / 2;
    hom_c = ngt - het - hom_r;
    for (het = mid; het <= nrare - 2; het += 2)
    {
        probs[het + 2] = probs[het] * 4.0 * hom_r * hom_c / ((het + 2.0) * (het + 1.0));
        sum += probs[het + 2];

        // add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
        hom_r--;
        hom_c--;
    }

    for (het=0; het<nrare+1; het++) probs[het] /= sum;

    double prob = probs[nhet];
    for (het = nhet + 1; het <= nrare; het++) prob += probs[het];
    *p_exc_het = prob;

    prob = 0;
    for (het=0; het <= nrare; het++)
    {
        if ( probs[het] > probs[nhet]) continue;
        prob += probs[het];
    }
    if ( prob > 1 ) prob = 1;
    *p_hwe = prob;
}

static inline void set_counts(pop_t *pop, int is_half, int is_hom, int is_hemi, kbitset_t *bset)
{
    kbitset_iter_t itr;
    int i;
    kbs_start(&itr);
    while ((i = kbs_next(bset, &itr)) >= 0)
    {
        if ( is_half ) pop->counts[i].nac++;
        else if ( !is_hom ) pop->counts[i].nhet++;
        else if ( !is_hemi ) pop->counts[i].nhom += 2;
        else pop->counts[i].nhemi++;
    }
    pop->ns++;
}
static void clean_counts(pop_t *pop, int nals)
{
    pop->ns = 0;
    memset(pop->counts,0,sizeof(counts_t)*nals);
}
static int cmpfloat_desc(const void *a, const void *b)
{
    float fa = *((float*)a);
    float fb = *((float*)b);
    if ( fa<fb ) return 1;
    if ( fa>fb ) return -1;
    return 0;
}
static void process_fmt(bcf1_t *rec)
{
    int i,j, nsmpl = bcf_hdr_nsamples(args->in_hdr);;

    bcf_fmt_t *fmt_gt = NULL;
    for (i=0; i<rec->n_fmt; i++)
        if ( rec->d.fmt[i].id==args->gt_id ) { fmt_gt = &rec->d.fmt[i]; break; }
    if ( !fmt_gt ) return;    // no GT tag

    hts_expand(int32_t,rec->n_allele, args->miarr, args->iarr);
    hts_expand(float,rec->n_allele*2, args->mfarr, args->farr);
    for (i=0; i<args->npop; i++)
        hts_expand(counts_t,rec->n_allele,args->pop[i].mcounts, args->pop[i].counts);

    for (i=0; i<args->npop; i++)
        clean_counts(&args->pop[i], rec->n_allele);

    if ( kbs_resize(&args->bset, rec->n_allele) < 0 ) error("kbs_resize: failed to store %d bits\n", rec->n_allele);

    #define BRANCH_INT(type_t,vector_end) \
    { \
        for (i=0; i<nsmpl; i++) \
        { \
            type_t *p = (type_t*) (fmt_gt->p + i*fmt_gt->size); \
            int ial, nbits = 0, nals = 0, is_half, is_hom, is_hemi; \
            kbs_clear(args->bset); \
            for (ial=0; ial<fmt_gt->n; ial++) \
            { \
                if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                if ( bcf_gt_is_missing(p[ial]) ) continue; /* missing allele */ \
                int idx = bcf_gt_allele(p[ial]); \
                nals++; \
                \
                if ( idx >= rec->n_allele ) \
                    error("Incorrect allele (\"%d\") in %s at %s:%"PRId64"\n",idx,args->in_hdr->samples[i],bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1); \
                if ( !kbs_exists(args->bset, idx) ) nbits++; \
                kbs_insert(args->bset, idx); \
            } \
            if ( nals==0 ) continue; /* missing genotype */ \
            is_hom = nbits==1 ? 1 : 0; /* only one bit is set for homs */ \
            if ( nals!=ial ) \
            { \
                if ( args->drop_missing ) is_hemi = 0, is_half = 1; \
                else is_hemi = 1, is_half = 0; \
            } \
            else if ( nals==1 ) is_hemi = 1, is_half = 0; \
            else is_hemi = 0, is_half = 0; \
            pop_t **pop = &args->smpl2pop[i*(args->npop+1)]; \
            while ( *pop ) { set_counts(*pop,is_half,is_hom,is_hemi,args->bset); pop++; } \
        } \
    }
    switch (fmt_gt->type) {
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
        default: error("The GT type is not recognised: %d at %s:%"PRId64"\n",fmt_gt->type, bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1); break;
    }
    #undef BRANCH_INT

    if ( args->tags & SET_NS )
    {
        for (i=0; i<args->npop; i++)
        {
            args->str.l = 0;
            ksprintf(&args->str, "NS%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,&args->pop[i].ns,1)!=0 )
                error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
        }
    }
    if ( args->tags & SET_AN )
    {
        for (i=0; i<args->npop; i++)
        {
            pop_t *pop = &args->pop[i];
            int32_t an = 0;
            for (j=0; j<rec->n_allele; j++)
                an += pop->counts[j].nhet + pop->counts[j].nhom + pop->counts[j].nhemi + pop->counts[j].nac;

            args->str.l = 0;
            ksprintf(&args->str, "AN%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,&an,1)!=0 )
                error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
        }
    }
    if ( args->tags & (SET_AF | SET_MAF) )
    {
        for (i=0; i<args->npop; i++)
        {
            int32_t an = 0;
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                for (j=0; j<rec->n_allele; j++)
                {
                    args->farr[j] = pop->counts[j].nhet + pop->counts[j].nhom + pop->counts[j].nhemi + pop->counts[j].nac;
                    an += args->farr[j];
                }
                if ( an )
                    for (j=0; j<rec->n_allele; j++) args->farr[j] /= an;
                else
                    for (j=0; j<rec->n_allele; j++) bcf_float_set_missing(args->farr[j]);
            }
            if ( args->tags & SET_AF )
            {
                args->str.l = 0;
                ksprintf(&args->str, "AF%s", args->pop[i].suffix);
                if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,args->farr+1,rec->n_allele-1)!=0 )
                    error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
            }
            if ( rec->n_allele > 1 && args->tags & SET_MAF )
            {
                if ( an ) qsort(args->farr,rec->n_allele,sizeof(float),cmpfloat_desc);
                args->str.l = 0;
                ksprintf(&args->str, "MAF%s", args->pop[i].suffix);
                if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,args->farr+1,1)!=0 )
                    error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
            }
        }
    }
    if ( args->tags & SET_AC )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                for (j=0; j<rec->n_allele; j++)
                    args->iarr[j] = pop->counts[j].nhet + pop->counts[j].nhom + pop->counts[j].nhemi + pop->counts[j].nac;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr+1,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
        }
    }
    if ( args->tags & SET_AC_Het )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->iarr, 0, sizeof(*args->iarr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++)
                    args->iarr[j-1] += pop->counts[j].nhet;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC_Het%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
        }
    }
    if ( args->tags & SET_AC_Hom )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->iarr, 0, sizeof(*args->iarr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++)
                    args->iarr[j-1] += pop->counts[j].nhom;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC_Hom%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
        }
    }
    if ( args->tags & SET_AC_Hemi && rec->n_allele > 1 )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->iarr, 0, sizeof(*args->iarr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++)
                    args->iarr[j-1] += pop->counts[j].nhemi;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC_Hemi%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
        }
    }
    if ( args->tags & (SET_HWE|SET_ExcHet) )
    {
        for (i=0; i<args->npop; i++)
        {
            float *fhwe = args->farr;
            float *fexc_het = args->farr + rec->n_allele;
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->farr,  0, sizeof(*args->farr)*(2*rec->n_allele));
                int nref_tot = pop->counts[0].nhom;
                for (j=0; j<rec->n_allele; j++) nref_tot += pop->counts[j].nhet;   // NB this neglects multiallelic genotypes
                for (j=1; j<rec->n_allele; j++)
                {
                    int nref = nref_tot - pop->counts[j].nhet;
                    int nalt = pop->counts[j].nhet + pop->counts[j].nhom;
                    int nhet = pop->counts[j].nhet;
                    if ( nref>0 && nalt>0 )
                        calc_hwe(args, nref, nalt, nhet, &fhwe[j-1], &fexc_het[j-1]);
                    else
                        fhwe[j-1] = fexc_het[j-1] = 1;
                }
            }
            if ( args->tags & SET_HWE )
            {
                args->str.l = 0;
                ksprintf(&args->str, "HWE%s", args->pop[i].suffix);
                if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,fhwe,rec->n_allele-1)!=0 )
                    error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
            }
            if ( args->tags & SET_ExcHet )
            {
                args->str.l = 0;
                ksprintf(&args->str, "ExcHet%s", args->pop[i].suffix);
                if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,fexc_het,rec->n_allele-1)!=0 )
                    error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
            }
        }
    }
}
static void process_info_af(bcf1_t *rec)
{
    if ( !(args->tags & SET_AF) ) return;
    if ( bcf_hdr_nsamples(args->in_hdr) ) return;

    int n = bcf_get_info_int32(args->in_hdr,rec,"AN",&args->iarr,&args->miarr);
    if ( n!=1 ) return;
    int an = args->iarr[0];
    if ( !an ) return;

    n = bcf_get_info_int32(args->in_hdr,rec,"AC",&args->iarr,&args->miarr);
    if ( n!=rec->n_allele-1 ) return;

    hts_expand(float,n,args->mfarr,args->farr);
    int i;
    for (i=0; i<n; i++) args->farr[i] = (double)args->iarr[i]/an;

    if ( bcf_update_info_float(args->out_hdr,rec,"AF", args->farr, n)!=0 )
        error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
}
static void process_vaf(bcf1_t *rec, int mode)
{
    int nsmpl = bcf_hdr_nsamples(args->in_hdr);
    int nval  = args->niarr / nsmpl;
    int nval1 = (mode & SET_VAF) ? rec->n_allele - 1 : 1;
    int nfarr = nval1 * nsmpl;
    hts_expand(float,nfarr,args->mfarr,args->farr);
    int i,j;
    for (i=0; i<nsmpl; i++)
    {
        int32_t *src = args->iarr + i*nval;
        float *dst = args->farr + i*nval1;
        float sum = 0;
        for (j=0; j<nval; j++)
        {
            if ( src[j]==bcf_int32_missing || src[j]==bcf_int32_vector_end ) break;
            sum += src[j];
        }
        if ( j!=nval )
        {
            bcf_float_set_missing(dst[0]);
            for (j=1; j<nval1; j++) bcf_float_set_vector_end(dst[j]);
            continue;
        }
        if ( mode & SET_VAF1 )
        {
            *dst = sum ? (sum - src[0])/sum : 0;
            continue;
        }
        for (j=0; j<nval1; j++)
            dst[j] = sum ? src[j+1]/sum : 0;
    }
    if ( bcf_update_format_float(args->out_hdr,rec,(mode & SET_VAF) ? "VAF" : "VAF1", args->farr, nfarr)!=0 )
        error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
}
static void process_vaf_vaf1(bcf1_t *rec)
{
    if ( !(args->tags & (SET_VAF|SET_VAF1)) ) return;
    if ( rec->n_allele <= 1 ) return;

    args->niarr = bcf_get_format_int32(args->in_hdr, rec, "AD", &args->iarr, &args->miarr);
    if ( args->niarr <= 0 )
    {
        if ( !(args->warned&(SET_VAF|SET_VAF1)) )
            fprintf(stderr,"Warning: cannot add the VAF/VAF1 annotations, the required FORMAT/AD tag is missing at %s:%"PRIhts_pos".\n"
                           "         (This warning is printed only once.)\n",bcf_seqname(args->in_hdr,rec),rec->pos+1);
        args->warned |= SET_VAF|SET_VAF1;
        return;
    }

    int nsmpl = bcf_hdr_nsamples(args->in_hdr);
    if ( args->niarr != nsmpl*rec->n_allele ) return;   // incorrect number of values (possibly all missing)

    if ( args->tags & SET_VAF ) process_vaf(rec, SET_VAF);
    if ( args->tags & SET_VAF1 ) process_vaf(rec, SET_VAF1);
}

bcf1_t *process(bcf1_t *rec)
{
    int i,j;

    bcf_unpack(rec, args->unpack);
    for (i=0; i<args->npop; i++)
        for (j=0; j<args->pop[i].nftf; j++)
            args->pop[i].ftf[j].func(args, rec, &args->pop[i], &args->pop[i].ftf[j]);

    if ( args->unpack & BCF_UN_FMT )
    {
        process_fmt(rec);
        process_info_af(rec);
        process_vaf_vaf1(rec);
    }

    if ( args->tags & SET_END )
    {
        int32_t end = rec->pos + rec->rlen;
        if ( bcf_update_info_int32(args->out_hdr,rec,"END",&end,1)!=0 )
            error("Error occurred while updating INFO/END at %s:%"PRId64"\n", bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
    }
    if ( args->tags & SET_TYPE )
    {
        int type = bcf_get_variant_types(rec);
        args->str.l = 0;
        if ( type == VCF_REF ) kputs("REF,", &args->str);
        if ( type & VCF_SNP ) kputs("SNP,", &args->str);
        if ( type & VCF_MNP ) kputs("MNP,", &args->str);
        if ( type & VCF_INDEL ) kputs("INDEL,", &args->str);
        if ( type & VCF_OTHER ) kputs("OTHER,", &args->str);
        if ( type & VCF_BND ) kputs("BND,", &args->str);
        if ( type & VCF_OVERLAP ) kputs("OVERLAP,", &args->str);
        if ( args->str.l )
        {
            args->str.s[args->str.l-1] = 0;
            if ( bcf_update_info_string(args->out_hdr,rec,"TYPE",args->str.s)!=0 )
                error("Error occurred while updating INFO/TYPE at %s:%"PRId64"\n", bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
        }
    }
    return rec;
}

void destroy(void)
{
    int i;
    for (i=0; i<args->npop; i++)
    {
        free(args->pop[i].name);
        free(args->pop[i].suffix);
        free(args->pop[i].smpl);
        free(args->pop[i].counts);
        ftf_destroy(&args->pop[i]);
    }
    kbs_destroy(args->bset);
    free(args->str.s);
    free(args->pop);
    free(args->smpl2pop);
    free(args->iarr);
    free(args->farr);
    free(args->hwe_probs);
    free(args);
}



