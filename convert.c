/*  convert.c -- functions for converting between VCF/BCF and related formats.

    Copyright (C) 2013-2014 Genome Research Ltd.

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
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "convert.h"

#define T_CHROM   1
#define T_POS     2
#define T_ID      3
#define T_REF     4
#define T_ALT     5
#define T_QUAL    6
#define T_FILTER  7
#define T_INFO    8
#define T_FORMAT  9
#define T_SAMPLE  10
#define T_SEP     11
#define T_IS_TS   12
#define T_TYPE    13
#define T_MASK    14
#define T_GT      15
#define T_TGT     16
#define T_LINE    17
#define T_CHROM_POS_ID 18   // not publicly advertised
#define T_GT_TO_PROB3  19   // not publicly advertised
#define T_PL_TO_PROB3  20   // not publicly advertised
#define T_FIRST_ALT    21   // not publicly advertised

typedef struct _fmt_t
{
    int type, id, is_gt_field, ready, subscript;
    char *key;
    bcf_fmt_t *fmt;
    void (*handler)(convert_t *, bcf1_t *, struct _fmt_t *, int, kstring_t *);
}
fmt_t;

struct _convert_t
{
    fmt_t *fmt;
    int nfmt, mfmt;
    int nsamples, *samples;
    bcf_hdr_t *header;
    int max_unpack;
    char *format_str;
    bcf_srs_t *readers; // required only for %MASK
    int nreaders;
    void *dat;
    int ndat;
};


static void process_chrom(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(convert->header->id[BCF_DT_CTG][line->rid].key, str); }
static void process_pos(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputw(line->pos+1, str); }
static void process_id(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(line->d.id, str); }
static void process_ref(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(line->d.allele[0], str); }
static void process_alt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i;
    if ( line->n_allele==1 )
    {
        kputc('.', str);
        return;
    }
    for (i=1; i<line->n_allele; i++)
    {
        if ( i>1 ) kputc(',', str);
        kputs(line->d.allele[i], str);
    }
}
static void process_first_alt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( line->n_allele==1 )
        kputc('.', str);
    else
        kputs(line->d.allele[1], str);
}
static void process_qual(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( bcf_float_is_missing(line->qual) ) kputc('.', str);
    else ksprintf(str, "%g", line->qual);
}
static void process_filter(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i;
    if ( line->d.n_flt )
    {
        for (i=0; i<line->d.n_flt; i++)
        {
            if (i) kputc(';', str);
            kputs(convert->header->id[BCF_DT_ID][line->d.flt[i]].key, str);
        }
    }
    else kputc('.', str);
}
static inline int bcf_array_ivalue(void *bcf_array, int type, int idx)
{
    if ( type==BCF_BT_INT8 ) return ((int8_t*)bcf_array)[idx];
    if ( type==BCF_BT_INT16 ) return ((int16_t*)bcf_array)[idx];
    return ((int32_t*)bcf_array)[idx];
}
static void process_info(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i;
    for (i=0; i<line->n_info; i++)
        if ( line->d.info[i].key == fmt->id ) break;

    // output "." if the tag is not present
    if ( i==line->n_info )
    {
        kputc('.', str);
        return;
    }

    bcf_info_t *info = &line->d.info[i];

    // if this is a flag, output 1
    if ( info->len <=0 )
    {
        kputc('1', str);
        return;
    }

    if ( info->len == 1 )
    {
        if ( info->type == BCF_BT_FLOAT ) ksprintf(str, "%g", info->v1.f);
        else if ( info->type != BCF_BT_CHAR ) kputw(info->v1.i, str);
        else kputc(info->v1.i, str);
    }
    else if ( fmt->subscript >=0 )
    {
        if ( info->len <= fmt->subscript )
        {
            kputc('.', str);
            return;
        }
        if ( info->type == BCF_BT_FLOAT ) ksprintf(str, "%g", ((float*)info->vptr)[fmt->subscript]);
        else if ( info->type != BCF_BT_CHAR ) kputw(bcf_array_ivalue(info->vptr,info->type,fmt->subscript), str);
        else error("TODO: %s:%d .. info->type=%d\n", __FILE__,__LINE__, info->type);
    }
    else
        bcf_fmt_array(str, info->len, info->type, info->vptr);
}
static void init_format(convert_t *convert, bcf1_t *line, fmt_t *fmt)
{
    fmt->id = bcf_hdr_id2int(convert->header, BCF_DT_ID, fmt->key);
    if ( fmt->id==-1 ) error("Error: no such tag defined in the VCF header: FORMAT/%s\n", fmt->key);
    fmt->fmt = NULL;
    int i;
    for (i=0; i<(int)line->n_fmt; i++)
        if ( line->d.fmt[i].id==fmt->id ) { fmt->fmt = &line->d.fmt[i]; break; }
    fmt->ready = 1;
}
static void process_format(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(convert, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }
    else if ( fmt->subscript >=0 )
    {
        if ( fmt->fmt->n <= fmt->subscript )
        {
            kputc('.', str);
            return;
        }
        if ( fmt->fmt->type == BCF_BT_FLOAT ) ksprintf(str, "%g", ((float*)(fmt->fmt->p + isample*fmt->fmt->size))[fmt->subscript]);
        else if ( fmt->fmt->type != BCF_BT_CHAR ) kputw(bcf_array_ivalue(fmt->fmt->p+isample*fmt->fmt->size,fmt->fmt->type,fmt->subscript), str);
        else error("TODO: %s:%d .. fmt->type=%d\n", __FILE__,__LINE__, fmt->fmt->type);
    }
    else
        bcf_fmt_array(str, fmt->fmt->n, fmt->fmt->type, fmt->fmt->p + isample*fmt->fmt->size);
}
static void process_gt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(convert, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }
    bcf_format_gt(fmt->fmt, isample, str);
}
static void process_tgt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(convert, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }

    assert( fmt->fmt->type==BCF_BT_INT8 );

    int l;
    int8_t *x = (int8_t*)(fmt->fmt->p + isample*fmt->fmt->size); // FIXME: does not work with n_alt >= 64
    for (l = 0; l < fmt->fmt->n && x[l] != bcf_int8_vector_end; ++l)
    {
        if (l) kputc("/|"[x[l]&1], str);
        if (x[l]>>1)
        {
            int ial = (x[l]>>1) - 1;
            kputs(line->d.allele[ial], str);
        }
        else
            kputc('.', str);
    }
    if (l == 0) kputc('.', str);
}
static void process_sample(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    kputs(convert->header->samples[isample], str);
}
static void process_sep(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { if (fmt->key) kputs(fmt->key, str); }
static void process_is_ts(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int is_ts = 0;
    if ( bcf_get_variant_types(line) & (VCF_SNP|VCF_MNP) )
        is_ts = abs(bcf_acgt2int(*line->d.allele[0])-bcf_acgt2int(*line->d.allele[1])) == 2 ? 1 : 0;
    kputc(is_ts ? '1' : '0', str);
}
static void process_type(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int line_type = bcf_get_variant_types(line);
    int i = 0;
    if ( line_type == VCF_REF ) { kputs("REF", str); i++; }
    if ( line_type & VCF_SNP ) { if (i) kputc(',',str); kputs("SNP", str); i++; }
    if ( line_type & VCF_MNP ) { if (i) kputc(',',str); kputs("MNP", str); i++; }
    if ( line_type & VCF_INDEL ) { if (i) kputc(',',str); kputs("INDEL", str); i++; }
    if ( line_type & VCF_OTHER ) { if (i) kputc(',',str); kputs("OTHER", str); i++; }
}
static void process_line(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    vcf_format1(convert->header, line, str);
}
static void process_chrom_pos_id(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( line->d.id[0]!='.' || line->d.id[1] )
    {
        // ID is present
        kputs(line->d.id, str);
    }
    else
    {
        // use CHROM:POS instead of ID
        kputs(convert->header->id[BCF_DT_CTG][line->rid].key, str);
        kputc(':', str);
        kputw(line->pos+1, str);
    }
}
static void process_gt_to_prob3(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int m,n,i;

    m = convert->ndat / sizeof(int32_t);
    n = bcf_get_genotypes(convert->header,line,&convert->dat,&m);
    convert->ndat = m * sizeof(int32_t);

    if ( n<=0 )
    {
        // Throw an error or silently proceed?
        //
        // for (i=0; i<convert->nsamples; i++) kputs(" 0.33 0.33 0.33", str);
        // return;

        error("Error parsing GT tag at %s:%d\n", bcf_seqname(convert->header,line),line->pos+1);
    }

    n /= convert->nsamples;
    for (i=0; i<convert->nsamples; i++)
    {
        int32_t *ptr = (int32_t*)convert->dat + i*n;
        int j;
        for (j=0; j<n; j++)
            if ( ptr[j]==bcf_int32_vector_end ) break;

        if ( j==2 )
        {
            // diploid
            if ( bcf_gt_allele(ptr[0])!=bcf_gt_allele(ptr[1]) )
                kputs(" 0 1 0", str);       // HET
            else if ( bcf_gt_allele(ptr[0])==1 )
                kputs(" 0 0 1", str);       // ALT HOM, first ALT allele
            else
                kputs(" 1 0 0", str);       // REF HOM or something else than first ALT
        }
        else if ( j==1 )
        {
            // haploid
            if ( ptr[0]==bcf_gt_missing )
                kputs(" 0.5 0.0 0.5", str);
            else if ( bcf_gt_allele(ptr[0])==1 )
                kputs(" 0 0 1", str);       // first ALT allele
            else
                kputs(" 1 0 0", str);       // REF or something else than first ALT
        }
        else error("FIXME: not ready for ploidy %d\n", j);
    }
}
static void process_pl_to_prob3(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int m,n,i;

    m = convert->ndat / sizeof(int32_t);
    n = bcf_get_format_int32(convert->header,line,"PL",&convert->dat,&m);
    convert->ndat = m * sizeof(int32_t);

    if ( n<=0 )
    {
        // Throw an error or silently proceed?
        //
        // for (i=0; i<convert->nsamples; i++) kputs(" 0.33 0.33 0.33", str);
        // return;

        error("Error parsing PL tag at %s:%d\n", bcf_seqname(convert->header,line),line->pos+1);
    }

    n /= convert->nsamples;
    for (i=0; i<convert->nsamples; i++)
    {
        int32_t *ptr = (int32_t*)convert->dat + i*n;
        int j;
        float sum = 0;
        for (j=0; j<n; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            sum += pow(10,-0.1*ptr[j]);
        }
        if ( j==line->n_allele )
        {
            // haploid
            kputc(' ',str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[0])/sum);
            kputs(" 0 ", str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[1])/sum);
        }
        else
        {
            // diploid
            kputc(' ',str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[0])/sum);
            kputc(' ',str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[1])/sum);
            kputc(' ',str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[2])/sum);
        }
    }
}

static fmt_t *register_tag(convert_t *convert, int type, char *key, int is_gtf)
{
    convert->nfmt++;
    if ( convert->nfmt > convert->mfmt )
    {
        convert->mfmt += 10;
        convert->fmt   = (fmt_t*) realloc(convert->fmt, convert->mfmt*sizeof(fmt_t));
    }
    fmt_t *fmt = &convert->fmt[ convert->nfmt-1 ];
    fmt->type  = type;
    fmt->key   = key ? strdup(key) : NULL;
    fmt->is_gt_field = is_gtf;
    fmt->subscript = -1;

    // Allow non-format tags, such as CHROM, INFO, etc., to appear amongst the format tags.
    if ( key )
    {
        int id = bcf_hdr_id2int(convert->header, BCF_DT_ID, key);
        if ( fmt->type==T_FORMAT && !bcf_hdr_idinfo_exists(convert->header,BCF_HL_FMT,id) )
        {
            if ( !strcmp("CHROM",key) ) { fmt->type = T_CHROM; }
            else if ( !strcmp("POS",key) ) { fmt->type = T_POS; }
            else if ( !strcmp("ID",key) ) { fmt->type = T_ID; }
            else if ( !strcmp("REF",key) ) { fmt->type = T_REF; }
            else if ( !strcmp("ALT",key) ) { fmt->type = T_ALT; }
            else if ( !strcmp("FIRST_ALT",key) ) { fmt->type = T_FIRST_ALT; }
            else if ( !strcmp("QUAL",key) ) { fmt->type = T_QUAL; }
            else if ( !strcmp("FILTER",key) ) { fmt->type = T_FILTER; }
            else if ( !strcmp("_CHROM_POS_ID",key) ) { fmt->type = T_CHROM_POS_ID; }
            else if ( id>=0 && bcf_hdr_idinfo_exists(convert->header,BCF_HL_INFO,id) )
            {
                fmt->type = T_INFO;
                fprintf(stderr,"Warning: Assuming INFO/%s\n", key);
            }
        }
    }

    switch (fmt->type)
    {
        case T_FIRST_ALT: fmt->handler = &process_first_alt; break;
        case T_CHROM_POS_ID: fmt->handler = &process_chrom_pos_id; break;
        case T_GT_TO_PROB3: fmt->handler = &process_gt_to_prob3; break;
        case T_PL_TO_PROB3: fmt->handler = &process_pl_to_prob3; break;
        case T_CHROM: fmt->handler = &process_chrom; break;
        case T_POS: fmt->handler = &process_pos; break;
        case T_ID: fmt->handler = &process_id; break;
        case T_REF: fmt->handler = &process_ref; break;
        case T_ALT: fmt->handler = &process_alt; break;
        case T_QUAL: fmt->handler = &process_qual; break;
        case T_FILTER: fmt->handler = &process_filter; convert->max_unpack |= BCF_UN_FLT; break;
        case T_INFO: fmt->handler = &process_info; convert->max_unpack |= BCF_UN_INFO; break;
        case T_FORMAT: fmt->handler = &process_format; convert->max_unpack |= BCF_UN_FMT; break;
        case T_SAMPLE: fmt->handler = &process_sample; break;
        case T_SEP: fmt->handler = &process_sep; break;
        case T_IS_TS: fmt->handler = &process_is_ts; break;
        case T_TYPE: fmt->handler = &process_type; break;
        case T_MASK: fmt->handler = NULL; break;
        case T_GT: fmt->handler = &process_gt; convert->max_unpack |= BCF_UN_FMT; break;
        case T_TGT: fmt->handler = &process_tgt; convert->max_unpack |= BCF_UN_FMT; break;
        case T_LINE: fmt->handler = &process_line; break;
        default: error("TODO: handler for type %d\n", fmt->type);
    }
    if ( key )
    {
        if ( fmt->type==T_INFO )
        {
            fmt->id = bcf_hdr_id2int(convert->header, BCF_DT_ID, key);
            if ( fmt->id==-1 ) error("Error: no such tag defined in the VCF header: INFO/%s\n", key);
        }
    }
    return fmt;
}

static int parse_subscript(char **p)
{
    char *q = *p;
    if ( *q!='{' ) return -1;
    q++;
    while ( *q && *q!='}' && isdigit(*q) ) q++;
    if ( *q!='}' ) return -1;
    int idx = atoi((*p)+1);
    *p = q+1;
    return idx;
}

static char *parse_tag(convert_t *convert, char *p, int is_gtf)
{
    char *q = ++p;
    while ( *q && (isalnum(*q) || *q=='_') ) q++;
    kstring_t str = {0,0,0};
    if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
    kputsn(p, q-p, &str);
    if ( is_gtf )
    {
        if ( !strcmp(str.s, "SAMPLE") ) register_tag(convert, T_SAMPLE, "SAMPLE", is_gtf);
        else if ( !strcmp(str.s, "GT") ) register_tag(convert, T_GT, "GT", is_gtf);
        else if ( !strcmp(str.s, "TGT") ) register_tag(convert, T_TGT, "GT", is_gtf);
        else
        {
            fmt_t *fmt = register_tag(convert, T_FORMAT, str.s, is_gtf);
            fmt->subscript = parse_subscript(&q);
        }
    }
    else
    {
        if ( !strcmp(str.s, "CHROM") ) register_tag(convert, T_CHROM, str.s, is_gtf);
        else if ( !strcmp(str.s, "POS") ) register_tag(convert, T_POS, str.s, is_gtf);
        else if ( !strcmp(str.s, "ID") ) register_tag(convert, T_ID, str.s, is_gtf);
        else if ( !strcmp(str.s, "REF") ) register_tag(convert, T_REF, str.s, is_gtf);
        else if ( !strcmp(str.s, "ALT") ) register_tag(convert, T_ALT, str.s, is_gtf);
        else if ( !strcmp(str.s, "FIRST_ALT") ) register_tag(convert, T_FIRST_ALT, str.s, is_gtf);
        else if ( !strcmp(str.s, "QUAL") ) register_tag(convert, T_QUAL, str.s, is_gtf);
        else if ( !strcmp(str.s, "FILTER") ) register_tag(convert, T_FILTER, str.s, is_gtf);
        else if ( !strcmp(str.s, "QUAL") ) register_tag(convert, T_QUAL, str.s, is_gtf);
        else if ( !strcmp(str.s, "IS_TS") ) register_tag(convert, T_IS_TS, str.s, is_gtf);
        else if ( !strcmp(str.s, "TYPE") ) register_tag(convert, T_TYPE, str.s, is_gtf);
        else if ( !strcmp(str.s, "MASK") ) register_tag(convert, T_MASK, str.s, is_gtf);
        else if ( !strcmp(str.s, "LINE") ) register_tag(convert, T_LINE, str.s, is_gtf);
        else if ( !strcmp(str.s, "_CHROM_POS_ID") ) register_tag(convert, T_CHROM_POS_ID, str.s, is_gtf);
        else if ( !strcmp(str.s, "_GT_TO_PROB3") ) register_tag(convert, T_GT_TO_PROB3, str.s, is_gtf);
        else if ( !strcmp(str.s, "_PL_TO_PROB3") ) register_tag(convert, T_PL_TO_PROB3, str.s, is_gtf);
        else if ( !strcmp(str.s, "INFO") )
        {
            if ( *q!='/' ) error("Could not parse format string: %s\n", convert->format_str);
            p = ++q;
            str.l = 0;
            while ( *q && (isalnum(*q) || *q=='_') ) q++;
            if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
            kputsn(p, q-p, &str);
            fmt_t *fmt = register_tag(convert, T_INFO, str.s, is_gtf);
            fmt->subscript = parse_subscript(&q);
        }
        else
        {
            fmt_t *fmt = register_tag(convert, T_INFO, str.s, is_gtf);
            fmt->subscript = parse_subscript(&q);
        }
    }
    free(str.s);
    return q;
}

static char *parse_sep(convert_t *convert, char *p, int is_gtf)
{
    char *q = p;
    kstring_t str = {0,0,0};
    while ( *q && *q!='[' && *q!=']' && *q!='%' )
    {
        if ( *q=='\\' )
        {
            q++;
            if ( *q=='n' ) kputc('\n', &str);
            else if ( *q=='t' ) kputc('\t', &str);
            else kputc(*q, &str);
        }
        else kputc(*q, &str);
        q++;
    }
    if ( !str.l ) error("Could not parse format string: %s\n", convert->format_str);
    register_tag(convert, T_SEP, str.s, is_gtf);
    free(str.s);
    return q;
}

convert_t *convert_init(bcf_hdr_t *hdr, int *samples, int nsamples, const char *format_str)
{
    convert_t *convert = (convert_t*) calloc(1,sizeof(convert_t));
    convert->header = hdr;
    convert->format_str = strdup(format_str);

    int i, is_gtf = 0;
    char *p = convert->format_str;
    while ( *p )
    {
        //fprintf(stderr,"<%s>\n", p);
        switch (*p)
        {
            case '[': is_gtf = 1; p++; break;
            case ']': is_gtf = 0; register_tag(convert, T_SEP, NULL, 0); p++; break;
            case '%': p = parse_tag(convert, p, is_gtf); break;
            default:  p = parse_sep(convert, p, is_gtf); break;
        }
    }

    if ( nsamples )
    {
        convert->nsamples = nsamples;
        convert->samples = (int*) malloc(sizeof(int)*nsamples);
        for (i=0; i<convert->nsamples; i++) convert->samples[i] = samples[i];
    }
    else
    {
        convert->nsamples = bcf_hdr_nsamples(convert->header);
        convert->samples = (int*) malloc(sizeof(int)*convert->nsamples);
        for (i=0; i<convert->nsamples; i++) convert->samples[i] = i;
    }
    return convert;
}

void convert_destroy(convert_t *convert)
{
    int i;
    for (i=0; i<convert->nfmt; i++)
        if ( convert->fmt[i].key ) free(convert->fmt[i].key);
    if ( convert->mfmt ) free(convert->fmt);
    free(convert->dat);
    free(convert->samples);
    free(convert->format_str);
    free(convert);
}


int convert_header(convert_t *convert, kstring_t *str)
{
    int i, icol = 0, l_ori = str->l;
    bcf_hdr_t *hdr = convert->header;

    // Supress the header output if LINE is present
    for (i=0; i<convert->nfmt; i++)
        if ( convert->fmt[i].type == T_LINE ) break;
    if ( i!=convert->nfmt )
        return str->l - l_ori;

    kputs("# ", str);
    for (i=0; i<convert->nfmt; i++)
    {
        // Genotype fields
        if ( convert->fmt[i].is_gt_field )
        {
            int j = i, js, k;
            while ( convert->fmt[j].is_gt_field ) j++;
            for (js=0; js<convert->nsamples; js++)
            {
                int ks = convert->samples[js];
                for (k=i; k<j; k++)
                {
                    if ( convert->fmt[k].type == T_SEP )
                    {
                        if ( convert->fmt[k].key ) kputs(convert->fmt[k].key, str);
                    }
                    else if ( convert->fmt[k].type == T_SAMPLE )
                        ksprintf(str, "[%d]%s", ++icol, convert->fmt[k].key);
                    else
                        ksprintf(str, "[%d]%s:%s", ++icol, hdr->samples[ks], convert->fmt[k].key);
                }
            }
            i = j-1;
            continue;
        }
        // Fixed fields
        if ( convert->fmt[i].type == T_SEP )
        {
            if ( convert->fmt[i].key ) kputs(convert->fmt[i].key, str);
            continue;
        }
        ksprintf(str, "[%d]%s", ++icol, convert->fmt[i].key);
    }
    return str->l - l_ori;
}

int convert_line(convert_t *convert, bcf1_t *line, kstring_t *str)
{
    int l_ori = str->l;
    bcf_unpack(line, convert->max_unpack);

    int i, ir;
    str->l = 0;
    for (i=0; i<convert->nfmt; i++)
    {
        // Genotype fields
        if ( convert->fmt[i].is_gt_field )
        {
            int j = i, js, k;
            while ( convert->fmt[j].is_gt_field )
            {
                convert->fmt[j].ready = 0;
                j++;
            }
            for (js=0; js<convert->nsamples; js++)
            {
                int ks = convert->samples[js];
                for (k=i; k<j; k++)
                {
                    if ( convert->fmt[k].type == T_MASK )
                    {
                        for (ir=0; ir<convert->nreaders; ir++)
                            kputc(bcf_sr_has_line(convert->readers,ir)?'1':'0', str);
                    }
                    else if ( convert->fmt[k].handler )
                        convert->fmt[k].handler(convert, line, &convert->fmt[k], ks, str);
                }
            }
            i = j-1;
            continue;
        }
        // Fixed fields
        if ( convert->fmt[i].type == T_MASK )
        {
            for (ir=0; ir<convert->nreaders; ir++)
                kputc(bcf_sr_has_line(convert->readers,ir)?'1':'0', str);
        }
        else if ( convert->fmt[i].handler )
            convert->fmt[i].handler(convert, line, &convert->fmt[i], -1, str);
    }
    return str->l - l_ori;
}

