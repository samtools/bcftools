/*  vcfquery.c -- Extracts fields from VCF/BCF file.

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
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"

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

struct _args_t;
typedef struct _args_t args_t;

typedef struct _fmt_t
{
    int type, id, is_gt_field, ready, subscript;
    int nsamples, *samples;
    char *key;
    bcf_fmt_t *fmt;
    void (*handler)(args_t *, bcf1_t *, struct _fmt_t *, int, kstring_t *);
}
fmt_t;

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

struct _args_t
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    fmt_t *fmt;
    int nfmt, mfmt;
    bcf_srs_t *files;
    bcf_hdr_t *header;
    int nsamples, *samples, sample_is_file;
    char **argv, *format, *sample_list, *targets_list, *regions_list, *vcf_list;
    int argc, list_columns, print_header;
};

/**
 *  ks_getline() - Read next line from $fp, appending it to $str.  The newline
 *  is stripped and \0 appended. Returns the number of characters read
 *  excluding the null byte.
 */
size_t ks_getline(FILE *fp, kstring_t *str)
{
    size_t nread=0;
    int c;
    while ((c=getc(fp))!= EOF && c!='\n')
    {
        nread++;
        if ( str->l+nread > str->m )
        {
            str->m += 1024;
            str->s = (char*) realloc(str->s, sizeof(char)*str->m);
        }
        str->s[str->l+nread-1] = c;
    }
    if ( str->l >= str->m )
    {
        str->m += 1024;
        str->s = (char*) realloc(str->s, sizeof(char)*str->m);
    }
    str->l += nread;
    str->s[ str->l ] = 0;
    return nread;
}
char **read_list(char *fname, int *_n)
{
    int n = 0;
    char **list = NULL;

    FILE *fp = fopen(fname,"r");
    if ( !fp ) error("%s: %s\n", fname, strerror(errno));

    kstring_t str = {0,0,0};
    while ( ks_getline(fp, &str) )
    {
        list = (char**) realloc(list, sizeof(char*)*(++n));
        list[n-1] = strdup(str.s);
        str.l = 0;
    }
    fclose(fp);
    if ( str.m ) free(str.s);
    *_n = n;
    return list;
}
void destroy_list(char **list, int n)
{
    int i;
    for (i=0; i<n; i++)
        free(list[i]);
    free(list);
}

static void process_chrom(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(args->header->id[BCF_DT_CTG][line->rid].key, str); }
static void process_pos(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputw(line->pos+1, str); }
static void process_id(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(line->d.id, str); }
static void process_ref(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(line->d.allele[0], str); }
static void process_alt(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
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
static void process_qual(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( bcf_float_is_missing(line->qual) ) kputc('.', str);
    else ksprintf(str, "%g", line->qual);
}
static void process_filter(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i;
    if ( line->d.n_flt )
    {
        for (i=0; i<line->d.n_flt; i++)
        {
            if (i) kputc(';', str);
            kputs(args->header->id[BCF_DT_ID][line->d.flt[i]].key, str);
        }
    }
    else kputc('.', str);
}
static int bcf_array_ivalue(void *bcf_array, int type, int idx)
{
    if ( type==BCF_BT_INT8 ) return ((int8_t*)bcf_array)[idx];
    if ( type==BCF_BT_INT16 ) return ((int16_t*)bcf_array)[idx];
    return ((int32_t*)bcf_array)[idx];
}
static void process_info(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
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
static void init_format(args_t *args, bcf1_t *line, fmt_t *fmt)
{
    fmt->id = bcf_hdr_id2int(args->header, BCF_DT_ID, fmt->key);
    if ( fmt->id==-1 ) error("Error: no such tag defined in the VCF header: FORMAT/%s\n", fmt->key);
    fmt->fmt = NULL;
    int i;
    for (i=0; i<(int)line->n_fmt; i++)
        if ( line->d.fmt[i].id==fmt->id ) { fmt->fmt = &line->d.fmt[i]; break; }
    fmt->ready = 1;
}
static void process_format(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(args, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }
    bcf_fmt_array(str, fmt->fmt->n, fmt->fmt->type, fmt->fmt->p + isample*fmt->fmt->size);
}
static void process_gt(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(args, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }
    bcf_format_gt(fmt->fmt, isample, str);
}
static void process_tgt(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(args, line, fmt);

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
static void process_sample(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    kputs(args->header->samples[isample], str);
}
static void process_sep(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { if (fmt->key) kputs(fmt->key, str); }
static void process_is_ts(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int is_ts = 0;
    if ( bcf_get_variant_types(line) & (VCF_SNP|VCF_MNP) )
        is_ts = abs(bcf_acgt2int(*line->d.allele[0])-bcf_acgt2int(*line->d.allele[1])) == 2 ? 1 : 0;
    kputc(is_ts ? '1' : '0', str);
}
static void process_type(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int line_type = bcf_get_variant_types(line);
    int i = 0;
    if ( line_type == VCF_REF ) { kputs("REF", str); i++; }
    if ( line_type & VCF_SNP ) { if (i) kputc(',',str); kputs("SNP", str); i++; }
    if ( line_type & VCF_MNP ) { if (i) kputc(',',str); kputs("MNP", str); i++; }
    if ( line_type & VCF_INDEL ) { if (i) kputc(',',str); kputs("INDEL", str); i++; }
    if ( line_type & VCF_OTHER ) { if (i) kputc(',',str); kputs("OTHER", str); i++; }
}
static void process_line(args_t *args, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    vcf_format1(args->header, line, str);
}

static fmt_t *register_tag(args_t *args, int type, char *key, int is_gtf)
{
    args->nfmt++;
    if ( args->nfmt > args->mfmt )
    {
        args->mfmt += 10;
        args->fmt   = (fmt_t*) realloc(args->fmt, args->mfmt*sizeof(fmt_t));

    }
    fmt_t *fmt = &args->fmt[ args->nfmt-1 ];
    fmt->type  = type;
    fmt->key   = key ? strdup(key) : NULL;
    fmt->is_gt_field = is_gtf;
    fmt->subscript = -1;

    // Allow non-format tags, such as CHROM, INFO, etc., to appear amongst the format tags.
    if ( key )
    {
        int id = bcf_hdr_id2int(args->header, BCF_DT_ID, key);
        if ( fmt->type==T_FORMAT && !bcf_hdr_idinfo_exists(args->header,BCF_HL_FMT,id) )
        {
            if ( !strcmp("CHROM",key) ) { fmt->type = T_CHROM; }
            else if ( !strcmp("POS",key) ) { fmt->type = T_POS; }
            else if ( !strcmp("ID",key) ) { fmt->type = T_ID; }
            else if ( !strcmp("REF",key) ) { fmt->type = T_REF; }
            else if ( !strcmp("ALT",key) ) { fmt->type = T_ALT; }
            else if ( !strcmp("QUAL",key) ) { fmt->type = T_QUAL; }
            else if ( !strcmp("FILTER",key) ) { fmt->type = T_FILTER; }
            else if ( id>=0 && bcf_hdr_idinfo_exists(args->header,BCF_HL_INFO,id) )
            {
                fmt->type = T_INFO;
                fprintf(stderr,"Warning: Assuming INFO/%s\n", key);
            }
        }
    }

    switch (fmt->type)
    {
        case T_CHROM: fmt->handler = &process_chrom; break;
        case T_POS: fmt->handler = &process_pos; break;
        case T_ID: fmt->handler = &process_id; break;
        case T_REF: fmt->handler = &process_ref; break;
        case T_ALT: fmt->handler = &process_alt; break;
        case T_QUAL: fmt->handler = &process_qual; break;
        case T_FILTER: fmt->handler = &process_filter; args->files->max_unpack |= BCF_UN_FLT; break;
        case T_INFO: fmt->handler = &process_info; args->files->max_unpack |= BCF_UN_INFO; break;
        case T_FORMAT: fmt->handler = &process_format; args->files->max_unpack |= BCF_UN_FMT; break;
        case T_SAMPLE: fmt->handler = &process_sample; break;
        case T_SEP: fmt->handler = &process_sep; break;
        case T_IS_TS: fmt->handler = &process_is_ts; break;
        case T_TYPE: fmt->handler = &process_type; break;
        case T_MASK: fmt->handler = NULL; break;
        case T_GT: fmt->handler = &process_gt; args->files->max_unpack |= BCF_UN_FMT; break;
        case T_TGT: fmt->handler = &process_tgt; args->files->max_unpack |= BCF_UN_FMT; break;
        case T_LINE: fmt->handler = &process_line; break;
        default: error("TODO: handler for type %d\n", fmt->type);
    }
    if ( key )
    {
        if ( fmt->type==T_INFO )
        {
            fmt->id = bcf_hdr_id2int(args->header, BCF_DT_ID, key);
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

static char *parse_tag(args_t *args, char *p, int is_gtf)
{
    char *q = ++p;
    while ( *q && (isalnum(*q) || *q=='_') ) q++;
    kstring_t str = {0,0,0};
    if ( q-p==0 ) error("Could not parse format string: %s\n", args->format);
    kputsn(p, q-p, &str);
    if ( is_gtf )
    {
        if ( !strcmp(str.s, "SAMPLE") ) register_tag(args, T_SAMPLE, "SAMPLE", is_gtf);
        else if ( !strcmp(str.s, "GT") ) register_tag(args, T_GT, "GT", is_gtf);
        else if ( !strcmp(str.s, "TGT") ) register_tag(args, T_TGT, "GT", is_gtf);
        else
        {
            fmt_t *fmt = register_tag(args, T_FORMAT, str.s, is_gtf);
            fmt->subscript = parse_subscript(&q);
        }
    }
    else
    {
        if ( !strcmp(str.s, "CHROM") ) register_tag(args, T_CHROM, str.s, is_gtf);
        else if ( !strcmp(str.s, "POS") ) register_tag(args, T_POS, str.s, is_gtf);
        else if ( !strcmp(str.s, "ID") ) register_tag(args, T_ID, str.s, is_gtf);
        else if ( !strcmp(str.s, "REF") ) register_tag(args, T_REF, str.s, is_gtf);
        else if ( !strcmp(str.s, "ALT") ) register_tag(args, T_ALT, str.s, is_gtf);
        else if ( !strcmp(str.s, "QUAL") ) register_tag(args, T_QUAL, str.s, is_gtf);
        else if ( !strcmp(str.s, "FILTER") ) register_tag(args, T_FILTER, str.s, is_gtf);
        else if ( !strcmp(str.s, "QUAL") ) register_tag(args, T_QUAL, str.s, is_gtf);
        else if ( !strcmp(str.s, "IS_TS") ) register_tag(args, T_IS_TS, str.s, is_gtf);
        else if ( !strcmp(str.s, "TYPE") ) register_tag(args, T_TYPE, str.s, is_gtf);
        else if ( !strcmp(str.s, "MASK") ) register_tag(args, T_MASK, str.s, is_gtf);
        else if ( !strcmp(str.s, "LINE") ) register_tag(args, T_LINE, str.s, is_gtf);
        else if ( !strcmp(str.s, "INFO") )
        {
            if ( *q!='/' ) error("Could not parse format string: %s\n", args->format);
            p = ++q;
            str.l = 0;
            while ( *q && (isalnum(*q) || *q=='_') ) q++;
            if ( q-p==0 ) error("Could not parse format string: %s\n", args->format);
            kputsn(p, q-p, &str);
            fmt_t *fmt = register_tag(args, T_INFO, str.s, is_gtf);
            fmt->subscript = parse_subscript(&q);
        }
        else
        {
            fmt_t *fmt = register_tag(args, T_INFO, str.s, is_gtf);
            fmt->subscript = parse_subscript(&q);
        }
    }
    free(str.s);
    return q;
}

static char *parse_sep(args_t *args, char *p, int is_gtf)
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
    if ( !str.l ) error("Could not parse format string: %s\n", args->format);
    register_tag(args, T_SEP, str.s, is_gtf);
    free(str.s);
    return q;
}

static void init_data(args_t *args)
{
    args->header = args->files->readers[0].header;
    int is_gtf = 0;
    char *p = args->format;
    while ( *p )
    {
        //fprintf(stderr,"<%s>\n", p);
        switch (*p)
        {
            case '[': is_gtf = 1; p++; break;
            case ']': is_gtf = 0; register_tag(args, T_SEP, NULL, 0); p++; break;
            case '%': p = parse_tag(args, p, is_gtf); break;
            default:  p = parse_sep(args, p, is_gtf); break;
        }
    }
    int i;
    if ( args->sample_list && strcmp("-",args->sample_list) )
    {
        for (i=0; i<args->files->nreaders; i++)
        {
            int ret = bcf_hdr_set_samples(args->files->readers[i].header,args->sample_list,args->sample_is_file);
            if ( ret<0 ) error("Error parsing the sample list\n");
            else if ( ret>0 ) error("Sample name mismatch: sample #%d not found in the header\n", ret);
        }

        if ( args->sample_list[0]!='^' )
        {
            // the sample ordering may be different if not negated
            int n;
            char **smpls = hts_readlist(args->sample_list, args->sample_is_file, &n);
            if ( !smpls ) error("Could not parse %s\n", args->sample_list);
            if ( n!=bcf_hdr_nsamples(args->files->readers[0].header) )
                error("The number of samples does not match, perhaps some are present multiple times?\n");
            args->nsamples = bcf_hdr_nsamples(args->files->readers[0].header);
            args->samples = (int*) malloc(sizeof(int)*args->nsamples);
            for (i=0; i<n; i++)
            {
                args->samples[i] = bcf_hdr_id2int(args->files->readers[0].header, BCF_DT_SAMPLE,smpls[i]);
                free(smpls[i]);
            }
            free(smpls);
        }
        else
        {
            args->nsamples = bcf_hdr_nsamples(args->files->readers[0].header);
            args->samples = (int*) malloc(sizeof(int)*args->nsamples);
            for (i=0; i<args->nsamples; i++) args->samples[i] = i;
        }
    }
    else
    {
        args->nsamples = bcf_hdr_nsamples(args->files->readers[0].header);
        args->samples = (int*) malloc(sizeof(int)*args->nsamples);
        for (i=0; i<args->nsamples; i++) args->samples[i] = i;
    }

    if ( args->filter_str )
        args->filter = filter_init(args->header, args->filter_str);
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nfmt; i++)
        if ( args->fmt[i].key ) free(args->fmt[i].key);
    if ( args->mfmt ) free(args->fmt);
    args->nfmt = args->mfmt = 0;
    args->fmt = NULL;
    if ( args->filter )
        filter_destroy(args->filter);
    free(args->samples);
}


static void print_header(args_t *args, kstring_t *str)
{
    int i, icol = 0;
    bcf_hdr_t *hdr = args->files->readers[0].header;

    // Supress the header output if LINE is present
    for (i=0; i<args->nfmt; i++)
        if ( args->fmt[i].type == T_LINE ) break;
    if ( i!=args->nfmt )
        return;

    kputs("# ", str);
    for (i=0; i<args->nfmt; i++)
    {
        // Genotype fields
        if ( args->fmt[i].is_gt_field )
        {
            int j = i, js, k;
            while ( args->fmt[j].is_gt_field ) j++;
            for (js=0; js<args->nsamples; js++)
            {
                int ks = args->samples[js];
                for (k=i; k<j; k++)
                {
                    if ( args->fmt[k].type == T_SEP )
                    {
                        if ( args->fmt[k].key ) kputs(args->fmt[k].key, str);
                    }
                    else if ( args->fmt[k].type == T_SAMPLE )
                        ksprintf(str, "[%d]%s", ++icol, args->fmt[k].key);
                    else
                        ksprintf(str, "[%d]%s:%s", ++icol, hdr->samples[ks], args->fmt[k].key);
                }
            }
            i = j-1;
            continue;
        }
        // Fixed fields
        if ( args->fmt[i].type == T_SEP )
        {
            if ( args->fmt[i].key ) kputs(args->fmt[i].key, str);
            continue;
        }
        ksprintf(str, "[%d]%s", ++icol, args->fmt[i].key);
    }
    fwrite(str->s, str->l, 1, stdout);
    str->l = 0;
}

static void query_vcf(args_t *args)
{
    kstring_t str = {0,0,0};

    args->files->max_unpack |= BCF_UN_STR;
    if ( args->print_header ) print_header(args, &str);

    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files,0) ) continue;
        bcf1_t *line = args->files->readers[0].buffer[0];
        bcf_unpack(line, args->files->max_unpack);

        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }

        int i, ir;
        str.l = 0;
        for (i=0; i<args->nfmt; i++)
        {
            // Genotype fields
            if ( args->fmt[i].is_gt_field )
            {
                int j = i, js, k;
                while ( args->fmt[j].is_gt_field )
                {
                    args->fmt[j].ready = 0;
                    j++;
                }
                for (js=0; js<args->nsamples; js++)
                {
                    int ks = args->samples[js];
                    for (k=i; k<j; k++)
                    {
                        if ( args->fmt[k].type == T_MASK )
                        {
                            for (ir=0; ir<args->files->nreaders; ir++)
                                kputc(bcf_sr_has_line(args->files,ir)?'1':'0', &str);
                        }
                        else if ( args->fmt[k].handler )
                            args->fmt[k].handler(args, line, &args->fmt[k], ks, &str);
                    }
                }
                i = j-1;
                continue;
            }
            // Fixed fields
            if ( args->fmt[i].type == T_MASK )
            {
                for (ir=0; ir<args->files->nreaders; ir++)
                    kputc(bcf_sr_has_line(args->files,ir)?'1':'0', &str);
            }
            else if ( args->fmt[i].handler )
                args->fmt[i].handler(args, line, &args->fmt[i], -1, &str);
        }
        if ( str.l )
            fwrite(str.s, str.l, 1, stdout);
    }
    if ( str.m ) free(str.s);
}

static void list_columns(args_t *args)
{
    int i;
    bcf_sr_t *reader = &args->files->readers[0];
    for (i=0; i<bcf_hdr_nsamples(reader->header); i++)
        printf("%s\n", reader->header->samples[i]);
}

static char **copy_header(bcf_hdr_t *hdr, char **src, int nsrc)
{
    char **dst = (char**) malloc(sizeof(char*)*nsrc);
    int i;
    for (i=0; i<nsrc; i++) dst[i] = strdup(src[i]);
    return dst;
}
static int compare_header(bcf_hdr_t *hdr, char **a, int na, char **b, int nb)
{
    if ( na!=nb ) return na-nb;
    int i;
    for (i=0; i<na; i++)
        if ( strcmp(a[i],b[i]) ) return 1;
    return 0;
}


static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Extracts fields from VCF/BCF file and prints them in user-defined format\n");
    fprintf(stderr, "Usage:   bcftools query [options] <A.vcf.gz> [<B.vcf.gz> [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -a, --annots <list>               alias for -f '%%CHROM\\t%%POS\\t%%MASK\\t%%REF\\t%%ALT\\t%%TYPE\\t' + tab-separated <list> of tags\n");
    fprintf(stderr, "    -c, --collapse <string>           collapse lines with duplicate positions for <snps|indels|both|all|some|none>, see man page [none]\n");
    fprintf(stderr, "    -e, --exclude <expr>              exclude sites for which the expression is true (see below for details)\n");
    fprintf(stderr, "    -f, --format <string>             learn by example, see below\n");
    fprintf(stderr, "    -H, --print-header                print header\n");
    fprintf(stderr, "    -i, --include <expr>              select sites for which the expression is true (see below for details)\n");
    fprintf(stderr, "    -l, --list-samples                print the list of samples and exit\n");
    fprintf(stderr, "    -r, --regions <region>            restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>         restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --samples <list>              list of samples to include\n");
    fprintf(stderr, "    -S, --samples-file <file>         file of samples to include\n");
    fprintf(stderr, "    -t, --targets <region>            similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>         similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "    -v, --vcf-list <file>             process multiple VCFs listed in the file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Format expressions:\n");
    fprintf(stderr, "\t%%CHROM          The CHROM column (similarly also other columns, such as POS, ID, QUAL, etc.)\n");
    fprintf(stderr, "\t%%INFO/TAG       Any tag in the INFO column\n");
    fprintf(stderr, "\t%%TYPE           Variant type (REF, SNP, MNP, INDEL, OTHER)\n");
    fprintf(stderr, "\t%%MASK           Indicates presence of the site in other files (with multiple files)\n");
    fprintf(stderr, "\t%%TAG{INT}       Curly brackets to subscript vectors (0-based)\n");
    fprintf(stderr, "\t[]              The brackets loop over all samples\n");
    fprintf(stderr, "\t%%GT             Genotype (e.g. 0/1)\n");
    fprintf(stderr, "\t%%TGT            Translated genotype (e.g. C/A)\n");
    fprintf(stderr, "\t%%LINE           Prints the whole line\n");
    fprintf(stderr, "\t%%SAMPLE         Sample name\n");
    //fprintf(stderr, "\t%*<A><B>        All format fields printed as KEY<A>VALUE<B>\n");
    fprintf(stderr, "\n");
    filter_expression_info(stderr);
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "\tbcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE=%%GT]\\n' file.vcf.gz\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfquery(int argc, char *argv[])
{
    int c, collapse = 0;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    int regions_is_file = 0, targets_is_file = 0;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"list-samples",0,0,'l'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"format",1,0,'f'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"annots",1,0,'a'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"print-header",0,0,'H'},
        {"collapse",1,0,'c'},
        {"vcf-list",1,0,'v'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "hlr:R:f:a:s:S:Ht:T:c:v:i:e:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'f': args->format = strdup(optarg); break;
            case 'H': args->print_header = 1; break;
            case 'v': args->vcf_list = optarg; break;
            case 'c':
                if ( !strcmp(optarg,"snps") ) collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"any") ) collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"all") ) collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"some") ) args->files->collapse |= COLLAPSE_SOME;
                else error("The --collapse string \"%s\" not recognised.\n", optarg);
                break;
            case 'a':
                {
                    kstring_t str = {0,0,0};
                    kputs("%CHROM\t%POS\t%MASK\t%REF\t%ALT\t%", &str);
                    char *p = optarg;
                    while ( *p )
                    {
                        if ( *p==',' )
                            kputs("\t%", &str);
                        else
                            kputc(*p, &str);
                        p++;
                    }
                    kputc('\n', &str);
                    args->format = str.s;
                    break;
                }
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'l': args->list_columns = 1; break;
            case 's': args->sample_list = optarg; break;
            case 'S': args->sample_list = optarg; args->sample_is_file = 1; break;
            case 'h':
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";
    }
    else fname = argv[optind];

    if ( args->list_columns )
    {
        if ( !fname ) error("Missing the VCF file name\n");
        args->files = bcf_sr_init();
        if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open or the file not indexed: %s\n", fname);
        list_columns(args);
        bcf_sr_destroy(args->files);
        free(args);
        return 0;
    }

    if ( !args->format ) usage();
    if ( !args->vcf_list )
    {
        if ( !fname ) usage();
        args->files = bcf_sr_init();
        args->files->collapse = collapse;
        if ( optind+1 < argc ) args->files->require_index = 1;
        if ( args->regions_list && bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
        if ( args->targets_list )
        {
            if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
                error("Failed to read the targets: %s\n", args->targets_list);
        }
        while ( fname )
        {
            if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open or the file not indexed: %s\n", fname);
            fname = ++optind < argc ? argv[optind] : NULL;
        }
        init_data(args);
        query_vcf(args);
        free(args->format);
        destroy_data(args);
        bcf_sr_destroy(args->files);
        free(args);
        return 0;
    }

    // multiple VCFs
    int i, k, nfiles, prev_nsamples = 0;
    char **fnames, **prev_samples = NULL;
    fnames = read_list(args->vcf_list, &nfiles);
    if ( !nfiles ) error("No files in %s?\n", args->vcf_list);
    for (i=0; i<nfiles; i++)
    {
        args->files = bcf_sr_init();
        args->files->collapse = collapse;
        if ( args->regions_list && bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
        if ( optind < argc ) args->files->require_index = 1;
        if ( args->targets_list )
        {
            if ( bcf_sr_set_targets(args->files, args->targets_list,targets_is_file, 0)<0 )
                error("Failed to read the targets: %s\n", args->targets_list);
        }
        if ( !bcf_sr_add_reader(args->files, fnames[i]) ) error("Failed to open or the file not indexed: %s\n", fnames[i]);
        for (k=optind; k<argc; k++)
            if ( !bcf_sr_add_reader(args->files, argv[k]) ) error("Failed to open or the file not indexed: %s\n", argv[k]);
        init_data(args);
        if ( i==0 )
            prev_samples = copy_header(args->header, args->files->readers[0].header->samples, bcf_hdr_nsamples(args->files->readers[0].header));
        else
        {
            args->print_header = 0;
            if ( compare_header(args->header, args->files->readers[0].header->samples, bcf_hdr_nsamples(args->files->readers[0].header), prev_samples, prev_nsamples) )
                error("Different samples in %s and %s\n", fnames[i-1],fnames[i]);
        }
        query_vcf(args);
        destroy_data(args);
        bcf_sr_destroy(args->files);
    }
    destroy_list(fnames, nfiles);
    destroy_list(prev_samples, prev_nsamples);
    free(args->format);
    free(args);
    return 0;
}


