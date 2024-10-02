/*
    Copyright (C) 2017-2023 Genome Research Ltd.

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
/*
    Split VCF by sample(s)
*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"
#include "filter.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    char **rename;      // use a new sample name (rename samples)
    int nsmpl, *smpl;   // number of samples to keep and their indices in the input header
    htsFile *fh;        // output file handle
    char *fname;        // output file name
    filter_t *filter;
    bcf_hdr_t *hdr;
    char *index_fn;
}
subset_t;

typedef struct
{
    char *filter_str;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)
    uint8_t *info_tags, *fmt_tags;
    int ninfo_tags, minfo_tags, nfmt_tags, mfmt_tags, keep_info, keep_fmt;
    int argc, region_is_file, target_is_file, output_type, clevel;
    int regions_overlap, targets_overlap;
    char **argv, *region, *target, *fname, *output_dir, *keep_tags, *samples_fname, *groups_fname;
    void *unique_fnames;
    bcf_hdr_t *hdr_in, *hdr_out;
    bcf_srs_t *sr;
    subset_t *sets;
    int nsets, nhts_opts;
    char **hts_opts;
    int write_index;
}
args_t;

const char *about(void)
{
    return "Split VCF by sample, creating single- or multi-sample VCFs\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Split VCF by sample, creating single- or multi-sample VCFs. The output files are named\n"
        "       by sample names whenever possible, with the characters from the set [ \\t:/\\] replaced\n"
        "       with \"_\", and a unique numeric suffix added in case of name clashes.\n"
        "Usage: bcftools +split [Options]\n"
        "Plugin options:\n"
        "   -e, --exclude EXPR              Exclude sites for which the expression is true (applied on the outputs)\n"
        "   -G, --groups-file FILE          Similar to -S, but the samples are split by group:\n"
        "                                       \n"
        "                                       # Create two output files (third column) with the second sample appearing\n"
        "                                       # in both. The second column is for optional renaming of the samples, use\n"
        "                                       # dash \"-\" to keep sample names unchanged\n"
        "                                       sample1   -          file1\n"
        "                                       sample2   -          file1,file2\n"
        "                                       sample3   new-name3  file2\n"
        "                                       \n"
        "   -i, --include EXPR              Include only sites for which the expression is true (applied on the outputs)\n"
        "   -k, --keep-tags LIST            List of tags to keep. By default all tags are preserved\n"
        "   -o, --output DIR                Write output to the directory DIR\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -r, --regions REGION            Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n"
        "   -S, --samples-file FILE         List of samples to keep with up to three columns, one line per output file:\n"
        "                                       \n"
        "                                       # Create two output files, the first sample is the basename\n"
        "                                       # of the new file\n"
        "                                       sample1\n"
        "                                       sample2,sample3\n"
        "                                       \n"
        "                                       # Optional second column to rename the samples\n"
        "                                       sample1           new-name2\n"
        "                                       sample2,sample3   new-name2,new-name3\n"
        "                                       \n"
        "                                       # Optional third column to provide output file base name, use dash \"-\"\n"
        "                                       # to keep sample names unchanged\n"
        "                                       sample1           new-name1   file1\n"
        "                                       sample2,sample3   -           file2\n"
        "                                       \n"
        "   -t, --targets REGION            Similar to -r but streams rather than index-jumps\n"
        "       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n"
        "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n"
        "       --hts-opts LIST             Low-level options to pass to HTSlib, e.g. block_size=32768\n"
        "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
        "\n"
        "Examples:\n"
        "   # Split a VCF file\n"
        "   bcftools +split input.bcf -Ob -o dir\n"
        "\n"
        "   # Exclude sites with missing or hom-ref genotypes\n"
        "   bcftools +split input.bcf -Ob -o dir -i'GT=\"alt\"'\n"
        "\n"
        "   # Keep all INFO tags but only GT and PL in FORMAT\n"
        "   bcftools +split input.bcf -Ob -o dir -k INFO,FMT/GT,PL\n"
        "\n"
        "   # Keep all FORMAT tags but drop all INFO tags\n"
        "   bcftools +split input.bcf -Ob -o dir -k FMT\n"
        "\n";
}

void mkdir_p(const char *fmt, ...) HTS_FORMAT(HTS_PRINTF_FMT, 1, 2);

static char *create_unique_file_name(args_t *args, const char *template)
{
    kstring_t str = {0,0,0};
    kputs(template, &str);
    char *ptr = str.s;
    while ( *ptr )
    {
        if ( *ptr==':' || *ptr=='\\' || *ptr=='/' || *ptr==' ' || *ptr=='\t' ) *ptr = '_';
        ptr++;
    }
    size_t ori_len = str.l;
    int id = 0;
    if ( !args->unique_fnames ) args->unique_fnames = khash_str2int_init();
    while ( khash_str2int_has_key(args->unique_fnames,str.s) )
    {
        str.l = ori_len;
        kputc('-', &str);
        kputw(++id, &str);
    }
    khash_str2int_inc(args->unique_fnames, strdup(str.s));
    return str.s;
}

void init_subsets(args_t *args)
{
    int i,j, nsmpl = bcf_hdr_nsamples(args->hdr_in);
    if ( !args->samples_fname && !args->groups_fname )
    {
        args->nsets = nsmpl;
        args->sets  = (subset_t*) calloc(nsmpl, sizeof(subset_t));
        for (i=0; i<nsmpl; i++)
        {
            subset_t *set = &args->sets[i];
            set->nsmpl = 1;
            set->smpl  = (int*) calloc(1, sizeof(*set->smpl));
            set->smpl[0] = i;
            set->fname   = create_unique_file_name(args, args->hdr_in->samples[i]);
        }
    }
    else if ( args->samples_fname )
    {
        kstring_t str = {0,0,0};
        int nfiles = 0;
        char **files = hts_readlines(args->samples_fname, &nfiles);
        if ( !nfiles || !files ) error("Failed to parse %s\n", args->samples_fname);
        args->nsets = 0;
        args->sets = (subset_t*) calloc(nfiles, sizeof(subset_t));
        for (i=0,j=0; i<nfiles; i++)
        {
            subset_t *set = &args->sets[args->nsets];
            set->nsmpl = 1;
            str.l = 0;
            int escaped = 0;
            char *ptr = files[i];
            while ( *ptr )
            {
                if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
                if ( isspace(*ptr) && !escaped ) break;
                if ( *ptr==',' ) set->nsmpl++;      // todo: allow commas in sample names
                kputc(*ptr, &str);
                escaped = 0;
                ptr++;
            }
            set->smpl = (int*) calloc(set->nsmpl, sizeof(*set->smpl));
            char *beg = str.s;
            j = 0;
            while ( *beg )
            {
                char *end = beg;
                while ( *end && *end!=',' ) end++;
                char tmp = *end;
                *end = 0;
                int idx = bcf_hdr_id2int(args->hdr_in, BCF_DT_SAMPLE, beg);
                if ( idx>=0 )
                    set->smpl[j++] = idx;
                else
                    fprintf(stderr,"Warning: The sample \"%s\" is not present in %s\n", beg,args->fname);
                if ( !tmp ) break;
                beg = end + 1;
            }
            if ( !j ) continue;

            while ( *ptr && isspace(*ptr) ) ptr++;
            j = 0;
            if ( *ptr )     // optional second column with new sample names
            {
                set->rename = (char**) calloc(set->nsmpl, sizeof(*set->rename));
                beg = ptr;
                while ( *beg && !isspace(*beg) )
                {
                    ptr = beg;
                    while ( *ptr && *ptr!=',' && !isspace(*ptr) ) ptr++;
                    char tmp = *ptr;
                    *ptr = 0;
                    if ( !strcmp("-",beg) )
                    {
                        if ( j ) error("Error: cannot combine sample names with \"-\"\n");
                        *ptr = tmp;
                        break;
                    }
                    set->rename[j++] = strdup(beg);
                    *ptr = tmp;
                    if ( !tmp || isspace(tmp) ) break;
                    beg = ptr + 1;
                    if ( j >= set->nsmpl )
                        error("Expected the same number of samples in the first and second column: %s\n",files[i]);
                }
                if ( j )
                    set->fname = create_unique_file_name(args, set->rename[0]);
                else
                {
                    free(set->rename);
                    set->rename = NULL;
                }
            }

            while ( *ptr && isspace(*ptr) ) ptr++;
            if ( *ptr )     // optional third column with file name
            {
                free(set->fname);
                set->fname = create_unique_file_name(args, ptr);
            }

            if ( !set->fname )
                set->fname = create_unique_file_name(args, args->hdr_in->samples[set->smpl[0]]);

            args->nsets++;
        }
        for (i=0; i<nfiles; i++) free(files[i]);
        free(files);
        free(str.s);
    }
    else // -G, args->groups_fname is set
    {
        void *fname2set = khash_str2int_init();
        kstring_t str = {0,0,0};
        int nsamples = 0;
        char **samples = hts_readlines(args->groups_fname, &nsamples);
        if ( !nsamples || !samples ) error("Failed to parse %s\n", args->groups_fname);
        args->nsets = 0;
        for (i=0; i<nsamples; i++)
        {
            char *ptr = samples[i];
            int escaped = 0, iset;
            while ( *ptr )
            {
                if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
                if ( isspace(*ptr) && !escaped ) break;
                escaped = 0;
                ptr++;
            }
            char tmp = *ptr;
            *ptr = 0;

            int idx = bcf_hdr_id2int(args->hdr_in, BCF_DT_SAMPLE, samples[i]);
            if ( idx < 0  )
            {
                fprintf(stderr,"Warning: The sample \"%s\" is not present in %s\n", samples[i],args->fname);
                continue;
            }

            char *rename = NULL;
            if ( tmp )      // two columns: new sample name
            {
                rename = ptr + 1;
                while ( *rename && isspace(*rename) ) rename++;
                if ( !*rename ) rename = NULL;  // trailing space
                else
                {
                    ptr = rename;
                    while ( *ptr && !isspace(*ptr) ) ptr++;
                    tmp = *ptr;
                    *ptr = 0;
                    if ( !strcmp("-",rename) ) rename = NULL;
                    if ( tmp ) ptr++;
                }
                while ( *ptr && isspace(*ptr) ) ptr++;
            }

            if ( !*ptr )    // no third column, use sample name as file name
                ptr = samples[i];

            char *beg = ptr;
            while ( *beg )
            {
                ptr = beg;
                while ( *ptr && *ptr!=',' ) ptr++;
                tmp = *ptr;
                *ptr = 0;
                if ( khash_str2int_get(fname2set,beg,&iset) == -1 )
                {
                    args->nsets++;
                    args->sets = (subset_t*) realloc(args->sets,args->nsets*sizeof(subset_t));
                    iset = args->nsets - 1;
                    if ( khash_str2int_set(fname2set,strdup(beg),iset) < 0 ) error("Failed to insert the key [%s]\n",beg);
                    memset(args->sets+iset,0,sizeof(*args->sets));
                }
                subset_t *set = &args->sets[iset];
                set->nsmpl++;
                set->smpl   = (int*) realloc(set->smpl,set->nsmpl*sizeof(*set->smpl));
                set->rename = (char**) realloc(set->rename,set->nsmpl*sizeof(*set->rename));
                set->smpl[set->nsmpl-1]   = idx;
                set->rename[set->nsmpl-1] = strdup(rename?rename:samples[i]);
                if ( !set->fname) set->fname = create_unique_file_name(args, beg);
                if ( !tmp ) break;
                beg = ptr + 1;
            }
            free(samples[i]);
        }
        khash_str2int_destroy_free(fname2set);
        free(str.s);
        free(samples);
    }
    if ( args->unique_fnames ) khash_str2int_destroy_free(args->unique_fnames);
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->region )
    {
        args->sr->require_index = 1;
        bcf_sr_set_opt(args->sr,BCF_SR_REGIONS_OVERLAP,args->regions_overlap);
        if ( bcf_sr_set_regions(args->sr, args->region, args->region_is_file)<0 ) error("Failed to read the regions: %s\n",args->region);
    }
    if ( args->target )
    {
        bcf_sr_set_opt(args->sr,BCF_SR_TARGETS_OVERLAP,args->targets_overlap);
        if ( bcf_sr_set_targets(args->sr, args->target, args->target_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->target);
    }
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr_in  = bcf_sr_get_header(args->sr,0);

    mkdir_p("%s/",args->output_dir);

    int i,j, nsmpl = bcf_hdr_nsamples(args->hdr_in);
    if ( !nsmpl ) error("No samples to split: %s\n", args->fname);
    init_subsets(args);

    // parse tags
    int is_info = 0, is_fmt = 0;
    char *beg = args->keep_tags;
    while ( beg && *beg )
    {
        if ( !strncasecmp("INFO/",beg,5) ) { is_info = 1; is_fmt = 0; beg += 5; }
        else if ( !strcasecmp("INFO",beg) ) { args->keep_info = 1; break; }
        else if ( !strncasecmp("INFO,",beg,5) ) { args->keep_info = 1; beg += 5; continue; }
        else if ( !strncasecmp("FMT/",beg,4) ) { is_info = 0; is_fmt = 1; beg += 4; }
        else if ( !strncasecmp("FORMAT/",beg,7) ) { is_info = 0; is_fmt = 1; beg += 7; }
        else if ( !strcasecmp("FMT",beg) ) { args->keep_fmt = 1; break; }
        else if ( !strcasecmp("FORMAT",beg) ) { args->keep_fmt = 1; break; }
        else if ( !strncasecmp("FMT,",beg,4) ) { args->keep_fmt = 1; beg += 4; continue; }
        else if ( !strncasecmp("FORMAT,",beg,7) ) { args->keep_fmt = 1; beg += 7; continue; }
        char *end = beg;
        while ( *end && *end!=',' ) end++;
        char tmp = *end; *end = 0;
        int id = bcf_hdr_id2int(args->hdr_in, BCF_DT_ID, beg);
        beg = tmp ? end + 1 : end;
        if ( is_info && bcf_hdr_idinfo_exists(args->hdr_in,BCF_HL_INFO,id) )
        {
            if ( id >= args->ninfo_tags ) args->ninfo_tags = id + 1;
            hts_expand0(uint8_t, args->ninfo_tags, args->minfo_tags, args->info_tags);
            args->info_tags[id] = 1;
        }
        if ( is_fmt && bcf_hdr_idinfo_exists(args->hdr_in,BCF_HL_FMT,id) )
        {
            if ( id >= args->nfmt_tags ) args->nfmt_tags = id + 1;
            hts_expand0(uint8_t, args->nfmt_tags, args->mfmt_tags, args->fmt_tags);
            args->fmt_tags[id] = 1;
        }
    }
    if ( !args->keep_info && !args->keep_fmt && !args->ninfo_tags && !args->nfmt_tags )
    {
        args->keep_info = args->keep_fmt = 1;
    }
    if ( !args->keep_fmt && !args->nfmt_tags ) args->keep_fmt = 1;

    bcf_hdr_t *tmp_hdr = bcf_hdr_dup(args->hdr_in);
    if ( !args->keep_info || args->ninfo_tags || args->nfmt_tags )
    {
        int j;
        for (j=tmp_hdr->nhrec-1; j>=0; j--)
        {
            bcf_hrec_t *hrec = tmp_hdr->hrec[j];
            if ( hrec->type!=BCF_HL_INFO && hrec->type!=BCF_HL_FMT ) continue;
            int k = bcf_hrec_find_key(hrec,"ID");
            assert( k>=0 ); // this should always be true for valid VCFs
            int remove = 0;
            if ( hrec->type==BCF_HL_INFO && (!args->keep_info || args->ninfo_tags) )
            {
                int id = bcf_hdr_id2int(tmp_hdr, BCF_DT_ID, hrec->vals[k]);
                if ( !args->keep_info || id >= args->ninfo_tags || !args->info_tags[id] ) remove = 1;
            }
            if ( hrec->type==BCF_HL_FMT && args->nfmt_tags )
            {
                int id = bcf_hdr_id2int(tmp_hdr, BCF_DT_ID, hrec->vals[k]);
                if ( id >= args->nfmt_tags || !args->fmt_tags[id] ) remove = 1;
            }
            if ( remove )
            {
                char *str = strdup(hrec->vals[k]);
                bcf_hdr_remove(tmp_hdr,hrec->type,str);
                free(str);
            }
        }
        if ( bcf_hdr_sync(tmp_hdr)!=0 ) error("Failed to update the VCF header\n");
    }

    kstring_t str = {0,0,0};
    for (i=0; i<args->nsets; i++)
    {
        subset_t *set = &args->sets[i];
        str.l = 0;
        kputs(args->output_dir, &str);
        if ( str.s[str.l-1] != '/' ) kputc('/', &str);
        kputs(set->fname, &str);
        char *suffix = NULL;
        if ( args->output_type & FT_BCF ) suffix = ".bcf";
        else if ( args->output_type & FT_GZ ) suffix = ".vcf.gz";
        else suffix = ".vcf";
        int len = strlen(set->fname);
        if ( len >= 4 && !strcasecmp(".bcf",set->fname+len-4) ) suffix = NULL;
        if ( len >= 4 && !strcasecmp(".vcf",set->fname+len-4) ) suffix = NULL;
        if ( len >= 7 && !strcasecmp(".vcf.gz",set->fname+len-7) ) suffix = NULL;
        if ( len >= 8 && !strcasecmp(".vcf.bgz",set->fname+len-8) ) suffix = NULL;
        if ( suffix ) kputs(suffix, &str);
        char wmode[8];
        set_wmode(wmode,args->output_type,str.s,args->clevel);
        set->fh = hts_open(str.s, wmode);
        if ( set->fh == NULL ) error("[%s] Error: cannot write to \"%s\": %s\n", __func__, str.s, strerror(errno));
        if ( args->hts_opts )
        {
            hts_opt *opts = NULL;
            for (j=0; j<args->nhts_opts; j++)
                if ( hts_opt_add(&opts,args->hts_opts[j]) ) error("Could not set the HTS option \"%s\"\n",args->hts_opts[j]);
            if ( hts_opt_apply(set->fh,opts) ) error("Could not apply the HTS options\n");
            hts_opt_free(opts);
        }
        set->hdr = tmp_hdr;     // dirty: reuse the same header to lower memory for large datasets
        bcf_hdr_nsamples(set->hdr) = set->nsmpl;
        for (j=0; j<set->nsmpl; j++)
            set->hdr->samples[j] = set->rename ? set->rename[j] : args->hdr_in->samples[set->smpl[j]];
        if ( bcf_hdr_write(set->fh, set->hdr)!=0 ) error("[%s] Error: cannot write the header to %s\n", __func__,str.s);
        if ( init_index2(set->fh,set->hdr,str.s,&set->index_fn,
                         args->write_index)<0 )
            error("Error: failed to initialise index for %s\n",str.s);
        if ( args->filter_str )
            set->filter = filter_init(set->hdr, args->filter_str);
    }
    args->hdr_out = tmp_hdr;
    free(str.s);
}
static void destroy_data(args_t *args)
{
    bcf_hdr_destroy(args->hdr_out);
    free(args->info_tags);
    free(args->fmt_tags);
    int i,j;
    for (i=0; i<args->nsets; i++)
    {
        subset_t *set = &args->sets[i];
        if ( args->write_index )
        {
            if ( bcf_idx_save(set->fh)<0 )
            {
                if ( hts_close(set->fh)!=0 ) error("Error: close failed .. %s\n", set->fname);
                error("Error: cannot write to index %s\n", set->index_fn);
            }
            free(set->index_fn);
        }
        if ( hts_close(set->fh)!=0 ) error("Error: close failed .. %s\n",set->fname);
        free(set->fname);
        free(set->smpl);
        if ( set->filter )
            filter_destroy(set->filter);
        //bcf_hdr_destroy(set->hdr);
        if ( set->rename )
        {
            for (j=0; j<set->nsmpl; j++) free(set->rename[j]);
            free(set->rename);
        }
    }
    for (i=0; i<args->nhts_opts; i++) free(args->hts_opts[i]);
    free(args->hts_opts);
    free(args->sets);
    bcf_sr_destroy(args->sr);
    free(args);
}

static bcf1_t *rec_set_info(args_t *args, subset_t *set, bcf1_t *rec, bcf1_t *out)
{
    if ( out )
    {
        out->n_sample = set->nsmpl;
        return out;
    }
    out = bcf_init1();
    out->rid  = rec->rid;
    out->pos  = rec->pos;
    out->rlen = rec->rlen;
    out->qual = rec->qual;
    out->n_allele = rec->n_allele;
    out->n_sample = set->nsmpl;
    if ( args->keep_info )
    {
        out->n_info = rec->n_info;
        out->shared.m = out->shared.l = rec->shared.l;
        out->shared.s = (char*) malloc(out->shared.l);
        memcpy(out->shared.s,rec->shared.s,out->shared.l);
        return out;
    }

    // build the BCF record
    kstring_t tmp = {0,0,0};
    char *ptr = rec->shared.s;
    kputsn_(ptr, rec->unpack_size[0], &tmp); ptr += rec->unpack_size[0]; // ID
    kputsn_(ptr, rec->unpack_size[1], &tmp); ptr += rec->unpack_size[1]; // REF+ALT
    kputsn_(ptr, rec->unpack_size[2], &tmp);                             // FILTER
    if ( args->ninfo_tags )
    {
        int i;
        for (i=0; i<rec->n_info; i++)
        {
            bcf_info_t *info = &rec->d.info[i];
            int id = info->key;
            if ( !args->info_tags[id] ) continue;
            kputsn_(info->vptr - info->vptr_off, info->vptr_len + info->vptr_off, &tmp);
            out->n_info++;
        }
    }
    out->shared.m = tmp.m;
    out->shared.s = tmp.s;
    out->shared.l = tmp.l;
    out->unpacked = 0;
    return out;
}

static bcf1_t *rec_set_format(args_t *args, subset_t *set, bcf1_t *src, bcf1_t *dst)
{
    dst->unpacked &= ~BCF_UN_FMT;
    dst->n_fmt = 0;
    kstring_t tmp = dst->indiv; tmp.l = 0;
    int i,j;
    for (i=0; i<src->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &src->d.fmt[i];
        int id = fmt->id;
        if ( !args->keep_fmt && (id>=args->nfmt_tags || !args->fmt_tags[id]) ) continue;

        bcf_enc_int1(&tmp, id);
        bcf_enc_size(&tmp, fmt->n, fmt->type);
        for (j=0; j<set->nsmpl; j++)
            kputsn_(fmt->p + set->smpl[j]*fmt->size, fmt->size, &tmp);

        dst->n_fmt++;
    }
    dst->indiv = tmp;
    return dst;
}

static void process(args_t *args)
{
    bcf1_t *rec = bcf_sr_get_line(args->sr,0);
    bcf_unpack(rec, BCF_UN_ALL);

    int i;
    bcf1_t *out = NULL;
    for (i=0; i<args->nsets; i++)
    {
        subset_t *set = &args->sets[i];
        bcf_hdr_nsamples(set->hdr) = set->nsmpl;    // dirty: prevent memory duplication for very large sets

        out = rec_set_info(args, set, rec, out);
        rec_set_format(args, set, rec, out);

        int pass = 1;
        if ( set->filter )
        {
            pass = filter_test(set->filter, out, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
        }
        if ( !pass ) continue;
        if ( bcf_write(set->fh, set->hdr, out)!=0 ) error("[%s] Error: failed to write the record\n", __func__);
    }
    if ( out ) bcf_destroy(out);
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_type  = FT_VCF;
    args->regions_overlap = 1;
    args->targets_overlap = 0;
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"hts-opts",required_argument,NULL,1},
        {"keep-tags",required_argument,NULL,'k'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,2},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"targets-overlap",required_argument,NULL,3},
        {"samples-file",required_argument,NULL,'S'},
        {"groups-file",required_argument,NULL,'G'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "vr:R:t:T:o:O:i:e:k:S:G:W::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->hts_opts = hts_readlist(optarg,0,&args->nhts_opts); break;
            case 'k': args->keep_tags = optarg; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'T': args->target = optarg; args->target_is_file = 1; break;
            case 't': args->target = optarg; break;
            case 'R': args->region = optarg; args->region_is_file = 1;  break;
            case 'S': args->samples_fname = optarg; break;
            case 'G': args->groups_fname = optarg; break;
            case 'r': args->region = optarg; break;
            case 'o': args->output_dir = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      }
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case  2 :
                args->regions_overlap = parse_overlap_option(optarg);
                if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  3 :
                args->targets_overlap = parse_overlap_option(optarg);
                if ( args->targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( args->samples_fname && args->groups_fname ) error("Error: only one of the options --groups-file or --samples-file can be given\n");
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s", usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s", usage_text());
    else args->fname = argv[optind];

    if ( !args->output_dir ) error("Missing the -o option\n");
    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");

    init_data(args);

    while ( bcf_sr_next_line(args->sr) ) process(args);

    destroy_data(args);
    return 0;
}


