/*  vcfmerge.c -- Merge multiple VCF/BCF files to create one multi-sample file.

    Copyright (C) 2012-2014 Genome Research Ltd.

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
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <math.h>
#include "bcftools.h"
#include "vcmp.h"

#include <htslib/khash.h>
KHASH_MAP_INIT_STR(strdict, int)
typedef khash_t(strdict) strdict_t;

#define SKIP_DONE 1
#define SKIP_DIFF 2

#define IS_VL_G(hdr,id) (bcf_hdr_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_G)
#define IS_VL_A(hdr,id) (bcf_hdr_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_A)

// For merging INFO Number=A,G,R tags
typedef struct
{
    const char *hdr_tag;
    int type, nvals;
    int nbuf, mbuf;
    uint8_t *buf;
}
AGR_info_t;

// Rules for merging arbitrary INFO tags
typedef struct _info_rule_t
{
    char *hdr_tag;
    void (*merger)(bcf_hdr_t *hdr, bcf1_t *line, struct _info_rule_t *rule);
    int type;           // one of BCF_HT_*
    int block_size;     // number of values in a block
    int nblocks;        // number of blocks in nvals (the number of merged files)
    int nvals, mvals;   // used and total size of vals array
    void *vals;         // the info tag values
}
info_rule_t;

// Auxiliary merge data for selecting the right combination
//  of buffered records across multiple readers. maux1_t
//  corresponds to one buffered line.
typedef struct
{
    int skip;
    int *map;   // mapping from input alleles to the output array
    int mmap;   // size of map array (only buffer[i].n_allele is actually used)
    int als_differ;
}
maux1_t;
typedef struct
{
    int n;  // number of readers
    char **als, **out_als;  // merged alleles (temp, may contain empty records) and merged alleles ready for output
    int nals, mals, nout_als, mout_als; // size of the output array
    int *cnt, ncnt; // number of records that refer to the alleles
    int *nbuf;      // readers have buffers of varying lengths
    int *smpl_ploidy, *smpl_nGsize; // ploidy and derived number of values in Number=G tags, updated for each line (todo: cache for missing cases)
    int *flt, mflt, minf;
    bcf_info_t *inf;// out_line's INFO fields
    bcf_fmt_t **fmt_map; // i-th output FORMAT field corresponds in j-th reader to i*nreader+j, first row is reserved for GT
    int nfmt_map;        // number of rows in the fmt_map array
    int *agr_map, nagr_map, magr_map;   // mapping between Number=AGR element indexes
    void *tmp_arr;
    int ntmp_arr;
    maux1_t **d;    // d[i][j] i-th reader, j-th buffer line
    AGR_info_t *AGR_info;
    int nAGR_info, mAGR_info;
    bcf_srs_t *files;
    int *has_line;  // which files are being merged
}
maux_t;

typedef struct
{
    vcmp_t *vcmp;
    maux_t *maux;
    int header_only, collapse, output_type, force_samples, merge_by_id;
    char *header_fname, *output_fname, *regions_list, *info_rules, *file_list;
    info_rule_t *rules;
    int nrules;
    strdict_t *tmph;
    kstring_t tmps;
    bcf_srs_t *files;
    bcf1_t *out_line;
    htsFile *out_fh;
    bcf_hdr_t *out_hdr;
    char **argv;
    int argc;
}
args_t;

static void info_rules_merge_sum(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = 0; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) ptr[j] += ptr[j+i*ndim]; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_avg(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = 0; \
        for (j=0; j<ndim; j++) \
        { \
            double sum = 0; \
            for (i=0; i<rule->nblocks; i++) sum += ptr[j+i*ndim]; \
            ptr[j] = sum / rule->nblocks; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_min(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing,set_missing,huge_val) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = huge_val; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) if ( ptr[j] > ptr[j+i*ndim] ) ptr[j] = ptr[j+i*ndim]; \
        } \
        for (i=0; i<rule->nvals; i++) if ( ptr[i]==huge_val ) set_missing; \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing, ptr[i]=bcf_int32_missing, INT32_MAX); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i]), bcf_float_set_missing(ptr[i]), HUGE_VAL); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_max(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing,set_missing,huge_val) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = huge_val; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) if ( ptr[j] < ptr[j+i*ndim] ) ptr[j] = ptr[j+i*ndim]; \
        } \
        for (i=0; i<rule->nvals; i++) if ( ptr[i]==huge_val ) set_missing; \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing, ptr[i]=bcf_int32_missing, INT32_MIN); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i]), bcf_float_set_missing(ptr[i]), -HUGE_VAL); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_join(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    if ( rule->type==BCF_HT_STR )
    {
        ((char*)rule->vals)[rule->nvals] = 0;
        bcf_update_info_string(hdr,line,rule->hdr_tag,rule->vals);
    }
    else
        bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,rule->nvals,rule->type);
}

static int info_rules_comp_key2(const void *a, const void *b)
{
    info_rule_t *rule1 = (info_rule_t*) a;
    info_rule_t *rule2 = (info_rule_t*) b;
    return strcmp(rule1->hdr_tag, rule2->hdr_tag);
}
static int info_rules_comp_key(const void *a, const void *b)
{
    char *key = (char*) a;
    info_rule_t *rule = (info_rule_t*) b;
    return strcmp(key, rule->hdr_tag);
}
static void info_rules_init(args_t *args)
{
    if ( args->info_rules && !strcmp("-",args->info_rules) ) return;

    kstring_t str = {0,0,0};
    if ( !args->info_rules )
    {
        if ( bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "DP")) ) kputs("DP:sum",&str);
        if ( bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "DP4")) )
        {
            if ( str.l ) kputc(',',&str);
            kputs("DP4:sum",&str);
        }
        if ( !str.l ) return;
        args->info_rules = str.s;
    }

    args->nrules = 1;
    char *ss = strdup(args->info_rules), *tmp = ss;
    int n = 0;
    while ( *ss )
    {
        if ( *ss==':' ) { *ss = 0; n++; if ( n%2==0 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules); }
        else if ( *ss==',' ) { *ss = 0; args->nrules++; n++; if ( n%2==1 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules); }
        ss++;
    }
    if ( n%2==0 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules);
    args->rules = (info_rule_t*) calloc(args->nrules,sizeof(info_rule_t));

    n = 0;
    ss = tmp;
    while ( n < args->nrules )
    {
        info_rule_t *rule = &args->rules[n];
        rule->hdr_tag = strdup(ss);
        int id = bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, rule->hdr_tag);
        if ( !bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,id) ) error("The tag is not defined in the header: \"%s\"\n", rule->hdr_tag);
        rule->type = bcf_hdr_id2type(args->out_hdr,BCF_HL_INFO,id);
        if ( rule->type!=BCF_HT_INT && rule->type!=BCF_HT_REAL && rule->type!=BCF_HT_STR ) error("The type is not supported: \"%s\"\n", rule->hdr_tag);

        while ( *ss ) ss++; ss++;
        if ( !*ss ) error("Could not parse INFO rules, missing logic of \"%s\"\n", rule->hdr_tag);

        int is_join = 0;
        if ( !strcasecmp(ss,"sum") ) rule->merger = info_rules_merge_sum;
        else if ( !strcasecmp(ss,"avg") ) rule->merger = info_rules_merge_avg;
        else if ( !strcasecmp(ss,"min") ) rule->merger = info_rules_merge_min;
        else if ( !strcasecmp(ss,"max") ) rule->merger = info_rules_merge_max;
        else if ( !strcasecmp(ss,"join") ) { rule->merger = info_rules_merge_join; is_join = 1; }
        else error("The rule logic \"%s\" not recognised\n", ss);

        if ( !is_join && rule->type==BCF_HT_STR )
            error("Numeric operation \"%s\" requested on non-numeric field: %s\n", ss, rule->hdr_tag);
        if ( bcf_hdr_id2number(args->out_hdr,BCF_HL_INFO,id)==0xfffff )
        {
            int is_agr = (
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_A ||
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_G ||
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_R
                    ) ? 1 : 0;
            if ( is_join && is_agr )
                error("Cannot -i %s:join on Number=[AGR] tags is not supported.\n", rule->hdr_tag);
            if ( !is_join && !is_agr )
                error("Only fixed-length vectors are supported with -i %s:%s\n", ss, rule->hdr_tag);
        }

        while ( *ss ) ss++; ss++; n++;
    }
    free(str.s);
    free(tmp);

    qsort(args->rules, args->nrules, sizeof(*args->rules), info_rules_comp_key2);
}
static void info_rules_destroy(args_t *args)
{
    int i;
    for (i=0; i<args->nrules; i++)
    {
        info_rule_t *rule = &args->rules[i];
        free(rule->hdr_tag);
        free(rule->vals);
    }
    free(args->rules);
}
static void info_rules_reset(args_t *args)
{
    int i;
    for (i=0; i<args->nrules; i++)
        args->rules[i].nblocks = args->rules[i].nvals = args->rules[i].block_size = 0;
}
static int info_rules_add_values(args_t *args, bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule, maux1_t *als, int var_len)
{
    int ret = bcf_get_info_values(hdr, line, rule->hdr_tag, &args->maux->tmp_arr, &args->maux->ntmp_arr, rule->type);
    if ( ret<=0 ) error("FIXME: error parsing %s at %s:%d .. %d\n", rule->hdr_tag,bcf_seqname(hdr,line),line->pos+1,ret);

    rule->nblocks++;

    if ( rule->type==BCF_HT_STR )
    {
        int need_comma = rule->nblocks==1 ? 0 : 1;
        hts_expand(char,rule->nvals+ret+need_comma+1,rule->mvals,rule->vals);            // 1 for null-termination
        char *tmp = (char*) rule->vals + rule->nvals;
        if ( rule->nvals>0 ) { *tmp = ','; tmp++; }
        strncpy(tmp,(char*)args->maux->tmp_arr,ret);
        rule->nvals += ret + need_comma;
        return 1;
    }

    int i, j;
    if ( var_len==BCF_VL_A )
    {
        assert( ret==line->n_allele-1 );
        args->maux->nagr_map = ret;
        hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
        // create mapping from source file ALT indexes to dst file indexes
        for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i+1] - 1;
        rule->block_size = args->maux->nout_als - 1;
    }
    else if ( var_len==BCF_VL_R )
    {
        assert( ret==line->n_allele );
        args->maux->nagr_map = ret;
        hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
        for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i];
        rule->block_size = args->maux->nout_als;
    }
    else if ( var_len==BCF_VL_G )
    {
        args->maux->nagr_map = bcf_alleles2gt(line->n_allele-1,line->n_allele-1)+1;
        assert( ret==line->n_allele || ret==args->maux->nagr_map );
        if ( ret==line->n_allele ) // haploid
        {
            args->maux->nagr_map = line->n_allele;
            hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
            for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i];
            rule->block_size = args->maux->nout_als;
        }
        else
        {
            hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
            int k_src = 0;
            for (i=0; i<line->n_allele; i++)
            {
                for (j=0; j<=i; j++)
                {
                    args->maux->agr_map[k_src] = bcf_alleles2gt(als->map[i],als->map[j]);
                    k_src++;
                }
            }
            rule->block_size = bcf_alleles2gt(args->maux->nout_als-1,args->maux->nout_als-1)+1;
        }
    }
    else
    {
        if ( rule->nblocks>1 && ret!=rule->block_size )
            error("Mismatch in number of values for INFO/%s at %s:%d\n", rule->hdr_tag,bcf_seqname(hdr,line),line->pos+1);
        rule->block_size = ret;
        args->maux->nagr_map = 0;
    }

    #define BRANCH(src_type_t,dst_type_t,set_missing) { \
        src_type_t *src = (src_type_t *) args->maux->tmp_arr; \
        hts_expand0(dst_type_t,(rule->nvals+rule->block_size),rule->mvals,rule->vals); \
        dst_type_t *dst = (dst_type_t *) rule->vals + rule->nvals; \
        rule->nvals += rule->block_size; \
        if ( !args->maux->nagr_map ) \
        { \
            for (i=0; i<ret; i++) dst[i] = src[i]; \
        } \
        else \
        { \
            for (i=0; i<rule->block_size; i++) set_missing; \
            for (i=0; i<ret; i++) dst[args->maux->agr_map[i]] = src[i]; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int, int32_t, dst[i] = bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, float, bcf_float_set_missing(dst[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    return 1;
}

int bcf_hdr_sync(bcf_hdr_t *h);

void bcf_hdr_merge(bcf_hdr_t *hw, const bcf_hdr_t *hr, const char *clash_prefix, int force_samples)
{
    // header lines
    int ret = bcf_hdr_combine(hw, hr);
    if ( ret!=0 ) error("Error occurred while merging the headers.\n");

    // samples
    int i;
    for (i=0; i<bcf_hdr_nsamples(hr); i++)
    {
        char *name = hr->samples[i];
        if ( bcf_hdr_id2int(hw, BCF_DT_SAMPLE, name)!=-1 )
        {
            // there is a sample with the same name
            if ( !force_samples ) error("Error: Duplicate sample names (%s), use --force-samples to proceed anyway.\n", name);

            int len = strlen(hr->samples[i]) + strlen(clash_prefix) + 1;
            name = (char*) malloc(sizeof(char)*(len+1));
            sprintf(name,"%s:%s",clash_prefix,hr->samples[i]);
            bcf_hdr_add_sample(hw,name);
            free(name);
        }
        else
            bcf_hdr_add_sample(hw,name);
    }
}

void debug_als(char **als, int nals)
{
    int k; for (k=0; k<nals; k++) fprintf(stderr,"%s ", als[k]);
    fprintf(stderr,"\n");
}

/**
 * normalize_alleles() - create smallest possible representation of the alleles
 * @als:    alleles to be merged, first is REF (rw)
 * @nals:   number of $a alleles
 *
 * Best explained on an example:
 *      In:  REF=GTTT  ALT=GTT
 *      Out: REF=GT    ALT=G
 *
 * Note: the als array will be modified
 */
void normalize_alleles(char **als, int nals)
{
    int j, i = 1, done = 0, rlen = strlen(als[0]);
    while ( i<rlen )
    {
        for (j=1; j<nals; j++)
        {
            int len = strlen(als[j]);
            if ( i>=len ) done = 1;
            if ( als[j][len-i] != als[0][rlen-i] ) { done = 1; break; }
        }
        if ( done ) break;
        i++;
    }
    if ( i>1 )
    {
        i--;
        als[0][rlen-i] = 0;
        for (j=1; j<nals; j++) als[j][strlen(als[j])-i] = 0;
    }
}

 /**
 * merge_alleles() - merge two REF,ALT records, $a and $b into $b.
 * @a:      alleles to be merged, first is REF
 * @na:     number of $a alleles
 * @map:    map from the original $a indexes to new $b indexes (0-based)
 * @b:      alleles to be merged, the array will be expanded as required
 * @nb:     number of $b alleles
 * @mb:     size of $b
 *
 * Returns NULL on error or $b expanded to incorporate $a alleles and sets
 * $map. Best explained on an example:
 *      In:     REF   ALT
 *           a: ACG,  AC,A    (1bp and 2bp deletion)
 *           b: ACGT, A       (3bp deletion)
 *      Out:
 *           b: ACGT, A,ACT,AT (3bp, 1bp and 2bp deletion)
 *           map: 0,2,3
 * Here the mapping from the original $a alleles to the new $b alleles is 0->0,
 * 1->2, and 2->3.
 */
char **merge_alleles(char **a, int na, int *map, char **b, int *nb, int *mb)
{
    // reference allele never changes
    map[0] = 0;

    int rla = !a[0][1] ? 1 : strlen(a[0]);
    int rlb = !b[0][1] ? 1 : strlen(b[0]);

    // the most common case: same SNPs
    if ( na==2 && *nb==2 && rla==1 && rlb==1 && a[1][0]==b[1][0] && !a[1][1] && !b[1][1] )
    {
        map[1] = 1;
        return b;
    }

    // Sanity check: reference prefixes must be identical
    if ( strncmp(a[0],b[0],rla<rlb?rla:rlb) )
    {
        fprintf(stderr, "The REF prefixes differ: %s vs %s (%d,%d)\n", a[0],b[0],rla,rlb);
        return NULL;
    }

    int n = *nb + na;
    hts_expand0(char*,n,*mb,b);

    // $b alleles need expanding
    int i,j;
    if ( rla>rlb )
    {
        for (i=0; i<*nb; i++)
        {
            int l = strlen(b[i]);
            b[i] = (char*) realloc(b[i],l+rla-rlb+1);
            memcpy(b[i]+l,a[0]+rlb,rla-rlb+1);
        }
    }

    // now check if the $a alleles are present and if not add them
    for (i=1; i<na; i++)
    {
        char *ai;
        if ( rlb>rla )  // $a alleles need expanding
        {
            int l = strlen(a[i]);
            ai = (char*) malloc(l+rlb-rla+1);
            memcpy(ai,a[i],l);
            memcpy(ai+l,b[0]+rla,rlb-rla+1);
        }
        else
            ai = a[i];

        for (j=1; j<*nb; j++)
            if ( !strcmp(ai,b[j]) ) break;

        if ( j<*nb ) // $b already has the same allele
        {
            map[i] = j;
            if ( rlb>rla ) free(ai);
            continue;
        }
        // new allele
        map[i] = *nb;
        b[*nb] = rlb>rla ? ai : strdup(ai);
        (*nb)++;
    }
    return b;
}

maux_t *maux_init(bcf_srs_t *files)
{
    maux_t *ma = (maux_t*) calloc(1,sizeof(maux_t));
    ma->n      = files->nreaders;
    ma->nbuf   = (int *) calloc(ma->n,sizeof(int));
    ma->d      = (maux1_t**) calloc(ma->n,sizeof(maux1_t*));
    ma->files  = files;
    int i, n_smpl = 0;
    for (i=0; i<ma->n; i++)
        n_smpl += bcf_hdr_nsamples(files->readers[i].header);
    ma->smpl_ploidy = (int*) calloc(n_smpl,sizeof(int));
    ma->smpl_nGsize = (int*) malloc(n_smpl*sizeof(int));
    ma->has_line = (int*) malloc(ma->n*sizeof(int));
    return ma;
}
void maux_destroy(maux_t *ma)
{
    int i;
    for (i=0; i<ma->n; i++) // for each reader
    {
        if ( !ma->d[i] ) continue;
        int j;
        for (j=0; j<ma->nbuf[i]; j++)  // for each buffered line
            if ( ma->d[i][j].map ) free(ma->d[i][j].map);
        free(ma->d[i]);
    }
    for (i=0; i<ma->mAGR_info; i++)
        free(ma->AGR_info[i].buf);
    free(ma->agr_map);
    free(ma->AGR_info);
    if (ma->ntmp_arr) free(ma->tmp_arr);
    if (ma->nfmt_map) free(ma->fmt_map);
    // ma->inf freed in bcf_destroy1
    free(ma->d);
    free(ma->nbuf);
    for (i=0; i<ma->mals; i++) free(ma->als[i]);
    if (ma->mout_als) free(ma->out_als);
    free(ma->als);
    free(ma->cnt);
    free(ma->smpl_ploidy);
    free(ma->smpl_nGsize);
    free(ma->has_line);
    free(ma);
}
void maux_expand1(maux_t *ma, int i)
{
    if ( ma->nbuf[i] <= ma->files->readers[i].nbuffer )
    {
        int n = ma->files->readers[i].nbuffer + 1;
        ma->d[i] = (maux1_t*) realloc(ma->d[i], sizeof(maux1_t)*n);
        memset(ma->d[i]+ma->nbuf[i],0,sizeof(maux1_t)*(n-ma->nbuf[i]));
        ma->nbuf[i] = n;
    }
}
void maux_reset(maux_t *ma)
{
    int i;
    for (i=0; i<ma->n; i++) maux_expand1(ma, i);
    for (i=1; i<ma->ncnt; i++) ma->cnt[i] = 0;
}
void maux_debug(maux_t *ma, int ir, int ib)
{
    printf("[%d,%d]\t", ir,ib);
    int i;
    for (i=0; i<ma->nals; i++)
    {
        printf(" %s [%d]", ma->als[i], ma->cnt[i]);
    }
    printf("\n");
}

void merge_chrom2qual(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
    kstring_t *tmps = &args->tmps;
    tmps->l = 0;

    maux_t *ma = args->maux;
    int *al_idxs = (int*) calloc(ma->nals,sizeof(int));
    bcf_float_set_missing(out->qual);

    // CHROM, POS, ID, QUAL
    out->pos = -1;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;

        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;

        // alleles
        int j;
        for (j=1; j<line->n_allele; j++)
            al_idxs[ ma->d[i][0].map[j] ] = 1;

        // position
        if ( out->pos==-1 )
        {
            const char *chr = hdr->id[BCF_DT_CTG][line->rid].key;
            out->rid = bcf_hdr_name2id(out_hdr, chr);
            if ( strcmp(chr,out_hdr->id[BCF_DT_CTG][out->rid].key) ) error("Uh\n");
            out->pos = line->pos;
        }

        // ID
        if ( line->d.id[0]!='.' || line->d.id[1] )
        {
            kitr = kh_get(strdict, tmph, line->d.id);
            if ( kitr == kh_end(tmph) )
            {
                if ( tmps->l ) kputc(';', tmps);
                kputs(line->d.id, tmps);
                kh_put(strdict, tmph, line->d.id, &ret);
            }
        }

        // set QUAL to the max qual value. Not exactly correct, but good enough for now
        if ( !bcf_float_is_missing(files->readers[i].buffer[0]->qual) )
        {
            if ( bcf_float_is_missing(out->qual) || out->qual < files->readers[i].buffer[0]->qual ) out->qual = files->readers[i].buffer[0]->qual;
        }
    }

    // set ID
    if ( !tmps->l ) kputs(".", tmps);
    if ( out->d.id ) free(out->d.id);
    out->d.id = strdup(tmps->s);

    // set alleles
    ma->nout_als = 0;
    for (i=1; i<ma->nals; i++)
    {
        if ( !al_idxs[i] ) continue;
        ma->nout_als++;

        // Adjust the indexes, the allele map could be created for multiple collapsed records,
        //  some of which might be unused for this output line
        int ir, j;
        for (ir=0; ir<files->nreaders; ir++)
        {
            if ( !ma->has_line[ir] ) continue;
            bcf1_t *line = files->readers[ir].buffer[0];
            for (j=1; j<line->n_allele; j++)
                if ( ma->d[ir][0].map[j]==i ) ma->d[ir][0].map[j] = ma->nout_als;
        }
    }
    // Expand the arrays and realloc the alleles string. Note that all alleles are in a single allocated block.
    ma->nout_als++;
    hts_expand0(char*, ma->nout_als, ma->mout_als, ma->out_als);
    int k = 0;
    for (i=0; i<ma->nals; i++)
        if ( i==0 || al_idxs[i] ) ma->out_als[k++] = strdup(ma->als[i]);
    assert( k==ma->nout_als );
    normalize_alleles(ma->out_als, ma->nout_als);
    bcf_update_alleles(out_hdr, out, (const char**) ma->out_als, ma->nout_als);
    free(al_idxs);
    for (i=0; i<ma->nout_als; i++) free(ma->out_als[i]);
}

void merge_filter(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    maux_t *ma = args->maux;
    out->d.n_flt = 0;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i]) continue;

        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        bcf_unpack(line, BCF_UN_ALL);

        int k;
        for (k=0; k<line->d.n_flt; k++)
        {
            const char *flt = hdr->id[BCF_DT_ID][line->d.flt[k]].key;
            kitr = kh_get(strdict, tmph, flt);
            if ( kitr == kh_end(tmph) )
            {
                int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, flt);
                if ( id==-1 ) error("The filter not defined: %s\n", flt);
                hts_expand(int,out->d.n_flt+1,ma->mflt,ma->flt);
                ma->flt[out->d.n_flt] = id;
                out->d.n_flt++;
                kh_put(strdict, tmph, flt, &ret);
            }
        }
    }
    // Check if PASS is not mixed with other filters
    if ( out->d.n_flt>1 )
    {
        int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "PASS");
        for (i=0; i<out->d.n_flt; i++)
            if ( ma->flt[i]==id ) break;
        if ( i<out->d.n_flt )
        {
            out->d.n_flt--;
            for (; i<out->d.n_flt; i++) ma->flt[i] = ma->flt[i+1];
        }
    }
    out->d.flt = ma->flt;
}

static void bcf_info_set_id(bcf1_t *line, bcf_info_t *info, int id, kstring_t *tmp_str)
{
    assert( !info->vptr_free );

    uint8_t *ptr = info->vptr - info->vptr_off;
    bcf_dec_typed_int1(ptr, &ptr);

    tmp_str->l = 0;
    bcf_enc_int1(tmp_str, id);

    if ( tmp_str->l == ptr - info->vptr + info->vptr_off )
    {
        // the new id is represented with the same number of bytes
        memcpy(info->vptr - info->vptr_off, tmp_str->s, tmp_str->l);
        return;
    }

    kputsn_(ptr, info->vptr - ptr, tmp_str);
    info->vptr_off = tmp_str->l;
    kputsn_(info->vptr, info->len << bcf_type_shift[info->type], tmp_str);

    info->vptr = (uint8_t*) tmp_str->s + info->vptr_off;
    info->vptr_free = 1;
    line->d.shared_dirty |= BCF1_DIRTY_INF;
    tmp_str->s = NULL;
    tmp_str->m = 0;
    tmp_str->l = 0;
}

void copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst)
{
    int ith_src = 0, start_src = 0;    // i-th field in src string
    while ( ith_src<isrc && start_src<src_len )
    {
        if ( src[start_src]==',' ) { ith_src++; }
        start_src++;
    }
    assert( ith_src==isrc ); // if ( ith_src<isrc ) return; // requested field not found
    int end_src = start_src;
    while ( end_src<src_len && src[end_src]!=',' ) end_src++;

    int nsrc_cpy = end_src - start_src;
    if ( nsrc_cpy==1 && src[start_src]=='.' ) return;   // don't write missing values, dst is already initialized

    int ith_dst = 0, start_dst = 0;
    while ( ith_dst<idst && start_dst<dst->l )
    {
        if ( dst->s[start_dst]==',' ) { ith_dst++; }
        start_dst++;
    }
    assert( ith_dst==idst ); // if ( ith_dst<idst ) return;
    int end_dst = start_dst;
    while ( end_dst<dst->l && dst->s[end_dst]!=',' ) end_dst++;

    if ( end_dst - start_dst>1 || dst->s[start_dst]!='.' ) return;   // do not overwrite non-empty values

    // Now start_dst and end_dst are indexes to the destination memory area
    // which needs to be replaced with nsrc_cpy
    // source bytes, end_dst points just after.
    int ndst_shift = nsrc_cpy - (end_dst - start_dst);
    int ndst_move  = dst->l - end_dst + 1;  // how many bytes must be moved (including \0)
    if ( ndst_shift )
    {
        ks_resize(dst, dst->l + ndst_shift + 1);    // plus \0
        memmove(dst->s+end_dst+ndst_shift, dst->s+end_dst, ndst_move);
    }
    memcpy(dst->s+start_dst, src+start_src, nsrc_cpy);
    dst->l += ndst_shift;
}

static void merge_AGR_info_tag(bcf1_t *line, bcf_info_t *info, int len, maux1_t *als, AGR_info_t *agr)
{
    int i;
    if ( !agr->nbuf )
    {
        if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 || info->type==BCF_BT_FLOAT )
        {
            agr->nbuf = 4 * agr->nvals;
            hts_expand(uint8_t,agr->nbuf,agr->mbuf,agr->buf);
            if ( info->type!=BCF_BT_FLOAT )
            {
                int32_t *tmp = (int32_t*) agr->buf;
                for (i=0; i<agr->nvals; i++) tmp[i] = bcf_int32_missing;
            }
            else
            {
                float *tmp = (float*) agr->buf;
                for (i=0; i<agr->nvals; i++) bcf_float_set_missing(tmp[i]);
            }
        }
        else if ( info->type==BCF_BT_CHAR )
        {
            kstring_t tmp; tmp.l = 0; tmp.m = agr->mbuf; tmp.s = (char*)agr->buf;
            kputc('.',&tmp);
            for (i=1; i<agr->nvals; i++) kputs(",.",&tmp);
            agr->mbuf = tmp.m; agr->nbuf = tmp.l; agr->buf = (uint8_t*)tmp.s;
        }
        else
            error("Not ready for type [%d]: %s at %d\n", info->type,agr->hdr_tag,line->pos+1);
    }

    if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 || info->type==BCF_BT_FLOAT )
    {
        if ( len==BCF_VL_A || len==BCF_VL_R )
        {
            int ifrom = len==BCF_VL_A ? 1 : 0;
            #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
                type_t *src = (type_t *) info->vptr; \
                out_type_t *tgt = (out_type_t *) agr->buf; \
                int iori, inew; \
                for (iori=ifrom; iori<line->n_allele; iori++) \
                { \
                    if ( is_vector_end ) break; \
                    if ( is_missing ) continue; \
                    inew = als->map[iori] - ifrom; \
                    tgt[inew] = *src; \
                    src++; \
                } \
            }
            switch (info->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  *src==bcf_int8_missing,  *src==bcf_int8_vector_end,  int); break;
                case BCF_BT_INT16: BRANCH(int16_t, *src==bcf_int16_missing, *src==bcf_int16_vector_end, int); break;
                case BCF_BT_INT32: BRANCH(int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, int); break;
                case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), float); break;
                default: fprintf(stderr,"TODO: %s:%d .. info->type=%d\n", __FILE__,__LINE__, info->type); exit(1);
            }
            #undef BRANCH
        }
        else
        {
            #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
                type_t *src = (type_t *) info->vptr; \
                out_type_t *tgt = (out_type_t *) agr->buf; \
                int iori,jori, inew,jnew; \
                for (iori=0; iori<line->n_allele; iori++) \
                { \
                    inew = als->map[iori]; \
                    for (jori=0; jori<=iori; jori++) \
                    { \
                        jnew = als->map[jori]; \
                        int kori = iori*(iori+1)/2 + jori; \
                        if ( is_vector_end ) break; \
                        if ( is_missing ) continue; \
                        int knew = inew>jnew ? inew*(inew+1)/2 + jnew : jnew*(jnew+1)/2 + inew; \
                        tgt[knew] = src[kori]; \
                    } \
                    if ( jori<=iori ) break; \
                } \
            }
            switch (info->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  src[kori]==bcf_int8_missing,  src[kori]==bcf_int8_vector_end,  int); break;
                case BCF_BT_INT16: BRANCH(int16_t, src[kori]==bcf_int16_missing, src[kori]==bcf_int16_vector_end, int); break;
                case BCF_BT_INT32: BRANCH(int32_t, src[kori]==bcf_int32_missing, src[kori]==bcf_int32_vector_end, int); break;
                case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(src[kori]), bcf_float_is_vector_end(src[kori]), float); break;
                default: fprintf(stderr,"TODO: %s:%d .. info->type=%d\n", __FILE__,__LINE__, info->type); exit(1);
            }
            #undef BRANCH
        }
    }
    else
    {
        kstring_t tmp; tmp.l = agr->nbuf; tmp.m = agr->mbuf; tmp.s = (char*)agr->buf;
        if ( len==BCF_VL_A || len==BCF_VL_R )
        {
            int iori, ifrom = len==BCF_VL_A ? 1 : 0;
            for (iori=ifrom; iori<line->n_allele; iori++)
                copy_string_field((char*)info->vptr, iori-ifrom, info->len, &tmp, als->map[iori]-ifrom);
        }
        else
        {
            int iori,jori, inew,jnew;
            for (iori=0; iori<line->n_allele; iori++)
            {
                inew = als->map[iori];
                for (jori=0; jori<=iori; jori++)
                {
                    jnew = als->map[jori];
                    int kori = iori*(iori+1)/2 + jori;
                    int knew = bcf_alleles2gt(inew,jnew);
                    copy_string_field((char*)info->vptr, kori, info->len, &tmp, knew);
                }
            }
        }
        agr->mbuf = tmp.m; agr->nbuf = tmp.l; agr->buf = (uint8_t*)tmp.s;
    }
}

void merge_info(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, j, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    maux_t *ma = args->maux;
    ma->nAGR_info = 0;
    out->n_info   = 0;
    info_rules_reset(args);
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        for (j=0; j<line->n_info; j++)
        {
            bcf_info_t *inf = &line->d.info[j];

            const char *key = hdr->id[BCF_DT_ID][inf->key].key;
            if ( !strcmp("AC",key) || !strcmp("AN",key) ) continue;  // AC and AN are done in merge_format() after genotypes are done

            int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, key);
            if ( id==-1 ) error("Error: The INFO field is not defined in the header: %s\n", key);

            kitr = kh_get(strdict, tmph, key);  // have we seen the tag in one of the readers?
            int len = bcf_hdr_id2length(hdr,BCF_HL_INFO,inf->key);
            if ( args->nrules )
            {
                info_rule_t *rule = (info_rule_t*) bsearch(key, args->rules, args->nrules, sizeof(*args->rules), info_rules_comp_key);
                if ( rule )
                {
                    maux1_t *als = ( len==BCF_VL_A || len==BCF_VL_G || len==BCF_VL_R ) ? &ma->d[i][0] : NULL;
                    if ( info_rules_add_values(args, hdr, line, rule, als, len) ) continue;
                }
            }

            // Todo: Number=AGR tags should use the newer info_rules_* functions (info_rules_merge_first to be added)
            // and merge_AGR_info_tag to be made obsolete.
            if ( len==BCF_VL_A || len==BCF_VL_G || len==BCF_VL_R  ) // Number=R,G,A requires special treatment
            {
                if ( kitr == kh_end(tmph) )
                {
                    // first occurance in this reader, alloc arrays
                    ma->nAGR_info++;
                    hts_expand0(AGR_info_t,ma->nAGR_info,ma->mAGR_info,ma->AGR_info);
                    kitr = kh_put(strdict, tmph, key, &ret);
                    kh_val(tmph,kitr) = ma->nAGR_info - 1;
                    ma->AGR_info[ma->nAGR_info-1].hdr_tag = key;
                    ma->AGR_info[ma->nAGR_info-1].type  = bcf_hdr_id2type(hdr,BCF_HL_INFO,inf->key);
                    ma->AGR_info[ma->nAGR_info-1].nbuf  = 0;    // size of the buffer
                    switch (len)
                    {
                        case BCF_VL_A: ma->AGR_info[ma->nAGR_info-1].nvals = ma->nout_als - 1; break;
                        case BCF_VL_G: ma->AGR_info[ma->nAGR_info-1].nvals = bcf_alleles2gt(ma->nout_als-1,ma->nout_als-1)+1; break;
                        case BCF_VL_R: ma->AGR_info[ma->nAGR_info-1].nvals = ma->nout_als; break;
                    }
                }
                kitr = kh_get(strdict, tmph, key);
                int idx = kh_val(tmph, kitr);
                if ( idx<0 ) error("Error occurred while processing INFO tag \"%s\" at %s:%d\n", key,bcf_seqname(hdr,line),line->pos+1);
                merge_AGR_info_tag(line,inf,len,&ma->d[i][0],&ma->AGR_info[idx]);
                continue;
            }

            if ( kitr == kh_end(tmph) )
            {
                hts_expand0(bcf_info_t,out->n_info+1,ma->minf,ma->inf);
                ma->inf[out->n_info].key  = id;
                ma->inf[out->n_info].type = inf->type;
                ma->inf[out->n_info].len  = inf->len;
                ma->inf[out->n_info].vptr = inf->vptr;
                ma->inf[out->n_info].v1.i = inf->v1.i;
                ma->inf[out->n_info].v1.f = inf->v1.f;
                ma->inf[out->n_info].vptr_off  = inf->vptr_off;
                ma->inf[out->n_info].vptr_len  = inf->vptr_len;
                ma->inf[out->n_info].vptr_free = inf->vptr_free;
                if ( (args->output_type & FT_BCF) && id!=bcf_hdr_id2int(hdr, BCF_DT_ID, key) )
                {
                    // The existing packed info cannot be reused. Change the id.
                    // Although quite hacky, it's faster than anything else given
                    // the data structures
                    bcf_info_set_id(out, &ma->inf[out->n_info], id, &args->tmps);
                }
                out->n_info++;
                kitr = kh_put(strdict, tmph, key, &ret);
                kh_val(tmph,kitr) = -(out->n_info-1);   // arbitrary negative value
            }
        }
    }
    out->d.info = ma->inf;
    out->d.m_info = ma->minf;
    for (i=0; i<args->nrules; i++)
        args->rules[i].merger(args->out_hdr, out, &args->rules[i]);
    for (i=0; i<ma->nAGR_info; i++)
    {
        AGR_info_t *agr = &ma->AGR_info[i];
        bcf_update_info(out_hdr,out,agr->hdr_tag,agr->buf,agr->nvals,agr->type);
    }
}

void update_AN_AC(bcf_hdr_t *hdr, bcf1_t *line)
{
    int32_t an = 0, *tmp = (int32_t*) malloc(sizeof(int)*line->n_allele);
    int ret = bcf_calc_ac(hdr, line, tmp, BCF_UN_FMT);
    if ( ret>0 )
    {
        int i;
        for (i=0; i<line->n_allele; i++) an += tmp[i];
        bcf_update_info_int32(hdr, line, "AN", &an, 1);
        bcf_update_info_int32(hdr, line, "AC", tmp+1, line->n_allele-1);
    }
    free(tmp);
}

void merge_GT(args_t *args, bcf_fmt_t **fmt_map, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = bcf_hdr_nsamples(out_hdr);

    int nsize = 0, msize = sizeof(int32_t);
    for (i=0; i<files->nreaders; i++)
    {
        if ( !fmt_map[i] ) continue;
        if ( fmt_map[i]->n > nsize ) nsize = fmt_map[i]->n;
    }

    if ( ma->ntmp_arr < nsamples*nsize*msize )
    {
        ma->ntmp_arr = nsamples*nsize*msize;
        ma->tmp_arr  = realloc(ma->tmp_arr, ma->ntmp_arr);
    }
    memset(ma->smpl_ploidy,0,nsamples*sizeof(int));

    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;
        bcf_fmt_t *fmt_ori = fmt_map[i];
        int32_t *tmp  = (int32_t *) ma->tmp_arr + ismpl*nsize;

        int j, k;
        if ( !fmt_ori )
        {
            // missing values: assume maximum ploidy
            for (j=0; j<bcf_hdr_nsamples(hdr); j++)
            {
                for (k=0; k<nsize; k++) { tmp[k] = 0; ma->smpl_ploidy[ismpl+j]++; }
                tmp += nsize;
            }
            ismpl += bcf_hdr_nsamples(hdr);
            continue;
        }

        #define BRANCH(type_t, missing, vector_end) { \
            type_t *p_ori  = (type_t*) fmt_ori->p; \
            if ( !ma->d[i][0].als_differ ) \
            { \
                /* the allele numbering is unchanged */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    for (k=0; k<fmt_ori->n; k++) \
                    { \
                        if ( p_ori[k]==vector_end ) break; /* smaller ploidy */ \
                        ma->smpl_ploidy[ismpl+j]++; \
                        if ( p_ori[k]==missing ) tmp[k] = 0; /* missing allele */ \
                        else tmp[k] = p_ori[k]; \
                    } \
                    for (; k<nsize; k++) tmp[k] = bcf_int32_vector_end; \
                    tmp += nsize; \
                    p_ori += fmt_ori->n; \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            /* allele numbering needs to be changed */ \
            for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
            { \
                for (k=0; k<fmt_ori->n; k++) \
                { \
                    if ( p_ori[k]==vector_end ) break; /* smaller ploidy */ \
                    ma->smpl_ploidy[ismpl+j]++; \
                    if ( !(p_ori[k]>>1) || p_ori[k]==missing ) tmp[k] = 0; /* missing allele */ \
                    else \
                    { \
                        int al = (p_ori[k]>>1) - 1; \
                        al = al<=0 ? al + 1 : ma->d[i][0].map[al] + 1; \
                        tmp[k] = (al << 1) | ((p_ori[k])&1); \
                    } \
                } \
                for (; k<nsize; k++) tmp[k] = bcf_int32_vector_end; \
                tmp += nsize; \
                p_ori += fmt_ori->n; \
            } \
            ismpl += bcf_hdr_nsamples(hdr); \
        }
        switch (fmt_ori->type)
        {
            case BCF_BT_INT8: BRANCH(int8_t,   bcf_int8_missing,  bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
            default: error("Unexpected case: %d\n", fmt_ori->type);
        }
        #undef BRANCH
    }
    bcf_update_format_int32(out_hdr, out, "GT", (int32_t*)ma->tmp_arr, nsamples*nsize);
}

void merge_format_field(args_t *args, bcf_fmt_t **fmt_map, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = bcf_hdr_nsamples(out_hdr);

    const char *key = NULL;
    int nsize = 0, length = BCF_VL_FIXED, type = -1;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        if ( !fmt_map[i] ) continue;
        if ( !key ) key = files->readers[i].header->id[BCF_DT_ID][fmt_map[i]->id].key;
        type = fmt_map[i]->type;
        if ( IS_VL_G(files->readers[i].header, fmt_map[i]->id) )
        {
            length = BCF_VL_G;
            nsize = out->n_allele*(out->n_allele + 1)/2;
            int nals_ori = files->readers[i].buffer[0]->n_allele;
            if ( fmt_map[i]->n != nals_ori*(nals_ori+1)/2 && fmt_map[i]->n != nals_ori )
                error("Incorrect number of %s fields at %s:%d, cannot merge.\n", key,bcf_seqname(args->out_hdr,out),out->pos+1);
            break;
        }
        if ( IS_VL_A(files->readers[i].header, fmt_map[i]->id) )
        {
            length = BCF_VL_A;
            nsize = out->n_allele - 1;
            int nals_ori = files->readers[i].buffer[0]->n_allele;
            if ( fmt_map[i]->n != nals_ori-1 )
                error("Incorrect number of %s fields at %s:%d, cannot merge.\n", key,bcf_seqname(args->out_hdr,out),out->pos+1);
            break;
        }
        if ( fmt_map[i]->n > nsize ) nsize = fmt_map[i]->n;
    }

    int msize = sizeof(float)>sizeof(int32_t) ? sizeof(float) : sizeof(int32_t);
    if ( ma->ntmp_arr < nsamples*nsize*msize )
    {
        ma->ntmp_arr = nsamples*nsize*msize;
        ma->tmp_arr  = realloc(ma->tmp_arr, ma->ntmp_arr);
    }

    // Fill the temp array for all samples by collecting values from all files
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;
        bcf_fmt_t *fmt_ori = fmt_map[i];
        if ( fmt_ori ) type = fmt_ori->type;

        // set the values
        #define BRANCH(tgt_type_t, src_type_t, src_is_missing, src_is_vector_end, tgt_set_missing, tgt_set_vector_end) { \
            int j, l, k; \
            tgt_type_t *tgt = (tgt_type_t *) ma->tmp_arr + ismpl*nsize; \
            if ( !fmt_ori ) \
            { \
                /* the field is not present in this file, set missing values */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt_set_missing; tgt++; for (l=1; l<nsize; l++) { tgt_set_vector_end; tgt++; } \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            assert( ma->has_line[i] ); \
            bcf1_t *line    = reader->buffer[0]; \
            src_type_t *src = (src_type_t*) fmt_ori->p; \
            if ( (length!=BCF_VL_G && length!=BCF_VL_A) || (line->n_allele==out->n_allele && !ma->d[i][0].als_differ) ) \
            { \
                /* alleles unchanged, copy over */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    for (l=0; l<fmt_ori->n; l++) \
                    { \
                        if ( src_is_vector_end ) break; \
                        else if ( src_is_missing ) tgt_set_missing; \
                        else *tgt = *src; \
                        tgt++; src++; \
                    } \
                    for (k=l; k<nsize; k++) { tgt_set_vector_end; tgt++; } \
                    src += fmt_ori->n - l; \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            /* allele numbering needs to be changed */ \
            if ( length==BCF_VL_G ) \
            { \
                /* Number=G tags */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize; \
                    int ngsize = ma->smpl_ploidy[ismpl+j]==1 ? out->n_allele : out->n_allele*(out->n_allele + 1)/2; \
                    for (l=0; l<ngsize; l++) { tgt_set_missing; tgt++; } \
                    for (; l<nsize; l++) { tgt_set_vector_end; tgt++; } \
                    int iori,jori, inew,jnew; \
                    for (iori=0; iori<line->n_allele; iori++) \
                    { \
                        inew = ma->d[i][0].map[iori]; \
                        for (jori=0; jori<=iori; jori++) \
                        { \
                            jnew = ma->d[i][0].map[jori]; \
                            int kori = iori*(iori+1)/2 + jori; \
                            int knew = inew>jnew ? inew*(inew+1)/2 + jnew : jnew*(jnew+1)/2 + inew; \
                            src = (src_type_t*) fmt_ori->p + j*fmt_ori->n + kori; \
                            tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + knew; \
                            if ( src_is_vector_end ) \
                            { \
                                iori = line->n_allele; \
                                break; \
                            } \
                            if ( src_is_missing ) tgt_set_missing; \
                            else *tgt = *src; \
                        } \
                    } \
                } \
            } \
            else \
            { \
                /* Number=A tags */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize; \
                    for (l=0; l<nsize; l++) { tgt_set_missing; tgt++; } \
                    int iori,inew; \
                    for (iori=1; iori<line->n_allele; iori++) \
                    { \
                        inew = ma->d[i][0].map[iori] - 1; \
                        tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + inew; \
                        if ( src_is_vector_end ) break; \
                        if ( src_is_missing ) tgt_set_missing; \
                        else *tgt = *src; \
                        src++; \
                    } \
                } \
            } \
            ismpl += bcf_hdr_nsamples(hdr); \
        }
        switch (type)
        {
            case BCF_BT_INT8:  BRANCH(int32_t,  int8_t, *src==bcf_int8_missing,  *src==bcf_int8_vector_end,  *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_INT16: BRANCH(int32_t, int16_t, *src==bcf_int16_missing, *src==bcf_int16_vector_end, *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_FLOAT: BRANCH(float, float, bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), bcf_float_set_missing(*tgt), bcf_float_set_vector_end(*tgt)); break;
            case BCF_BT_CHAR:  BRANCH(uint8_t, uint8_t, *src==bcf_str_missing, *src==bcf_str_vector_end, *tgt=bcf_str_missing, *tgt=bcf_str_vector_end); break;
            default: error("Unexpected case: %d, %s\n", type, key);
        }
        #undef BRANCH
    }
    if ( type==BCF_BT_FLOAT )
        bcf_update_format_float(out_hdr, out, key, (float*)ma->tmp_arr, nsamples*nsize);
    else if ( type==BCF_BT_CHAR )
        bcf_update_format_char(out_hdr, out, key, (float*)ma->tmp_arr, nsamples*nsize);
    else
        bcf_update_format_int32(out_hdr, out, key, (int32_t*)ma->tmp_arr, nsamples*nsize);
}

void merge_format(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    if ( !ma->nfmt_map )
    {
        ma->nfmt_map = 2;
        ma->fmt_map  = (bcf_fmt_t**) calloc(ma->nfmt_map*files->nreaders, sizeof(bcf_fmt_t*));
    }
    else
        memset(ma->fmt_map, 0, ma->nfmt_map*files->nreaders*sizeof(bcf_fmt_t**));

    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
    int i, j, ret, has_GT = 0, max_ifmt = 0; // max fmt index
    for (i=0; i<files->nreaders; i++)
    {
        if ( !ma->has_line[i] ) continue;
        bcf_sr_t *reader = &files->readers[i];
        bcf1_t *line = reader->buffer[0];
        bcf_hdr_t *hdr = reader->header;
        for (j=0; j<line->n_fmt; j++)
        {
            // Wat this tag already seen?
            bcf_fmt_t *fmt = &line->d.fmt[j];
            const char *key = hdr->id[BCF_DT_ID][fmt->id].key;
            kitr = kh_get(strdict, tmph, key);

            int ifmt;
            if ( kitr != kh_end(tmph) )
                ifmt = kh_value(tmph, kitr);    // seen
            else
            {
                // new FORMAT tag
                if ( key[0]=='G' && key[1]=='T' && key[2]==0 ) { has_GT = 1; ifmt = 0; }
                else
                {
                    ifmt = ++max_ifmt;
                    if ( max_ifmt >= ma->nfmt_map )
                    {
                        ma->fmt_map = (bcf_fmt_t**) realloc(ma->fmt_map, sizeof(bcf_fmt_t*)*(max_ifmt+1)*files->nreaders);
                        memset(ma->fmt_map+ma->nfmt_map*files->nreaders, 0, (max_ifmt-ma->nfmt_map+1)*files->nreaders*sizeof(bcf_fmt_t*));
                        ma->nfmt_map = max_ifmt+1;
                    }
                }
                kitr = kh_put(strdict, tmph, key, &ret);
                kh_value(tmph, kitr) = ifmt;
            }
            ma->fmt_map[ifmt*files->nreaders+i] = fmt;
        }
        // Check if the allele numbering must be changed
        for (j=1; j<reader->buffer[0]->n_allele; j++)
            if ( ma->d[i][0].map[j]!=j ) break;
        ma->d[i][0].als_differ = j==reader->buffer[0]->n_allele ? 0 : 1;
    }

    out->n_sample = bcf_hdr_nsamples(out_hdr);
    if ( has_GT )
        merge_GT(args, ma->fmt_map, out);
    update_AN_AC(out_hdr, out);

    if ( out->d.info!=ma->inf )
    {
        // hacky, we rely on htslib internals: bcf_update_info() reallocated the info
        ma->inf  = out->d.info;
        ma->minf = out->d.m_info;
    }

    for (i=1; i<=max_ifmt; i++)
        merge_format_field(args, &ma->fmt_map[i*files->nreaders], out);
    out->d.indiv_dirty = 1;
}

// The core merging function, one or none line from each reader
void merge_line(args_t *args)
{
    bcf1_t *out = args->out_line;
    bcf_clear1(out);
    out->unpacked = BCF_UN_ALL;

    merge_chrom2qual(args, out);
    merge_filter(args, out);
    merge_info(args, out);
    merge_format(args, out);

    bcf_write1(args->out_fh, args->out_hdr, out);
}


void debug_buffers(FILE *fp, bcf_srs_t *files);
void debug_buffer(FILE *fp, bcf_sr_t *reader);

#define SWAP(type_t,a,b) { type_t tmp = (a); (a) = (b); (b) = tmp; }

// Clean the reader's buffer to and make it ready for the next next_line() call.
// Moves finished records (SKIP_DONE flag set) at the end of the buffer and put
// the rest to the beggining. Then shorten the buffer so that the last element
// points to the last unfinished record. There are two special cases: the last
// line of the buffer typically has a different position and must stay at the
// end; next, the first record of the buffer must be one of those already
// printed, as it will be discarded by next_line().
//
void shake_buffer(maux_t *maux, int ir, int pos)
{
    bcf_sr_t *reader = &maux->files->readers[ir];
    maux1_t *m = maux->d[ir];

    if ( !reader->buffer ) return;

    int i;
    // FILE *fp = stdout;
    // fprintf(fp,"<going to shake> nbuf=%d\t", reader->nbuffer); for (i=0; i<reader->nbuffer; i++) fprintf(fp," %d", skip[i]); fprintf(fp,"\n");
    // debug_buffer(fp,reader);
    // fprintf(fp,"--\n");

    int a = 1, b = reader->nbuffer;
    if ( reader->buffer[b]->pos != pos ) b--;   // move the last line separately afterwards

    while ( a<b )
    {
        if ( !(m[a].skip&SKIP_DONE) ) { a++; continue; }
        if ( m[b].skip&SKIP_DONE ) { b--; continue; }
        SWAP(bcf1_t*, reader->buffer[a], reader->buffer[b]);
        SWAP(maux1_t, m[a], m[b]);
        a++;
        b--;
    }

    // position $a to the after the first unfinished record
    while ( a<=reader->nbuffer && !(m[a].skip&SKIP_DONE) ) a++;

    if ( a<reader->nbuffer )
    {
        // there is a gap between the unfinished lines at the beggining and the
        // last line. The last line must be brought forward to fill the gap
        if ( reader->buffer[reader->nbuffer]->pos != pos )
        {
            SWAP(bcf1_t*, reader->buffer[a], reader->buffer[reader->nbuffer]);
            SWAP(maux1_t, m[a], m[reader->nbuffer]);
            reader->nbuffer = a;
        }
    }

    if ( !(m[0].skip&SKIP_DONE) && reader->buffer[0]->pos==pos )
    {
        // the first record is unfinished, replace it with an empty line
        // from the end of the buffer or else next_line will remove it
        if ( reader->nbuffer + 1 >= maux->nbuf[ir] )
        {
            reader->nbuffer++;
            maux_expand1(maux, ir);
            reader->nbuffer--;
            m = maux->d[ir];
        }
        if ( reader->nbuffer+1 >= reader->mbuffer )
            error("Uh, did not expect this: %d vs %d\n", reader->nbuffer,reader->mbuffer);

        if ( reader->buffer[reader->nbuffer]->pos!=pos )
        {
            // 4way swap
            bcf1_t *tmp = reader->buffer[0];
            reader->buffer[0] = reader->buffer[reader->nbuffer+1];
            reader->buffer[reader->nbuffer+1] = reader->buffer[reader->nbuffer];
            reader->buffer[reader->nbuffer] = tmp;
            m[reader->nbuffer].skip   = m[0].skip;
            m[reader->nbuffer+1].skip = SKIP_DIFF;
            reader->nbuffer++;
        }
        else
        {
            SWAP(bcf1_t*, reader->buffer[0], reader->buffer[reader->nbuffer+1]);
            SWAP(maux1_t, m[0], m[reader->nbuffer+1]);
        }
    }

    // debug_buffer(fp,reader);
    // fprintf(fp,"<shaken>\t"); for (i=0; i<reader->nbuffer; i++) fprintf(fp," %d", skip[i]);
    // fprintf(fp,"\n\n");

    // set position of finished buffer[0] line to -1, otherwise swapping may
    // bring it back after next_line()
    reader->buffer[0]->pos = -1;

    // trim the buffer, remove finished lines from the end
    i = reader->nbuffer;
    while ( i>=1 && m[i--].skip&SKIP_DONE )
        reader->nbuffer--;
}

void debug_maux(args_t *args, int pos, int var_type)
{
    bcf_srs_t *files = args->files;
    maux_t *maux = args->maux;
    int j,k,l;

    fprintf(stderr,"Alleles to merge at %d\n", pos+1);
    for (j=0; j<files->nreaders; j++)
    {
        bcf_sr_t *reader = &files->readers[j];
        fprintf(stderr," reader %d: ", j);
        for (k=0; k<=reader->nbuffer; k++)
        {
            if ( maux->d[j][k].skip==SKIP_DONE ) continue;
            bcf1_t *line = reader->buffer[k];
            if ( line->pos!=pos ) continue;
            fprintf(stderr,"\t");
            if ( maux->d[j][k].skip ) fprintf(stderr,"[");  // this record will not be merged in this round
            for (l=0; l<line->n_allele; l++)
                fprintf(stderr,"%s%s", l==0?"":",", line->d.allele[l]);
            if ( maux->d[j][k].skip ) fprintf(stderr,"]");
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr," counts: ");
    for (j=0; j<maux->nals; j++) fprintf(stderr,"%s   %dx %s", j==0?"":",",maux->cnt[j], maux->als[j]); fprintf(stderr,"\n");
    for (j=0; j<files->nreaders; j++)
    {
        bcf_sr_t *reader = &files->readers[j];
        fprintf(stderr," out %d: ", j);
        for (k=0; k<=reader->nbuffer; k++)
        {
            if ( maux->d[j][k].skip==SKIP_DONE ) continue;
            bcf1_t *line = reader->buffer[k];
            if ( line->pos!=pos ) continue;
            if ( maux->d[j][k].skip ) continue;
            fprintf(stderr,"\t");
            for (l=0; l<line->n_allele; l++)
                fprintf(stderr,"%s%s", l==0?"":",", maux->als[maux->d[j][k].map[l]]);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
}

// Determine which line should be merged from which reader: go through all
// readers and all buffered lines, expand REF,ALT and try to match lines with
// the same ALTs. A step towards output independent on input ordering of the
// lines.
void merge_buffer(args_t *args)
{
    bcf_srs_t *files = args->files;
    int i, pos = -1, var_type = 0;
    char *id = NULL;
    maux_t *maux = args->maux;
    maux_reset(maux);

    // set the current position
    for (i=0; i<files->nreaders; i++)
    {
        if ( bcf_sr_has_line(files,i) )
        {
            bcf1_t *line = bcf_sr_get_line(files,i);
            pos = line->pos;
            var_type = bcf_get_variant_types(line);
            id = line->d.id;
            break;
        }
    }

    // In this loop we select from each reader compatible candidate lines.
    // (i.e. SNPs or indels). Go through all files and all lines at this
    // position and normalize relevant alleles.
    // REF-only sites may be associated with both SNPs and indels.
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        if ( !reader->buffer ) continue;
        int j, k;
        for (j=0; j<=reader->nbuffer; j++)
        {
            bcf1_t *line = reader->buffer[j];
            int line_type = bcf_get_variant_types(line);

            // select relevant lines
            maux->d[i][j].skip = SKIP_DIFF;
            if ( pos!=line->pos )
            {
                if ( j==0 ) maux->d[i][j].skip |= SKIP_DONE; // left from previous run, force to ignore
                continue;
            }
            if ( args->merge_by_id )
            {
                if ( strcmp(id,line->d.id) ) continue;
            }
            else
            {
                if ( args->collapse==COLLAPSE_NONE && maux->nals )
                {
                    // All alleles of the tested record must be present in the
                    // selected maux record plus variant types must be the same
                    if ( var_type!=line->d.var_type ) continue;
                    if ( vcmp_set_ref(args->vcmp,maux->als[0],line->d.allele[0]) < 0 ) continue;   // refs not compatible
                    for (k=1; k<line->n_allele; k++)
                    {
                        if ( vcmp_find_allele(args->vcmp,maux->als+1,maux->nals-1,line->d.allele[k])>=0 ) break;
                    }
                    if ( k==line->n_allele ) continue;  // no matching allele
                }
                if ( var_type&VCF_SNP && !(line_type&VCF_SNP) && !(args->collapse&COLLAPSE_ANY) && line_type!=VCF_REF ) continue;
                if ( var_type&VCF_INDEL && !(line_type&VCF_INDEL) && !(args->collapse&COLLAPSE_ANY) && line_type!=VCF_REF ) continue;
            }
            maux->d[i][j].skip = 0;

            hts_expand(int, line->n_allele, maux->d[i][j].mmap, maux->d[i][j].map);
            if ( !maux->nals )    // first record, copy the alleles to the output
            {
                maux->nals = line->n_allele;
                hts_expand0(char*, maux->nals, maux->mals, maux->als);
                hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
                for (k=0; k<maux->nals; k++)
                {
                    maux->als[k] = strdup(line->d.allele[k]);
                    maux->d[i][j].map[k] = k;
                    maux->cnt[k] = 1;
                }
                pos = line->pos;
                continue;
            }

            // normalize alleles
            maux->als = merge_alleles(line->d.allele, line->n_allele, maux->d[i][j].map, maux->als, &maux->nals, &maux->mals);
            if ( !maux->als ) error("Failed to merge alleles at %s:%d\n",bcf_seqname(args->out_hdr,line),line->pos+1);
            hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
            for (k=1; k<line->n_allele; k++)
                maux->cnt[ maux->d[i][j].map[k] ]++;    // how many times an allele appears in the files
            maux->cnt[0]++;
        }
    }

    // debug_maux(args, pos, var_type);

    // Select records that have the same alleles; the input ordering of indels
    // must not matter. Multiple VCF lines can be emitted from this loop.
    // We expect only very few alleles and not many records with the same
    // position in the buffers, therefore the nested loops should not slow us
    // much.
    while (1)
    {
        // take the most frequent allele present in multiple files
        int icnt = 0;
        for (i=1; i<maux->nals; i++)
            if ( maux->cnt[i] > maux->cnt[icnt] ) icnt = i;
        if ( maux->cnt[icnt]<0 ) break;

        int nmask = 0;
        for (i=0; i<files->nreaders; i++)
        {
            maux->has_line[i] = 0;

            bcf_sr_t *reader = &files->readers[i];
            if ( !reader->buffer ) continue;

            // find lines with the same allele
            int j;
            for (j=0; j<=reader->nbuffer; j++)
            {
                if ( maux->d[i][j].skip ) continue;
                int k;
                for (k=0; k<reader->buffer[j]->n_allele; k++)
                    if ( icnt==maux->d[i][j].map[k] ) break;
                if ( k<reader->buffer[j]->n_allele ) break;
            }
            if ( j>reader->nbuffer )
            {
                // no matching allele found in this file
                if ( args->collapse==COLLAPSE_NONE ) continue;

                for (j=0; j<=reader->nbuffer; j++)
                {
                    if ( maux->d[i][j].skip ) continue;
                    if ( args->collapse&COLLAPSE_ANY ) break;
                    int line_type = bcf_get_variant_types(reader->buffer[j]);
                    if ( var_type&VCF_SNP && line_type&VCF_SNP && (args->collapse&COLLAPSE_SNPS) ) break;
                    if ( var_type&VCF_INDEL && line_type&VCF_INDEL && (args->collapse&COLLAPSE_INDELS) ) break;
                    if ( line_type==VCF_REF )
                    {
                        if ( var_type&VCF_SNP && (args->collapse&COLLAPSE_SNPS) ) break;
                        if ( var_type&VCF_INDEL && (args->collapse&COLLAPSE_INDELS) ) break;
                    }
                    else if ( var_type==VCF_REF )
                    {
                        if ( line_type&VCF_SNP && (args->collapse&COLLAPSE_SNPS) ) break;
                        if ( line_type&VCF_INDEL && (args->collapse&COLLAPSE_INDELS) ) break;
                    }
                }
            }
            if ( j<=reader->nbuffer )
            {
                // found a suitable line for merging, place it at the beggining
                if ( j>0 )
                {
                    SWAP(bcf1_t*, reader->buffer[0], reader->buffer[j]);
                    SWAP(maux1_t, maux->d[i][0], maux->d[i][j]);
                }
                // mark as finished so that it's ignored next time
                maux->d[i][0].skip |= SKIP_DONE;
                maux->has_line[i] = 1;
                nmask++;
            }
        }
        if ( !nmask ) break;    // done, no more lines suitable for merging found
        merge_line(args);       // merge and output the line
        maux->cnt[icnt] = -1;   // do not pick this allele again, mark it as finished
    }

    // clean the alleles
    for (i=0; i<maux->nals; i++)
    {
        free(maux->als[i]);
        maux->als[i] = 0;
    }
    maux->nals = 0;

    // get the buffers ready for the next next_line() call
    for (i=0; i<files->nreaders; i++)
        shake_buffer(maux, i, pos);
}

void bcf_hdr_append_version(bcf_hdr_t *hdr, int argc, char **argv, const char *cmd)
{
    kstring_t str = {0,0,0};
    ksprintf(&str,"##%sVersion=%s+htslib-%s\n", cmd, bcftools_version(), hts_version());
    bcf_hdr_append(hdr,str.s);

    str.l = 0;
    ksprintf(&str,"##%sCommand=%s", cmd, argv[0]);
    int i;
    for (i=1; i<argc; i++)
    {
        if ( strchr(argv[i],' ') )
            ksprintf(&str, " '%s'", argv[i]);
        else
            ksprintf(&str, " %s", argv[i]);
    }
    kputc('\n', &str);
    bcf_hdr_append(hdr,str.s);
    free(str.s);

    bcf_hdr_sync(hdr);
}

void merge_vcf(args_t *args)
{
    args->out_fh  = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    args->out_hdr = bcf_hdr_init("w");

    if ( args->header_fname )
    {
        if ( bcf_hdr_set(args->out_hdr,args->header_fname) ) error("Could not read/parse the header: %s\n", args->header_fname);
    }
    else
    {
        int i;
        for (i=0; i<args->files->nreaders; i++)
        {
            char buf[10]; snprintf(buf,10,"%d",i+1);
            bcf_hdr_merge(args->out_hdr, args->files->readers[i].header,buf,args->force_samples);
        }
        bcf_hdr_append_version(args->out_hdr, args->argc, args->argv, "bcftools_merge");
        bcf_hdr_sync(args->out_hdr);
    }
    info_rules_init(args);

    bcf_hdr_set_version(args->out_hdr, bcf_hdr_get_version(args->files->readers[0].header));
    bcf_hdr_write(args->out_fh, args->out_hdr);
    if ( args->header_only )
    {
        bcf_hdr_destroy(args->out_hdr);
        hts_close(args->out_fh);
        return;
    }

    if ( args->collapse==COLLAPSE_NONE ) args->vcmp = vcmp_init();
    args->maux = maux_init(args->files);
    args->out_line = bcf_init1();
    args->tmph = kh_init(strdict);
    int ret;
    while ( (ret=bcf_sr_next_line(args->files)) )
    {
        merge_buffer(args);
    }
    info_rules_destroy(args);
    maux_destroy(args->maux);
    bcf_hdr_destroy(args->out_hdr);
    hts_close(args->out_fh);
    bcf_destroy1(args->out_line);
    kh_destroy(strdict, args->tmph);
    if ( args->tmps.m ) free(args->tmps.s);
    if ( args->vcmp ) vcmp_destroy(args->vcmp);
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.\n");
    fprintf(stderr, "         Compatible records are combined into one according to the -m option.\n");
    fprintf(stderr, "Usage:   bcftools merge [options] <A.vcf.gz> <B.vcf.gz> [...]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "        --force-samples                resolve duplicate sample names\n");
    fprintf(stderr, "        --print-header                 print only the merged header and exit\n");
    fprintf(stderr, "        --use-header <file>            use the provided header\n");
    fprintf(stderr, "    -f, --apply-filters <list>         require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -i, --info-rules <tag:method,..>   rules for merging INFO fields (method is one of sum,avg,min,max,join) or \"-\" to turn off the default [DP:sum,DP4:sum]\n");
    fprintf(stderr, "    -l, --file-list <file>             read file names from the file\n");
    fprintf(stderr, "    -m, --merge <string>               merge sites with differing alleles for <snps|indels|both|all|none|id>, see man page for details [both]\n");
    fprintf(stderr, "    -o, --output <file>                write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type <b|u|z|v>        'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <region>             restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>          restrict to regions listed in a file\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfmerge(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->collapse = COLLAPSE_BOTH;
    int regions_is_file = 0;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"merge",1,0,'m'},
        {"file-list",1,0,'l'},
        {"apply-filters",1,0,'f'},
        {"use-header",1,0,1},
        {"print-header",0,0,2},
        {"force-samples",0,0,3},
        {"output",1,0,'o'},
        {"output-type",1,0,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"info-rules",1,0,'i'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "hm:f:r:R:o:O:i:l:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'l': args->file_list = optarg; break;
            case 'i': args->info_rules = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 'm':
                args->collapse = COLLAPSE_NONE;
                if ( !strcmp(optarg,"snps") ) args->collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) args->collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) args->collapse |= COLLAPSE_BOTH;
                else if ( !strcmp(optarg,"any") ) args->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"all") ) args->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"none") ) args->collapse = COLLAPSE_NONE;
                else if ( !strcmp(optarg,"id") ) { args->collapse = COLLAPSE_NONE; args->merge_by_id = 1; }
                else error("The -m type \"%s\" is not recognised.\n", optarg);
                break;
            case 'f': args->files->apply_filters = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case  1 : args->header_fname = optarg; break;
            case  2 : args->header_only = 1; break;
            case  3 : args->force_samples = 1; break;
            case 'h':
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( argc==optind && !args->file_list ) usage();
    if ( argc-optind<2 && !args->file_list ) usage();

    args->files->require_index = 1;
    if ( args->regions_list && bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
        error("Failed to read the regions: %s\n", args->regions_list);

    while (optind<argc)
    {
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open: %s\n", argv[optind]);
        optind++;
    }
    if ( args->file_list )
    {
        int nfiles, i;
        char **files = hts_readlines(args->file_list, &nfiles);
        if ( !files ) error("Failed to read from %s\n", args->file_list);
        for (i=0;i<nfiles; i++)
            if ( !bcf_sr_add_reader(args->files, files[i]) ) error("Failed to open: %s\n", files[i]);
        for (i=0; i<nfiles; i++) free(files[i]);
        free(files);
    }
    merge_vcf(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

