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
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/regidx.h>
#include <htslib/khash_str2int.h>
#include <htslib/kseq.h>

typedef struct
{
    int sex, ploidy;
}
sex_ploidy_t;

static kstring_t tmp_str = {0,0,0};
static bcf_hdr_t *in_hdr = NULL, *out_hdr = NULL;
static regidx_t *idx = NULL;
static int *sample2sex = NULL;
static int n_sample = 0, nsex = 0, *sex2ploidy = NULL;
static int32_t ngt_arr = 0, *gt_arr = NULL, *gt_arr2 = NULL, ngt_arr2 = 0; 

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

const char *about(void)
{
    return "Fix ploidy.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Fix ploidy\n"
        "Usage: bcftools +fixploidy [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -p, --ploidy <file>   space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY\n"
        "   -s, --sex <files>     list of samples, \"NAME SEX\"\n"
        "   -t, --tags <list>     VCF tags to fix [GT]\n"
        "\n"
        "Example:\n"
        "   # Default ploidy, if -p not given. Unlisted regions have ploidy 2\n"
        "   X 1 60000 M 1\n"
        "   X 2699521 154931043 M 1\n"
        "   Y 1 59373566 M 1\n"
        "   Y 1 59373566 F 0\n"
        "   MT 1 16569 M 1\n"
        "   MT 1 16569 F 1\n"
        "   \n"
        "   # Example of -s file, sex of unlisted samples is \"F\"\n"
        "   sampleName1 M\n"
        "   \n"
        "   bcftools +fixploidy in.vcf -- -s samples.txt\n"
        "\n";
}

int ploidy_parse(const char *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *sex2id)
{
    // Fill CHR,FROM,TO
    int i, ret = regidx_parse_tab(line,chr_beg,chr_end,reg,NULL,NULL);
    if ( ret!=0 ) return ret;

    // Skip the fields already parsed by regidx_parse_tab
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    for (i=0; i<3; i++)
    {
        while ( *ss && !isspace(*ss) ) ss++;
        if ( !*ss ) return -2;  // wrong number of fields
        while ( *ss && isspace(*ss) ) ss++;
    }
    if ( !*ss ) return -2;

    // Parse the payload
    char *se = ss;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se || se==ss ) error("Could not parse: %s\n", line);
    tmp_str.l = 0;
    kputsn(ss,se-ss,&tmp_str);

    sex_ploidy_t *sp = (sex_ploidy_t*) payload;
    if ( khash_str2int_get(sex2id, tmp_str.s, &sp->sex) != 0 )
        sp->sex = khash_str2int_inc(sex2id,strdup(tmp_str.s));

    ss = se;
    while ( *se && isspace(*se) ) se++;
    if ( !*se ) error("Could not parse: %s\n", line);
    sp->ploidy = strtol(ss,&se,10);
    if ( ss==se ) error("Could not parse: %s\n", line);

    return 0;
}

void set_samples(char *fname, bcf_hdr_t *hdr, void *sex2id, int *sample2sex)
{
    kstring_t tmp = {0,0,0};

    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);
    while ( hts_getline(fp, KS_SEP_LINE, &tmp) > 0 )
    {
        char *ss = tmp.s;
        while ( *ss && isspace(*ss) ) ss++;
        if ( !*ss ) error("Could not parse: %s\n", tmp.s);
        if ( *ss=='#' ) continue;
        char *se = ss;
        while ( *se && !isspace(*se) ) se++;
        char x = *se; *se = 0;

        int ismpl = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, ss);
        if ( ismpl < 0 ) { fprintf(stderr,"Warning: No such sample in the VCF: %s\n",ss); continue; }

        *se = x;
        ss = se+1;
        while ( *ss && isspace(*ss) ) ss++;
        if ( !*ss )  error("Could not parse: %s\n", tmp.s);
        se = ss;
        while ( *se && !isspace(*se) ) se++;
        if ( se==ss ) error("Could not parse: %s\n", tmp.s);

        int sex_id;
        if ( khash_str2int_get(sex2id, ss, &sex_id)!=0 )
            sex_id = khash_str2int_inc(sex2id, ss);
        sample2sex[ismpl] = sex_id;
    }
    if ( hts_close(fp) ) error("Close failed: %s\n", fname);
    free(tmp.s);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    int c;
    char *tags_str = "GT";
    char *ploidy_fname = NULL, *sex_fname = NULL;

    static struct option loptions[] =
    {
        {"ploidy",1,0,'p'},
        {"sex",1,0,'s'},
        {"tags",1,0,'t'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "?ht:s:p:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'p': ploidy_fname = optarg; break;
            case 's': sex_fname = optarg; break;
            case 't': tags_str = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( strcasecmp("GT",tags_str) ) error("Only -t GT is currently supported, sorry\n");

    void *sex2id = khash_str2int_init();
    if ( ploidy_fname )
        idx = regidx_init(ploidy_fname,ploidy_parse,NULL,sizeof(sex_ploidy_t),sex2id);
    else
    {
        idx = regidx_init(NULL,ploidy_parse,NULL,sizeof(sex_ploidy_t),sex2id);
        regidx_insert(idx,"X 1 60000 M 1");
        regidx_insert(idx,"X 2699521 154931043 M 1");
        regidx_insert(idx,"Y 1 59373566 M 1");
        regidx_insert(idx,"Y 1 59373566 F 0");
        regidx_insert(idx,"MT 1 16569 M 1");
        regidx_insert(idx,"MT 1 16569 F 1");
        regidx_insert(idx,NULL);
    }

    int dflt_sex_id = 0;
    if ( khash_str2int_get(sex2id,"F",&dflt_sex_id)!=0 )
        dflt_sex_id = khash_str2int_inc(sex2id, strdup("F"));

    int i;

    nsex = khash_str2int_size(sex2id);
    sex2ploidy = (int*) malloc(sizeof(int)*nsex);

    in_hdr  = in;
    out_hdr = out;
    n_sample = bcf_hdr_nsamples(in);
    sample2sex = (int*) calloc(n_sample,sizeof(int));
    for (i=0; i<n_sample; i++) sample2sex[i] = dflt_sex_id;
    if ( sex_fname ) set_samples(sex_fname, in, sex2id, sample2sex);

    khash_str2int_destroy_free(sex2id);
    return 0;
}


bcf1_t *process(bcf1_t *rec)
{
    int i,j, max_ploidy = -1, min_ploidy = 99;
    for (i=0; i<nsex; i++) sex2ploidy[i] = 2;

    regitr_t itr;
    if ( regidx_overlap(idx, (char*)bcf_seqname(in_hdr,rec), rec->pos,rec->pos, &itr) )
    {
        while ( REGITR_OVERLAP(itr,rec->pos,rec->pos) )
        {
            int sex    = REGITR_PAYLOAD(itr,sex_ploidy_t).sex;
            int ploidy = REGITR_PAYLOAD(itr,sex_ploidy_t).ploidy;
            if ( ploidy!=2 ) 
            {
                sex2ploidy[ sex ] = ploidy;
                if ( min_ploidy > ploidy ) min_ploidy = ploidy;
                if ( max_ploidy < ploidy ) max_ploidy = ploidy;
            }
            itr.i++;
        }
        if ( max_ploidy==-1 ) max_ploidy = 2;
        if ( min_ploidy==99 ) min_ploidy = 2;
    }
    else
        min_ploidy = max_ploidy = 2;

    int ngts = bcf_get_genotypes(in_hdr, rec, &gt_arr, &ngt_arr);
    if ( ngts % n_sample )
        error("Error at %s:%d: wrong number of GT fields\n",bcf_seqname(in_hdr,rec),rec->pos+1);

    ngts /= n_sample;
    if ( ngts < max_ploidy )
    {
        hts_expand(int32_t,max_ploidy*n_sample,ngt_arr2,gt_arr2);
        for (i=0; i<n_sample; i++)
        {
            int ploidy = sex2ploidy[ sample2sex[i] ];
            int32_t *src = &gt_arr[i*ngts];
            int32_t *dst = &gt_arr2[i*max_ploidy];
            j = 0;
            if ( !ploidy ) { dst[j] = bcf_gt_missing; j++; }
            else
                while ( j<ngts && j<ploidy && src[j]!=bcf_int32_vector_end ) { dst[j] = src[j]; j++; }
            assert( j );
            while ( j<ploidy ) { dst[j] = dst[j-1]; j++; } // expand "." to "./." and "0" to "0/0"
            while ( j<max_ploidy ) { dst[j] = bcf_int32_vector_end; j++; }
        }
        if ( bcf_update_genotypes(out_hdr,rec,gt_arr2,n_sample*max_ploidy) )
            error("Could not update GT field at %s:%d\n", bcf_seqname(in_hdr,rec),rec->pos+1);
    }
    else if ( ngts!=1 || max_ploidy!=1 )
    {
        for (i=0; i<n_sample; i++)
        {
            int ploidy = sex2ploidy[ sample2sex[i] ];
            int32_t *gts = &gt_arr[i*ngts];
            j = 0;
            if ( !ploidy ) { gts[j] = bcf_gt_missing; j++; }
            else 
                while ( j<ngts && j<ploidy && gts[j]!=bcf_int32_vector_end ) j++;
            assert( j );
            while ( j<ploidy ) { gts[j] = gts[j-1]; j++; } // expand "." to "./." and "0" to "0/0"
            while ( j<ngts ) { gts[j] = bcf_int32_vector_end; j++; }
        }
        if ( bcf_update_genotypes(out_hdr,rec,gt_arr,n_sample*ngts) )
            error("Could not update GT field at %s:%d\n", bcf_seqname(in_hdr,rec),rec->pos+1);
    }
    return rec;
}


void destroy(void)
{
    free(gt_arr);
    free(gt_arr2);
    regidx_destroy(idx);
    free(tmp_str.s);
    free(sample2sex);
    free(sex2ploidy);
}


