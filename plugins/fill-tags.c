/* The MIT License

   Copyright (c) 2015 Genome Research Ltd.

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
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include "bcftools.h"
#include "htslib/khash_str2int.h"

#define SET_AN      (1<<0)
#define SET_AC      (1<<1)
#define SET_AC_Hom  (1<<2)
#define SET_AC_Het  (1<<3)
#define SET_AC_Hemi (1<<4)
#define SET_AF      (1<<5)
#define SET_NS      (1<<6)
#define SET_MAF     (1<<7)

static const int POS_NS = 0;
static const int POS_AN = 1;
static const int POS_AC = 2;
static const int POS_AF = 3;
static const int POS_Hom = 4;
static const int POS_Het = 5;
static const int POS_Hemi = 6;
static const int POS_MAF = 7;

typedef struct
{
    int nhom, nhet, nhemi, nac;
}
counts_t;

typedef struct
{
    bcf_hdr_t *in_hdr, *out_hdr;
    int marr, mfarr, mcounts, gt_id, drop_missing, *n_samples, **sample_pos, n_sample_files;
    long tags;
    int32_t *arr;
    float *farr;
    counts_t *counts;
    char ***ids, **sample_files, ***samples;
}
args_t;

static args_t args;

const char *about(void)
{
    return "Set INFO tags MAF, AF, AN, AC, NS, AC_Hom, AC_Het, AC_Hemi.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Set INFO tags MAF, AF, AN, AC, NS, AC_Hom, AC_Het, AC_Hemi.\n"
        "Usage: bcftools +fill-tags [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -d, --drop-missing      do not count half-missing genotypes \"./1\" as hemizygous\n"
        "   -t, --tags LIST         list of output tags. By default, all tags are filled.\n"
        "   -S, --samples [PREFIX,]FILE      Prefix with a file containing a list of samples to use (one id per line) for AC calculations, separated by ','.\n"
        "\n"
        "Example:\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t AN,AC\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -d\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -S GRP_A_,my_a_samples.txt -S GRP_B_,my_b_samples.txt -S NOT_GRP_B,^my_b_samples.txt\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -S my_samples.txt\n"
        "\n";
}

char* concat_prefix(const char * prefix, const char *s2)
{
    char *result = (char *) malloc(strlen(prefix)+strlen(s2)+1);
    strcpy(result, prefix);
    strcat(result, s2);
    return result;
}

int hdr_append(bcf_hdr_t *hdr, const char *str, args_t *args, char * prefix)
{
    const size_t len1 = strlen(prefix);
    const size_t len2 = strlen(str);
    char *result = (char*)malloc(len1+len2+1);
    sprintf(result, str, prefix);
    bcf_hdr_append(args->out_hdr, result);
    free(result);
    return 0;
}

void addHeader(args_t *args, char * prefix)
{
    if ( args->tags&SET_AN ) hdr_append(args->out_hdr, "##INFO=<ID=%sAN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">", args, prefix);
    if ( args->tags&SET_AC ) hdr_append(args->out_hdr, "##INFO=<ID=%sAC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">", args, prefix);
    if ( args->tags&SET_NS ) hdr_append(args->out_hdr, "##INFO=<ID=%sNS,Number=1,Type=Integer,Description=\"Number of samples with data\">", args, prefix);
    if ( args->tags&SET_AC_Hom ) hdr_append(args->out_hdr, "##INFO=<ID=%sAC_Hom,Number=A,Type=Integer,Description=\"Allele counts in homozygous genotypes\">", args, prefix);
    if ( args->tags&SET_AC_Het ) hdr_append(args->out_hdr, "##INFO=<ID=%sAC_Het,Number=A,Type=Integer,Description=\"Allele counts in heterozygous genotypes\">", args, prefix);
    if ( args->tags&SET_AC_Hemi ) hdr_append(args->out_hdr, "##INFO=<ID=%sAC_Hemi,Number=A,Type=Integer,Description=\"Allele counts in hemizygous genotypes\">", args, prefix);
    if ( args->tags&SET_AF ) hdr_append(args->out_hdr, "##INFO=<ID=%sAF,Number=A,Type=Float,Description=\"Allele frequency\">", args, prefix);
    if ( args->tags&SET_MAF ) hdr_append(args->out_hdr, "##INFO=<ID=%sMAF,Number=A,Type=Float,Description=\"Minor Allele frequency\">", args, prefix);
}

void load_samples(args_t *args)
{
    int i, ifile, n_samples_local;

    void *hdr_samples = khash_str2int_init();
    n_samples_local = bcf_hdr_nsamples(args->in_hdr);
    if (args->n_sample_files > 0) {
        args->ids = (char ***) malloc(args->n_sample_files * sizeof(const char**));
        args->n_samples = (int *)malloc(args->n_sample_files * sizeof(int));
        memset(args->n_samples,0,args->n_sample_files * sizeof(int));
        args->samples = (char ***) malloc(args->n_sample_files*sizeof(const char***));
        args->sample_pos = (int **) malloc(args->n_sample_files*sizeof(const int*));
    } else {
        args->ids = (char ***) malloc(1 * sizeof(char**));
        args->n_samples = (int *)malloc(sizeof(int));
        args->n_samples[0] = n_samples_local;
        args->samples = NULL;
        char * empty = "";
        args->ids[0] = (char **)malloc(8 * sizeof(char*)); // array for 7 elements NS,AN,AC,AF,Het,Hom,Hemi
        args->ids[0][POS_NS] = concat_prefix(empty, "NS");
        args->ids[0][POS_AN] = concat_prefix(empty, "AN");
        args->ids[0][POS_AC] = concat_prefix(empty, "AC");
        args->ids[0][POS_AF] = concat_prefix(empty, "AF");
        args->ids[0][POS_MAF] = concat_prefix(empty, "MAF");
        args->ids[0][POS_Het] = concat_prefix(empty, "AC_Het");
        args->ids[0][POS_Hom] = concat_prefix(empty, "AC_Hom");
        args->ids[0][POS_Hemi] = concat_prefix(empty, "AC_Hemi");
        addHeader(args, empty);
    }
    for (i=0; i< n_samples_local; i++) {
        khash_str2int_inc(hdr_samples, bcf_hdr_int2id(args->in_hdr,BCF_DT_SAMPLE,i));
    }
    if (args->n_sample_files == 0) { // load all samples
        args->sample_pos = (int **)  malloc(sizeof(int));
        args->sample_pos[0] = (int *)  malloc((n_samples_local + 1) * sizeof(int));
        args->samples    = (char***) malloc(sizeof(const char**));
        args->samples[0]    = (char**) malloc((n_samples_local + 1) * sizeof(const char*));
        for (i=0; i < args->n_samples[0]; i++)
        {
            args->samples[0][i] = strdup(bcf_hdr_int2id(args->in_hdr,BCF_DT_SAMPLE,i));
            args->sample_pos[0][i] = i;
        }
        args->n_sample_files = 1; // to make it easier to iterate
    } else { // only load specific samples
        for (ifile = 0; ifile < args->n_sample_files; ++ifile) {
            int moff = 0, *off = NULL;
            int ncols = ksplit_core(args->sample_files[ifile],',',&moff,&off);
            if ( ncols>2 ) error("Could not parse the paramteter: %s\n", args->sample_files[ifile]);
            char* pref = "";
            char* file = strdup(args->sample_files[ifile]);
            if (ncols == 2) {
                free(file);
                pref = strdup(&(args->sample_files[ifile][off[0]]));
                file = strdup(&(args->sample_files[ifile][off[1]]));
            }
            args->ids[ifile] = (char **)malloc(8 * sizeof(char*)); // array for 7 elements NS,AN,AC,AF,Het,Hom,Hemi
            args->ids[ifile][POS_NS] = concat_prefix(pref, "NS");
            args->ids[ifile][POS_AN] = concat_prefix(pref, "AN");
            args->ids[ifile][POS_AC] = concat_prefix(pref, "AC");
            args->ids[ifile][POS_AF] = concat_prefix(pref, "AF");
            args->ids[ifile][POS_MAF] = concat_prefix(pref, "MAF");
            args->ids[ifile][POS_Het] = concat_prefix(pref, "AC_Het");
            args->ids[ifile][POS_Hom] = concat_prefix(pref, "AC_Hom");
            args->ids[ifile][POS_Hemi] = concat_prefix(pref, "AC_Hemi");
            addHeader(args, pref);

            void *exclude = (file[0]=='^') ? khash_str2int_init() : NULL;
            int nsmpl;
            char **smpl = NULL;
            args->n_samples[ifile] = 0;
            smpl = hts_readlist(exclude ? &file[1] : file, 1, &nsmpl);
            if ( !smpl )
            {
                error("Could not read the list: \"%s\"\n", exclude ? &file[1] : file);
            }
            int n = 0;
            if ( exclude )
            {
                for (i=0; i<nsmpl; i++) {
                    if (!khash_str2int_has_key(hdr_samples,smpl[i])) {
                        error("Error: exclude called for sample that does not exist in header: \"%s\". Use \"--force-samples\" to ignore this error.\n", smpl[i]);
                    }
                    khash_str2int_inc(exclude, smpl[i]);
                    ++n;
                }

                for (i=0; i<bcf_hdr_nsamples(args->in_hdr); i++)
                {
                    if ( exclude && khash_str2int_has_key(exclude,bcf_hdr_int2id(args->in_hdr,BCF_DT_SAMPLE,i))  ) continue;
                    args->samples[ifile] = (char**) realloc(args->samples[ifile], (args->n_samples[ifile]+1)*sizeof(const char*));
                    args->samples[ifile][args->n_samples[ifile]++] = strdup(bcf_hdr_int2id(args->in_hdr,BCF_DT_SAMPLE,i));
                }
                khash_str2int_destroy(exclude);
            }
            else
            {
                for (i=0; i<nsmpl; i++) {
                    if (!khash_str2int_has_key(hdr_samples,smpl[i])) {
                        error("Error: subset called for sample that does not exist in header: \"%s\". Use \"--force-samples\" to ignore this error.\n", smpl[i]);
                    }
                    if (i == 0) {
                        args->samples[ifile] = (char**) malloc((args->n_samples[ifile]+1)*sizeof(const char**));
                    } else {
                        args->samples[ifile] = (char**) realloc(args->samples[ifile], (args->n_samples[ifile]+1)*sizeof(const char*));
                    }
                    args->samples[ifile][args->n_samples[ifile]++] = strdup(smpl[i]);
                }
            }
            args->sample_pos[ifile] = (int *) malloc( (args->n_samples[ifile] + 1) * sizeof(int));
            for (i=0; i < args->n_samples[ifile]; ++i) {
                args->sample_pos[ifile][i] = khash_str2int_inc(hdr_samples, args->samples[ifile][i]); // build up position index
            }

            for (i=0; i<nsmpl; i++) free(smpl[i]);
            free(smpl);
	    if (args->n_samples[ifile] == 0) {
	      fprintf(stderr, "Warn: subsetting has removed all samples for prefix '%s'\n", pref);
	    }
            if (strlen(pref) > 0) {
                free(pref);
            }
            free(file);
        }
    }
    khash_str2int_destroy(hdr_samples);
}


int parse_tags(args_t *args, const char *str)
{
    int i = 0, n_tags;
    long flag = 0l;
    char **tags = hts_readlist(str, 0, &n_tags);
    for(i=0; i<n_tags; i++)
    {
        if ( !strcasecmp(tags[i],"AN") ) flag |= SET_AN;
        else if ( !strcasecmp(tags[i],"AC") ) flag |= SET_AC;
        else if ( !strcasecmp(tags[i],"NS") ) flag |= SET_NS;
        else if ( !strcasecmp(tags[i],"AC_Hom") ) flag |= SET_AC_Hom;
        else if ( !strcasecmp(tags[i],"AC_Het") ) flag |= SET_AC_Het;
        else if ( !strcasecmp(tags[i],"AC_Hemi") ) flag |= SET_AC_Hemi;
        else if ( !strcasecmp(tags[i],"AF") ) flag |= SET_AF;
        else if ( !strcasecmp(tags[i],"MAF") ) flag |= SET_MAF;
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

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.in_hdr  = in;
    args.out_hdr = out;

    static struct option loptions[] =
    {
        {"drop-missing",0,0,'d'},
        {"tags",1,0,'t'},
        {"samples-file",1,0,'S'},
        {0,0,0,0}
    };
    args.n_sample_files = 0;
    args.sample_files = NULL;
    int c;
    while ((c = getopt_long(argc, argv, "?ht:dS:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'd': args.drop_missing = 1; break;
            case 't': args.tags |= parse_tags(&args,optarg); break;
            case 'S': 
                args.n_sample_files++;
                args.sample_files = (char**) realloc(args.sample_files, args.n_sample_files*sizeof(char*));
                args.sample_files[args.n_sample_files-1] = strdup(optarg);
                break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }

    if ( optind != argc ) error(usage());
    args.gt_id = bcf_hdr_id2int(args.in_hdr,BCF_DT_ID,"GT");
    if ( args.gt_id<0 ) error("Error: GT field is not present\n");
    if ( !args.tags ) args.tags |= SET_AN|SET_AC|SET_NS|SET_AC_Hom|SET_AC_Het|SET_AC_Hemi|SET_AF|SET_MAF;
    load_samples(&args);
    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    int nsf, i;

    bcf_unpack(rec, BCF_UN_FMT);
    bcf_fmt_t *fmt_gt = NULL;
    for (i=0; i<rec->n_fmt; i++)
        if ( rec->d.fmt[i].id==args.gt_id ) { fmt_gt = &rec->d.fmt[i]; break; }
    if ( !fmt_gt ) return rec;    // no GT tag

    hts_expand(int32_t,rec->n_allele, args.marr, args.arr);
    hts_expand(float,rec->n_allele, args.mfarr, args.farr);
    hts_expand(counts_t,rec->n_allele,args.mcounts, args.counts);

    for (nsf=0; nsf < args.n_sample_files; ++nsf) {
        int pos, ns = 0, is_hom, is_hemi, is_half;
        memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
        memset(args.counts,0,sizeof(*args.counts)*rec->n_allele);

       #define BRANCH_INT(type_t,vector_end) { \
         for (pos=0; pos<args.n_samples[nsf]; pos++) \
         { \
            i = args.sample_pos[nsf][pos]; \
            type_t *p = (type_t*) (fmt_gt->p + i*fmt_gt->size); \
            int ial, als = 0, nals = 0; \
            for (ial=0; ial<fmt_gt->n; ial++) \
            { \
                if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                if ( bcf_gt_is_missing(p[ial]) ) continue; /* missing allele */ \
                int idx = bcf_gt_allele(p[ial]); \
                nals++; \
                \
                if ( idx >= rec->n_allele ) \
                    error("Incorrect allele (\"%d\") in %s at %s:%d\n",idx,args.in_hdr->samples[i],bcf_seqname(args.in_hdr,rec),rec->pos+1); \
                als |= (1<<idx);  /* this breaks with too many alleles */ \
            } \
            if ( nals==0 ) continue; /* missing genotype */ \
            ns++; \
            is_hom = als && !(als & (als-1)); /* only one bit is set */ \
            if ( nals!=ial ) \
            { \
                if ( args.drop_missing ) is_hemi = 0, is_half = 1; \
                else is_hemi = 1, is_half = 0; \
            } \
            else if ( nals==1 ) is_hemi = 1, is_half = 0; \
            else is_hemi = 0, is_half = 0; \
            for (ial=0; als; ial++) \
            { \
                if ( als&1 ) \
                { \
                    if ( is_half ) \
                        args.counts[ial].nac++; \
                    else if ( !is_hom ) \
                        args.counts[ial].nhet++; \
                    else if ( !is_hemi ) \
                        args.counts[ial].nhom += 2; \
                    else \
                        args.counts[ial].nhemi++; \
                } \
                als >>= 1; \
            } \
          } \
        }
        switch (fmt_gt->type) {
            case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
            default: error("The GT type is not recognised: %d at %s:%d\n",fmt_gt->type, bcf_seqname(args.in_hdr,rec),rec->pos+1); break;
        }

        #undef BRANCH_INT

        if ( args.tags&SET_NS )
        {
            char * id = args.ids[nsf][POS_NS];
            if ( bcf_update_info_int32(args.out_hdr,rec,id,&ns,1)!=0 )
                error("Error occurred while updating NS at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
        if ( args.tags&SET_AN )
        {
            char * id = args.ids[nsf][POS_AN];
            args.arr[0] = 0;
            for (i=0; i<rec->n_allele; i++)
                args.arr[0] += args.counts[i].nhet + args.counts[i].nhom + args.counts[i].nhemi + args.counts[i].nac;
            if ( bcf_update_info_int32(args.out_hdr,rec,id,args.arr,1)!=0 )
                error("Error occurred while updating AN at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
        if ( args.tags&SET_AF || args.tags&SET_MAF )
        {
            int n = rec->n_allele-1;
            if ( n>0 )
            {
                args.arr[0] = 0;
                for (i=0; i<rec->n_allele; i++)
                    args.arr[0] += args.counts[i].nhet + args.counts[i].nhom + args.counts[i].nhemi + args.counts[i].nac;
                for (i=1; i<rec->n_allele; i++)
                    args.farr[i] = (args.counts[i].nhet + args.counts[i].nhom + args.counts[i].nhemi + args.counts[i].nac)*1.0/args.arr[0];
            }
            if ( args.arr[0] )
            {
                if (args.tags&SET_AF) {
                char * id = args.ids[nsf][POS_AF];
                if ( bcf_update_info_float(args.out_hdr,rec,id,args.farr+1,n)!=0 )
                    error("Error occurred while updating AF at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
                }
                if (args.tags&SET_MAF) {
                    for (i=1; i<rec->n_allele; i++) {
                        if (args.farr[i] > 0.5) {
                            args.farr[i] = 1 - args.farr[i];
                        }
                    }
                    char * id = args.ids[nsf][POS_MAF];
                    if ( bcf_update_info_float(args.out_hdr,rec,id,args.farr+1,n)!=0 )
                        error("Error occurred while updating MAF at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
                }
            }
        }
        if ( args.tags&SET_AC )
        {
            int n = rec->n_allele-1;
            if ( n>0 )
            {
                memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
                for (i=1; i<rec->n_allele; i++)
                    args.arr[i] = args.counts[i].nhet + args.counts[i].nhom + args.counts[i].nhemi + args.counts[i].nac;
            }
            char * id = args.ids[nsf][POS_AC];
            if ( bcf_update_info_int32(args.out_hdr,rec,id,args.arr+1,n)!=0 )
                error("Error occurred while updating AC at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
        if ( args.tags&SET_AC_Het )
        {
            int n = rec->n_allele-1;
            if ( n>0 )
            {
                memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
                for (i=1; i<rec->n_allele; i++)
                    args.arr[i] += args.counts[i].nhet;
            }
            char * id = args.ids[nsf][POS_Het];
            if ( bcf_update_info_int32(args.out_hdr,rec,id,args.arr+1,n)!=0 )
                error("Error occurred while updating AC_Het at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
        if ( args.tags&SET_AC_Hom )
        {
            int n = rec->n_allele-1;
            if ( n>0 )
            {
                memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
                for (i=1; i<rec->n_allele; i++)
                    args.arr[i] += args.counts[i].nhom;
            }
            char * id = args.ids[nsf][POS_Hom];
            if ( bcf_update_info_int32(args.out_hdr,rec,id,args.arr+1,n)!=0 )
                error("Error occurred while updating AC_Hom at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
        if ( args.tags&SET_AC_Hemi )
        {
            int n = rec->n_allele-1;
            if ( n>0 )
            {
                memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
                for (i=1; i<rec->n_allele; i++)
                    args.arr[i] += args.counts[i].nhemi;
            }
            char * id = args.ids[nsf][POS_Hemi];
            if ( bcf_update_info_int32(args.out_hdr,rec,id,args.arr+1,n)!=0 )
                error("Error occurred while updating AC_Hemi at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
    }
    return rec;
}

void destroy(void)
{
    int i,j;
    for (i=0; i < args.n_sample_files; ++i) {
        for (j = 0; j < args.n_samples[i]; ++j) {
            free(args.samples[i][j]);
        }
        free(args.samples[i]);
        free(args.sample_pos[i]);
        free(args.ids[i][POS_NS]);
        free(args.ids[i][POS_AN]);
        free(args.ids[i][POS_AC]);
        free(args.ids[i][POS_AF]);
        free(args.ids[i][POS_MAF]);
        free(args.ids[i][POS_Het]);
        free(args.ids[i][POS_Hom]);
        free(args.ids[i][POS_Hemi]);
        free(args.ids[i]);
    }

    free(args.n_samples);
    free(args.samples);
    free(args.sample_pos);
    free(args.ids);

    free(args.counts);
    free(args.arr);
    free(args.farr);
}



