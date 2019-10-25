/* The MIT License

   Copyright (c) 2015-2019 Genome Research Ltd.

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

typedef struct _args_t args_t;
typedef struct _ftf_t ftf_t;
typedef int (*fill_tag_f)(args_t *, bcf1_t *, ftf_t *);
struct _ftf_t
{
    char *src_tag, *dst_tag;
    fill_tag_f func;
    int *pop_vals;      // for now assuming only 1 integer value per annotation
};

typedef struct
{
    int nhom, nhet, nhemi, nac;
}
counts_t;

typedef struct
{
    int ns;
    int ncounts, mcounts;
    counts_t *counts;
    char *name, *suffix;
    int nsmpl, *smpl;
}
pop_t;

struct _args_t
{
    bcf_hdr_t *in_hdr, *out_hdr;
    int npop, tags, drop_missing, gt_id;
    pop_t *pop, **smpl2pop;
    float *farr;
    int32_t *iarr, niarr, miarr, nfarr, mfarr;
    double *hwe_probs;
    int mhwe_probs;
    kstring_t str;
    kbitset_t *bset;
    ftf_t *ftf;
    int nftf;
};

static args_t *args;

const char *about(void)
{
    return "Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS and more.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS\n"
        "   or custom INFO/TAG=func(FMT/TAG), use -l for detailed description\n"
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
        "   # Fill all available tags\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t all\n"
        "\n"
        "   # Calculate HWE for sample groups (possibly multiple) read from a file\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -S sample-group.txt -t HWE\n"
        "\n"
        "   # Calculate total read depth (INFO/DP) from per-sample depths (FORMAT/DP)\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t 'DP=sum(DP)'\n"
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

void ftf_destroy(args_t *args)
{
    int i;
    for (i=0; i<args->nftf; i++)
    {
        ftf_t *ftf = &args->ftf[i];
        free(ftf->src_tag);
        free(ftf->dst_tag);
        free(ftf->pop_vals);
    }
    free(args->ftf);
}
int ftf_sum(args_t *args, bcf1_t *rec, ftf_t *ftf)
{
    int nsmpl = bcf_hdr_nsamples(args->in_hdr);
    int nval = bcf_get_format_int32(args->in_hdr, rec, ftf->src_tag, &args->iarr, &args->miarr);
    if ( nval<=0 ) return 0;
    nval /= nsmpl;

    int i;
    for (i=0; i<args->npop; i++)
        ftf->pop_vals[i] = -1;

    for (i=0; i<nsmpl; i++)
    {
        if ( args->iarr[i*nval]==bcf_int32_missing || args->iarr[i*nval]==bcf_int32_vector_end ) continue;

        pop_t **pop = &args->smpl2pop[i*(args->npop+1)];
        while ( *pop )
        {
            int ipop = (int)(*pop - args->pop);
            if ( ftf->pop_vals[ipop]<0 ) ftf->pop_vals[ipop] = 0;
            ftf->pop_vals[ipop] += args->iarr[i*nval];
            pop++;
        }
    }

    for (i=0; i<args->npop; i++)
    {
        if ( ftf->pop_vals[i]<0 ) continue;
        args->str.l = 0;
        ksprintf(&args->str, "%s%s", ftf->dst_tag,args->pop[i].suffix);
        if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,ftf->pop_vals+i,1)!=0 )
            error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
    }

    return 0;
}

void hdr_append(args_t *args, char *fmt)
{
    int i;
    for (i=0; i<args->npop; i++)
        bcf_hdr_printf(args->out_hdr, fmt, args->pop[i].suffix,*args->pop[i].name ? " in " : "",args->pop[i].name);
}

int parse_func(args_t *args, char *tag, char *expr)
{
    args->nftf++;
    args->ftf = (ftf_t *)realloc(args->ftf,sizeof(*args->ftf)*args->nftf);
    ftf_t *ftf = &args->ftf[ args->nftf - 1 ];

    ftf->pop_vals = (int*)calloc(args->npop,sizeof(*ftf->pop_vals));
    ftf->dst_tag = (char*)calloc(expr-tag,1);
    memcpy(ftf->dst_tag, tag, expr-tag-1);

    if ( !strncasecmp(expr,"sum(",4) ) { ftf->func  = ftf_sum; expr += 4; }
    else error("Error: the expression not recognised: %s\n",tag);

    char *tmp = expr; 
    while ( *tmp && *tmp!=')' ) tmp++;
    if ( !*tmp ) error("Error: could not parse: %s\n",tag);

    ftf->src_tag = (char*)calloc(tmp-expr+2,1);
    memcpy(ftf->src_tag, expr, tmp-expr);

    int id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,ftf->src_tag);
    if ( !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,id) ) error("Error: the field FORMAT/%s is not present\n",ftf->src_tag);

    int i = 0;
    for (i=0; i<args->npop; i++)
    {
        args->str.l = 0;
        ksprintf(&args->str, "%s%s", ftf->dst_tag,args->pop[i].suffix);
        id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,args->str.s);
        if ( bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,id) )
        {
            if ( bcf_hdr_id2length(args->in_hdr,BCF_HL_FMT,id)!=BCF_VL_FIXED )
                error("Error: the field INFO/%s already exists with a definition different from Number=1\n",args->str.s);
            if ( bcf_hdr_id2number(args->in_hdr,BCF_HL_FMT,id)!=1 )
                error("Error: the field INFO/%s already exists with a definition different from Number=1\n",args->str.s);
            if ( bcf_hdr_id2type(args->in_hdr,BCF_HT_INT,id)!=BCF_HT_INT )
                error("Error: the field INFO/%s already exists with a definition different from Type=Integer\n",args->str.s);
        }
        else
            bcf_hdr_printf(args->out_hdr, "##INFO=<ID=%s,Number=1,Type=Integer,Description=\"%s%s%s\">",args->str.s,tag,*args->pop[i].name ? " in " : "",args->pop[i].name);
    }
    return SET_FUNC;
}
int parse_tags(args_t *args, const char *str)
{
    if ( !args->in_hdr ) error("%s", usage());

    int i,j, flag = 0, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags), *ptr;
    for(i=0; i<n_tags; i++)
    {
        if ( !strcasecmp(tags[i],"all") )
        {
            for (j=0; j<=10; j++) flag |= 1<<j;
        }
        else if ( !strcasecmp(tags[i],"AN") ) flag |= SET_AN;
        else if ( !strcasecmp(tags[i],"AC") ) flag |= SET_AC;
        else if ( !strcasecmp(tags[i],"NS") ) flag |= SET_NS;
        else if ( !strcasecmp(tags[i],"AC_Hom") ) flag |= SET_AC_Hom;
        else if ( !strcasecmp(tags[i],"AC_Het") ) flag |= SET_AC_Het;
        else if ( !strcasecmp(tags[i],"AC_Hemi") ) flag |= SET_AC_Hemi;
        else if ( !strcasecmp(tags[i],"AF") ) flag |= SET_AF;
        else if ( !strcasecmp(tags[i],"MAF") ) flag |= SET_MAF;
        else if ( !strcasecmp(tags[i],"HWE") ) flag |= SET_HWE;
        else if ( !strcasecmp(tags[i],"ExcHet") ) flag |= SET_ExcHet;
        else if ( (ptr=strchr(tags[i],'=')) ) flag |= parse_func(args,tags[i],ptr+1);
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
        "INFO/AN       Number:1  Type:Integer  ..  Total number of alleles in called genotypes\n"
        "INFO/AC       Number:A  Type:Integer  ..  Allele count in genotypes\n"
        "INFO/NS       Number:1  Type:Integer  ..  Number of samples with data\n"
        "INFO/AC_Hom   Number:A  Type:Integer  ..  Allele counts in homozygous genotypes\n"
        "INFO/AC_Het   Number:A  Type:Integer  ..  Allele counts in heterozygous genotypes\n"
        "INFO/AC_Hemi  Number:A  Type:Integer  ..  Allele counts in hemizygous genotypes\n"
        "INFO/AF       Number:A  Type:Float    ..  Allele frequency\n"
        "INFO/MAF      Number:A  Type:Float    ..  Minor Allele frequency\n"
        "INFO/HWE      Number:A  Type:Float    ..  HWE test (PMID:15789306); 1=good, 0=bad\n"
        "INFO/ExcHet   Number:A  Type:Float    ..  Test excess heterozygosity; 1=good, 0=bad\n"
        "TAG=func(TAG) Number:1  Type:Integer  ..  Experimental support for user-defined\n"
        "    expressions such as \"DP=sum(DP)\". This is currently very basic, to be extended.\n"
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

    args->gt_id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"GT");
    if ( args->gt_id<0 ) error("Error: GT field is not present\n");

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
    if ( args->tags & SET_MAF ) hdr_append(args, "##INFO=<ID=MAF%s,Number=A,Type=Float,Description=\"Minor Allele frequency%s%s\">");
    if ( args->tags & SET_HWE ) hdr_append(args, "##INFO=<ID=HWE%s,Number=A,Type=Float,Description=\"HWE test%s%s (PMID:15789306); 1=good, 0=bad\">");
    if ( args->tags & SET_ExcHet ) hdr_append(args, "##INFO=<ID=ExcHet%s,Number=A,Type=Float,Description=\"Test excess heterozygosity%s%s; 1=good, 0=bad\">");

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

bcf1_t *process(bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_FMT);

    int i,j, nsmpl = bcf_hdr_nsamples(args->in_hdr);;

    for (i=0; i<args->nftf; i++)
        args->ftf[i].func(args, rec, &args->ftf[i]);

    bcf_fmt_t *fmt_gt = NULL;
    for (i=0; i<rec->n_fmt; i++)
        if ( rec->d.fmt[i].id==args->gt_id ) { fmt_gt = &rec->d.fmt[i]; break; }
    if ( !fmt_gt ) return rec;    // no GT tag

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
                memset(args->farr, 0, sizeof(*args->farr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++) 
                    args->farr[j-1] += pop->counts[j].nhet + pop->counts[j].nhom + pop->counts[j].nhemi + pop->counts[j].nac;
                an = pop->counts[0].nhet + pop->counts[0].nhom + pop->counts[0].nhemi + pop->counts[0].nac;
                for (j=1; j<rec->n_allele; j++) an += args->farr[j-1];
                if ( an )
                    for (j=1; j<rec->n_allele; j++) args->farr[j-1] /= an;
                else
                    for (j=1; j<rec->n_allele; j++) bcf_float_set_missing(args->farr[j-1]);
            }
            if ( args->tags & SET_AF )
            {
                args->str.l = 0;
                ksprintf(&args->str, "AF%s", args->pop[i].suffix);
                if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,args->farr,rec->n_allele-1)!=0 )
                    error("Error occurred while updating %s at %s:%"PRId64"\n", args->str.s,bcf_seqname(args->in_hdr,rec),(int64_t) rec->pos+1);
            }
            if ( args->tags & SET_MAF )
            {
                if ( an )
                {
                    for (j=1; j<rec->n_allele; j++)
                        if ( args->farr[j-1] > 0.5 ) args->farr[j-1] = 1 - args->farr[j-1];     // todo: this is incorrect for multiallelic sites
                }
                args->str.l = 0;
                ksprintf(&args->str, "MAF%s", args->pop[i].suffix);
                if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,args->farr,rec->n_allele-1)!=0 )
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
                memset(args->iarr, 0, sizeof(*args->iarr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++) 
                    args->iarr[j-1] += pop->counts[j].nhet + pop->counts[j].nhom + pop->counts[j].nhemi + pop->counts[j].nac;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr,rec->n_allele-1)!=0 )
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
    }
    kbs_destroy(args->bset);
    free(args->str.s);
    free(args->pop);
    free(args->smpl2pop);
    free(args->iarr);
    free(args->farr);
    free(args->hwe_probs);
    ftf_destroy(args);
    free(args);
}



