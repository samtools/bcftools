/* The MIT License

   Copyright (c) 2016 Genome Research Ltd.

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
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"


#define SET_AN  (1<<0)
#define SET_AC  (1<<1)
#define SET_AF  (1<<2)
#define SET_HWE (1<<3)

typedef struct
{
    char *name;
    int nsmpl, *smpl;
}
pop_t;

typedef struct
{
    bcf_hdr_t *in_hdr, *out_hdr;
    int npop, tags;
    pop_t *pop;
    int32_t *gt_arr, *ac, *tot_ac;
    float *af;
    int mgt_arr, mac, maf, mtot_ac;
    double *hwe_probs;
    int mhwe_probs;
    kstring_t str;
}
args_t;

args_t args;

const char *about(void)
{
    return "Fill population INFO tags AN, AC, AF, HWE.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Fill AN_*,AC_*,AF_*,HWE_* tags, for each population and across all samples.\n"
        "Usage: bcftools +aggregate-tags [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -S, --samples-file <file>    List of samples (first column) and comma-separated list of populations (second column)\n"
        "   -t, --tags <list>            Comma-separated list of tags to fill. By default, all tags are filled\n"
        "\n"
        "Example:\n"
        "   bcftools +aggregate-tags file.bcf -- -S samples.txt\n"
        "\n";
}

int cmpint(const void *a, const void *b)
{
    if ( *((const int*)a) < *((const int*)b) ) return -1;
    if ( *((const int*)a) == *((const int*)b) ) return 0;
    return 1;
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
        // NA12400 CEU
        // NA18507 YRI
        int ncols = ksplit_core(str.s,'\t',&moff,&off);
        if ( ncols!=2 ) error("Could not parse the ped file: %s\n", str.s);

        char *smpl = strdup(&str.s[off[0]]);
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
        khash_str2int_inc(smpli,smpl);

        char *pop_names = &str.s[off[1]];
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

    int i;
    if ( nsmpl != bcf_hdr_nsamples(args->in_hdr) )
    {
        fprintf(stderr,"Warning: %d samples in the list, %d samples in the VCF. Unlisted samples marked as UNKN.\n", nsmpl,bcf_hdr_nsamples(args->in_hdr));

        // unlisted samples
        for (i=0; i<bcf_hdr_nsamples(args->in_hdr); i++)
        {   
            if ( khash_str2int_has_key(smpli,bcf_hdr_int2id(args->in_hdr,BCF_DT_SAMPLE,i)) ) continue;  // sample already assigned to a population

            char *pop_name = "UNKN";
            if ( !khash_str2int_has_key(pop2i,pop_name) )
            {
                pop_name = strdup(pop_name);
                khash_str2int_set(pop2i,pop_name,args->npop);
                args->npop++;
                args->pop = (pop_t*) realloc(args->pop,args->npop*sizeof(*args->pop));
                memset(args->pop+args->npop-1,0,sizeof(*args->pop));
                args->pop[args->npop-1].name = pop_name;
            }
            int ipop = 0;
            khash_str2int_get(pop2i,pop_name,&ipop);
            pop_t *pop = &args->pop[ipop];
            pop->nsmpl++;
            pop->smpl = (int*) realloc(pop->smpl,pop->nsmpl*sizeof(*pop->smpl));
            pop->smpl[pop->nsmpl-1] = i;
            nsmpl++;
        }
        if ( nsmpl != bcf_hdr_nsamples(args->in_hdr) ) error("fixme, this should not happen: %d %d\n",nsmpl,bcf_hdr_nsamples(args->in_hdr));
    }

    for (i=0; i<args->npop; i++)
        qsort(args->pop[i].smpl, args->pop[i].nsmpl, sizeof(*args->pop[i].smpl), cmpint);

    khash_str2int_destroy(pop2i);
    khash_str2int_destroy_free(smpli);
    free(str.s);
    free(off);
    hts_close(fp);
}

int parse_tags(args_t *args, const char *str)
{
    int i, flag = 0, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags);
    for(i=0; i<n_tags; i++)
    {
        if ( !strcasecmp(tags[i],"AN") ) flag |= SET_AN;
        else if ( !strcasecmp(tags[i],"AC") ) flag |= SET_AC;
        else if ( !strcasecmp(tags[i],"AF") ) flag |= SET_AF;
        else if ( !strcasecmp(tags[i],"HWE") ) flag |= SET_HWE;
        else
        {
            fprintf(stderr,"Unknown tag \"%s\" in \"%s\"\n", tags[i], str);
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
    char *samples_fname = NULL;
    static struct option loptions[] =
    {
        {"samples-file",1,0,'S'},
        {"tags",1,0,'t'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?hS:t:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 't': args.tags |= parse_tags(&args,optarg); break;
            case 'S': samples_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !samples_fname ) error("Expected the -S option\n");
    parse_samples(&args, samples_fname);

    if ( !args.tags ) args.tags = SET_AN|SET_AC|SET_AF|SET_HWE;

    int i;
    if ( args.tags&SET_AN )
    {
        bcf_hdr_append(args.out_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
        for (i=0; i<args.npop; i++)
            bcf_hdr_printf(args.out_hdr, "##INFO=<ID=AN_%s,Number=1,Type=Integer,Description=\"Total number of alleles in %s\">",args.pop[i].name,args.pop[i].name);
    }

    if ( args.tags&SET_AC )
    {
        bcf_hdr_append(args.out_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">");
        for (i=0; i<args.npop; i++)
            bcf_hdr_printf(args.out_hdr, "##INFO=<ID=AC_%s,Number=A,Type=Integer,Description=\"Allele count in %s\">",args.pop[i].name,args.pop[i].name);
    }

    if ( args.tags&SET_AF )
    {
        bcf_hdr_append(args.out_hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">");
        for (i=0; i<args.npop; i++)
            bcf_hdr_printf(args.out_hdr, "##INFO=<ID=AF_%s,Number=A,Type=Float,Description=\"Allele frequency in %s\">",args.pop[i].name,args.pop[i].name);
    }

    if ( args.tags&SET_HWE )
    {
        bcf_hdr_append(args.out_hdr, "##INFO=<ID=HWE,Number=1,Type=Float,Description=\"HWE test (PMID:15789306)\">");
        for (i=0; i<args.npop; i++)
            bcf_hdr_printf(args.out_hdr, "##INFO=<ID=HWE_%s,Number=1,Type=Float,Description=\"HWE test in %s\">",args.pop[i].name,args.pop[i].name);
    }

    return 0;
}

/* 
    Wigginton 2005, PMID: 15789306 

    nref .. number of reference alleles
    nalt .. number of alt alleles
    nhet .. number of het genotypes, assuming number of genotypes = (nref+nalt)*2

*/
float calc_hwe(args_t *args, int nref, int nalt, int nhet)
{
    int ngt   = (nref+nalt) / 2;
    int nrare = nref < nalt ? nref : nalt;

    // sanity check: there is odd/even number of rare alleles iff there is odd/even number of hets
    if ( (nrare & 1) ^ (nhet & 1) ) error("nrare/nhet should be both odd or even: nrare=%d nref=%d nalt=%d nhet=%d\n",nrare,nref,nalt,nhet);
    if ( nrare < nhet ) error("More rare alleles than hets? nrare=%d nref=%d nalt=%d nhet=%d\n",nrare,nref,nalt,nhet);
    if ( (nref+nalt) & 1 ) error("Expected diploid genotypes: nref=%d nalt=%d\n",nref,nalt);

    // initialize het probs
    hts_expand(double,nrare+1,args->mhwe_probs,args->hwe_probs);
    memset(args->hwe_probs, 0, sizeof(*args->hwe_probs)*(nrare+1));
    double *probs = args->hwe_probs;

    // start at midpoint
    int mid = nrare * (nref + nalt - nrare) / (nref + nalt);

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

    double p_rank = 0.0;
    for (het=0; het <= nrare; het++)
    {
        if ( probs[het] > probs[nhet]) continue;
        p_rank += probs[het];
    }

    return p_rank > 1 ? 1.0 : p_rank;
}

bcf1_t *process(bcf1_t *rec)
{
    int ngt = bcf_get_genotypes(args.in_hdr, rec, &args.gt_arr, &args.mgt_arr);
    if ( ngt<0 ) return rec;

    int nsmpl = bcf_hdr_nsamples(args.in_hdr);
    if ( ngt % nsmpl != 0 ) return rec;  // never happens
    ngt /= nsmpl;
    if ( ngt!=2 ) return rec;   // not diploid

    hts_expand(int32_t,rec->n_allele,args.mac,args.ac);
    hts_expand(int32_t,rec->n_allele,args.mtot_ac,args.tot_ac);
    hts_expand(float,rec->n_allele,args.maf,args.af);

    int i;
    for (i=0; i<rec->n_allele; i++) args.tot_ac[i] = 0;

    int ipop, tot_an = 0, tot_nhet = 0;
    for (ipop=0; ipop<args.npop; ipop++)
    {
        pop_t *pop = &args.pop[ipop];

        int nhet = 0;
        for (i=0; i<rec->n_allele; i++) args.ac[i] = 0;
        for (i=0; i<pop->nsmpl; i++)
        {
            int32_t *ptr = args.gt_arr + ngt*pop->smpl[i];
            if ( ptr[1]==bcf_int32_vector_end ) continue;
            if ( bcf_gt_is_missing(ptr[0]) ) continue;
            if ( bcf_gt_is_missing(ptr[1]) ) continue;

            int a = bcf_gt_allele(ptr[0]);
            if ( a >= rec->n_allele ) error("Wrong allele index at %s:%d\n",bcf_seqname(args.in_hdr,rec),rec->pos+1);

            int b = bcf_gt_allele(ptr[1]);
            if ( b >= rec->n_allele ) error("Wrong allele index at %s:%d\n",bcf_seqname(args.in_hdr,rec),rec->pos+1);

            args.ac[a]++;
            args.ac[b]++;
            if ( (!a && b) || (a && !b) ) nhet++;     // approximation: collapse all alt alleles into one
        }

        int32_t ac = 0, an = 0;
        for (i=0; i<rec->n_allele; i++)
        {
            args.tot_ac[i] += args.ac[i];
            an += args.ac[i];
        }
        for (i=1; i<rec->n_allele; i++) ac += args.ac[i];   // collapsed alt ac
        tot_an   += an;
        tot_nhet += nhet;

        if ( args.tags & SET_AN )
        {
            args.str.l = 0;
            ksprintf(&args.str, "AN_%s", pop->name);
            if ( bcf_update_info_int32(args.out_hdr,rec,args.str.s,&an,1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args.str.s,bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
        if ( args.tags & SET_AC )
        {
            args.str.l = 0;
            ksprintf(&args.str, "AC_%s", pop->name);
            if ( bcf_update_info_int32(args.out_hdr,rec,args.str.s,args.ac+1,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args.str.s,bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
        if ( args.tags & SET_AF )
        {
            for (i=1; i<rec->n_allele; i++) args.af[i] = an ? (float)args.ac[i]/an : 0;
            args.str.l = 0;
            ksprintf(&args.str, "AF_%s", pop->name);
            if ( bcf_update_info_float(args.out_hdr,rec,args.str.s,args.af+1,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args.str.s,bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
        if ( args.tags & SET_HWE )
        {
            float hwe = (args.ac[0]>0 && ac>0) ? calc_hwe(&args, args.ac[0], ac, nhet) : -1;
            args.str.l = 0;
            ksprintf(&args.str, "HWE_%s", pop->name);
            if ( bcf_update_info_float(args.out_hdr,rec,args.str.s,&hwe,1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args.str.s,bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
    }

    if ( args.tags & SET_AN )
    {
        if ( bcf_update_info_int32(args.out_hdr,rec,"AN",&tot_an,1)!=0 )
            error("Error occurred while updating AN at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    if ( args.tags & SET_AC )
    {
        if ( bcf_update_info_int32(args.out_hdr,rec,"AC",args.tot_ac+1,rec->n_allele-1)!=0 )
            error("Error occurred while updating AC at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    if ( args.tags & SET_AF )
    {
        for (i=1; i<rec->n_allele; i++) args.af[i] = tot_an ? (float)args.tot_ac[i]/tot_an : 0;
        if ( bcf_update_info_float(args.out_hdr,rec,"AF",args.af+1,rec->n_allele-1)!=0 )
            error("Error occurred while updating AF at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    if ( args.tags & SET_HWE )
    {
        int32_t tot_ac = 0;
        for (i=1; i<rec->n_allele; i++) tot_ac += args.tot_ac[i];
        float hwe = (args.tot_ac[0]>0 && tot_ac>0) ? calc_hwe(&args, args.tot_ac[0], tot_ac, tot_nhet) : -1;
        if ( bcf_update_info_float(args.out_hdr,rec,"HWE",&hwe,1)!=0 )
            error("Error occurred while updating HWE at %s:%d\n", args.str.s,bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    return rec;
}

void destroy(void)
{
    int i;
    for (i=0; i<args.npop; i++) 
    {
        free(args.pop[i].smpl);
        free(args.pop[i].name);
    }
    free(args.pop);
    free(args.gt_arr);
    free(args.hwe_probs);
    free(args.ac);
    free(args.af);
    free(args.tot_ac);
    free(args.str.s);
}
