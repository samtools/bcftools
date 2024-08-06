/* The MIT License

   Copyright (c) 2024 Genome Research Ltd.

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
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>     // for isatty
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <assert.h>
#include <errno.h>
#include "mpileup2/mpileup.h"
#include "bcftools.h"
#include "regidx.h"


// VAF profile across many samples at a single site. Used as a regidx payload
typedef struct
{
    char ref, alt;      // todo: indels
    uint32_t no_good:1, // is the site suitable for training?
             nval:31,   // number of values
            *dist;      // histogram of VAF frequencies
}
site_t;

typedef struct
{
    char *fname;
    int nbins;          // number of VAF bins
    int nval;           // number of values contributing to the computed profile (ie number of good sites)
    double *mean, *var; // mean and variance of each VAF bin across good sites indexed by regidx
    regidx_t *idx;      // sites in the batch; all batches should have the same set of sites
    site_t good_calls;  // cumulative VAF histogram from good calls, ie binom(VAF)<=binom_th
}
batch_t;

typedef struct
{
    int argc, output_type, record_cmd_line, clevel, use_bam_idx;
    char **argv, *output_fname, *sites_fname, *aln_fname, *fasta_fname, *batch_fname, *batch;
    char **bams;
    int nbams;
    double binom_th;    // include in the profile only calls with binom(VAF)<=binom_th
    BGZF *out_fh;
    regidx_t *sites_idx;
    regitr_t *sites_itr;
    mpileup_t *mplp;
    batch_t profile;
}
args_t;

const char *about(void)
{
    return "Localised assessment of sequencing artefacts, estimate site noisiness\n";
}
static const char *usage_text(void)
{
    return
        "\n"
        "About: Assess site noisiness from a large number of unaffected parental samples\n"
        "Usage: bcftools +noise-profile [OPTIONS]\n"
        "\n"
        "Required options:\n"
        "   -a, --alns FILE               List of BAM/CRAM files\n"
        "   -f, --fasta-ref FILE          Reference file in fasta format\n"
        "   -s, --sites FILE              A tab-delimited file name of sites to assess (chr,pos,ref,alt)\n"
        "\n"
        "Model options:\n"
        "   -B, --binom-th FLOAT          Maximum p-value of calls included in the expected profile [1e-4]\n"
        "   -n, --nbins INT               Number of VAF bins [20]\n"
        "\n"
        "Other options:\n"
        "   -b, --batch I/N               Run I-th batch out of N, 1-based\n"
        "   -i, --use-index               Use index to jump, rather than stream, the alignments\n"
        "   -m, --merge-batches FILE      Merge files produced with -b, --batch\n"
        "   -o, --output FILE             Output file name [stdout]\n"
        "   -O, --output-type t|z[0-9]    t/z: un/compressed text file, 0-9: compression level [t]\n"
        "\n"
        "Example:\n"
        "   # Typical run\n"
        "   bcftools +noise-profile -f ref.fa -a bams.txt -s sites.txt -o scores.txt\n"
        "\n"
        "   # Run in batches. Let's say one batch can have at most 3 bams and there are 5 bams in total\n"
        "   #   1) find out the number of required batches with `--batch k=3` (the program outputs: 2)\n"
        "   #   2) run all batches individually with `--batch 1/2` and `--batch 2/2`\n"
        "   #   3) create a list of outputs and merge with `--merge-batches list.txt`\n"
        "   bcftools +noise-profile -a bams.txt --batch k=3     # prints 2\n"
        "   bcftools +noise-profile -f ref.fa -a bams.txt -s sites.txt -o scores1.txt --batch 1/2\n"
        "   bcftools +noise-profile -f ref.fa -a bams.txt -s sites.txt -o scores2.txt --batch 2/2\n"
        "   bcftools +noise-profile --merge-batches list.txt\n"
        "\n";
}

static int parse_sites(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    // parses space-delimited lines like this:
    //      chr	3	C	A
    args_t *args = (args_t*) usr;

    // CHR part
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && !isspace(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se ) error("Could not parse the CHR part of the line: %s\n",line);

    // POS part
    ss = se+1;
    *beg = strtod(ss, &se);
    if ( ss==se ) error("Could not parse the POS part of the line: %s\n",line);
    (*beg)--;

    // REF part and REF length
    ss = ++se;
    while ( *se && !isspace(*se) ) se++;
    *end = *beg + se-ss-1;

    site_t *site = (site_t*)payload;
    memset(site,0,sizeof(site_t));
    site->ref = ss[0];
    site->dist = calloc(args->profile.nbins,sizeof(*site->dist));

    // ALT part
    ss = ++se;
    while ( *se && !isspace(*se) ) se++;
    site->alt = ss[0];

    return 0;
}
static void free_sites(void *payload)
{
    site_t *site = (site_t*)payload;
    free(site->dist);
}
static void batch_profile_destroy(args_t *args)
{
    if ( args->out_fh && bgzf_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    int i;
    for (i=0; i<args->nbams; i++) free(args->bams[i]);
    free(args->bams);
    if ( args->mplp ) mpileup_destroy(args->mplp);
    if ( args->sites_idx ) regidx_destroy(args->sites_idx);
    if ( args->sites_itr ) regitr_destroy(args->sites_itr);
    free(args->profile.good_calls.dist);
    free(args->profile.mean);
    free(args->profile.var);
    free(args);
}
static int batch_profile_init(args_t *args)
{
    args->bams = hts_readlist(args->aln_fname, 1, &args->nbams);
    if ( args->batch )
    {
        char *tmp;
        if ( !strncmp("k=",args->batch,2) )
        {
            // find out the number of required batches
            int k = strtol(args->batch+2,&tmp,10);
            if ( *tmp || k<=0 ) error("Error: could not parse --batch %s\n",args->batch);
            printf("# Number of required batches with %d files total and max %d files per batch:\n",args->nbams,k);
            printf("%.0f\n",ceil((double)args->nbams/k));
            batch_profile_destroy(args);
            return 1;
        }
        int ith = strtol(args->batch,&tmp,10);
        if ( !*tmp || *tmp!='/' || ith<=0 ) error("Error: could not parse --batch %s\n",args->batch);
        if ( ith > args->nbams ) error("Error: asked for %d-th batch in a list of %d files\n",ith,args->nbams);
        int nbatches = strtol(tmp+1,&tmp,10);
        if ( *tmp || nbatches<=0 ) error("Error: could not parse --batch %s\n",args->batch);
        if ( nbatches > args->nbams ) error("Error: cannot create %d batches from a list of %d files\n",nbatches,args->nbams);
        int nper_batch = ceil((double)args->nbams/nbatches);
        int isrc = (ith-1)*nper_batch;
        if ( isrc + nper_batch > args->nbams ) nper_batch = args->nbams - isrc;
        int i;
        for (i=0; i<isrc; i++)
        {
            free(args->bams[i]);
            args->bams[i] = NULL;
        }
        for (i=isrc+nper_batch; i<args->nbams; i++)
        {
            free(args->bams[i]);
            args->bams[i] = NULL;
        }
        if ( isrc ) memmove(args->bams,args->bams+isrc,(args->nbams-isrc)*sizeof(*args->bams));
        args->nbams = nper_batch;
    }

    int has_args = args->fasta_fname && args->sites_fname ? 1 : 0;
    if ( !has_args ) error("%s", usage_text());

    args->sites_idx = regidx_init(args->sites_fname,parse_sites,free_sites,sizeof(site_t),args);
    args->sites_itr = regitr_init(args->sites_idx);
    args->out_fh = bgzf_open(args->output_fname, args->output_type&FT_GZ ? "wg" : "wu");
    args->profile.good_calls.dist = calloc(args->profile.nbins,sizeof(*args->profile.good_calls.dist));
    return 0;
}
inline static int nn2bin(int nbin, int nref, int nalt)
{
    if ( !nalt && !nref ) return -1;
    if ( !nalt ) return 0;
    return (int)((double)(nbin-1)*nalt/(nref+nalt));
}

static int batch_profile_run1(args_t *args, char *aln_fname)
{
    // clean everything from the previous run
    if ( args->mplp ) mpileup_destroy(args->mplp);

    // init mpileup for next sample
    args->mplp = mpileup_alloc();
    mpileup_set(args->mplp, MAX_DP_PER_SAMPLE, 250);
    mpileup_set(args->mplp, MIN_MQ, 0);
    mpileup_set(args->mplp, MAX_BQ, 60);
    mpileup_set(args->mplp, DELTA_BQ, 30);
    mpileup_set(args->mplp, MIN_REALN_FRAC, 0.05);
    mpileup_set(args->mplp, MIN_REALN_DP, 2);
    mpileup_set(args->mplp, MAX_REALN_DP, 250);
    mpileup_set(args->mplp, MAX_REALN_LEN, 500);
    mpileup_set(args->mplp, SKIP_ANY_SET, BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);
    if ( args->use_bam_idx )
        mpileup_set(args->mplp, REGIONS_FNAME, args->sites_fname);
    else
        mpileup_set(args->mplp, TARGETS_FNAME, args->sites_fname);
    mpileup_set(args->mplp, LEGACY_MODE, 1);
    if ( mpileup_set(args->mplp, FASTA_REF, args->fasta_fname)!=0 ) error("Error: could not read the reference %s\n",args->fasta_fname);
    if ( mpileup_set(args->mplp, BAM, aln_fname)!=0 ) error("Error: could not reat %s\n",aln_fname);
    if ( mpileup_init(args->mplp)!=0 ) error("Error: could not initialize mpileup2\n");
    int nsmpl  = mpileup_get_val(args->mplp,int,N_SAMPLES);
    int *n_plp = mpileup_get_val(args->mplp,int*,N_READS);
    int ret,i,j,len;

    // process the entire bam
    while ( (ret=mpileup_next(args->mplp))==1 )
    {
        char *chr = mpileup_get_val(args->mplp,char*,CHROM);
        hts_pos_t pos = mpileup_get_val(args->mplp,int,POS);

        if ( !regidx_overlap(args->sites_idx, chr,pos,pos, args->sites_itr) )
        {
            fprintf(stderr,"No overlap: %s:%"PRIhts_pos"\n",chr,pos+1);
            continue;
        }

        bam_pileup1_t **plp = mpileup_get_val(args->mplp,bam_pileup1_t**,LEGACY_PILEUP);
        char *ref = mpileup_get(args->mplp,REF,&len);

        site_t *site = NULL;
        while ( regitr_overlap(args->sites_itr) )
        {
            site = &regitr_payload(args->sites_itr,site_t);
            if ( ref[0]==site->ref ) break;
            fprintf(stderr,"No ref match at %s:%"PRIhts_pos" ... %c vs %c\n",chr,pos+1,ref[0],site->ref);
            site = NULL;
        }
        if ( !site ) continue;
        for (i=0; i<nsmpl; i++)
        {
            int nref = 0, nalt = 0;
            for (j=0; j<n_plp[i]; j++)
            {
                const bam_pileup1_t *plp1 = plp[i] + j;
                int bi = bam_seqi(bam_get_seq(plp1->b), plp1->qpos);
                char bc = bi ? seq_nt16_str[bi] : 'x';
                if ( bc==site->ref ) nref++;
                else if ( bc==site->alt ) nalt++;
            }
            int ifreq = nn2bin(args->profile.nbins,nref,nalt);
            if ( ifreq<0 ) continue;
            if ( nref > nalt && calc_binom_two_sided(nref,nalt,0.5) <= args->binom_th )
            {
                args->profile.good_calls.nval++;
                args->profile.good_calls.dist[ifreq]++;
            }
            else
                site->no_good = 1;  // this site in this sample is not good for training: either low coverage or high VAF
            site->nval++;
            site->dist[ifreq]++;
        }
    }
    return 0;
}
static void batch_profile_set_mean_var(batch_t *batch)
{
    // this is to calculate the mean and variance for each VAF bin
    free(batch->mean);
    free(batch->var);
    batch->mean = calloc(batch->nbins,sizeof(*batch->mean));
    batch->var  = calloc(batch->nbins,sizeof(*batch->var));
    batch->nval = 0;

    // calculate profile from all good sites
    int i;
    regitr_t *itr = regitr_init(batch->idx);
    while ( regitr_loop(itr) )
    {
        site_t *site = &regitr_payload(itr,site_t);
        if ( site->no_good ) continue;

        // normalize the site and add to the mean and variance calculation
        double max_val = site->dist[0];
        for (i=1; i<batch->nbins; i++)
            if ( max_val < site->dist[i] ) max_val = site->dist[i];
        for (i=0; i<batch->nbins; i++)
        {
            double val = site->dist[i]/max_val;
            batch->mean[i] += val;
            batch->var[i]  += val*val;
        }
        batch->nval++;
    }
    assert(batch->nval);
    double min_nonzero_var = HUGE_VAL;
    for (i=0; i<batch->nbins; i++)
    {
        batch->mean[i] = batch->mean[i]/batch->nval;
        batch->var[i]  = batch->var[i]/batch->nval - batch->mean[i]*batch->mean[i];
        if ( batch->var[i]>0 && batch->var[i] < min_nonzero_var ) min_nonzero_var = batch->var[i];
    }
    // to avoid infinite scores, make sure we never see zero variance
    for (i=0; i<batch->nbins; i++)
        if ( batch->var[i]==0 ) batch->var[i] = min_nonzero_var;
    regitr_destroy(itr);
}
static double score_site(batch_t *batch, site_t *site)
{
    int i;
    double max_val = site->dist[0];
    for (i=1; i<batch->nbins; i++)
        if ( max_val < site->dist[i] ) max_val = site->dist[i];

    double score = 0;
    for (i=0; i<batch->nbins; i++)
    {
        double tmp = site->dist[i]/max_val - batch->mean[i];
        score += tmp*tmp/batch->var[i];
    }
    return score;
}
static int write_batch(args_t *args, batch_t *batch)
{
    // output the site profiles
    int i;
    kstring_t str = {0,0,0};
    if ( batch->idx )
    {
        regitr_t *itr = regitr_init(batch->idx);
        while ( regitr_loop(itr) )
        {
            site_t *site = &regitr_payload(itr,site_t);
            double score = score_site(batch,site);
            str.l = 0;
            ksprintf(&str,"SITE\t%s\t%d\t%c\t%c\t%d\t%e\t", itr->seq, itr->beg+1,site->ref,site->alt,(int)site->no_good,score);
            for (i=0; i<batch->nbins; i++) ksprintf(&str," %d",site->dist[i]);
            ksprintf(&str,"\n");
            if ( bgzf_write(args->out_fh,str.s,str.l)!=str.l ) error("Failed to write to %s\n",args->output_fname);
        }
        regitr_destroy(itr);
    }

    // output the good calls' profile
    if ( batch->good_calls.dist )
    {
        str.l = 0;
        ksprintf(&str,"COUNTS\t");
        for (i=0; i<batch->nbins; i++) ksprintf(&str," %d",batch->good_calls.dist[i]);
        ksprintf(&str,"\n");
        if ( bgzf_write(args->out_fh,str.s,str.l)!=str.l ) error("Failed to write to %s\n",args->output_fname);
    }

    // output the mean and variation at good calls
    if ( batch->mean )
    {
        str.l = 0;
        ksprintf(&str,"MEAN\t");
        for (i=0; i<batch->nbins; i++) ksprintf(&str," %e",batch->mean[i]);
        ksprintf(&str,"\n");
        ksprintf(&str,"VAR2\t");
        for (i=0; i<batch->nbins; i++) ksprintf(&str," %e",batch->var[i]);
        ksprintf(&str,"\n");
        if ( bgzf_write(args->out_fh,str.s,str.l)!=str.l ) error("Failed to write to %s\n",args->output_fname);
    }

    free(str.s);
    return 0;
}
static int batch_profile_run(args_t *args)
{
    // collect the profiles across all bams. This is the I/O intensive part
    int i;
    for (i=0; i<args->nbams; i++)
        batch_profile_run1(args, args->bams[i]);

    args->profile.idx = args->sites_idx;
    batch_profile_set_mean_var(&args->profile);
    write_batch(args, &args->profile);
    args->profile.idx = NULL;   // otherwise it would be destroyed twice

    return 0;
}
static uint32_t *parse_bins(const char *line, uint32_t *nbins)
{
    int i, n = 0;
    const char *ptr = line;
    while ( *ptr )
    {
        while ( *ptr && !isspace(*ptr) ) ptr++;
        n++;
        if ( *ptr ) ptr++;
    }
    ptr = line;
    uint32_t *bins = calloc(n,sizeof(*bins));
    for (i=0; i<n; i++)
    {
        char *tmp;
        bins[i] = strtol(ptr,&tmp,10);
        if ( *tmp && *tmp!=' ' ) error("Could not parse the DIST part of the line: %s\n",line);
        ptr = tmp+1;
    }
    *nbins = n;
    return bins;
}
static int parse_batch(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    // parses lines like this:
    //  SITE	chr	3	A	C	1	7.194245e-02	 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    batch_t *batch = (batch_t*) usr;

    uint32_t nbins, *bins;
    if ( !strncmp(line,"COUNTS\t",7) )
    {
        int i;
        batch->good_calls.dist = parse_bins(line+8, &nbins);
        if ( batch->nbins && batch->nbins!=nbins ) error("Different number of bins, %d vs %d: %s\n",batch->nbins,nbins,line+8);
        for (i=0; i<nbins; i++)
            batch->good_calls.nval += batch->good_calls.dist[i];
        return -1;
    }
    if ( strncmp(line,"SITE\t",5) ) return -1;     // skip

    // CHR part
    char *ss = (char*) line + 5;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -2;      // unexpected format

    char *se = ss;
    while ( *se && !isspace(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se ) error("Could not parse the CHR part of the line: %s\n",line);

    // POS part
    ss = se+1;
    *beg = strtod(ss, &se);
    if ( ss==se ) error("Could not parse the POS part of the line: %s\n",line);
    (*beg)--;

    // REF part and REF length
    ss = ++se;
    while ( *se && !isspace(*se) ) se++;
    *end = *beg + se-ss-1;

    site_t *site = (site_t*)payload;
    memset(site,0,sizeof(site_t));
    site->ref = ss[0];

    // ALT part
    ss = ++se;
    while ( *se && !isspace(*se) ) se++;
    site->alt = ss[0];

    // NO_GOOD part
    if ( *se!='\t' ) error("Could not parse the NO_GOOD part of the line: %s\n",line);
    se++;
    if ( *se!='1' && *se!='0' ) error("Could not parse the NO_GOOD part of the line: %s\n",line);
    site->no_good = *se=='0' ? 0 : 1;
    se++;
    while ( *se && isspace(*se) ) se++;

    // skip the SCORE part
    ss = se;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se ) error("Could not parse the SCORE part of the line: %s\n",line);
    while ( *se && isspace(*se) ) se++;
    if ( !*se ) error("Could not parse the SCORE part of the line: %s\n",line);

    // read the PROFILE part
    bins = parse_bins(se, &nbins);
    if ( !batch->nbins ) batch->nbins = nbins;
    if ( batch->nbins!=nbins ) error("Different number of bins, %d vs %d: %s\n",batch->nbins,nbins,line);
    site->dist = bins;

    return 0;
}
static batch_t *batch_read(char *fname)
{
    batch_t *batch = calloc(1,sizeof(batch_t));
    batch->idx = regidx_init(fname,parse_batch,free_sites,sizeof(site_t),batch);
    batch->fname = strdup(fname);
    return batch;
}
static int batch_merge(batch_t *tgt, batch_t *src)
{
    if ( tgt->nbins!=src->nbins ) error("Different bin size in %s and %s\n",tgt->fname,src->fname);
    regitr_t *src_itr = regitr_init(src->idx);
    regitr_t *tgt_itr = regitr_init(tgt->idx);
    int i;
    regitr_reset(src->idx, src_itr);
    while ( regitr_loop(src_itr) )
    {
        site_t *src_site = &regitr_payload(src_itr,site_t);
        if ( !regidx_overlap(tgt->idx, src_itr->seq,src_itr->beg,src_itr->end, tgt_itr) ) error("uh: fixme\n");
        while ( regitr_overlap(tgt_itr) )
        {
            site_t *tgt_site = &regitr_payload(tgt_itr,site_t);
            if ( src_site->no_good ) tgt_site->no_good = 1;
            for (i=0; i<tgt->nbins; i++) tgt_site->dist[i] += src_site->dist[i];
            tgt_site->nval += src_site->nval;
            break;
        }
    }
    regitr_destroy(src_itr);
    regitr_destroy(tgt_itr);
    tgt->nval += src->nval;
    for (i=0; i<tgt->nbins; i++)
        tgt->good_calls.dist[i] += src->good_calls.dist[i];
    return 0;
}
static void batch_destroy(batch_t *batch)
{
    free(batch->good_calls.dist);
    if ( batch->idx ) regidx_destroy(batch->idx);
    free(batch->fname);
    free(batch->mean);
    free(batch->var);
    free(batch);
}
static int merge(args_t *args)
{
    args->out_fh = bgzf_open(args->output_fname, args->output_type&FT_GZ ? "wg" : "wu");

    batch_t *batch = NULL, *tmp;
    int i,nfile = 0;
    char **file = hts_readlist(args->batch_fname, 1, &nfile);
    if ( !file ) error("No files to merge: %s\n",args->batch_fname);
    for (i=0; i<nfile; i++)
    {
        tmp = batch_read(file[i]);
        free(file[i]);
        if ( !i ) { batch = tmp; continue; }
        batch_merge(batch,tmp);
        batch_destroy(tmp);
    }

    batch_profile_set_mean_var(batch);
    write_batch(args,batch);

    if ( args->out_fh && bgzf_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    batch_destroy(batch);
    free(file);
    free(args);
    return 0;
}
int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->record_cmd_line = 1;
    args->clevel = -1;
    args->binom_th = 1e-4;
    args->profile.nbins = 20;
    static struct option loptions[] =
    {
        {"batch",required_argument,NULL,'b'},
        {"merge-batches",required_argument,NULL,'m'},
        {"use-index",required_argument,NULL,'i'},
        {"nbins",required_argument,NULL,'n'},
        {"binom-th",required_argument,NULL,'B'},
        {"fasta-ref",required_argument,NULL,'f'},
        {"alns",required_argument,NULL,'a'},
        {"sites",required_argument,NULL,'s'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "o:O:s:t:a:f:B:n:ib:m:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'B':
                args->binom_th = strtod(optarg,&tmp);
                if ( *tmp || args->binom_th<=0 || args->binom_th>1 ) error("Could not parse argument: --binom-th %s, expecting value (0,1]\n", optarg);
                break;
            case 'n':
                args->profile.nbins = strtol(optarg,&tmp,10);
                if ( *tmp || args->profile.nbins<10 ) error("Could not parse argument: --nbins %s; the minimum value is 10\n", optarg);
                break;
            case 'i': args->use_bam_idx = 1; break;
            case 'f': args->fasta_fname = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 't': args->output_type = FT_TAB_TEXT; break;
                          case 'z': args->output_type = FT_GZ; break;
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      };
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 's': args->sites_fname = optarg; break;
            case 'a': args->aln_fname = optarg; break;
            case 'b': args->batch = optarg; break;
            case 'm': args->batch_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    int do_merge = args->batch_fname ? 1 : 0;
    int do_profile = args->aln_fname ? 1 : 0;
    if ( !do_merge && !do_profile ) error("%s", usage_text());
    if ( do_merge && do_profile ) error("%s", usage_text());

    if ( do_profile )
    {
        if ( batch_profile_init(args) > 0 ) return 0;
        batch_profile_run(args);
        batch_profile_destroy(args);
    }
    if ( do_merge ) merge(args);

    return 0;
}
