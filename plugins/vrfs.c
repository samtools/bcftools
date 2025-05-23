/* The MIT License

   Copyright (c) 2024-2025 Genome Research Ltd.

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
#include <sys/time.h>
#include "mpileup2/mpileup.h"
#include "bcftools.h"
#include "regidx.h"


// Variance at good sites as found in real data (~2000 samples). Originally, the code would attempt to
// identify good sites with an automatic heuristic, but that did not always work. A more robust way is
// to use a pre-computed profile and allow the user to override (todo)
#define N_BINS 20
const double var2[N_BINS] = {1,1.441327e-03,4.382657e-05,6.160600e-07,1.414270e-08,2.828540e-09,2.357117e-09,
    2.020386e-09,1.767838e-09,1.571411e-09,1.414270e-09,1.285700e-09,1.178558e-09,1.087900e-09,1.010193e-09,
    9.428468e-10,8.839188e-10,8.319236e-10,7.857056e-10,7.443527e-10};

#define VAR2_DATA 1
#define VAR2_FILE 2
#define VAR2_HC   3

// VAF profile across many samples at a single site. Used as a regidx payload
typedef struct
{
    char *ref, *alt;    // note indels are preserved, but internally all treated as a single type
    uint32_t nval,      // number of values
            *dist;      // histogram of VAF frequencies
}
site_t;

typedef struct
{
    char *fname;
    int nbins;              // number of VAF bins
    int nval;               // number of values contributing to the computed profile
    double *mean, *var2;    // mean and variance of each VAF bin across all sites indexed by regidx
    regidx_t *idx;          // sites in the batch; all batches should have the same set of sites
}
batch_t;

typedef struct
{
    char *seq;
    int beg;
    double score;
    site_t *site;
}
prn_site_t;

typedef struct
{
    int argc, output_type, record_cmd_line, clevel, use_bam_idx, verbose;
    char **argv, *output_fname, *sites_fname, *aln_fname, *fasta_fname, *batch_fname, *batch;
    char **bams, **batch_fnames, *recalc_type_str;
    int recalc_type;    // one of VAR2_* types
    int nbams, nbatch_fnames;
    int min_dp;    // minimum read depth to consider
    prn_site_t *prn_site_buf;
    int nprn_site_buf;
    BGZF *out_fh;
    regidx_t *sites_idx;
    regitr_t *sites_itr;
    mpileup_t *mplp;
    batch_t profile;
}
args_t;

const char *about(void)
{
    return "Localised assessment of sequencing artefacts, estimate site noisiness (variant read frequency score)\n";
}
static const char *usage_text(void)
{
    return
        "\n"
        "About: Assess site noisiness (variant read frequency score) from a large number of unaffected parental samples\n"
        "Usage: bcftools +vrfs [OPTIONS]\n"
        "\n"
        "Required options:\n"
        "   -a, --alns FILE               List of BAM/CRAM files\n"
        "   -f, --fasta-ref FILE          Reference file in fasta format\n"
        "   -s, --sites FILE              A tab-delimited file name of sites to assess (chr,pos,ref,alt)\n"
        "\n"
        "Model options:\n"
        "   -d, --min-depth INT           Minimum depth of calls to include in the profile [10]\n"
        "   -n, --nbins INT               Number of VAF bins [20]\n"
        "   -r, --recalc TYPE             Recalculate scores based on provided variances [hc]\n"
        "                                   data .. as observed in the data\n"
        "                                   file .. read profile from a file, one value per line (-r file:/path/to/file.txt)\n"
        "                                   hc   .. hard-coded profile\n"
        "\n"
        "Other options:\n"
        "   -b, --batch I/N               Run I-th batch out of N, 1-based\n"
        "   -i, --use-index               Use index to jump the alignments, rather than stream (faster with few sites)\n"
        "   -m, --merge-batches FILE      Merge files produced with -b, --batch, FILE is a list of batch files\n"
        "   -M, --merge-files FILE ...    Same as -m, FILE is one or more batch files\n"
        "   -o, --output FILE             Output file name [stdout]\n"
        "   -O, --output-type t|z[0-9]    t/z: un/compressed text file, 0-9: compression level [t]\n"
        "   -v, --verbosity INT           Verbosity level, eg can print the elapsed and estimated total running time\n"
        "\n"
        "Example:\n"
        "   # Typical run\n"
        "   bcftools +vrfs -f ref.fa -a bams.txt -s sites.txt -o scores.txt\n"
        "\n"
        "   # Run in batches. Let's say one batch can have at most 3 bams and there are 5 bams in total\n"
        "   #   1) find out the number of required batches with `--batch k=3` (the program outputs: 2)\n"
        "   #   2) run all batches individually with `--batch 1/2` and `--batch 2/2`\n"
        "   #   3) create a list of outputs and merge with `--merge-batches list.txt`\n"
        "   bcftools +vrfs -a bams.txt --batch k=3     # prints 2\n"
        "   bcftools +vrfs -f ref.fa -a bams.txt -s sites.txt -o scores1.txt --batch 1/2\n"
        "   bcftools +vrfs -f ref.fa -a bams.txt -s sites.txt -o scores2.txt --batch 2/2\n"
        "   bcftools +vrfs --merge-batches list.txt\n"
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
    while ( *se && isspace(*se) ) se++;
    ss = se;
    while ( *se && !isspace(*se) ) se++;
    int ref_len = se - ss;
    if ( !ref_len ) error("Could not parse the REF part of the line: %s\n",line);
    *end = *beg;    // we are interested in overlaps at the POS only, not variant length

    site_t *site = (site_t*)payload;
    memset(site,0,sizeof(site_t));
    site->ref = malloc(ref_len+1);
    strncpy(site->ref,ss,ref_len);
    site->ref[ref_len] = 0;
    site->dist = calloc(args->profile.nbins,sizeof(*site->dist));

    // ALT part
    while ( *se && isspace(*se) ) se++;
    ss = se;
    while ( *se && !isspace(*se) ) se++;
    int alt_len = se - ss;
    if ( !alt_len ) error("Could not parse the ALT part of the line: %s\n",line);
    site->alt = malloc(alt_len+1);
    strncpy(site->alt,ss,alt_len);
    site->alt[alt_len] = 0;

    if ( ref_len==alt_len )     // whatever this is, treat it as SNV, not MNP, should have been split
        site->ref[1] = site->alt[1] = 0;
    return 0;
}
static void free_sites(void *payload)
{
    site_t *site = (site_t*)payload;
    free(site->ref);
    free(site->alt);
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
    free(args->profile.mean);
    free(args->profile.var2);
    free(args->prn_site_buf);
    args->prn_site_buf  = NULL;
    args->nprn_site_buf = 0;
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
        if ( ith > nbatches ) error("Error: the batch index is outside the permitted range [1,%d]\n",nbatches);
        if ( nbatches > args->nbams ) error("Error: cannot create %d batches from a list of %d files\n",nbatches,args->nbams);
        int nper_batch = ceil((double)args->nbams/nbatches);
        int isrc = (ith-1)*nper_batch;
        if ( isrc > args->nbams )
        {
            args->out_fh = bgzf_open(args->output_fname, args->output_type&FT_GZ ? "wg" : "wu");
            kstring_t str = {0,0,0};
            ksprintf(&str,
                "# This is the %d-th chunk out of %d requested. As you can see, it is empty: there are %d files per batch\n"
                "# and %d files in total. Don't worry, it can still be used with the options --merge-batches and --merge-files.\n",
                ith,nbatches,nper_batch,args->nbams);
            if ( bgzf_write(args->out_fh,str.s,str.l)!=str.l ) error("Failed to write to %s\n",args->output_fname);
            free(str.s);
            batch_profile_destroy(args);
            return 1;
        }
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

    int ret,i,j,len;
    if ( args->use_bam_idx )
        ret = mpileup_set(args->mplp, REGIONS_FNAME, args->sites_fname);
    else
    {
        ret = mpileup_set(args->mplp, TARGETS_FNAME, args->sites_fname);
        if ( args->verbose )
            fprintf(stderr,"Note: the -i, --use-index option is not given, streaming the alignment files\n");
    }
    if ( ret ) error("Error: could not initialize site list %s\n",args->sites_fname);

    mpileup_set(args->mplp, LEGACY_MODE, 1);
    if ( mpileup_set(args->mplp, FASTA_REF, args->fasta_fname)!=0 ) error("Error: could not read the reference %s\n",args->fasta_fname);
    if ( mpileup_set(args->mplp, BAM, aln_fname)!=0 ) error("Error: could not reat %s\n",aln_fname);
    if ( mpileup_init(args->mplp)!=0 ) error("Error: could not initialize mpileup2\n");
    int nsmpl  = mpileup_get_val(args->mplp,int,N_SAMPLES);
    int *n_plp = mpileup_get_val(args->mplp,int*,N_READS);

    if ( args->verbose>=3 )
    {
        static int hdr_printed = 0;
        if ( !hdr_printed )
            fprintf(stderr,"# [1]DEBUG_SITE\t[2]file\t[3]chr:pos\t[4]REF\t[5]ALT\t[6]nREF\t[7]nALT\t[8]iAF_Bin\n");
        hdr_printed = 1;
    }

    // process the entire bam
    while ( (ret=mpileup_next(args->mplp))==1 )
    {
        char *chr = mpileup_get_val(args->mplp,char*,CHROM);
        hts_pos_t pos = mpileup_get_val(args->mplp,hts_pos_t,POS);

        if ( !regidx_overlap(args->sites_idx, chr,pos,pos, args->sites_itr) )
        {
            fprintf(stderr,"No overlap: %s:%"PRIhts_pos"\n",chr,pos+1);
            continue;
        }

        bam_pileup1_t **plp = mpileup_get_val(args->mplp,bam_pileup1_t**,LEGACY_PILEUP);
        char *ref = mpileup_get(args->mplp,REF,&len);

        // There can be duplicate positions with different ALT alleles, we make each of these accessible
        // through ialt2site array, the index uses the encoding 0/1/2/3 for A/C/G/T and 4 for I/D
        site_t *ialt2site[5] = {NULL,NULL,NULL,NULL,NULL};
        int has_site = 0;
        while ( regitr_overlap(args->sites_itr) )
        {
            site_t *site = &regitr_payload(args->sites_itr,site_t);
            if ( ref[0]!=site->ref[0] )
            {
                fprintf(stderr,"No ref match at %s:%"PRIhts_pos" ... %c vs \"%s\"\n",chr,pos+1,ref[0],site->ref);
                continue;
            }
            int ialt = site->ref[1] || site->alt[1] ? 4 : seq_nt16_int[seq_nt16_table[(int)site->alt[0]]];  // indel or SNV?
            if ( !ialt2site[ialt] ) ialt2site[ialt] = site;     // keep only one duplicate position of the same type
            has_site = 1;
        }
        if ( !has_site ) continue;
        for (i=0; i<nsmpl; i++)
        {
            // Collect counts: the number of reads in total, alt alleles total, and each alt separately
            int ntot = 0, nalt[5] = {0,0,0,0,0};
            for (j=0; j<n_plp[i]; j++)
            {
                const bam_pileup1_t *plp1 = plp[i] + j;
                int ialt = -1;
                if ( plp1->indel ) ialt = 4;
                else
                {
                    int bi = bam_seqi(bam_get_seq(plp1->b), plp1->qpos);
                    assert( bi );   // when does this happen?
                    char bc = seq_nt16_str[bi];
                    if ( bc!=ref[0] )
                        ialt = seq_nt16_int[seq_nt16_table[(int)bc]];
                    else
                        ntot++;
                }
                if ( ialt>=0 )
                {
                    nalt[ialt]++;
                    ntot++;
                }
            }
            if ( ntot < args->min_dp ) continue;

            // Increment site counters
            for (j=0; j<5; j++)
            {
                site_t *site = ialt2site[j];
                if ( !site ) continue;
                int ifreq = nn2bin(args->profile.nbins,ntot-nalt[j],nalt[j]);
                assert( ifreq>=0 );
                site->nval++;
                site->dist[ifreq]++;
                if ( args->verbose >= 3 )
                    fprintf(stderr,"DEBUG_SITE\t%s\t%s:%"PRIhts_pos"\t%s\t%s\t%d\t%d\t%d\n",aln_fname,chr,pos+1,site->ref,site->alt,ntot-nalt[j],nalt[j],ifreq);
            }
        }
    }
    return 0;
}

/*
    In the mode "hc" and "file" the function returns NULL, otherwise returns a newly allocated array
    of size batch->nbins, to be filled with VAR2 values collected from the data.
    The modes "hc" and "file" may need rescaling, if the number of bins is different from --nbins.
*/
static double *init_var2(args_t *args, batch_t *batch)
{
    assert( !batch->var2 );

    double *ori = NULL;
    int i, nnew = batch->nbins , nori = 0;
    if ( !strncmp(args->recalc_type_str,"file:",5)  )
    {
        int n;
        char *tmp, **list = hts_readlist(args->recalc_type_str+5, 1, &n);
        if ( !list || !n ) error("Error: could not parse %s\n",args->recalc_type_str+5);
        ori = malloc(sizeof(*batch->var2)*n);
        for (i=0; i<n; i++)
        {
            ori[i] = strtod(list[i], &tmp);
            if ( tmp==list[i] ) error("Error: could not parse the %d-th value in %s\n",i+1,args->recalc_type_str+5);
            free(list[i]);
        }
        free(list);
        nori = n;
        args->recalc_type = VAR2_FILE;
    }
    else if ( !strcmp(args->recalc_type_str,"hc") )
    {
        ori = malloc(sizeof(*ori)*N_BINS);
        for (i=0; i<N_BINS; i++) ori[i] = var2[i];
        nori = N_BINS;
        args->recalc_type = VAR2_HC;
    }

    if ( ori && nori==nnew )
    {
        batch->var2 = ori;
        return NULL;
    }
    if ( ori && nori!=nnew )
    {
        // intrapolate data points
        double *tmp = malloc(sizeof(*tmp)*nnew);

        tmp[0] = ori[0];
        tmp[nnew - 1] = ori[nori - 1];

        double dx = (double)1/(nnew-1);
        double DX = (double)1/(nori-1);
        int J=1;
        for (i=1; i<nnew - 1; i++)
        {
            double x = i*dx;
            while ( J*DX < x ) J++;
            assert( J < nori );
            double X = (J-1)*DX;
            assert( X <= x );
            if ( x==X )
                tmp[i] = ori[J-1];
            else
                tmp[i] = ori[J-1] + (ori[J]-ori[J-1])*(x-X)/DX;
        }
        free(ori);
        batch->var2 = tmp;
        return NULL;
    }

    if ( !strcmp(args->recalc_type_str,"data") ) args->recalc_type = VAR2_DATA;
    else error("Error: the mode --recalc %s is not recognized\n",args->recalc_type_str);

    batch->var2 = calloc(batch->nbins,sizeof(*batch->var2));
    return batch->var2;
}
static void batch_profile_set_mean_var2(args_t *args, batch_t *batch)
{
    // this is to calculate the mean and variance for each VAF bin
    free(batch->mean);
    batch->mean = calloc(batch->nbins,sizeof(*batch->mean));

    // if a pointer is returned, collect VAR2 from the data
    double *var2_ptr = init_var2(args, batch);
    batch->nval = 0;

    // calculate profile from all sites
    int i;
    regitr_t *itr = regitr_init(batch->idx);
    while ( regitr_loop(itr) )
    {
        site_t *site = &regitr_payload(itr,site_t);
        if ( !site->nval ) continue;

        // normalize the site and add to the mean and variance calculation
        double max_val = site->dist[0];
        for (i=0; i<batch->nbins; i++)
            if ( max_val < site->dist[i] ) max_val = site->dist[i];
        for (i=0; i<batch->nbins; i++)
        {
            double val = site->dist[i]/max_val;
            batch->mean[i] += val;
            if ( var2_ptr ) var2_ptr[i] += val*val;
        }
        batch->nval++;
    }
    for (i=0; i<batch->nbins; i++)
        batch->mean[i] = batch->mean[i]/batch->nval;

    if ( var2_ptr && batch->nval && args->recalc_type==VAR2_DATA )
    {
        double min_nonzero_var2 = 1;
        for (i=0; i<batch->nbins; i++)
        {
            var2_ptr[i] = var2_ptr[i]/batch->nval - batch->mean[i]*batch->mean[i];
            if ( var2_ptr[i]>0 && var2_ptr[i] < min_nonzero_var2 ) min_nonzero_var2 = var2_ptr[i];
        }
        // to avoid infinite scores, make sure we never see zero variance,
        // but make it ever decreasing to penalize higher VAF bins
        for (i=0; i<batch->nbins; i++)
            if ( var2_ptr[i]==0 ) var2_ptr[i] = min_nonzero_var2/(i?i+1:1);
    }
    regitr_destroy(itr);
}
static double score_site(batch_t *batch, site_t *site)
{
    int i;
    double max_val = site->dist[0];
    for (i=0; i<batch->nbins; i++)
        if ( max_val < site->dist[i] ) max_val = site->dist[i];
    if ( !max_val ) max_val = 1;

    // This naturally counts excesses in non-zero VAF bins, cleaner sites are not penalized.
    double score = 0;
    for (i=1; i<batch->nbins; i++)
    {
        double tmp = site->dist[i]/max_val;
        score += tmp*tmp / batch->var2[i];
    }
    return 10*log(1 + score);       // multiply by 10 for a better range of [0,250) or so
}


// How to treat duplicate positions: one alt may look good, but we want to consider this
// as a site - one bad alt spoils it for all. With indels, it's even more complicated,
// we cannot detect the alleles quite accurately (reference bias, realignment issues),
// therefore we decided to see just a generic 'indel', not individual indel types. We
// are aware this can lead to missed calls at multiallelic sites, but reducing FDR is
// a lesser of the two evils.
static void print_buffered_sites(args_t *args, batch_t *batch, prn_site_t *buf, int nbuf, kstring_t *str)
{
    double max_score = 0;
    int i,j;
    for (i=0; i<nbuf; i++)
        if ( max_score < buf[i].score ) max_score = buf[i].score;

    // For technical reasons we keep stats on a generic indel, specifically, the first indel record.
    // We use its distribution and score for all indel types
    prn_site_t *fst_indel = NULL;
    for (i=0; i<nbuf; i++)
    {
        // ad-hoc rule: increase the lower score by 75% of the difference to the max score
        double score = (max_score - buf[i].score)*0.75 + buf[i].score;

        site_t *site = buf[i].site;
        uint32_t *dist = buf[i].site->dist;

        if ( site->ref[1] || site->alt[1] )
        {
            if ( !fst_indel ) fst_indel = &buf[i];
            dist  = fst_indel->site->dist;
            score = fst_indel->score;
        }

        str->l = 0;
        ksprintf(str,"SITE\t%s\t%d\t%s\t%s\t%e\t", buf[i].seq, buf[i].beg+1, site->ref, site->alt, score);
        for (j=0; j<batch->nbins; j++) ksprintf(str,"%s%d",j==0?"":"-",dist[j]);
        ksprintf(str,"\n");
        if ( bgzf_write(args->out_fh,str->s,str->l)!=str->l ) error("Failed to write to %s\n",args->output_fname);
    }
}

static int write_batch(args_t *args, batch_t *batch)
{
    // output the site profiles
    int i;
    kstring_t str = {0,0,0};
    if ( batch->idx )
    {
        // Duplicate positions need to be re-scored, when one alt looks good while another bad, it is still a noisy site.
        int nbuf = 0;
        regitr_t *itr = regitr_init(batch->idx);
        while ( regitr_loop(itr) )
        {
            site_t *site = &regitr_payload(itr,site_t);
            if ( nbuf && (args->prn_site_buf[nbuf-1].seq!=itr->seq || args->prn_site_buf[nbuf-1].beg!=itr->beg) )
            {
                print_buffered_sites(args,batch,args->prn_site_buf,nbuf,&str);
                nbuf = 0;
            }
            nbuf++;
            hts_resize(prn_site_t, nbuf, &args->nprn_site_buf, &args->prn_site_buf, 0);
            prn_site_t *tmp = &args->prn_site_buf[nbuf-1];
            tmp->beg = itr->beg;
            tmp->seq = itr->seq;
            tmp->score = score_site(batch,site);
            tmp->site  = site;
        }
        print_buffered_sites(args,batch,args->prn_site_buf,nbuf,&str);
        regitr_destroy(itr);
    }

    // output the mean and variance
    if ( batch->mean )
    {
        str.l = 0;
        ksprintf(&str,"MEAN\t");
        for (i=0; i<batch->nbins; i++) ksprintf(&str," %e",batch->mean[i]);
        ksprintf(&str,"\n");
        ksprintf(&str,"VAR2\t");
        for (i=0; i<batch->nbins; i++) ksprintf(&str," %e",batch->var2[i]);
        ksprintf(&str,"\n");
        if ( bgzf_write(args->out_fh,str.s,str.l)!=str.l ) error("Failed to write to %s\n",args->output_fname);
    }

    free(str.s);
    return 0;
}
static void ksprint_time(kstring_t *str, double delta)
{
    if ( delta > 60*60*24 )
    {
        ksprintf(str,"%.0fd",ceil(delta/60./60./24));
        return;
    }
    if ( delta > 60*60 )
    {
        ksprintf(str,"%.0fh",ceil(delta/60./60.));
        return;
    }
    if ( delta > 60 )
    {
        ksprintf(str,"%.0fm",ceil(delta/60.));
        return;
    }
    ksprintf(str,"%.0fs",delta>0 ? ceil(delta) : 1);
}
static int batch_profile_run(args_t *args)
{
    kstring_t str = {0,0,0};
    struct timeval t0, t1;
    gettimeofday(&t0, NULL);
    double delta_prev = 0;
    int i;

    // collect the profiles across all bams. This is the I/O intensive part
    for (i=0; i<args->nbams; i++)
    {
        batch_profile_run1(args, args->bams[i]);

        if ( args->verbose )
        {
            // Report time
            gettimeofday(&t1, NULL);
            double delta = (t1.tv_sec - t0.tv_sec) * 1e6 + (t1.tv_usec - t0.tv_usec);
            double avg = delta / (i+1);
            double tot = avg * args->nbams;
            str.l = 0;
            ksprintf(&str,"Time required to process %s .. ",args->bams[i]);
            ksprint_time(&str,(delta-delta_prev)/1e6);
            kputs("  (avg/elapsed/est tot/eta: ",&str);
            ksprint_time(&str,avg/1e6);
            kputc('/',&str);
            ksprint_time(&str,delta/1e6);
            kputc('/',&str);
            ksprint_time(&str,tot/1e6);
            kputc('/',&str);
            ksprint_time(&str,(tot-delta)/1e6);
            fprintf(stderr,"%s)\n",str.s);
            delta_prev = delta;
        }
    }
    free(str.s);

    args->profile.idx = args->sites_idx;
    batch_profile_set_mean_var2(args, &args->profile);
    write_batch(args, &args->profile);
    args->profile.idx = NULL;   // otherwise it would be destroyed twice

    return 0;
}
static uint32_t *parse_bins(const char *line, uint32_t *nbins, uint32_t *nval)
{
    if ( nval ) *nval = 0;
    int i, n = 0;
    const char *ptr = line;
    while ( *ptr )
    {
        while ( *ptr && *ptr!='-' ) ptr++;
        n++;
        if ( *ptr ) ptr++;
    }
    ptr = line;
    uint32_t *bins = calloc(n,sizeof(*bins));
    for (i=0; i<n; i++)
    {
        char *tmp;
        bins[i] = strtol(ptr,&tmp,10);
        if ( *tmp && *tmp!='-' ) error("Could not parse the DIST part of the line: %s\n",line);
        if ( nval ) *nval += bins[i];
        ptr = tmp+1;
    }
    *nbins = n;
    return bins;
}
static double *parse_float_array(const char *line, int *narray)
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
    double *array = calloc(n,sizeof(*array));
    for (i=0; i<n; i++)
    {
        char *tmp;
        array[i] = strtod(ptr,&tmp);
        if ( *tmp && !isspace(*tmp) ) error("Could not parse the float array: %s\n",line);
        ptr = tmp+1;
    }
    *narray = n;
    return array;
}
static int parse_batch(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    // parses lines like this:
    //  SITE	chr	3	A	C	1	7.194245e-02	 1-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0
    batch_t *batch = (batch_t*) usr;

    if ( !strncmp(line,"MEAN\t",5) )
    {
        int nbins;
        batch->mean = parse_float_array(line+6, &nbins);
        if ( batch->nbins && batch->nbins!=nbins ) error("Different number of bins, %d vs %d: %s\n",batch->nbins,nbins,line);
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

    // REF part
    ss = ++se;
    while ( *se && !isspace(*se) ) se++;
    int ref_len = se - ss;
    *end = *beg;

    site_t *site = (site_t*)payload;
    memset(site,0,sizeof(site_t));
    site->ref = malloc(ref_len+1);
    strncpy(site->ref,ss,ref_len);
    site->ref[ref_len] = 0;

    // ALT part
    ss = ++se;
    while ( *se && !isspace(*se) ) se++;
    int alt_len = se - ss;
    site->alt = malloc(alt_len+1);
    strncpy(site->alt,ss,alt_len);
    site->alt[alt_len] = 0;

    // skip the SCORE part
    while ( *se && isspace(*se) ) se++;
    ss = se;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se ) error("Could not parse the SCORE part of the line: %s\n",line);
    while ( *se && isspace(*se) ) se++;
    if ( !*se ) error("Could not parse the SCORE part of the line: %s\n",line);

    // read the PROFILE part
    uint32_t nbins, nval, *bins = parse_bins(se, &nbins,&nval);
    if ( !batch->nbins ) batch->nbins = nbins;
    if ( batch->nbins!=nbins ) error("Different number of bins, %d vs %d: %s\n",batch->nbins,nbins,line);
    site->dist = bins;
    site->nval = nval;

    return 0;
}
static batch_t *batch_read(char *fname)
{
    batch_t *batch = calloc(1,sizeof(batch_t));
    batch->idx = regidx_init(fname,parse_batch,free_sites,sizeof(site_t),batch);
    if ( !batch->idx ) error("Could not read the batch file: %s\n",fname);
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
            if ( strcmp(src_site->ref,tgt_site->ref) ) continue;
            if ( strcmp(src_site->alt,tgt_site->alt) ) continue;
            for (i=0; i<tgt->nbins; i++) tgt_site->dist[i] += src_site->dist[i];
            tgt_site->nval += src_site->nval;
            break;
        }
    }
    regitr_destroy(src_itr);
    regitr_destroy(tgt_itr);
    tgt->nval += src->nval;
    return 0;
}
static void batch_destroy(batch_t *batch)
{
    if ( batch->idx ) regidx_destroy(batch->idx);
    free(batch->fname);
    free(batch->mean);
    free(batch->var2);
    free(batch);
}
static void merge_add_batch(args_t *args, const char *fname)
{
    args->nbatch_fnames++;
    args->batch_fnames = (char**)realloc(args->batch_fnames,sizeof(*args->batch_fnames)*args->nbatch_fnames);
    args->batch_fnames[args->nbatch_fnames - 1] = strdup(fname);
}
static int merge(args_t *args)
{
    args->out_fh = bgzf_open(args->output_fname, args->output_type&FT_GZ ? "wg" : "wu");

    if ( !args->nbatch_fnames )
    {
        args->batch_fnames = hts_readlist(args->batch_fname, 1, &args->nbatch_fnames);
        if ( !args->batch_fnames ) error("Could not read the file: %s\n",args->batch_fname);
    }

    batch_t *batch = NULL, *tmp;
    int i;
    for (i=0; i<args->nbatch_fnames; i++)
    {
        tmp = batch_read(args->batch_fnames[i]);
        free(args->batch_fnames[i]);
        if ( !tmp->nbins ) // no data
        {
            batch_destroy(tmp);
            continue;
        }
        if ( !batch )
        {
            batch = tmp;
            continue;
        }
        batch_merge(batch,tmp);
        batch_destroy(tmp);
    }
    if ( !batch ) error("Error: failed to merge the files, no usable data found\n");

    batch_profile_set_mean_var2(args,batch);
    write_batch(args,batch);

    free(args->batch_fnames);
    batch_destroy(batch);
    batch_profile_destroy(args);
    return 0;
}
int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type = -1;
    args->record_cmd_line = 1;
    args->clevel = -1;
    args->min_dp = 10;
    args->profile.nbins = N_BINS;
    args->recalc_type_str = "hc";
    static struct option loptions[] =
    {
        {"recalc",required_argument,NULL,'r'},
        {"batch",required_argument,NULL,'b'},
        {"merge-batches",required_argument,NULL,'m'},
        {"merge-files",required_argument,NULL,'M'},
        {"use-index",required_argument,NULL,'i'},
        {"nbins",required_argument,NULL,'n'},
        {"min-depth",required_argument,NULL,'d'},
        {"fasta-ref",required_argument,NULL,'f'},
        {"alns",required_argument,NULL,'a'},
        {"sites",required_argument,NULL,'s'},
        {"verbose",optional_argument,NULL,'v'},
        {"verbosity",optional_argument,NULL,'v'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "o:O:s:t:a:f:d:n:ib:m:M:v::r:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'd':
                args->min_dp = strtol(optarg,&tmp,10);
                if ( *tmp || args->min_dp<0 ) error("Could not parse argument: --min-depth %s\n", optarg);
                break;
            case 'n':
                args->profile.nbins = strtol(optarg,&tmp,10);
                if ( *tmp || args->profile.nbins<10 ) error("Could not parse argument: --nbins %s; the minimum value is 10\n", optarg);
                break;
            case 'i': args->use_bam_idx = 1; break;
            case 'r': args->recalc_type_str = optarg; break;
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
            case 'v':
                if (!optarg) args->verbose++;
                else
                {
                    args->verbose = strtol(optarg,&tmp,10);
                    if ( *tmp || args->verbose<0 ) error("Could not parse argument: --verbosity %s\n", optarg);
                    if ( args->verbose > 3 ) hts_verbose = args->verbose;
                }
                break;
            case 'm': args->batch_fname = optarg; break;
            case 'M': merge_add_batch(args,optarg); break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( args->output_type==-1 )
    {
        if ( strlen(args->output_fname) > 3 && !strcasecmp(".gz",args->output_fname+strlen(args->output_fname)-3) ) args->output_type = FT_GZ;
        else args->output_type = FT_TAB_TEXT;
    }
    int i;
    for (i=optind; i<argc; i++) merge_add_batch(args,argv[i]);
    int do_merge = args->batch_fname || args->nbatch_fnames ? 1 : 0;
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
