/*  vcfroh.c -- HMM model for detecting runs of autozygosity.

    Copyright (C) 2013-2015 Genome Research Ltd.

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
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "HMM.h"
#include "smpl_ilist.h"

#define STATE_HW 0        // normal state, follows Hardy-Weinberg allele frequencies
#define STATE_AZ 1        // autozygous state

/** Genetic map */
typedef struct
{
    int pos;
    double rate;
}
genmap_t;

/** HMM data for each sample */
typedef struct
{
    double *eprob;      // emission probs [2*nsites,msites]
    uint32_t *sites;    // positions [nsites,msites]
    int nsites, msites;
    int beg;            // overlap sites to ignore
    int igenmap;        // current position in genmap
    int nused;          // some stats to detect if things didn't go wrong
    int nrid, *rid, *rid_off;   // for viterbi training, keep all chromosomes
    void *snapshot;             // hmm snapshot
}
smpl_t;

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    double t2AZ, t2HW;      // P(AZ|HW) and P(HW|AZ) parameters
    double unseen_PL, dflt_AF;

    char *genmap_fname;
    genmap_t *genmap;
    int ngenmap, mgenmap, igenmap;
    double rec_rate;        // constant recombination rate if > 0

    hmm_t *hmm;
    double baum_welch_th;
    int nrids, *rids, *rid_offs;    // multiple chroms with vi_training
    int nbuf_max, nbuf_olap;

    float *AFs;
    int32_t *itmp;
    int mAFs, nitmp, mitmp, pl_hdr_id, gt_hdr_id;

    double pl2p[256], *pdg;
    int32_t skip_rid, prev_rid, prev_pos;

    int ntot;                   // some stats to detect if things didn't go wrong
    smpl_t *smpl;               // HMM data for each sample
    smpl_ilist_t *af_smpl;      // list of samples to estimate AF from (--estimate-AF)
    smpl_ilist_t *roh_smpl;     // list of samples to analyze (--samples, --samples-file)
    char *estimate_AF;          // list of samples for AF estimate and query sample
    char **argv, *targets_list, *regions_list, *af_fname, *af_tag, *samples, *buffer_size;
    int argc, fake_PLs, snps_only, vi_training, samples_is_file;
}
args_t;

void set_tprob_genmap(hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob);
void set_tprob_rrate(hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob);

void *smalloc(size_t size)
{
    void *mem = malloc(size);
    if ( !mem ) error("malloc: Could not allocate %d bytes\n", (int)size);
    return mem;
}

static void init_data(args_t *args)
{
    int i;

    args->prev_rid = args->skip_rid = -1;
    args->hdr = args->files->readers[0].header;

    if ( !bcf_hdr_nsamples(args->hdr) ) error("No samples in the VCF?\n");

    if ( !args->fake_PLs )
    {
        args->pl_hdr_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, "PL");
        if ( !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,args->pl_hdr_id) )
            error("Error: The FORMAT/PL tag not found in the header, consider running with -G\n");
        if ( bcf_hdr_id2type(args->hdr,BCF_HL_FMT,args->pl_hdr_id)!=BCF_HT_INT ) 
            error("Error: The FORMAT/PL tag not defined as Integer in the header\n");
    }

    if ( args->estimate_AF && strcmp("-",args->estimate_AF) )
        args->af_smpl = smpl_ilist_init(args->hdr, args->estimate_AF, 1, SMPL_NONE);

    if ( args->estimate_AF || args->fake_PLs )
    {
        args->gt_hdr_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, "GT");
        if ( !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,args->gt_hdr_id) )
            error("Error: The FORMAT/GT tag not found in the header\n");
    }

    args->roh_smpl = smpl_ilist_init(args->hdr, args->samples, args->samples_is_file, SMPL_NONE);
    if ( args->samples )
    {
        // we may be able to subset to a few samples, for a text VCF this can be a major speedup
        if ( (bcf_sr_get_reader(args->files,0))->file->format.format==vcf )
        {
            kstring_t str = {0,0,0};
            smpl_ilist_t *tmp = args->roh_smpl, *rmme = NULL;
            if ( args->af_smpl )
            {
                for (i=0; i<args->roh_smpl->n; i++)
                {
                    if ( str.l ) kputc(',', &str);
                    kputs(args->hdr->samples[args->roh_smpl->idx[i]], &str);
                }
                for (i=0; i<args->af_smpl->n; i++)
                {
                    kputc(',', &str);
                    kputs(args->hdr->samples[args->af_smpl->idx[i]], &str);
                }
                rmme = tmp = smpl_ilist_init(args->hdr, str.s, 0, SMPL_NONE);
            }
            if ( tmp->n < bcf_hdr_nsamples(args->hdr) )
            {
                str.l = 0;
                for (i=0; i<tmp->n; i++)
                {
                    if ( str.l ) kputc(',', &str);
                    kputs(args->hdr->samples[tmp->idx[i]], &str);
                }
                int ret = bcf_hdr_set_samples(args->hdr, str.s, 0);
                if ( ret<0 ) error("Error parsing the list of samples: %s\n", str.s);
                else if ( ret>0 ) error("The %d-th sample not found in the VCF: %s\n", ret,str.s);

                // update sample ids
                smpl_ilist_destroy(args->roh_smpl);
                args->roh_smpl = smpl_ilist_init(args->hdr, args->samples, args->samples_is_file, SMPL_NONE);

                if ( args->af_smpl )
                {
                    smpl_ilist_destroy(args->af_smpl);
                    args->af_smpl = smpl_ilist_init(args->hdr, args->estimate_AF, 1, SMPL_NONE);
                }
            }
            free(str.s);
            if ( rmme )
                smpl_ilist_destroy(rmme);
        }
    }

    // check whether all samples are in this list. If so, the lookup will not be needed
    if ( args->af_smpl && args->af_smpl->n == bcf_hdr_nsamples(args->hdr) )
    {
        // all samples are in this list
        smpl_ilist_destroy(args->af_smpl);
        args->af_smpl = NULL;
    }

    if ( args->buffer_size )
    {
        args->nbuf_olap = -1;
        char *end;
        double tmp = strtod(args->buffer_size,&end);
        if ( *end )
        {
            if ( *end!=',') error("Could not parse: --buffer-size %s\n", args->buffer_size);
            args->nbuf_olap = strtol(end+1,&end,10);
            if ( *end || args->nbuf_olap<0 ) error("Could not parse: --bufer-size %s\n", args->buffer_size);
        }
        if ( tmp<0 )
            args->nbuf_max = fabs(tmp)*1e6/(4+8*2)/args->roh_smpl->n;
        else
            args->nbuf_max = tmp;

        if ( args->nbuf_olap<0 )
            args->nbuf_olap = args->nbuf_max*0.01;
    }
    fprintf(stderr,"Number of target samples: %d\n", args->roh_smpl->n);
    fprintf(stderr,"Number of --estimate-AF samples: %d\n", args->af_smpl ? args->af_smpl->n : 0);
    fprintf(stderr,"Number of sites in the buffer/overlap: ");
    if ( args->nbuf_max ) fprintf(stderr,"%d/%d\n", args->nbuf_max,args->nbuf_olap);
    else fprintf(stderr,"unlimited\n");

    args->smpl = (smpl_t*) calloc(args->roh_smpl->n,sizeof(smpl_t));

    for (i=0; i<256; i++) args->pl2p[i] = pow(10., -i/10.);

    // Init transition matrix and HMM
    double tprob[4];
    MAT(tprob,2,STATE_HW,STATE_HW) = 1 - args->t2AZ;
    MAT(tprob,2,STATE_HW,STATE_AZ) = args->t2HW;
    MAT(tprob,2,STATE_AZ,STATE_HW) = args->t2AZ;
    MAT(tprob,2,STATE_AZ,STATE_AZ) = 1 - args->t2HW; 

    if ( args->genmap_fname ) 
    {
        args->hmm = hmm_init(2, tprob, 0);
        hmm_set_tprob_func(args->hmm, set_tprob_genmap, args);
    }
    else if ( args->rec_rate > 0 )
    {
        args->hmm = hmm_init(2, tprob, 0);
        hmm_set_tprob_func(args->hmm, set_tprob_rrate, args);

    }
    else
        args->hmm = hmm_init(2, tprob, 10000);

    // print header
    printf("# This file was produced by: bcftools roh(%s+htslib-%s)\n", bcftools_version(),hts_version());
    printf("# The command line was:\tbcftools %s", args->argv[0]);
    for (i=1; i<args->argc; i++)
        printf(" %s",args->argv[i]);
    printf("\n#\n");
    i = 2;
    printf("# ROH");
    printf("\t[%d]Sample", i++);
    printf("\t[%d]Chromosome", i++);
    printf("\t[%d]Position", i++);
    printf("\t[%d]State (0:HW, 1:AZ)", i++);
    printf("\t[%d]Quality (fwd-bwd phred score)", i++);
    printf("\n");

    if ( args->vi_training)
    {
        i = 2;
        printf("# VT, Viterbi Training");
        printf("\t[%d]Sample", i++);
        printf("\t[%d]Iteration", i++);
        printf("\t[%d]dAZ", i++);
        printf("\t[%d]dHW", i++);
        printf("\t[%d]1 - P(HW|HW)", i++);
        printf("\t[%d]P(AZ|HW)", i++);
        printf("\t[%d]1 - P(AZ|AZ)", i++);
        printf("\t[%d]P(HW|AZ)", i++);
        printf("\n");
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->roh_smpl->n; i++)
    {
        free(args->smpl[i].eprob);
        free(args->smpl[i].sites);
        free(args->smpl[i].rid);
        free(args->smpl[i].rid_off);
        free(args->smpl[i].snapshot);
    }
    free(args->smpl);
    if ( args->af_smpl ) smpl_ilist_destroy(args->af_smpl);
    smpl_ilist_destroy(args->roh_smpl);
    free(args->rids);
    free(args->rid_offs);
    hmm_destroy(args->hmm);
    bcf_sr_destroy(args->files);
    free(args->AFs); free(args->pdg);
    free(args->genmap);
    free(args->itmp);
    free(args->samples);
}

static int load_genmap(args_t *args, bcf1_t *line)
{
    if ( !args->genmap_fname ) { args->ngenmap = 0; return 0; }

    kstring_t str = {0,0,0};
    char *fname = strstr(args->genmap_fname,"{CHROM}");
    if ( fname )
    {
        kputsn(args->genmap_fname, fname - args->genmap_fname, &str);
        kputs(bcf_seqname(args->hdr,line), &str);
        kputs(fname+7,&str);
        fname = str.s;
    }
    else
        fname = args->genmap_fname;

    htsFile *fp = hts_open(fname, "rb");
    if ( !fp )
    {
        args->ngenmap = 0;
        return -1;
    }

    hts_getline(fp, KS_SEP_LINE, &str);
    if ( strcmp(str.s,"position COMBINED_rate(cM/Mb) Genetic_Map(cM)") )
        error("Unexpected header, found:\n\t[%s], but expected:\n\t[position COMBINED_rate(cM/Mb) Genetic_Map(cM)]\n", fname, str.s);

    args->ngenmap = args->igenmap = 0;
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        args->ngenmap++;
        hts_expand(genmap_t,args->ngenmap,args->mgenmap,args->genmap);
        genmap_t *gm = &args->genmap[args->ngenmap-1];

        char *tmp, *end;
        gm->pos = strtol(str.s, &tmp, 10);
        if ( str.s==tmp ) error("Could not parse %s: %s\n", fname, str.s);

        // skip second column
        tmp++;
        while ( *tmp && !isspace(*tmp) ) tmp++;

        // read the genetic map in cM
        gm->rate = strtod(tmp+1, &end);
        if ( tmp+1==end ) error("Could not parse %s: %s\n", fname, str.s);
    }
    if ( !args->ngenmap ) error("Genetic map empty?\n");
    int i;
    for (i=0; i<args->ngenmap; i++) args->genmap[i].rate /= args->genmap[args->ngenmap-1].rate; // scale to 1
    if ( hts_close(fp) ) error("Close failed\n");
    free(str.s);
    return 0;
}

static double get_genmap_rate(args_t *args, int start, int end)
{
    // position i to be equal to or smaller than start
    int i = args->igenmap;
    if ( args->genmap[i].pos > start )
    {
        while ( i>0 && args->genmap[i].pos > start ) i--;
    }
    else
    {
        while ( i+1<args->ngenmap && args->genmap[i+1].pos < start ) i++;
    }
    // position j to be equal or larger than end
    int j = i;
    while ( j+1<args->ngenmap && args->genmap[j].pos < end ) j++;
    if ( i==j )
    {
        args->igenmap = i;
        return 0;
    }

    if ( start <  args->genmap[i].pos ) start = args->genmap[i].pos;
    if ( end >  args->genmap[j].pos ) end = args->genmap[j].pos;
    double rate = (args->genmap[j].rate - args->genmap[i].rate)/(args->genmap[j].pos - args->genmap[i].pos) * (end-start);
    args->igenmap = j;
    return rate;
}

void set_tprob_genmap(hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob)
{
    args_t *args = (args_t*) data;
    double ci = get_genmap_rate(args, prev_pos, pos);
    MAT(tprob,2,STATE_HW,STATE_AZ) *= ci;
    MAT(tprob,2,STATE_AZ,STATE_HW) *= ci;
    MAT(tprob,2,STATE_AZ,STATE_AZ)  = 1 - MAT(tprob,2,STATE_HW,STATE_AZ);
    MAT(tprob,2,STATE_HW,STATE_HW)  = 1 - MAT(tprob,2,STATE_AZ,STATE_HW);
}

void set_tprob_rrate(hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob)
{
    args_t *args = (args_t*) data;
    double ci = (pos - prev_pos) * args->rec_rate;
    MAT(tprob,2,STATE_HW,STATE_AZ) *= ci;
    MAT(tprob,2,STATE_AZ,STATE_HW) *= ci;
    MAT(tprob,2,STATE_AZ,STATE_AZ)  = 1 - MAT(tprob,2,STATE_HW,STATE_AZ);
    MAT(tprob,2,STATE_HW,STATE_HW)  = 1 - MAT(tprob,2,STATE_AZ,STATE_HW);
}


/**
 *  This function implements the HMM model:
 *    D = Data, AZ = autozygosity, HW = Hardy-Weinberg (non-autozygosity),
 *    f = non-ref allele frequency
 *
 *  Emission probabilities:
 *    oAZ = P_i(D|AZ) = (1-f)*P(D|RR) + f*P(D|AA)
 *    oHW = P_i(D|HW) = (1-f)^2 * P(D|RR) + f^2 * P(D|AA) + 2*f*(1-f)*P(D|RA)
 *
 *  Transition probabilities:
 *    tAZ = P(AZ|HW)  .. parameter
 *    tHW = P(HW|AZ)  .. parameter
 *
 *    ci  = P_i(C)    .. probability of cross-over at site i, from genetic map
 *
 *    AZi = P_i(AZ)   .. probability of site i being AZ/non-AZ, scaled so that AZi+HWi = 1
 *    HWi = P_i(HW)
 *
 *    P_i(AZ|HW) = P(AZ|HW) * ci * HW{i-1}     = tAZ * ci * (1 - AZ{i-1})
 *    P_i(HW|AZ) = P(HW|AZ) * ci * AZ{i-1}     = tHW * ci * AZ{i-1}
 *    P_i(AZ|AZ) = 1 - P_i(HW|AZ)
 *    P_i(HW|HW) = 1 - P_i(AZ|HW)
 *
 */

static void flush_viterbi(args_t *args, int ismpl)
{
    smpl_t *smpl = &args->smpl[ismpl];
    if ( !smpl->nsites ) return;

    const char *name = args->hdr->samples[ args->roh_smpl->idx[ismpl] ];

    int i,j,k;

    if ( !args->vi_training ) // single viterbi pass
    {
        hmm_restore(args->hmm, smpl->snapshot); 
        int end = (args->nbuf_max && smpl->nsites >= args->nbuf_max && smpl->nsites > args->nbuf_olap) ? smpl->nsites - 0.5*args->nbuf_olap : smpl->nsites;
        if ( end < smpl->nsites )
            smpl->snapshot = hmm_snapshot(args->hmm, smpl->snapshot, smpl->nsites - args->nbuf_olap - 1);

        args->igenmap = smpl->igenmap;
        hmm_run_viterbi(args->hmm, smpl->nsites, smpl->eprob, smpl->sites);
        hmm_run_fwd_bwd(args->hmm, smpl->nsites, smpl->eprob, smpl->sites);
        double *fwd = hmm_get_fwd_bwd_prob(args->hmm);

        const char *chr  = bcf_hdr_id2name(args->hdr,args->prev_rid);
        uint8_t *vpath   = hmm_get_viterbi_path(args->hmm);

        for (i=smpl->beg; i<end; i++)
        {
            int state = vpath[i*2]==STATE_AZ ? 1 : 0;
            double *pval = fwd + i*2;
            printf("ROH\t%s\t%s\t%d\t%d\t%.1f\n", name,chr,smpl->sites[i]+1, state, phred_score(1.0-pval[state]));
        }

        if ( end < smpl->nsites )
        {
            end = smpl->nsites - args->nbuf_olap;
            memmove(smpl->sites, smpl->sites + end, sizeof(*smpl->sites)*args->nbuf_olap);
            memmove(smpl->eprob, smpl->eprob + end*2, sizeof(*smpl->eprob)*args->nbuf_olap*2);
            smpl->nsites  = args->nbuf_olap;
            smpl->beg     = 0.5*args->nbuf_olap;
            smpl->igenmap = args->igenmap;
        }
        else
        {
            smpl->nsites  = 0;
            smpl->beg     = 0;
            smpl->igenmap = 0;
        }

        return;
    }


    // viterbi training, multiple chromosomes
    double t2az_prev, t2hw_prev;
    double deltaz, delthw;

    double *tprob_arr = hmm_get_tprob(args->hmm);
    MAT(tprob_arr,2,STATE_HW,STATE_HW) = 1 - args->t2AZ;
    MAT(tprob_arr,2,STATE_HW,STATE_AZ) = args->t2HW;
    MAT(tprob_arr,2,STATE_AZ,STATE_HW) = args->t2AZ;
    MAT(tprob_arr,2,STATE_AZ,STATE_AZ) = 1 - args->t2HW; 
    if ( args->genmap_fname || args->rec_rate > 0 )
        hmm_set_tprob(args->hmm, tprob_arr, 0);
    else
        hmm_set_tprob(args->hmm, tprob_arr, 10000);

    int niter = 0;
    do
    {
        tprob_arr = hmm_get_tprob(args->hmm);
        t2az_prev = MAT(tprob_arr,2,STATE_AZ,STATE_HW); //args->t2AZ;
        t2hw_prev = MAT(tprob_arr,2,STATE_HW,STATE_AZ); //args->t2HW;
        double tprob_new[] = { 0,0,0,0 };
        for (i=0; i<smpl->nrid; i++)
        {
            int ioff = smpl->rid_off[i];
            int nsites = (i+1==smpl->nrid ? smpl->nsites : smpl->rid_off[i+1]) - ioff;
            args->igenmap = 0;
            tprob_arr = hmm_run_baum_welch(args->hmm, nsites, smpl->eprob+ioff*2, smpl->sites+ioff);
            for (j=0; j<2; j++)
                for (k=0; k<2; k++) MAT(tprob_new,2,j,k) += MAT(tprob_arr,2,j,k);
        }
        for (j=0; j<2; j++)
            for (k=0; k<2; k++) MAT(tprob_new,2,j,k) /= smpl->nrid;

        if ( args->genmap_fname || args->rec_rate > 0 )
            hmm_set_tprob(args->hmm, tprob_new, 0);
        else
            hmm_set_tprob(args->hmm, tprob_new, 10000);

        deltaz = fabs(MAT(tprob_new,2,STATE_AZ,STATE_HW)-t2az_prev);
        delthw = fabs(MAT(tprob_new,2,STATE_HW,STATE_AZ)-t2hw_prev);
        niter++;
        printf("VT\t%s\t%d\t%e\t%e\t%e\t%e\t%e\t%e\n", 
            name,niter,deltaz,delthw,
            1-MAT(tprob_new,2,STATE_HW,STATE_HW),MAT(tprob_new,2,STATE_AZ,STATE_HW),
            1-MAT(tprob_new,2,STATE_AZ,STATE_AZ),MAT(tprob_new,2,STATE_HW,STATE_AZ));
    }
    while ( deltaz > args->baum_welch_th || delthw > args->baum_welch_th );
    
    // output the results
    for (i=0; i<smpl->nrid; i++)
    {
        int ioff = smpl->rid_off[i];
        int nsites = (i+1==smpl->nrid ? smpl->nsites : smpl->rid_off[i+1]) - ioff;
        args->igenmap = 0;
        hmm_run_viterbi(args->hmm, nsites, smpl->eprob+ioff*2, smpl->sites+ioff);
        hmm_run_fwd_bwd(args->hmm, nsites, smpl->eprob+ioff*2, smpl->sites+ioff);
        uint8_t *vpath = hmm_get_viterbi_path(args->hmm);
        double  *fwd   = hmm_get_fwd_bwd_prob(args->hmm);

        const char *chr = bcf_hdr_id2name(args->hdr,smpl->rid[i]);
        for (j=0; j<nsites; j++)
        {
            int state = vpath[j*2]==STATE_AZ ? 1 : 0;
            double *pval = fwd + j*2;
            printf("ROH\t%s\t%s\t%d\t%d\t%.1f\n", name,chr,smpl->sites[ioff+j]+1, state, phred_score(1.0-pval[state]));
        }
    }
}


int read_AF(bcf_sr_regions_t *tgt, bcf1_t *line, double *alt_freq)
{
    if ( tgt->nals != line->n_allele ) return -1;    // number of alleles does not match

    int i;
    for (i=0; i<tgt->nals; i++)
        if ( strcmp(line->d.allele[i],tgt->als[i]) ) break; // we could be smarter, see vcmp
    if ( i<tgt->nals ) return -1;

    char *tmp, *str = tgt->line.s;
    i = 0;
    while ( *str && i<3 ) 
    {
        if ( *str=='\t' ) i++;
        str++;
    }
    *alt_freq = strtod(str, &tmp);
    if ( *tmp && !isspace(*tmp) )
    {
        if ( str[0]=='.' && (!str[1] || isspace(str[1])) ) return -1; // missing value
        error("Could not parse: [%s]\n", tgt->line.s);
    }
    if ( *alt_freq<0 || *alt_freq>1 ) error("Could not parse AF: [%s]\n", tgt->line.s);
    return 0;
}

int8_t *get_GT(args_t *args, bcf1_t *line)
{
    int i;
    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==args->gt_hdr_id ) break;
    if ( i==line->n_fmt ) return NULL;        // the tag is not present in this record

    bcf_fmt_t *fmt = &line->d.fmt[i];
    if ( fmt->n!=2 ) return NULL;             // not diploid

    if ( fmt->type!=BCF_BT_INT8 ) error("This is unexpected, GT type is %d\n", fmt->type);
    return (int8_t*) fmt->p;
}

int estimate_AF_from_GT(args_t *args, int8_t *gt, double *alt_freq)
{
    int i, nalt = 0, nref = 0;
    if ( args->af_smpl )        // subset samples for AF estimate
    {
        for (i=0; i<args->af_smpl->n; i++)
        {
            int ismpl = args->af_smpl->idx[i];
            if ( bcf_gt_is_missing(gt[2*ismpl]) || bcf_gt_is_missing(gt[2*ismpl+1]) ) continue;

            if ( bcf_gt_allele(gt[2*ismpl]) ) nalt++;
            else nref++;

            if ( bcf_gt_allele(gt[2*ismpl+1]) ) nalt++;
            else nref++;
        }
    }
    else                        // all samples used in AF estimate
    {
        int8_t *end = gt + 2*bcf_hdr_nsamples(args->hdr);
        while ( gt < end )
        {
            if ( bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1]) ) continue;

            if ( bcf_gt_allele(gt[0]) ) nalt++;
            else nref++;

            if ( bcf_gt_allele(gt[1]) ) nalt++;
            else nref++;

            gt += 2;
        }
    }
    if ( !nalt && !nref ) return -1;

    *alt_freq = (double)nalt / (nalt + nref);
    return 0;
}


int process_line(args_t *args, bcf1_t *line)
{
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    double alt_freq;
    int8_t *GTs = NULL;

    // Set allele frequency
    int ret = 0;
    if ( args->af_tag )
    {
        // Use an INFO tag provided by the user
        ret = bcf_get_info_float(args->hdr, line, args->af_tag, &args->AFs, &args->mAFs);
        if ( ret==1 )
            alt_freq = args->AFs[0];
        if ( ret==-2 )
            error("Type mismatch for INFO/%s tag at %s:%d\n", args->af_tag, bcf_seqname(args->hdr,line), line->pos+1);
    }
    else if ( args->af_fname ) 
    {
        // Read AF from a file
        ret = read_AF(args->files->targets, line, &alt_freq);
    }
    else if ( args->dflt_AF > 0 )
    {
        alt_freq = args->dflt_AF;
    }
    else if ( args->estimate_AF )
    {
        // Estimate AF from GTs of all samples or samples listed in a file
        GTs = get_GT(args, line);
        if ( !GTs ) return -1;
        ret = estimate_AF_from_GT(args, GTs, &alt_freq);
    }
    else
    {
        // Use AC/AN
        int AC = -1, AN = 0;
        ret = bcf_get_info_int32(args->hdr, line, "AN", &args->itmp, &args->mitmp);
        if ( ret==1 )
        {
            AN = args->itmp[0];
            ret = bcf_get_info_int32(args->hdr, line, "AC", &args->itmp, &args->mitmp);
            if ( ret>0 )
                AC = args->itmp[0];
        }
        if ( AN<=0 || AC<0 ) 
            ret = -1;
        else 
            alt_freq = (double) AC/AN;
    }

    if ( ret<0 ) return ret;
    if ( alt_freq==0.0 ) return -1;

    int i,j;
    bcf_fmt_t *fmt_pl = NULL;
    if ( args->fake_PLs )
    {
        if ( !GTs ) GTs = get_GT(args, line);
    }
    else
    {
        for (i=0; i<line->n_fmt; i++)
            if ( line->d.fmt[i].id==args->pl_hdr_id ) break;
        if ( i==line->n_fmt ) return -1;                    // the PL tag is not present in this record

        fmt_pl = &line->d.fmt[i];
        if ( fmt_pl->n != 3 ) return -1;                    // not biallelic & diploid
    }

    for (i=0; i<args->roh_smpl->n; i++)
    {
        int ismpl = args->roh_smpl->idx[i];

        // set P(D|G)
        double pdg[3];
        if ( args->fake_PLs )
        {
            int8_t *gt = GTs + 2*ismpl;
            if ( bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1]) ) continue;

            int a = bcf_gt_allele(gt[0]);
            int b = bcf_gt_allele(gt[1]);
            if ( a!=b )
            {
                pdg[0] = pdg[2] = args->unseen_PL;
                pdg[1] = 1 - 2*args->unseen_PL;
            }
            else if ( a==0 )
            {
                pdg[0] = 1 - 2*args->unseen_PL;
                pdg[1] = pdg[2] = args->unseen_PL;
            }
            else
            {
                pdg[0] = pdg[1] = args->unseen_PL;
                pdg[2] = 1 - 2*args->unseen_PL;
            }
        }
        else
        {
            #define BRANCH(type_t) \
            { \
                type_t *p = (type_t*)fmt_pl->p + 3*ismpl; \
                if ( p[0]<0 || p[1]<0 || p[2]<0 ) continue;    /* missing value */ \
                if ( p[0]==p[1] && p[0]==p[2] ) continue;      /* all values are the same */ \
                for (j=0; j<3; j++) \
                    pdg[j] = p[j] < 256 ? args->pl2p[ p[j] ] : args->pl2p[255]; \
            }
            switch (fmt_pl->type) {
                case BCF_BT_INT8:  BRANCH(int8_t); break;
                case BCF_BT_INT16: BRANCH(int16_t); break;
                case BCF_BT_INT32: BRANCH(int32_t); break;
                default: fprintf(stderr,"Unknown format type for PL: %s:%d .. fmt->type=%d\n", __FILE__,__LINE__, fmt_pl->type); exit(1);
            }
            #undef BRANCH
        }

        double sum = pdg[0] + pdg[1] + pdg[2];
        if ( !sum ) continue;
        for (j=0; j<3; j++) pdg[j] /= sum;

        smpl_t *smpl = &args->smpl[i];
        smpl->nused++;

        if ( smpl->nsites >= smpl->msites )
        {
            hts_expand(uint32_t,smpl->nsites+1,smpl->msites,smpl->sites);
            smpl->eprob = (double*) realloc(smpl->eprob,sizeof(*smpl->eprob)*smpl->msites*2);
        }
        
        // Calculate emission probabilities P(D|AZ) and P(D|HW)
        double *eprob = &smpl->eprob[2*smpl->nsites];
        eprob[STATE_AZ] = pdg[0]*(1-alt_freq) + pdg[2]*alt_freq;
        eprob[STATE_HW] = pdg[0]*(1-alt_freq)*(1-alt_freq) + 2*pdg[1]*(1-alt_freq)*alt_freq + pdg[2]*alt_freq*alt_freq;
        
        smpl->sites[smpl->nsites] = line->pos;
        smpl->nsites++;

        if ( args->vi_training )
        {
            if ( !smpl->nrid || line->rid!=smpl->rid[smpl->nrid-1] )
            {
                smpl->nrid++;
                smpl->rid = (int*) realloc(smpl->rid,sizeof(*smpl->rid)*smpl->nrid);
                smpl->rid[smpl->nrid-1] = line->rid;
                smpl->rid_off = (int*) realloc(smpl->rid_off,sizeof(*smpl->rid_off)*smpl->nrid);
                smpl->rid_off[smpl->nrid-1] = smpl->nsites - 1;
            }
        }
        else if ( args->nbuf_max && smpl->nsites >= args->nbuf_max ) flush_viterbi(args, i);
    }

    return 0;
}

static void vcfroh(args_t *args, bcf1_t *line)
{
    int i;

    // Are we done?
    if ( !line )
    { 
        for (i=0; i<args->roh_smpl->n; i++) flush_viterbi(args, i);
        return; 
    }
    args->ntot++;

    // Skip unwanted lines
    if ( line->rid == args->skip_rid ) return;
    if ( line->n_allele==1 ) return;    // no ALT allele
    if ( line->n_allele!=2 ) return;    // only biallelic sites
    if ( args->snps_only && !bcf_is_snp(line) ) return;

    // Initialize genetic map
    int skip_rid = 0;
    if ( args->prev_rid<0 )
    {
        args->prev_rid = line->rid;
        args->prev_pos = line->pos;
        skip_rid = load_genmap(args, line);
    }

    // New chromosome?
    if ( args->prev_rid!=line->rid )
    {
        skip_rid = load_genmap(args, line);
        if ( !args->vi_training )
        {
            for (i=0; i<args->roh_smpl->n; i++) flush_viterbi(args, i);
        }
        args->prev_rid = line->rid;
        args->prev_pos = line->pos;
    }

    if ( skip_rid )
    {
        fprintf(stderr,"Skipping the sequence, no genmap for %s\n", bcf_seqname(args->hdr,line));
        args->skip_rid = line->rid;
        return;
    }
    if ( args->prev_pos > line->pos ) error("The file is not sorted?!\n");

    args->prev_rid = line->rid;
    args->prev_pos = line->pos;


    // parse the new line
    process_line(args, line);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   HMM model for detecting runs of autozygosity.\n");
    fprintf(stderr, "Usage:   bcftools roh [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "        --AF-dflt <float>              if AF is not known, use this allele frequency [skip]\n");
    fprintf(stderr, "        --AF-tag <TAG>                 use TAG for allele frequency\n");
    fprintf(stderr, "        --AF-file <file>               read allele frequencies from file (CHR\\tPOS\\tREF,ALT\\tAF)\n");
    fprintf(stderr, "    -b  --buffer-size <int[,int]>      buffer size and the number of overlapping sites, 0 for unlimited [0]\n");
    fprintf(stderr, "                                           If the first number is negative, it is interpreted as the maximum memory to\n");
    fprintf(stderr, "                                           use, in MB. The default overlap is set to roughly 1%% of the buffer size.\n");
    fprintf(stderr, "    -e, --estimate-AF <file>           estimate AF from GTs of all samples (\"-\") or samples listed in <file>\n");
    fprintf(stderr, "    -G, --GTs-only <float>             use GTs and ignore PLs, instead using <float> for PL of the two least likely genotypes.\n");
    fprintf(stderr, "                                           Safe value to use is 30 to account for GT errors.\n");
    fprintf(stderr, "    -I, --skip-indels                  skip indels as their genotypes are enriched for errors\n");
    fprintf(stderr, "    -m, --genetic-map <file>           genetic map in IMPUTE2 format, single file or mask, where string \"{CHROM}\"\n");
    fprintf(stderr, "                                           is replaced with chromosome name\n");
    fprintf(stderr, "    -M, --rec-rate <float>             constant recombination rate per bp\n");
    fprintf(stderr, "    -r, --regions <region>             restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>          restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --samples <list>               list of samples to analyze [all samples]\n");
    fprintf(stderr, "    -S, --samples-file <file>          file of samples to analyze [all samples]\n");
    fprintf(stderr, "    -t, --targets <region>             similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>          similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "HMM Options:\n");
    fprintf(stderr, "    -a, --hw-to-az <float>             P(AZ|HW) transition probability from HW (Hardy-Weinberg) to AZ (autozygous) state [6.7e-8]\n");
    fprintf(stderr, "    -H, --az-to-hw <float>             P(HW|AZ) transition probability from AZ to HW state [5e-9]\n");
    fprintf(stderr, "    -V, --viterbi-training <float>     estimate HMM parameters, <float> is the convergence threshold, e.g. 1e-10 (experimental)\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfroh(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->t2AZ    = 6.7e-8;
    args->t2HW    = 5e-9;
    args->rec_rate = 0;
    int regions_is_file = 0, targets_is_file = 0;

    static struct option loptions[] =
    {
        {"AF-tag",1,0,0},
        {"AF-file",1,0,1},
        {"AF-dflt",1,0,2},
        {"buffer-size",1,0,'b'},
        {"estimate-AF",1,0,'e'},
        {"GTs-only",1,0,'G'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"hw-to-az",1,0,'a'},
        {"az-to-hw",1,0,'H'},
        {"viterbi-training",1,0,'V'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"genetic-map",1,0,'m'},
        {"rec-rate",1,0,'M'},
        {"skip-indels",0,0,'I'},
        {0,0,0,0}
    };

    int naf_opts = 0;
    char *tmp;
    while ((c = getopt_long(argc, argv, "h?r:R:t:T:H:a:s:S:m:M:G:Ia:e:V:b:",loptions,NULL)) >= 0) {
        switch (c) {
            case 0: args->af_tag = optarg; naf_opts++; break;
            case 1: args->af_fname = optarg; naf_opts++; break;
            case 2: 
                args->dflt_AF = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: --AF-dflt %s\n", optarg);
                break;
            case 'e': args->estimate_AF = optarg; naf_opts++; break;
            case 'b': args->buffer_size = optarg; break;
            case 'I': args->snps_only = 1; break;
            case 'G':
                args->fake_PLs = 1; 
                args->unseen_PL = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -G %s\n", optarg);
                args->unseen_PL = pow(10,-args->unseen_PL/10.); 
                break;
            case 'm': args->genmap_fname = optarg; break;
            case 'M':
                args->rec_rate = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -M %s\n", optarg);
                break;
            case 's': args->samples = strdup(optarg); break;
            case 'S': args->samples = strdup(optarg); args->samples_is_file = 1; break;
            case 'a':
                args->t2AZ = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -a %s\n", optarg);
                break;
            case 'H':
                args->t2HW = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -H %s\n", optarg);
                break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 'V': 
                args->vi_training = 1; 
                args->baum_welch_th = strtod(optarg,&tmp); 
                if ( *tmp ) error("Could not parse: --viterbi-training %s\n", optarg);
                break;
            case 'h': 
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( argc<optind+1 ) usage(args);
    if ( args->vi_training && args->buffer_size ) error("Error: cannot use -b with -V\n");
    if ( args->t2AZ<0 || args->t2AZ>1 ) error("Error: The parameter --hw-to-az is not in [0,1]\n", args->t2AZ);
    if ( args->t2HW<0 || args->t2HW>1 ) error("Error: The parameter --az-to-hw is not in [0,1]\n", args->t2HW);
    if ( naf_opts>1 ) error("Error: The options --AF-tag, --AF-file and -e are mutually exclusive\n");
    if ( args->af_fname && args->targets_list ) error("Error: The options --AF-file and -t are mutually exclusive\n");
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( args->af_fname )
    {
        if ( bcf_sr_set_targets(args->files, args->af_fname, 1, 3)<0 )
            error("Failed to read the targets: %s\n", args->af_fname);
    }
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open %s: %s\n", argv[optind],bcf_sr_strerror(args->files->errnum));

    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        vcfroh(args, args->files->readers[0].buffer[0]);
    }
    vcfroh(args, NULL);
    int i, nmin = 0;
    for (i=0; i<args->roh_smpl->n; i++)
        if ( !i || args->smpl[i].nused < nmin ) nmin = args->smpl[i].nused;
    fprintf(stderr,"Number of lines total/processed: %d/%d\n", args->ntot,nmin);
    if ( nmin==0 )
    {
        fprintf(stderr,"No usable sites were found.");
        if ( !naf_opts && !args->dflt_AF ) fprintf(stderr, " Consider using one of the AF options.\n");
    }
    destroy_data(args);
    free(args);
    return 0;
}


