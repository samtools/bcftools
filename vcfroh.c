#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "rbuf.h"

/** Buffered sites, separate for each sample to allow missing genotypes */
typedef struct
{
    rbuf_t rbuf;
    double *ohw, *oaz;  // P(D|HW) and P(D|AZ)
    double last_az;
    uint32_t *pos;
}
smpl_t;

/** Genetic map */
typedef struct
{
    int pos;
    double rate;
}
genmap_t;

/** Viterbi path element for two-state HMM */
typedef struct
{
    double pAZ;
    char ptr;
}
path_t;

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    double tAZ, tHW;        // P(AZ|HW), P(HW|AZ)
    double *fwd, *bwd;      // HMM forward and backward autozygosity probs scaled to hw+az=1
    path_t *viterbi;

    genmap_t *genmap;
    int ngenmap, mgenmap, igenmap;

    int nsmpl, *ismpl, *als;
    int mwin;
    smpl_t *smpl;
    int32_t *PLs, *AN, *ACs;
    double *AFs;
    double pl2p[256], *pdg;
    int mPLs, mAFs, mAN, mACs, mpdg;
    int ntot, nused;
    int prev_rid, skip_rid;
    double unseen_PL;

    char **argv, *targets_fname, *regions_fname, *samples_fname;
    char *genmap_fname;
    int argc, counts_only, fwd_bwd, fake_PLs, biallelic_only, snps_only, estimate_AF;
}
args_t;

void *smalloc(size_t size)
{
    void *mem = malloc(size);
    if ( !mem ) error("malloc: Could not allocate %d bytes\n", (int)size);
    return mem;
}

static void init_data(args_t *args)
{
    args->prev_rid = args->skip_rid = -1;
    args->hdr = args->files->readers[0].header;
    if ( !args->samples_fname ) args->samples_fname = "-";
    if ( !bcf_sr_set_samples(args->files, args->samples_fname) )
        error("Error: could not set the samples %s\n", args->samples_fname);
    args->nsmpl = args->files->n_smpl;
    args->ismpl = args->files->readers[0].samples;

    args->als   = (int*) smalloc(sizeof(int)*args->nsmpl);
    args->smpl  = (smpl_t*) smalloc(sizeof(smpl_t)*args->nsmpl);
    if ( args->fwd_bwd )
    {
        args->fwd   = (double*) smalloc(sizeof(double)*(args->mwin+1));
        args->bwd   = (double*) smalloc(sizeof(double)*(args->mwin+1));
    }
    args->viterbi = (path_t*) smalloc(sizeof(path_t)*(args->mwin+1));
    args->pdg = (double*) smalloc(sizeof(double)*args->nsmpl*3);
    int i;
    for (i=0; i<args->nsmpl; i++)
    {
        smpl_t *smpl = &args->smpl[i];
        rbuf_init(&smpl->rbuf, args->mwin);
        smpl->ohw = (double*) smalloc(sizeof(double)*args->mwin);
        smpl->oaz = (double*) smalloc(sizeof(double)*args->mwin);
        smpl->pos = (uint32_t*) smalloc(sizeof(uint32_t)*args->mwin);
        smpl->last_az = -1;
    }
    for (i=0; i<256; i++)
        args->pl2p[i] = pow(10., -i/10.);

    // print header
    printf("# This file was produced by: bcftools roh(%s)\n", bcftools_version());
    printf("# The command line was:\tbcftools %s", args->argv[0]);
    for (i=1; i<args->argc; i++)
        printf(" %s",args->argv[i]);
    printf("\n#\n");
    if ( args->counts_only ) 
        printf("# [1]Sample\t[2]Chromosome\t[3]Position\t[4]HOM rate\t[5]HET rate\n");
    else if ( args->fwd_bwd )
        printf("# [1]Sample\t[2]Chromosome\t[3]Position\t[4]ROH p-value\n");
    else
        printf("# [1]Sample\t[2]Chromosome\t[3]Position\t[4]p-value\t[5]ROH\n");
}

static void destroy_data(args_t *args)
{
    bcf_sr_destroy(args->files);
    int i;
    for (i=0; i<args->nsmpl; i++)
    {
        smpl_t *smpl = &args->smpl[i];
        free(smpl->ohw); free(smpl->oaz); free(smpl->pos);
    }
    free(args->als);
    free(args->smpl);
    free(args->fwd); free(args->bwd); free(args->viterbi);
    free(args->PLs); free(args->AFs); free(args->pdg);
    free(args->AN); free(args->ACs);
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
        gm->rate = strtod(tmp+2, &end);
        if ( tmp+2==end ) error("Could not parse %s: %s\n", fname, str.s);
        gm->rate *= 0.01;
    }
    if ( hts_close(fp) ) error("Close failed\n");
    free(str.s);
    return 0;
}

static double get_genmap_rate(args_t *args, int pos)
{
    int i = args->igenmap;
    if ( args->genmap[i].pos > pos )
    {
        while ( i>0 && args->genmap[i].pos > pos ) i--;
    }
    else
    {
        while ( i+1<args->ngenmap && args->genmap[i+1].pos <= pos ) i++;
    }
    args->igenmap = i;
    return args->genmap[i].rate;
}


static void flush_counts(args_t *args, int ismpl, int n)
{
    smpl_t *smpl = &args->smpl[ismpl];
    
    int pos = smpl->pos[ smpl->rbuf.f ];
    int nhw = 0, naz = 0;
    int ir;

    for (ir=-1; rbuf_next(&smpl->rbuf,&ir); )
    {
        if ( smpl->ohw[ir] > smpl->oaz[ir] ) nhw++;
        else naz++;
        if ( smpl->pos[ir] - pos + 1 >= args->mwin  )
        {
            double az_rate = (double)naz*args->mwin/(smpl->pos[ir] - pos + 1);
            double hw_rate = (double)nhw*args->mwin/(smpl->pos[ir] - pos + 1);
            printf("%s\t%s\t%d\t%e\t%e\n", args->hdr->samples[args->ismpl[ismpl]], args->hdr->id[BCF_DT_CTG][args->prev_rid].key, pos+1, az_rate, hw_rate);
            pos = smpl->pos[ir];
            nhw = naz = 0;
        }
    }
    rbuf_shift_n(&smpl->rbuf, n);
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
 *    P(AZ|AZ) = 1 - P(HW|AZ) = 1 - tHW
 *    P(HW|HW) = 1 - P(AZ|HW) = 1 - tAZ
 *
 *    ci  = P_i(C)    .. probability of cross-over at site i, from genetic map
 *
 *    AZi = P_i(AZ)   .. probability of site i being AZ/non-AZ, scaled so that AZi+HWi = 1
 *    HWi = P_i(HW) 
 *
 *    P_i(AZ|AZ) = P(AZ|AZ) * (1-ci) * AZ{i-1} = (1-tHW) * (1-ci) * AZ{i-1}
 *    P_i(AZ|HW) = P(AZ|HW) * ci * HW{i-1}     = tAZ * ci * (1 - AZ{i-1})
 *
 *    P_i(HW|HW) = P(HW|HW) * (1-ci) * HW{i-1} = (1-tAZ) * (1-ci) * (1 - AZ{i-1})
 *    P_i(HW|AZ) = P(HW|AZ) * ci * AZ{i-1}     = tHW * ci * AZ{i-1}
 *
 *    ------------------------------------------------------
 *
 *    P_{i+1}(AZ) = P_{i+1}(D|AZ) * [ P_i(AZ|AZ) + P_i(AZ|HW) ] 
 *                = oAZ * [ (1-tHW) * (1-ci) * AZ{i-1} + tAZ * ci * (1-AZ{i-1})]
 *    P_{i+1}(HW) = P_{i+1}(D|HW) * [ P_i(HW|HW) + P_i(HW|AZ) ] 
 *                = oHW * [ (1-tAZ) * (1-ci) * (1-AZ{i-1}) + tHW * ci * AZ{i-1} ]
 */
static void flush_buffer_fwd_bwd(args_t *args, int ismpl, int n)
{
    smpl_t *smpl = &args->smpl[ismpl];
    args->fwd[0] = smpl->last_az < 0 ?  0.5 : smpl->last_az;

    int i, ir;
    double pAZ, pHW;
    for (i=1; i<=n; i++)
    {
        ir = rbuf_kth(&smpl->rbuf, i-1);
        double ci = args->ngenmap ? get_genmap_rate(args, smpl->pos[ir]) : 0.5;

        // P_{i+1}(AZ) = oAZ * [(1-tHW) * (1-ci) * AZ{i-1} + tAZ * ci * (1-AZ{i-1})]
        pAZ = smpl->oaz[ir] * ( (1-args->tHW) * (1-ci) * args->fwd[i-1] + args->tAZ * ci * (1-args->fwd[i-1]) );

        // P_{i+1}(HW) = oHW * [(1-tAZ) * (1-ci) * (1-AZ{i-1}) + tHW * ci * AZ{i-1}]
        pHW = smpl->ohw[ir] * ( (1-args->tAZ) * (1-ci) * (1-args->fwd[i-1]) + args->tHW * ci * args->fwd[i-1]);

        args->fwd[i] = pAZ / (pAZ+pHW);
        // printf("fwd\ti=%d ir=%d  %d\t oaz=%e ohw=%e \t pAZ=%e  pHW=%e  fwd=%e\n", i,ir,smpl->pos[ir]+1, smpl->oaz[ir],smpl->ohw[ir],pAZ,pHW,args->fwd[i]);
    }

    args->bwd[0] = 0.5;
    for (i=1; i<=n; i++)
    {
        ir  = rbuf_kth(&smpl->rbuf, n-i);
        double ci = args->ngenmap ? get_genmap_rate(args, smpl->pos[ir]) : 0.5;

        pAZ = smpl->oaz[ir] * ( (1-args->tHW) * (1-ci) * args->bwd[i-1] + args->tAZ * ci * (1-args->bwd[i-1]) );
        pHW = smpl->ohw[ir] * ( (1-args->tAZ) * (1-ci) * (1-args->bwd[i-1]) + args->tHW * ci * args->bwd[i-1]);

        args->bwd[i] = pAZ / (pAZ+pHW);
        //printf("bwd\ti=%d ir=%d  %d\t oaz=%e ohw=%e \t pAZ=%e  pHW=%e  bwd=%e\n", i,ir,smpl->pos[ir]+1, smpl->oaz[ir],smpl->ohw[ir],pAZ,pHW,args->bwd[i]);
    }

    for (i=1; i<=n; i++)
    {
        ir = rbuf_kth(&smpl->rbuf, i);
        printf("%s\t%s\t%d\t%f\n", args->hdr->samples[args->ismpl[ismpl]], args->hdr->id[BCF_DT_CTG][args->prev_rid].key, smpl->pos[ir]+1, args->fwd[i]*args->bwd[i]);
    }

    smpl->last_az = args->fwd[n]*args->bwd[n];
    rbuf_shift_n(&smpl->rbuf, n);
}

static void flush_buffer_viterbi(args_t *args, int ismpl, int n)
{
    smpl_t *smpl = &args->smpl[ismpl];
    args->viterbi[0].pAZ = smpl->last_az < 0 ? 0.5 : smpl->last_az;

    int i, ir;
    double pAZ, pHW, fromAZ, fromHW;
    for (i=1; i<=n; i++)
    {
        ir = rbuf_kth(&smpl->rbuf, i-1);
        double ci = args->ngenmap ? get_genmap_rate(args, smpl->pos[ir]) : 0.5;

        // P_{i+1}(AZ) = oAZ * max[(1-tHW) * (1-ci) * AZ{i-1} , tAZ * ci * (1-AZ{i-1})]
        fromAZ = (1-args->tHW) * (1-ci) * args->viterbi[i-1].pAZ;
        fromHW = args->tAZ * ci * (1-args->viterbi[i-1].pAZ);
        if ( fromAZ > fromHW )
        {
            pAZ = smpl->oaz[ir] * fromAZ;
            args->viterbi[i].ptr = 0;
        }
        else
        {
            pAZ = smpl->oaz[ir] * fromHW;
            args->viterbi[i].ptr = 1;
        }

        // P_{i+1}(HW) = oHW * max[(1-tAZ) * (1-ci) * (1-AZ{i-1}) , tHW * ci * AZ{i-1}]
        fromHW = (1-args->tAZ) * (1-ci) * (1-args->viterbi[i-1].pAZ);
        fromAZ = args->tHW * ci * args->viterbi[i-1].pAZ;
        if ( fromAZ > fromHW )
        {
            pHW = smpl->ohw[ir] * fromAZ;
        }
        else
        {
            pHW = smpl->ohw[ir] * fromHW;
            args->viterbi[i].ptr |= 2;
        }

        args->viterbi[i].pAZ = pAZ / (pAZ+pHW);
        // printf("viterbi\ti=%d ir=%d  %d\t oaz=%e ohw=%e \t pAZ=%e  pHW=%e  fwd=%e\n", i,ir,smpl->pos[ir]+1, smpl->oaz[ir],smpl->ohw[ir],pAZ,pHW,args->viterbi[i].pAZ);
    }
    // traceback: set ptr to 0 for AZ or 1 for HW
    int mask = args->viterbi[n].pAZ > 0.5 ? 1 : 2;
    for (i=n-1; i>0; i--)
    {
        int ptr = mask & args->viterbi[i+1].ptr;
        args->viterbi[i+1].ptr = ptr;
        mask = args->viterbi[i].ptr & mask ? 2 : 1;
    }
    for (i=1; i<n; i++)
    {
        int AZ = args->viterbi[i].ptr ? 0 : 1;
        double pAZ = AZ ? args->viterbi[i].pAZ : 1 - args->viterbi[i].pAZ;
        ir = rbuf_kth(&smpl->rbuf, i);
        printf("%s\t%s\t%d\t%f\t%d\n", args->hdr->samples[args->ismpl[ismpl]], args->hdr->id[BCF_DT_CTG][args->prev_rid].key, smpl->pos[ir]+1, pAZ,AZ);
    }

    smpl->last_az = args->viterbi[n].pAZ;
    rbuf_shift_n(&smpl->rbuf, n);
    if ( !smpl->rbuf.n ) smpl->last_az = 0.5;
}

static void flush_buffer(args_t *args, int ismpl, int n)
{
    if ( args->counts_only ) flush_counts(args, ismpl, n);
    else if ( args->fwd_bwd ) flush_buffer_fwd_bwd(args, ismpl, n);
    else flush_buffer_viterbi(args, ismpl, n);
}

void set_AF(args_t *args, bcf1_t *line, int32_t *GTs, int nGTs)
{
    int i, j;
    hts_expand(double, line->n_allele, args->mAFs, args->AFs);

    // Get the allele frequencies
    if ( !args->estimate_AF )
    {
        if ( bcf_get_info_int(args->hdr, line, "AN", &args->AN, &args->mAN) != 1 ) 
            error("No AN tag at %s:%d? Use -e to calculate AC,AN on the fly.\n", bcf_seqname(args->hdr,line), line->pos+1);
        int nAC = bcf_get_info_int(args->hdr, line, "AC", &args->ACs, &args->mACs);
        if ( nAC <= 0 ) 
            error("No AC tag at %s:%d? Use -e to calculate AC,AN on the fly.\n", bcf_seqname(args->hdr,line), line->pos+1);

        int nalt = 0; for (i=0; i<line->n_allele-1; i++) nalt += args->ACs[i];    // number of non-ref alleles total
        args->AFs[0] = (double) (args->AN[0] - nalt)/args->AN[0];    // REF frequency
        for (i=1; i<line->n_allele; i++) args->AFs[i] = (double)args->ACs[i-1] / args->AN[0];  // ALT frequencies
    }
    else
    {
        if ( !GTs )
        {
            nGTs = bcf_get_genotypes(args->hdr, line, &args->PLs, &args->mPLs);
            if ( nGTs < 0 ) error("Cannot recalculate AC,AN, GT is not present at %s:%d\n", bcf_seqname(args->hdr,line), line->pos+1);
            if ( nGTs != 2*bcf_hdr_nsamples(args->hdr) ) error("Not diploid at %s:%d?\n", bcf_seqname(args->hdr,line), line->pos+1);
            GTs = args->PLs;
            nGTs /= bcf_hdr_nsamples(args->hdr);
        }

        hts_expand(int32_t, line->n_allele, args->mACs, args->ACs);
        for (i=0; i<line->n_allele; i++) args->ACs[i] = 0;
        if ( args->estimate_AF==1 ) // all samples
        {
            for (i=0; i<bcf_hdr_nsamples(args->hdr); i++)
            {
                int32_t *gt = &args->PLs[i*nGTs];
                for (j=0; j<2; j++)
                {
                    if ( gt[j]==bcf_gt_missing ) continue;
                    if ( gt[j]==bcf_int32_vector_end ) break;
                    args->ACs[bcf_gt_allele(gt[j])]++;
                }
            }
        }
        else                        // subset samples
        {
            for (i=0; i<args->nsmpl; i++)
            {
                int32_t *gt = &args->PLs[args->ismpl[i]*nGTs];
                for (j=0; j<2; j++)
                {
                    if ( gt[0]==bcf_gt_missing || gt[1]==bcf_gt_missing ) continue;
                    args->ACs[bcf_gt_allele(gt[0])]++;
                }
            }
        }
        int ntot = 0; for (i=0; i<line->n_allele; i++) ntot += args->ACs[i]; 
        for (i=0; i<line->n_allele; i++) args->AFs[i] = (double)args->ACs[i] / ntot;  // ALT frequencies
    }
}

int set_pdg_from_PLs(args_t *args, bcf1_t *line)
{
    int nPLs = bcf_get_format_int(args->hdr, line, "PL", &args->PLs, &args->mPLs);
    if ( nPLs!=bcf_hdr_nsamples(args->hdr)*line->n_allele*(line->n_allele+1)/2 ) return -1;
    nPLs /= bcf_hdr_nsamples(args->hdr);

    // Convert PLs to probabilities
    int i, j, k;
    for (i=0; i<args->nsmpl; i++)
    {
        int32_t *pl = &args->PLs[args->ismpl[i]*nPLs], min = pl[0];
        int jmin = 0, kmin = 0, pl_sum = 0;
        for (k=1; k<line->n_allele; k++)
        {
            for (j=0; j<=k; j++)
            {
                pl++;
                if ( *pl==bcf_int32_vector_end ) { k = line->n_allele; break; }
                if ( *pl==bcf_int32_missing ) { continue; }
                pl_sum += *pl;
                if ( *pl < min ) { min = *pl; jmin = j; kmin = k; }
            }
        }
        double *pdg = &args->pdg[i*3];
        if ( !pl_sum ) { pdg[0] = -1; continue; }   // missing data

        pl = &args->PLs[args->ismpl[i]*nPLs];
        double sum = 0;
        for (j=0; j<nPLs; j++)
        {
            assert( pl[j]<256 );
            sum += args->pl2p[ pl[j] ];
        }
        if ( sum==0 ) { pdg[0] = -1; continue; }

        int *als = &args->als[i];
        if ( line->n_allele==2 )
        {
            *als = 1;
            pdg[0] = args->pl2p[ pl[0] ] / sum;
            pdg[1] = args->pl2p[ pl[1] ] / sum;
            pdg[2] = args->pl2p[ pl[2] ] / sum;
            continue;
        }

        if ( jmin==kmin )   // RR or AA
        {
            if ( jmin!=0 )  // AA
                jmin = 0;
            else            // RR
            {
                // Choose ALT at multiallelic sites. Typically, the more frequent ALT will be selected as it comes first 
                pl = &args->PLs[args->ismpl[i]*nPLs]; min = pl[2]; kmin = 1;
                for (k=1; k<line->n_allele; k++)
                {
                    for (j=0; j<=k; j++)
                    {
                        pl++;
                        if ( *pl==bcf_int32_vector_end ) { k = line->n_allele; break; }
                        if ( *pl==bcf_int32_missing ) { continue; }
                    }
                    if ( j==k && *pl < min ) { min = *pl; kmin = k; }
                }
                pl = &args->PLs[args->ismpl[i]*nPLs];
            }
        }

        int idx;
        idx = bcf_alleles2gt(jmin,jmin); pdg[0] = args->pl2p[ pl[idx] ] / sum;
        idx = bcf_alleles2gt(jmin,kmin); pdg[1] = args->pl2p[ pl[idx] ] / sum;
        idx = bcf_alleles2gt(kmin,kmin); pdg[2] = args->pl2p[ pl[idx] ] / sum;
        *als = jmin<<4 | kmin;
    }

    set_AF(args, line, NULL, 0);
    return 0;
}

int set_pdg_from_GTs(args_t *args, bcf1_t *line)
{
    // misusing PL array for GT, which is OK, only the naming is confusing
    int nGTs = bcf_get_genotypes(args->hdr, line, &args->PLs, &args->mPLs);
    if ( nGTs != 2*bcf_hdr_nsamples(args->hdr) ) return -1; // not diploid?

    nGTs /= bcf_hdr_nsamples(args->hdr);

    // Convert GTs to fake probabilities
    int i;
    for (i=0; i<args->nsmpl; i++)
    {
        double *pdg = &args->pdg[i*3];
        int32_t *gt = &args->PLs[args->ismpl[i]*nGTs];

        if (  gt[0]==bcf_gt_missing || gt[1]==bcf_gt_missing ) { pdg[0] = -1; continue; }

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
        int *als = &args->als[i];
        *als = a<<4 | b;
    }
    set_AF(args, line, args->PLs, nGTs);
    return 0;
}

static void vcfroh(args_t *args, bcf1_t *line)
{
    int i;
    if ( !line )
    { 
        for (i=0; i<args->nsmpl; i++)
            flush_buffer(args, i, args->smpl[i].rbuf.n); 
        return; 
    }
    if ( line->rid == args->skip_rid ) return;
    if ( line->n_allele==1 ) return;    // no ALT allele
    if ( args->biallelic_only && line->n_allele!=2 ) return;
    if ( args->snps_only && !bcf_is_snp(line) ) return;

    int skip_rid = 0;
    if ( args->prev_rid<0 ) 
    {
        args->prev_rid = line->rid;
        skip_rid = load_genmap(args, line);
    }
    if ( args->prev_rid!=line->rid )
    {
        for (i=0; i<args->nsmpl; i++)
            flush_buffer(args, i, args->smpl[i].rbuf.n);
        skip_rid = load_genmap(args, line);
    }
    if ( skip_rid )
    {
        fprintf(stderr,"Skipping the sequence: %s\n", bcf_seqname(args->hdr,line));
        args->skip_rid = line->rid;
        return;
    }
    args->prev_rid = line->rid;
    args->ntot++;
    
    int ret;
    if ( !args->fake_PLs ) 
        ret = set_pdg_from_PLs(args, line);
    else
        ret = set_pdg_from_GTs(args, line);

    if ( ret )
    {
        if ( !args->fake_PLs ) error("Could not parse PL field at %s:%d, please run with -G option\n", bcf_seqname(args->hdr,line), line->pos+1);
        error("Could not parse GT field at %s:%d\n", bcf_seqname(args->hdr,line), line->pos+1);
    }
    args->nused++;

    // Calculate emission probabilities P(D|AZ) and P(D|HW)
    for (i=0; i<args->nsmpl; i++)
    {
        double *pdg = &args->pdg[i*3];
        if ( pdg[0]<0 ) continue;   // missing values;

        int *als = &args->als[i];
        int ira  = (*als)>>4;
        int irb  = (*als)&0xf;

        smpl_t *smpl = &args->smpl[i];
        int idx = rbuf_add(&smpl->rbuf);
        double raf = args->AFs[ira];
        double aaf = args->AFs[irb];
        smpl->oaz[idx] = pdg[0]*raf + pdg[2]*aaf;
        smpl->ohw[idx] = pdg[0]*raf*raf + pdg[2]*aaf*aaf + pdg[1]*raf*aaf*2;
        smpl->pos[idx] = line->pos;
        // printf("%d .. als:%d,%d  oaz=%e ohw=%e  raf=%e aaf=%e pdg=%e,%e,%e\n", line->pos+1,ira,irb,smpl->oaz[idx],smpl->ohw[idx],raf,aaf,pdg[0],pdg[1],pdg[2]);
        // printf("%s\t%d\toaz=%e\tohw=%e\t%d\n", args->hdr->id[BCF_DT_CTG][args->prev_rid].key, line->pos+1,smpl->oaz[idx],smpl->ohw[idx], idx);
        if ( smpl->rbuf.n >= args->mwin ) flush_buffer(args, i, smpl->rbuf.n);
    }
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   HMM model for detecting runs of autozygosity.\n");
    fprintf(stderr, "Usage:   bcftools roh [OPTIONS] <in.bcf>|<in.vcf>|<in.vcf.gz>|-\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "    -b, --biallelic-sites              consider only bi-allelic sites\n");
    fprintf(stderr, "    -e, --estimate-AF <all|subset>     calculate AC,AN counts on the fly, using either all samples or samples given via -s\n");
    fprintf(stderr, "    -f, --fwd-bwd                      run forward-backward algorithm instead of Viterbi\n");
    fprintf(stderr, "    -G, --GTs-only <float>             use GTs, ignore PLs, set PL of unseen genotypes to <float>. Safe value to use is 30 to account for GT errors.\n");
    fprintf(stderr, "    -I, --skip-indels                  skip indels as their genotypes are enriched for errors\n");
    fprintf(stderr, "    -m, --genetic-map <file>           genetic map in IMPUTE2 format, single file or mask, where string \"{CHROM}\" is replaced with chromosome name\n");
    fprintf(stderr, "    -r, --regions <reg|file>           restrict to comma-separated list of regions or regions listed in a file, see man page for details\n");
    fprintf(stderr, "    -s, --samples <list|file>          list of samples (file or comma separated list) [null]\n");
    fprintf(stderr, "    -t, --targets <reg|file>           similar to -r but streams rather than index-jumps, see man page for details\n");
    fprintf(stderr, "    -w, --win <int>                    maximum window length [100_000]\n");
    fprintf(stderr, "HMM Options:\n");
    fprintf(stderr, "    -a, --hw-to-az <float>             P(AZ|HW) transition probability from AZ (autozygous) to HW (Hardy-Weinberg) state [1e-4]\n");
    fprintf(stderr, "    -H, --az-to-hw <float>             P(HW|AZ) transition probability from HW to AZ state [1e-3]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfroh(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->tAZ     = 1e-4;
    args->tHW     = 1e-3;
    args->mwin    = (int)1e5;   // maximum number of sites that can be processed in one go

    static struct option loptions[] = 
    {
        {"estimate-AF",1,0,'e'},
        {"GTs-only",1,0,'G'},
        {"counts-only",0,0,'c'},
        {"win",1,0,'w'},
        {"samples",1,0,'s'},
        {"hw-to-az",1,0,'a'},
        {"az-to-hw",1,0,'H'},
        {"targets",1,0,'t'},
        {"regions",1,0,'r'},
        {"genetic-map",1,0,'m'},
        {"fwd-bwd",1,0,'f'},
        {"biallelic-sites",0,0,'b'},
        {"skip-indels",0,0,'I'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h?r:t:H:a:w:s:cm:fG:bIa:e:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'e': 
                if (!strcmp("all",optarg)) args->estimate_AF = 1; 
                else if  (!strcmp("subset",optarg)) args->estimate_AF = 2;
                else error("Expected 'all' or 'subset' with -s\n");
            case 'b': args->biallelic_only = 1; break;
            case 'I': args->snps_only = 1; break;
            case 'G': args->fake_PLs = 1; args->unseen_PL = pow(10,-atof(optarg)/10.); break;
            case 'm': args->genmap_fname = optarg; break;
            case 'c': args->counts_only = 1; break;
            case 'f': args->fwd_bwd = 1; break;
            case 's': args->samples_fname = optarg; break;
            case 'w': args->mwin = (int)atof(optarg); break;
            case 'a': args->tAZ = atof(optarg); break;
            case 'H': args->tHW = atof(optarg); break;
            case 't': args->targets_fname = optarg; break;
            case 'r': args->regions_fname = optarg; break;
            case 'h': 
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( argc<optind+1 ) usage(args);
    if ( args->regions_fname )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_fname)<0 )
            error("Failed to read the regions: %s\n", args->regions_fname);
    }
    if ( args->targets_fname )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_fname, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_fname);
    }
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    
    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        vcfroh(args, args->files->readers[0].buffer[0]);
    }
    vcfroh(args, NULL);
    fprintf(stderr,"Number of lines: total/processed: %d/%d\n", args->ntot,args->nused);
    destroy_data(args);
    free(args);
    return 0;
}
