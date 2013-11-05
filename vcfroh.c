#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"
#include "rbuf.h"


typedef struct
{
    rbuf_t rbuf;
    double *ohw, *oaz;  // P(D|HW) and P(D|AZ)
    uint32_t *pos;
    int j_al, k_al;
}
smpl_t;

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    double az_hw, hw_az;    // P(AZ|HW), P(HW|AZ)
    double *fwd, *bwd;      // HMM forward and backward autozygosity probs scaled to hw+az=1

    int nsmpl, *ismpl;
    int mwin;
    smpl_t *smpl;
    int32_t *PLs, *AN, *ACs;
    double *AFs;
    double pl2p[256], *pdg;
    int mPLs, mAFs, mAN, mACs, mpdg;

    char **argv, *targets_fname, *regions_fname, *samples_fname;
    int argc, counts_only;
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
    args->hdr = args->files->readers[0].header;
    if ( !args->samples_fname ) args->samples_fname = "-";
    if ( !bcf_sr_set_samples(args->files, args->samples_fname) )
        error("Error: could not set the samples %s\n", args->samples_fname);
    args->nsmpl = args->files->n_smpl;
    args->ismpl = args->files->readers[0].samples;

    args->smpl  = (smpl_t*) smalloc(sizeof(smpl_t)*args->nsmpl);
    args->fwd   = (double*) smalloc(sizeof(double)*args->mwin);
    args->bwd   = (double*) smalloc(sizeof(double)*args->mwin);
    int i;
    for (i=0; i<args->nsmpl; i++)
    {
        smpl_t *smpl = &args->smpl[i];
        rbuf_init(&smpl->rbuf, args->mwin);
        smpl->ohw = (double*) smalloc(sizeof(double)*args->mwin);
        smpl->oaz = (double*) smalloc(sizeof(double)*args->mwin);
        smpl->pos = (uint32_t*) smalloc(sizeof(uint32_t)*args->mwin);
    }
    for (i=0; i<256; i++)
        args->pl2p[i] = pow(10., -i/10.);
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
    free(args->smpl);
    free(args->fwd); free(args->bwd);
    free(args->PLs); free(args->AFs); free(args->pdg);
    free(args->AN); free(args->ACs);
}

static void flush_counts(args_t *args, int ismpl, int n)
{
    smpl_t *smpl = &args->smpl[ismpl];
    
    int ir  = smpl->rbuf.f;
    int pos = smpl->pos[ smpl->rbuf.f ];
    int nhw = 0, naz = 0;

    for (ir=-1; rbuf_next(&smpl->rbuf,&ir); )
    {
        if ( smpl->ohw[ir] > smpl->oaz[ir] ) nhw++;
        else naz++;
        if ( smpl->pos[ir] - pos > args->mwin  )
        {
            double az_rate = (double)naz*args->mwin/(smpl->pos[ir] - pos);
            double hw_rate = (double)nhw*args->mwin/(smpl->pos[ir] - pos);
            printf("%s\t%d\t%e\t%e\n", args->hdr->samples[args->ismpl[ismpl]], pos+1, az_rate, hw_rate);
            pos = smpl->pos[ir];
            nhw = naz = 0;
        }
    }
    rbuf_shift_n(&smpl->rbuf, n);
}

/**
 *  This function implements the HMM model:
 *      D = Data, AZ = autozygosity, HW = Hardy-Weinberg
 *      oaz = P(D|AZ)
 *      ohw = P(D|HW)
 *      hw_az = P(HW|AZ)
 *      az_hw = P(AZ|HW)
 *      P(HW|HW) = 1 - az_hw
 *      P(AZ|AZ) = 1 - hw_az
 */
static void flush_buffer(args_t *args, int ismpl, int n)
{
    if ( args->counts_only )
    {
        flush_counts(args, ismpl, n);
        return;
    }

    smpl_t *smpl = &args->smpl[ismpl];

    int ir = rbuf_kth(&smpl->rbuf, 0);
    double faz = 0.5 * (args->az_hw - args->hw_az + 1.0) * smpl->oaz[ir];
    double fhw = 0.5 * (args->hw_az - args->az_hw + 1.0) * smpl->ohw[ir];
    args->fwd[0] = faz / (faz + fhw);

    int i;
    for (i=1; i<n; i++)
    {
        ir  = rbuf_kth(&smpl->rbuf, i);
        faz = (args->fwd[i-1]*(1 - args->hw_az) + (1 - args->fwd[i-1])*args->az_hw) * smpl->oaz[ir];
        fhw = (args->fwd[i-1]*args->hw_az + (1-args->fwd[i-1])*(1 - args->az_hw)) * smpl->ohw[ir];
        args->fwd[i] = faz / (faz + fhw);
    }

    ir = rbuf_kth(&smpl->rbuf, n-1);
    faz = 0.5 * (args->az_hw - args->hw_az + 1.0) * smpl->oaz[ir];
    fhw = 0.5 * (args->hw_az - args->az_hw + 1.0) * smpl->ohw[ir];
    args->bwd[0] = faz / (faz + fhw);
    for (i=1; i<n; i++)
    {
        ir  = rbuf_kth(&smpl->rbuf, n-1-i);
        faz = (args->bwd[i-1]*(1 - args->hw_az) + (1 - args->bwd[i-1])*args->az_hw) * smpl->oaz[ir];
        fhw = (args->bwd[i-1]*args->hw_az + (1-args->bwd[i-1])*(1 - args->az_hw)) * smpl->ohw[ir];
        args->bwd[i] = faz / (faz + fhw);
    }

    for (i=0; i<n; i++)
        printf("%s\t%d\t%f\n", args->hdr->samples[args->ismpl[ismpl]], smpl->pos[i]+1, args->fwd[i]*args->bwd[i]);

    rbuf_shift_n(&smpl->rbuf, n);
}

static void vcfroh(args_t *args, bcf1_t *line)
{
    int i, j, k;
    if ( !line )
    { 
        for (i=0; i<args->nsmpl; i++)
            flush_buffer(args, i, args->smpl[i].rbuf.n); 
        return; 
    }
    if ( line->n_allele==1 || line->d.allele[1][0]=='X' ) return;    // skip ref-only sites
    
    // Get the PLs
    int nPLs = bcf_get_format_int(args->hdr, line, "PL", &args->PLs, &args->mPLs);
    if ( nPLs!=bcf_hdr_nsamples(args->hdr)*line->n_allele*(line->n_allele+1)/2 ) return;
    nPLs /= bcf_hdr_nsamples(args->hdr);

    // Get the allele frequencies and convert PLs to probabilities
    hts_expand(double, nPLs*args->nsmpl, args->mpdg, args->pdg);
    if ( bcf_get_info_int(args->hdr, line, "AN", &args->AN, &args->mAN) != 1 ) error("No AN tag at %s:%d?\n", bcf_seqname(args->hdr,line), line->pos+1);
    int nAC = bcf_get_info_int(args->hdr, line, "AC", &args->ACs, &args->mACs);
    if ( nAC <= 0 ) error("No AC tag at %s:%d?\n", bcf_seqname(args->hdr,line), line->pos+1);
    if ( nAC < 1 ) return;  // skip
    hts_expand(double, line->n_allele, args->mAFs, args->AFs);
    int nalt = 0; for (i=0; i<line->n_allele-1; i++) nalt += args->ACs[i];    // number of non-ref alleles total
    args->AFs[0] = (double) (args->AN[0] - nalt)/args->AN[0];    // REF frequency
    for (i=1; i<line->n_allele; i++) args->AFs[i] = (double)args->ACs[i] / args->AN[0];  // ALT frequencies

    for (i=0; i<args->nsmpl; i++)
    {
        int32_t *pl = &args->PLs[args->ismpl[i]*nPLs], min = pl[0];
        int jmin = 0, kmin = 0;
        for (k=1; k<line->n_allele; k++)
        {
            for (j=0; j<=k; j++)
            {
                pl++;
                if ( *pl==bcf_int32_vector_end ) { j = line->n_allele; break; }
                if ( *pl==bcf_int32_missing ) { continue; }
                if ( *pl < min ) { min = *pl; jmin = j; kmin = k; }
            }
        }
        if ( jmin==kmin )
        {
            if ( jmin==0 ) 
            {
                // skip ref call
                args->pdg[i*nPLs] = -1;
                continue;
            }
            jmin = 0;   // homALT is most likely; leave k at ALT and set jmin to REF
        }
        args->smpl[i].j_al = jmin;
        args->smpl[i].k_al = kmin;

        double *pdg = &args->pdg[i*nPLs];
        pl = &args->PLs[args->ismpl[i]*nPLs];
        double sum = 0;
        for (j=0; j<nPLs; j++)
        {
            assert( pl[j]<256 );
            pdg[j] = args->pl2p[ pl[j] ];
            sum += pdg[j];
        }
        // Normalize: sum_i pdg_i = 1
        if ( sum!=nPLs )
            for (j=0; j<nPLs; j++) pdg[j] /= sum;
        else
            pdg[0] = -1;
    }

    // Calculate emission probabilities P(D|AZ) and P(D|HW)
    for (i=0; i<args->nsmpl; i++)
    {
        double *pdg = &args->pdg[i*nPLs];
        if ( pdg[0]<0 ) continue;   // missing values;

        smpl_t *smpl = &args->smpl[i];
        int idx = rbuf_add(&smpl->rbuf);
        j = smpl->j_al;
        k = smpl->k_al; 
        smpl->oaz[idx] = pdg[(j+1)*(j+2)/2-1]*args->AFs[j] + pdg[(k+1)*(k+2)/2-1]*args->AFs[k];
        smpl->ohw[idx] = pdg[(j+1)*(j+2)/2-1]*args->AFs[j]*args->AFs[j] + pdg[(k+1)*(k+2)/2-1]*args->AFs[k]*args->AFs[k] + pdg[k*(k+1)/2+j]*args->AFs[j]*args->AFs[k]*2;
        smpl->pos[idx] = line->pos;
        //printf("%d\t%e\t%e\t QHet=%f  frac=%e %e (%d,%d)\n", line->pos+1,smpl->oaz[idx],smpl->ohw[idx], -4.3429*log(pdg[k*(k+1)/2+j]), args->AFs[j],args->AFs[k],j,k);
        if ( smpl->rbuf.n + 1 >= args->mwin ) flush_buffer(args, i, smpl->rbuf.n);
    }
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   HMM model for detecting runs of autozygosity.\n");
    fprintf(stderr, "Usage:   bcftools roh [OPTIONS] <in.bcf>|<in.vcf>|<in.vcf.gz>|-\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "    -c, --counts-only              no HMM, simply report counts of HETs and HOMs per win\n");
    fprintf(stderr, "    -r, --regions <reg|file>       same as -t but index-jumps rather than streams to a region (requires indexed VCF/BCF)\n");
    fprintf(stderr, "    -s, --samples <list|file>      list of samples (file or comma separated list) [null]\n");
    fprintf(stderr, "    -t, --targets <reg|file>       restrict to positions in tab-delimited tabix indexed file <chr,pos> or <chr,from,to>, 1-based, inclusive\n");
    fprintf(stderr, "    -w, --win <int>                maximum window length [100_000]\n");
    fprintf(stderr, "HMM Options:\n");
    fprintf(stderr, "    -a, --autozygosity <float>     P(AZ|HW) transition probability [3.8e-9]\n");
    fprintf(stderr, "    -H, --Hardy-Weinberg <float>   P(HW|AZ) transition probability [6.0e-8]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfroh(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->az_hw   = 3.8e-9;
    args->hw_az   = 6.0e-8;
    args->mwin    = (int)1e5;   // maximum number of sites that can be processed in one go

    static struct option loptions[] = 
    {
        {"counts-only",0,0,'c'},
        {"win",1,0,'w'},
        {"samples",1,0,'s'},
        {"autozygosity",1,0,'a'},
        {"Hardy-Weinberg",1,0,'H'},
        {"targets",1,0,'t'},
        {"regions",1,0,'r'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h?r:t:H:a:w:s:c",loptions,NULL)) >= 0) {
        switch (c) {
            case 'c': args->counts_only = 1; break;
            case 's': args->samples_fname = optarg; break;
            case 'w': args->mwin = (int)atof(optarg); break;
            case 'a': args->az_hw = atof(optarg); break;
            case 'H': args->hw_az = atof(optarg); break;
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
    if ( args->counts_only ) printf("# [1]Sample\t[2]Position\t[3]HOM rate\t[4]HET rate\n");
    while ( bcf_sr_next_line(args->files) )
    {
        vcfroh(args, args->files->readers[0].buffer[0]);
    }
    vcfroh(args, NULL);
    destroy_data(args);
    free(args);
    return 0;
}
