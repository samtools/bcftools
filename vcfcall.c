#include <stdarg.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <htslib/vcf.h>
#include <time.h>
#include <zlib.h>
#include <stdarg.h>
#include <htslib/kfunc.h>
#include <htslib/synced_bcf_reader.h>
#include <sys/stat.h>
#include "bcftools.h"
#include "call.h"
#include "prob1.h"

void error(const char *format, ...);

#ifdef _WIN32
#define srand48(x) srand(x)
#define lrand48() rand()
#endif

#include <htslib/kseq.h>
KSTREAM_INIT(gzFile, gzread, 16384)

#define CF_NO_GENO      1
//                      (1<<1)
#define CF_CCALL        (1<<2)
//                      (1<<3)
//                      (1<<4)
//                      (1<<5)
#define CF_ACGT_ONLY    (1<<6)
#define CF_QCALL        (1<<7)
#define CF_ADJLD        (1<<8)
#define CF_NO_INDEL     (1<<9)
#define CF_ANNO_MAX     (1<<10)
#define CF_MCALL        (1<<11)
#define CF_PAIRCALL     (1<<12)
#define CF_QCNT         (1<<13)
#define CF_INDEL_ONLY   (1<<14)

typedef struct 
{
    int flag;   // combination of CF_* flags above
    int output_type;
    htsFile *bcf_in, *out_fh;
    char *bcf_fname;
    char **samples;             // for subsampling and ploidy
    int nsamples, *samples_map;
    char *regions, *targets;    // regions to process

    call_t aux;     // parameters and temporary data

    int argc;
    char **argv;

	//  int flag, prior_type, n1, n_sub, *sublist, n_perm;
	//  uint32_t *trio_aux;
	//  char *prior_file, **subsam;
	//  uint8_t *ploidy;
	//  double theta, pref, indel_frac, min_smpl_frac, min_lrt, min_ma_lrt;
    // Permutation tests
    //  int n_perm, *seeds;
    //  double min_perm_p;
	//  void *bed;
}
args_t;


static family_t *get_family(family_t *fam, int nfam, char *name)
{
    int i;
    for (i=0; i<nfam; i++)
    {
        if ( !strcmp(fam[i].name, name) ) return &fam[i];
    }
    return NULL;
}

static char **add_sample(char **sam, int *n, int *m, char *name, int ploidy, int *ith)
{
    int i;
    for (i=0; i<*n; i++)
    {
        if ( !strcmp(sam[i], name) ) 
        {
            *ith = i;
            return sam;
        }
    }
    hts_expand(char*,(*n+1),*m,sam);
    int len = strlen(name);
    sam[*n] = (char*) malloc(len+2);
    memcpy(sam[*n],name,len+1);
    sam[*n][len+1] = ploidy;
    *ith = *n;
    (*n)++;
    return sam;
}

static char **read_ped_samples(call_t *call, const char *fn, int *_n)
{
	char **sam = 0;
    int dret, n = 0, max = 0, i;
    kstream_t *ks;
    kstring_t s = {0,0,0};
    gzFile fp;
    fp = gzopen(fn, "r");
    if ( !fp ) error("Could not read the file: %s\n", fn);
    ks = ks_init(fp);
    while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) 
    {
        char *col_ends[5], *tmp = s.s;
        i = 0;
        while ( *tmp && i<5 )
        {
            if ( *tmp=='\t' || *tmp==' ' ) 
            { 
                col_ends[i] = tmp;
                *tmp = 0;
                i++; 
            }
            tmp++;
        }
        if ( i!=5 ) break;

        family_t *fam = get_family(call->fams, call->nfams, s.s);
        if ( !fam )
        {
            call->nfams++;
            hts_expand(family_t, call->nfams, call->mfams, call->fams);
            fam = &call->fams[call->nfams-1];
            fam->name = strdup(s.s);
            for (i=0; i<3; i++) fam->sample[i] = -1;
        } 
        
        int ploidy = 2;
        if ( call->flag & (CALL_CHR_X|CALL_CHR_Y) )
        {
            if ( col_ends[3][1]=='1' ) ploidy = 1; // male: one chrX and one chrY copy
            else
                ploidy = call->flag & CALL_CHR_X ? 2 : 0; // female: two chrX copies, no chrY
        }
        sam = add_sample(sam, &n, &max, col_ends[0]+1, ploidy, &i);
        if ( col_ends[2]-col_ends[1]!=2 || col_ends[1][1]!='0' )    // father
        {
            if ( fam->sample[CHILD]>=0 ) error("Multiple childs in %s\n", s.s);
            fam->sample[CHILD] = i;
            if ( fam->sample[FATHER]>=0 ) error("Two fathers in %s?\n", s.s);
            sam = add_sample(sam, &n, &max, col_ends[1]+1, call->flag & (CALL_CHR_X|CALL_CHR_Y) ? 1 : 2, &fam->sample[FATHER]); 
        }
        if ( col_ends[3]-col_ends[2]!=2 || col_ends[2][1]!='0' )    // mother
        {
            if ( fam->sample[MOTHER]>=0 ) error("Two mothers in %s?\n", s.s);
            sam = add_sample(sam, &n, &max, col_ends[2]+1, call->flag & CALL_CHR_Y ? 0 : 2, &fam->sample[MOTHER]); 
        }
    }
    for (i=0; i<call->nfams; i++)
        assert( call->fams[i].sample[0]>=0 && call->fams[i].sample[1]>=0 && call->fams[i].sample[2]>=0 ); //todo

    ks_destroy(ks);
    gzclose(fp);
    free(s.s);
    *_n = n;
    return sam;
}


/*
 *  Reads sample names and their ploidy (optional) from a file.
 *  Alternatively, if no such file exists, the file name is interpreted
 *  as a comma-separated list of samples. When ploidy is not present,
 *  the default ploidy 2 is assumed.
 *
 *  Returns an array of sample names, where the byte value just after \0
 *  indicates the ploidy.
 */
static char **read_samples(call_t *call, const char *fn, int *_n)
{
	int dret, n = 0, max = 0;
	char **sam = 0;
	*_n = 0;

    struct stat sbuf;
    if ( stat(fn, &sbuf) != 0  )
    {
        // it is not a file, interpret as list of sample names
        const char *t = fn, *p = t;
        while (*t)
        {
            t++;
            if ( *t==',' || !*t )
            { 
                sam = (char**) realloc(sam, sizeof(char*)*(n+1));
                sam[n] = (char*) malloc(sizeof(char)*(t-p+2));
                memcpy(sam[n], p, t-p);
                sam[n][t-p]   = 0;
                sam[n][t-p+1] = 2;    // assume diploid
                p = t+1;
                n++; 
            }
        }
        *_n = n;
        return sam;
    }

    sam = read_ped_samples(call, fn, _n);
    if ( sam ) return sam;

	kstream_t *ks;
	kstring_t s = {0,0,0};
	gzFile fp;
	fp = gzopen(fn, "r");
    if ( !fp ) error("Could not read the file: %s\n", fn);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_SPACE, &s, &dret) >= 0) 
    {
        hts_expand(char*,(n+1),max,sam);
		int l = s.l;
		sam[n] = (char*) malloc(sizeof(char)*(s.l+2));
		strcpy(sam[n], s.s);
		sam[n][l+1] = 2; // by default, diploid
		if (dret != '\n') {
			if (ks_getuntil(ks, KS_SEP_SPACE, &s, &dret) >= 0) { // read ploidy, 1 or 2
				int x = (int)s.s[0] - '0'; // Convert ASCII digit to decimal
				if (x == 0 || x == 1 || x == 2) sam[n][l+1] = x;
				else fprintf(stderr, "(%s) ploidy can only be 0, 1 or 2; assuming diploid\n", __func__);
			}
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
		}
		++n;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(s.s);
	*_n = n;
	return sam;
}


static void init_data(args_t *args)
{
    args->aux.srs = bcf_sr_init();

    // Open files for input and output, initialize structures
    if ( args->targets )
    {
        if ( bcf_sr_set_targets(args->aux.srs, args->targets, args->aux.flag&CALL_CONSTR_ALLELES ? 3 : 0)<0 )
            error("Failed to read the targets: %s\n", args->targets);
    }
    if ( args->regions )
    {
        if ( bcf_sr_set_regions(args->aux.srs, args->regions)<0 )
            error("Failed to read the targets: %s\n", args->regions);
    }
    
    int i;
    if ( !bcf_sr_add_reader(args->aux.srs, args->bcf_fname) ) error("Failed to open: %s\n", args->bcf_fname);

    if ( args->nsamples && args->nsamples != args->aux.srs->readers[0].header->n[BCF_DT_SAMPLE] )
    {
        args->samples_map = (int *) malloc(sizeof(int)*args->nsamples);
        args->aux.hdr = bcf_hdr_subset(args->aux.srs->readers[0].header, args->nsamples, args->samples, args->samples_map);
        for (i=0; i<args->nsamples; i++)
            if ( args->samples_map[i]<0 ) fprintf(stderr,"Warning: no such sample: \"%s\"\n", args->samples[i]);
    }
    else
        args->aux.hdr = bcf_hdr_dup(args->aux.srs->readers[0].header);

    // Reorder ploidy and family indexes to match mpileup's output and exclude samples which are not available
    if ( args->aux.ploidy )
    {
        for (i=0; i<args->aux.nfams; i++)
        {
            int j;
            for (j=0; j<3; j++)
            {
                int k = bcf_hdr_id2int(args->aux.hdr, BCF_DT_SAMPLE, args->samples[ args->aux.fams[i].sample[j] ]);
                if ( k<0 ) error("No such sample: %s\n", args->samples[ args->aux.fams[i].sample[j] ]);
                args->aux.fams[i].sample[j] = k;
            }
        }
        uint8_t *ploidy = (uint8_t*) calloc(args->aux.hdr->n[BCF_DT_SAMPLE], 1);
        for (i=0; i<args->nsamples; i++)    // i index in -s sample list
        {
            int j = bcf_hdr_id2int(args->aux.hdr, BCF_DT_SAMPLE, args->samples[i]);     // j index in the output VCF / subset VCF
            if ( j<0 ) 
            {
                fprintf(stderr,"Warning: no such sample: \"%s\"\n", args->samples[i]);
                continue;
            }
            ploidy[j] = args->aux.ploidy[i];
        }
        args->nsamples = args->aux.hdr->n[BCF_DT_SAMPLE];
        for (i=0; i<args->nsamples; i++)
            assert( ploidy[i]==0 || ploidy[i]==1 || ploidy[i]==2 );
        free(args->aux.ploidy);
        args->aux.ploidy = ploidy;
    }

    args->out_fh = hts_open("-", hts_bcf_wmode(args->output_type));

    if ( args->flag & CF_QCALL ) 
        return;

    if ( args->flag & CF_MCALL ) 
        mcall_init(&args->aux);

    if ( args->flag & CF_CCALL )
        ccall_init(&args->aux);

    bcf_hdr_append_version(args->aux.hdr, args->argc, args->argv, "bcftools_call");
    bcf_hdr_fmt_text(args->aux.hdr);
    bcf_hdr_write(args->out_fh, args->aux.hdr);
}

static void destroy_data(args_t *args)
{
    if ( args->flag & CF_CCALL ) ccall_destroy(&args->aux);
    else if ( args->flag & CF_MCALL ) mcall_destroy(&args->aux);
    else if ( args->flag & CF_QCALL ) qcall_destroy(&args->aux);
    int i;
    for (i=0; i<args->nsamples; i++) free(args->samples[i]);
    if ( args->aux.fams )
    {
        for (i=0; i<args->aux.nfams; i++) free(args->aux.fams[i].name);
        free(args->aux.fams);
    }
    free(args->samples);
    free(args->samples_map);
    free(args->aux.ploidy);
    bcf_hdr_destroy(args->aux.hdr);
    hts_close(args->out_fh);
    bcf_sr_destroy(args->aux.srs);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About: This command replaces the former \"bcftools view\" caller. Some of the original functionality has been\n");
    fprintf(stderr, "       temporarily lost in the process of transition under htslib, but will be added back on popular demand. The original\n");
    fprintf(stderr, "       calling model can be invoked with the -c option. Note that we use the new multiallelic -m caller by default,\n");
    fprintf(stderr, "       therefore -c is not as well tested as -m. If you encounter bugs, please do let us know.\n");
    fprintf(stderr, "Usage: bcftools call [options] <in.bcf|in.vcf|in.vcf.gz> [reg]\n");
    fprintf(stderr, "File format options:\n");
    fprintf(stderr, "   -o, --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "   -r, --region <reg|file>         restrict to comma-separated list of regions or regions listed in tab-delimited indexed file\n");
    fprintf(stderr, "   -s, --samples <list|file>       sample list, PED file or a file with optional second column for ploidy (0, 1 or 2) [all samples]\n");
    fprintf(stderr, "   -t, --targets <reg|file>        same as -r but streams rather than index-jumps to it. Coordinates are 1-based, inclusive\n");
    fprintf(stderr, "\nInput/output options:\n");
    fprintf(stderr, "   -A, --keep-alts                 keep all possible alternate alleles at variant sites\n");
    fprintf(stderr, "   -N, --skip-Ns                   skip sites where REF is not A/C/G/T\n");
    fprintf(stderr, "   -S, --skip <snps|indels>        skip indels/snps\n");
    fprintf(stderr, "   -v, --variants-only             output variant sites only\n");
    fprintf(stderr, "\nConsensus/variant calling options:\n");
    fprintf(stderr, "   -c, --consensus-caller          the original calling method (conflicts with -m)\n");
    fprintf(stderr, "   -C, --constrain <str>           one of: alleles, trio (see manual)\n");
    fprintf(stderr, "   -m, --multiallelic-caller       alternative model for multiallelic and rare-variant calling (conflicts with -c)\n");
    fprintf(stderr, "   -p, --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5] or another allele accepted if P(chi^2)>=1-FLOAT with -m [1e-2]\n");
    fprintf(stderr, "   -X, --chromosome-X              haploid output for male samples (requires PED file with -s)\n");
    fprintf(stderr, "   -Y, --chromosome-Y              haploid output for males and skips females (requires PED file with -s)\n");

    // todo (and more)
    // fprintf(stderr, "\nContrast calling and association test options:\n");
    // fprintf(stderr, "       -1 INT    number of group-1 samples [0]\n");
    // fprintf(stderr, "       -C FLOAT  posterior constrast for LRT<FLOAT and P(ref|D)<0.5 [%g]\n", args->aux.min_lrt);
    // fprintf(stderr, "       -U INT    number of permutations for association testing (effective with -1) [0]\n");
    // fprintf(stderr, "       -X FLOAT  only perform permutations for P(chi^2)<FLOAT [%g]\n", args->aux.min_perm_p);
    fprintf(stderr, "\n");
    exit(-1);
}

int main_vcfcall(int argc, char *argv[])
{
    char *samples_fname = NULL;
    args_t args;
	memset(&args, 0, sizeof(args_t));
    args.argc = argc; args.argv = argv;
    args.aux.prior_type = -1;
    args.aux.indel_frac = -1;
    args.aux.theta      = 1e-3;
    args.aux.pref       = 0.5;
    args.aux.min_perm_p = 0.01;
    args.aux.min_lrt    = 1;
    args.aux.min_ma_lrt = 1 - 1e-2;

    float p_arg = -1;
    int i, c;

    static struct option loptions[] = 
    {
        {"help",0,0,'h'},
        {"output-type",1,0,'o'},
        {"region",1,0,'r'},
        {"samples",1,0,'s'},
        {"targets",1,0,'t'},
        {"keep-alts",0,0,'A'},
        {"skip-Ns",0,0,'N'},
        {"skip",1,0,'S'},
        {"variants-only",0,0,'v'},
        {"consensus-caller",0,0,'c'},
        {"constrain",1,0,'C'},
        {"multiallelic-caller",0,0,'m'},
        {"pval-threshold",1,0,'p'},
        {"chromosome-X",0,0,'X'},
        {"chromosome-Y",0,0,'Y'},
        {0,0,0,0}
    };

	while ((c = getopt_long(argc, argv, "h?o:r:s:t:ANS:vcmp:C:XY", loptions, NULL)) >= 0) 
    {
		switch (c) 
        {
            case 'N': args.flag |= CF_ACGT_ONLY; break;                 // omit sites where first base in REF is N
            case 'A': args.aux.flag |= CALL_KEEPALT; break;
            case 'c': args.flag |= CF_CCALL; break;          // the original EM based calling method
            case 'v': args.aux.flag |= CALL_VARONLY; break;
            case 'o': 
                      switch (optarg[0]) {
                          case 'b': args.output_type = FT_BCF_GZ; break;
                          case 'u': args.output_type = FT_BCF; break;
                          case 'z': args.output_type = FT_VCF_GZ; break;
                          case 'v': args.output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'C': 
                      if ( !strcasecmp(optarg,"alleles") ) args.aux.flag |= CALL_CONSTR_ALLELES;
                      else if ( !strcasecmp(optarg,"trio") ) args.aux.flag |= CALL_CONSTR_TRIO;
                      else error("Unknown argument to -C: \"%s\"\n", optarg);
                      break;
            case 'X': args.aux.flag |= CALL_CHR_X; break;
            case 'Y': args.aux.flag |= CALL_CHR_Y; break;
            case 'S': 
                      if ( !strcasecmp(optarg,"snps") ) args.flag |= CF_INDEL_ONLY;
                      else if ( !strcasecmp(optarg,"indels") ) args.flag |= CF_NO_INDEL;
                      else error("Unknown argument to -I: \"%s\"\n", optarg);
            case 'm': args.flag |= CF_MCALL; break;         // multiallelic calling method
            case 'p': p_arg = atof(optarg); break;
            case 'r': args.regions = optarg; break;
            case 't': args.targets = optarg; break;
            case 's': samples_fname = optarg; break;
            default: usage(&args);
        }
    }
	if (argc == optind) usage(&args); 
    args.bcf_fname = argv[optind++];

    // Sanity check options and initialize
    if ( samples_fname )
    {
        args.samples = read_samples(&args.aux, samples_fname, &args.nsamples);
        args.aux.ploidy = (uint8_t*) calloc(args.nsamples+1, 1);
        for (i=0; i<args.nsamples; i++) args.aux.ploidy[i] = args.samples[i][strlen(args.samples[i]) + 1];
    }
    if ( (args.flag & CF_CCALL ? 1 : 0) + (args.flag & CF_MCALL ? 1 : 0) + (args.flag & CF_QCALL ? 1 : 0) > 1 ) error("Only one of -c or -m options can be given\n");
    if ( !(args.flag & CF_CCALL) && !(args.flag & CF_MCALL) && !(args.flag & CF_QCALL) ) error("Expected -c or -m option\n");
	if ( args.aux.n_perm && args.aux.ngrp1_samples<=0 ) error("Expected -1 with -U\n");    // not sure about this, please fix
    if ( p_arg!=-1 ) { args.aux.pref = p_arg; args.aux.min_ma_lrt = 1 - p_arg; }  // only one is actually used
    if ( args.aux.flag & CALL_CONSTR_ALLELES )
    {
        if ( !args.targets ) error("Expected -t with \"-C alleles\"\n");
        if ( !(args.flag & CF_MCALL) ) error("The \"-C alleles\" mode requires -m\n");
    }
    if ( args.aux.flag & CALL_CHR_X && args.aux.flag & CALL_CHR_Y ) error("Only one of -X or -Y should be given\n");

    init_data(&args);

    while ( bcf_sr_next_line(args.aux.srs) )
    {
        bcf1_t *bcf_rec = args.aux.srs->readers[0].buffer[0];
        if ( args.samples_map ) bcf_subset(args.aux.hdr, bcf_rec, args.nsamples, args.samples_map);
        bcf_unpack(bcf_rec, BCF_UN_STR);

        // Skip unwanted sites
        if ( args.aux.flag & CALL_VARONLY )
        {
            if ( bcf_rec->n_allele==1 ) continue;                                       // not a variant
            if ( bcf_rec->n_allele==2 && bcf_rec->d.allele[1][0]=='X' ) continue;       // second allele is mpileup's X, not a variant
        }
        if ( (args.flag & CF_INDEL_ONLY) && bcf_is_snp(bcf_rec) ) continue;    // not an indel
        if ( (args.flag & CF_NO_INDEL) && !bcf_is_snp(bcf_rec) ) continue;     // not a SNP
        if ( (args.flag & CF_ACGT_ONLY) && (bcf_rec->d.allele[0][0]=='N' || bcf_rec->d.allele[0][0]=='n') ) continue;   // REF[0] is 'N'

        bcf_unpack(bcf_rec, BCF_UN_ALL);

        // Various output modes: QCall output (todo)
        if ( args.flag & CF_QCALL ) 
        {
            qcall(&args.aux, bcf_rec);
            continue;
        }

        // Calling modes which output VCFs
        int ret;
        if ( args.flag & CF_MCALL )
            ret = mcall(&args.aux, bcf_rec);
        else
            ret = ccall(&args.aux, bcf_rec);

        if ( ret==-1 ) error("Something is wrong\n");
        if ( (args.aux.flag & CALL_VARONLY) && ret==0 ) continue;     // not a variant

        // Output
        bcf_write1(args.out_fh, args.aux.hdr, bcf_rec);
    }
    destroy_data(&args);
	return 0;
}

