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
#define CF_VCFIN        (1<<4)
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
    bcf_srs_t *files;           // synced reader
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


/*
 *  Reads sample names and their ploidy (optional) from a file.
 *  Alternatively, if no such file exists, the file name is interpreted
 *  as a comma-separated list of samples. When ploidy is not present,
 *  the default ploidy 2 is assumed.
 *
 *  Returns an array of sample names, where the byte value just after \0
 *  indicates the ploidy.
 */
static char **read_samples(const char *fn, int *_n)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t s = {0,0,0};
	int dret, n = 0, max = 0;
	char **sam = 0;
	*_n = 0;
	fp = gzopen(fn, "r");
	if (fp == 0) 
    {
        // interpret as sample names, not as a file name
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
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, &s, &dret) >= 0) {
		int l;
		if (max == n) {
			max = max? max<<1 : 4;
			sam = (char**) realloc(sam, sizeof(char*)*max);
		}
		l = s.l;
		sam[n] = (char*) malloc(sizeof(char)*(s.l+2));
		strcpy(sam[n], s.s);
		sam[n][l+1] = 2; // by default, diploid
		if (dret != '\n') {
			if (ks_getuntil(ks, 0, &s, &dret) >= 0) { // read ploidy, 1 or 2
				int x = (int)s.s[0] - '0'; // Convert ASCII digit to decimal
				if (x == 1 || x == 2) sam[n][l+1] = x;
				else fprintf(stderr, "(%s) ploidy can only be 1 or 2; assume diploid\n", __func__);
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
    args->files = bcf_sr_init();

    // Open files for input and output, initialize structures
    if ( args->targets )
    {
        if ( bcf_sr_set_targets(args->files, args->targets,0)<0 )
            error("Failed to read the targets: %s\n", args->targets);
    }
    if ( args->regions )
    {
        args->files->require_index = 1;
        if ( bcf_sr_set_targets(args->files, args->regions,0)<0 )
            error("Failed to read the targets: %s\n", args->regions);
    }
    
    int mode = FT_UNKN;
    if ( !strcmp(args->bcf_fname,"-") )
    {
        mode |= args->flag & CF_VCFIN ? FT_VCF : FT_BCF;        // VCF or BCF on input?
    }
    if ( !bcf_sr_open_reader(args->files, args->bcf_fname, mode) ) error("Failed to open: %s\n", args->bcf_fname);

    if ( args->nsamples && args->nsamples != args->files->readers[0].header->n[BCF_DT_SAMPLE] )
    {
        args->samples_map = (int *) malloc(sizeof(int)*args->nsamples);
        args->aux.hdr = bcf_hdr_subset(args->files->readers[0].header, args->nsamples, args->samples, args->samples_map);

        // Check if the supplied samples are present in the VCF and print a warning if not.
        if ( args->nsamples != args->aux.hdr->n[BCF_DT_SAMPLE] )
        {
            int i;
            for (i=0; i<args->nsamples; i++)
                if ( args->samples_map[i]<0 ) fprintf(stderr,"Warning: no such sample: \"%s\"\n", args->samples[i]);
        }
    }
    else
        args->aux.hdr = bcf_hdr_dup(args->files->readers[0].header);

    args->out_fh = hts_open("-", hts_bcf_wmode(args->output_type), 0);

    if ( args->flag & CF_QCALL ) 
        return;

    if ( args->flag & CF_MCALL ) 
        mcall_init(&args->aux);

    if ( args->flag & CF_CCALL )
        ccall_init(&args->aux);

    bcf_hdr_append_version(args->aux.hdr, args->argc, args->argv, "bcftools_call");
    bcf_hdr_fmt_text(args->aux.hdr);
    vcf_hdr_write(args->out_fh, args->aux.hdr);
}

static void destroy_data(args_t *args)
{
    if ( args->flag & CF_CCALL ) ccall_destroy(&args->aux);
    else if ( args->flag & CF_MCALL ) mcall_destroy(&args->aux);
    else if ( args->flag & CF_QCALL ) qcall_destroy(&args->aux);
    int i;
    for (i=0; i<args->nsamples; i++) free(args->samples[i]);
    free(args->samples);
    free(args->samples_map);
    free(args->aux.ploidy);
    bcf_hdr_destroy(args->aux.hdr);
    hts_close(args->out_fh);
    bcf_sr_destroy(args->files);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: bcftools call [options] <in.bcf> [reg]\n");
    fprintf(stderr, "File format options:\n");
    fprintf(stderr, "       -o TYPE   output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "       -Q        output the QCALL likelihood format\n");
    fprintf(stderr, "       -r REG    restrict to comma-separated list of regions or regions listed in tab-delimited indexed file\n");
    fprintf(stderr, "       -s STR    comma-separated list or file name with a list of samples to use [all samples]\n");
    fprintf(stderr, "       -l REG    same as -r but streams rather than index-jumps to it. Coordinates are 1-based, inclusive\n");
    fprintf(stderr, "       -V        input is VCF\n");
    fprintf(stderr, "\nInput/output options:\n");
    fprintf(stderr, "       -A        keep all possible alternate alleles at variant sites\n");
    fprintf(stderr, "       -G        output sites only, drop genotype fields\n");
    fprintf(stderr, "       -L        calculate LD for adjacent sites\n");
    fprintf(stderr, "       -N        skip sites where REF is not A/C/G/T\n");
    fprintf(stderr, "       -v        output potential variant sites only\n");
    fprintf(stderr, "\nConsensus/variant calling options:\n");
    fprintf(stderr, "       -c        the original calling method (conflicts with -m)\n");
    fprintf(stderr, "       -d FLOAT  skip loci where less than FLOAT fraction of samples covered [0]\n");
    fprintf(stderr, "       -i FLOAT  indel-to-substitution ratio (-c only) [%.4g]\n", args->aux.indel_frac);
    fprintf(stderr, "       -S STR    skip the given variant type <snps|indels>\n");
    fprintf(stderr, "       -m        alternative model for multiallelic and rare-variant calling (conflicts with -c)\n");
    fprintf(stderr, "       -p FLOAT  variant if P(ref|D)<FLOAT with -c [0.5] or another allele accepted if P(chi^2)>=1-FLOAT with -m [1e-2]\n");
    fprintf(stderr, "       -P STR    type of prior: full, cond2, flat (-c only) [full]\n");
    fprintf(stderr, "       -t FLOAT  scaled substitution mutation rate (-c only) [%.4g]\n", args->aux.theta);
    fprintf(stderr, "       -T STR    constrained calling; STR can be: pair, trioauto, trioxd and trioxs (see manual) [null]\n");
    fprintf(stderr, "\nContrast calling and association test options:\n");
    fprintf(stderr, "       -1 INT    number of group-1 samples [0]\n");
    fprintf(stderr, "       -C FLOAT  posterior constrast for LRT<FLOAT and P(ref|D)<0.5 [%g]\n", args->aux.min_lrt);
    fprintf(stderr, "       -U INT    number of permutations for association testing (effective with -1) [0]\n");
    fprintf(stderr, "       -X FLOAT  only perform permutations for P(chi^2)<FLOAT [%g]\n", args->aux.min_perm_p);
    fprintf(stderr, "\n");
    exit(-1);
}

int main_vcfcall(int argc, char *argv[])
{
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

    // Note: Some of the functionality was lost in the transition but will be put back on demand (pd3 todo)

	while ((c = getopt(argc, argv, "N1:l:cC:AGvVP:t:p:QLi:S:s:U:X:d:T:mr:l:o:")) >= 0) 
    {
		switch (c) 
        {
            case 'N': args.flag |= CF_ACGT_ONLY; break;                 // omit sites where first base in REF is N
            case '1': args.aux.ngrp1_samples = atoi(optarg); break;     // is anyone using this? Please expand docs
            //  case 'l': args.bed = bed_read(optarg);          // todo: with few sites, use index to jump rather than stream; genotypes given alleles
            //            if (!args.bed) 
            //                error("Could not read \"%s\"\n", optarg); 
            //            break;
            case 'G': args.flag |= CF_NO_GENO; break;       // output only sites, no genotype fields
            case 'A': args.aux.flag |= CALL_KEEPALT; break;
            case 'V': args.flag |= CF_VCFIN; break;
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
            case 'S': 
                      if ( !strcasecmp(optarg,"snps") ) args.flag |= CF_INDEL_ONLY;
                      else if ( !strcasecmp(optarg,"indels") ) args.flag |= CF_NO_INDEL;
                      else error("Unknown argument to -I: \"%s\"\n", optarg);
            case 'm': args.flag |= CF_MCALL; break;         // multiallelic calling method
            case 't': args.aux.theta = atof(optarg); break;
            case 'p': p_arg = atof(optarg); break;
            case 'i': args.aux.indel_frac = atof(optarg); break;
            case 'Q': args.flag |= CF_QCALL; break;
            case 'L': args.flag |= CF_ADJLD; break;
            case 'U': args.aux.n_perm = atoi(optarg); break;
            case 'C': args.aux.min_lrt = atof(optarg); break;
            case 'X': args.aux.min_perm_p = atof(optarg); break;
            // case 'd': args.aux.min_smpl_frac = atof(optarg); break; // todo
            case 'r': args.regions = optarg; break;
            case 'l': args.targets = optarg; break;
            case 's': 
                      args.samples = read_samples(optarg, &args.nsamples);
                      args.aux.nsamples = args.nsamples;
                      args.aux.ploidy = (uint8_t*) calloc(args.nsamples+1, 1);
                      for (i=0; i<args.nsamples; i++) args.aux.ploidy[i] = args.samples[i][strlen(args.samples[i]) + 1];
                      break;
            // todo
            // case 'T': 
            //           if (strcmp(optarg, "trioauto") == 0) args.aux.trio = bcf_trio_prep(0, 0);
            //           else if (strcmp(optarg, "trioxd") == 0) args.aux.trio = bcf_trio_prep(1, 0);
            //           else if (strcmp(optarg, "trioxs") == 0) args.aux.trio = bcf_trio_prep(1, 1);
            //           else if (strcmp(optarg, "pair") == 0) args.flag |= CF_PAIRCALL;
            //           else error("[%s] Option '-T' can only take value trioauto, trioxd or trioxs.\n", __func__);
            // case 'P':
            //           if (strcmp(optarg, "full") == 0) args.aux.prior_type = MC_PTYPE_FULL;
            //           else if (strcmp(optarg, "cond2") == 0) args.aux.prior_type = MC_PTYPE_COND2;
            //           else if (strcmp(optarg, "flat") == 0) args.aux.prior_type = MC_PTYPE_FLAT;
            //           else args.aux.prior_file = strdup(optarg);
            //           break;

            //  case 'e': args.flag |= CF_EM; break;        // is anyone using -e without anything else? Disabling for now
            //  case 'M': args.flag |= CF_ANNO_MAX; break;  // undocumented -> not supported for now
            //  case 'Y': args.flag |= CF_QCNT; break;      // undocumented -> not supported for now
            //  case 'K': bcf_p1_fp_lk = gzopen(optarg, "w"); break;    // undocumented -> not supported for now
            default: usage(&args);
        }
    }
	if (argc == optind) usage(&args); 
    args.bcf_fname = argv[optind++];

    // Sanity check options and initialize
    if ( (args.flag & CF_CCALL ? 1 : 0) + (args.flag & CF_MCALL ? 1 : 0) + (args.flag & CF_QCALL ? 1 : 0) > 1 ) error("Only one of -c, -m, or -Q options can be given\n");
    if ( !(args.flag & CF_CCALL) && !(args.flag & CF_MCALL) && !(args.flag & CF_QCALL) ) error("Expected one of -c, -m, or -Q options\n");
	if ( args.aux.n_perm && args.aux.ngrp1_samples<=0 ) error("Expected -1 with -U\n");    // not sure about this, please fix
    if ( p_arg!=-1 ) { args.aux.pref = p_arg; args.aux.min_ma_lrt = 1 - p_arg; }  // only one is actually used

    init_data(&args);

    while ( bcf_sr_next_line(args.files) )
    {
        bcf1_t *bcf_rec = args.files->readers[0].buffer[0];
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

        // todo: subsample, bed overlaps

        // QCall output (todo)
        if ( args.flag & CF_QCALL ) 
        {
            qcall(&args.aux, bcf_rec);
            continue;
        }

        // Output from different calling models
        int ret;
        if ( args.flag & CF_MCALL )
            ret = mcall(&args.aux, bcf_rec);
        else
            ret = ccall(&args.aux, bcf_rec);

        if ( ret==-1 ) error("Something is wrong\n");
        if ( (args.aux.flag & CALL_VARONLY) && ret==0 ) continue;     // not a variant

        if ( args.output_type & FT_BCF ) bcf1_sync(bcf_rec);
        vcf_write1(args.out_fh, args.aux.hdr, bcf_rec);
    }
    destroy_data(&args);
	return 0;
}

