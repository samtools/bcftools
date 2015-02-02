/*  vcfcall.c -- SNP/indel variant calling from VCF/BCF.

    Copyright (C) 2013-2014 Genome Research Ltd.

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

#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <htslib/vcf.h>
#include <time.h>
#include <zlib.h>
#include <stdarg.h>
#include <htslib/kfunc.h>
#include <htslib/synced_bcf_reader.h>
#include <ctype.h>
#include "bcftools.h"
#include "call.h"
#include "prob1.h"

void error(const char *format, ...);

#ifdef _WIN32
#define srand48(x) srand(x)
#define lrand48() rand()
#endif

#define CF_NO_GENO      1
#define CF_INS_MISSED   (1<<1)
#define CF_CCALL        (1<<2)
#define CF_GVCF         (1<<3)
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
    char *bcf_fname, *output_fname;
    char **samples;             // for subsampling and ploidy
    int nsamples, *samples_map;
    char *regions, *targets;    // regions to process
    int regions_is_file, targets_is_file;

    bcf1_t *missed_line;
    call_t aux;     // parameters and temporary data
    gvcf_t gvcf;

    int argc;
    char **argv;

    //  int flag, prior_type, n1, n_sub, *sublist, n_perm;
    //  uint32_t *trio_aux;
    //  char *prior_file, **subsam;
    //  uint8_t *ploidy;
    //  double theta, pref, indel_frac, min_smpl_frac, min_lrt;
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

static char **parse_ped_samples(call_t *call, char **vals, int _n)
{
    int i, j, max = 0, n = 0;
    kstring_t str = {0,0,0};
    char **sam = NULL;
    for (i=0; i<_n; i++)
    {
        str.l = 0;
        kputs(vals[i], &str);
        char *col_ends[5], *tmp = str.s;
        j = 0;
        while ( *tmp && j<5 )
        {
            if ( isspace(*tmp) )
            {
                *tmp = 0;
                ++tmp;
                while ( isspace(*tmp) ) tmp++;  // allow multiple spaces
                col_ends[j] = tmp-1;
                j++;
                continue;
            }
            tmp++;
        }
        if ( j!=5 ) break;

        family_t *fam = get_family(call->fams, call->nfams, str.s);
        if ( !fam )
        {
            call->nfams++;
            hts_expand(family_t, call->nfams, call->mfams, call->fams);
            fam = &call->fams[call->nfams-1];
            fam->name = strdup(str.s);
            for (j=0; j<3; j++) fam->sample[j] = -1;
        }

        int ploidy = 2;
        if ( call->flag & (CALL_CHR_X|CALL_CHR_Y) )
        {
            if ( col_ends[3][1]=='1' ) ploidy = 1; // male: one chrX and one chrY copy
            else
                ploidy = call->flag & CALL_CHR_X ? 2 : 0; // female: two chrX copies, no chrY
        }
        sam = add_sample(sam, &n, &max, col_ends[0]+1, ploidy, &j);
        if ( strcmp(col_ends[1]+1,"0") )    // father
        {
            if ( fam->sample[CHILD]>=0 ) error("Multiple childs in %s [%s,%s]\n", str.s, sam[j],sam[fam->sample[CHILD]]);
            fam->sample[CHILD] = j;
            if ( fam->sample[FATHER]>=0 ) error("Two fathers in %s?\n", str.s);
            sam = add_sample(sam, &n, &max, col_ends[1]+1, call->flag & (CALL_CHR_X|CALL_CHR_Y) ? 1 : 2, &fam->sample[FATHER]);
        }
        if ( strcmp(col_ends[2]+1,"0") )    // mother
        {
            if ( fam->sample[MOTHER]>=0 ) error("Two mothers in %s?\n", str.s);
            sam = add_sample(sam, &n, &max, col_ends[2]+1, call->flag & CALL_CHR_Y ? 0 : 2, &fam->sample[MOTHER]);
        }
    }
    free(str.s);

    if ( i!=_n ) // not a ped file
    {
        if ( i>0 ) error("Could not parse the samples, thought it was PED format, some rows have 5 columns?!\n");
        return NULL;
    }
    assert( n==_n );
    for (i=0; i<call->nfams; i++)
        assert( call->fams[i].sample[0]>=0 && call->fams[i].sample[1]>=0 && call->fams[i].sample[2]>=0 ); // multiple childs, not a trio

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
static char **read_samples(call_t *call, const char *fn, int is_file, int *_n)
{
    int i, n;
    char **vals = hts_readlist(fn, is_file, &n);
    if ( !vals ) error("Could not read the file: %s\n", fn);

    char **smpls = parse_ped_samples(call, vals, n);
    if ( !smpls )
    {
        smpls = (char**) malloc(sizeof(char*)*n);
        for (i=0; i<n; i++)
        {
            char *s = vals[i];
            while ( *s && !isspace(*s) ) s++;
            int len = s-vals[i];
            smpls[i] = (char*) malloc(len+2);
            strncpy(smpls[i],vals[i],len);
            smpls[i][len] = 0;
            while ( *s && isspace(*s) ) s++;
            int x = 2;
            if ( *s )
            {
                x = (int)s[0] - '0'; // Convert ASCII digit to decimal
                if (x != 0 && x != 1 && x != 2) error("Ploidy can only be 0, 1 or 2: %s\n", vals[i]);
            }
            smpls[i][len+1] = x;
        }
    }

    for (i=0; i<n; i++) free(vals[i]);
    free(vals);

    *_n = n;
    return smpls;
}

static void init_missed_line(args_t *args)
{
    int i;
    for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++)
    {
        args->aux.gts[i*2]   = bcf_gt_missing;
        args->aux.gts[i*2+1] = bcf_int32_vector_end;
    }
    args->missed_line = bcf_init1();
    bcf_update_genotypes(args->aux.hdr, args->missed_line, args->aux.gts, 2*bcf_hdr_nsamples(args->aux.hdr));
    bcf_float_set_missing(args->missed_line->qual);
}

static void print_missed_line(bcf_sr_regions_t *regs, void *data)
{
    args_t *args = (args_t*) data;
    call_t *call = &args->aux;
    bcf1_t *missed = args->missed_line;

    if ( args->flag & CF_GVCF ) error("todo: Combine --gvcf and --insert-missed\n");

    char *ss = regs->line.s;
    int i = 0;
    while ( i<args->aux.srs->targets_als-1 && *ss )
    {
        if ( *ss=='\t' ) i++;
        ss++;
    }
    if ( !*ss ) error("Could not parse: [%s] (%d)\n", regs->line.s,args->aux.srs->targets_als);

    missed->rid  = bcf_hdr_name2id(call->hdr,regs->seq_names[regs->prev_seq]);
    missed->pos  = regs->start;
    bcf_update_alleles_str(call->hdr, missed,ss);

    bcf_write1(args->out_fh, call->hdr, missed);
}

static void init_data(args_t *args)
{
    args->aux.srs = bcf_sr_init();

    // Open files for input and output, initialize structures
    if ( args->targets )
    {
        if ( bcf_sr_set_targets(args->aux.srs, args->targets, args->targets_is_file, args->aux.flag&CALL_CONSTR_ALLELES ? 3 : 0)<0 )
            error("Failed to read the targets: %s\n", args->targets);

        if ( args->aux.flag&CALL_CONSTR_ALLELES && args->flag&CF_INS_MISSED )
        {
            args->aux.srs->targets->missed_reg_handler = print_missed_line;
            args->aux.srs->targets->missed_reg_data = args;
        }
    }
    if ( args->regions )
    {
        if ( bcf_sr_set_regions(args->aux.srs, args->regions, args->regions_is_file)<0 )
            error("Failed to read the targets: %s\n", args->regions);
    }

    int i;
    if ( !bcf_sr_add_reader(args->aux.srs, args->bcf_fname) ) error("Failed to open %s: %s\n", args->bcf_fname,bcf_sr_strerror(args->aux.srs->errnum));

    if ( args->nsamples && args->nsamples != bcf_hdr_nsamples(args->aux.srs->readers[0].header) )
    {
        args->samples_map = (int *) malloc(sizeof(int)*args->nsamples);
        args->aux.hdr = bcf_hdr_subset(args->aux.srs->readers[0].header, args->nsamples, args->samples, args->samples_map);
        if ( !args->aux.hdr ) error("Error occurred while subsetting samples\n");
        for (i=0; i<args->nsamples; i++)
            if ( args->samples_map[i]<0 ) error("No such sample: %s\n", args->samples[i]);
        if ( !bcf_hdr_nsamples(args->aux.hdr) ) error("No matching sample found\n");
    }
    else
    {
        args->aux.hdr = bcf_hdr_dup(args->aux.srs->readers[0].header);
        for (i=0; i<args->nsamples; i++)
            if ( bcf_hdr_id2int(args->aux.hdr,BCF_DT_SAMPLE,args->samples[i])<0 )
                error("No such sample: %s\n", args->samples[i]);
    }

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
        uint8_t *ploidy = (uint8_t*) calloc(bcf_hdr_nsamples(args->aux.hdr), 1);
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
        args->nsamples = bcf_hdr_nsamples(args->aux.hdr);
        for (i=0; i<args->nsamples; i++)
            assert( ploidy[i]==0 || ploidy[i]==1 || ploidy[i]==2 );
        free(args->aux.ploidy);
        args->aux.ploidy = ploidy;
    }

    args->out_fh = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));

    if ( args->flag & CF_QCALL )
        return;

    if ( args->flag & CF_MCALL )
        mcall_init(&args->aux);

    if ( args->flag & CF_CCALL )
        ccall_init(&args->aux);

    if ( args->flag&CF_GVCF )
    {
        bcf_hdr_append(args->aux.hdr,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
        args->gvcf.rid  = -1;
        args->gvcf.line = bcf_init1();
        args->gvcf.gt   = (int32_t*) malloc(2*sizeof(int32_t)*bcf_hdr_nsamples(args->aux.hdr));
        for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++)
        {
            args->gvcf.gt[2*i+0] = bcf_gt_unphased(0);
            args->gvcf.gt[2*i+1] = bcf_gt_unphased(0);
        }
    }

    bcf_hdr_remove(args->aux.hdr, BCF_HL_INFO, "QS");
    bcf_hdr_remove(args->aux.hdr, BCF_HL_INFO, "I16");

    bcf_hdr_append_version(args->aux.hdr, args->argc, args->argv, "bcftools_call");
    bcf_hdr_write(args->out_fh, args->aux.hdr);

    if ( args->flag&CF_INS_MISSED ) init_missed_line(args);
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
    if ( args->missed_line ) bcf_destroy(args->missed_line);
    if ( args->gvcf.line ) bcf_destroy(args->gvcf.line);
    free(args->gvcf.gt);
    free(args->gvcf.dp);
    free(args->samples);
    free(args->samples_map);
    free(args->aux.ploidy);
    bcf_hdr_destroy(args->aux.hdr);
    hts_close(args->out_fh);
    bcf_sr_destroy(args->aux.srs);
}

void parse_novel_rate(args_t *args, const char *str)
{
    if ( sscanf(str,"%le,%le,%le",&args->aux.trio_Pm_SNPs,&args->aux.trio_Pm_del,&args->aux.trio_Pm_ins)==3 )  // explicit for all
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_del  = 1 - args->aux.trio_Pm_del;
        args->aux.trio_Pm_ins  = 1 - args->aux.trio_Pm_ins;
    }
    else if ( sscanf(str,"%le,%le",&args->aux.trio_Pm_SNPs,&args->aux.trio_Pm_del)==2 )   // dynamic for indels
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_ins  = -1;    // negative value for dynamic calculation
    }
    else if ( sscanf(str,"%le",&args->aux.trio_Pm_SNPs)==1 )  // same for all
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_del  = -1;
        args->aux.trio_Pm_ins  = -1;
    }
    else error("Could not parse --novel-rate %s\n", str);
}

static int parse_format_flag(const char *str)
{
    int flag = 0;
    const char *ss = str;
    while ( *ss )
    {
        const char *se = ss;
        while ( *se && *se!=',' ) se++;
        if ( !strncasecmp(ss,"GQ",se-ss) ) flag |= CALL_FMT_GQ;
        else if ( !strncasecmp(ss,"GP",se-ss) ) flag |= CALL_FMT_GP;
        else
        {
            fprintf(stderr,"Could not parse \"%s\"\n", str);
            exit(1);
        }
        if ( !*se ) break;
        ss = se + 1;
    }
    return flag;
}


static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   SNP/indel variant calling from VCF/BCF. To be used in conjunction with samtools mpileup.\n");
    fprintf(stderr, "         This command replaces the former \"bcftools view\" caller. Some of the original\n");
    fprintf(stderr, "         functionality has been temporarily lost in the process of transition to htslib,\n");
    fprintf(stderr, "         but will be added back on popular demand. The original calling model can be\n");
    fprintf(stderr, "         invoked with the -c option.\n");
    fprintf(stderr, "Usage:   bcftools call [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "File format options:\n");
    fprintf(stderr, "   -o, --output <file>             write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "   -r, --regions <region>          restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file <file>       restrict to regions listed in a file\n");
    fprintf(stderr, "   -s, --samples <list>            list of samples to include [all samples]\n");
    fprintf(stderr, "   -S, --samples-file <file>       PED file or a file with optional second column for ploidy (0, 1 or 2) [all samples]\n");
    fprintf(stderr, "   -t, --targets <region>          similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "   -T, --targets-file <file>       similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input/output options:\n");
    fprintf(stderr, "   -A, --keep-alts                 keep all possible alternate alleles at variant sites\n");
    fprintf(stderr, "   -f, --format-fields <list>      output format fields: GQ,GP (lowercase allowed) []\n");
    fprintf(stderr, "   -g, --gvcf <minDP>              output gVCF blocks of homozygous REF calls\n");
    fprintf(stderr, "   -i, --insert-missed             output also sites missed by mpileup but present in -T\n");
    fprintf(stderr, "   -M, --keep-masked-ref           keep sites with masked reference allele (REF=N)\n");
    fprintf(stderr, "   -V, --skip-variants <type>      skip indels/snps\n");
    fprintf(stderr, "   -v, --variants-only             output variant sites only\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus/variant calling options:\n");
    fprintf(stderr, "   -c, --consensus-caller          the original calling method (conflicts with -m)\n");
    fprintf(stderr, "   -C, --constrain <str>           one of: alleles, trio (see manual)\n");
    fprintf(stderr, "   -m, --multiallelic-caller       alternative model for multiallelic and rare-variant calling (conflicts with -c)\n");
    fprintf(stderr, "   -n, --novel-rate <float>,[...]  likelihood of novel mutation for constrained trio calling, see man page for details [1e-8,1e-9,1e-9]\n");
    fprintf(stderr, "   -p, --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5]\n");
    fprintf(stderr, "   -P, --prior <float>             mutation rate (use bigger for greater sensitivity) [1.1e-3]\n");
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
    args.aux.theta      = 1.1e-3;
    args.aux.pref       = 0.5;
    args.aux.min_perm_p = 0.01;
    args.aux.min_lrt    = 1;
    args.flag           = CF_ACGT_ONLY;
    args.output_fname   = "-";
    args.output_type    = FT_VCF;
    args.aux.trio_Pm_SNPs = 1 - 1e-8;
    args.aux.trio_Pm_ins  = args.aux.trio_Pm_del  = 1 - 1e-9;

    int i, c, samples_is_file = 0;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"gvcf",1,0,'g'},
        {"format-fields",1,0,'f'},
        {"output",1,0,'o'},
        {"output-type",1,0,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"keep-alts",0,0,'A'},
        {"insert-missed",0,0,'i'},
        {"skip-Ns",0,0,'N'},            // now the new default
        {"keep-masked-refs",0,0,'M'},
        {"skip-variants",1,0,'V'},
        {"variants-only",0,0,'v'},
        {"consensus-caller",0,0,'c'},
        {"constrain",1,0,'C'},
        {"multiallelic-caller",0,0,'m'},
        {"pval-threshold",1,0,'p'},
        {"prior",1,0,'P'},
        {"chromosome-X",0,0,'X'},
        {"chromosome-Y",0,0,'Y'},
        {"novel-rate",1,0,'n'},
        {0,0,0,0}
    };

    char *tmp = NULL;
    while ((c = getopt_long(argc, argv, "h?o:O:r:R:s:S:t:T:ANMV:vcmp:C:XYn:P:f:ig:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'g':
                args.flag |= CF_GVCF;
                args.gvcf.min_dp = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse, expected integer argument: -g %s\n", optarg);
                break;
            case 'f': args.aux.output_tags |= parse_format_flag(optarg); break;
            case 'M': args.flag &= ~CF_ACGT_ONLY; break;     // keep sites where REF is N
            case 'N': args.flag |= CF_ACGT_ONLY; break;      // omit sites where first base in REF is N (the new default)
            case 'A': args.aux.flag |= CALL_KEEPALT; break;
            case 'c': args.flag |= CF_CCALL; break;          // the original EM based calling method
            case 'i': args.flag |= CF_INS_MISSED; break;
            case 'v': args.aux.flag |= CALL_VARONLY; break;
            case 'o': args.output_fname = optarg; break;
            case 'O':
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
            case 'V':
                      if ( !strcasecmp(optarg,"snps") ) args.flag |= CF_INDEL_ONLY;
                      else if ( !strcasecmp(optarg,"indels") ) args.flag |= CF_NO_INDEL;
                      else error("Unknown skip category \"%s\" (-S argument must be \"snps\" or \"indels\")\n", optarg);
                      break;
            case 'm': args.flag |= CF_MCALL; break;         // multiallelic calling method
            case 'p':
                args.aux.pref = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: --pval-threshold %s\n", optarg);
                break;
            case 'P': args.aux.theta = strtod(optarg,&tmp);
                      if ( *tmp ) error("Could not parse, expected float argument: -P %s\n", optarg);
                      break;
            case 'n': parse_novel_rate(&args,optarg); break;
            case 'r': args.regions = optarg; break;
            case 'R': args.regions = optarg; args.regions_is_file = 1; break;
            case 't': args.targets = optarg; break;
            case 'T': args.targets = optarg; args.targets_is_file = 1; break;
            case 's': samples_fname = optarg; break;
            case 'S': samples_fname = optarg; samples_is_file = 1; break;
            default: usage(&args);
        }
    }
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args.bcf_fname = "-";  // reading from stdin
        else usage(&args);
    }
    else args.bcf_fname = argv[optind++];

    // Sanity check options and initialize
    if ( samples_fname )
    {
        args.samples = read_samples(&args.aux, samples_fname, samples_is_file, &args.nsamples);
        args.aux.ploidy = (uint8_t*) calloc(args.nsamples+1, 1);
        args.aux.all_diploid = 1;
        for (i=0; i<args.nsamples; i++)
        {
            args.aux.ploidy[i] = args.samples[i][strlen(args.samples[i]) + 1];
            if ( args.aux.ploidy[i]!=2 ) args.aux.all_diploid = 0;
        }
    }
    if ( args.flag & CF_GVCF )
    {
        // Force some flags to avoid unnecessary branching
        args.aux.flag &= ~CALL_KEEPALT;
        args.aux.flag |= CALL_VARONLY;
    }
    if ( (args.flag & CF_CCALL ? 1 : 0) + (args.flag & CF_MCALL ? 1 : 0) + (args.flag & CF_QCALL ? 1 : 0) > 1 ) error("Only one of -c or -m options can be given\n");
    if ( !(args.flag & CF_CCALL) && !(args.flag & CF_MCALL) && !(args.flag & CF_QCALL) ) error("Expected -c or -m option\n");
    if ( args.aux.n_perm && args.aux.ngrp1_samples<=0 ) error("Expected -1 with -U\n");    // not sure about this, please fix
    if ( args.aux.flag & CALL_CONSTR_ALLELES )
    {
        if ( !args.targets ) error("Expected -t or -T with \"-C alleles\"\n");
        if ( !(args.flag & CF_MCALL) ) error("The \"-C alleles\" mode requires -m\n");
    }
    if ( args.aux.flag & CALL_CHR_X && args.aux.flag & CALL_CHR_Y ) error("Only one of -X or -Y should be given\n");
    if ( args.flag & CF_INS_MISSED && !(args.aux.flag&CALL_CONSTR_ALLELES) ) error("The -i option requires -C alleles\n");
    init_data(&args);

    while ( bcf_sr_next_line(args.aux.srs) )
    {
        bcf1_t *bcf_rec = args.aux.srs->readers[0].buffer[0];
        if ( args.samples_map ) bcf_subset(args.aux.hdr, bcf_rec, args.nsamples, args.samples_map);
        bcf_unpack(bcf_rec, BCF_UN_STR);

        // Skip unwanted sites
        if ( args.aux.flag & CALL_VARONLY )
        {
            int is_ref = 0;
            if ( bcf_rec->n_allele==1 ) is_ref = 1;     // not a variant
            else if ( bcf_rec->n_allele==2 )
            {
                // second allele is mpileup's X, not a variant
                if ( bcf_rec->d.allele[1][0]=='X' ) is_ref = 1;
                else if ( bcf_rec->d.allele[1][0]=='<' && bcf_rec->d.allele[1][1]=='X' && bcf_rec->d.allele[1][2]=='>' ) is_ref = 1;
                else if ( bcf_rec->d.allele[1][0]=='<' && bcf_rec->d.allele[1][1]=='*' && bcf_rec->d.allele[1][2]=='>' ) is_ref = 1;
            }
            if ( is_ref )
            {
                // gVCF output
                if ( args.flag & CF_GVCF ) gvcf_write(args.out_fh, &args.gvcf, args.aux.hdr, bcf_rec, 1);
                continue;
            }
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

        // gVCF output
        if ( args.flag & CF_GVCF )
        {
            gvcf_write(args.out_fh, &args.gvcf, args.aux.hdr, bcf_rec, ret?0:1);
            continue;
        }

        // Normal output
        if ( (args.aux.flag & CALL_VARONLY) && ret==0 ) continue;     // not a variant
        bcf_write1(args.out_fh, args.aux.hdr, bcf_rec);
    }
    if ( args.flag & CF_GVCF ) gvcf_write(args.out_fh, &args.gvcf, args.aux.hdr, NULL, 0);
    if ( args.flag & CF_INS_MISSED ) bcf_sr_regions_flush(args.aux.srs->targets);
    destroy_data(&args);
    return 0;
}

