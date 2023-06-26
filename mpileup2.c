/*  mpileup2.c -- mpileup2, new version of mpileup

    Copyright (C) 2022-2023 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include <inttypes.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>
#include <htslib/hts_os.h>
#include <assert.h>
#include "bcftools.h"
#include "mpileup2/mpileup.h"
#include "bam2bcf.h"
#include "gvcf.h"

#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)
#define MPLP_REALN_PARTIAL  (1<<13)

typedef struct
{
    int argc;
    char **argv, *output_fname;
    int flag;
    int output_type, clevel, max_read_len;
    int record_cmd_line;
    mpileup_t *mplp;
    bcf_callaux_t bca;
    gvcf_t *gvcf;
}
args_t;

static int run_mpileup(args_t *args)
{
    int ret;
    while ( (ret=mpileup_next(args->mplp))==1 )
    {
    }
    return ret;
}

#define SET_FMT_FLAG(str,bit,msg) \
    if (!strcasecmp(tag,str) || !strcasecmp(tag,"FMT/"str) || !strcasecmp(tag,"FORMAT/"str)) \
    { \
        if ( *msg ) fprintf(stderr,"%s",msg); \
        if ( exclude ) \
            *flag &= ~bit; \
        else \
            *flag |= bit; \
        free(tags[i]); \
        continue; \
    }
#define SET_INFO_FLAG(str,bit,msg) if (!strcasecmp(tag,"INFO/"str)) \
    { \
        if ( exclude ) \
            *flag &= ~bit; \
        else \
            *flag |= bit; \
        free(tags[i]); \
        continue; \
    }

static void parse_format_flag(uint32_t *flag, const char *str)
{
    int i, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags);
    for(i=0; i<n_tags; i++)
    {
        int exclude = tags[i][0]=='-' ? 1 : 0;
        char *tag = exclude ? tags[i]+1 : tags[i];
        SET_FMT_FLAG("AD", B2B_FMT_AD, "");
        SET_FMT_FLAG("ADF", B2B_FMT_ADF, "");
        SET_FMT_FLAG("ADR", B2B_FMT_ADR, "");
        SET_FMT_FLAG("DP", B2B_FMT_DP, "");
        SET_FMT_FLAG("DP4", B2B_FMT_DP4, "[warning] tag DP4 functional, but deprecated. Please switch to `ADF` and `ADR` in future.\n");
        SET_FMT_FLAG("DPR", B2B_FMT_DPR, "[warning] tag DPR functional, but deprecated. Please switch to `AD` in future.\n");
        SET_FMT_FLAG("DV", B2B_FMT_DV, "[warning] tag DV functional, but deprecated. Please switch to `AD` in future.\n");
        SET_FMT_FLAG("NMBZ", B2B_FMT_NMBZ, "");
        SET_FMT_FLAG("QS", B2B_FMT_QS, "");
        SET_FMT_FLAG("SP", B2B_FMT_SP, "");
        SET_FMT_FLAG("SCR", B2B_FMT_SCR, "");
        SET_INFO_FLAG("DPR", B2B_INFO_DPR, "[warning] tag INFO/DPR functional, but deprecated. Please switch to `INFO/AD` in future.\n");
        SET_INFO_FLAG("AD", B2B_INFO_AD, "");
        SET_INFO_FLAG("ADF", B2B_INFO_ADF, "");
        SET_INFO_FLAG("ADR", B2B_INFO_ADR, "");
        SET_INFO_FLAG("BQBZ", B2B_INFO_BQBZ, "");
        SET_INFO_FLAG("FS", B2B_INFO_FS, "");
        SET_INFO_FLAG("IDV", B2B_INFO_IDV, "");
        SET_INFO_FLAG("IMF", B2B_INFO_IMF, "");
        SET_INFO_FLAG("MIN_PL_SUM", B2B_INFO_MIN_PL_SUM, "");
        SET_INFO_FLAG("MQ0F", B2B_INFO_MQ0F, "");
        SET_INFO_FLAG("MQBZ", B2B_INFO_MQBZ, "");
        SET_INFO_FLAG("NM", B2B_INFO_NM, "");
        SET_INFO_FLAG("NMBZ", B2B_INFO_NMBZ, "");
        SET_INFO_FLAG("RPBZ", B2B_INFO_RPBZ, "");
        SET_INFO_FLAG("SCBZ", B2B_INFO_SCBZ, "");
        SET_INFO_FLAG("SCR", B2B_INFO_SCR, "");
        SET_INFO_FLAG("SGB", B2B_INFO_SGB, "");
        SET_INFO_FLAG("VDB", B2B_INFO_VDB, "");
        fprintf(stderr,"Could not parse tag \"%s\" in \"%s\"\n", tag, str);
        exit(EXIT_FAILURE);
    }
    if (n_tags) free(tags);
}

// todo: make it possible to turn off some annotations or change the defaults,
//      specifically RPB, VDB, MWU, SGB tests. It would be good to do some
//      benchmarking first to see if it's worth it.
static void list_annotations(FILE *fp)
{
    fprintf(fp,
        "Annotations added by default are in this list prefixed with \"*\". To suppress their output, run with\n"
        "e.g. \"-a -FORMAT/AD\".\n"
        "\n"
        "FORMAT annotation tags available (\"FORMAT/\" prefix is optional):\n"
        "\n"
        "  FORMAT/AD   .. Allelic depth (Number=R,Type=Integer)\n"
        "  FORMAT/ADF  .. Allelic depths on the forward strand (Number=R,Type=Integer)\n"
        "  FORMAT/ADR  .. Allelic depths on the reverse strand (Number=R,Type=Integer)\n"
        "  FORMAT/DP   .. Number of high-quality bases (Number=1,Type=Integer)\n"
        "  FORMAT/NMBZ .. Mann-Whitney U-z test of Number of Mismatches within supporting reads (Number=1,Type=Float)\n"
        "  FORMAT/QS   .. Allele phred-score quality sum for use with `call -mG` and +trio-dnm (Number=R,Type=Integer)\n"
        "  FORMAT/SP   .. Phred-scaled strand bias P-value (Number=1,Type=Integer)\n"
        "  FORMAT/SCR  .. Number of soft-clipped reads (Number=1,Type=Integer)\n"
        "\n"
        "INFO annotation tags available:\n"
        "\n"
        "  INFO/AD    .. Total allelic depth (Number=R,Type=Integer)\n"
        "  INFO/ADF   .. Total allelic depths on the forward strand (Number=R,Type=Integer)\n"
        "  INFO/ADR   .. Total allelic depths on the reverse strand (Number=R,Type=Integer)\n"
        "* INFO/BQBZ  .. Mann-Whitney U test of Base Quality Bias (Number=1,Type=Float)\n"
        "* INFO/FS    .. Phred-scaled p-value using Fisher's exact test to detect strand bias (Number=1,Type=Float)\n"
        "* INFO/IDV   .. Maximum number of raw reads supporting an indel (Number=1,Type=Integer)\n"
        "* INFO/IMF   .. Maximum fraction of raw reads supporting an indel (Number=1,Type=Float)\n"
        "  INFO/MIN_PL_SUM\n"
        "             .. Sum of min PL across all samples before normalization, experimental (Number=1,Type=Integer)\n"
        "* INFO/MQ0F  .. Fraction of reads with zero mapping quality (Number=1,Type=Float)\n"
        "* INFO/MQBZ  .. Mann-Whitney U test of Mapping Quality Bias (Number=1,Type=Float)\n"
        "* INFO/MQSBZ .. Mann-Whitney U-z test of Mapping Quality vs Strand Bias (Number=1,Type=Float)\n"
        "  INFO/NM    .. Approximate average number of mismatches in ref and alt reads, experimental (Number=2,Type=Float)\n"
        "* INFO/NMBZ  .. Mann-Whitney U-z test of Number of Mismatches within supporting reads (Number=1,Type=Float)\n"
        "* INFO/RPBZ  .. Mann-Whitney U test of Read Position Bias (Number=1,Type=Float)\n"
        "* INFO/SCBZ  .. Mann-Whitney U-z test of Soft-Clip Length Bias (Number=1,Type=Float)\n"
        "  INFO/SCR   .. Number of soft-clipped reads (Number=1,Type=Integer)\n"
        "* INFO/SGB   .. Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf (Number=1,Type=Float)\n"
        "* INFO/VDB   .. Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (Number=1,Type=Float)\n"
        "\n");
}

static void print_usage(args_t *args, FILE *fp)
{
    char *tmp_skip_all_set = bam_flag2str(mpileup_get_opt(args->mplp, int, SKIP_ALL_SET));
    char *tmp_skip_any_unset = bam_flag2str(mpileup_get_opt(args->mplp, int, SKIP_ANY_UNSET));
    char *tmp_skip_all_unset = bam_flag2str(mpileup_get_opt(args->mplp, int, SKIP_ALL_UNSET));
    char *tmp_skip_any_set = bam_flag2str(mpileup_get_opt(args->mplp, int, SKIP_ANY_SET));

    fprintf(fp,
        "\n"
        "Usage: bcftools mpileup2 [OPTIONS] in1.bam [in2.bam [...]]\n"
        "\n"
        "Input options:\n"
        "  -6, --illumina1.3+      Quality is in the Illumina-1.3+ encoding\n"
        "  -A, --count-orphans     Do not discard anomalous read pairs\n"
        "  -b, --bam-list FILE     List of input BAM filenames, one per line\n"
        "  -B, --no-BAQ            Disable BAQ (per-Base Alignment Quality)\n"
        "  -C, --adjust-MQ INT     Adjust mapping quality [0]\n"
        "  -D, --full-BAQ          Apply BAQ everywhere, not just in problematic regions\n"
        "  -d, --max-depth INT     Max raw per-sample depth; avoids excessive memory usage [%d]\n", mpileup_get_opt(args->mplp, int, MAX_DP_PER_SAMPLE));
    fprintf(fp,
        "  -E, --redo-BAQ          Recalculate BAQ on the fly, ignore existing BQs\n"
        "  -f, --fasta-ref FILE    Faidx indexed reference sequence file\n"
        "      --no-reference      Do not require fasta reference file\n"
        "  -G, --read-groups FILE  Select or exclude read groups listed in the file\n"
        "  -q, --min-MQ INT        Skip alignments with mapQ smaller than INT [%d]\n", mpileup_get_opt(args->mplp, int, MIN_MQ));
    fprintf(fp,
        "  -Q, --min-BQ INT        Skip bases with baseQ/BAQ smaller than INT [%d]\n", mpileup_get_opt(args->mplp, int, MIN_BQ));
    fprintf(fp,
        "      --max-BQ INT        Limit baseQ/BAQ to no more than INT [%d]\n", mpileup_get_opt(args->mplp, int, MAX_BQ));
    fprintf(fp,
        "      --delta-BQ INT      Use neighbour_qual + INT if less than qual [%d]\n", mpileup_get_opt(args->mplp, int, DELTA_BQ));
    fprintf(fp,
        "  -r, --regions REG[,...] Comma separated list of regions in which pileup is generated\n"
        "  -R, --regions-file FILE Restrict to regions listed in a file\n"
        "      --ignore-RG         Ignore RG tags (one BAM = one sample)\n"
        "  --ls, --skip-all-set STR|INT  Skip reads with all of the bits set []\n");
    fprintf(fp,
        "  --ns, --skip-any-set STR|INT  Skip reads with any of the bits set [%s]\n", tmp_skip_any_set);
    fprintf(fp,
        "  --lu, --skip-all-unset STR|INT  Skip reads with all of the bits unset []\n"
        "  --nu, --skip-any-unset STR|INT  Skip reads with any of the bits unset []\n");
    fprintf(fp,
        "  -s, --samples LIST      Comma separated list of samples to include\n"
        "  -S, --samples-file FILE File of samples to include\n"
        "  -t, --targets REG[,...] Similar to -r but streams rather than index-jumps\n"
        "  -T, --targets-file FILE Similar to -R but streams rather than index-jumps\n"
        "  -x, --ignore-overlaps   Disable read-pair overlap detection\n"
        "      --seed INT          Random number seed used for sampling deep regions [0]\n"
        "\n"
        "Output options:\n"
        "  -a, --annotate LIST     Optional tags to output; '\\?' to list available tags []\n"
        "  -g, --gvcf INT[,...]    Group non-variant sites into gVCF blocks according\n"
        "                          To minimum per-sample DP\n"
        "      --no-version        Do not append version and command line to the header\n"
        "  -o, --output FILE       Write output to FILE [standard output]\n"
        "  -O, --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;\n"
        "                          'z' compressed VCF; 'v' uncompressed VCF; 0-9 compression level [v]\n"
        "      --threads INT       Use multithreading with INT worker threads [0]\n"
        "\n"
        "SNP/INDEL genotype likelihoods options:\n"
        "  -X, --config STR        Specify platform specific profiles (see below)\n"
        "  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [%d]\n", args->bca.extQ);
    fprintf(fp,
        "  -F, --gap-frac FLOAT    Minimum fraction of gapped reads [%f]\n", mpileup_get_opt(args->mplp, float, MIN_REALN_FRAC));
    fprintf(fp,
        "  -h, --tandem-qual INT   Coefficient for homopolymer errors [%d]\n", args->bca.tandemQ);
    fprintf(fp,
        "  -I, --skip-indels       Do not perform indel calling\n"
        "  -L, --max-idepth INT    Subsample to maximum per-sample depth for INDEL calling [%d]\n", mpileup_get_opt(args->mplp, int, MAX_REALN_DP));
    fprintf(fp,
        "  -m, --min-ireads INT    Minimum number gapped reads for indel candidates [%d]\n", mpileup_get_opt(args->mplp, int, MIN_REALN_DP));
    fprintf(fp,
        "  -M, --max-read-len INT  Maximum length of read to pass to BAQ algorithm [%d]\n", mpileup_get_opt(args->mplp, int, MAX_REALN_LEN));
    fprintf(fp,
        "  -o, --open-prob INT     Phred-scaled gap open seq error probability [%d]\n", args->bca.openQ);
    fprintf(fp,
        "  -p, --per-sample-mF     Apply -m and -F per-sample for increased sensitivity\n"
        "  -P, --platforms STR     Comma separated list of platforms for indels [all]\n"
        "  --ar, --ambig-reads STR   What to do with ambiguous indel reads: drop,incAD,incAD0 [drop]\n");
    fprintf(fp,
        "      --indel-bias FLOAT  Raise to favour recall over precision [%.2f]\n", args->bca.indel_bias);
    fprintf(fp,
        "      --indel-size INT    Approximate maximum indel size considered [%d]\n", args->bca.indel_win_size);
    fprintf(fp,"\n");
    fprintf(fp,
        "Configuration profiles activated with -X, --config:\n"
        "    1.12:        -Q13 -h100 -m1 -F0.002\n"
        "    illumina:    [ default values ]\n"
        "    ont:         -B -Q5 --max-BQ 30 -I [also try eg |bcftools call -P0.01]\n"
        "    pacbio-ccs:  -D -Q5 --max-BQ 50 -F0.1 -o25 -e1 --delta-BQ 10 -M99999\n"
        "\n"
        "Notes: Assuming diploid individuals.\n"
        "\n"
        "Example:\n"
        "   # See also http://samtools.github.io/bcftools/howtos/variant-calling.html\n"
        "   bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf\n"
        "\n");

    free(tmp_skip_all_set);
    free(tmp_skip_any_unset);
    free(tmp_skip_all_unset);
    free(tmp_skip_any_set);
}

static void clean_up(args_t *args)
{
    mpileup_destroy(args->mplp);
    free(args);
}

int main_mpileup2(int argc, char *argv[])
{
    int c, ret = 0;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->flag   = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_REALN_PARTIAL | MPLP_SMART_OVERLAPS;
    args->bca.openQ   = 40;
    args->bca.extQ    = 20;
    args->bca.tandemQ = 500;

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

    hts_srand48(0);

#if 0
    const char *file_list = NULL;
    char **fn = NULL;
    int nfiles = 0, use_orphan = 0, noref = 0;
    mplp_conf_t mplp;
    memset(&mplp, 0, sizeof(mplp_conf_t));
    mplp.max_baseQ = 60;
    mplp.delta_baseQ = 30;
    mplp.capQ_thres = 0;
    mplp.max_depth = 250; mplp.max_indel_depth = 250;
    mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 500;
    mplp.min_frac = 0.05; mplp.indel_bias = 1.0; mplp.min_support = 2;
    mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_REALN_PARTIAL
              | MPLP_SMART_OVERLAPS;
    mplp.argc = argc; mplp.argv = argv;
    mplp.rflag_skip_any_set = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    mplp.output_fname = NULL;
    mplp.output_type = FT_VCF;
    mplp.record_cmd_line = 1;
    mplp.n_threads = 0;
    // the default to be changed in future, see also parse_format_flag()
    //mplp.fmt_flag = B2B_INFO_VDB|B2B_INFO_RPB|B2B_INFO_SCB|B2B_INFO_ZSCORE;
    //mplp.ambig_reads = B2B_DROP;
    mplp.max_read_len = 500;
    mplp.indel_win_size = 110;
    mplp.clevel = -1;
#endif

    static const struct option lopts[] =
    {
        {"nu", required_argument, NULL, 16},
        {"lu", required_argument, NULL, 17},
        {"rf", required_argument, NULL, 17},   // old --rf, --incl-flags = --lu, --skip-all-unset
        {"ns", required_argument, NULL, 18},
        {"ff", required_argument, NULL, 18},   // old --ff, --excl-flags = --ns, --skip-any-set
        {"ls", required_argument, NULL, 19},
        {"skip-any-unset", required_argument, NULL, 16},
        {"skip-all-unset", required_argument, NULL, 17},
        {"skip-any-set", required_argument, NULL, 18},
        {"skip-all-set", required_argument, NULL, 19},
        {"output", required_argument, NULL, 3},
        {"open-prob", required_argument, NULL, 4},
        {"ignore-RG", no_argument, NULL, 5},
        {"ignore-rg", no_argument, NULL, 5},
        {"gvcf", required_argument, NULL, 'g'},
        {"no-reference", no_argument, NULL, 7},
        {"no-version", no_argument, NULL, 8},
        {"threads",required_argument,NULL,9},
        {"illumina1.3+", no_argument, NULL, '6'},
        {"count-orphans", no_argument, NULL, 'A'},
        {"bam-list", required_argument, NULL, 'b'},
        {"no-BAQ", no_argument, NULL, 'B'},
        {"no-baq", no_argument, NULL, 'B'},
        {"full-BAQ", no_argument, NULL, 'D'},
        {"full-baq", no_argument, NULL, 'D'},
        {"adjust-MQ", required_argument, NULL, 'C'},
        {"adjust-mq", required_argument, NULL, 'C'},
        {"max-depth", required_argument, NULL, 'd'},
        {"redo-BAQ", no_argument, NULL, 'E'},
        {"redo-baq", no_argument, NULL, 'E'},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"read-groups", required_argument, NULL, 'G'},
        {"region", required_argument, NULL, 'r'},
        {"regions", required_argument, NULL, 'r'},
        {"regions-file", required_argument, NULL, 'R'},
        {"targets", required_argument, NULL, 't'},
        {"targets-file", required_argument, NULL, 'T'},
        {"min-MQ", required_argument, NULL, 'q'},
        {"min-mq", required_argument, NULL, 'q'},
        {"min-BQ", required_argument, NULL, 'Q'},
        {"min-bq", required_argument, NULL, 'Q'},
        {"max-bq", required_argument, NULL, 11},
        {"max-BQ", required_argument, NULL, 11},
        {"delta-BQ", required_argument, NULL, 12},
        {"ignore-overlaps", no_argument, NULL, 'x'},
        {"output-type", required_argument, NULL, 'O'},
        {"samples", required_argument, NULL, 's'},
        {"samples-file", required_argument, NULL, 'S'},
        {"annotate", required_argument, NULL, 'a'},
        {"ext-prob", required_argument, NULL, 'e'},
        {"gap-frac", required_argument, NULL, 'F'},
        {"indel-bias", required_argument, NULL, 10},
        {"indel-size", required_argument, NULL, 15},
        {"tandem-qual", required_argument, NULL, 'h'},
        {"skip-indels", no_argument, NULL, 'I'},
        {"max-idepth", required_argument, NULL, 'L'},
        {"min-ireads", required_argument, NULL, 'm'},
        {"per-sample-mF", no_argument, NULL, 'p'},
        {"per-sample-mf", no_argument, NULL, 'p'},
        {"platforms", required_argument, NULL, 'P'},
        {"max-read-len", required_argument, NULL, 'M'},
        {"config", required_argument, NULL, 'X'},
        {"mwu-u", no_argument, NULL, 'U'},
        {"seed", required_argument, NULL, 13},
        {"ambig-reads", required_argument, NULL, 14},
        {"ar", required_argument, NULL, 14},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, "Ag:f:r:R:q:Q:C:BDd:L:b:P:po:e:h:Im:F:EG:6O:xa:s:S:t:T:M:X:U",lopts,NULL)) >= 0) {
        switch (c) {
        case 'x': mpileup_set(args->mplp, SMART_OVERLAPS, 0); break;
        case  16 :
        {
            int flag = bam_str2flag(optarg);
            if ( flag<0 ) { fprintf(stderr,"Could not parse --nu %s\n", optarg); return 1; }
            mpileup_set(args->mplp, SKIP_ANY_UNSET, flag);
            break;
        }
        case  17 :
        {
            int flag = bam_str2flag(optarg);
            if ( flag<0 ) { fprintf(stderr,"Could not parse --lu %s\n", optarg); return 1; }
            mpileup_set(args->mplp, SKIP_ALL_UNSET, flag);
            break;
        }
        case  18 :
        {
            int flag = bam_str2flag(optarg);
            if ( flag<0 ) { fprintf(stderr,"Could not parse --ns %s\n", optarg); return 1; }
            mpileup_set(args->mplp, SKIP_ANY_SET, flag);
            break;
        }
        case  19 :
        {
            int flag = bam_str2flag(optarg);
            if ( flag<0 ) { fprintf(stderr,"Could not parse --ls %s\n", optarg); return 1; }
            mpileup_set(args->mplp, SKIP_ALL_SET, flag);
            break;
        }
        case  3 : args->output_fname = optarg; break;
        case  4 : args->bca.openQ = atoi(optarg); break;
        case  5 : mpileup_set(args->mplp, SMPL_IGNORE_RG, 1); break;
        case 'g':
            args->gvcf = gvcf_init(optarg);
            if ( !args->gvcf ) error("Could not parse: --gvcf %s\n", optarg);
            break;
        case 'f': mpileup_set(args->mplp, FASTA_REF, optarg); break;
        case  7 : error("todo: --no-reference"); /*noref = 1;*/ break;
        case  8 : args->record_cmd_line = 0; break;
        case  9 : error("todo: --threads\n"); /* note this was used *only* for output VCF compression, practically useless, moreover -Ou is recommended... */ break;
        case 'd': mpileup_set(args->mplp, MAX_DP_PER_SAMPLE, atoi(optarg)); break;
        case 'r': mpileup_set(args->mplp, REGIONS, optarg); break;
        case 'R': mpileup_set(args->mplp, REGIONS_FNAME, optarg); break;
        case 't': mpileup_set(args->mplp, TARGETS, optarg); break;
        case 'T': mpileup_set(args->mplp, TARGETS_FNAME, optarg); break;
        case 'P': error("todo: --platforms\n"); /* this was never used: mplp.pl_list = strdup(optarg); */ break;
        case 'p': args->flag |= MPLP_PER_SAMPLE; break;
        case 'B': args->flag &= ~MPLP_REALN; break;
        case 'D': args->flag &= ~MPLP_REALN_PARTIAL; break;
        case 'I': args->flag |= MPLP_NO_INDEL; break;
        case 'E': args->flag |= MPLP_REDO_BAQ; break;
        case '6': args->flag |= MPLP_ILLUMINA13; break;
        case 's': mpileup_set(args->mplp, SAMPLES, optarg); break;
        case 'S': mpileup_set(args->mplp, SAMPLES_FNAME, optarg); break;
        case 'O':
            switch (optarg[0]) {
                case 'b': args->output_type = FT_BCF_GZ; break;
                case 'u': args->output_type = FT_BCF; break;
                case 'z': args->output_type = FT_VCF_GZ; break;
                case 'v': args->output_type = FT_VCF; break;
                default:
                {
                    char *tmp;
                    args->clevel = strtol(optarg,&tmp,10);
                    if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                }
            }
            if ( optarg[1] )
            {
                char *tmp;
                args->clevel = strtol(optarg+1,&tmp,10);
                if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --output-type %s\n", optarg+1);
            }
            break;
        case 'C': mpileup_set(args->mplp, ADJUST_MQ, atoi(optarg)); break;
        case 'q': mpileup_set(args->mplp, MIN_MQ, atoi(optarg)); break;
        case 'Q': mpileup_set(args->mplp, MIN_BQ, atoi(optarg)); break;
        case  11: mpileup_set(args->mplp, MAX_BQ, atoi(optarg)); break;
        case  12: mpileup_set(args->mplp, DELTA_BQ, atoi(optarg)); break;
        case 'b': if ( mpileup_set(args->mplp, BAM_FNAME, optarg)!=0 ) clean_up(args); break;
        case 'o':
            /* this option used to be -o INT and -o FILE. In the new code disable -o as an alias of --open-prob as its use is very rare */
            args->output_fname = optarg;
            break;
        case 'e': args->bca.extQ = atoi(optarg); break;
        case 'h': args->bca.tandemQ = atoi(optarg); break;
        case 10: // --indel-bias (inverted so higher => more indels called)
            if (atof(optarg) < 1e-2)
                args->bca.indel_bias = 1/1e2;
            else
                args->bca.indel_bias = 1/atof(optarg);
            break;
        case  15: {
                char *tmp;
                args->bca.indel_win_size = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --indel-size %s\n", optarg);
                if ( args->bca.indel_win_size < 110 )
                {
                    args->bca.indel_win_size = 110;
                    fprintf(stderr,"Warning: running with --indel-size %d, the requested value is too small\n",args->bca.indel_win_size);
                }
            }
            break;
        case 'A': args->flag &= ~MPLP_NO_ORPHAN; break;
        case 'F': mpileup_set(args->mplp, MIN_REALN_FRAC, atof(optarg)); break;
        case 'm': mpileup_set(args->mplp, MIN_REALN_DP, atoi(optarg)); break;
        case 'L': mpileup_set(args->mplp, MAX_REALN_DP, atoi(optarg)); break;
        case 'G': mpileup_set(args->mplp, READ_GROUPS_FNAME, optarg); break;
        case 'a':
            if (optarg[0]=='?') {
                list_annotations(stderr);
                return 1;
            }
            parse_format_flag(&args->bca.fmt_flag,optarg);
        break;
        case 'M': args->max_read_len = atoi(optarg); break;
        case 'X':
            if (strcasecmp(optarg, "pacbio-ccs") == 0) {
                mpileup_set(args->mplp, MIN_REALN_FRAC, 0.1);
                mpileup_set(args->mplp, MIN_BQ, 5);
                mpileup_set(args->mplp, MAX_BQ, 50);
                mpileup_set(args->mplp, DELTA_BQ, 10);
                args->bca.openQ = 25;
                args->bca.extQ = 1;
                args->flag |= MPLP_REALN_PARTIAL;
                mpileup_set(args->mplp, MAX_REALN_LEN, 99999);
            } else if (strcasecmp(optarg, "ont") == 0) {
                fprintf(stderr, "For ONT it may be beneficial to also run bcftools call with "
                        "a higher -P, eg -P0.01 or -P 0.1\n");
                mpileup_set(args->mplp, MIN_BQ, 5);
                mpileup_set(args->mplp, MAX_BQ, 30);
                args->flag &= ~MPLP_REALN;
                args->flag |= MPLP_NO_INDEL;
            } else if (strcasecmp(optarg, "1.12") == 0) {
                // 1.12 and earlier
                mpileup_set(args->mplp, MIN_REALN_FRAC, 0.002);
                args->bca.min_support = 1;
                mpileup_set(args->mplp, MIN_BQ, 13);
                args->bca.tandemQ = 100;
                args->flag &= ~MPLP_REALN_PARTIAL;
                args->flag |= MPLP_REALN;
            } else if (strcasecmp(optarg, "illumina") == 0) {
                args->flag |= MPLP_REALN_PARTIAL;
            } else {
                fprintf(stderr, "Unknown configuration name '%s'\n"
                        "Please choose from 1.12, illumina, pacbio-ccs or ont\n",
                        optarg);
                return 1;
            }
            break;
        case 13: hts_srand48(atoi(optarg)); break;
        case 14:
            error("todo: --ambig-reads\n");
            // if ( !strcasecmp(optarg,"drop") ) .ambig_reads = B2B_DROP;
            // else if ( !strcasecmp(optarg,"incAD") ) mplp.ambig_reads = B2B_INC_AD;
            // else if ( !strcasecmp(optarg,"incAD0") ) mplp.ambig_reads = B2B_INC_AD0;
            // else error("The option to --ambig-reads not recognised: %s\n",optarg);
            break;
        default:
            fprintf(stderr,"Invalid option: '%c'\n", c);
            return 1;
        }
    }

    if ( args->gvcf && !(args->bca.fmt_flag&B2B_FMT_DP) ) args->bca.fmt_flag |= B2B_FMT_DP;
    if ( argc==1 )
    {
        print_usage(args, stderr);
        clean_up(args);
        return 1;
    }
    int i;
    for (i=0; i<argc-optind; i++)
        if ( mpileup_set(args->mplp, BAM, argv[optind+i])!=0 ) clean_up(args);
    mpileup_init(args->mplp);

    ret = run_mpileup(args);

    clean_up(args);
    return ret;
}
