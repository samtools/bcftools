#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct _args_t
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // one of FLT_IN/EXCLUDE (-i or -e)

    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hnull, *hsub; // original header, sites-only header, subset header
    char **argv, *format, *sample_names, *subset_fname, *targets_fname, *regions_fname;
    int argc, clevel, output_type, print_header, update_info, header_only, n_samples, *imap;
    int trim_alts, sites_only, known, novel, multiallelic, biallelic, exclude_ref, private_vars, exclude_uncalled, min_ac, max_ac, calc_ac;
    char *fn_ref, *fn_out, **samples;
    char *include_types, *exclude_types;
    int include, exclude;
    htsFile *out;
}
args_t;

static void init_data(args_t *args)
{
    int i;
    args->hdr = args->files->readers[0].header;
    
    if (args->calc_ac && args->update_info)
    {
        bcf_hdr_append(args->hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
        bcf_hdr_append(args->hdr,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    }
    bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_subset");

    // setup sample data    
    if (args->sample_names)
    {
        struct stat sbuf;
        if ( stat(args->sample_names, &sbuf) == 0  )
            args->samples = hts_readlines(args->sample_names, &args->n_samples);
        else
        {
            int m = 0, n = 0;
            char **s = 0;
            const char *q, *p;
            for (q = p = args->sample_names;; ++p) {
                if (*p == ',' || *p == 0) {
                    if (m == n) {
                        m = m? m<<1 : 16;
                        s = (char**)realloc(s, m * sizeof(void*));
                    }
                    s[n] = (char*)calloc(p - q + 1, 1);
                    strncpy(s[n++], q, p - q);
                    q = p + 1;
                    if (*p == 0) break;
                }
            }
            s = (char**)realloc(s, n * sizeof(void*));
            args->samples = s;
            args->n_samples = n;
        }
    }
    
    if (args->n_samples)
        args->imap = (int*)malloc(args->n_samples * sizeof(int));
    
    // determine variant types to include/exclude
    if (args->include_types || args->exclude_types) {
        if (args->include_types && args->exclude_types) {
            fprintf(stderr, "Error: only supply one of --include_types --exclude-types options\n");
            exit(1);
        }
        char **type_list = 0;
        int m = 0, n = 0;
        const char *q, *p;
        for (q = p = args->include_types ? args->include_types : args->exclude_types;; ++p) {
            if (*p == ',' || *p == 0) {
                if (m == n) {
                    m = m? m<<1 : 16;
                    type_list = (char**)realloc(type_list, m * sizeof(void*));
                }
                type_list[n] = (char*)calloc(p - q + 1, 1);
                strncpy(type_list[n++], q, p - q);
                q = p + 1;
                if (*p == 0) break;
            }
        }
        type_list = (char**)realloc(type_list, n * sizeof(void*));

        if (args->include_types) {
            args->include = 0;
            for (i = 0; i < n; ++i) {
                if (strcmp(type_list[i], "snps") == 0) args->include |= VCF_SNP;
                else if (strcmp(type_list[i], "indels") == 0) args->include |= VCF_INDEL;
                else if (strcmp(type_list[i], "mnps") == 0) args->include |= VCF_MNP;
                else if (strcmp(type_list[i], "other") == 0) args->include |= VCF_OTHER;
                else {
                    fprintf(stderr, "[E::%s] unknown type\n", type_list[i]);
                    exit(1);
                }
            }
        }
        if (args->exclude_types) {
            args->exclude = 0;
            for (i = 0; i < n; ++i) {
                if (strcmp(type_list[i], "snps") == 0) args->exclude |= VCF_SNP;
                else if (strcmp(type_list[i], "indels") == 0) args->exclude |= VCF_INDEL;
                else if (strcmp(type_list[i], "mnps") == 0) args->exclude |= VCF_MNP;
                else if (strcmp(type_list[i], "other") == 0) args->exclude |= VCF_OTHER;
                else {
                    fprintf(stderr, "[E::%s] unknown type\n", type_list[i]);
                    exit(1);
                }
            }
        }
        for (i = 0; i < n; ++i)
            free(type_list[i]);
        free(type_list);
    }

    // setup output
    char modew[8];
    strcpy(modew, "w");
    if (args->clevel >= 0 && args->clevel <= 9) sprintf(modew + 1, "%d", args->clevel);
    if (args->output_type==FT_BCF) strcat(modew, "bu");         // uncompressed BCF
    else if (args->output_type & FT_BCF) strcat(modew, "b");    // compressed BCF
    else if (args->output_type & FT_GZ) strcat(modew,"z");      // compressed VCF
    args->out = hts_open(args->fn_out ? args->fn_out : "-", modew);
 
    // headers: hdr=full header, hsub=subset header, hnull=sites only header
    if (args->sites_only)
        args->hnull = bcf_hdr_subset(args->hdr, 0, 0, 0);
    if (args->n_samples > 0)
        args->hsub = bcf_hdr_subset(args->hdr, args->n_samples, args->samples, args->imap);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    bcf_hdr_fmt_text(args->hdr);
}

static void destroy_data(args_t *args)
{
    int i;
    if ( args->imap ) {
        for (i = 0; i < args->n_samples; ++i)
            free(args->samples[i]);
        free(args->samples);
        free(args->imap);
    }
    if (args->hnull) bcf_hdr_destroy(args->hnull);
    if (args->hsub) bcf_hdr_destroy(args->hsub);
    if ( args->filter )
        filter_destroy(args->filter);
}

int subset_vcf(args_t *args, bcf1_t *line)
{
    if ( args->multiallelic && !(line->n_allele>2) ) return 0; // select multiallelic sites
    if ( args->biallelic && !(line->n_allele==2) ) return 0; // skip multiallelic sites
    if (args->novel || args->known)
    {
        if ( args->novel && (line->d.id[0]!='.' || line->d.id[1]!=0) ) return 0; // skip sites which are known, ID != '.'
        if ( args->known && line->d.id[0]=='.' && line->d.id[1]==0 ) return 0;  // skip sites which are novel, ID == '.'
    }
                
    if (args->include || args->exclude)
    {
        int line_type = bcf_get_variant_types(line);
        if ( args->include && !(line_type&args->include) ) return 0; // include only given variant types
        if ( args->exclude &&   line_type&args->exclude  ) return 0; // exclude given variant types
    }

    if ( args->filter )
    {
        int ret = filter_test(args->filter, line);
        if ( args->filter_logic==FLT_INCLUDE ) { if ( !ret ) return 0; }
        else if ( ret ) return 0;
    }

    int i, an = 0, n_ac = 0, *ac = (int*) calloc(line->n_allele+1,sizeof(int));
    if (args->calc_ac) {
        bcf_calc_ac(args->hdr, line, ac, BCF_UN_INFO|BCF_UN_FMT); // get original AC and AN values from INFO field if available, otherwise calculate
        for (i=1; i<=line->n_allele; i++)
            n_ac += ac[i];
        for (i=0; i<line->n_allele; i++)
            an+=ac[i];
    }

    if (args->n_samples)
    {
        int n_ac_sub = 0, *ac_sub = (int*) calloc(line->n_allele+1,sizeof(int));
        bcf_subset(args->hdr, line, args->n_samples, args->imap);
        if (args->calc_ac) {
            bcf_calc_ac(args->hsub, line, ac_sub, BCF_UN_FMT); // recalculate AC and AN
            an = 0;
            for (i=0; i<line->n_allele; i++)
                an+=ac_sub[i];
            for (i=1; i<=line->n_allele; i++)
                n_ac_sub += ac_sub[i];
            if (args->private_vars && !(n_ac_sub > 0 && n_ac == n_ac_sub)) { free(ac); free(ac_sub); return 0; }
            n_ac = n_ac_sub;
            for (i=0; i<=line->n_allele; i++)
                ac[i] = ac_sub[i];
        }
        free(ac_sub);
    }
    if (args->min_ac && args->min_ac>n_ac) { free(ac); return 0; }
    if (args->max_ac && args->max_ac<n_ac) { free(ac); return 0; }
    if (args->exclude_uncalled && n_ac == 0 && ac[0] == 0) { free(ac); return 0; }
    if (args->calc_ac && args->update_info) {
        bcf_update_info_int32(args->hdr, line, "AC", &ac[1], line->n_allele-1);
        bcf_update_info_int32(args->hdr, line, "AN", &an, 1);
    }
    free(ac);
    if (args->exclude_ref && n_ac == 0) return 0;
    if (args->trim_alts) bcf_trim_alleles(args->hsub ? args->hsub : args->hdr, line);
    if (args->sites_only) bcf_subset(args->hsub ? args->hsub : args->hdr, line, 0, 0);
    return 1;
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   View, subset and filter VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools subset [options] <in.bcf>|<in.vcf>|<in.vcf.gz> [region1 [...]]\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -O, --out FILE                 output file name [stdout]\n");
    fprintf(stderr, "    -l INT                         compression level [%d]\n", args->clevel);
    fprintf(stderr, "    -h                             suppress the header in VCF output\n");
    fprintf(stderr, "    -H                             print the header only\n");
    fprintf(stderr, "    -G,                            drop individual genotype information (after subsetting if -s option set)\n");
	fprintf(stderr, "    -o, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <reg|file>       same as -t but index-jumps rather than streams to a region (requires indexed VCF/BCF)\n");
    fprintf(stderr, "    -t, --targets <reg|file>       restrict to positions in tab-delimited tabix indexed file <chr,pos> or <chr,from,to>, 1-based, inclusive\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Subset options:\n");
    fprintf(stderr, "    -a, --trim-alt-alleles      trim alternate alleles not seen in subset\n");
    fprintf(stderr, "    -I, --no-update             do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)\n");
    fprintf(stderr, "    -s, --samples STR/FILE      list of samples (FILE or comma separated list STR) [null]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    -R,   --exclude-ref                            exclude sites without a non-reference genotype\n");
    fprintf(stderr, "    -U,   --exclude-uncalled                       exclude sites without a called genotype\n");
    fprintf(stderr, "    -p,   --private                                print sites where only the subset samples carry an non-reference allele\n");
    fprintf(stderr, "    -i/e, --include/exclude-filters <expr>         include/exclude sites for which the expression is true (see vcffilter for details)\n");
    fprintf(stderr, "    -f,   --apply-filters <list>                   require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -k/n, --known/--novel                          print known/novel sites only (ID is/not '.')\n");
    fprintf(stderr, "    -m/M, --multiallelic/--biallelic               print multiallelic/biallelic sites only\n");
    fprintf(stderr, "    -c/C  --min-ac/--max-ac                        minimum/maximum allele count (INFO/AC) of sites to be printed\n");
    fprintf(stderr, "    -1/2  --singletons/--doubletons                print singleton/doubleton sites only (shortcut for -c1 -C1/-c2 -C2)\n");
    fprintf(stderr, "    -v/V  --include-types/--exclude-types STR      comma-separated list of variant types to include/exclude: snps,indels,mnps,other [null]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfsubset(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->clevel  = -1;
    args->print_header = 1;
    args->update_info = 1;
    args->output_type = FT_VCF;

    static struct option loptions[] = 
    {
        {"exclude-filters",1,0,'e'},
        {"include-filters",1,0,'i'},
        {"trim-alt-alleles",0,0,'a'},
        {"exclude-ref",0,0,'R'},
        {"no-update",0,0,'I'},
        {"private",0,0,'p'},
        {"exclude-uncalled",0,0,'U'},
        {"apply-filters",1,0,'f'},
        {"known",0,0,'k'},
        {"novel",0,0,'n'},
        {"multiallelic",0,0,'m'},
        {"biallelic",0,0,'M'},
        {"samples",1,0,'s'},
        {"output-type",1,0,'o'},
        {"out",1,0,'O'},
        {"file-name",1,0,'n'},
        {"include-types",1,0,'v'},
        {"exclude-types",1,0,'V'},
        {"targets",1,0,'t'},
        {"regions",1,0,'r'},
        {"min-ac",1,0,'c'},
        {"max-ac",1,0,'C'},
        {"singletons",0,0,'1'},
        {"doubletons",0,0,'2'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "l:St:r:o:O:s:Gf:knv:V:mMaRpUhHc:C:12Ie:i:",loptions,NULL)) >= 0) {
        switch (c) {
    	    case 'o': 
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'l': args->clevel = atoi(optarg); args->output_type |= FT_GZ; break;
            case 'O': args->fn_out = optarg; break;
            case 'h': args->print_header = 0; break;
            case 'H': args->header_only = 1; break;
            
            case 't': args->targets_fname = optarg; break;
            case 'r': args->regions_fname = optarg; break;
            
            case 's': args->sample_names = optarg; break;
            case 'a': args->trim_alts = 1; args->calc_ac = 1; break;
            case 'I': args->update_info = 0; break;
            case 'G': args->sites_only = 1; break;
            
            case 'f': args->files->apply_filters = optarg; break;
            case 'k': args->known = 1; break;
            case 'n': args->novel = 1; break;
            case 'm': args->multiallelic = 1; break;
            case 'M': args->biallelic = 1; break;
            case 'v': args->include_types = optarg; break;
            case 'V': args->exclude_types = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;

            case 'c': args->min_ac = atoi(optarg); args->calc_ac = 1; break;
            case 'C': args->max_ac = atoi(optarg); args->calc_ac = 1; break;
            case '1': args->min_ac = 1; args->max_ac = 1; args->calc_ac = 1; break;
            case '2': args->min_ac = 2; args->max_ac = 2; args->calc_ac = 1; break;

            case 'R': args->exclude_ref = 1; args->calc_ac = 1; break;
            case 'p': args->private_vars = 1; args->calc_ac = 1; break;
            case 'U': args->exclude_uncalled = 1; args->calc_ac = 1; break;
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( argc<optind+1 ) usage(args);
    // read in the regions from the command line
    if ( args->regions_fname )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_fname)<0 )
            error("Failed to read the regions: %s\n", args->regions_fname);
    }
    else if ( optind+1 < argc ) 
    {
        int i;
        kstring_t tmp = {0,0,0};
        kputs(argv[optind+1],&tmp);
        for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
        if ( bcf_sr_set_regions(args->files, tmp.s)<0 )
            error("Failed to read the regions: %s\n", tmp.s);
        free(tmp.s);
    }
    if ( args->targets_fname )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_fname,0)<0 )
            error("Failed to read the targets: %s\n", args->targets_fname);
    }

    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    
    init_data(args);
    bcf_hdr_t *out_hdr = args->hnull ? args->hnull : (args->hsub ? args->hsub : args->hdr);
    if (args->print_header)
        bcf_hdr_write(args->out, out_hdr);
    if (!args->header_only) 
    {
        while ( bcf_sr_next_line(args->files) )
        {
            bcf1_t *line = args->files->readers[0].buffer[0];
            if ( line->errcode && out_hdr!=args->hdr ) error("Undefined tags in the header, cannot proceed in the sample subset mode.\n");
            if ( subset_vcf(args, line) )
                bcf_write1(args->out, out_hdr, line);
        }
    }
    hts_close(args->out);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
