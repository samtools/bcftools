#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <htslib/hts.h>
#include "bcftools.h"

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

int main_tabix(int argc, char *argv[]);
int main_vcfindex(int argc, char *argv[]);
int main_vcfstats(int argc, char *argv[]);
int main_vcfisec(int argc, char *argv[]);
int main_vcfmerge(int argc, char *argv[]);
int main_vcfquery(int argc, char *argv[]);
int main_vcffilter(int argc, char *argv[]);
int main_vcfsom(int argc, char *argv[]);
int main_vcfnorm(int argc, char *argv[]);
int main_vcfgtcheck(int argc, char *argv[]);
int main_vcfview(int argc, char *argv[]);
int main_vcfcall(int argc, char *argv[]);
int main_vcfannotate(int argc, char *argv[]);
int main_vcfroh(int argc, char *argv[]);
int main_vcfconcat(int argc, char *argv[]);

typedef struct
{
    int (*func)(int, char*[]);
    const char *alias, *help;
}
cmd_t;

static cmd_t cmds[] =
{
    { .func  = NULL, 
      .alias = "Indexing:",
      .help  = NULL
    },
    { .func = main_vcfindex,   
      .alias = "index",
      .help = "index VCF/BCF files"
    },
    { .func = main_tabix,   
      .alias = "tabix",
      .help = "-tabix for BGZF'd BED, GFF, SAM, VCF and more" // do not advertise; only keep here for testing
    },
    { .func  = NULL, 
      .alias = "Core VCF/BCF tools:",
      .help  = NULL
    },
    { .func  = main_vcfannotate,  
      .alias = "annotate", 
      .help  = "annotate and edit VCF/BCF files",
    },
    { .func  = main_vcfcall,  
      .alias = "call", 
      .help  = "SNP/indel calling (former \"view\")"
    },
    { .func  = main_vcfconcat,  
      .alias = "concat", 
      .help  = "-combine VCF/BCF files (one-sample files yield one-sample file)"    // do not advertise yet
    },
    { .func  = main_vcffilter, 
      .alias = "filter",
      .help  = "filter VCF/BCF files using fixed thresholds"
    },
    { .func  = main_vcfgtcheck, 
      .alias = "gtcheck",
      .help  = "check sample concordance, detect sample swaps and contamination"
    },
    { .func  = main_vcfisec,  
      .alias = "isec", 
      .help  = "intersections of VCF/BCF files"
    },
    { .func  = main_vcfmerge, 
      .alias = "merge",
      .help  = "merge VCF/BCF files (one-sample files yield multi-sample file)"
    },
    { .func  = main_vcfnorm, 
      .alias = "norm",
      .help  = "left-align normalize indels"
    },
    { .func  = main_vcfquery, 
      .alias = "query",
      .help  = "transform VCF/BCF into user-defined formats"
    },
    { .func  = main_vcfstats, 
      .alias = "stats",
      .help  = "produce VCF/BCF stats (former vcfcheck)"
    },
    { .func  = main_vcfview,  
      .alias = "view", 
      .help  = "VCF/BCF conversion, view, subset and filter VCF/BCF files"
    },
    { .func  = NULL, 
      .alias = "Other/Experimental tools:" ,
      .help  = NULL
    },
    { .func  = main_vcfroh, 
      .alias = "roh",
      .help  = "identify runs of autozygosity (HMM)",
    },
    { .func  = main_vcfsom, 
      .alias = "som",
      .help  = "filter using Self-Organized Maps (experimental)"
    },
    { .func  = NULL,
      .alias = NULL,
      .help  = NULL
    }
};

char *bcftools_version(void)
{
    return BCFTOOLS_VERSION;
}

static void usage(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)\n");
    fprintf(fp, "Version: %s (using htslib %s)\n", bcftools_version(), hts_version());
    fprintf(fp, "\n");
    fprintf(fp, "Usage:   bcftools <command> <argument>\n");
    fprintf(fp, "\n");
    fprintf(fp, "Commands:\n");

    int i = 0;
    const char *sep = NULL;
    while (cmds[i].alias)
    {
        if ( !cmds[i].func ) sep = cmds[i].alias;
        if ( sep )
        {
            fprintf(fp, "\n -- %s\n", sep);
            sep = NULL;
        }
        if ( cmds[i].func && cmds[i].help[0]!='-' ) fprintf(fp, "\t%-15s %s\n", cmds[i].alias, cmds[i].help);
        i++;
    }

    fprintf(fp,"\n");
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(stderr); return 1; }

    if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0) {
        printf("bcftools %s\nUsing htslib %s\nCopyright (C) 2014 Genome Research Ltd.\n", bcftools_version(), hts_version());
        return 0;
    }
    else if (strcmp(argv[1], "--version-only") == 0) {
        printf("%s+htslib-%s\n", bcftools_version(), hts_version());
        return 0;
    }
    else if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        if (argc == 2) { usage(stdout); return 0; }
        // Otherwise change "bcftools help COMMAND [...]" to "bcftools COMMAND";
        // main_xyz() functions by convention display the subcommand's usage
        // when invoked without any arguments.
        argv++;
        argc = 2;
    }

    int i = 0;
    while (cmds[i].alias)
    {
        if ( !strcmp(argv[1],cmds[i].alias) ) 
        {
            return cmds[i].func(argc-1,argv+1);
        }
        i++;
    }
    fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
    return 1;
}

