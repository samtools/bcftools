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

int main_samview(int argc, char *argv[]);
int main_vcfview(int argc, char *argv[]);
int main_bcfidx(int argc, char *argv[]);
int main_tabix(int argc, char *argv[]);
int main_vcfcheck(int argc, char *argv[]);
int main_vcfisec(int argc, char *argv[]);
int main_vcfmerge(int argc, char *argv[]);
int main_vcfquery(int argc, char *argv[]);
int main_vcffilter(int argc, char *argv[]);
int main_vcfnorm(int argc, char *argv[]);
int main_vcfgtcheck(int argc, char *argv[]);
int main_vcfsubset(int argc, char *argv[]);
int main_vcfcall(int argc, char *argv[]);

typedef struct
{
    int (*func)(int, char*[]);
    const char *alias, *help, *sep;
}
cmd_t;

static cmd_t cmds[] =
{
    { .func  = main_tabix,    
      .alias = "tabix",
      .help  = "tabix for BGZF'd BED, GFF, SAM, VCF and more",
      .sep   = NULL
    },
    { .func = main_bcfidx,   
      .alias = "index",
      .help = "index BCF",
      .sep   = NULL
    },
    { .func  = main_vcfcall,  
      .alias = "call", 
      .help  = "SNP/indel calling (former \"view\"; this version is broken)",
      .sep   = "VCF/BCF tools:"
    },
    { .func  = main_vcfcheck, 
      .alias = "check",
      .help  = "produce VCF stats",
      .sep   = NULL
    },
    { .func  = main_vcffilter, 
      .alias = "filter",
      .help  = "filter VCF files (this version is broken)",
      .sep   = NULL
    },
    { .func  = main_vcfgtcheck, 
      .alias = "gtcheck",
      .help  = "tool for detecting swaps and contaminations",
      .sep   = NULL
    },
    { .func  = main_vcfisec,  
      .alias = "isec", 
      .help  = "intersections of VCF files",
      .sep   = NULL
    },
    { .func  = main_vcfmerge, 
      .alias = "merge",
      .help  = "merge VCF files",
      .sep   = NULL
    },
    { .func  = main_vcfnorm, 
      .alias = "norm",
      .help  = "normalize indels",
      .sep   = NULL
    },
    { .func  = main_vcfquery, 
      .alias = "query",
      .help  = "transform VCF into user-defined formats",
      .sep   = NULL
    },
    { .func  = main_vcfsubset, 
      .alias = "subset",
      .help  = "subset and filter vcf and bcf",
      .sep   = NULL
    },
    { .func  = main_vcfview,  
      .alias = "view", 
      .help  = "VCF<->BCF conversion",
      .sep   = NULL
    },
    { .func  = NULL,
      .alias = NULL,
      .help  = NULL,
      .sep   = NULL
    }
};

char *bcftools_version_string = NULL;

char *bcftools_version(void)
{
    if ( !bcftools_version_string )
    {
        int len = strlen(hts_version()) + strlen(BCFTOOLS_VERSION) + 2;
        bcftools_version_string = (char*) malloc(len);
        snprintf(bcftools_version_string,len,"%s:%s", BCFTOOLS_VERSION,hts_version());
    }
    return bcftools_version_string;
}

static int usage(char *argv0)
{
	fprintf(stderr, "\n");
    fprintf(stderr, "Version: %s\n", bcftools_version());
	fprintf(stderr, "Usage:   %s <command> <argument>\n", argv0);
	fprintf(stderr, "Commands:\n");

    int i = 0;
    const char *sep = NULL;
    while (cmds[i].func)
    {
        if ( cmds[i].sep ) sep = cmds[i].sep;
        if ( sep )
        {
            printf("\n -- %s\n", sep);
            sep = NULL;
        }
        printf("\t%-15s %s\n", cmds[i].alias, cmds[i].help);
        i++;
    }

    fprintf(stderr,"\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage(argv[0]);

    int i = 0;
    while (cmds[i].func)
    {
        if ( !strcmp(argv[1],cmds[i].alias) ) 
        {
            int ret = cmds[i].func(argc-1,argv+1);
            free(bcftools_version_string);
            return ret;
        }
        i++;
    }
    fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
    return 1;
}

