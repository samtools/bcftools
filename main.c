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

int main_bcfidx(int argc, char *argv[]);
int main_tabix(int argc, char *argv[]);
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

typedef struct
{
    int (*func)(int, char*[]);
    const char *alias, *help;
}
cmd_t;

static cmd_t cmds[] =
{
    { .func  = main_tabix,    
      .alias = "tabix",
      .help  = "tabix for BGZF'd BED, GFF, SAM, VCF and more"
    },
    { .func = main_bcfidx,   
      .alias = "index",
      .help = "index BCF"
    },
    { .func  = NULL, 
      .alias = "Core VCF/BCF tools:",
      .help  = NULL
    },
    { .func  = main_vcfcall,  
      .alias = "call", 
      .help  = "SNP/indel calling (former \"view\")"
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
      .help  = "merge VCF/BCF files"
    },
    { .func  = main_vcfnorm, 
      .alias = "norm",
      .help  = "normalize indels"
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
    { .func  = main_vcfannotate,  
      .alias = "annotate", 
      .help  = "-annotate and edit VCF/BCF files",  // do not advertise yet
    },
    { .func  = main_vcfroh, 
      .alias = "roh",
      .help  = "-identify runs of autozygosity (HMM)",  // do not advertise yet
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

char *bcftools_version_string = NULL;

char *bcftools_version(void)
{
    if ( !bcftools_version_string )
    {
        int len = strlen(hts_version()) + strlen(BCFTOOLS_VERSION) + 9;
        bcftools_version_string = (char*) malloc(len);
        snprintf(bcftools_version_string,len,"%s+htslib-%s", BCFTOOLS_VERSION,hts_version());
    }
    return bcftools_version_string;
}

static int usage(void)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)\n");
  fprintf(stderr, "Version: %s\n", bcftools_version());
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   bcftools <command> <argument>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Commands:\n");

    int i = 0;
    const char *sep = NULL;
    while (cmds[i].alias)
    {
        if ( !cmds[i].func ) sep = cmds[i].alias;
        if ( sep )
        {
            printf("\n -- %s\n", sep);
            sep = NULL;
        }
        if ( cmds[i].func && cmds[i].help[0]!='-' ) printf("\t%-15s %s\n", cmds[i].alias, cmds[i].help);
        i++;
    }

    fprintf(stderr,"\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();

    int i = 0;
    while (cmds[i].alias)
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

