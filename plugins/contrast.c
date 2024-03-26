/* The MIT License

   Copyright (c) 2018-2023 Genome Research Ltd.

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
   THE SOFTWARE.

 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <strings.h>
#include <errno.h>
#include <unistd.h>     // for isatty
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/kfunc.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"
#include "filter.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define PRINT_PASSOC  (1<<0)
#define PRINT_FASSOC  (1<<1)
#define PRINT_NASSOC  (1<<2)
#define PRINT_NOVELAL (1<<3)
#define PRINT_NOVELGT (1<<4)

typedef struct
{
    int argc, filter_logic, regions_is_file, targets_is_file, output_type, force_samples, clevel;
    int regions_overlap, targets_overlap;
    uint32_t annots;
    char **argv, *output_fname, *fname, *regions, *targets, *filter_str, *annots_str;
    char *control_samples_str, *case_samples_str, *max_AC_str;
    int *control_smpl, *case_smpl, ncontrol_smpl, ncase_smpl;
    filter_t *filter;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fh;
    int32_t *gts;
    int mgts;
    uint32_t *control_gts;
    int ncontrol_gts, mcontrol_gts, ntotal, nskipped, ntested, ncase_al, ncase_gt;
    kstring_t case_als_smpl, case_gts_smpl;
    int max_AC, nals[4];    // nals: number of control-ref, control-alt, case-ref and case-alt alleles in the region
    char *index_fn;
    int write_index;
}
args_t;

args_t args;

const char *about(void)
{
    return "Find novel alleles and genotypes in two groups of samples.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Runs a basic association test, per-site or in a region, and checks for novel alleles and\n"
        "       genotypes in two groups of samples. Adds the following INFO annotations:\n"
        "       - PASSOC  .. Fisher's exact test probability of genotypic association (REF vs non-REF allele)\n"
        "       - FASSOC  .. proportion of non-REF allele in controls and cases\n"
        "       - NASSOC  .. number of control-ref, control-alt, case-ref and case-alt alleles\n"
        "       - NOVELAL .. lists samples with a novel allele not observed in the control group\n"
        "       - NOVELGT .. lists samples with a novel genotype not observed in the control group\n"
        "Usage: bcftools +contrast [Plugin Options]\n"
        "Plugin options:\n"
        "   -a, --annots LIST                List of annotations to output [PASSOC,FASSOC,NOVELAL]\n"
        "   -0, --control-samples LIST|FILE  File or comma-separated list of control (background) samples\n"
        "   -1, --case-samples LIST|FILE     File or comma-separated list of samples where novel allele or genotype is expected\n"
        "   -e, --exclude EXPR               Exclude sites and samples for which the expression is true\n"
        "   -f, --max-allele-freq NUM        Calculate enrichment of rare alleles. Floating point numbers between 0 and 1 are\n"
        "                                        interpreted as ALT allele frequencies, integers as ALT allele counts\n"
        "       --force-samples              Continue even if some samples listed in the -0,-1 files are missing from the VCF\n"
        "   -i, --include EXPR               Include sites and samples for which the expression is true\n"
        "   -o, --output FILE                Output file name [stdout]\n"
        "   -O, --output-type u|b|v|z[0-9]   u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -r, --regions REG                Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE          Restrict to regions listed in a file\n"
        "       --regions-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n"
        "   -t, --targets REG                Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE          Similar to -R but streams rather than index-jumps\n"
        "       --targets-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n"
        "   -W, --write-index[=FMT]          Automatically index the output files [off]\n"
        "\n"
        "Example:\n"
        "   # Test if any of the samples a,b is different from the samples c,d,e\n"
        "   bcftools +contrast -0 c,d,e -1 a,b file.bcf\n"
        "\n"
        "   # Same as above, but read samples from a file. In case of a name collision, the sample name\n"
        "   # has precedence: the existence of a file with a list of samples is not checked unless no such\n"
        "   # sample exists in the VCF. Use a full path (e.g. \"./string\" instead of \"string\") to avoid\n"
        "   # name clashes\n"
        "   bcftools +contrast -0 samples0.txt -1 samples1.txt file.bcf\n"
        "\n"
        "   # The same as above but checks for enrichment of rare alleles, AF<0.001 in this example, in a region\n"
        "   bcftools +contrast -r 20:1000-2000 -f 0.001 -0 samples0.txt -1 samples1.txt file.bcf\n"
        "\n";
}

static int cmp_int(const void *a, const void *b)
{
    if ( *((int*)a) < *((int*)b) ) return -1;
    if ( *((int*)a) > *((int*)b) ) return -1;
    return 0;
}
static void read_sample_list_or_file(bcf_hdr_t *hdr, const char *str, int **smpl, int *nsmpl, int force_samples)
{
    char **str_list = NULL;
    int i,j, *list = NULL, nlist = 0, is_file, nskipped = 0;

    for (is_file=0; is_file<=1; is_file++)
    {
        if ( str_list )
        {
            for (i=0; i<nlist; i++) free(str_list[i]);
            free(str_list);
            free(list);
            str_list = NULL;
            list = NULL;
            nlist = 0;
        }

        str_list = hts_readlist(str, is_file, &nlist);
        if ( !str_list )
        {
            if ( force_samples ) continue;
            error("The sample \"%s\", is not present in the VCF\n", str);
        }

        list = (int*) malloc(sizeof(int)*nlist);
        for (i=0,j=0; i<nlist; i++,j++)
        {
            list[j] = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, str_list[i]);
            if ( list[j] >= 0 ) continue;
            if ( is_file )
            {
                if ( !force_samples ) error("The sample \"%s\" is not present in the VCF. Use --force-samples to proceed anyway.\n", str_list[i]);
                j--;
                nskipped++;
                continue;
            }
            break;
        }
        if ( i==nlist ) break;
    }
    for (i=0; i<nlist; i++) free(str_list[i]);
    nlist -= nskipped;
    if ( !nlist && !force_samples ) error("None of the samples are present in the VCF: %s\n", str);
    if ( nskipped ) fprintf(stderr,"Warning: using %d sample%s, %d from %s %s not present in the VCF\n", nlist,nlist>1?"s":"",nskipped,str,nskipped>1?"are":"is");
    free(str_list);
    qsort(list,nlist,sizeof(*list),cmp_int);
    *smpl = list;
    *nsmpl = nlist;
}

static void init_data(args_t *args)
{
    int ntmp, i;
    char **tmp = hts_readlist(args->annots_str, 0, &ntmp);
    for (i=0; i<ntmp; i++)
    {
        if ( !strcasecmp("PASSOC",tmp[i]) ) args->annots |= PRINT_PASSOC;
        else if ( !strcasecmp("FASSOC",tmp[i]) ) args->annots |= PRINT_FASSOC;
        else if ( !strcasecmp("NASSOC",tmp[i]) ) args->annots |= PRINT_NASSOC;
        else if ( !strcasecmp("NOVELAL",tmp[i]) ) args->annots |= PRINT_NOVELAL;
        else if ( !strcasecmp("NOVELGT",tmp[i]) ) args->annots |= PRINT_NOVELGT;
        else error("The annotation is not recognised: %s\n", tmp[i]);
        free(tmp[i]);
    }
    free(tmp);

    args->sr = bcf_sr_init();
    if ( args->regions )
    {
        args->sr->require_index = 1;
        bcf_sr_set_opt(args->sr,BCF_SR_REGIONS_OVERLAP,args->regions_overlap);
        if ( bcf_sr_set_regions(args->sr, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n",args->regions);
    }
    if ( args->targets )
    {
        bcf_sr_set_opt(args->sr,BCF_SR_TARGETS_OVERLAP,args->targets_overlap);
        if ( bcf_sr_set_targets(args->sr, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->targets);
    }
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);
    args->hdr_out = bcf_hdr_dup(args->hdr);
    if ( args->annots & PRINT_PASSOC )
        bcf_hdr_append(args->hdr_out, "##INFO=<ID=PASSOC,Number=1,Type=Float,Description=\"Fisher's exact test probability of genotypic association (REF vs non-REF allele)\">");
    if ( args->annots & PRINT_FASSOC )
        bcf_hdr_append(args->hdr_out, "##INFO=<ID=FASSOC,Number=2,Type=Float,Description=\"Proportion of non-REF allele in controls and cases\">");
    if ( args->annots & PRINT_NASSOC )
        bcf_hdr_append(args->hdr_out, "##INFO=<ID=NASSOC,Number=4,Type=Integer,Description=\"Number of control-ref, control-alt, case-ref and case-alt alleles\">");
    if ( args->annots & PRINT_NOVELAL )
        bcf_hdr_append(args->hdr_out, "##INFO=<ID=NOVELAL,Number=.,Type=String,Description=\"List of samples with novel alleles. Note that samples listed here are not listed in NOVELGT again.\">");
    if ( args->annots & PRINT_NOVELGT )
        bcf_hdr_append(args->hdr_out, "##INFO=<ID=NOVELGT,Number=.,Type=String,Description=\"List of samples with novel genotypes\">");

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    read_sample_list_or_file(args->hdr, args->control_samples_str, &args->control_smpl, &args->ncontrol_smpl, args->force_samples);
    read_sample_list_or_file(args->hdr, args->case_samples_str, &args->case_smpl, &args->ncase_smpl, args->force_samples);

    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( bcf_hdr_write(args->out_fh, args->hdr_out)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    if ( init_index2(args->out_fh,args->hdr_out,args->output_fname,
                     &args->index_fn, args->write_index)<0 )
        error("Error: failed to initialise index for %s\n",args->output_fname);

    if ( args->max_AC_str )
    {
        char *tmp;
        args->max_AC = strtol(args->max_AC_str, &tmp, 10);
        if ( tmp==args->max_AC_str || *tmp )
        {
            double val = strtod(args->max_AC_str, &tmp);
            if ( tmp==args->max_AC_str || *tmp ) error("Could not parse the argument: -f, --max-allele-freq %s\n", args->max_AC_str);
            if ( val<0 || val>1 ) error("Expected integer or float from the range [0,1]: -f, --max-allele-freq %s\n", args->max_AC_str);
            args->max_AC = val * bcf_hdr_nsamples(args->hdr);
            if ( !args->max_AC ) args->max_AC = 1;
        }
    }
}
static void destroy_data(args_t *args)
{
    bcf_hdr_destroy(args->hdr_out);
    if ( args->write_index )
    {
        if ( bcf_idx_save(args->out_fh)<0 )
        {
            if ( hts_close(args->out_fh)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
            error("Error: cannot write to index %s\n", args->index_fn);
        }
        free(args->index_fn);
    }
    if ( hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    free(args->case_als_smpl.s);
    free(args->case_gts_smpl.s);
    free(args->gts);
    free(args->control_gts);
    free(args->control_smpl);
    free(args->case_smpl);
    if ( args->filter ) filter_destroy(args->filter);
    bcf_sr_destroy(args->sr);
    free(args);
}
static inline int binary_search(uint32_t val, uint32_t *dat, int ndat)
{
    int i = -1, imin = 0, imax = ndat - 1;
    while ( imin<=imax )
    {
        i = (imin+imax)/2;
        if ( dat[i] < val ) imin = i + 1;
        else if ( dat[i] > val ) imax = i - 1;
        else return 1;
    }
    return 0;
}
static inline void binary_insert(uint32_t val, uint32_t **dat, int *ndat, int *mdat)
{
    int i = -1, imin = 0, imax = *ndat - 1;
    while ( imin<=imax )
    {
        i = (imin+imax)/2;
        if ( (*dat)[i] < val ) imin = i + 1;
        else if ( (*dat)[i] > val ) imax = i - 1;
        else return;
    }
    while ( i>=0 && (*dat)[i]>val ) i--;

    (*ndat)++;
    hts_expand(uint32_t, (*ndat), (*mdat), (*dat));

    if ( *ndat > 1 )
        memmove(*dat + i + 1, *dat + i, sizeof(uint32_t)*(*ndat - i - 1));

    (*dat)[i+1] = val;
}
static int process_record(args_t *args, bcf1_t *rec)
{
    args->ntotal++;

    static int warned = 0;
    int ngts = bcf_get_genotypes(args->hdr, rec, &args->gts, &args->mgts);
    ngts /= rec->n_sample;
    if ( ngts>2 ) error("todo: ploidy=%d\n", ngts);

    args->ncontrol_gts = 0;
    uint32_t control_als = 0;
    int32_t nals[4] = {0,0,0,0};    // ctrl-ref, ctrl-alt, case-ref, case-alt
    int i,j;
    for (i=0; i<args->ncontrol_smpl; i++)
    {
        uint32_t gt  = 0;
        int32_t *ptr = args->gts + args->control_smpl[i]*ngts;
        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(ptr[j]) ) continue;
            int ial = bcf_gt_allele(ptr[j]);
            if ( ial > 31 )
            {
                if ( !warned )
                {
                    fprintf(stderr,"Too many alleles (>32) at %s:%"PRId64", skipping the site.\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                    warned = 1;
                }
                args->nskipped++;
                return -1;
            }
            control_als |= 1<<ial;
            gt |= 1<<ial;
            if ( ial ) nals[1]++;
            else nals[0]++;
        }
        if ( args->annots & PRINT_NOVELGT )
            binary_insert(gt, &args->control_gts, &args->ncontrol_gts, &args->mcontrol_gts);
    }
    if ( !control_als && args->ncontrol_smpl )
    {
        // all are missing
        args->nskipped++;
        return -1;
    }

    args->case_als_smpl.l = 0;
    args->case_gts_smpl.l = 0;

    int has_gt = 0;
    for (i=0; i<args->ncase_smpl; i++)
    {
        int case_al = 0;
        uint32_t gt  = 0;
        int32_t *ptr = args->gts + args->case_smpl[i]*ngts;
        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(ptr[j]) ) continue;
            int ial = bcf_gt_allele(ptr[j]);
            if ( ial > 31 )
            {
                if ( !warned )
                {
                    fprintf(stderr,"Too many alleles (>32) at %s:%"PRId64", skipping. (todo?)\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                    warned = 1;
                }
                args->nskipped++;
                return -1;
            }
            if ( !(control_als & (1<<ial)) ) case_al = 1;
            gt |= 1<<ial;
            if ( ial ) nals[3]++;
            else nals[2]++;
        }
        if ( !gt ) continue;
        has_gt = 1;

        char *smpl = args->hdr->samples[ args->case_smpl[i] ];
        if ( case_al && (args->annots & PRINT_NOVELAL) )
        {
            if ( args->case_als_smpl.l ) kputc(',', &args->case_als_smpl);
            kputs(smpl, &args->case_als_smpl);
        }
        else if ( (args->annots & PRINT_NOVELGT) && !binary_search(gt, args->control_gts, args->ncontrol_gts) )
        {
            if ( args->case_gts_smpl.l ) kputc(',', &args->case_gts_smpl);
            kputs(smpl, &args->case_gts_smpl);
        }
    }
    if ( !has_gt && args->ncase_smpl )
    {
        // all are missing
        args->nskipped++;
        return -1;
    }

    if ( args->max_AC )
    {
        if ( nals[0]+nals[2] > nals[1]+nals[3] )
        {
            if ( nals[1]+nals[3] <= args->max_AC )
                for (i=0; i<4; i++) args->nals[i] += nals[i];
        }
        else
        {
            if ( nals[0]+nals[2] <= args->max_AC )
            {
                args->nals[0] += nals[1];
                args->nals[1] += nals[0];
                args->nals[2] += nals[3];
                args->nals[3] += nals[2];
            }
        }
    }

    float vals[2];
    if ( (args->annots & PRINT_PASSOC) && args->ncontrol_smpl && args->ncase_smpl )
    {
        double left, right, fisher;
        kt_fisher_exact(nals[0],nals[1],nals[2],nals[3], &left,&right,&fisher);
        vals[0] = fisher;
        bcf_update_info_float(args->hdr_out, rec, "PASSOC", vals, 1);
    }
    if ( (args->annots & PRINT_FASSOC) && args->ncontrol_smpl && args->ncase_smpl )
    {
        if ( nals[0]+nals[1] ) vals[0] = (float)nals[1]/(nals[0]+nals[1]);
        else bcf_float_set_missing(vals[0]);
        if ( nals[2]+nals[3] ) vals[1] = (float)nals[3]/(nals[2]+nals[3]);
        else bcf_float_set_missing(vals[1]);
        bcf_update_info_float(args->hdr_out, rec, "FASSOC", vals, 2);
    }
    if ( args->annots & PRINT_NASSOC )
        bcf_update_info_int32(args->hdr_out, rec, "NASSOC", nals, 4);

    if ( args->case_als_smpl.l )
    {
        bcf_update_info_string(args->hdr_out, rec, "NOVELAL", args->case_als_smpl.s);
        args->ncase_al++;
    }
    if ( args->case_gts_smpl.l )
    {
        bcf_update_info_string(args->hdr_out, rec, "NOVELGT", args->case_gts_smpl.s);
        args->ncase_gt++;
    }
    args->ntested++;
    return 0;
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->annots_str = "PASSOC,FASSOC";
    args->regions_overlap = 1;
    args->targets_overlap = 0;
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"max-allele-freq",required_argument,0,'f'},
        {"annots",required_argument,0,'a'},
        {"force-samples",no_argument,0,1},
        {"bg-samples",required_argument,0,'0'},     // renamed to --control-samples, leaving it in for backward compatibility
        {"control-samples",required_argument,0,'0'},
        {"novel-samples",required_argument,0,'1'},  // renamed to --case-samples, leaving it in for backward compatibility
        {"case-samples",required_argument,0,'1'},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,4},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "O:o:i:e:r:R:t:T:0:1:a:f:W::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->force_samples = 1; break;
            case 'f': args->max_AC_str = optarg; break;
            case 'a': args->annots_str = optarg; break;
            case '0': args->control_samples_str = optarg; break;
            case '1': args->case_samples_str = optarg; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; args->targets_is_file = 1; break;
            case 'r': args->regions = optarg; break;
            case 'R': args->regions = optarg; args->regions_is_file = 1; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      };
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case  3 :
                args->regions_overlap = parse_overlap_option(optarg);
                if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  4 :
                args->targets_overlap = parse_overlap_option(optarg);
                if ( args->targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s",usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s",usage_text());
    else args->fname = argv[optind];

    if ( !args->control_samples_str ) error("Error: missing the -0, --control-samples option\n");
    if ( !args->case_samples_str ) error("Error: missing the -1, --case-samples option\n");

    init_data(args);

    while ( bcf_sr_next_line(args->sr) )
    {
        bcf1_t *rec = bcf_sr_get_line(args->sr,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, rec, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }
        process_record(args, rec);
        if ( bcf_write(args->out_fh, args->hdr_out, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    }

    fprintf(stderr,"Total/processed/skipped/case_allele/case_gt:\t%d\t%d\t%d\t%d\t%d\n", args->ntotal, args->ntested, args->nskipped, args->ncase_al, args->ncase_gt);
    if ( args->max_AC )
    {
        double val1, val2, fisher;
        kt_fisher_exact(args->nals[0],args->nals[1],args->nals[2],args->nals[3], &val1,&val2,&fisher);
        val1 = args->nals[0]+args->nals[1] ? (float)args->nals[1]/(args->nals[0]+args->nals[1]) : 0;
        val2 = args->nals[2]+args->nals[3] ? (float)args->nals[3]/(args->nals[2]+args->nals[3]) : 0;
        fprintf(stderr,"max_AC/PASSOC/FASSOC/NASSOC:\t%d\t%e\t%f,%f\t%d,%d,%d,%d\n",args->max_AC,fisher,val1,val2,args->nals[0],args->nals[1],args->nals[2],args->nals[3]);
    }
    destroy_data(args);

    return 0;
}
