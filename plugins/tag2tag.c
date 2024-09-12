/*  plugins/tag2tag.c -- convert between similar tags

    Copyright (C) 2014-2024 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include "bcftools.h"

typedef struct
{
    int type;
    const char *str, *hdr;
}
tag_t;
enum tag { UNKN, GP, GL, PL, GT, QRQA, QS, XX, LXX, LAA, LPL, LAD, AD };
const tag_t tags[] =
{
    [UNKN] = { .type = 0 ,          .str = NULL,  .hdr = NULL },
    [GP]   = { .type = BCF_HT_REAL, .str = "GP",  .hdr = "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype probabilities\">" },
    [GL]   = { .type = BCF_HT_REAL, .str = "GL",  .hdr = "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods\">" },
    [PL]   = { .type = BCF_HT_INT,  .str = "PL",  .hdr = "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">" },
    [GT]   = { .type = BCF_HT_STR,  .str = "GT",  .hdr = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" },
    [QRQA] = { .type = BCF_HT_INT,  .str = NULL,  .hdr = NULL },
    [QS]   = { .type = BCF_HT_INT,  .str = "QS",  .hdr = "##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Phred-score allele quality sum\">" },

    // conversion from/to localized allelles get special treatment
    [XX]   = { .type = 0 ,          .str = NULL,  .hdr = NULL },
    [LXX]  = { .type = 0 ,          .str = NULL,  .hdr = NULL },
    [LAA]  = { .type = BCF_HT_INT,  .str = "LAA", .hdr = "##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"Localized alleles, 1-based index to ALT\">" },
    [LPL]  = { .type = BCF_HT_INT,  .str = "LPL", .hdr = "##FORMAT=<ID=LPL,Number=.,Type=Integer,Description=\"Localized phred-scaled genotype likelihoods\">" },
    [LAD]  = { .type = BCF_HT_INT,  .str = "LAD", .hdr = "##FORMAT=<ID=LAD,Number=.,Type=Integer,Description=\"Localized allelic depths\">" },
    [AD]   = { .type = BCF_HT_INT,  .str = "AD",  .hdr = "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">" },
};
const int tags_XX[]  = { PL, AD, 0 };
const int tags_LXX[] = { LPL, LAD, LAA, 0 };

typedef struct
{
    enum tag src, dst;
    int loc_src, loc_dst;   // bitmask of localized tags marked for conversion (eg 1<<LPL)
    int drop_src, skip_nalt;
    bcf_hdr_t *in_hdr, *out_hdr;
    float *farr, gp_th;
    int32_t dflt_AD, dflt_PL;
    int32_t *iarr, *iarr2, *iarr3, *iarr4;
    int mfarr, miarr, miarr2, miarr3, miarr4;
}
args_t;

static args_t *args;

const char *about(void)
{
    return "Convert between similar tags, such as GL,PL,GP or QR,QA,QS or localized alleles, eg PL,LPL.\n";
}

const char *usage(void)
{
    return
        "\n"
        "About: Convert between similar tags such as GL,PL,GP or QR,QA,QS or localized alleles, eg PL,LPL.\n"
        "Usage: bcftools +tag2tag [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "       --ORI-to-NEW           Convert from source tags FORMAT/GL,PL,GP to FORMAT/GL,PL,GP,GT\n"
        "       --QR-QA-to-QS          Convert FORMAT/QR,QA to FORMAT/QS\n"
        "       --XX-to-LXX            Convert from normal to localized tags (from PL,AD to LAA,LPL,LAD)\n"
        "       --LXX-to-XX            Convert from localized to normal tags (from LAA,LPL,LAD to PL,AD)\n"
        "   -d, --defaults LIST        Values to use in place of missing --LXX-to-XX fields [AD:.,PL:.]\n"
        "   -s, --skip-nalt INT        Do not modify sites with fewer (--XX-to-LXX) or more (--LXX-to-XX)\n"
        "                              than INT alternate alleles, 0=all sites [0]\n"
        "   -r, --replace              Drop the source tag\n"
        "   -t, --threshold FLOAT      Threshold for GP to GT hard-call [0.1]\n"
        "\n"
        "Examples:\n"
        "   bcftools +tag2tag in.vcf -- -r --GP-to-GL\n"
        "   bcftools +tag2tag in.vcf --    --PL-to-GT\n"
        "\n"
        "   # Expand the family of localized tags (LPL,LAD,LAF) to normal tags (PL,AD,AF). With the -s 3\n"
        "   # option, only sites with fewer than three alternate alleles will be expanded.\n"
        "   # Note that FORMAT/LAA must be present.\n"
        "   bcftools +tag2tag in.vcf -- --LXX-to-XX -s 3\n"
        "\n"
        "   # The same as above, but convert only the LAD tag\n"
        "   bcftools +tag2tag in.vcf -- --LAD-to-AD -s 3\n"
        "\n"
        "   # Replace the family of normal tags (PL,AD,AF) with localized tags (LPL,LAD,LAF). With the -s 3\n"
        "   # option, only sites with more than three alternate alleles will be modified. Note that in order\n"
        "   # to remove the original tags, -r must be given.\n"
        "   bcftools +tag2tag in.vcf -- --XX-to-LXX -s 3 -r\n"
        "\n"
        "   # The same as above, but convert only the PL tag\n"
        "   bcftools +tag2tag in.vcf -- --PL-to-LPL -s 3\n"
        "\n";
}

static int parse_ori2new_option(args_t *args, char *optarg)
{
    if ( args->src!=UNKN && args->src!=XX && args->src!=LXX ) error("Multiple modes selected: %s\n",optarg);
    if ( strncmp("--",optarg,2) ) return -1;

    if ( !strcasecmp("--LXX-to-XX",optarg) )
    {
        args->src = LXX;
        args->dst = XX;
        args->loc_src = (1<<LPL)|(1<<LAD)|(1<<LAA);
        args->loc_dst = (1<<PL)|(1<<AD);
    }
    else if ( !strcasecmp("--LPL-to-PL",optarg) )
    {
        args->src = LXX;
        args->dst = XX;
        args->loc_src = (1<<LPL)|(1<<LAA);
        args->loc_dst = (1<<PL);
    }
    else if ( !strcasecmp("--LAD-to-AD",optarg) )
    {
        args->src = LXX;
        args->dst = XX;
        args->loc_src = (1<<LAD)|(1<<LAA);
        args->loc_dst = (1<<AD);
    }
    else if ( !strcasecmp("--XX-to-LXX",optarg) )
    {
        args->src = XX;
        args->dst = LXX;
        args->loc_src = (1<<PL)|(1<<AD);
        args->loc_dst = (1<<LPL)|(1<<LAD)|(1<<LAA);
    }
    else if ( !strcasecmp("--PL-to-LPL",optarg) )
    {
        args->src = XX;
        args->dst = LXX;
        args->loc_src = (1<<PL);
        args->loc_dst = (1<<LPL)|(1<<LAA);
    }
    else if ( !strcasecmp("--AD-to-LAD",optarg) )
    {
        args->src = XX;
        args->dst = LXX;
        args->loc_src = (1<<AD);
        args->loc_dst = (1<<LAD)|(1<<LAA);
    }
    else if ( !strcasecmp("--QR-QA-to-QS",optarg) )
    {
        args->src = QRQA;
        args->dst = QS;
    }
    else
    {
        char *tmp = optarg+2;
        if ( !strncasecmp("GL-to-",tmp,6) ) args->src = GL, tmp += 6;
        else if ( !strncasecmp("PL-to-",tmp,6) ) args->src = PL, tmp += 6;
        else if ( !strncasecmp("GP-to-",tmp,6) ) args->src = GP, tmp += 6;
        if ( args->src==UNKN ) error("The conversion is not supported: %s\n",optarg);

        if ( !strcasecmp("GL",tmp) ) args->dst = GL;
        else if ( !strcasecmp("PL",tmp) ) args->dst = PL;
        else if ( !strcasecmp("GP",tmp) ) args->dst = GP;
        else if ( !strcasecmp("GT",tmp) ) args->dst = GT;
        if ( args->dst==UNKN || args->src==args->dst ) error("The conversion is not supported: %s\n",optarg);
    }

    return 0;
}

void parse_defaults(args_t *args, char *optarg)
{
    char *tmp, *ptr = optarg;
    while ( *ptr )
    {
        if ( !strncasecmp("AD:",ptr,3) )
        {
            ptr += 3;
            args->dflt_AD = strtol(ptr,&tmp,10);
            if ( *tmp && *tmp!=',' && *tmp!='.' ) error("Could not parse: --defaults %s\n",optarg);
            if ( *tmp=='.' ) { args->dflt_AD = bcf_int32_missing; tmp++; }
            if ( *tmp ) tmp++;
            ptr = tmp;
        }
        else if ( !strncasecmp("PL:",ptr,3) )
        {
            ptr += 3;
            args->dflt_PL = strtol(ptr,&tmp,10);
            if ( *tmp && *tmp!=',' && *tmp!='.' ) error("Could not parse: --defaults %s\n",optarg);
            if ( *tmp=='.' ) { args->dflt_PL = bcf_int32_missing; tmp++; }
            if ( *tmp ) tmp++;
            ptr = tmp;
        }
    }
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t*) calloc(1,sizeof(args_t));
    args->in_hdr  = in;
    args->out_hdr = out;
    args->src = args->dst = UNKN;
    args->dflt_AD = bcf_int32_missing;
    args->dflt_PL = bcf_int32_missing;

    static struct option loptions[] =
    {
        {"replace",no_argument,NULL,'r'},
        {"threshold",required_argument,NULL,'t'},
        {"skip-nalt",required_argument,NULL,'s'},
        {"defaults",required_argument,NULL,'d'},
        {NULL,0,NULL,0}
    };
    opterr = 0;
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "?hrt:s:d:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'r': args->drop_src = 1; break;
            case 'd': parse_defaults(args,optarg); break;
            case 's':
                args->skip_nalt = strtol(optarg,&tmp,10);
                if ( *tmp || args->skip_nalt<0 ) error("Could not parse: --skip-nalt %s\n", optarg);
                break;
            case 't':
                args->gp_th = strtod(optarg,&tmp);
                if ( *tmp || args->gp_th<0 || args->gp_th>1 ) error("Expected value between [0-1] for -t, --threshold, got %s\n",optarg);
                break;
            case 'h': error("%s", usage()); break;
            default:
                if ( optind<=0 || parse_ori2new_option(args,argv[optind-1])<0 ) error("%s", usage());
                break;
        }
    }
    if ( args->src==UNKN ) error("%s", usage());

    // Check if the source tags are defined in the header
    int tag_id,i;
    if ( args->src==QRQA )
    {
        if ( (tag_id=bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"QR"))<0 || !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The source tag does not exist: QR\n");
        if ( bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id)!=tags[args->src].type ) error("The source tag is of different type than required by the VCF specification\n");
        if ( (tag_id=bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"QA"))<0 || !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The source tag does not exist: QA\n");
        if ( bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id)!=tags[args->src].type ) error("The source tag is of different type than required by the VCF specification\n");
    }
    else if ( args->src==LXX )
    {
        for (i=0; tags_LXX[i]; i++)
        {
            int j = tags_LXX[i];
            if ( !(args->loc_src & (1<<j)) ) continue;
            if ( (tag_id=bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,tags[j].str))<0 || !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The source tag does not exist: %s\n",tags[j].str);
            if ( bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id)!=tags[j].type )
                error("The source tag %s is of different type than required by the VCF specification (%d vs %d)\n",tags[j].str,bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id),tags[j].type);
        }
    }
    else if ( args->src==XX )
    {
        for (i=0; tags_XX[i]; i++)
        {
            int j = tags_XX[i];
            if ( !(args->loc_src & (1<<j)) ) continue;
            if ( (tag_id=bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,tags[j].str))<0 || !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The source tag does not exist: %s\n",tags[j].str);
            if ( bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id)!=tags[j].type )
                error("The source tag %s is of different type than required by the VCF specification (%d vs %d)\n",tags[j].str,bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id),tags[j].type);
        }
    }
    else
    {
        if ( (tag_id=bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,tags[args->src].str))<0 || !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The source tag does not exist: %s\n",tags[args->src].str);
        if ( bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id)!=tags[args->src].type )
            error("The source tag %s is of different type than required by the VCF specification (%d vs %d)\n",tags[args->src].str,bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id),tags[args->src].type);
    }

    // Remove tags from the header if -r, --replace was given. However, do not remove if -s, --skip-nalt was given,
    // in that case some of the records may retain the tags.
    if ( args->drop_src )
    {
        if ( args->src==QRQA )
        {
            bcf_hdr_remove(args->out_hdr,BCF_HL_FMT,"QR");
            bcf_hdr_remove(args->out_hdr,BCF_HL_FMT,"QA");
        }
        else if ( args->src==LXX )
        {
            for (i=0; tags_LXX[i]; i++)
            {
                int j = tags_LXX[i];
                if ( args->skip_nalt || !(args->loc_src & (1<<j)) ) continue;
                bcf_hdr_remove(args->out_hdr,BCF_HL_FMT,tags[j].str);
            }
        }
        else if ( args->src==XX )
        {
            for (i=0; tags_XX[i]; i++)
            {
                int j = tags_XX[i];
                if ( args->skip_nalt || !(args->loc_src & (1<<j)) ) continue;
                bcf_hdr_remove(args->out_hdr,BCF_HL_FMT,tags[j].str);
            }
        }
        else
            bcf_hdr_remove(args->out_hdr,BCF_HL_FMT,tags[args->src].str);
    }

    // Add the new tags to the header
    if ( tags[args->dst].hdr )
        bcf_hdr_append(args->out_hdr,tags[args->dst].hdr);
    else if ( args->dst==LXX )
    {
        for (i=0; tags_LXX[i]; i++)
        {
            int j = tags_LXX[i];
            if ( !(args->loc_dst & (1<<j)) ) continue;
            bcf_hdr_append(args->out_hdr,tags[j].hdr);
        }
    }
    else if ( args->dst==XX )
    {
        for (i=0; tags_XX[i]; i++)
        {
            int j = tags_XX[i];
            if ( !(args->loc_dst & (1<<j)) ) continue;
            bcf_hdr_append(args->out_hdr,tags[j].hdr);
        }
    }

    return 0;
}

bcf1_t *process_LXX(bcf1_t *rec)
{
    if ( args->skip_nalt && rec->n_allele > args->skip_nalt ) return rec;

    // Use LAA to determine mapping between LPL->PL and LAD->AD indices for each sample
    int nlaa = bcf_get_format_int32(args->in_hdr,rec,"LAA",&args->iarr,&args->miarr);
    if ( nlaa <=0 ) return rec;

    int nals  = rec->n_allele;
    int nsmpl = bcf_hdr_nsamples(args->in_hdr);
    int i,j,k, nlaa1 = nlaa / nsmpl;
    int nsrc;
    int drop_laa = args->loc_src;

    if ( (args->loc_src & (1<<LAD)) && (nsrc=bcf_get_format_int32(args->in_hdr,rec,"LAD",&args->iarr2,&args->miarr2))>0 )
    {
        int nsrc1 = nsrc / nsmpl;
        int ndst1 = nals;
        if ( hts_resize(int32_t,nsmpl*ndst1,&args->miarr3,&args->iarr3,0)!=0 ) return rec;   // too big
        for (i=0; i<nsmpl; i++)
        {
            // This code is agnostic to the allele ordering in LAA, although VCF specification requires
            // them to be in ascending order. This is violated by some programs, eg DRAGEN.
            int32_t *laa = args->iarr  + i*nlaa1;
            int32_t *dst = args->iarr3 + i*ndst1;
            int32_t *src = args->iarr2 + i*nsrc1;
            dst[0] = src[0];
            for (j=1; j<nals; j++) dst[j] = args->dflt_AD;
            for (j=1; j<nsrc1; j++)
                if ( laa[j-1]>=0 && laa[j-1]<nals ) dst[laa[j-1]] = src[j];
        }
        if ( bcf_update_format_int32(args->out_hdr,rec,"AD",args->iarr3,nsmpl*nals)!=0 )
            error("Error: Failed to set FORMAT/AD at %s:%"PRIhts_pos"\n",bcf_seqname(args->in_hdr,rec),rec->pos+1);

        if ( args->drop_src )
        {
            bcf_update_format_int32(args->out_hdr,rec,"LAD",NULL,0);
            drop_laa &= ~(1<<LAD);
        }
    }

    if ( (args->loc_src & (1<<LPL)) && (nsrc=bcf_get_format_int32(args->in_hdr,rec,"LPL",&args->iarr2,&args->miarr2))>0 )
    {
        // for convenience to keep LAA alleles which don't include 0
        hts_expand(int32_t,nlaa1+1,args->miarr4,args->iarr4);
        int32_t *tmp_laa = args->iarr4;
        tmp_laa[0] = 0;

        int nsrc1 = nsrc / nsmpl;
        int ndst1 = nals*(nals+1)/2;
        if ( hts_resize(int32_t,nsmpl*ndst1,&args->miarr3,&args->iarr3,0)!=0 ) return rec;   // too big
        for (i=0; i<nsmpl; i++)
        {
            int32_t *laa = args->iarr  + i*nlaa1;
            int32_t *dst = args->iarr3 + i*ndst1;
            int32_t *src = args->iarr2 + i*nsrc1;

            // assuming the simple case: diploid
            for (j=0; j<ndst1; j++) dst[j] = args->dflt_PL;
            for (j=0; j<nlaa1; j++) tmp_laa[j+1] = laa[j];
            for (j=0; j<=nlaa1 && tmp_laa[j]>=0 && tmp_laa[j]<nals; j++)
                for (k=0; k<=j; k++)
                {
                    int idx = tmp_laa[j]*(tmp_laa[j]+1)/2 + tmp_laa[k];
                    dst[idx] = *src;
                    src++;
                }

        }
        bcf_update_format_int32(args->out_hdr,rec,"PL",args->iarr3,nsmpl*ndst1);
        if ( args->drop_src )
        {
            bcf_update_format_int32(args->out_hdr,rec,"LPL",NULL,0);
            drop_laa &= ~(1<<LPL);
        }
    }

    if ( args->drop_src && drop_laa==1<<LAA )
        bcf_update_format_int32(args->out_hdr,rec,"LAA",NULL,0);

    return rec;
}

bcf1_t *process_XX(bcf1_t *rec)
{
    error("todo: --XX-to-LXX\n");
    if ( args->skip_nalt && rec->n_allele < args->skip_nalt ) return rec;
    return rec;
}

bcf1_t *process(bcf1_t *rec)
{
    int i,j,n;

    if ( args->src==LXX ) return process_LXX(rec);
    if ( args->src==XX ) return process_XX(rec);

    if ( args->src==QRQA )
    {
        int nals  = rec->n_allele;
        int nsmpl = bcf_hdr_nsamples(args->in_hdr);

        n = bcf_get_format_int32(args->in_hdr,rec,"QR",&args->iarr,&args->miarr);
        if ( n<=0 ) return rec;
        if ( n!=nsmpl ) error("Unexpected number of QR values at %s:%"PRIhts_pos"\n",bcf_seqname(args->in_hdr,rec),rec->pos+1);
        if ( nals==1 )
            bcf_update_format_int32(args->out_hdr,rec,"QS",args->iarr,nsmpl);
        else
        {
            int n2 = bcf_get_format_int32(args->in_hdr,rec,"QA",&args->iarr2,&args->miarr2);
            if ( n2<=0 ) return rec;
            if ( n*(nals-1) != n2 ) error("Unexpected number of QR vs QA values at %s:%"PRIhts_pos"\n",bcf_seqname(args->in_hdr,rec),rec->pos+1);
            hts_expand(int32_t,nsmpl*nals,args->miarr3,args->iarr3);
            for (i=0; i<n; i++)
            {
                args->iarr3[i*nals] = args->iarr[i];
                for (j=1; j<nals; j++) args->iarr3[i*nals+j] = args->iarr2[i*(nals-1)+j-1];
            }
            bcf_update_format_int32(args->out_hdr,rec,"QS",args->iarr3,nals*nsmpl);
        }
        if ( args->drop_src )
        {
            bcf_update_format_int32(args->out_hdr,rec,"QR",NULL,0);
            bcf_update_format_int32(args->out_hdr,rec,"QA",NULL,0);
        }
        return rec;
    }

    // convert source tags to GL. If performance is an issue, can be smarter
    if ( tags[args->src].type==BCF_HT_REAL )
    {
        // from GL,GP to something
        n = bcf_get_format_float(args->in_hdr,rec,tags[args->src].str,&args->farr,&args->mfarr);
        if ( n<=0 ) return rec;

        // convert to GL
        if ( args->src==GP )
        {
            for (i=0; i<n; i++)
            {
                if ( bcf_float_is_missing(args->farr[i]) || bcf_float_is_vector_end(args->farr[i]) ) continue;
                args->farr[i] = args->farr[i] ? log10(args->farr[i]) : -99;
            }
        }
    }
    else if ( tags[args->src].type==BCF_HT_INT )
    {
        // from PL to something
        n = bcf_get_format_int32(args->in_hdr,rec,tags[args->src].str,&args->iarr,&args->miarr);
        if ( n<=0 ) return rec;
        hts_expand(float, n, args->mfarr, args->farr);

        // convert to GL
        for (i=0; i<n; i++)
        {
            if ( args->iarr[i]==bcf_int32_missing )
                bcf_float_set_missing(args->farr[i]);
            else if ( args->iarr[i]==bcf_int32_vector_end )
                bcf_float_set_vector_end(args->farr[i]);
            else
                args->farr[i] = -0.1 * args->iarr[i];
        }
    }
    else error("fixme: expected BCF_HT_REAL or BCF_HT_INT for source tag\n");

    if ( args->dst==GL )
        bcf_update_format_float(args->out_hdr,rec,"GL",args->farr,n);
    else if ( args->dst==PL )
    {
        hts_expand(int32_t, n, args->miarr, args->iarr);
        for (i=0; i<n; i++)
        {
            if ( bcf_float_is_missing(args->farr[i]) ) args->iarr[i] = bcf_int32_missing;
            else if ( bcf_float_is_vector_end(args->farr[i]) ) args->iarr[i] = bcf_int32_vector_end;
            else args->iarr[i] = lroundf(-10*args->farr[i]);
        }
        bcf_update_format_int32(args->out_hdr,rec,"PL",args->iarr,n);
    }
    else if ( args->dst==GP || args->dst==GT )
    {
        int nsmpl = bcf_hdr_nsamples(args->in_hdr);
        int n1 = n / nsmpl;
        for (i=0; i<nsmpl; i++)
        {
            float *ptr = args->farr + i*n1;
            float sum = 0;
            for (j=0; j<n1; j++)
            {
                if ( bcf_float_is_missing(ptr[j]) ) continue;
                if ( bcf_float_is_vector_end(ptr[j]) ) break;
                ptr[j] = pow(10, ptr[j]);
                sum += ptr[j];
            }
            if ( sum<=0 ) continue;
            for (j=0; j<n1; j++)
            {
                if ( bcf_float_is_missing(ptr[j]) ) continue;
                if ( bcf_float_is_vector_end(ptr[j]) ) break;
                ptr[j] /= sum;
            }
        }
        if ( args->dst==GP )
            bcf_update_format_float(args->out_hdr,rec,"GP",args->farr,n);
    }
    if ( args->dst==GT )
    {
        // farr contains GP. If performance is an issue, note that if the src was GP, we made an unnecessary conversion GP->GL->GP
        int nsmpl = bcf_hdr_nsamples(args->in_hdr);
        int n1   = n / nsmpl;
        int nals = rec->n_allele;
        hts_expand(int32_t,nsmpl*2,args->miarr,args->iarr);
        for (i=0; i<nsmpl; i++)
        {
            float *ptr = args->farr + i*n1;
            if ( bcf_float_is_missing(ptr[0]) )
            {
                args->iarr[2*i] = args->iarr[2*i+1] = bcf_gt_missing;
                continue;
            }
            int jmax = 0;
            for (j=1; j<n1; j++)
            {
                if ( bcf_float_is_missing(ptr[j]) || bcf_float_is_vector_end(ptr[j]) ) break;
                if ( ptr[j] > ptr[jmax] ) jmax = j;
            }

            // haploid genotype
            if ( j==nals )
            {
                args->iarr[2*i]   = ptr[jmax] < 1 - args->gp_th ? bcf_gt_missing : bcf_gt_unphased(jmax);
                args->iarr[2*i+1] = bcf_int32_vector_end;
                continue;
            }

            if ( j!=nals*(nals+1)/2 )
                error("Wrong number of GP values for diploid genotype at %s:%"PRId64", expected %d, found %d\n",
                        bcf_seqname(args->in_hdr,rec), rec->pos+1, nals*(nals+1)/2,j);

            if ( ptr[jmax] < 1 - args->gp_th )
            {
                args->iarr[2*i] = args->iarr[2*i+1] = bcf_gt_missing;
                continue;
            }

            // most common case: RR
            if ( jmax==0 )
            {
                args->iarr[2*i] = args->iarr[2*i+1] = bcf_gt_unphased(0);
                continue;
            }

            int a,b;
            bcf_gt2alleles(jmax,&a,&b);
            args->iarr[2*i]   = bcf_gt_unphased(a);
            args->iarr[2*i+1] = bcf_gt_unphased(b);
        }
        bcf_update_genotypes(args->out_hdr,rec,args->iarr,nsmpl*2);
    }

    if ( args->drop_src )
    {
        if ( tags[args->src].type==BCF_HT_REAL )
            bcf_update_format_float(args->out_hdr,rec,tags[args->src].str,NULL,0);
        else if ( tags[args->src].type==BCF_HT_INT )
            bcf_update_format_int32(args->out_hdr,rec,tags[args->src].str,NULL,0);
    }
    return rec;
}

void destroy(void)
{
    free(args->farr);
    free(args->iarr);
    free(args->iarr2);
    free(args->iarr3);
    free(args->iarr4);
    free(args);
}


