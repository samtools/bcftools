/*  plugins/tag2tag.c -- convert between similar tags

    Copyright (C) 2014-2023 Genome Research Ltd.

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
enum tag { UNKN, GP, GL, PL, GT, QRQA, QS };
tag_t tags[] =
{
    [UNKN] = { .type = 0 ,          .str = NULL, .hdr = NULL },
    [GP]   = { .type = BCF_HT_REAL, .str = "GP", .hdr = "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype probabilities\">" },
    [GL]   = { .type = BCF_HT_REAL, .str = "GL", .hdr = "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods\">" },
    [PL]   = { .type = BCF_HT_INT,  .str = "PL", .hdr = "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">" },
    [GT]   = { .type = BCF_HT_STR,  .str = "GT", .hdr = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" },
    [QRQA] = { .type = BCF_HT_INT,  .str = NULL, .hdr = NULL },
    [QS]   = { .type = BCF_HT_INT,  .str = "QS", .hdr = "##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Phred-score allele quality sum\">" },
};

typedef struct
{
    enum tag src, dst;
    int drop_src;
    bcf_hdr_t *in_hdr, *out_hdr;
    float *farr, gp_th;
    int32_t *iarr, *iarr2, *iarr3;
    int mfarr, miarr, miarr2, miarr3;
}
args_t;

static args_t *args;

const char *about(void)
{
    return "Convert between similar tags, such as GL,PL,GP or QR,QA,QS.\n";
}

const char *usage(void)
{
    return
        "\n"
        "About: Convert between similar tags such as GL,PL,GP or QR,QA,QS.\n"
        "Usage: bcftools +tag2tag [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "       --ORI-to-NEW           Convert from source tags FORMAT/GL,PL,GP to FORMAT/GL,PL,GP,GT\n"
        "       --QR-QA-to-QS          Convert FORMAT/QR,QA to FORMAT/QS\n"
        "   -r, --replace              Drop the source tag\n"
        "   -t, --threshold FLOAT      Threshold for GP to GT hard-call [0.1]\n"
        "\n"
        "Example:\n"
        "   bcftools +tag2tag in.vcf -- -r --GP-to-GL\n"
        "   bcftools +tag2tag in.vcf --    --PL-to-GT\n"
        "\n";
}

static int parse_ori2new_option(args_t *args, char *optarg)
{
    if ( args->src!=UNKN ) error("Multiple modes selected: %s\n",optarg);
    if ( strncmp("--",optarg,2) ) return -1;

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
    return 0;
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t*) calloc(1,sizeof(args_t));
    args->in_hdr  = in;
    args->out_hdr = out;
    args->src = args->dst = UNKN;

    static struct option loptions[] =
    {
        {"replace",no_argument,NULL,'r'},
        {"qr-qa-to-qs",no_argument,NULL,1},
        {"QR-QA-to-QS",no_argument,NULL,1},
        {"threshold",required_argument,NULL,'t'},
        {NULL,0,NULL,0}
    };
    opterr = 0;
    int c;
    while ((c = getopt_long(argc, argv, "?hrt:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->src = QRQA; args->dst = QS; break;
            case 'r': args->drop_src = 1; break;
            case 't':
                args->gp_th = atof(optarg);
                if ( args->gp_th<0 || args->gp_th>1 ) error("Expected value between [0-1] for -t, --threshold, got %f\n",args->gp_th);
                break;
            case 'h': error("%s", usage()); break;
            default:
                if ( optind<=0 || parse_ori2new_option(args,argv[optind-1])<0 ) error("%s", usage());
                break;
        }
    }
    if ( args->src==UNKN ) error("%s", usage());
    int tag_id;
    if ( args->src==QRQA )
    {
        if ( (tag_id=bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"QR"))<0 || !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The source tag does not exist: QR\n");
        if ( bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id)!=tags[args->src].type ) error("The source tag is of different type than required by the VCF specification\n");
        if ( (tag_id=bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"QA"))<0 || !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The source tag does not exist: QA\n");
        if ( bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id)!=tags[args->src].type ) error("The source tag is of different type than required by the VCF specification\n");
    }
    else
    {
        if ( (tag_id=bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,tags[args->src].str))<0 || !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The source tag does not exist: %s\n",tags[args->src].str);
        if ( bcf_hdr_id2type(args->in_hdr,BCF_HL_FMT,tag_id)!=tags[args->src].type ) error("The source tag is of different type than required by the VCF specification\n");
    }
    if ( args->drop_src )
    {
        if ( args->src==QRQA )
        {
            bcf_hdr_remove(args->out_hdr,BCF_HL_FMT,"QR");
            bcf_hdr_remove(args->out_hdr,BCF_HL_FMT,"QA");
        }
        else
            bcf_hdr_remove(args->out_hdr,BCF_HL_FMT,tags[args->src].str);
    }
    assert( tags[args->dst].hdr );
    bcf_hdr_append(args->out_hdr,tags[args->dst].hdr);

    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    int i,j,n;
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
    free(args);
}


