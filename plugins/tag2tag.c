/*  plugins/tag2tag.c -- convert between similar tags

    Copyright (C) 2014-2021 Genome Research Ltd.

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


#define GP_TO_GL 1
#define GL_TO_PL 2
#define GP_TO_GT 3
#define PL_TO_GL 4
#define QRQA_TO_QS 5

static int mode = 0, drop_source_tag = 0;
static bcf_hdr_t *in_hdr, *out_hdr;
static float *farr = NULL, thresh = 0.1;
static int32_t *iarr = NULL;
static int32_t *iarr2 = NULL;
static int32_t *iarr3 = NULL;
static int mfarr = 0, miarr = 0, miarr2 = 0, miarr3 = 0;

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
        "       --GP-to-GL           convert FORMAT/GP to FORMAT/GL\n"
        "       --GP-to-GT           convert FORMAT/GP to FORMAT/GT by taking argmax of GP\n"
        "       --GL-to-PL           convert FORMAT/GL to FORMAT/PL\n"
        "       --PL-to-GL           convert FORMAT/PL to FORMAT/GL\n"
        "       --QR-QA-to-QS        convert FORMAT/QR,QA to FORMAT/QS\n"
        "   -r, --replace            drop the source tag\n"
        "   -t, --threshold <float>  threshold for GP to GT hard-call [0.1]\n"
        "\n"
        "Example:\n"
        "   bcftools +tag2tag in.vcf -- -r --gp-to-gl\n"
        "\n";
}


static void init_header(bcf_hdr_t *hdr, const char *ori, int ori_type, const char *new_hdr_line)
{
    if ( ori )
        bcf_hdr_remove(hdr,ori_type,ori);

    bcf_hdr_append(hdr, new_hdr_line);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    static struct option loptions[] =
    {
        {"replace",no_argument,NULL,'r'},
        {"gp-to-gl",no_argument,NULL,1},
        {"GP-to-GL",no_argument,NULL,1},
        {"gl-to-pl",no_argument,NULL,2},
        {"GL-to-PL",no_argument,NULL,2},
        {"gp-to-gt",no_argument,NULL,3},
        {"GP-to-GT",no_argument,NULL,3},
        {"pl-to-gl",no_argument,NULL,4},
        {"PL-to-GL",no_argument,NULL,4},
        {"qr-qa-to-qs",no_argument,NULL,5},
        {"QR-QA-to-QS",no_argument,NULL,5},
        {"threshold",required_argument,NULL,'t'},
        {NULL,0,NULL,0}
    };
    int c;
    char *src_tag = "GP";
    int src_type  = BCF_HT_REAL;
    while ((c = getopt_long(argc, argv, "?hrt:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case  1 : src_tag = "GP"; mode = GP_TO_GL; src_type = BCF_HT_REAL; break;
            case  2 : src_tag = "GL"; mode = GL_TO_PL; src_type = BCF_HT_REAL; break;
            case  3 : src_tag = "GP"; mode = GP_TO_GT; src_type = BCF_HT_REAL; break;
            case  4 : src_tag = "PL"; mode = PL_TO_GL; src_type = BCF_HT_INT; break;
            case  5 : src_tag = "QR"; mode = QRQA_TO_QS; src_type = BCF_HT_INT; break;
            case 'r': drop_source_tag = 1; break;
            case 't': thresh = atof(optarg); break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !mode ) mode = GP_TO_GL;

    in_hdr  = in;
    out_hdr = out;

    if ( mode==GP_TO_GL )
        init_header(out_hdr,drop_source_tag?"GP":NULL,BCF_HL_FMT,"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihoods\">");
    else if ( mode==GL_TO_PL )
        init_header(out_hdr,drop_source_tag?"GL":NULL,BCF_HL_FMT,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">");
    else if ( mode==PL_TO_GL )
        init_header(out_hdr,drop_source_tag?"PL":NULL,BCF_HL_FMT,"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods\">");
    else if ( mode==GP_TO_GT )
    {
        if (thresh<0||thresh>1) error("--threshold must be in the range [0,1]: %f\n", thresh);
        init_header(out_hdr,drop_source_tag?"GP":NULL,BCF_HL_FMT,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    }
    else if ( mode==QRQA_TO_QS )
    {
        if ( drop_source_tag )
        {
            bcf_hdr_remove(out_hdr,BCF_HL_FMT,"QR");
            bcf_hdr_remove(out_hdr,BCF_HL_FMT,"QA");
        }
        bcf_hdr_append(out_hdr, "##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Phred-score allele quality sum\">");
    }

    int tag_id;
    if ( (tag_id=bcf_hdr_id2int(in_hdr,BCF_DT_ID,src_tag))<0 || !bcf_hdr_idinfo_exists(in_hdr,BCF_HL_FMT,tag_id) )
        error("The source tag does not exist: %s\n", src_tag);
    if ( bcf_hdr_id2type(in_hdr,BCF_HL_FMT,tag_id) != src_type )
        error("The source tag type does not match the VCF specification, expected Type=%s. Use `bcftools reheader` to fix.\n",src_type==BCF_HT_REAL?"Float":"Integer");
    if ( mode==QRQA_TO_QS )
    {
        src_tag = "QA";
        if ( (tag_id=bcf_hdr_id2int(in_hdr,BCF_DT_ID,src_tag))<0 || !bcf_hdr_idinfo_exists(in_hdr,BCF_HL_FMT,tag_id) )
            error("The source tag does not exist: %s\n", src_tag);
        if ( bcf_hdr_id2type(in_hdr,BCF_HL_FMT,tag_id) != src_type )
            error("The source tag type does not match the VCF specification, expected Type=%s. Use `bcftools reheader` to fix.\n",src_type==BCF_HT_REAL?"Float":"Integer");
    }

    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    int i, n;
    if ( mode==GP_TO_GL )
    {
        n = bcf_get_format_float(in_hdr,rec,"GP",&farr,&mfarr);
        if ( n<=0 ) return rec;
        for (i=0; i<n; i++)
        {
            if ( bcf_float_is_missing(farr[i]) || bcf_float_is_vector_end(farr[i]) ) continue;
            farr[i] = farr[i] ? log10(farr[i]) : -99;
        }
        bcf_update_format_float(out_hdr,rec,"GL",farr,n);
        if ( drop_source_tag )
            bcf_update_format_float(out_hdr,rec,"GP",NULL,0);
    }
    else if ( mode==PL_TO_GL )
    {
        n = bcf_get_format_int32(in_hdr,rec,"PL",&iarr,&miarr);
        if ( n<=0 ) return rec;
        hts_expand(float, n, mfarr, farr);
        for (i=0; i<n; i++)
        {
            if ( iarr[i]==bcf_int32_missing )
                bcf_float_set_missing(farr[i]);
            else if ( iarr[i]==bcf_int32_vector_end )
                bcf_float_set_vector_end(farr[i]);
            else
                farr[i] = -0.1 * iarr[i];
        }
        bcf_update_format_float(out_hdr,rec,"GL",farr,n);
        if ( drop_source_tag )
            bcf_update_format_int32(out_hdr,rec,"PL",NULL,0);
    }
    else if ( mode==GL_TO_PL )
    {
        n = bcf_get_format_float(in_hdr,rec,"GL",&farr,&mfarr);
        if ( n<=0 ) return rec;
        hts_expand(int32_t, n, miarr, iarr);
        for (i=0; i<n; i++)
        {
            if ( bcf_float_is_missing(farr[i]) )
                iarr[i] = bcf_int32_missing;
            else if ( bcf_float_is_vector_end(farr[i]) )
                iarr[i] = bcf_int32_vector_end;
            else
                iarr[i] = lroundf(-10 * farr[i]);
        }
        bcf_update_format_int32(out_hdr,rec,"PL",iarr,n);
        if ( drop_source_tag )
            bcf_update_format_float(out_hdr,rec,"GL",NULL,0);
    }
    else if ( mode==GP_TO_GT )
    {
        int nals  = rec->n_allele;
        int nsmpl = bcf_hdr_nsamples(in_hdr);
        hts_expand(int32_t,nsmpl*2,miarr,iarr);

        n = bcf_get_format_float(in_hdr,rec,"GP",&farr,&mfarr);
        if ( n<=0 ) return rec;

        n /= nsmpl;
        for (i=0; i<nsmpl; i++)
        {
            float *ptr = farr + i*n;
            if ( bcf_float_is_missing(ptr[0]) )
            {
                iarr[2*i] = iarr[2*i+1] = bcf_gt_missing;
                continue;
            }

            int j, jmax = 0;
            for (j=1; j<n; j++)
            {
                if ( bcf_float_is_missing(ptr[j]) || bcf_float_is_vector_end(ptr[j]) ) break;
                if ( ptr[j] > ptr[jmax] ) jmax = j;
            }

            // haploid genotype
            if ( j==nals )
            {
                iarr[2*i]   = ptr[jmax] < 1-thresh ? bcf_gt_missing : bcf_gt_unphased(jmax);
                iarr[2*i+1] = bcf_int32_vector_end;
                continue;
            }

            if ( j!=nals*(nals+1)/2 )
                error("Wrong number of GP values for diploid genotype at %s:%"PRId64", expected %d, found %d\n",
                    bcf_seqname(in_hdr,rec),(int64_t) rec->pos+1, nals*(nals+1)/2,j);

            if (ptr[jmax] < 1-thresh)
            {
                iarr[2*i] = iarr[2*i+1] = bcf_gt_missing;
                continue;
            }

            // most common case: RR
            if ( jmax==0 )
            {
                iarr[2*i] = iarr[2*i+1] = bcf_gt_unphased(0);
                continue;
            }

            int a,b;
            bcf_gt2alleles(jmax,&a,&b);
            iarr[2*i]   = bcf_gt_unphased(a);
            iarr[2*i+1] = bcf_gt_unphased(b);
        }
        bcf_update_genotypes(out_hdr,rec,iarr,nsmpl*2);
        if ( drop_source_tag )
            bcf_update_format_float(out_hdr,rec,"GP",NULL,0);
    }
    else if ( mode==QRQA_TO_QS )
    {
        int nals  = rec->n_allele;
        int nsmpl = bcf_hdr_nsamples(in_hdr);

        n = bcf_get_format_int32(in_hdr,rec,"QR",&iarr,&miarr);
        if ( n<=0 ) return rec;
        if ( n!=nsmpl ) error("Unexpected number of QR values at %s:%"PRIhts_pos"\n",bcf_seqname(in_hdr,rec),rec->pos+1);
        if ( nals==1 )
            bcf_update_format_int32(out_hdr,rec,"QS",iarr,nsmpl);
        else
        {
            int n2 = bcf_get_format_int32(in_hdr,rec,"QA",&iarr2,&miarr2);
            if ( n2<=0 ) return rec;
            if ( n*(nals-1) != n2 ) error("Unexpected number of QR vs QA values at %s:%"PRIhts_pos"\n",bcf_seqname(in_hdr,rec),rec->pos+1);
            hts_expand(int32_t,nsmpl*nals,miarr3,iarr3);
            for (i=0; i<n; i++)
            {
                int j;
                iarr3[i*nals] = iarr[i];
                for (j=1; j<nals; j++) iarr3[i*nals+j] = iarr2[i*(nals-1)+j-1];
            }
            bcf_update_format_int32(out_hdr,rec,"QS",iarr3,nals*nsmpl);
        }
        if ( drop_source_tag )
        {
            bcf_update_format_int32(out_hdr,rec,"QR",NULL,0);
            bcf_update_format_int32(out_hdr,rec,"QA",NULL,0);
        }
    }
    return rec;
}

void destroy(void)
{
    free(farr);
    free(iarr);
    free(iarr2);
    free(iarr3);
}


