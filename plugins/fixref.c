/* The MIT License

   Copyright (c) 2016 Genome Research Ltd.

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

/*
    Illumina TOP/BOT strand convention causes a lot of pain. This tool
    attempts to determine the strand convention and convert it to the
    forward reference strand.

    On TOP strand, we can encounter
        unambiguous SNPs:
            A/G
            A/C
        ambiguous (context-dependent) SNPs:
            C/G
            A/T

    On BOT strand:
        unambiguous SNPs:
            T/G
            T/C
        ambiguous (context-dependent) SNPs:
            T/A
            G/C


    For unambiguous pairs (A/C, A/G, T/C, T/G), the knowledge of reference base
    at the SNP position is enough to determine the strand:

         TOP      REF   ->  ALLELES   TOP_ON_STRAND
         -------------------------------------------
         A/C     A or C      A/C         1
          "      T or G      T/G        -1
         A/G     A or G      A/G         1
          "      T or C      T/C        -1


    For ambiguous pairs (A/T, C/G), a sequence walking must be performed
    (simultaneously upstream and downstream) until the first unambiguous pair
    is encountered. The 5' base determines the strand:

         TOP    5'REF_BASE   ->  ALLELES   TOP_ON_STRAND
         ------------------------------------------------
         A/T    A or T             A/T          1
          "     C or G             T/A         -1
         C/G    A or T             C/G          1
          "     C or G             G/C         -1

 */


#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/kfunc.h>
#include <htslib/faidx.h>
#include "bcftools.h"

#define MODE_STATS    1
#define MODE_TOP2FWD  2
#define MODE_FLIP2FWD 3

typedef struct
{
    int mode, discard;
    bcf_hdr_t *hdr;
    faidx_t *fai;
    int32_t *gts, ngts;
    uint32_t nsite,nok,nflip,nerr,nswap,nflip_swap,nonSNP,nonACGT,nonbiallelic;
    uint32_t count[4][4];
}
args_t;

args_t args;

const char *about(void)
{
    return "Fix reference strand orientation, e.g. from Illumina/TOP to fwd.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: This tool helps to determine and fix strand orientation.\n"
        "       Currently the following modes are recognised:\n"
        "           flip  .. flips non-ambiguous SNPs and ignores the rest\n"
        "           stats .. converts from Illumina TOP strand to fwd\n"
        "           top   .. converts from Illumina TOP strand to fwd\n"
        "\n"
        "       WARNING: Do not use the program blindly, make an effort to\n"
        "       understand what strand convention your data uses! Make sure\n"
        "       the reason for mismatching REF alleles is not a different\n"
        "       reference build!!\n"
        "\n"
        "Usage: bcftools +fixref [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -d, --discard               Discard sites which could not be resolved\n"
        "   -f, --fasta-ref <file.fa>   Reference sequence\n"
        "   -m, --mode <string>         Collect stats (\"stats\") or convert (\"top\", \"flip\") [stats]\n"
        "\n"
        "Examples:\n"
        "   # run stats\n"
        "   bcftools +fixref file.bcf -- -f ref.fa\n"
        "\n"
        "   # convert from TOP to fwd\n"
        "   bcftools +fixref file.bcf -Ob -o out.bcf -- -f ref.fa -m top\n"
        "\n"
        "   # assuming the reference build is correct, just flip to fwd, discarding the rest\n"
        "   bcftools +fixref file.bcf -Ob -o out.bcf -- -d -f ref.fa -m flip\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.hdr = in;
    args.mode = MODE_STATS;
    char *ref_fname = NULL;
    static struct option loptions[] =
    {
        {"mode",required_argument,NULL,'m'},
        {"discard",no_argument,NULL,'d'},
        {"fasta-ref",required_argument,NULL,'f'},
        {NULL,0,NULL,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?hf:m:d",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'm': 
                if ( !strcasecmp(optarg,"top") ) args.mode = MODE_TOP2FWD; 
                else if ( !strcasecmp(optarg,"flip") ) args.mode = MODE_FLIP2FWD; 
                else if ( !strcasecmp(optarg,"stats") ) args.mode = MODE_STATS; 
                else error("The source strand convention not recognised: %s\n", optarg);
                break;
            case 'd': args.discard = 1; break;
            case 'f': ref_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !ref_fname ) error("Expected the -f option\n");
    args.fai = fai_load(ref_fname);
    if ( !args.fai ) error("Failed to load the fai index: %s\n", ref_fname);

    if ( args.mode==MODE_STATS ) return 1;
    return 0;
}

bcf1_t *set_ref_alt(args_t *args, bcf1_t *rec, const char ref, const char alt, int swap)
{
    rec->d.allele[0][0] = ref;
    rec->d.allele[1][0] = alt;
    rec->d.shared_dirty |= BCF1_DIRTY_ALS;

    if ( !swap ) return rec;    // only fix the alleles, leaving GTs unchanged

    int ngts = bcf_get_genotypes(args->hdr, rec, &args->gts, &args->ngts);
    int i, j, nsmpl = bcf_hdr_nsamples(args->hdr);
    ngts /= nsmpl;
    for (i=0; i<nsmpl; i++)
    {
        int32_t *ptr = args->gts + i*ngts;
        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_gt_unphased(0) ) ptr[j] = bcf_gt_unphased(1);
            else if ( ptr[j]==bcf_gt_phased(0) ) ptr[j] = bcf_gt_phased(1);
            else if ( ptr[j]==bcf_gt_unphased(1) ) ptr[j] = bcf_gt_unphased(0);
            else if ( ptr[j]==bcf_gt_phased(1) ) ptr[j] = bcf_gt_phased(0);
        }
    }
    bcf_update_genotypes(args->hdr,rec,args->gts,args->ngts);
    
    return rec;
}

static inline int nt2int(char nt)
{
    nt = toupper(nt);
    if ( nt=='A' ) return 0;
    if ( nt=='C' ) return 1;
    if ( nt=='G' ) return 2;
    if ( nt=='T' ) return 3;
    return -1;
}
#define int2nt(x) "ACGT"[x]
#define revint(x) ("3210"[x]-'0')

bcf1_t *process(bcf1_t *rec)
{
    bcf1_t *ret = args.mode==MODE_STATS ? NULL : rec;
    args.nsite++;

    // Skip non-SNPs
    if ( bcf_get_variant_types(rec)!=VCF_SNP )
    {
        args.nonSNP++;
        return args.discard ? NULL : ret;
    }

    // Get the reference allele
    int len;
    char *ref = faidx_fetch_seq(args.fai, (char*)bcf_seqname(args.hdr,rec), rec->pos, rec->pos, &len);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", bcf_seqname(args.hdr,rec),rec->pos+1);
    int ir = nt2int(*ref);
    free(ref);
    if ( ir<0 )
    {
        args.nonACGT++;
        return args.discard ? NULL : ret;     // not A,C,G,T
    }

    if ( rec->n_allele!=2 )
    {
        // not a biallelic site
        args.nonbiallelic++;
        return args.discard ? NULL : ret;
    }

    int ia = nt2int(rec->d.allele[0][0]);
    if ( ia<0 )
    {
        // not A,C,G,T
        args.nonACGT++;
        return args.discard ? NULL : ret;
    }

    int ib = nt2int(rec->d.allele[1][0]);
    if ( ib<0 )
    {
        // not A,C,G,T
        args.nonACGT++;
        return args.discard ? NULL : ret;
    }

    if ( ia==ib )
    {
        // should not happen in well-formed VCF
        args.nonSNP++;
        return args.discard ? NULL : ret;
    }
    args.count[ia][ib]++;

    if ( ir==ia ) args.nok++;

    if ( args.mode==MODE_FLIP2FWD )
    {
        int pair = 1 << ia | 1 << ib;
        if ( pair==0x9 || pair==0x6 )   // skip ambiguous pairs: A/T or C/G
        {
            args.nerr++;
            return args.discard ? NULL : ret;
        }
        if ( ir==ia ) return ret;
        if ( ir==ib ) { args.nswap++; return set_ref_alt(&args,rec,int2nt(ib),int2nt(ia),1); }
        if ( ir==revint(ia) ) { args.nflip++; return set_ref_alt(&args,rec,int2nt(revint(ia)),int2nt(revint(ib)),0); }
        if ( ir==revint(ib) ) { args.nflip_swap++; return set_ref_alt(&args,rec,int2nt(revint(ib)),int2nt(revint(ia)),1); }
        error("FIXME: this should not happen %s:%d\n", bcf_seqname(args.hdr,rec),rec->pos+1);
    }
    if ( args.mode==MODE_TOP2FWD )
    {
        if ( ia==0 && (ib==1 || ib==2) )    // unambiguous pair: A/C or A/G
        {
            if ( ir==ia ) return ret;

            int ia_rev = revint(ia);
            if ( ir==ia_rev )               // vcfref is A, faref is T, flip
            {
                args.nflip++;
                return set_ref_alt(&args,rec,int2nt(ia_rev),int2nt(revint(ib)),0);
            }
            if ( ir==ib )                   // vcfalt is faref, swap
            {
                args.nswap++;
                return set_ref_alt(&args,rec,int2nt(ib),int2nt(ia),1);
            }
            assert( ib==revint(ir) );

            args.nflip_swap++;
            return set_ref_alt(&args,rec,int2nt(revint(ib)),int2nt(revint(ia)),1);
        }
        else    // ambiguous pair, sequence walking must be performed
        {
            int win = rec->pos > 100 ? 100 : rec->pos, beg = rec->pos - win, end = rec->pos + win;
            ref = faidx_fetch_seq(args.fai, (char*)bcf_seqname(args.hdr,rec), beg,end, &len);
            if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", bcf_seqname(args.hdr,rec),rec->pos+1);
            if ( end - beg + 1 != len ) error("FIXME: check win=%d,len=%d at %s:%d  (%d %d %d)\n", win,len, bcf_seqname(args.hdr,rec),rec->pos+1);

            int i, mid = rec->pos - beg, strand = 0;
            for (i=1; i<=win; i++)
            {
                int ra = nt2int(ref[mid-i]);
                int rb = nt2int(ref[mid+i]);
                if ( ra<0 || rb<0 || ra==rb ) continue;     // skip N's and non-infomative pairs: A/A, C/C, G/G, T/T
                int pair = 1 << ra | 1 << rb;
                if ( pair==0x9 || pair==0x6 ) continue;     // skip ambiguous pairs: A/T or C/G
                strand = ra & 0x9 ? 1 : -1;
                break;
            }
            free(ref);
            
            if ( strand==1 )
            {
                if ( ir==ia ) return ret;
                if ( ir==ib )
                {
                    args.nswap++;
                    return set_ref_alt(&args,rec,int2nt(ib),int2nt(ia),1);
                }
            }
            else if ( strand==-1 )
            {
                int ia_rev = revint(ia);
                int ib_rev = revint(ib);
                if ( ir==ia_rev )
                {
                    args.nflip++;
                    return set_ref_alt(&args,rec,int2nt(ia_rev),int2nt(ib_rev),0);
                }
                if ( ir==ib_rev )
                {
                    args.nflip_swap++;
                    return set_ref_alt(&args,rec,int2nt(ib_rev),int2nt(ia_rev),1);
                }
            }

            args.nerr++;
            return args.discard ? NULL : ret;
        }
    }
    return ret;
}

int top_mask[4][4] = 
{ 
    {0,1,1,1},
    {0,0,1,0},
    {0,0,0,0},
    {0,0,0,0},
};
int bot_mask[4][4] = 
{ 
    {0,0,0,0},
    {0,0,0,0},
    {0,1,0,0},
    {1,1,1,0},
};

void destroy(void)
{
    uint32_t i,j,tot = 0;
    uint32_t top_err = 0, bot_err = 0;
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            tot += args.count[i][j];
            if ( !top_mask[i][j] && args.count[i][j] ) top_err++;
            if ( !bot_mask[i][j] && args.count[i][j] ) bot_err++;
        }
    }
    uint32_t nskip = args.nonACGT+args.nonSNP+args.nonbiallelic;
    uint32_t ncmp  = args.nsite - nskip;

    fprintf(stderr,"# SC, guessed strand convention\n");
    fprintf(stderr,"SC\tTOP-compatible\t%d\n",top_err?0:1);
    fprintf(stderr,"SC\tBOT-compatible\t%d\n",bot_err?0:1);

    fprintf(stderr,"# ST, substitution types\n");
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            if ( i==j ) continue;
            fprintf(stderr,"ST\t%c>%c\t%u\t%.1f%%\n", int2nt(i),int2nt(j),args.count[i][j], args.count[i][j]*100./tot);
        }
    }
    fprintf(stderr,"# NS, Number of sites:\n");
    fprintf(stderr,"NS\ttotal        \t%u\n", args.nsite);
    fprintf(stderr,"NS\tref match    \t%u\t%.1f%%\n", args.nok,100.*args.nok/ncmp);
    fprintf(stderr,"NS\tref mismatch \t%u\t%.1f%%\n", ncmp-args.nok,100.*(ncmp-args.nok)/ncmp);
    if ( args.mode!=MODE_STATS )
    {
        fprintf(stderr,"NS\tflipped      \t%u\t%.1f%%\n", args.nflip,100.*args.nflip/(args.nsite-nskip));
        fprintf(stderr,"NS\tswapped      \t%u\t%.1f%%\n", args.nswap,100.*args.nswap/(args.nsite-nskip));
        fprintf(stderr,"NS\tflip+swap    \t%u\t%.1f%%\n", args.nflip_swap,100.*args.nflip_swap/(args.nsite-nskip));
        fprintf(stderr,"NS\tunresolved   \t%u\t%.1f%%\n", args.nerr,100.*args.nerr/(args.nsite-nskip));
    }
    fprintf(stderr,"NS\tskipped      \t%u\n", nskip);
    fprintf(stderr,"NS\tnon-ACGT     \t%u\n", args.nonACGT);
    fprintf(stderr,"NS\tnon-SNP      \t%u\n", args.nonSNP);
    fprintf(stderr,"NS\tnon-biallelic\t%u\n", args.nonbiallelic);
    
    //fprintf(stderr,"Number of correct/swapped/flipped/unresolved sites:\t%u\t%u\t%u\t%u\n", args.nok,args.nswap,args.nflip,args.nerr);
    //fprintf(stderr,"Percentage of correct/swapped/flipped/unresolved sites:\t%.1f\t%.1f\t%.1f\t%.1f\n", 
    //    100.*args.nok/args.ncmp,100.*args.nswap/args.ncmp,100.*args.nflip/args.ncmp,100.*args.nerr/args.ncmp);

    free(args.gts);
    if ( args.fai ) fai_destroy(args.fai);
}
