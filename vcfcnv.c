/* The MIT License

   Copyright (c) 2014 Genome Research Ltd.

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
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kstring.h>
#include <htslib/kfunc.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"
#include "HMM.h"

#define DBG_HMM_PRN 0

#define N_STATES 5
#define CN0 0
#define CN1 1
#define CN2 2
#define CN3 3
#define CNx 4

typedef struct
{
    char *name;
    int idx;
    double pobs[N_STATES];
    FILE *dat_fh, *cn_fh, *summary_fh;
    char *dat_fname, *cn_fname, *summary_fname;
}
sample_t;

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    int prev_rid, ntot, nused;
    sample_t query_sample, control_sample;

    int nstates;    // number of states: N_STATES for one sample, N_STATES^2 for two samples
    double baf_sigma2, lrr_sigma2;          // squared std dev of B-allele frequency and LRR distribution
    double lrr_bias, baf_bias;              // LRR/BAF weights
    double same_prob, ij_prob;              // prior of both samples being the same and the transition probability P(i|j)
    double err_prob;                        // constant probability of erroneous measurement
    double pRR, pRA, pAA, pRR_dflt, pRA_dflt, pAA_dflt;

    double *tprob, *tprob_arr;  // array of transition matrices, precalculated up to ntprob_arr positions
    int ntprob_arr;

    hmm_t *hmm;
    double *eprob;          // emission probs [nstates*nsites,msites]
    uint32_t *sites;        // positions [nsites,msites]
    int nsites, msites;

    double baum_welch_th;
    float plot_th;
    FILE *summary_fh;
    char **argv, *regions_list, *summary_fname, *output_dir;
    char *targets_list, *af_fname;
    int argc;
}
args_t;

FILE *open_file(char **fname, const char *mode, const char *fmt, ...);


static inline void hmm2cn_state(int nstates, int i, int *a, int *b)
{
    *a = i / N_STATES;
    *b = i - (*a)*N_STATES;
}
static double *init_tprob_matrix(int ndim, double ij_prob, double same_prob)
{
    int i,j;
    double *mat = (double*) malloc(sizeof(double)*ndim*ndim);

    assert( ndim==N_STATES || ndim==N_STATES*N_STATES);
    double pii = 1 - ij_prob*(N_STATES-1);

    if ( ndim==N_STATES )   // one sample
    {
        if ( pii < ij_prob ) error("Error: -x set a bit too high, P(x|x) < P(x|y): %e vs %e\n", pii,ij_prob);
        for (j=0; j<ndim; j++)
        {
            double sum = 0;
            for (i=0; i<ndim; i++)
            {
                // transition from j-th to i-th state
                if ( i==j )
                    MAT(mat,ndim,i,j) = pii;
                else
                    MAT(mat,ndim,i,j) = ij_prob;
                sum += MAT(mat,ndim,i,j);
            }
            assert( fabs(sum - 1.0)<1e-15 );
        }
    }
    else    // two samples
    {
        for (j=0; j<ndim; j++)
        {
            int ja,jb;
            hmm2cn_state(ndim, j, &ja, &jb);

            double sum = 0;
            for (i=0; i<ndim; i++)
            {
                int ia,ib;
                hmm2cn_state(ndim, i, &ia, &ib);

                // transition from (ja,jb)-th to (ia,ib)-th state
                double pa = ja==ia ? pii : ij_prob;
                double pb = jb==ib ? pii : ij_prob;

                if ( ia==ib && ja==jb )
                    MAT(mat,ndim,i,j) = pa*pb - pa*pb*same_prob + sqrt(pa*pb)*same_prob;
                else if ( ia==ib )
                    MAT(mat,ndim,i,j) = pa*pb;
                else
                    MAT(mat,ndim,i,j) = pa*pb*(1-same_prob);

                sum += MAT(mat,ndim,i,j);
            }
            for (i=0; i<ndim; i++) MAT(mat,ndim,i,j) /= sum;
        }
    }
    return mat;
}

static void init_sample_files(sample_t *smpl, char *dir)
{
    smpl->dat_fh = open_file(&smpl->dat_fname,"w","%s/dat.%s.tab",dir,smpl->name);
    smpl->cn_fh  = open_file(&smpl->cn_fname,"w","%s/cn.%s.tab",dir,smpl->name);
    smpl->summary_fh = open_file(&smpl->summary_fname,"w","%s/summary.%s.tab",dir,smpl->name);
    fprintf(smpl->dat_fh,"# [1]Chromosome\t[2]Position\t[3]BAF\t[4]LRR\n");
    fprintf(smpl->cn_fh,"# [1]Chromosome\t[2]Position\t[3]CN\t[4]P(CN0)\t[5]P(CN1)\t[6]P(CN2)\t[7]P(CN3)\t[8]P(CNx)\n");
    fprintf(smpl->summary_fh,"# RG, Regions [2]Chromosome\t[3]Start\t[4]End\t[5]Copy Number state\t[6]Quality\n");
}
static void close_sample_files(sample_t *smpl)
{
    fclose(smpl->dat_fh);
    fclose(smpl->cn_fh);
    fclose(smpl->summary_fh);
}

static void init_data(args_t *args)
{
    args->prev_rid = -1;
    args->hdr = args->files->readers[0].header;

    if ( !args->query_sample.name )
    {
        if ( bcf_hdr_nsamples(args->hdr)>1 ) error("Multi-sample VCF, missing the -s option\n");
        args->query_sample.name = strdup(args->hdr->samples[0]);
    }
    else 
        if ( bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,args->query_sample.name)<0 ) error("The sample \"%s\" not found\n", args->query_sample.name);
    if ( !args->files->readers[0].file->is_bin )
    {
        int ret;
        kstring_t tmp = {0,0,0};
        if ( args->control_sample.name )
        {
            ksprintf(&tmp, "%s,%s", args->query_sample.name,args->control_sample.name);
            ret = bcf_hdr_set_samples(args->hdr, tmp.s, 0);
        }
        else
        {
            ret = bcf_hdr_set_samples(args->hdr, args->query_sample.name, 0);
            tmp.s = args->query_sample.name;
        }
        if ( ret<0 ) error("Error parsing the list of samples: %s\n", tmp.s);
        else if ( ret>0 ) error("The sample not found in the VCF: %s\n", ret==1 ? args->query_sample.name : args->control_sample.name);

        if ( args->control_sample.name ) free(tmp.s);
    }
    args->query_sample.idx = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,args->query_sample.name);
    args->control_sample.idx = args->control_sample.name ? bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,args->control_sample.name) : -1;
    args->nstates = args->control_sample.name ? N_STATES*N_STATES : N_STATES;
    args->tprob = init_tprob_matrix(args->nstates, args->ij_prob, args->same_prob);
    args->hmm = hmm_init(args->nstates, args->tprob, 10000);

    args->summary_fh = stdout;
    if ( args->output_dir )
    {
        init_sample_files(&args->query_sample, args->output_dir);
        if ( args->control_sample.name )
        {
            init_sample_files(&args->control_sample, args->output_dir);
            args->summary_fh = open_file(&args->summary_fname,"w","%s/summary.tab",args->output_dir);
        }
        else
            args->summary_fh = NULL;    // one sample only, no two-file summary
    }

    int i;
    FILE *fh = args->summary_fh ? args->summary_fh : args->query_sample.summary_fh;

    fprintf(fh, "# This file was produced by: bcftools cnv(%s+htslib-%s)\n", bcftools_version(),hts_version());
    fprintf(fh, "# The command line was:\tbcftools %s", args->argv[0]);
    for (i=1; i<args->argc; i++) fprintf(fh, " %s",args->argv[i]);
    if ( args->control_sample.name )
        fprintf(fh, "\n#\n"
                "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:%s\t[6]Copy number:%s\t[7]Quality\n",
                args->query_sample.name,args->control_sample.name
               );
    else
        fprintf(fh, "\n#\n"
                "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:%s\t[6]Quality\n",
                args->query_sample.name
               );
}

char *msprintf(const char *fmt, ...);
static void py_plot_cnv(char *script, float th)
{
    if ( th>100 ) return;   // create no plots

    char *cmd = msprintf("python %s -p %f", script, th);
    int ret = system(cmd);
    if ( ret) fprintf(stderr, "The command returned non-zero status %d: %s\n", ret, cmd);
    free(cmd);
}

static void plot_sample(args_t *args, sample_t *smpl)
{
    char *fname;
    FILE *fp = open_file(&fname,"w","%s/plot.%s.py",args->output_dir,smpl->name);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import csv\n"
            "import numpy as np\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "\n"
            "dat = {}\n"
            "with open('%s', 'rb') as f:\n"
            "    reader = csv.reader(f, 'tab')\n"
            "    for row in reader:\n"
            "        chr = row[0]\n"
            "        if chr[0]=='#': continue\n"
            "        if chr not in dat: dat[chr] = []\n"
            "        dat[chr].append([row[1], float(row[2]), float(row[3])])\n"
            "\n"
            "cnv = {}\n"
            "with open('%s', 'rb') as f:\n"
            "    reader = csv.reader(f, 'tab')\n"
            "    for row in reader:\n"
            "        chr = row[0]\n"
            "        if chr[0]=='#': continue\n"
            "        if chr not in cnv: cnv[chr] = []\n"
            "        row[2] = int(row[2]) + 0.5\n"
            "        cnv[chr].append(row[1:])\n"
            "\n"
            "for chr in dat:\n"
            "    fig,(ax1, ax2, ax3) = plt.subplots(3,1,figsize=(10,8),sharex=True)\n"
            "    ax1.plot([x[0] for x in dat[chr]],[x[2] for x in dat[chr]], '.', ms=3)\n"
            "    ax2.plot([x[0] for x in dat[chr]],[x[1] for x in dat[chr]], '.', ms=3)\n"
            "    cn_dat = cnv[chr]\n"
            "    xgrid = [float(x[0]) for x in cn_dat]\n"
            "    ygrid = np.linspace(0,5,6)\n"
            "    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n"
            "    heat = np.zeros_like(xgrid)\n"
            "    for x in range(len(heat[0])-1):\n"
            "       heat[0][x] = cn_dat[x][2]\n"
            "       heat[1][x] = cn_dat[x][3]\n"
            "       heat[2][x] = cn_dat[x][4]\n"
            "       heat[3][x] = cn_dat[x][5]\n"
            "       heat[4][x] = cn_dat[x][6]\n"
            "    mesh = ax3.pcolormesh(xgrid, ygrid, heat, cmap='bwr_r')\n"
            "    mesh.set_clim(vmin=-1,vmax=1)\n"
            "    ax3.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'.-',ms=3,color='black')\n"
            "    fig.suptitle('%s (chr '+chr+')')\n"
            "    ax1.set_title('Log-R intensities Ratio',fontsize=10)\n"
            "    ax2.set_title('B-Allele Frequency',fontsize=10)\n"
            "    ax3.set_title('Copy Number Variation',fontsize=10)\n"
            "    ax1.set_ylabel('LRR')\n"
            "    ax2.set_ylabel('BAF')\n"
            "    ax3.set_ylabel('CN')\n"
            "    ax3.set_xlabel('Coordinate (chrom '+chr+')',fontsize=10)\n"
            "    ax3.set_ylim(-0.1,5.1)\n"
            "    ax3.set_yticks([0.5,1.5,2.5,3.5,4.5])\n"
            "    ax3.set_yticklabels(['CN0','CN1','CN2','CN3','CN4'])\n"
            "    plt.subplots_adjust(left=0.08,right=0.95,bottom=0.08,top=0.92)\n"
            "    plt.savefig('%s/plot.%s.chr'+chr+'.png')\n"
            "    plt.close()\n"
            "\n", 
            smpl->dat_fname,smpl->cn_fname,smpl->name,args->output_dir,smpl->name
    );
    fclose(fp);

    py_plot_cnv(fname, args->plot_th);
    free(fname);
}

static void create_plots(args_t *args)
{
    close_sample_files(&args->query_sample);
    if ( args->control_sample.name ) close_sample_files(&args->control_sample);
    if ( args->summary_fh ) fclose(args->summary_fh);

    if ( !args->control_sample.name )
    {
        plot_sample(args, &args->query_sample);
        return;
    }

    char *fname;
    FILE *fp = open_file(&fname,"w","%s/plot.%s.%s.py",args->output_dir,args->control_sample.name,args->query_sample.name);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import csv,argparse\n"
            "import numpy as np\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "\n"
            "control_sample = '%s'\n"
            "query_sample   = '%s'\n"
            "\n"
            "parser = argparse.ArgumentParser()\n"
            "parser.add_argument('-p', '--plot-threshold', type=float)\n"
            "parser.add_argument('-c', '--chromosome')\n"
            "args = parser.parse_args()\n"
            "if args.plot_threshold==None: args.plot_threshold = 0\n"
            "\n"
            "def chroms_to_plot(th):\n"
            "   dat = {}\n"
            "   with open('%s/summary.tab', 'rb') as f:\n"
            "       reader = csv.reader(f, 'tab')\n"
            "       for row in reader:\n"
            "           if row[0]!='RG': continue\n"
            "           chr   = row[1]\n"
            "           start = row[2]\n"
            "           end   = row[3]\n"
            "           qual  = float(row[6])\n"
            "           if row[4]==row[5] and args.plot_threshold!=0: continue\n"
            "           if chr not in dat: dat[chr] = 0.0\n"
            "           if qual > dat[chr]: dat[chr] = qual\n"
            "   out = {}\n"
            "   for chr in dat:\n"
            "       if (chr not in dat) or dat[chr]<th: continue\n"
            "       out[chr] = 1\n"
            "   return out\n"
            "if args.chromosome!=None:\n"
            "   plot_chroms = { args.chromosome:1 }\n"
            "else:\n"
            "   plot_chroms = chroms_to_plot(args.plot_threshold)\n"
            "\n"
            "def read_dat(file,dat):\n"
            "   with open(file, 'rb') as f:\n"
            "       reader = csv.reader(f, 'tab')\n"
            "       for row in reader:\n"
            "           chr = row[0]\n"
            "           if chr[0]=='#': continue\n"
            "           if chr not in plot_chroms: continue\n"
            "           if chr not in dat: dat[chr] = []\n"
            "           dat[chr].append([row[1], float(row[2]), float(row[3])])\n"
            "def read_cnv(file,cnv):\n"
            "   with open(file, 'rb') as f:\n"
            "       reader = csv.reader(f, 'tab')\n"
            "       for row in reader:\n"
            "           chr = row[0]\n"
            "           if chr[0]=='#': continue\n"
            "           if chr not in cnv: cnv[chr] = []\n"
            "           row[2] = int(row[2]) + 0.5\n"
            "           cnv[chr].append(row[1:])\n"
            "def find_diffs(a,b):\n"
            "    out = []\n"
            "    diff = []\n"
            "    for i in range(len(a)):\n"
            "        if a[i][1]!=b[i][1]:\n"
            "            if i>0: diff.append([b[i-1][0],b[i-1][1],a[i-1][1]])\n"
            "            diff.append([b[i][0],b[i][1],a[i][1]])\n"
            "        elif len(diff):\n"
            "            diff.append([b[i][0],b[i][1],a[i][1]])\n"
            "            out.append(diff)\n"
            "            diff = []\n"
            "    if len(diff): out.append(diff)\n"
            "    return out\n"
            "\n"
            "control_dat = {}\n"
            "control_cnv = {}\n"
            "query_dat   = {}\n"
            "query_cnv   = {}\n"
            "read_dat('%s',control_dat)\n"
            "read_dat('%s',query_dat)\n"
            "read_cnv('%s',control_cnv)\n"
            "read_cnv('%s',query_cnv)\n"
            "\n"
            "for chr in query_dat:\n"
            "    fig,(ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6,1,figsize=(10,8),sharex=True)\n"
            "    ax1.plot([x[0] for x in control_dat[chr]],[x[2] for x in control_dat[chr]], '.', ms=3,color='red')\n"
            "    ax2.plot([x[0] for x in control_dat[chr]],[x[1] for x in control_dat[chr]], '.', ms=3,color='red')\n"
            "    cn_dat = control_cnv[chr]\n"
            "    xgrid = [float(x[0]) for x in cn_dat]\n"
            "    ygrid = np.linspace(0,5,6)\n"
            "    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n"
            "    heat = np.zeros_like(xgrid)\n"
            "    for x in range(len(heat[0])-1):\n"
            "       heat[0][x] = cn_dat[x][2]\n"
            "       heat[1][x] = cn_dat[x][3]\n"
            "       heat[2][x] = cn_dat[x][4]\n"
            "       heat[3][x] = cn_dat[x][5]\n"
            "       heat[4][x] = cn_dat[x][6]\n"
            "    mesh = ax3.pcolormesh(xgrid, ygrid, heat, cmap='bwr')\n"
            "    mesh.set_clim(vmin=-1,vmax=1)\n"
            "    ax3.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'-',ms=3,color='black',lw=1.7)\n"
            "\n"
            "    ax6.plot([x[0] for x in query_dat[chr]],[x[2] for x in query_dat[chr]], '.', ms=3)\n"
            "    ax5.plot([x[0] for x in query_dat[chr]],[x[1] for x in query_dat[chr]], '.', ms=3)\n"
            "    cn_dat = query_cnv[chr]\n"
            "    xgrid = [float(x[0]) for x in cn_dat]\n"
            "    ygrid = np.linspace(0,5,6)\n"
            "    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n"
            "    heat = np.zeros_like(xgrid)\n"
            "    for x in range(len(heat[0])-1):\n"
            "       heat[0][x] = cn_dat[x][2]\n"
            "       heat[1][x] = cn_dat[x][3]\n"
            "       heat[2][x] = cn_dat[x][4]\n"
            "       heat[3][x] = cn_dat[x][5]\n"
            "       heat[4][x] = cn_dat[x][6]\n"
            "    mesh = ax4.pcolormesh(xgrid, ygrid, heat, cmap='bwr_r')\n"
            "    mesh.set_clim(vmin=-1,vmax=1)\n"
            "    ax4.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'-',ms=3,color='black',lw=1.7)\n"
            "    ax3.annotate(control_sample, xy=(0.02,0.1), xycoords='axes fraction', color='red',fontsize=12, va='bottom',ha='left')\n"
            "    ax4.annotate(query_sample, xy=(0.02,0.9), xycoords='axes fraction', color='blue',fontsize=12, va='top',ha='left')\n"
            "\n"
            "    diffs = find_diffs(control_cnv[chr],query_cnv[chr])\n"
            "    for diff in diffs:\n"
            "        ax3.plot([x[0] for x in diff],[x[1] for x in diff],'-',ms=3,color='blue',lw=1.7)\n"
            "        ax4.plot([x[0] for x in diff],[x[2] for x in diff],'-',ms=3,color='red',lw=1.7)\n"
            "\n"
            "    fig.suptitle('chr '+chr+', '+control_sample+' vs '+query_sample)\n"
            "    ax1.tick_params(axis='both', labelsize=8)\n"
            "    ax2.tick_params(axis='both', labelsize=8)\n"
            "    ax3.tick_params(axis='both', labelsize=8)\n"
            "    ax4.tick_params(axis='both', labelsize=8)\n"
            "    ax5.tick_params(axis='both', labelsize=8)\n"
            "    ax6.tick_params(axis='both', labelsize=8)\n"
            "    ax6.set_xlabel('Coordinate (chrom '+chr+')',fontsize=8)\n"
            "    ax1.set_ylabel('LRR')\n"
            "    ax2.set_ylabel('BAF')\n"
            "    ax3.set_ylabel('CN')\n"
            "    ax6.set_ylabel('LRR')\n"
            "    ax5.set_ylabel('BAF')\n"
            "    ax4.set_ylabel('CN')\n"
            "    ax3.set_ylim(-0.1,5.1)\n"
            "    ax3.set_yticks([0.5,1.5,2.5,3.5,4.5])\n"
            "    ax3.set_yticklabels(['CN0','CN1','CN2','CN3','CN4'])\n"
            "    ax4.set_ylim(-0.1,5.1)\n"
            "    ax4.set_yticks([0.5,1.5,2.5,3.5,4.5])\n"
            "    ax4.set_yticklabels(['CN0','CN1','CN2','CN3','CN4'])\n"
            "    plt.subplots_adjust(left=0.08,right=0.95,bottom=0.08,top=0.92,hspace=0)\n"
            "    plt.savefig('%s/plot.%s.%s.chr'+chr+'.png')\n"
            "    plt.close()\n"
            "\n", 
            args->control_sample.name,args->query_sample.name,
            args->output_dir,
            args->control_sample.dat_fname,args->query_sample.dat_fname,
            args->control_sample.cn_fname,args->query_sample.cn_fname,
            args->output_dir,args->control_sample.name,args->query_sample.name
        );
        fclose(fp);

    py_plot_cnv(fname,args->plot_th);
    free(fname);
}

static void destroy_data(args_t *args)
{
    bcf_sr_destroy(args->files);
    hmm_destroy(args->hmm);
    free(args->sites);
    free(args->eprob);
    free(args->tprob);
    free(args->summary_fname);
    free(args->query_sample.name);
    free(args->query_sample.dat_fname);
    free(args->query_sample.cn_fname);
    free(args->query_sample.summary_fname);
    free(args->control_sample.dat_fname);
    free(args->control_sample.cn_fname);
    free(args->control_sample.summary_fname);
}

static inline char copy_number_state(args_t *args, int istate, int ismpl)
{
    char code[] = "01234";
    if ( !args->control_sample.name ) return code[istate];
    int idx = ismpl ? istate - (istate/N_STATES)*N_STATES : istate/N_STATES;
    return code[idx];
}

static double phred_score(double prob)
{
    if ( prob==0 ) return 99;
    prob = -4.3429*log(prob);
    return prob>99 ? 99 : prob;
}

static double avg_ii_prob(int n, double *mat)
{
    int i;
    double avg = 0;
    for (i=0; i<n; i++) avg += MAT(mat,n,i,i);
    return avg/n;
}

static void cnv_flush_viterbi(args_t *args)
{
    if ( !args->nsites ) return;

    hmm_t *hmm = args->hmm;
    hmm_set_tprob(args->hmm, args->tprob, 10000);
    while ( args->baum_welch_th!=0 )
    {
        double ori_ii = avg_ii_prob(hmm->nstates,hmm->tprob_arr);
        hmm_run_baum_welch(hmm, args->nsites, args->eprob, args->sites);
        double new_ii = avg_ii_prob(hmm->nstates,hmm->tprob_arr);
        fprintf(stderr,"%e\t%e\t%e\n", ori_ii,new_ii,new_ii-ori_ii);
        double *tprob = init_tprob_matrix(hmm->nstates, 1-new_ii, args->same_prob);
        hmm_set_tprob(args->hmm, tprob, 10000);
        free(tprob);
        if ( fabs(new_ii - ori_ii) < args->baum_welch_th )
        {
            int i,j;
            for (i=0; i<hmm->nstates; i++)
            {
                for (j=0; j<hmm->nstates; j++)
                {
                    printf(" %.15f", MAT(hmm->tprob_arr,hmm->nstates,j,i));
                }
                printf("\n");
            }
            break;
        }
    }
    hmm_run_viterbi(hmm, args->nsites, args->eprob, args->sites);
    hmm_run_fwd_bwd(hmm, args->nsites, args->eprob, args->sites);


    // Output the results
    double qual = 0;
    int i,j, isite, start_cn = hmm->vpath[0], start_pos = args->sites[0], istart_pos = 0;
    for (isite=0; isite<args->nsites; isite++)
    {
        int state = hmm->vpath[args->nstates*isite];
        double *pval = hmm->fwd + isite*args->nstates;

        qual += pval[start_cn];

        // output CN and fwd-bwd likelihood for each site
        if ( args->query_sample.cn_fh )
        {
            fprintf(args->query_sample.cn_fh, "%s\t%d\t%c", bcf_hdr_id2name(args->hdr,args->prev_rid), args->sites[isite]+1, copy_number_state(args,state,0));
            if ( !args->control_sample.cn_fh )
                for (i=0; i<args->nstates; i++) fprintf(args->query_sample.cn_fh, "\t%f", pval[i]);
            else
                for (i=0; i<N_STATES; i++)
                {
                    double sum = 0;
                    for (j=0; j<N_STATES; j++) sum += pval[i*N_STATES+j];
                    fprintf(args->query_sample.cn_fh, "\t%f", sum);
                }
            fprintf(args->query_sample.cn_fh, "\n");
        }
        if ( args->control_sample.cn_fh )
        {
            fprintf(args->control_sample.cn_fh, "%s\t%d\t%c", bcf_hdr_id2name(args->hdr,args->prev_rid), args->sites[isite]+1, copy_number_state(args,state,1));
            for (i=0; i<N_STATES; i++)
            {
                double sum = 0;
                for (j=0; j<N_STATES; j++) sum += pval[i+N_STATES*j];
                fprintf(args->control_sample.cn_fh, "\t%f", sum);
            }
            fprintf(args->control_sample.cn_fh, "\n");
        }

        if ( start_cn != state )
        {
            char start_cn_query = copy_number_state(args,start_cn,0);
            qual = phred_score(1 - qual/(isite - istart_pos));
            fprintf(args->query_sample.summary_fh,"RG\t%s\t%d\t%d\t%c\t%.1f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite],start_cn_query,qual);

            if ( args->control_sample.name )
            {
                // regions 0-based, half-open
                char start_cn_ctrl = copy_number_state(args,start_cn,1);
                fprintf(args->control_sample.summary_fh,"RG\t%s\t%d\t%d\t%c\t%.1f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite],start_cn_ctrl,qual);
                fprintf(args->summary_fh,"RG\t%s\t%d\t%d\t%c\t%c\t%.1f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite],start_cn_query,start_cn_ctrl,qual);
            }

            istart_pos = isite;
            start_pos = args->sites[isite];
            start_cn = state;
            qual = 0;
        }
    }
    qual = phred_score(1 - qual/(isite - istart_pos));
    char start_cn_query = copy_number_state(args,start_cn,0);
    fprintf(args->query_sample.summary_fh,"RG\t%s\t%d\t%d\t%c\t%.1f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite-1]+1,start_cn_query,qual);
    if ( args->control_sample.name )
    {
        char start_cn_ctrl = copy_number_state(args,start_cn,1);
        fprintf(args->control_sample.summary_fh,"RG\t%s\t%d\t%d\t%c\t%.1f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite-1]+1,start_cn_ctrl,qual);
        fprintf(args->summary_fh,"RG\t%s\t%d\t%d\t%c\t%c\t%.1f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite-1]+1,start_cn_query,start_cn_ctrl,qual);
    }
}

static int set_observed_prob(args_t *args, bcf_fmt_t *baf_fmt, bcf_fmt_t *lrr_fmt, sample_t *smpl)
{
    float baf, lrr;
    baf = ((float*)(baf_fmt->p + baf_fmt->size*smpl->idx))[0];
    if ( bcf_float_is_missing(baf) || isnan(baf) ) baf = -0.1;    // arbitrary negative value == missing value

    lrr = ((float*)(lrr_fmt->p + lrr_fmt->size*smpl->idx))[0];
    if ( bcf_float_is_missing(lrr) || isnan(lrr) ) baf = -0.1;

    if ( baf>=0 )    // skip missing values
        fprintf(smpl->dat_fh,"%s\t%d\t%.3f\t%.3f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), args->sites[args->nsites-1]+1,baf,lrr);

    if ( baf<0 )
    {
        // no call: either some technical issue or the call could not be made because it is CN0
        int i;
        smpl->pobs[CN0] = 0.5;
        for (i=1; i<N_STATES; i++) smpl->pobs[i] = (1.0-smpl->pobs[CN0])/(N_STATES-1);
        return 0;
    }

    double pk0, pk14, pk13, pk12, pk23, pk34, pk1;
    pk0  = exp(-baf*baf/args->baf_sigma2);
    pk14 = exp(-(baf-1/4.)*(baf-1/4.)/args->baf_sigma2);
    pk13 = exp(-(baf-1/3.)*(baf-1/3.)/args->baf_sigma2);
    pk12 = exp(-(baf-1/2.)*(baf-1/2.)/args->baf_sigma2);
    pk23 = exp(-(baf-2/3.)*(baf-2/3.)/args->baf_sigma2);
    pk34 = exp(-(baf-3/4.)*(baf-3/4.)/args->baf_sigma2);
    pk1  = exp(-(baf-1.0)*(baf-1.0)/args->baf_sigma2);

    double cn1_baf, cn2_baf, cn3_baf, cn4_baf;
    cn1_baf = pk0*(args->pRR+args->pRA/2.)  + pk1*(args->pAA+args->pRA/2.);
    cn2_baf = pk0*args->pRR + pk1*args->pAA + pk12*args->pRA;
    cn3_baf = pk0*args->pRR + pk1*args->pAA + (pk13 + pk23)*args->pRA/2.;
    cn4_baf = pk0*args->pRR + pk1*args->pAA + (pk14 + pk12 + pk34)*args->pRA/3.;

    double cn1_lrr, cn2_lrr, cn3_lrr, cn4_lrr;
    cn1_lrr = exp(-(lrr + 0.45)*(lrr + 0.45)/args->lrr_sigma2);
    cn2_lrr = exp(-(lrr - 0.00)*(lrr - 0.00)/args->lrr_sigma2);
    cn3_lrr = exp(-(lrr - 0.30)*(lrr - 0.30)/args->lrr_sigma2);
    cn4_lrr = exp(-(lrr - 0.75)*(lrr - 0.75)/args->lrr_sigma2);

    smpl->pobs[CN0] = 0;
    smpl->pobs[CN1] = args->err_prob + (1 - args->baf_bias + args->baf_bias*cn1_baf)*(1 - args->lrr_bias + args->lrr_bias*cn1_lrr);
    smpl->pobs[CN2] = args->err_prob + (1 - args->baf_bias + args->baf_bias*cn2_baf)*(1 - args->lrr_bias + args->lrr_bias*cn2_lrr);
    smpl->pobs[CN3] = args->err_prob + (1 - args->baf_bias + args->baf_bias*cn3_baf)*(1 - args->lrr_bias + args->lrr_bias*cn3_lrr);
    smpl->pobs[CNx] = args->err_prob + (1 - args->baf_bias + args->baf_bias*cn4_baf)*(1 - args->lrr_bias + args->lrr_bias*cn4_lrr);

    return 0;
}

static void set_emission_prob(args_t *args)
{
    double *eprob = &args->eprob[args->nstates*(args->nsites-1)];
    int i;
    for (i=0; i<N_STATES; i++)
        eprob[i] = args->query_sample.pobs[i];
}

static void set_emission_prob2(args_t *args)
{
    double *eprob = &args->eprob[args->nstates*(args->nsites-1)];
    int i, j;
    for (i=0; i<N_STATES; i++)
    {
        for (j=0; j<N_STATES; j++)
        {
            eprob[i*N_STATES+j] = args->query_sample.pobs[i]*args->control_sample.pobs[j];
        }
    }
}

int read_AF(bcf_sr_regions_t *tgt, bcf1_t *line, double *alt_freq);

static void cnv_next_line(args_t *args, bcf1_t *line)
{
    if ( !line ) 
    {
        // Done, flush viterbi
        cnv_flush_viterbi(args);
        return;
    }

    if ( line->rid!=args->prev_rid )
    {
        // New chromosome
        cnv_flush_viterbi(args);
        args->prev_rid = line->rid;
        args->nsites = 0;
    }

    // Process line
    args->ntot++;

    bcf_fmt_t *baf_fmt, *lrr_fmt;
    if ( !(baf_fmt = bcf_get_fmt(args->hdr, line, "BAF")) ) return; 
    if ( !(lrr_fmt = bcf_get_fmt(args->hdr, line, "LRR")) ) return;

    // Realloc buffers needed by viterbi and fwd-bwd
    args->nsites++;
    int m = args->msites;
    hts_expand(uint32_t,args->nsites,args->msites,args->sites);
    if ( args->msites!=m )
        args->eprob = (double*) realloc(args->eprob,sizeof(double)*args->msites*args->nstates);
    args->sites[args->nsites-1] = line->pos;

    double alt_freq;
    if ( !args->af_fname || read_AF(args->files->targets, line, &alt_freq) < 0 )
    {
        args->pRR = args->pRR_dflt;
        args->pRA = args->pRA_dflt;
        args->pAA = args->pAA_dflt;
    }
    else
    {
        args->pRR = (1 - alt_freq)*(1 - alt_freq);
        args->pRA = 2*(1 - alt_freq)*alt_freq;
        args->pAA = alt_freq*alt_freq;
    }

    int ret = set_observed_prob(args, baf_fmt,lrr_fmt, &args->query_sample);
    if ( ret<0 ) 
    {
        args->nsites--;
        return;
    }
    if ( args->control_sample.name )
    {
        ret = set_observed_prob(args, baf_fmt,lrr_fmt, &args->control_sample);
        if ( ret<0 )
        {
            args->nsites--;
            return;
        }
    }

    if ( args->control_sample.name )
        set_emission_prob2(args);
    else
        set_emission_prob(args);

    args->nused++;
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Copy number variation caller, requires Illumina's B-allele frequency (BAF) and Log R\n");
    fprintf(stderr, "         Ratio intensity (LRR). The HMM considers the following copy number states: CN 2\n");
    fprintf(stderr, "         (normal), 1 (single-copy loss), 0 (complete loss), 3 (single-copy gain), x (other)\n");
    fprintf(stderr, "Usage:   bcftools cnv [OPTIONS] <file.vcf>\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "    -c, --control-sample <string>      optional control sample name to highlight differences\n");
    fprintf(stderr, "    -f, --AF-file <file>               read allele frequencies from file (CHR\\tPOS\\tREF,ALT\\tAF)\n");
    fprintf(stderr, "    -o, --output-dir <path>            \n");
    fprintf(stderr, "    -p, --plot-threshold <float>       plot aberrant chromosomes with quality at least 'float'\n");
    fprintf(stderr, "    -r, --regions <region>             restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>          restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --query-sample <string>        query samply name\n");
    fprintf(stderr, "    -t, --targets <region>             similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>          similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "HMM Options:\n");
    fprintf(stderr, "    -b, --BAF-weight <float>           relative contribution from BAF [1]\n");
    fprintf(stderr, "    -e, --err-prob <float>             probability of error [1e-4]\n");
    fprintf(stderr, "    -l, --LRR-weight <float>           relative contribution from LRR [0.2]\n");
    fprintf(stderr, "    -P, --same-prob <float>            prior probability of -s/-c being same [1e-1]\n");
    fprintf(stderr, "    -x, --xy-prob <float>              P(x|y) transition probability [1e-8]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfcnv(int argc, char *argv[])
{
    int c;
    args_t *args    = (args_t*) calloc(1,sizeof(args_t));
    args->argc      = argc; args->argv = argv;
    args->files     = bcf_sr_init();
    args->plot_th   = 1e9;   // by default plot none

    // How much FORMAT/LRR and FORMAT/BAF matter
    args->lrr_bias  = 0.2;
    args->baf_bias  = 1.0;
    args->err_prob  = 1e-4;

    // Transition probability to a different state and the prior of both samples being the same
    args->ij_prob   = 1e-8;
    args->same_prob = 1e-1;

    // Squared std dev of BAF and LRR values (gaussian noise), estimated from real data (hets, one sample, one chr)
    args->baf_sigma2 = 0.08*0.08;   // illumina: 0.03
    args->lrr_sigma2 = 0.4*0.4; //0.20*0.20;   // illumina: 0.18

    // Priors for RR, RA, AA genotypes
    args->pRR_dflt = 0.76;
    args->pRA_dflt = 0.14;
    args->pAA_dflt = 0.098;
    // args->pRR = 0.69;
    // args->pRA = 0.18;
    // args->pAA = 0.11;

    int regions_is_file = 0, targets_is_file = 0;
    static struct option loptions[] = 
    {
        {"AF-file",1,0,'f'},
        {"baum-welch",1,0,'W'},
        {"err-prob",1,0,'e'},
        {"BAF-weight",1,0,'b'},
        {"LRR-weight",1,0,'l'},
        {"same-prob",1,0,'P'},
        {"xy-prob",1,0,'x'},
        {"sample",1,0,'s'},
        {"control",1,0,'c'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"plot",1,0,'p'},
        {"output-dir",1,0,'o'},
        {0,0,0,0}
    };
    char *tmp = NULL;
    while ((c = getopt_long(argc, argv, "h?r:R:t:T:s:o:p:l:T:c:b:P:x:e:W:f:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'f': args->af_fname = optarg; break;
            case 'W':
                args->baum_welch_th = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -W %s\n", optarg);
                break;
            case 'e': 
                args->err_prob = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -e %s\n", optarg);
                break;
            case 'b': 
                args->baf_bias = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -b %s\n", optarg);
                break;
            case 'x': 
                args->ij_prob = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -x %s\n", optarg);
                break;
            case 'P': 
                args->same_prob = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -P %s\n", optarg);
                break;
            case 'l': 
                args->lrr_bias = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -l %s\n", optarg);
                break;
            case 'p': 
                args->plot_th = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -p %s\n", optarg);
                break;
            case 'o': args->output_dir = optarg; break;
            case 's': args->query_sample.name = strdup(optarg); break;
            case 'c': args->control_sample.name = optarg; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 'h': 
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";
    }
    else fname = argv[optind];
    if ( !fname ) usage(args);

    if ( args->plot_th<=100 && !args->output_dir ) error("Expected -o option with -p\n");
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( args->af_fname )
    {
        if ( bcf_sr_set_targets(args->files, args->af_fname, 1, 3)<0 )
            error("Failed to read the targets: %s\n", args->af_fname);
    }
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));
    
    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        cnv_next_line(args, line);
    }
    cnv_next_line(args, NULL);
    create_plots(args);
    fprintf(stderr,"Number of lines: total/processed: %d/%d\n", args->ntot,args->nused);
    destroy_data(args);
    free(args);
    return 0;
}


