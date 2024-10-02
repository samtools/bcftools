/* The MIT License

   Copyright (c) 2013-2021 Genome Research Ltd.

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
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/hts_defs.h>
#include "bcftools.h"
#include "peakfit.h"

typedef struct
{
    int nvals;          // all values, including RR,AA peaks
    double *xvals;      // pointer to args_t.xvals
    double *yvals;
    int copy_number;    // heuristics to skip futile CN1 fits when no het peak is detected
    int irr, ira, iaa;  // chop off RR and AA peaks
    char *chr;
}
dist_t;

typedef struct
{
    int ndist, nbins, ra_rr_scaling, smooth;
    double *xvals;
    dist_t *dist;
    char **argv, *output_dir;
    double fit_th, peak_symmetry, cn_penalty, min_peak_size, min_fraction;
    int argc, plot, verbose, regions_is_file, targets_is_file, include_aa, force_cn;
    int regions_overlap, targets_overlap;
    char *dat_fname, *fname, *regions_list, *targets_list, *sample;
    FILE *dat_fp;
}
args_t;

FILE *open_file(char **fname, const char *mode, const char *fmt, ...) HTS_FORMAT(HTS_PRINTF_FMT, 3, 4);

static void init_dist(args_t *args, dist_t *dist, int verbose)
{
    // isolate RR and AA peaks and rescale so that they are comparable to hets
    int i, irr, iaa, n = dist->nvals;

    // smooth the distribution, this is just to find the peaks
    double *tmp = (double*) malloc(sizeof(double)*n);
    int win  = args->smooth ? abs(args->smooth)*2 + 1 : 7;   // must be an odd number
    int hwin = win/2;
    double avg = tmp[0] = dist->yvals[0];
    for (i=1; i<hwin; i++)
    {
        avg += dist->yvals[2*i-1]; 
        tmp[i] = avg/(2*i+1);
    }
    avg = 0;
    for (i=0; i<n; i++)
    {
        avg += dist->yvals[i];
        if ( i>=win-1 )
        {
            tmp[i-hwin] = avg/win;
            avg -= dist->yvals[i-win+1];
        }
    }
    for (i=n-hwin; i<n; i++)
    {
        avg -= dist->yvals[i-hwin];
        hwin--;
        tmp[i] = avg/(2*hwin+1);
        avg -= dist->yvals[i-hwin];
    }

    // find the extremes; first a simple approach: find a gap
    for (irr=0,i=0; i<n/2; i++) if ( tmp[i] < tmp[irr] ) irr = i;
    for (iaa=n-1,i=n-1; i>=n/2; i--) if ( tmp[i] < tmp[iaa] ) iaa = i;
    irr += win*0.5;
    iaa += win*0.5;
    if ( iaa>=n ) iaa = n-1;
    if ( irr>=iaa ) error("FIXME: oops, dist normalization failed for %s: %d vs %d\n", dist->chr,irr,iaa); // we may need to be smarter
    if ( args->smooth>0 ) for (i=0; i<n; i++) dist->yvals[i] = tmp[i];
    free(tmp);

    // clean the data: the AA peak is occasionally not centered at 1.0 but is closer to the center, chop off
    int imax_aa = iaa;
    for (i=iaa; i<n; i++)
        if ( dist->yvals[imax_aa] < dist->yvals[i] ) imax_aa = i;
    dist->nvals = imax_aa+1;
    if ( iaa>=dist->nvals ) iaa = dist->nvals-1;

    // find the maximum and scale the peaks (first draft: no attempt to join the segments smootly)
    double max_rr = 0, max_aa = 0, max_ra = 0, srr = 0, saa = 0, sra = 0;
    for (i=0; i<irr; i++)
    {
        srr += dist->yvals[i];
        if ( max_rr < dist->yvals[i] ) max_rr = dist->yvals[i];
    }
    for (i=irr; i<=iaa; i++)
    {
        sra += dist->yvals[i];
        if ( max_ra < dist->yvals[i] ) max_ra = dist->yvals[i];
    }
    for (i=iaa+1; i<n; i++)
    {
        saa += dist->yvals[i];
        if ( max_aa < dist->yvals[i] ) max_aa = dist->yvals[i];
    }


    // Does the het peak exist at all? Usually the numbers are as follows:
    //  1:  cn=0     ra/rr=0.205730      aa/ra=0.674922      nra=7066
    //  20: cn=0     ra/rr=0.258019      aa/ra=0.548929      nra=2381
    //  X:  cn=0     ra/rr=0.005976      aa/ra=44.116667     nra=60
    //  Y:  cn=0     ra/rr=0.008316      aa/ra=7.250000      nra=12
    //  MT: cn=0     ra/rr=0.013699      aa/ra=0.666667      nra=3

    if ( !args->ra_rr_scaling ) max_ra = max_aa = max_rr;
    if ( !sra || (sra/srr<0.1 && saa/sra>1.0) ) // too few hets, CN1
    {
        max_ra = max_aa;
        dist->copy_number = 1;
    }
    else if ( sra/srr<0.1 || saa/sra>1.0 )
    {
        max_ra = max_aa;
        dist->copy_number = -1;     // unknown copy number
    }
    if ( max_rr ) for (i=0; i<irr; i++) dist->yvals[i] /= max_rr;
    if ( max_ra ) for (i=irr; i<=iaa; i++) dist->yvals[i] /= max_ra;
    if ( max_aa ) for (i=iaa+1; i<n; i++) dist->yvals[i] /= max_aa;

    dist->irr = irr;
    dist->iaa = iaa;
    dist->ira = n*0.5;

    if ( verbose )
        fprintf(stderr,"%s:\t irr,ira,iaa=%.2f,%.2f,%.2f \t cn=%2d \t ra/rr=%f \t aa/ra=%f \t nra=%d\n", 
            dist->chr, dist->xvals[irr],dist->xvals[dist->ira],dist->xvals[iaa],
            dist->copy_number,sra/srr,saa/sra, (int)sra);
}

static void init_data(args_t *args)
{
    bcf_srs_t *files = bcf_sr_init();
    if ( args->regions_list )
    {
        bcf_sr_set_opt(files,BCF_SR_REGIONS_OVERLAP,args->regions_overlap);
        if ( bcf_sr_set_regions(files, args->regions_list, args->regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        bcf_sr_set_opt(files,BCF_SR_TARGETS_OVERLAP,args->targets_overlap);
        if ( bcf_sr_set_targets(files, args->targets_list, args->targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( !bcf_sr_add_reader(files, args->fname) ) error("Failed to read from %s: %s\n", !strcmp("-",args->fname)?"standard input":args->fname,bcf_sr_strerror(files->errnum));
    bcf_hdr_t *hdr = files->readers[0].header;
    if ( !args->sample )
    {
        if ( bcf_hdr_nsamples(hdr)>1 ) error("Missing the option -s, --sample\n");
        args->sample = hdr->samples[0];
    }
    else if ( bcf_hdr_id2int(hdr,BCF_DT_SAMPLE,args->sample)<0 ) error("No such sample: %s\n", args->sample);
    int ret = bcf_hdr_set_samples(hdr, args->sample, 0);
    if ( ret<0 ) error("Error setting the sample: %s\n", args->sample);

    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,bcf_hdr_id2int(hdr,BCF_DT_ID,"BAF")) )
        error("The tag FORMAT/BAF is not present in the VCF: %s\n", args->fname);

    int i;
    args->xvals = (double*) calloc(args->nbins,sizeof(double));
    for (i=0; i<args->nbins; i++) args->xvals[i] = 1.0*i/(args->nbins-1);

    // collect BAF distributions for all chromosomes
    int idist = -1, nbaf = 0, nprocessed = 0, ntotal = 0, prev_chr = -1;
    float *baf = NULL;
    while ( bcf_sr_next_line(files) )
    {
        ntotal++;

        bcf1_t *line = bcf_sr_get_line(files,0);
        if ( bcf_get_format_float(hdr,line,"BAF",&baf,&nbaf) != 1 ) continue;
        if ( bcf_float_is_missing(baf[0]) ) continue;

        nprocessed++;

        if ( prev_chr==-1 || prev_chr!=line->rid )
        {
            // new chromosome
            idist = args->ndist++;
            args->dist = (dist_t*) realloc(args->dist, sizeof(dist_t)*args->ndist);
            memset(&args->dist[idist],0,sizeof(dist_t));
            args->dist[idist].chr   = strdup(bcf_seqname(hdr,line));
            args->dist[idist].yvals = (double*) calloc(args->nbins,sizeof(double));
            args->dist[idist].xvals = args->xvals;
            args->dist[idist].nvals = args->nbins;
            prev_chr = line->rid;
        }
        int bin = baf[0]*(args->nbins-1);
        args->dist[idist].yvals[bin]++;   // the distribution
    }
    free(baf);
    bcf_sr_destroy(files);

    for (idist=0; idist<args->ndist; idist++)
    {
        #if 0
            int j;
            for (j=0; j<args->nbins; j++)
            {
                double x = args->dist[idist].xvals[j];
                args->dist[idist].yvals[j] = exp(-(x-0.5)*(x-0.5)/1e-3);
            }
        #endif
        init_dist(args, &args->dist[idist],args->verbose);
    }

    args->dat_fp = open_file(&args->dat_fname,"w","%s/dist.dat", args->output_dir);
    fprintf(args->dat_fp, "# This file was produced by: bcftools polysomy(%s+htslib-%s), the command line was:\n", bcftools_version(),hts_version());
    fprintf(args->dat_fp, "# \t bcftools %s ", args->argv[0]);
    for (i=1; i<args->argc; i++)
        fprintf(args->dat_fp, " %s",args->argv[i]);
    fprintf(args->dat_fp,"\n#\n");
    fprintf(args->dat_fp,"# DIST\t[2]Chrom\t[3]BAF\t[4]Normalized Count\n");
    fprintf(args->dat_fp,"# FIT\t[2]Goodness of Fit\t[3]iFrom\t[4]iTo\t[5]The Fitted Function\n");
    fprintf(args->dat_fp,"# CN\t[2]Chrom\t[3]Estimated Copy Number\t[4]Absolute fit deviation\n");

    char *fname = NULL;
    FILE *fp = open_file(&fname,"w","%s/dist.py", args->output_dir);
//-------- matplotlib script --------------
    fprintf(fp,
        "#!/usr/bin/env python\n"
        "#\n"
        "import matplotlib as mpl\n"
        "mpl.use('Agg')\n"
        "import matplotlib.pyplot as plt\n"
        "import csv,sys,argparse\n"
        "from math import exp\n"
        "\n"
        "outdir = '%s'\n"
        "\n"
        "def read_dat(dat,fit,cn):\n"
        "   csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)\n"
        "   with open(outdir+'/dist.dat', 'r') as f:\n"
        "      reader = csv.reader(f, 'tab')\n"
        "      for row in reader:\n"
        "          if row[0][0]=='#': continue\n"
        "          type = row[0]\n"
        "          chr  = row[1]\n"
        "          if type=='DIST':\n"
        "              if chr not in dat: dat[chr] = []\n"
        "              dat[chr].append(row)\n"
        "          elif type=='FIT':\n"
        "              if chr not in fit: fit[chr] = []\n"
        "              fit[chr].append(row)\n"
        "          elif type=='CN':\n"
        "              cn[chr] = row[2]\n"
        "\n"
        "def plot_dist(dat,fit,chr):\n"
        "   fig, ax = plt.subplots(1, 1, figsize=(7,5))\n"
        "   ax.plot([x[2] for x in dat[chr]],[x[3] for x in dat[chr]],'k-',label='Distribution')\n"
        "   if chr in fit:\n"
        "       for i in range(len(fit[chr])):\n"
        "           pfit = fit[chr][i]\n"
        "           exec('def xfit(x): return '+pfit[5])\n"
        "           istart = int(pfit[3])\n"
        "           iend   = int(pfit[4])+1\n"
        "           vals   = dat[chr][istart:iend]\n"
        "           args   = {}\n"
        "           if i==0: args = {'label':'Target to Fit'}\n"
        "           ax.plot([x[2] for x in vals],[x[3] for x in vals],'r-',**args)\n"
        "           if i==0: args = {'label':'Best Fit'}\n"
        "           ax.plot([x[2] for x in vals],[xfit(float(x[2])) for x in vals],'g-',**args)\n"
        "   ax.set_title('BAF distribution, chr'+chr)\n"
        "   ax.set_xlabel('BAF')\n"
        "   ax.set_ylabel('Frequency')\n"
        "   ax.legend(loc='best',prop={'size':7},frameon=False)\n"
        "   plt.savefig(outdir+'/dist.chr'+chr+'.png')\n"
        "   plt.close()\n"
        "\n"
        "def plot_copy_number(cn):\n"
        "   fig, ax = plt.subplots(1, 1, figsize=(7,5))\n"
        "   xlabels = sorted(cn.keys())\n"
        "   xvals = range(len(xlabels))\n"
        "   yvals = [float(cn[x]) for x in xlabels]\n"
        "   ax.plot(xvals,yvals,'o',color='red')\n"
        "   for i in range(len(xvals)):\n"
        "       if yvals[i]==-1: ax.annotate('?', xy=(xvals[i],0.5),va='center',ha='center',color='red',fontweight='bold')\n"
        "   ax.tick_params(axis='both', which='major', labelsize=9)\n"
        "   ax.set_xticks(xvals)\n"
        "   ax.set_xticklabels(xlabels,rotation=45)\n"
        "   ax.set_xlim(-1,len(xlabels))\n"
        "   ax.set_ylim(0,5.0)\n"
        "   ax.set_yticks([1.0,2.0,3.0,4.0])\n"
        "   ax.set_xlabel('Chromosome')\n"
        "   ax.set_ylabel('Copy Number')\n"
        "   plt.savefig(outdir+'/copy-number.png')\n"
        "   plt.close()\n"
        "\n"
        "class myParser(argparse.ArgumentParser):\n"
        "   def error(self, message):\n"
        "       self.print_help()\n"
        "       sys.stderr.write('error: %%s\\n' %% message)\n"
        "       sys.exit(2)\n"
        "\n"
        "def main():\n"
        "   parser = myParser()\n"
        "   parser.add_argument('-a', '--all', action='store_true', help='Create all plots')\n"
        "   parser.add_argument('-c', '--copy-number', action='store_true', help='Create copy-number plot')\n"
        "   parser.add_argument('-d', '--distrib', metavar='CHR', help='Plot BAF distribution of a single chromosome')\n"
        "   args = parser.parse_args()\n"
        "   dat = {}; fit = {}; cn = {}\n"
        "   read_dat(dat,fit,cn)\n"
        "   if args.distrib!=None:\n"
        "       plot_dist(dat,fit,args.distrib)\n"
        "   if args.all:\n"
        "       for chr in dat: plot_dist(dat,fit,chr)\n"
        "       plot_copy_number(cn)\n"
        "   elif args.copy_number:\n"
        "       plot_copy_number(cn)\n"
        "   else:\n"
        "       for chr in dat: plot_dist(dat,fit,chr)\n"
        "\n"
        "if __name__ == '__main__':\n"
        "   main()\n",
        args->output_dir);
//---------------------------------------
    chmod(fname, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH|S_IXUSR|S_IXGRP|S_IXOTH);
    if ( fclose(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,fname);
    free(fname);
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->ndist; i++)
    {
        free(args->dist[i].chr);
        free(args->dist[i].yvals);
    }
    free(args->dist);
    free(args->xvals);
    if ( fclose(args->dat_fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->dat_fname);
    free(args->dat_fname);
}

static void save_dist(args_t *args, dist_t *dist)
{
    int i;
    for (i=0; i<args->nbins; i++)
        fprintf(args->dat_fp,"DIST\t%s\t%f\t%f\n",dist->chr,dist->xvals[i],dist->yvals[i]);
}
static void fit_curves(args_t *args)
{
    peakfit_t *pkf = peakfit_init();
    peakfit_verbose(pkf,args->verbose);

    int i, nmc = 50;
    for (i=0; i<args->ndist; i++)
    {
        dist_t *dist = &args->dist[i];
        save_dist(args, &args->dist[i]);

        if ( dist->copy_number!=0 )
        {
            fprintf(args->dat_fp,"CN\t%s\t%.2f\n", dist->chr,(float)dist->copy_number);
            continue;
        }

        if ( args->verbose )
            fprintf(stderr,"%s:\n", dist->chr);

        int nrr_aa  = dist->iaa - dist->irr + 1;
        int nrr_ra  = dist->ira - dist->irr + 1;
        int naa_max = dist->nvals - dist->iaa;
        double xrr  = dist->xvals[dist->irr], *xrr_vals = &dist->xvals[dist->irr], *yrr_vals = &dist->yvals[dist->irr];
        double xaa  = dist->xvals[dist->iaa], *xaa_vals = &dist->xvals[dist->iaa], *yaa_vals = &dist->yvals[dist->iaa];
        double xra  = dist->xvals[dist->ira];
        double xmax = dist->xvals[dist->nvals-1];


        // CN2
        double cn2aa_fit = 0, cn2ra_fit, cn2_fit;
        char *cn2aa_func = 0, *cn2ra_func;
        double cn2aa_params[3] = {1,1,1} ,cn2ra_params[3];
        if ( args->include_aa )
        {
            peakfit_reset(pkf);
            peakfit_add_exp(pkf, 1.0,1.0,0.2, 5);
            peakfit_set_mc(pkf, 0.01,0.3,2,nmc);
            peakfit_set_mc(pkf, 0.05,1.0,0,nmc);
            cn2aa_fit  = peakfit_run(pkf, naa_max, xaa_vals, yaa_vals);
            cn2aa_func = strdup(peakfit_sprint_func(pkf));
            peakfit_get_params(pkf,0,cn2aa_params,3);
        }
        peakfit_reset(pkf);
        peakfit_add_bounded_gaussian(pkf, 1.0,0.5,0.03, 0.45,0.55, 7);
        peakfit_set_mc(pkf, 0.01,0.3,2,nmc);
        peakfit_set_mc(pkf, 0.05,1.0,0,nmc);
        cn2ra_fit  = peakfit_run(pkf, nrr_aa,xrr_vals,yrr_vals);
        cn2ra_func = strdup(peakfit_sprint_func(pkf));
        cn2_fit    = cn2ra_fit + cn2aa_fit;
        peakfit_get_params(pkf,0,cn2ra_params,3);


        // CN3: fit two peaks, then enforce the symmetry and fit again
        double cn3rra_params[5], cn3raa_params[5], *cn3aa_params = cn2aa_params;
        double cn3aa_fit = cn2aa_fit, cn3ra_fit;
        char *cn3aa_func = cn2aa_func, *cn3ra_func;
        double min_dx3   = 0.5 - 1./(args->min_fraction+2);
        peakfit_reset(pkf);
        peakfit_add_bounded_gaussian(pkf, 1.0,1/3.,0.03, xrr,xra-min_dx3, 7);
        peakfit_set_mc(pkf, xrr,xra-min_dx3, 1,nmc);
        peakfit_add_bounded_gaussian(pkf, 1.0,2/3.,0.03, xra+min_dx3,xaa, 7);
        peakfit_set_mc(pkf, xra+min_dx3,xaa, 1,nmc);
        peakfit_run(pkf, nrr_aa, xrr_vals, yrr_vals);
        // force symmetry around x=0.5
        peakfit_get_params(pkf,0,cn3rra_params,5);
        peakfit_get_params(pkf,1,cn3raa_params,5);
        double cn3_dx = (0.5-cn3rra_params[1] + cn3raa_params[1]-0.5)*0.5;
        if ( cn3_dx > 0.5/3 ) cn3_dx = 0.5/3;   // CN3 peaks should not be separated by more than 1/3
        peakfit_reset(pkf);
        peakfit_add_gaussian(pkf, cn3rra_params[0],0.5-cn3_dx,cn3rra_params[2], 5);
        peakfit_add_gaussian(pkf, cn3raa_params[0],0.5+cn3_dx,cn3raa_params[2], 5);
        cn3ra_fit  = peakfit_run(pkf, nrr_aa, xrr_vals, yrr_vals);
        cn3ra_func = strdup(peakfit_sprint_func(pkf));
        // compare peak sizes
        peakfit_get_params(pkf,0,cn3rra_params,3);
        peakfit_get_params(pkf,1,cn3raa_params,3);
        double cn3rra_size = cn3rra_params[0]*cn3rra_params[0];
        double cn3raa_size = cn3raa_params[0]*cn3raa_params[0];
        double cn3_dy      = cn3rra_size > cn3raa_size ? cn3raa_size/cn3rra_size : cn3rra_size/cn3raa_size;
        double cn3_frac    = (1 - 2*cn3rra_params[1]) / cn3rra_params[1];
        double cn3_fit     = cn3ra_fit + cn3aa_fit;
        // A very reasonable heuristics: check if the peak's width converged, exclude far too broad or far too narrow peaks
        if ( cn3rra_params[2]>0.3  || cn3raa_params[2]>0.3 ) cn3_fit = HUGE_VAL;
        if ( cn3rra_params[2]<1e-2 || cn3raa_params[2]<1e-2 ) cn3_fit = HUGE_VAL;


        // CN4 (contaminations)
        // - first fit only the [0,0.5] part of the data, then enforce the symmetry and fit again
        // - min_frac=1 (resp. 0.5) is interpreted as 50:50% (rep. 75:25%) contamination
        double cn4AAaa_params[3] = {1,1,1} ,cn4AAra_params[3] = {1,1,1}, cn4RAra_params[3], cn4RArr_params[5], cn4RAaa_params[5];
        double cn4aa_fit = 0, cn4ra_fit;
        char *cn4aa_func = 0, *cn4ra_func;
        double min_dx4   = 0.25*args->min_fraction;
        if ( args->include_aa )
        {
            peakfit_reset(pkf);
            peakfit_add_exp(pkf, 0.5,1.0,0.2, 5);
            peakfit_set_mc(pkf, 0.01,0.3,2,nmc);
            peakfit_add_bounded_gaussian(pkf, 0.4,(xaa+xmax)*0.5,2e-2, xaa,xmax, 7);
            peakfit_set_mc(pkf, xaa,xmax, 1,nmc);
            cn4aa_fit  = peakfit_run(pkf, naa_max, xaa_vals,yaa_vals);
            cn4aa_func = strdup(peakfit_sprint_func(pkf));
            peakfit_get_params(pkf,0,cn4AAaa_params,3);
            peakfit_get_params(pkf,1,cn4AAra_params,5);
        }
        peakfit_reset(pkf);
        // first fit only the [0,0.5] part of the data
        peakfit_add_gaussian(pkf, 1.0,0.5,0.03, 5);
        peakfit_add_bounded_gaussian(pkf, 0.6,0.3,0.03, xrr,xra-min_dx4, 7);
        peakfit_set_mc(pkf, xrr,xra-min_dx4,2,nmc);
        peakfit_run(pkf, nrr_ra , xrr_vals, yrr_vals);
        // now forcet symmetry around x=0.5
        peakfit_get_params(pkf,0,cn4RAra_params,3);
        peakfit_get_params(pkf,1,cn4RArr_params,5);
        double cn4_dx = 0.5-cn4RArr_params[1];
        if ( cn4_dx > 0.25 ) cn4_dx = 0.25;   // CN4 peaks should not be separated by more than 0.5
        peakfit_reset(pkf);
        peakfit_add_gaussian(pkf, cn4RAra_params[0],0.5,cn4RAra_params[2], 5);
        peakfit_add_gaussian(pkf, cn4RArr_params[0],0.5-cn4_dx,cn4RArr_params[2], 5);
        peakfit_add_gaussian(pkf, cn4RArr_params[0],0.5+cn4_dx,cn4RArr_params[2], 5);
        peakfit_set_mc(pkf, 0.1,cn4RAra_params[0],0,nmc);
        peakfit_set_mc(pkf, 0.01,0.1,2,nmc);
        cn4ra_fit  = peakfit_run(pkf, nrr_aa , xrr_vals, yrr_vals);
        cn4ra_func = strdup(peakfit_sprint_func(pkf));
        peakfit_get_params(pkf,0,cn4RAra_params,3);
        peakfit_get_params(pkf,1,cn4RArr_params,3);
        peakfit_get_params(pkf,2,cn4RAaa_params,3);
        double cn4RAra_size = cn4RAra_params[0]==0 ? HUGE_VAL : cn4RAra_params[0]*cn4RAra_params[0];
        double cn4RArr_size = cn4RArr_params[0]*cn4RArr_params[0];
        double cn4RAaa_size = cn4RAaa_params[0]*cn4RAaa_params[0];
        double cn4RArr_dy   = cn4RArr_size < cn4RAra_size ? cn4RArr_size/cn4RAra_size : cn4RAra_size/cn4RArr_size;
        double cn4RAaa_dy   = cn4RAaa_size < cn4RAra_size ? cn4RAaa_size/cn4RAra_size : cn4RAra_size/cn4RAaa_size;
        double cn4_dy       = cn4RArr_dy < cn4RAaa_dy ? cn4RArr_dy/cn4RAaa_dy : cn4RAaa_dy/cn4RArr_dy;
        double cn4_ymin     = cn4RArr_size < cn4RAaa_size ? cn4RArr_size/cn4RAra_size : cn4RAaa_size/cn4RAra_size;
        cn4_dx              = (cn4RAaa_params[1]-0.5) - (0.5-cn4RArr_params[1]);
        double cn4_frac     = cn4RAaa_params[1] - cn4RArr_params[1];
        double cn4_fit      = cn4ra_fit + cn4aa_fit;
        // A very reasonable heuristics: check if the peak's width converged, exclude far too broad or far too narrow peaks
        if ( cn4RAra_params[2]>0.3 || cn4RArr_params[2]>0.3 || cn4RAaa_params[2]>0.3 ) cn4_fit = HUGE_VAL;
        if ( cn4RAra_params[2]<1e-2 || cn4RArr_params[2]<1e-2 || cn4RAaa_params[2]<1e-2 ) cn4_fit = HUGE_VAL;


        // Choose the best match
        char cn2_fail = '*', cn3_fail = '*', cn4_fail = '*';
        if ( cn2_fit > args->fit_th ) cn2_fail = 'f';

        if ( cn3_fit > args->fit_th ) cn3_fail = 'f';
        else if ( cn3_dy < args->peak_symmetry ) cn3_fail = 'y';    // size difference is too big

        if ( cn4_fit > args->fit_th ) cn4_fail = 'f';
        else if ( cn4_ymin < args->min_peak_size ) cn4_fail = 'y';      // side peak is too small
        else if ( cn4_dy < args->peak_symmetry ) cn4_fail = 'Y';    // size difference is too big
        else if ( cn4_dx > 0.1 ) cn4_fail = 'x';                    // side peaks placed assymetrically

        double cn = -1, fit = cn2_fit;
        if ( cn2_fail == '*' ) { cn = 2; fit = cn2_fit; }
        if ( cn3_fail == '*' )
        {
            // Use cn_penalty as a tiebreaker. If set to 0.3, cn3_fit must be 30% smaller than cn2_fit.
            if ( cn<0 || cn3_fit < (1-args->cn_penalty) * fit )
            {
                cn = 2 + cn3_frac; 
                fit = cn3_fit; 
                if ( cn2_fail=='*' ) cn2_fail = 'p';
            }
            else cn3_fail = 'p';
        }
        if ( cn4_fail == '*' )
        {
            if ( cn<0 || cn4_fit < (1-args->cn_penalty) * fit )
            {
                cn = 3 + cn4_frac;
                fit = cn4_fit;
                if ( cn2_fail=='*' ) cn2_fail = 'p';
                if ( cn3_fail=='*' ) cn3_fail = 'p';
            }
            else cn4_fail = 'p';
        }

        if ( args->verbose )
        {
            fprintf(stderr,"\tcn2 %c fit=%e\n", cn2_fail, cn2_fit);
            fprintf(stderr,"\t       .. %e\n", cn2ra_fit);
            fprintf(stderr,"\t            RA:   %f %f %f\n", cn2ra_params[0],cn2ra_params[1],cn2ra_params[2]);
            fprintf(stderr,"\t       .. %e\n", cn2aa_fit);
            fprintf(stderr,"\t            AA:   %f %f %f\n", cn2aa_params[0],cn2aa_params[1],cn2aa_params[2]);
            fprintf(stderr,"\t      func:\n");
            fprintf(stderr,"\t            %s\n", cn2ra_func);
            fprintf(stderr,"\t            %s\n", cn2aa_func);
            fprintf(stderr,"\n");
            fprintf(stderr,"\tcn3 %c fit=%e  frac=%f  symmetry=%f\n", cn3_fail, cn3_fit, cn3_frac, cn3_dy);
            fprintf(stderr,"\t       .. %e\n", cn3ra_fit);
            fprintf(stderr,"\t            RRA:  %f %f %f\n", cn3rra_params[0],cn3rra_params[1],cn3rra_params[2]);
            fprintf(stderr,"\t            RAA:  %f %f %f\n", cn3raa_params[0],cn3raa_params[1],cn3raa_params[2]);
            fprintf(stderr,"\t       .. %e\n", cn3aa_fit);
            fprintf(stderr,"\t            AAA:  %f %f %f\n", cn3aa_params[0],cn3aa_params[1],cn3aa_params[2]);
            fprintf(stderr,"\t      func:\n");
            fprintf(stderr,"\t            %s\n", cn3ra_func);
            fprintf(stderr,"\t            %s\n", cn3aa_func);
            fprintf(stderr,"\n");
            fprintf(stderr,"\tcn4 %c fit=%e  frac=%f  symmetry=%f ymin=%f\n", cn4_fail, cn4_fit, cn4_frac, cn4_dy, cn4_ymin);
            fprintf(stderr,"\t       .. %e\n", cn4ra_fit);
            fprintf(stderr,"\t            RArr:  %f %f %f\n", cn4RArr_params[0],cn4RArr_params[1],cn4RArr_params[2]);
            fprintf(stderr,"\t            RAra:  %f %f %f\n", cn4RAra_params[0],cn4RAra_params[1],cn4RAra_params[2]);
            fprintf(stderr,"\t            RAaa:  %f %f %f\n", cn4RAaa_params[0],cn4RAaa_params[1],cn4RAaa_params[2]);
            fprintf(stderr,"\t       .. %e\n", cn4aa_fit);
            fprintf(stderr,"\t            AAaa:  %f %f %f\n", cn4AAaa_params[0],cn4AAaa_params[1],cn4AAaa_params[2]);
            fprintf(stderr,"\t      func:\n");
            fprintf(stderr,"\t            %s\n", cn4ra_func);
            fprintf(stderr,"\t            %s\n", cn4aa_func);
            fprintf(stderr,"\n");
        }

        if ( args->force_cn==2 || cn2_fail == '*' )
        {
            fprintf(args->dat_fp,"FIT\t%s\t%e\t%d\t%d\t%s\n", dist->chr,cn2ra_fit,dist->irr,dist->iaa,cn2ra_func);
            if ( cn2aa_func ) fprintf(args->dat_fp,"FIT\t%s\t%e\t%d\t%d\t%s\n", dist->chr,cn2aa_fit,dist->iaa,dist->nvals-1,cn2aa_func);
        }
        if ( args->force_cn==3 || cn3_fail == '*' )
        {
            fprintf(args->dat_fp,"FIT\t%s\t%e\t%d\t%d\t%s\n", dist->chr,cn3ra_fit,dist->irr,dist->iaa,cn3ra_func);
            if ( cn3aa_func ) fprintf(args->dat_fp,"FIT\t%s\t%e\t%d\t%d\t%s\n", dist->chr,cn3aa_fit,dist->iaa,dist->nvals-1,cn3aa_func);
        }
        if ( args->force_cn==4 || cn4_fail == '*' )
        {
            fprintf(args->dat_fp,"FIT\t%s\t%e\t%d\t%d\t%s\n", dist->chr,cn4ra_fit,dist->irr,dist->iaa,cn4ra_func);
            if ( cn4aa_func ) fprintf(args->dat_fp,"FIT\t%s\t%e\t%d\t%d\t%s\n", dist->chr,cn4aa_fit,dist->iaa,dist->nvals-1,cn4aa_func);
        }
        fprintf(args->dat_fp,"CN\t%s\t%.2f\t%f\n", dist->chr, cn, fit);

        free(cn2aa_func);
        free(cn2ra_func);
        free(cn3ra_func);
        free(cn4ra_func);
        free(cn4aa_func);
    }

    peakfit_destroy(pkf);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Detect number of chromosomal copies from Illumina's B-allele frequency (BAF)\n");
    fprintf(stderr, "Usage:   bcftools polysomy [OPTIONS] FILE.vcf\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr, "    -o, --output-dir PATH          \n");
    fprintf(stderr, "    -r, --regions REGION           Restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file FILE        Restrict to regions listed in a file\n");
    fprintf(stderr, "        --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(stderr, "    -s, --sample NAME              Sample to analyze\n");
    fprintf(stderr, "    -t, --targets REGION           Similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file FILE        Similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "        --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n");
    fprintf(stderr, "    -v, --verbose                  \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Algorithm options:\n");
    fprintf(stderr, "    -b, --peak-size FLOAT          Minimum peak size (0-1, larger is stricter) [0.1]\n");
    fprintf(stderr, "    -c, --cn-penalty FLOA>         Penalty for increasing CN (0-1, larger is stricter) [0.7]\n");
    fprintf(stderr, "    -f, --fit-th FLOAT             Goodness of fit threshold (>0, smaller is stricter) [3.3]\n");
    fprintf(stderr, "    -i, --include-aa               Include the AA peak in CN2 and CN3 evaluation\n");
    fprintf(stderr, "    -m, --min-fraction FLOAT       Minimum distinguishable fraction of aberrant cells [0.1]\n");
    fprintf(stderr, "    -p, --peak-symmetry FLOAT      Peak symmetry threshold (0-1, larger is stricter) [0.5]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_polysomy(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->nbins  = 150;
    args->fit_th = 3.3;
    args->cn_penalty = 0.7;
    args->peak_symmetry = 0.5;
    args->min_peak_size = 0.1;
    args->ra_rr_scaling = 1;
    args->min_fraction = 0.1;
    args->smooth = -3;
    args->regions_overlap = 1;
    args->targets_overlap = 0;

    static struct option loptions[] =
    {
        {"ra-rr-scaling",0,0,1},    // hidden option
        {"force-cn",1,0,2},         // hidden option
        {"smooth",1,0,'S'},         // hidden option
        {"nbins",1,0,'n'},          // hidden option
        {"include-aa",0,0,'i'},
        {"peak-size",1,0,'b'},
        {"min-fraction",1,0,'m'},
        {"verbose",0,0,'v'},
        {"fit-th",1,0,'f'},
        {"cn-penalty",1,0,'c'},
        {"peak-symmetry",1,0,'p'},
        {"output-dir",1,0,'o'},
        {"sample",1,0,'s'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,4},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {0,0,0,0}
    };
    char *tmp;
    int c;
    while ((c = getopt_long(argc, argv, "h?o:vt:T:r:R:s:f:p:c:im:b:n:S:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->ra_rr_scaling = 0; break;
            case  2 : args->force_cn = atoi(optarg); break;
            case  3 :
                args->regions_overlap = parse_overlap_option(optarg);
                if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  4 :
                args->targets_overlap = parse_overlap_option(optarg);
                if ( args->targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case 'n': args->nbins = atoi(optarg); break;
            case 'S': args->smooth = atoi(optarg); break;
            case 'i': args->include_aa = 1; break;
            case 'b':
                args->min_peak_size = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -b %s\n", optarg);
                if ( args->min_peak_size<0 || args->min_peak_size>1 ) error("Range error: -b %s\n", optarg);
                break;
            case 'm':
                args->min_fraction = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -n %s\n", optarg);
                if ( args->min_fraction<0 || args->min_fraction>1 ) error("Range error: -n %s\n", optarg);
                break;
            case 'f':
                args->fit_th = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -f %s\n", optarg);
                break;
            case 'p':
                args->peak_symmetry = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -p %s\n", optarg);
                break;
            case 'c':
                args->cn_penalty = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -c %s\n", optarg);
                break;
            case 's': args->sample = optarg; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; args->targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; args->regions_is_file = 1; break;
            case 'o': args->output_dir = optarg; break;
            case 'v': args->verbose++; break;
            default: usage(args); break;
        }
    }
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";
    }
    else args->fname = argv[optind];
    if ( !args->fname ) usage(args);
    if ( !args->output_dir ) error("Missing the -o option\n");

    init_data(args);
    fit_curves(args);
    destroy_data(args);
    free(args);

    return 0;
}


