/* The MIT License

   Copyright (c) 2013-2014 Genome Research Ltd.

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
#include "bcftools.h"

typedef struct
{
    int nvals;      // number of data points to fit against (excluding RR,AA peaks)
    double *xvals;  // xvalues, pointer to dist_t.xvals
    double *yvals;  // yvalues, pointer to dist_t.yvals
    int ngauss;     // number of gaussian functions
}
data_t;

typedef struct
{
    data_t dat;
    int nvals;          // all values, including RR,AA peaks
    double *xvals;      // pointer to args_t.xvals
    double *yvals;   
    int copy_number;    // heuristics to skip futile CN1 fits when no het peak is detected
    int irr, iaa;       // chop off RR and AA peaks
    char *chr;
}
dist_t;

typedef struct
{
    int ndist, nbins;
    double *xvals;
    dist_t *dist;
    char **argv, *output_dir;
    double fit_th, peak_symmetry, cn_penalty;
    int argc, plot, verbose, regions_is_file, targets_is_file;
    char *dat_fname, *fname, *regions_list, *targets_list, *sample;
    FILE *dat_fp;
}
args_t;

FILE *open_file(char **fname, const char *mode, const char *fmt, ...);

static void init_dist(dist_t *dist, int verbose)
{
    // isolate RR and AA peaks and rescale so that they are comparable to hets
    int i, irr, iaa, n = dist->nvals;

    // smooth the distribution
    double *tmp = (double*) malloc(sizeof(double)*n);
    int win = 0.02*n < 1 ? 1 : 0.02*n;
    double avg = 0;
    for (i=0; i<win; i++) avg += dist->yvals[i];
    for (i=0; i<n-win; i++)
    {
        tmp[i] = avg / win;
        avg += -dist->yvals[i] + dist->yvals[i+win];
    }
    for (; i<n; i++) tmp[i] = avg;

    // find the extremes; first a simple approach: find a gap
    for (irr=0,i=0; i<n/2; i++) if ( tmp[i] < tmp[irr] ) irr = i;
    for (iaa=n-1,i=n-1; i>=n/2; i--) if ( tmp[i] < tmp[iaa] ) iaa = i;
    irr += win*0.5;
    iaa += win*0.5;
    if ( irr>=iaa ) error("FIXME: oops, dist normalization failed for %s: %d vs %d\n", dist->chr,irr,iaa); // we may need to be smarter
    free(tmp);

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

    if ( !sra || (sra/srr<0.1 && saa/sra>1.0) )
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

    dist->dat.yvals = &dist->yvals[irr];
    dist->dat.xvals = &dist->xvals[irr];
    dist->dat.nvals = iaa - irr + 1;

    if ( verbose )
        fprintf(stderr,"%s:\t cn=%2d \t ra/rr=%f \t aa/ra=%f \t nra=%d\n", dist->chr,dist->copy_number,sra/srr,saa/sra, (int)sra);
}

static void init_data(args_t *args)
{
    bcf_srs_t *files = bcf_sr_init();
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(files, args->regions_list, args->regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(files, args->targets_list, args->targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( !bcf_sr_add_reader(files, args->fname) ) error("Failed to open %s: %s\n", args->fname,bcf_sr_strerror(files->errnum));
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
        init_dist(&args->dist[idist],args->verbose);

    args->dat_fp = open_file(&args->dat_fname,"w","%s/dist.dat", args->output_dir);
    fprintf(args->dat_fp, "# This file was produced by: bcftools cnv(%s+htslib-%s), the command line was:\n", bcftools_version(),hts_version());
    fprintf(args->dat_fp, "# \t bcftools %s ", args->argv[0]);
    for (i=1; i<args->argc; i++)
        fprintf(args->dat_fp, " %s",args->argv[i]);
    fprintf(args->dat_fp,"\n#\n");
    fprintf(args->dat_fp,"# DIST\t[2]Chrom\t[3]BAF\t[4]Normalized Count\n");
    fprintf(args->dat_fp,"# FIT\t[2]Chrom\t[3]Mean of fitted Gaussian\t[4]Scale\t[5]Sigma[6]\tMean etc.\n");
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
        "import csv,math,sys,argparse\n"
        "\n"
        "outdir = '%s'\n"
        "\n"
        "def read_dat(dat,fit,cn):\n"
        "   csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)\n"
        "   with open(outdir+'/dist.dat', 'rb') as f:\n"
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
        "              fit[chr] = row[2:]\n"
        "          elif type=='CN':\n"
        "              cn[chr] = row[2]\n"
        "\n"
        "def fitted_func(xvals,params):\n"
        "   n = len(params)/3\n"
        "   out = []\n"
        "   for x in xvals:\n"
        "       y = 0\n"
        "       for i in range(n):\n"
        "           mean  = float(params[i*3+0])\n"
        "           scale = float(params[i*3+1])\n"
        "           sigma = float(params[i*3+2])\n"
        "           y += scale * math.exp(-(float(x)-mean)**2/sigma**2)\n"
        "       out.append(y)\n"
        "   return out\n"
        "\n"
        "def plot_dist(dat,fit,chr):\n"
        "   fig, ax = plt.subplots(1, 1, figsize=(7,5))\n"
        "   ax.plot([x[2] for x in dat[chr]],[x[3] for x in dat[chr]],'-',label='Distribution')\n"
        "   ax.plot([x[2] for x in dat[chr]],fitted_func([x[2] for x in dat[chr]], fit[chr]),'-',label='Best Fit')\n"
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
        "def main():\n"
        "   parser = argparse.ArgumentParser()\n"
        "   parser.add_argument('-a', '--all', action='store_true', help='Create all plots')\n"
        "   parser.add_argument('-c', '--copy-number', action='store_true', help='Create copy-number plot')\n"
        "   parser.add_argument('-d', '--distrib', metavar='CHR', help='Plot BAF distribution of a single chromosome')\n"
        "   args = parser.parse_args()\n"
        "   if args.distrib==None and not args.all and not args.copy_number: parser.print_help()\n"
        "   dat = {}; fit = {}; cn = {}\n"
        "   read_dat(dat,fit,cn)\n"
        "   if args.distrib!=None:\n"
        "       plot_dist(dat,fit,args.distrib)\n"
        "   if args.all:\n"
        "       for chr in dat: plot_dist(dat,fit,chr)\n"
        "       plot_copy_number(cn)\n"
        "   elif args.copy_number:\n"
        "       plot_copy_number(cn)\n"
        "\n"
        "if __name__ == '__main__':\n"
        "   main()\n",
        args->output_dir);
//---------------------------------------
    chmod(fname, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH|S_IXUSR|S_IXGRP|S_IXOTH);
    free(fname);
    fclose(fp);
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
    free(args->dat_fname);
    fclose(args->dat_fp);
}

static void save_dist(args_t *args, int idist, int ngauss, double *params)
{
    int i;
    for (i=0; i<args->nbins; i++)
        fprintf(args->dat_fp,"DIST\t%s\t%f\t%f\n",args->dist[idist].chr,args->dist[idist].xvals[i],args->dist[idist].yvals[i]);
    fprintf(args->dat_fp,"FIT\t%s", args->dist[idist].chr);
    for (i=0; i<ngauss*3; i++) fprintf(args->dat_fp,"\t%f", params[i]);
    fprintf(args->dat_fp,"\n");
}

int func_f(const gsl_vector *params, void *data, gsl_vector *yvals)
{
    data_t *dat = (data_t *) data;

    int i, j;
    for (i=0; i<dat->nvals; i++)
    {
        double xi = dat->xvals[i];
        double yi = 0;
        for (j=0; j<dat->ngauss; j++)
        {
            double center = gsl_vector_get(params,j*3 + 0);
            double scale  = gsl_vector_get(params,j*3 + 1);
            double sigma  = gsl_vector_get(params,j*3 + 2);

            double zi = (xi - center) / sigma;
            yi += scale*scale * exp(-zi*zi);
        }
        gsl_vector_set(yvals, i, (yi - dat->yvals[i])/0.1);
    }
    return GSL_SUCCESS;
}

int func_df(const gsl_vector *params, void *data, gsl_matrix *jacobian)
{
    data_t *dat = (data_t *) data;

    int i, j;
    for (i=0; i<dat->nvals; i++)
    {
        // Jacobian matrix J(i,j) = dfi / dxj,
        // where fi = (Yi - yi),
        //       Yi = scale^2 * exp(-(center - xi)^2/sigma^2)
        //

        double xi = dat->xvals[i];
        for (j=0; j<dat->ngauss; j++)
        {
            double center = gsl_vector_get(params,j*3 + 0);
            double scale  = gsl_vector_get(params,j*3 + 1);
            double sigma  = gsl_vector_get(params,j*3 + 2);

            double zi = (xi - center) / sigma;
            double ei = exp(-zi*zi);

            gsl_matrix_set(jacobian, i, j*3 + 0, 2*scale*scale*(xi-center)/(sigma*sigma)*ei);
            gsl_matrix_set(jacobian, i, j*3 + 1, 2*scale*ei);
            gsl_matrix_set(jacobian, i, j*3 + 2, 2*scale*scale*(xi-center)*(xi-center)/(sigma*sigma*sigma)*ei);
        }
    }
    return GSL_SUCCESS;
}

int func_set(const gsl_vector *params, void *data, gsl_vector *yvals, gsl_matrix *jacobian)
{
    func_f(params, data, yvals);
    func_df(params, data, jacobian);
    return GSL_SUCCESS;
}
static double eval_fit(int nvals, double *xvals, double *yvals, int ngauss, double *params)
{
    double sum = 0;
    int i, j;
    for (i=0; i<nvals; i++)
    {
        double yval = 0;
        for (j=0; j<ngauss; j++)
        {
            double center = params[j*3 + 0];
            double scale  = params[j*3 + 1];
            double sigma  = params[j*3 + 2];
            
            double zi = (xvals[i] - center) / sigma;
            yval += scale * exp(-zi*zi);
        }
        sum += fabs(yval - yvals[i]);
    }
    return sum;
}

static int gauss_fit(dist_t *dist, int ngauss, double *params)
{
    data_t *dat = &dist->dat;
    dat->ngauss = ngauss;

    gsl_multifit_function_fdf mfunc;
    mfunc.f   = &func_f;
    mfunc.df  = &func_df;
    mfunc.fdf = &func_set;
    mfunc.n   = dat->nvals;
    mfunc.p   = ngauss*3;      // number of fitting parameters
    mfunc.params = dat;

    const gsl_multifit_fdfsolver_type *solver_type;
    gsl_multifit_fdfsolver *solver;
    gsl_vector_view vview = gsl_vector_view_array(params, mfunc.p);
    solver_type = gsl_multifit_fdfsolver_lmsder;
    solver = gsl_multifit_fdfsolver_alloc(solver_type, dat->nvals, mfunc.p);
    gsl_multifit_fdfsolver_set(solver, &mfunc, &vview.vector);

    int i, status;
    size_t iter = 0;
    do
    {
        status = gsl_multifit_fdfsolver_iterate(solver);
        if ( status ) break;
        status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-4, 1e-4);
    }
    while (status == GSL_CONTINUE && iter++ < 500);

    for (i=0; i<mfunc.p; i++)
        params[i] = gsl_vector_get(solver->x, i);

    gsl_multifit_fdfsolver_free(solver);
    return iter>500 ? -1 : 0;
}

static double best_fit(args_t *args, dist_t *dist, int ngauss, double *params)
{
    if ( ngauss==1 )
    {
        gauss_fit(dist,ngauss,params);
        params[1] *= params[1];
        return eval_fit(dist->dat.nvals, dist->dat.xvals, dist->dat.yvals, ngauss,params);
    }

    int i, j, n = 3;
    int ipk = 3*(ngauss-1);
    double delta = 0.5 * (params[ipk] - params[0]) / n;
    double best_params[9], tmp_params[9], best_fit = HUGE_VAL;
    for (i=0; i<n; i++)
    {
        memcpy(tmp_params,params,sizeof(double)*ngauss*3);
        tmp_params[0]   += delta*i;
        tmp_params[ipk] -= delta*i;
        if ( gauss_fit(dist,ngauss,tmp_params)<0 ) continue;    // did not converge

        for (j=0; j<ngauss; j++) tmp_params[j*3+1] *= tmp_params[j*3+1];

        // From the nature of the data, we can assume that in presence of
        // multiple peaks they will be placed symmetrically around 0.5. Also
        // their size should be about the same. We evaluate the fit with this
        // in mind.
        double dx = fabs(0.5 - tmp_params[0]) + fabs(tmp_params[ipk] - 0.5);
        tmp_params[0] = 0.5 - dx*0.5;
        tmp_params[ipk] = 0.5 + dx*0.5;
        double fit = eval_fit(dist->dat.nvals, dist->dat.xvals, dist->dat.yvals, ngauss, tmp_params);

        if ( best_fit < fit ) continue;     // worse than previous
        best_fit = fit;
        memcpy(best_params,tmp_params,sizeof(double)*ngauss*3);
    }
    memcpy(params,best_params,sizeof(double)*ngauss*3);
    return best_fit;
}

static void print_params(data_t *dat, int ngauss, double *params, float fit, float frac, char fail, char comment)
{
    int i, j;
    printf("\t%c%c fit=%f frac=%.2f .. center,scale,sigma = ", comment,fail?fail:'o',fit,frac);
    for (i=0; i<ngauss; i++)
    {
        if ( i!=0 ) printf("\t");
        for (j=0; j<3; j++) printf(" %f", params[i*3+j]);
    }
    printf("\n");
}
static void fit_curves(args_t *args)
{
    int i;
    for (i=0; i<args->ndist; i++)
    {
        dist_t *dist = &args->dist[i];
        if ( dist->copy_number!=0 )
        {
            fprintf(args->dat_fp,"CN\t%s\t%.2f\n", dist->chr,(float)dist->copy_number);
            save_dist(args, i, 0, NULL);
            continue;
        }

        // Parameters (center,scale,sigma) for gaussian peaks.
        double params_cn2[] = { 1/2.,0.5,0.05 };
        double params_cn3[] = { 1/3.,0.5,0.05, 2/3.,0.5,0.05 };
        double params_cn4[] = { 1/4.,0.5,0.05, 1/2.,0.5,0.05, 3/4.,0.5,0.05 };

        double fit_cn2 = best_fit(args,&args->dist[i],1,params_cn2);
        double fit_cn3 = best_fit(args,&args->dist[i],2,params_cn3);
        double fit_cn4 = best_fit(args,&args->dist[i],3,params_cn4);

        double dx_cn3  = fabs(params_cn3[0] - params_cn3[3]);
        double dx_cn4  = fabs(params_cn4[0] - params_cn4[6]);
        double dy_cn3  = params_cn3[1] > params_cn3[4] ? params_cn3[4]/params_cn3[1] : params_cn3[1]/params_cn3[4];
        double dy_cn4a = params_cn4[1] > params_cn4[7] ? params_cn4[7]/params_cn4[1] : params_cn4[1]/params_cn4[7]; // side peaks
        double ymax = params_cn4[1] > params_cn4[7] ? params_cn4[1] : params_cn4[7];
        double dy_cn4b = ymax > params_cn4[4] ? params_cn4[4]/ymax : ymax/params_cn4[4];    // middle peak

        // Three peaks (CN4) are always a better fit than two (CN3) or one (CN2). Therefore
        // check that peaks are well separated and that the peak sizes are reasonable
        char cn2_fail = 0, cn3_fail = 0, cn4_fail = 0;
        if ( fit_cn2 > args->fit_th ) cn2_fail = 'f';

        if ( fit_cn3 > args->fit_th ) cn3_fail = 'f';
        else if ( dx_cn3 < 0.05 ) cn3_fail = 'x';         // peak separation: at least ~10% of cells
        else if ( dy_cn3 < args->peak_symmetry ) cn3_fail = 'y';    

        if ( fit_cn4 > args->fit_th ) cn4_fail = 'f';
        else if ( dx_cn4 < 0.1 ) cn4_fail = 'x';            // peak separation
        else if ( dy_cn4a < args->peak_symmetry ) cn4_fail = 'y';
        else if ( dy_cn4b < args->peak_symmetry ) cn4_fail = 'Y';

        // Estimate fraction of affected cells. For CN4 we estimate
        // contamination (the fraction of foreign cells), which is more
        // common than CN4; hence the value is from the interval [0,0.5].
        //      CN3 .. f = 2*dx/(1-dx)
        //      CN4 .. f = dx
        dx_cn3 = 2*dx_cn3 / (1-dx_cn3);

        double cn = -1, fit = fit_cn2;
        if ( !cn2_fail ) { cn = 2; fit = fit_cn2; }
        if ( !cn3_fail && fit_cn3 < args->cn_penalty * fit ) { cn = 3; fit = fit_cn3; }
        if ( !cn4_fail && fit_cn4 < args->cn_penalty * fit ) { cn = 4; fit = fit_cn4; }

        if ( cn==-1 ) save_dist(args, i, 0, NULL);
        else if ( cn==2 ) save_dist(args, i, 1, params_cn2);
        else if ( cn==3 ) { save_dist(args, i, 2, params_cn3); cn = 2 + dx_cn3; }
        else if ( cn==4 ) { save_dist(args, i, 3, params_cn4); cn = 3 + dx_cn4; }

        if ( args->verbose )
        {
            printf("%s: \n", args->dist[i].chr);
            print_params(&args->dist[i].dat, 1, params_cn2, fit_cn2, 1.0,    cn2_fail, cn==2 ? '*' : ' ');
            print_params(&args->dist[i].dat, 2, params_cn3, fit_cn3, dx_cn3, cn3_fail, cn>2 && cn<=3 ? '*' : ' ');
            print_params(&args->dist[i].dat, 3, params_cn4, fit_cn4, dx_cn4, cn4_fail, cn>3 ? '*' : ' ');
            printf("\n");
        }
        fprintf(args->dat_fp,"CN\t%s\t%.2f\t%f\n", dist->chr, cn, fit);
    }
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Detect number of chromosomal copies from Illumina's B-allele frequency (BAF)\n");
    fprintf(stderr, "Usage:   bcftools polysomy [OPTIONS] <file.vcf>\n");
    fprintf(stderr, "General options:\n");
    fprintf(stderr, "    -o, --output-dir <path>        \n");
    fprintf(stderr, "    -r, --regions <region>         restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>      restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --sample <name>            sample to analyze\n");
    fprintf(stderr, "    -t, --targets <region>         similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>      similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "    -v, --verbose                  \n");
    fprintf(stderr, "Algorithm options:\n");
    fprintf(stderr, "    -c, --cn-penalty <float>       penalty for increasing CN (smaller more strict) [0.7]\n");
    fprintf(stderr, "    -f, --fit-th <float>           goodness of fit threshold (smaller more strict) [3.0]\n");
    fprintf(stderr, "    -p, --peak-symmetry <float>    peak symmetry threshold (bigger more strict) [0.7]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_polysomy(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->nbins  = 150;
    args->fit_th = 3.0;
    args->cn_penalty = 0.7;
    args->peak_symmetry = 0.7;

    static struct option loptions[] = 
    {
        {"verbose",0,0,'v'},
        {"fit-th",1,0,'f'},
        {"cn-penalty",1,0,'c'},
        {"peak-symmetry",1,0,'p'},
        {"output-dir",1,0,'o'},
        {"sample",1,0,'s'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {0,0,0,0}
    };
    char c, *tmp;
    while ((c = getopt_long(argc, argv, "h?o:vt:T:r:R:s:f:p:c:",loptions,NULL)) >= 0) 
    {
        switch (c) 
        {
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
            case 'v': args->verbose = 1; break;
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


