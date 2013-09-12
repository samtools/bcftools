#include <math.h>
#include <htslib/kfunc.h>
#include "call.h"

int ccall(call_t *call, bcf1_t *rec) { return 0; }
int qcall(call_t *call, bcf1_t *rec) 
{ 
    // QCall format: 
    //  chromosome, position, reference allele, depth, mapping quality, 0, ..
    error("TODO: qcall output\n");
    return 0; 
}

void mcall_init(call_t *call) 
{ 
    call->pl2p = (double*) malloc(sizeof(double)*256);
    int i;
    for (i=0; i<256; i++)
        call->pl2p[i] = pow(10., -i/10.);

    call->nqsum = 4;
    call->qsum  = (float*) malloc(sizeof(float)*call->nqsum); 
    call->nals_map = 4;
    call->als_map  = (int*) malloc(sizeof(int)*call->nals_map);
    call->npl_map  = 4*(4+1)/2;
    call->pl_map   = (int*) malloc(sizeof(int)*call->npl_map);
    call->gts = (int*) calloc(call->hdr->n[BCF_DT_SAMPLE]*2,sizeof(int));   // assuming at most diploid everywhere

    bcf_hdr_append(call->hdr,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=ICB,Number=1,Type=Float,Description=\"Inbreeding Coefficient Binomial test (bigger is better)\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">");

    return; 
}
void ccall_init(call_t *call) { return; }
void qcall_init(call_t *call) { return; }

void mcall_destroy(call_t *call) 
{ 
    free(call->pl2p);
    free(call->PLs);
    free(call->qsum);
    free(call->als_map);
    free(call->pl_map);
    free(call->gts);
    free(call->pdg);
    free(call->als);
    return; 
}
void ccall_destroy(call_t *call) { return; }
void qcall_destroy(call_t *call) { return; }


// Inits P(D|G): convert PLs from log space and normalize. In case of zero
// depth, missing PLs are all zero. In this case, pdg's are set to 0
// so that the corresponding genotypes can be set as missing and the
// qual calculation is not affected.
void set_pdg(call_t *call, int *PLs, double *pdg, int n_smpl, int n_gt)
{
    int i, j;
    for (i=0; i<n_smpl; i++)
    {
        double sum = 0;
        for (j=0; j<n_gt; j++)
        {
            assert( PLs[j]<256 );
            pdg[j] = call->pl2p[ PLs[j] ];
            sum += pdg[j];
        }
        // Normalize: sum_i pdg_i = 1
        if ( sum!=n_gt )
            for (j=0; j<n_gt; j++) pdg[j] /= sum;
        else
            for (j=0; j<n_gt; j++) pdg[j] = 0;

        PLs += n_gt;
        pdg += n_gt;
    }
}

// Create mapping between old and new (trimmed) alleles
void init_allele_trimming_maps(call_t *call, int als, int nals)
{
    int i, j;

    // als_map: old(i) -> new(j)
    for (i=0, j=0; i<nals; i++)
    {
        if ( als & 1<<i ) call->als_map[i] = j++;
        else call->als_map[i] = -1;
    }

    // pl_map: new(k) -> old(l)
    int k = 0, l = 0;
    for (i=0; i<nals; i++)
    {
        for (j=0; j<=i; j++)
        {
            if ( (als & 1<<i) && (als & 1<<j) ) call->pl_map[k++] = l;
            l++;
        }
    }
}

double binom_dist(int N, double p, int k)
{
    int mean = (int) (N*p);
    if ( mean==k ) return 1.0;

    double log_p = (k-mean)*log(p) + (mean-k)*log(1.0-p);
    if ( k > N - k ) k = N - k;
    if ( mean > N - mean ) mean = N - mean;
    
    if ( k < mean ) { int tmp = k; k = mean; mean = tmp; }
    double diff = k - mean;

    double val = 1.0;
    int i;
    for (i=0; i<diff; i++)
        val = val * (N-mean-i) / (k-i);

    return exp(log_p)/val;
}


// Inbreeding Coefficient, binomial test
float calc_ICB(int nref, int nalt, int nhets, int ndiploid)
{
    if ( !nref || !nalt || !ndiploid ) return 1.0;

    double fref = (double)nref/(nref+nalt); // fraction of reference allelels
    double falt = (double)nalt/(nref+nalt); // non-ref als
    double q = 2*fref*falt;                 // probability of a het, assuming HWE
    double mean = q*ndiploid;

    //fprintf(stderr,"\np=%e N=%d k=%d  .. nref=%d nalt=%d nhets=%d ndiploid=%d\n", q,ndiploid,nhets, nref,nalt,nhets,ndiploid);

    // Can we use normal approximation? The second condition is for performance only
    // and is not well justified. 
    if ( (mean>10 && (1-q)*ndiploid>10 ) || ndiploid>200 )
        return exp(-0.5*(nhets-mean)*(nhets-mean)/(mean*(1-q)));

    return binom_dist(ndiploid, q, nhets);
}

/**
  *  log(sum_i exp(a_i))
  */
inline double logsumexp(double *vals, int nvals)
{
    int i;
    double max_exp = vals[0];
    for (i=1; i<nvals; i++)
        if ( max_exp < vals[i] ) max_exp = vals[i];

    double sum = 0;
    for (i=0; i<nvals; i++)
        sum += exp(vals[i] - max_exp);

    return log(sum) + max_exp;
}
inline double logsumexp2(double a, double b)
{
    if ( a>b )
        return log(1 + exp(b-a)) + a;
    else
        return log(1 + exp(a-b)) + b;
}

// Macro to set the most likely and second most likely alleles
#define UPDATE_MAX_LKs(als) { \
     if ( max_lk<lk_tot ) { max_lk2 = max_lk; max_als2 = max_als; max_lk = lk_tot; max_als = (als); } \
     else if ( max_lk2<lk_tot ) { max_lk2 = lk_tot; max_als2 = (als); } \
     if ( lk_tot_set ) lk_sum = logsumexp2(lk_tot,lk_sum); \
}

#define SWAP(type_t,x,y) {type_t tmp; tmp = x; x = y; y = tmp; }


/**
  *  This function implements the multiallelic calling model. It has two major parts:
  *   1) determine the most likely set of alleles and init calculate the quality of ref/non-ref site
  *   2) determine and set the genotypes
  *  In various places in between, the BCF record gets updated.
  */
int mcall(call_t *call, bcf1_t *rec)
{
    int nals = rec->n_allele;
    if ( nals>4 )
        error("FIXME: Not ready for more than 4 alleles at %s:%d (%d)\n", call->hdr->id[BCF_DT_CTG][rec->rid].key,rec->pos+1, nals);

    // Get the genotype likelihoods
    int nsmpl = call->hdr->n[BCF_DT_SAMPLE];
    int npl = bcf_get_format_int(call->hdr, rec, "PL", &call->PLs, &call->nPLs);
    if ( npl!=nsmpl*nals*(nals+1)/2 && npl!=nsmpl*nals )
        error("Wrong number of PL fields? nals=%d npl=%d\n", nals,npl);

    // Convert PLs to probabilities
    int ngts = nals*(nals+1)/2;
    hts_expand(double, npl, call->npdg, call->pdg);
    set_pdg(call, call->PLs, call->pdg, call->hdr->n[BCF_DT_SAMPLE], ngts);

    // Get sum of qualities
    int i, nqs = bcf_get_info_float(call->hdr, rec, "QS", &call->qsum, &call->nqsum);
    assert( nals<=call->nqsum );
    for (i=nqs; i<nals; i++) call->qsum[i] = 0;

    // 1) Determine the most likely combination of alleles. In this implementation, at most three-allelic sites are considered.
    int ia,ib,ic;   // iterators over up to three alleles
    int max_als=0, max_als2=0;  // most likely and second-most likely combination of alleles
    double ref_lk = 0, max_lk = INT_MIN, max_lk2 = INT_MIN; // likelihood of the reference and of most likely combination of alleles
    double lk_sum = -HUGE_VAL;    // for normalizing the likelihoods

    // Single allele
    for (ia=0; ia<nals; ia++)
    {
        double lk_tot  = 0;
        int lk_tot_set = 0;
        int iaa = (ia+1)*(ia+2)/2-1;    // index in PL which corresponds to the homozygous "ia/ia" genotype
        int isample;
        double *pdg = call->pdg + iaa;
        for (isample=0; isample<nsmpl; isample++)
        {
            if ( *pdg ) { lk_tot += log(*pdg); lk_tot_set = 1; }
            pdg += ngts;
        }
        if ( ia==0 ) ref_lk = lk_tot;   // likelihood of 0/0 for all samples
        UPDATE_MAX_LKs(1<<ia);
    }

    // Two alleles
    if ( nals>1 )
    {
        for (ia=0; ia<nals; ia++)
        {
            if ( call->qsum[ia]==0 ) continue;
            int iaa = (ia+1)*(ia+2)/2-1;
            for (ib=0; ib<ia; ib++)
            {
                if ( call->qsum[ib]==0 ) continue;
                double lk_tot  = 0;
                int lk_tot_set = 0;
                double fa  = call->qsum[ia]/(call->qsum[ia]+call->qsum[ib]);
                double fb  = call->qsum[ib]/(call->qsum[ia]+call->qsum[ib]);
                double fab = 2*fa*fb; fa *= fa; fb *= fb;
                int isample, ibb = (ib+1)*(ib+2)/2-1, iab = iaa - ia + ib;
                double *pdg  = call->pdg;
                for (isample=0; isample<nsmpl; isample++)
                {
                    double val;
                    if ( call->ploidy && call->ploidy[isample]==1 )
                        val = fa*pdg[iaa] + fb*pdg[ibb];
                    else 
                        val = fa*pdg[iaa] + fb*pdg[ibb] + fab*pdg[iab];
                    if ( val ) { lk_tot += log(val); lk_tot_set = 1; }
                    pdg += ngts;
                }
                UPDATE_MAX_LKs(1<<ia|1<<ib);
            }
        }
    }

    // Three alleles
    if ( nals>2 )
    {
        for (ia=0; ia<nals; ia++)
        {
            if ( call->qsum[ia]==0 ) continue;
            int iaa = (ia+1)*(ia+2)/2-1;
            for (ib=0; ib<ia; ib++)
            {
                if ( call->qsum[ib]==0 ) continue;
                int ibb = (ib+1)*(ib+2)/2-1; 
                int iab = iaa - ia + ib;
                for (ic=0; ic<ib; ic++)
                {
                    if ( call->qsum[ic]==0 ) continue;
                    double lk_tot  = 0;
                    int lk_tot_set = 1;
                    double fa  = call->qsum[ia]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fb  = call->qsum[ib]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fc  = call->qsum[ic]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fab = 2*fa*fb, fac = 2*fa*fc, fbc = 2*fb*fc; fa *= fa; fb *= fb; fc *= fc;
                    int isample, icc = (ic+1)*(ic+2)/2-1;
                    int iac = iaa - ia + ic, ibc = ibb - ib + ic;
                    double *pdg = call->pdg;
                    for (isample=0; isample<nsmpl; isample++)
                    {
                        double val = 0;
                        if ( call->ploidy && call->ploidy[isample]==1 ) 
                            val = fa*pdg[iaa] + fb*pdg[ibb] + fc*pdg[icc];
                        else
                            val = fa*pdg[iaa] + fb*pdg[ibb] + fc*pdg[icc] + fab*pdg[iab] + fac*pdg[iac] + fbc*pdg[ibc];
                        if ( val ) { lk_tot += log(val); lk_tot_set = 1; }
                        pdg += ngts;
                    }
                    UPDATE_MAX_LKs(1<<ia|1<<ib|1<<ic);
                }
            }
        }
    }

    // How many alleles to call? Add a new allele only if it increases the likelihood significantly.
    // If the most likely set has more alleles than the second most likely set but the difference is
    // not big, go with fewer alleles. Here we use chi-squared distribution with 1 degree of freedom,
    // but is it correct? Is the number of degrees of freedom constant or does it depend on the
    // number of alleles? 
    int n1=0, n2=0;
    for (i=0; i<nals; i++) if ( max_als  & 1<<i) n1++;
    for (i=0; i<nals; i++) if ( max_als2 & 1<<i) n2++;

    //fprintf(stderr,"max_lk=%e max_lk2=%e ref_lk=%e lk_sum=%e lrt=%e (%e)  .. n1=%d n2=%d .. als1=%d als2=%d\n", max_lk,max_lk2,ref_lk,lk_sum,kf_gammap(0.5,max_lk-max_lk2),call->min_ma_lrt, n1,n2, max_als,max_als2);

    // Xi^2 = -2*ln(P0 / P1) and CDF(0.5*ndf,0.5*x)
    if ( n1>n2 && kf_gammap(0.5,max_lk-max_lk2)<call->min_ma_lrt )
    {
        // Going with fewer alleles, the likelihood is not significantly bigger
        SWAP(double, max_lk, max_lk2);
        SWAP(int, max_als, max_als2);
        SWAP(int, n1, n2);
    }

    // With -A, keep all ALTs except X
    if ( call->flag & CALL_KEEPALT )
    {
        n1 = 0;
        for (i=0; i<nals; i++)
            if ( rec->d.allele[i][0]!='X' ) { max_als |= 1<<i; n1++; }
    }
    // Make sure the REF allele is always present
    else if ( !(max_als&1) )
    {
        max_als |= 1;
        n1++;
    }


    // 2) Set the genotypes
    int ac[4] = {0,0,0,0};
    if ( max_als==1 )
    {
        // Reference only calls: return or set all genotypes to 0/0
        if ( call->flag & CALL_VARONLY ) return 0;
        init_allele_trimming_maps(call, max_als, nals);

        // Set the quality
        rec->qual = lk_sum==-HUGE_VAL ? 0 : -4.343*log(1 - exp(ref_lk - lk_sum));

        // Set all genotypes to 0/0 and remove PL vector
        int *gts    = call->gts;
        double *pdg = call->pdg;
        int isample;
        for (isample = 0; isample < nsmpl; isample++) 
        {
            int ploidy = call->ploidy ? call->ploidy[isample] : 2;
            for (i=0; i<ngts; i++) if ( pdg[i]!=0.0 ) break;
            if ( i==ngts ) 
            {
                gts[0] = bcf_gt_missing;
                gts[1] = ploidy==2 ? bcf_gt_missing : bcf_int32_vector_end;
            }
            else
            {
                gts[0] = bcf_gt_unphased(0);
                gts[1] = ploidy==2 ? bcf_gt_unphased(0) : bcf_int32_vector_end;
                ac[0] += ploidy;
            }
            gts += 2;
            pdg += ngts;
        }
        bcf1_update_format_int32(call->hdr, rec, "PL", NULL, 0);    // remove PL, useless now
    }
    else
    {
        // The most likely set of alleles includes non-reference allele (or was enforced), call genotypes.
        // Note that it is a valid outcome if the called genotypes exclude some of the ALTs.
        init_allele_trimming_maps(call, max_als, nals);

        int npls_src = ngts, npls_dst = n1*(n1+1)/2;     // number of PL values in diploid samples, ori and new
        double *pdg  = call->pdg - ngts;
        int *gts  = call->gts - 2;
        int *pls_src = call->PLs - npls_src, *pls_dst = call->PLs - npls_dst;

        int isample, ndiploid = 0, nhets = 0;
        for (isample = 0; isample < nsmpl; isample++) 
        {
            int ploidy = call->ploidy ? call->ploidy[isample] : 2;
            if ( ploidy==2 ) ndiploid++;

            pdg += ngts;
            gts += 2;
            pls_src += npls_src;
            pls_dst += npls_dst;

            // Skip samples with all pdg's equal to 1. These have zero depth.
            for (i=0; i<ngts; i++) if ( pdg[i]!=0.0 ) break;
            if ( i==ngts ) 
            {
                gts[0] = bcf_gt_missing;
                gts[1] = ploidy==2 ? bcf_gt_missing : bcf_int32_vector_end;
            }
            else
            {
                // Non-zero depth, determine the most likely genotype
                int ia;
                double best_lk = 0;
                for (ia=0; ia<nals; ia++)
                {
                    if ( !(max_als & 1<<ia) ) continue;     // ia-th allele not in the final selection, skip
                    int iaa = (ia+1)*(ia+2)/2-1;            // PL index of the ia/ia genotype
                    double lk = pdg[iaa]*call->qsum[ia]*call->qsum[ia];
                    if ( best_lk < lk ) 
                    { 
                        best_lk = lk; 
                        gts[0] = bcf_gt_unphased(call->als_map[ia]); 
                    }
                }
                if ( ploidy==2 ) 
                {
                    gts[1] = gts[0];
                    for (ia=0; ia<nals; ia++)
                    {
                        if ( !(max_als & 1<<ia) ) continue;
                        int iaa = (ia+1)*(ia+2)/2-1;
                        for (ib=0; ib<ia; ib++)
                        {
                            if ( !(max_als & 1<<ib) ) continue;
                            int iab = iaa - ia + ib;
                            double lk = 2*pdg[iab]*call->qsum[ia]*call->qsum[ib];
                            if ( best_lk < lk ) 
                            { 
                                best_lk = lk; 
                                gts[0] = bcf_gt_unphased(call->als_map[ia]); 
                                gts[1] = bcf_gt_unphased(call->als_map[ib]); 
                            }
                        }
                    }
                    if ( gts[0] != gts[1] ) nhets++;
                }
                else
                    gts[1] = bcf_int32_vector_end;
            }

            // Subset PL vector
            if ( npls_src!=npls_dst )
            {
                if ( ploidy==2 )
                {
                    for (ia=0; ia<npls_dst; ia++)
                        pls_dst[ia] =  pls_src[ call->pl_map[ia] ];
                }
                else
                {
                    for (ia=0; ia<n1; ia++)
                    {
                        int isrc = call->pl_map[ia]; 
                        isrc = (isrc+1)*(isrc+2)/2-1;
                        pls_dst[ia] = pls_src[isrc];
                    }
                    pls_dst[ia] = bcf_int32_vector_end;
                }
            }

            ac[ bcf_gt_allele(gts[0]) ]++;
            if ( gts[1]!=bcf_int32_vector_end ) ac[ bcf_gt_allele(gts[1]) ]++;
        }

        if ( npls_src!=npls_dst ) 
            bcf1_update_format_int32(call->hdr, rec, "PL", call->PLs, npls_dst*nsmpl);

        // Skip the site if all samples are 0/0. This can happen occasionally.
        if ( ac[1]+ac[2]+ac[3]==0 )
        {
            if ( call->flag & CALL_VARONLY ) return 0;
        }
        else
        {
            float icb = calc_ICB(ac[0],ac[1]+ac[2]+ac[3], nhets, ndiploid);
            bcf1_update_info_float(call->hdr, rec, "ICB", &icb, 1);
        }

        // Set the quality of a REF
        if ( lk_sum==-HUGE_VAL )
            rec->qual = 0;
        else 
            rec->qual = n1==1 ? -4.343*log(1 - exp(ref_lk - lk_sum)) : -4.343*(ref_lk - lk_sum);
    }

    hts_expand(char*,n1,call->nals,call->als);
    for (i=0; i<nals; i++)
        if ( call->als_map[i]>=0 ) call->als[call->als_map[i]] = rec->d.allele[i];  
    bcf1_update_alleles(call->hdr, rec, (const char**)call->als, n1);
    bcf1_update_genotypes(call->hdr, rec, call->gts, nsmpl*2);

    // DP4 tag
    int ndp = 16; float *anno16 = call->anno16;
    bcf_get_info_float(call->hdr, rec, "I16", &anno16, &ndp);
    assert( anno16==call->anno16 ); // this shouldn't happen unless the VCF is broken
    int32_t dp[4]; dp[0] = call->anno16[0]; dp[1] = call->anno16[1]; dp[2] = call->anno16[2]; dp[3] = call->anno16[3];
    bcf1_update_info_int32(call->hdr, rec, "DP4", dp, 4);

    bcf1_update_info_int32(call->hdr, rec, "I16", NULL, 0);     // remove I16 tag
    bcf1_update_info_int32(call->hdr, rec, "QS", NULL, 0);      // remove QS tag

    // AC, AN
    if ( n1>1 ) bcf1_update_info_int32(call->hdr, rec, "AC", ac+1, n1-1);
    ac[0] += ac[1] + ac[2] + ac[3];
    bcf1_update_info_int32(call->hdr, rec, "AN", ac, 1);

    if ( rec->qual>999 ) rec->qual = 999;
    if ( rec->qual>50 ) rec->qual = rint(rec->qual);

    return n1;
}

