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
    call->gts = (int*) calloc(call->hdr_out->n[BCF_DT_SAMPLE]*2,sizeof(int));   // assuming at most diploid everywhere

    bcf_hdr_append(call->hdr_out,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

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


// Inits P(D|G): convert PLs from log space and normalize
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
        assert( sum!=0 );
        for (j=0; j<n_gt; j++) pdg[j] /= sum;

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


// Macro to set the most likely and second most likely alleles
#define UPDATE_MAX_LKs(als) { \
     if ( max_lk<lk_tot ) { max_lk2 = max_lk; max_als2 = max_als; max_lk = lk_tot; max_als = (als); } \
     else if ( max_lk2<lk_tot ) { max_lk2 = lk_tot; max_als2 = (als); } \
     lk_sum = lk_tot>lk_sum ? lk_tot + log(1+exp(lk_sum-lk_tot)) : lk_sum + log(1+exp(lk_tot-lk_sum)); \
}

#define SWAP(type_t,x,y) {type_t tmp; tmp = x; x = y; y = tmp; }

int mcall(call_t *call, bcf1_t *rec)
{
    int nals = rec->n_allele;
    if ( nals>4 )
        error("FIXME: Not ready for more than 4 alleles at %s:%d (%d)\n", call->hdr_in->id[BCF_DT_CTG][rec->rid].key,rec->pos+1, nals);

    // Get the genotype likelihoods
    int nsmpl = call->hdr_in->n[BCF_DT_SAMPLE];
    int npl = bcf_get_format_int(call->hdr_in, rec, "PL", &call->PLs, &call->nPLs);
    if ( npl!=nsmpl*nals*(nals+1)/2 && npl!=nsmpl*nals )
        error("Wrong number of PL fields? nals=%d npl=%d\n", nals,npl);

    // Convert PLs to probabilities
    int ngts = nals*(nals+1)/2;
    hts_expand(double, npl, call->npdg, call->pdg);
    set_pdg(call, call->PLs, call->pdg, call->hdr_in->n[BCF_DT_SAMPLE], ngts);

    // Get sum of qualities
    int i, nqs = bcf_get_info_float(call->hdr_in, rec, "QS", &call->qsum, &call->nqsum);
    assert( nals<call->nqsum );
    for (i=nqs; i<nals; i++) call->qsum[i] = 0;

    // Determine the most likely combination of alleles. In this implementation, at most three-allelic sites are considered.
    int ia,ib,ic;   // iterators over up to three alleles
    int max_als=0, max_als2=0;  // most likely and second-most likely combination of alleles
    double ref_lk = 0, max_lk = INT_MIN, max_lk2 = INT_MIN; // likelihood of the reference and of most likely combination of alleles
    double lk_sum = INT_MIN, lk_sums[3];    // for normalizing the likelihoods

    // Single allele
    for (ia=0; ia<nals; ia++)
    {
        double lk_tot = 0;
        int iaa = (ia+1)*(ia+2)/2-1;    // index in PL which corresponds to the homozygous "ia/ia" genotype
        int isample;
        double *pdg = call->pdg + iaa;
        for (isample=0; isample<nsmpl; isample++)
        {
            lk_tot += log(*pdg);
            pdg += ngts;
        }
        if ( ia==0 ) ref_lk = lk_tot;   // likelihood of 0/0 for all samples
        UPDATE_MAX_LKs(1<<ia);
    }
    lk_sums[0] = lk_sum;    // normalizing term for single allele genotypes

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
                double lk_tot = 0;
                double fa  = call->qsum[ia]/(call->qsum[ia]+call->qsum[ib]);
                double fb  = call->qsum[ib]/(call->qsum[ia]+call->qsum[ib]);
                double fab = 2*fa*fb; fa *= fa; fb *= fb;
                int isample, ibb = (ib+1)*(ib+2)/2-1, iab = iaa - ia + ib;
                double *pdg  = call->pdg;
                for (isample=0; isample<nsmpl; isample++)
                {
                    if ( call->ploidy && call->ploidy[isample]==1 )
                        lk_tot +=  log(fa*pdg[iaa] + fb*pdg[ibb]);
                    else 
                        lk_tot +=  log(fa*pdg[iaa] + fb*pdg[ibb] + fab*pdg[iab]);
                    pdg += ngts;
                }
                UPDATE_MAX_LKs(1<<ia|1<<ib);
            }
        }
        lk_sums[1] = lk_sum;
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
                    double lk_tot = 0;
                    double fa  = call->qsum[ia]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fb  = call->qsum[ib]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fc  = call->qsum[ic]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fab = 2*fa*fb, fac = 2*fa*fc, fbc = 2*fb*fc; fa *= fa; fb *= fb; fc *= fc;
                    int isample, icc = (ic+1)*(ic+2)/2-1;
                    int iac = iaa - ia + ic, ibc = ibb - ib + ic;
                    double *pdg = call->pdg;
                    for (isample=0; isample<nsmpl; isample++)
                    {
                        if ( call->ploidy && call->ploidy[isample]==1 ) 
                            lk_tot += log(fa*pdg[iaa] + fb*pdg[ibb] + fc*pdg[icc]);
                        else
                            lk_tot += log(fa*pdg[iaa] + fb*pdg[ibb] + fc*pdg[icc] + fab*pdg[iab] + fac*pdg[iac] + fbc*pdg[ibc]);
                        pdg += ngts;
                    }
                    UPDATE_MAX_LKs(1<<ia|1<<ib|1<<ic);
                }
            }
        }
        lk_sums[2] = lk_sum;
    }

    // How many alleles to call? Add a new allele only if it increases the likelihood significantly.
    // If the most likely set has more alleles than the second most likely set but the difference is
    // not big, go with fewer alleles. Here we use chi-squared distribution with 1 degree of freedom,
    // but is it correct? Is the number of degrees of freedom constant or does it depend on the
    // number of alleles? 
    int n1=0, n2=0;
    for (i=0; i<nals; i++) if ( max_als  & 1<<i) n1++;
    for (i=0; i<nals; i++) if ( max_als2 & 1<<i) n2++;

//printf("%e %e (%e %e) .. %d %d .. %d %d\n", max_lk,max_lk2, max_lk - lk_sums[n1-1],max_lk2 - lk_sums[n2-1], n1,n2, max_als,max_als2);

    // Xi^2 = -2*ln(P0 / P1) and CDF(0.5*ndf,0.5*x)
    if ( n1>n2 && kf_gammap(0.5,max_lk-max_lk2)<call->min_ma_lrt )
    {
        // Going with fewer alleles, the likelihood is not significantly bigger
        SWAP(double, max_lk, max_lk2);
        SWAP(int, max_als, max_als2);
        SWAP(int, n1, n2);
    }

    // Call genotypes
    if ( max_als==1 )
    {
        if ( call->flag & CALL_VARONLY ) return 0;
        assert(0); // todo: set all gts to 0
    }
    else
    {
        init_allele_trimming_maps(call, max_als, nals);

        int npls_src = ngts, npls_dst = n1*(n1+1)/2;     // number of PL values in diploid samples, ori and new
        double *pdg  = call->pdg - ngts;
        int *gts  = call->gts - 2;
        int *pls_src = call->PLs - npls_src, *pls_dst = call->PLs - npls_dst;

        int isample, ac[4] = {0,0,0,0};
        int ndiploid = 0, nhets = 0;
        for (isample = 0; isample < nsmpl; isample++) 
        {
            int ploidy = call->ploidy ? call->ploidy[isample] : 2;
            if ( ploidy==2 ) ndiploid++;

            pdg += ngts;
            gts += 2;
            pls_src += npls_src;
            pls_dst += npls_dst;

            // Determine the most likely genotype
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

        // Skip the site if all samples are 0/0. It can happen occasionally.
        if ( !nhets && (call->flag & CALL_VARONLY) && ac[1]+ac[2]+ac[3]==0 ) return 0;
    }

    hts_expand(char*,n1,call->nals,call->als);
    for (i=0; i<nals; i++)
        if ( call->als_map[i]>=0 ) call->als[call->als_map[i]] = rec->d.allele[i];
    bcf1_update_alleles(call->hdr_out, rec, (const char**)call->als, n1);

    bcf1_update_genotypes(call->hdr_out, rec, call->gts, nsmpl*2);
    bcf1_update_info_int32(call->hdr_out, rec, "I16", NULL, 0);     // remove I16 tag
    bcf1_update_info_int32(call->hdr_out, rec, "QS", NULL, 0);      // remove QS tag

// ICF
    rec->qual = nals>1 ? -4.343*(ref_lk - max_lk) : -4.343*log(1-exp(ref_lk - max_lk));
    if ( rec->qual>999 ) rec->qual = 999;
    if ( rec->qual>50 ) rec->qual = rint(rec->qual);

    return n1;
}

