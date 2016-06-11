/*
 This file is part of program wsprd, a detector/demodulator/decoder
 for the Weak Signal Propagation Reporter (WSPR) mode.

 File name: wsprd.c

 Copyright 2001-2015, Joe Taylor, K1JT

 Much of the present code is based on work by Steven Franke, K9AN,
 which in turn was based on earlier work by K1JT.

 Copyright 2014-2015, Steven Franke, K9AN

 License: GNU GPL v3

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <fftw3.h>

#include "wsprd.h"
#include "fano.h"
#include "nhash.h"
#include "wsprd_utils.h"
#include "wsprsim_utils.h"

#define DT       1.0/375.0
#define DF       375.0/256.0
#define TWOPIDT  2*M_PI*DT

#define max(x,y) ((x) > (y) ? (x) : (y))
// Possible PATIENCE options: FFTW_ESTIMATE, FFTW_ESTIMATE_PATIENT,
// FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
#define PATIENCE FFTW_ESTIMATE
fftwf_plan PLAN1,PLAN2,PLAN3;


unsigned char pr3[162]= {
    1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
    0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
    0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
    1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
    0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
    0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
    0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
    0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
    0,0
};

unsigned long nr;

int printdata=0;


//***************************************************************************
void sync_and_demodulate(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1, float fstep,
                         int *shift1, int lagmin, int lagmax, int lagstep,
                         float *drift1, int symfac, float *sync, int mode) {
    /***********************************************************************
     * mode = 0: no frequency or drift search. find best time lag.          *
     *        1: no time lag or drift search. find best frequency.          *
     *        2: no frequency or time lag search. calculate soft-decision   *
     *           symbols using passed frequency and shift.                  *
     ************************************************************************/

    float fbest=0.0;
    float f0=0.0,fp,ss;
    int lag;
    static float fplast=-10000.0;
    float i0[162],q0[162],i1[162],q1[162],i2[162],q2[162],i3[162],q3[162];
    float p0,p1,p2,p3,cmet,totp,syncmax,fac;
    float c0[256],s0[256],c1[256],s1[256],c2[256],s2[256],c3[256],s3[256];
    float dphi0, cdphi0, sdphi0,
          dphi1, cdphi1, sdphi1,
          dphi2, cdphi2, sdphi2,
          dphi3, cdphi3, sdphi3;
    float fsum=0.0, f2sum=0.0, fsymb[162];
    int best_shift = 0;
    int ifmin=0, ifmax=0;

    syncmax=-1e30;
    if( mode == 0 ) {
        ifmin=0;
        ifmax=0;
        fstep=0.0;
        f0=*f1;
    }
    if( mode == 1 ) {
        lagmin=*shift1;
        lagmax=*shift1;
        ifmin=-5;
        ifmax=5;
        f0=*f1;
    }
    if( mode == 2 ) {
        lagmin=*shift1;
        lagmax=*shift1;
        ifmin=0;
        ifmax=0;
        f0=*f1;
    }

    for(int ifreq=ifmin; ifreq<=ifmax; ifreq++) {
        f0=*f1+ifreq*fstep;
        for(lag=lagmin; lag<=lagmax; lag=lag+lagstep) {
            ss=0.0;
            totp=0.0;
            for (int i=0; i<162; i++) {
                fp = f0 + ((float)*drift1/2.0)*((float)i-81.0)/81.0;
                if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
                    dphi0=TWOPIDT*(fp-1.5*DF);
                    cdphi0=cosf(dphi0);
                    sdphi0=sinf(dphi0);

                    dphi1=TWOPIDT*(fp-0.5*DF);
                    cdphi1=cosf(dphi1);
                    sdphi1=sinf(dphi1);

                    dphi2=TWOPIDT*(fp+0.5*DF);
                    cdphi2=cosf(dphi2);
                    sdphi2=sinf(dphi2);

                    dphi3=TWOPIDT*(fp+1.5*DF);
                    cdphi3=cosf(dphi3);
                    sdphi3=sinf(dphi3);

                    c0[0]=1;
                    s0[0]=0;
                    c1[0]=1;
                    s1[0]=0;
                    c2[0]=1;
                    s2[0]=0;
                    c3[0]=1;
                    s3[0]=0;

                    for (int j=1; j<256; j++) {
                        c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                        s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                        c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                        s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                        c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                        s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                        c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                        s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
                    }
                    fplast = fp;
                }

                i0[i]=0.0;
                q0[i]=0.0;
                i1[i]=0.0;
                q1[i]=0.0;
                i2[i]=0.0;
                q2[i]=0.0;
                i3[i]=0.0;
                q3[i]=0.0;

                for (int j=0; j<256; j++) {
                    int k=lag+i*256+j;
                    if( (k>0) & (k<np) ) {
                        i0[i]=i0[i] + id[k]*c0[j] + qd[k]*s0[j];
                        q0[i]=q0[i] - id[k]*s0[j] + qd[k]*c0[j];
                        i1[i]=i1[i] + id[k]*c1[j] + qd[k]*s1[j];
                        q1[i]=q1[i] - id[k]*s1[j] + qd[k]*c1[j];
                        i2[i]=i2[i] + id[k]*c2[j] + qd[k]*s2[j];
                        q2[i]=q2[i] - id[k]*s2[j] + qd[k]*c2[j];
                        i3[i]=i3[i] + id[k]*c3[j] + qd[k]*s3[j];
                        q3[i]=q3[i] - id[k]*s3[j] + qd[k]*c3[j];
                    }
                }
                p0=i0[i]*i0[i] + q0[i]*q0[i];
                p1=i1[i]*i1[i] + q1[i]*q1[i];
                p2=i2[i]*i2[i] + q2[i]*q2[i];
                p3=i3[i]*i3[i] + q3[i]*q3[i];

                p0=sqrtf(p0);
                p1=sqrtf(p1);
                p2=sqrtf(p2);
                p3=sqrtf(p3);

                totp=totp+p0+p1+p2+p3;
                cmet=(p1+p3)-(p0+p2);
                ss=ss+cmet*(2*pr3[i]-1);
                if( mode == 2) {                 //Compute soft symbols
                    if(pr3[i]) {
                        fsymb[i]=p3-p1;
                    } else {
                        fsymb[i]=p2-p0;
                    }
                }
            }

            if( ss/totp > syncmax ) {          //Save best parameters
                syncmax=ss/totp;
                best_shift=lag;
                fbest=f0;
            }
        } // lag loop
    } //freq loop

    if( mode <=1 ) {                       //Send best params back to caller
        *sync=syncmax;
        *shift1=best_shift;
        *f1=fbest;
        return;
    }

    if( mode == 2 ) {
        *sync=syncmax;
        for (int i=0; i<162; i++) {              //Normalize the soft symbols
            fsum=fsum+fsymb[i]/162.0;
            f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
        }
        fac=sqrtf(f2sum-fsum*fsum);
        for (int i=0; i<162; i++) {
            fsymb[i]=symfac*fsymb[i]/fac;
            if( fsymb[i] > 127) fsymb[i]=127.0;
            if( fsymb[i] < -128 ) fsymb[i]=-128.0;
            symbols[i]=fsymb[i] + 128;
        }
        return;
    }
    return;
}


/***************************************************************************
 symbol-by-symbol signal subtraction
 ****************************************************************************/
void subtract_signal(float *id, float *qd, long np,
                     float f0, int shift0, float drift0, unsigned char* channel_symbols) {

    float i0,q0;
    float c0[256],s0[256];
    float dphi, cdphi, sdphi;

    for (int i=0; i<162; i++) {
        float fp = f0 + ((float)drift0/2.0)*((float)i-81.0)/81.0;

        dphi=TWOPIDT*(fp+((float)channel_symbols[i]-1.5)*DF);
        cdphi=cosf(dphi);
        sdphi=sinf(dphi);

        c0[0]=1;
        s0[0]=0;

        for (int j=1; j<256; j++) {
            c0[j]=c0[j-1]*cdphi - s0[j-1]*sdphi;
            s0[j]=c0[j-1]*sdphi + s0[j-1]*cdphi;
        }

        i0=0.0;
        q0=0.0;

        for (int j=0; j<256; j++) {
            int k=shift0+i*256+j;
            if( (k>0) & (k<np) ) {
                i0=i0 + id[k]*c0[j] + qd[k]*s0[j];
                q0=q0 - id[k]*s0[j] + qd[k]*c0[j];
            }
        }


        // subtract the signal here.

        i0=i0/256.0; //will be wrong for partial symbols at the edges...
        q0=q0/256.0;

        for (int j=0; j<256; j++) {
            int k=shift0+i*256+j;
            if( (k>0) & (k<np) ) {
                id[k]=id[k]- (i0*c0[j] - q0*s0[j]);
                qd[k]=qd[k]- (q0*c0[j] + i0*s0[j]);
            }
        }
    }
    return;
}


/******************************************************************************
 Fully coherent signal subtraction
 *******************************************************************************/
void subtract_signal2(float *id, float *qd, long np,
                      float f0, int shift0, float drift0, unsigned char* channel_symbols) {

    float phi=0, dphi, cs;
    int nsym=162, nspersym=256,  nfilt=256; //nfilt must be even number.
    int nsig=nsym*nspersym;
    int nc2=45000;

    float *refi, *refq, *ci, *cq, *cfi, *cfq;

    refi=malloc(sizeof(float)*nc2);
    refq=malloc(sizeof(float)*nc2);
    ci=malloc(sizeof(float)*nc2);
    cq=malloc(sizeof(float)*nc2);
    cfi=malloc(sizeof(float)*nc2);
    cfq=malloc(sizeof(float)*nc2);

    memset(refi,0,sizeof(float)*nc2);
    memset(refq,0,sizeof(float)*nc2);
    memset(ci,0,sizeof(float)*nc2);
    memset(cq,0,sizeof(float)*nc2);
    memset(cfi,0,sizeof(float)*nc2);
    memset(cfq,0,sizeof(float)*nc2);


    /******************************************************************************
     Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
     Reference is:                       r(t) = exp( j*phi(t) )
     Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
     so c(t) has phase angle theta-phi
     Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
     *******************************************************************************/

    // create reference wspr signal vector, centered on f0.
    //
    for (int i=0; i<nsym; i++) {

        cs=(float)channel_symbols[i];

        dphi=TWOPIDT* ( f0 +
                        ((float)drift0/2.0)*((float)i-(float)nsym/2.0)/((float)nsym/2.0) +
                        (cs-1.5)*DF  );

        for (int j=0; j<nspersym; j++ ) {
            int ii=nspersym*i+j;
            refi[ii]=refi[ii]+cosf(phi); //cannot precompute sin/cos because dphi is changing
            refq[ii]=refq[ii]+sinf(phi);
            phi=phi+dphi;
        }
    }

    // s(t) * conjugate(r(t))
    // beginning of first symbol in reference signal is at i=0
    // beginning of first symbol in received data is at shift0.
    // filter transient lasts nfilt samples
    // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    for (int i=0; i<nsym*nspersym; i++) {
        int k=shift0+i;
        if( (k>0) & (k<np) ) {
            ci[i+nfilt] = id[k]*refi[i] + qd[k]*refq[i];
            cq[i+nfilt] = qd[k]*refi[i] - id[k]*refq[i];
        }
    }

    //quick and dirty filter - may want to do better
    float w[nfilt], norm=0, partialsum[nfilt];
    memset(partialsum,0,sizeof(float)*nfilt);
    for (int i=0; i<nfilt; i++) {
        w[i]=sinf(M_PI*(float)i/(float)(nfilt-1));
        norm=norm+w[i];
    }
    for (int i=0; i<nfilt; i++) {
        w[i]=w[i]/norm;
    }
    for (int i=1; i<nfilt; i++) {
        partialsum[i]=partialsum[i-1]+w[i];
    }

    // LPF
    for (int i=nfilt/2; i<45000-nfilt/2; i++) {
        cfi[i]=0.0;
        cfq[i]=0.0;
        for (int j=0; j<nfilt; j++) {
            cfi[i]=cfi[i]+w[j]*ci[i-nfilt/2+j];
            cfq[i]=cfq[i]+w[j]*cq[i-nfilt/2+j];
        }
    }

    // subtract c(t)*r(t) here
    // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq)+cq*refi)
    // beginning of first symbol in reference signal is at i=nfilt
    // beginning of first symbol in received data is at shift0.
    for (int i=0; i<nsig; i++) {
        if( i<nfilt/2 ) {        // take care of the end effect (LPF step response) here
            norm=partialsum[nfilt/2+i];
        } else if( i>(nsig-1-nfilt/2) ) {
            norm=partialsum[nfilt/2+nsig-1-i];
        } else {
            norm=1.0;
        }
        int k=shift0+i;
        int j=i+nfilt;
        if( (k>0) & (k<np) ) {
            id[k]=id[k] - (cfi[j]*refi[i]-cfq[j]*refq[i])/norm;
            qd[k]=qd[k] - (cfi[j]*refq[i]+cfq[j]*refi[i])/norm;
        }
    }

    free(refi);
    free(refq);
    free(ci);
    free(cq);
    free(cfi);
    free(cfq);

    return;
}



//***************************************************************************
int wspr_decode(float *idat, float *qdat, unsigned int npoints,
                struct decoder_options options, struct decoder_results *decodes, int *n_results) {

    int i,j,k;
    unsigned char *symbols, *decdata;
    signed char message[]= {-9,13,-35,123,57,-39,64,0,0,0,0};
    char *callsign, *call_loc_pow, *call, *loc, *pwr;
    char *data_dir=NULL;
    char wisdom_fname[200],all_fname[200],spots_fname[200];
    char timer_fname[200],hash_fname[200];
    int delta,verbose=0;
    int writenoise=0,wspr_type=2, ipass;
    int shift1, lagmin, lagmax, lagstep, worth_a_try, not_decoded;
    unsigned int nbits=81;
    unsigned int metric, maxcycles, cycles, maxnp;
    float freq0[200],snr0[200],drift0[200],sync0[200];
    int shift0[200];
    float dt_print;
    float freq_print;
    float dialfreq= (float)options.freq / 1e6; // check
    float dialfreq_error=0.0;
    float f1, fstep, sync1=0.0, drift1;
    int noprint=0;
    int uniques=0;
    float fmin=-110.0;
    float fmax=110.0;

    char *hashtab;
    hashtab=malloc(sizeof(char)*32768*13);
    memset(hashtab,0,sizeof(char)*32768*13);

    int nh;
    symbols=malloc(sizeof(char)*nbits*2);
    decdata=malloc((nbits+7)/8);
    callsign=malloc(sizeof(char)*13);
    call_loc_pow=malloc(sizeof(char)*23);
    call=malloc(sizeof(char)*13);
    loc=malloc(sizeof(char)*7);
    pwr=malloc(sizeof(char)*3);


    float allfreqs[100];
    char allcalls[100][13];
    memset(allfreqs,0,sizeof(float)*100);
    memset(allcalls,0,sizeof(char)*100*13);


    // Parameters used for performance-tuning:
    maxcycles=10000;                         //Fano timeout limit
    double minsync1=0.10;                    //First sync limit
    double minsync2=0.12;                    //Second sync limit
    int iifac=3;                             //Step size in final DT peakup
    int symfac=50;                           //Soft-symbol normalizing factor
    int maxdrift=4;                          //Maximum (+/-) drift
    double minrms=52.0 * (symfac/64.0);      //Final test for plausible decoding
    delta=60;                                //Fano threshold step

    fftwf_complex *fftin, *fftout;

#include "./metric_tables.c"

    int mettab[2][256];
    float bias=0.42;

    // setup metric table
    for(i=0; i<256; i++) {
        mettab[0][i]=round( 10*(metric_tables[2][i]-bias) );
        mettab[1][i]=round( 10*(metric_tables[2][255-i]-bias) );
    }

    FILE *fp_fftw_wisdom_file, *fhash;
    strcpy(wisdom_fname,".");
    strcpy(all_fname,".");
    strcpy(spots_fname,".");
    strcpy(timer_fname,".");
    strcpy(hash_fname,".");
    if(data_dir != NULL) {
        strcpy(wisdom_fname,data_dir);
        strcpy(all_fname,data_dir);
        strcpy(spots_fname,data_dir);
        strcpy(timer_fname,data_dir);
        strcpy(hash_fname,data_dir);
    }
    strncat(wisdom_fname,"/wspr_wisdom.dat",20);
    strncat(all_fname,"/ALL_WSPR.TXT",20);
    strncat(spots_fname,"/wspr_spots.txt",20);
    strncat(timer_fname,"/wspr_timer.out",20);
    strncat(hash_fname,"/hashtable.txt",20);
    if ((fp_fftw_wisdom_file = fopen(wisdom_fname, "r"))) {  //Open FFTW wisdom
        fftwf_import_wisdom_from_file(fp_fftw_wisdom_file);
        fclose(fp_fftw_wisdom_file);
    }

    // Do windowed ffts over 2 symbols, stepped by half symbols
    int nffts=4*floor(npoints/512)-1;
    fftin=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
    fftout=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
    PLAN3 = fftwf_plan_dft_1d(512, fftin, fftout, FFTW_FORWARD, PATIENCE);

    float ps[512][nffts];
    float w[512];
    for(i=0; i<512; i++) {
        w[i]=sin(0.006147931*i);
    }

    if( options.usehashtable ) {
        char line[80], hcall[12];
        if( (fhash=fopen(hash_fname,"r+")) ) {
            while (fgets(line, sizeof(line), fhash) != NULL) {
                sscanf(line,"%d %s",&nh,hcall);
                strcpy(hashtab+nh*13,hcall);
            }
        } else {
            fhash=fopen(hash_fname,"w+");
        }
        fclose(fhash);
    }

    //*************** main loop starts here *****************
    for (ipass=0; ipass<options.npasses; ipass++) {

        if( ipass == 1 && uniques == 0 ) break;
        if( ipass == 1 ) {  //otherwise we bog down on the second pass
            options.quickmode = 1;
        }

        memset(ps,0.0, sizeof(float)*512*nffts);
        for (i=0; i<nffts; i++) {
            for(j=0; j<512; j++ ) {
                k=i*128+j;
                fftin[j][0]=idat[k] * w[j];
                fftin[j][1]=qdat[k] * w[j];
            }
            fftwf_execute(PLAN3);
            for (j=0; j<512; j++ ) {
                k=j+256;
                if( k>511 )
                    k=k-512;
                ps[j][i]=fftout[k][0]*fftout[k][0]+fftout[k][1]*fftout[k][1];
            }
        }

        // Compute average spectrum
        float psavg[512]= {0};
        //memset(psavg,0.0, sizeof(float)*512);
        for (i=0; i<nffts; i++) {
            for (j=0; j<512; j++) {
                psavg[j]=psavg[j]+ps[j][i];
            }
        }

        // Smooth with 7-point window and limit spectrum to +/-150 Hz
        int window[7]= {1,1,1,1,1,1,1};
        float smspec[411];
        for (i=0; i<411; i++) {
            smspec[i]=0.0;
            for(j=-3; j<=3; j++) {
                k=256-205+i+j;
                smspec[i]=smspec[i]+window[j+3]*psavg[k];
            }
        }

        // Sort spectrum values, then pick off noise level as a percentile
        float tmpsort[411];
        for (j=0; j<411; j++) {
            tmpsort[j]=smspec[j];
        }
        qsort(tmpsort, 411, sizeof(float), floatcomp);

        // Noise level of spectrum is estimated as 123/411= 30'th percentile
        float noise_level = tmpsort[122];

        /* Renormalize spectrum so that (large) peaks represent an estimate of snr.
         * We know from experience that threshold snr is near -7dB in wspr bandwidth,
         * corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
         * The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15. */

        float min_snr, snr_scaling_factor;
        min_snr = pow(10.0,-7.0/10.0); //this is min snr in wspr bw
        if( wspr_type == 2 ) {
            snr_scaling_factor=26.3;
        } else {
            snr_scaling_factor=35.3;
        }
        for (j=0; j<411; j++) {
            smspec[j]=smspec[j]/noise_level - 1.0;
            if( smspec[j] < min_snr) smspec[j]=0.1;
            continue;
        }

        // Find all local maxima in smoothed spectrum.
        for (i=0; i<200; i++) {
            freq0[i]=0.0;
            snr0[i]=0.0;
            drift0[i]=0.0;
            shift0[i]=0;
            sync0[i]=0.0;
        }

        int npk=0;
        for(j=1; j<410; j++) {
            if((smspec[j]>smspec[j-1]) && (smspec[j]>smspec[j+1]) && (npk<200)) {
                freq0[npk]=(j-205)*(DF/2.0);
                snr0[npk]=10*log10(smspec[j])-snr_scaling_factor;
                npk++;
            }
        }

        // Compute corrected fmin, fmax, accounting for dial frequency error
        fmin += dialfreq_error;    // dialfreq_error is in units of Hz
        fmax += dialfreq_error;

        // Don't waste time on signals outside of the range [fmin,fmax].
        i=0;
        for( j=0; j<npk; j++) {
            if( freq0[j] >= fmin && freq0[j] <= fmax ) {
                freq0[i]=freq0[j];
                snr0[i]=snr0[j];
                i++;
            }
        }
        npk=i;

        // bubble sort on snr, bringing freq along for the ride
        int pass;
        float tmp;
        for (pass = 1; pass <= npk - 1; pass++) {
            for (k = 0; k < npk - pass ; k++) {
                if (snr0[k] < snr0[k+1]) {
                    tmp = snr0[k];
                    snr0[k] = snr0[k+1];
                    snr0[k+1] = tmp;
                    tmp = freq0[k];
                    freq0[k] = freq0[k+1];
                    freq0[k+1] = tmp;
                }
            }
        }

        /* Make coarse estimates of shift (DT), freq, and drift
         * Look for time offsets up to +/- 8 symbols (about +/- 5.4 s) relative
           to nominal start time, which is 2 seconds into the file
         * Calculates shift relative to the beginning of the file
         * Negative shifts mean that signal started before start of file
         * The program prints DT = shift-2 s
         * Shifts that cause sync vector to fall off of either end of the data
           vector are accommodated by "partial decoding", such that missing
           symbols produce a soft-decision symbol value of 128
         * The frequency drift model is linear, deviation of +/- drift/2 over the
           span of 162 symbols, with deviation equal to 0 at the center of the
           signal vector.
         */

        int idrift,ifr,if0,ifd,k0;
        int kindex;
        float smax,ss,pow,p0,p1,p2,p3;
        for(j=0; j<npk; j++) {                              //For each candidate...
            smax=-1e30;
            if0=freq0[j]/(DF/2.0)+256;
            for (ifr=if0-1; ifr<=if0+1; ifr++) {                      //Freq search
                for( k0=-10; k0<22; k0++) {                             //Time search
                    for (idrift=-maxdrift; idrift<=maxdrift; idrift++) {  //Drift search
                        ss=0.0;
                        pow=0.0;
                        for (k=0; k<162; k++) {                             //Sum over symbols
                            ifd=ifr+((float)k-81.0)/81.0*( (float)idrift )/DF;
                            kindex=k0+2*k;
                            if( kindex < nffts ) {
                                p0=ps[ifd-3][kindex];
                                p1=ps[ifd-1][kindex];
                                p2=ps[ifd+1][kindex];
                                p3=ps[ifd+3][kindex];

                                p0=sqrt(p0);
                                p1=sqrt(p1);
                                p2=sqrt(p2);
                                p3=sqrt(p3);

                                ss=ss+(2*pr3[k]-1)*((p1+p3)-(p0+p2));
                                pow=pow+p0+p1+p2+p3;
                                sync1=ss/pow;
                            }
                        }
                        if( sync1 > smax ) {                  //Save coarse parameters
                            smax=sync1;
                            shift0[j]=128*(k0+1);
                            drift0[j]=idrift;
                            freq0[j]=(ifr-256)*(DF/2.0);
                            sync0[j]=sync1;
                        }
                    }
                }
            }
        }

        /*
         Refine the estimates of freq, shift using sync as a metric.
         Sync is calculated such that it is a float taking values in the range
         [0.0,1.0].

         Function sync_and_demodulate has three modes of operation
         mode is the last argument:

         0 = no frequency or drift search. find best time lag.
         1 = no time lag or drift search. find best frequency.
         2 = no frequency or time lag search. Calculate soft-decision
         symbols using passed frequency and shift.

         NB: best possibility for OpenMP may be here: several worker threads
         could each work on one candidate at a time.
         */

        for (j=0; j<npk; j++) {
            memset(symbols,0,sizeof(char)*nbits*2);
            memset(callsign,0,sizeof(char)*13);
            memset(call_loc_pow,0,sizeof(char)*23);
            memset(call,0,sizeof(char)*13);
            memset(loc,0,sizeof(char)*7);
            memset(pwr,0,sizeof(char)*3);

            f1=freq0[j];
            drift1=drift0[j];
            shift1=shift0[j];
            sync1=sync0[j];

            // Fine search for best sync lag (mode 0)
            fstep=0.0;
            lagmin=shift1-144;
            lagmax=shift1+144;
            lagstep=8;
            if(options.quickmode) lagstep=16;

            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);

            // Fine search for frequency peak (mode 1)
            fstep=0.1;
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);

            if( sync1 > minsync1 ) {
                worth_a_try = 1;
            } else {
                worth_a_try = 0;
            }

            int idt=0, ii=0, jiggered_shift;
            double y,sq,rms;
            not_decoded=1;

            while ( worth_a_try && not_decoded && idt<=(128/iifac)) {
                ii=(idt+1)/2;
                if( idt%2 == 1 ) ii=-ii;
                ii=iifac*ii;
                jiggered_shift=shift1+ii;

                // Use mode 2 to get soft-decision symbols
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, fstep,
                                    &jiggered_shift, lagmin, lagmax, lagstep, &drift1, symfac,
                                    &sync1, 2);

                sq=0.0;
                for(i=0; i<162; i++) {
                    y=(double)symbols[i] - 128.0;
                    sq += y*y;
                }
                rms=sqrt(sq/162.0);

                if((sync1 > minsync2) && (rms > minrms)) {
                    deinterleave(symbols);

                    not_decoded = fano(&metric,&cycles,&maxnp,decdata,symbols,nbits,
                                       mettab,delta,maxcycles);
                }
                idt++;
                if( options.quickmode ) break;
            }

            if( worth_a_try && !not_decoded ) {
                for(i=0; i<11; i++) {
                    if( decdata[i]>127 ) {
                        message[i]=decdata[i]-256;
                    } else {
                        message[i]=decdata[i];
                    }
                }

                // Unpack the decoded message, update the hashtable, apply
                // sanity checks on grid and power, and return
                // call_loc_pow string and also callsign (for de-duping).
                noprint=unpk_(message,hashtab,call_loc_pow,call,loc,pwr,callsign);

                if( options.subtraction && (ipass == 0) && !noprint ) {
                    unsigned char channel_symbols[162];

                    if( get_wspr_channel_symbols(call_loc_pow, hashtab, channel_symbols) ) {
                        subtract_signal2(idat, qdat, npoints, f1, shift1, drift1, channel_symbols);
                    } else {
                        break;
                    }

                }

                // Remove dupes (same callsign and freq within 3 Hz)
                int dupe=0;
                for (i=0; i<uniques; i++) {
                    if(!strcmp(callsign,allcalls[i]) && (fabs(f1-allfreqs[i]) <3.0))
                        dupe=1;
                }

                if( (verbose || !dupe) && !noprint) {
                    strcpy(allcalls[uniques],callsign);
                    allfreqs[uniques]=f1;
                    uniques++;

                    // Add an extra space at the end of each line so that wspr-x doesn't
                    // truncate the power (TNX to DL8FCL!)
                    if( wspr_type == 15 ) {
                        freq_print=dialfreq+(1500+112.5+f1/8.0)/1e6;
                        dt_print=shift1*8*DT-2.0;
                    } else {
                        freq_print=dialfreq+(1500+f1)/1e6;
                        dt_print=shift1*DT-2.0;
                    }

                    decodes[uniques-1].sync=sync1;
                    decodes[uniques-1].snr=snr0[j];
                    decodes[uniques-1].dt=dt_print;
                    decodes[uniques-1].freq=freq_print;
                    decodes[uniques-1].drift=drift1;
                    decodes[uniques-1].cycles=cycles;
                    decodes[uniques-1].jitter=ii;
                    strcpy(decodes[uniques-1].message,call_loc_pow);
                    strcpy(decodes[uniques-1].call,call);
                    strcpy(decodes[uniques-1].loc,loc);
                    strcpy(decodes[uniques-1].pwr,pwr);
                }
            }
        }
    }

    // sort the result
    struct decoder_results temp;
    for (j = 1; j <= uniques - 1; j++) {
        for (k = 0; k < uniques - j ; k++) {
            if (decodes[k].snr < decodes[k+1].snr) {
                temp = decodes[k];
                decodes[k]=decodes[k+1];;
                decodes[k+1] = temp;
            }
        }
    }

    // Return number of spots to the calling fct
    *n_results = uniques;

    fftwf_free(fftin);
    fftwf_free(fftout);

    if ((fp_fftw_wisdom_file = fopen(wisdom_fname, "w"))) {
        fftwf_export_wisdom_to_file(fp_fftw_wisdom_file);
        fclose(fp_fftw_wisdom_file);
    }

    fftwf_destroy_plan(PLAN1);
    fftwf_destroy_plan(PLAN2);
    fftwf_destroy_plan(PLAN3);

    if( options.usehashtable ) {
        fhash=fopen(hash_fname,"w");
        for (i=0; i<32768; i++) {
            if( strncmp(hashtab+i*13,"\0",1) != 0 ) {
                fprintf(fhash,"%5d %s\n",i,hashtab+i*13);
            }
        }
        fclose(fhash);
    }

    if(writenoise == 999) return -1;  //Silence compiler warning
    return 0;
}
