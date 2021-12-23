/*
 This file is part of program wsprd, a detector/demodulator/decoder
 for the Weak Signal Propagation Reporter (WSPR) mode.

 File name: wsprd.c

 Copyright 2001-2015, Joe Taylor, K1JT

 Much of the present code is based on work by Steven Franke, K9AN,
 which in turn was based on earlier work by K1JT.

 Copyright 2014-2015, Steven Franke, K9AN

 Minor modifications

 Copyright 2016, Guenael Jouchet, VA2GKA

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

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>

#include "./wsprd.h"
#include "./fano.h"
#include "./nhash.h"
#include "./wsprd_utils.h"
#include "./wsprsim_utils.h"
#include "./metric_tables.h"

#define SIGNAL_LENGHT       120
#define SIGNAL_SAMPLE_RATE  375
#define SIGNAL_SAMPLES      SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE
#define NBITS               81
#define NSYM                162
#define NSPERSYM            256
#define DF                  375.0 / 256.0
#define DT                  1.0 / 375.0
#define DF05                DF * 0.5
#define DF15                DF * 1.5
#define TWOPIDT             2.0 * M_PI * DT


/* Possible PATIENCE options: F
   FTW_ESTIMATE,
   FFTW_ESTIMATE_PATIENT,
   FFTW_MEASURE,
   FFTW_PATIENT,
   FFTW_EXHAUSTIVE
*/
#define PATIENCE FFTW_ESTIMATE

fftwf_plan PLAN;
int32_t printdata = 0;

uint8_t pr3vector[NSYM] = {
    1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0,
    0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1,
    1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1,
    0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0,
    0, 0};


/* mode = 0: no frequency or drift search. find best time lag.
 *        1: no time lag or drift search. find best frequency.
 *        2: no frequency or time lag search. calculate soft-decision
 *           symbols using passed frequency and shift.
 */
void sync_and_demodulate(float *id,
                         float *qd,
                         long  np,
                         unsigned char *symbols,
                         float *freq,
                         int   ifmin,
                         int   ifmax,
                         float fstep,
                         int   *shift,
                         int   lagmin,
                         int   lagmax,
                         int   lagstep,
                         float *drift,
                         int   symfac,
                         float *sync,
                         int   mode) {

    float i0[NSYM], q0[NSYM], 
          i1[NSYM], q1[NSYM],
          i2[NSYM], q2[NSYM],
          i3[NSYM], q3[NSYM];
    float c0[NSPERSYM], s0[NSPERSYM],
          c1[NSPERSYM], s1[NSPERSYM],
          c2[NSPERSYM], s2[NSPERSYM],
          c3[NSPERSYM], s3[NSPERSYM];
    float fsymb[NSYM];

    float fbest = 0.0,
          fsum  = 0.0,
          f2sum = 0.0;

    int   best_shift = 0;
    static float fplast = -10000.0;
    float syncmax = -1e30;

    if (mode == 0) {
        ifmin = 0;
        ifmax = 0;
        fstep = 0.0;
    } else if (mode == 1) {
        lagmin = *shift;
        lagmax = *shift;
    } else if (mode == 2) {
        lagmin = *shift;
        lagmax = *shift;
        ifmin  = 0;
        ifmax  = 0;
    }

    for (int ifreq = ifmin; ifreq <= ifmax; ifreq++) {
        float f0 = *freq + ifreq * fstep;
        for (int lag = lagmin; lag <= lagmax; lag = lag + lagstep) {
            float ss = 0.0;
            float totp = 0.0;
            for (int i = 0; i < NSYM; i++) {
                float fp = f0 + (*drift / 2.0) * ((float)i - (float)NBITS) / (float)NBITS;
                if (i == 0 || (fp != fplast)) {  // only calculate sin/cos if necessary
                    float dphi0  = TWOPIDT * (fp - DF15);
                    float cdphi0 = cosf(dphi0);
                    float sdphi0 = sinf(dphi0);

                    float dphi1  = TWOPIDT * (fp - DF05);
                    float cdphi1 = cosf(dphi1);
                    float sdphi1 = sinf(dphi1);

                    float dphi2  = TWOPIDT * (fp + DF05);
                    float cdphi2 = cosf(dphi2);
                    float sdphi2 = sinf(dphi2);

                    float dphi3  = TWOPIDT * (fp + DF15);
                    float cdphi3 = cosf(dphi3);
                    float sdphi3 = sinf(dphi3);

                    c0[0] = 1; s0[0] = 0;
                    c1[0] = 1; s1[0] = 0;
                    c2[0] = 1; s2[0] = 0;
                    c3[0] = 1; s3[0] = 0;

                    for (int j = 1; j < NSPERSYM; j++) {
                        c0[j] = c0[j - 1] * cdphi0 - s0[j - 1] * sdphi0;
                        s0[j] = c0[j - 1] * sdphi0 + s0[j - 1] * cdphi0;
                        c1[j] = c1[j - 1] * cdphi1 - s1[j - 1] * sdphi1;
                        s1[j] = c1[j - 1] * sdphi1 + s1[j - 1] * cdphi1;
                        c2[j] = c2[j - 1] * cdphi2 - s2[j - 1] * sdphi2;
                        s2[j] = c2[j - 1] * sdphi2 + s2[j - 1] * cdphi2;
                        c3[j] = c3[j - 1] * cdphi3 - s3[j - 1] * sdphi3;
                        s3[j] = c3[j - 1] * sdphi3 + s3[j - 1] * cdphi3;
                    }
                    fplast = fp;
                }

                i0[i] = 0.0; q0[i] = 0.0;
                i1[i] = 0.0; q1[i] = 0.0;
                i2[i] = 0.0; q2[i] = 0.0;
                i3[i] = 0.0; q3[i] = 0.0;

                for (int j = 0; j < NSPERSYM; j++) {
                    int k = lag + i * NSPERSYM + j;
                    if ((k > 0) && (k < np)) {
                        i0[i] = i0[i] + id[k] * c0[j] + qd[k] * s0[j];
                        q0[i] = q0[i] - id[k] * s0[j] + qd[k] * c0[j];
                        i1[i] = i1[i] + id[k] * c1[j] + qd[k] * s1[j];
                        q1[i] = q1[i] - id[k] * s1[j] + qd[k] * c1[j];
                        i2[i] = i2[i] + id[k] * c2[j] + qd[k] * s2[j];
                        q2[i] = q2[i] - id[k] * s2[j] + qd[k] * c2[j];
                        i3[i] = i3[i] + id[k] * c3[j] + qd[k] * s3[j];
                        q3[i] = q3[i] - id[k] * s3[j] + qd[k] * c3[j];
                    }
                }

                float p0 = sqrt(i0[i] * i0[i] + q0[i] * q0[i]);
                float p1 = sqrt(i1[i] * i1[i] + q1[i] * q1[i]);
                float p2 = sqrt(i2[i] * i2[i] + q2[i] * q2[i]);
                float p3 = sqrt(i3[i] * i3[i] + q3[i] * q3[i]);

                totp = totp + p0 + p1 + p2 + p3;
                float cmet = (p1 + p3) - (p0 + p2);
                ss = (pr3vector[i] == 1) ? ss + cmet : ss - cmet;
                if (mode == 2) {  // Compute soft symbols
                    if (pr3vector[i] == 1) {
                        fsymb[i] = p3 - p1;
                    } else {
                        fsymb[i] = p2 - p0;
                    }
                }
            }
            ss = ss / totp;
            if (ss > syncmax) {  // Save best parameters
                syncmax = ss;
                best_shift = lag;
                fbest = f0;
            }
        }  // lag loop
    }  // freq loop

    if (mode <= 1) {  // Send best params back to caller
        *sync = syncmax;
        *shift = best_shift;
        *freq = fbest;
        return;
    }

    if (mode == 2) {
        *sync = syncmax;
        for (int i = 0; i < NSYM; i++) {  // Normalize the soft symbols
            fsum  += fsymb[i] / NSYM;
            f2sum += fsymb[i] * fsymb[i] / NSYM;
        }
        float fac = sqrt(f2sum - fsum * fsum);
        for (int i = 0; i < NSYM; i++) {
            fsymb[i] = symfac * fsymb[i] / fac;
            if (fsymb[i] > 127) fsymb[i] = 127.0;
            if (fsymb[i] < -128) fsymb[i] = -128.0;
            symbols[i] = fsymb[i] + 128;
        }
        return;
    }
    return;
}


/* symbol-by-symbol signal subtraction */
void subtract_signal(float *id,
                     float *qd,
                     long  np,
                     float f0,
                     int   shift,
                     float drift,
                     unsigned char *channel_symbols) {

    float c0[NSPERSYM], s0[NSPERSYM];

    for (int i = 0; i < NSYM; i++) {
        float fp = f0 + ((float)drift / 2.0) * ((float)i - (float)NBITS) / (float)NBITS;

        float dphi  = TWOPIDT * (fp + ((float)channel_symbols[i] - 1.5) * DF);
        float cdphi = cosf(dphi);
        float sdphi = sinf(dphi);

        c0[0] = 1;
        s0[0] = 0;

        for (int j = 1; j < NSPERSYM; j++) {
            c0[j] = c0[j - 1] * cdphi - s0[j - 1] * sdphi;
            s0[j] = c0[j - 1] * sdphi + s0[j - 1] * cdphi;
        }

        float i0 = 0.0;
        float q0 = 0.0;

        for (int j = 0; j < NSPERSYM; j++) {
            int k = shift + i * NSPERSYM + j;
            if ((k > 0) & (k < np)) {
                i0 = i0 + id[k] * c0[j] + qd[k] * s0[j];
                q0 = q0 - id[k] * s0[j] + qd[k] * c0[j];
            }
        }

        // subtract the signal here.
        i0 = i0 / (float)NSPERSYM;  // will be wrong for partial symbols at the edges...
        q0 = q0 / (float)NSPERSYM;

        for (int j = 0; j < NSPERSYM; j++) {
            int k = shift + i * NSPERSYM + j;
            if ((k > 0) & (k < np)) {
                id[k] = id[k] - (i0 * c0[j] - q0 * s0[j]);
                qd[k] = qd[k] - (q0 * c0[j] + i0 * s0[j]);
            }
        }
    }
    return;
}


/* Subtract the coherent component of a signal */
void subtract_signal2(float *id,
                      float *qd,
                      long np,
                      float f0,
                      int shift,
                      float drift,
                      unsigned char *channel_symbols) {

    float phi = 0.0;
    const int nfilt = 360;  // nfilt must be even number.

    float refi[SIGNAL_SAMPLES] = {0}, refq[SIGNAL_SAMPLES] = {0},
          ci[SIGNAL_SAMPLES]   = {0}, cq[SIGNAL_SAMPLES]   = {0},
          cfi[SIGNAL_SAMPLES]  = {0}, cfq[SIGNAL_SAMPLES]  = {0};

    /******************************************************************************
     Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
     Reference is:                       r(t) = exp( j*phi(t) )
     Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
     so c(t) has phase angle theta-phi
     Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
     *******************************************************************************/

    /* create reference wspr signal vector, centered on f0. */
    for (int i = 0; i < NSYM; i++) {
        float cs = (float)channel_symbols[i];

        float dphi = TWOPIDT * (f0 + (drift / 2.0) * ((float)i - (float)NSYM / 2.0) / ((float)NSYM / 2.0) + (cs - 1.5) * DF);

        for (int j = 0; j < NSPERSYM; j++) {
            int ii = NSPERSYM * i + j;
            refi[ii] = cosf(phi);  // cannot precompute sin/cos because dphi is changing
            refq[ii] = sinf(phi);
            phi = phi + dphi;
        }
    }

    float w[nfilt], norm = 0, partialsum[nfilt];

    /* lowpass filter and remove startup transient */
    for (int i = 0; i < nfilt; i++) {
        partialsum[i] = 0.0;
    }
    for (int i = 0; i < nfilt; i++) {
        w[i] = sinf(M_PI * (float)i / (float)(nfilt - 1));
        norm = norm + w[i];
    }
    for (int i = 0; i < nfilt; i++) {
        w[i] = w[i] / norm;
    }
    for (int i = 1; i < nfilt; i++) {
        partialsum[i] = partialsum[i - 1] + w[i];
    }

    // s(t) * conjugate(r(t))
    // beginning of first symbol in reference signal is at i=0
    // beginning of first symbol in received data is at shift value.
    // filter transient lasts nfilt samples
    // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    for (int i = 0; i < NSYM * NSPERSYM; i++) {
        int k = shift + i;
        if ((k > 0) && (k < np)) {
            ci[i + nfilt] = id[k] * refi[i] + qd[k] * refq[i];
            cq[i + nfilt] = qd[k] * refi[i] - id[k] * refq[i];
        }
    }

    // LPF
    for (int i = nfilt / 2; i < SIGNAL_SAMPLES - nfilt / 2; i++) {
        cfi[i] = 0.0;
        cfq[i] = 0.0;
        for (int j = 0; j < nfilt; j++) {
            cfi[i] = cfi[i] + w[j] * ci[i - nfilt / 2 + j];
            cfq[i] = cfq[i] + w[j] * cq[i - nfilt / 2 + j];
        }
    }

    // subtract c(t)*r(t) here
    // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq+cq*refi)
    // beginning of first symbol in reference signal is at i=nfilt
    // beginning of first symbol in received data is at shift value.
    for (int i = 0; i < NSYM * NSPERSYM; i++) {
        if (i < nfilt / 2) {  // take care of the end effect (LPF step response) here
            norm = partialsum[nfilt / 2 + i];
        } else if (i > (NSYM * NSPERSYM - 1 - nfilt / 2)) {
            norm = partialsum[nfilt / 2 + NSYM * NSPERSYM - 1 - i];
        } else {
            norm = 1.0;
        }
        int k = shift + i;
        int j = i + nfilt;
        if ((k > 0) && (k < np)) {
            id[k] = id[k] - (cfi[j] * refi[i] - cfq[j] * refq[i]) / norm;
            qd[k] = qd[k] - (cfi[j] * refq[i] + cfq[j] * refi[i]) / norm;
        }
    }
    return;
}


int wspr_decode(float  *idat, 
                float  *qdat, 
                int    samples,
                struct decoder_options options, 
                struct decoder_results *decodes,
                int    *n_results) {

    /* Parameters used for performance-tuning */
    float minsync1  = 0.10;                    // First sync limit
    float minsync2  = 0.12;                    // Second sync limit
    int   iifac     = 3;                       // Step size in final DT peakup
    int   symfac    = 50;                      // Soft-symbol normalizing factor
    int   maxdrift  = 4;                       // Maximum (+/-) drift
    float minrms    = 52.0 * (symfac / 64.0);  // Final test for plausible decoding
    int   delta     = 60;                      // Fano threshold step
    int   maxcycles = 10000;                   // Fano timeout limit
    float fmin      = -110.0;
    float fmax      =  110.0;

    /* Search live parameters */
    float fstep;
    int   lagmin;
    int   lagmax;
    int   lagstep;
    int   ifmin;
    int   ifmax;

    /* Decoder flags */
    int   worth_a_try;
    int   uniques = 0;

    /* CPU usage stats */
    uint32_t metric, cycles, maxnp;

    /* Candidates properties */
    struct cand candidates[200];

    /* Decoded candidate */
    uint8_t symbols[NBITS * 2] = {0};
    uint8_t decdata[(NBITS + 7) / 8] = {0};
    int8_t  message[12] = {0};

    /* Results */
    char  callsign[13] = {0};
    char  call_loc_pow[23] = {0};
    char  call[13] = {0};
    char  loc[7] = {0};
    char  pwr[3] = {0};
    float allfreqs[100] = {0};
    char  allcalls[100][13] = {0};

    /* Setup metric table */
    int32_t mettab[2][256];
    float bias = 0.45;
    for (int i = 0; i < 256; i++) {
        mettab[0][i] = roundf(10.0 * (metric_tables[2][i] - bias));
        mettab[1][i] = roundf(10.0 * (metric_tables[2][255 - i] - bias));
    }

    /* Setup/Load hash tables */
    FILE  *fhash;
    int   nh;
    char  hashtab[32768 * 13] = {0};
    char  loctab[32768 * 5] = {0};

    if (options.usehashtable) {
        char line[80], hcall[12], hgrid[5];;
        if ((fhash = fopen("hashtable.txt", "r+"))) {
            while (fgets(line, sizeof(line), fhash) != NULL) {
                hgrid[0] = '\0';
                sscanf(line, "%d %s %s", &nh, hcall, hgrid);
                strcpy(hashtab + nh * 13, hcall);
                if (strlen(hgrid) > 0) strcpy(loctab + nh * 5, hgrid);
            }
            fclose(fhash);
        }
    }

    /* FFT buffer (512 bins) */
    fftwf_complex *fftin, *fftout;
    fftin  = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 512);
    fftout = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 512);
    PLAN   = fftwf_plan_dft_1d(512, fftin, fftout, FFTW_FORWARD, PATIENCE);

    /* Recover FFTW optimization settings */
    FILE *fp_fftw_wisdom_file;
    if ((fp_fftw_wisdom_file = fopen("fftw_wisdom.dat", "r"))) {  // Open FFTW wisdom
        fftwf_import_wisdom_from_file(fp_fftw_wisdom_file);
        fclose(fp_fftw_wisdom_file);
    }

    /* Hann function */
    float hann[512];
    for (int i = 0; i < 512; i++) {
        hann[i] = sinf(0.006147931 * i);
    }

    /* FFT output alloc */
    const int blocks = 4 * floor(samples / 512) - 1;
    float ps[512][blocks];
    memset(ps, 0.0, sizeof(float) * 512 * blocks);

    /* Main loop starts here */
    for (int ipass = 0; ipass < options.npasses; ipass++) {
        if (ipass == 1 && uniques == 0)
            break;
        if (ipass < 2) {
            maxdrift = 4;
            minsync2 = 0.12;
        }
        if (ipass == 2) {
            maxdrift = 0;    // no drift for smaller frequency estimator variance
            minsync2 = 0.10;
        }

        /* Compute FFT
         * FFT over 2 symbols, stepped by half symbols
         */
        for (int i = 0; i < blocks; i++) {
            /* Load samples */
            for (int j = 0; j < 512; j++) {
                int k = i * 128 + j;
                fftin[j][0] = idat[k] * hann[j];
                fftin[j][1] = qdat[k] * hann[j];
            }

            fftwf_execute(PLAN);

            /* Recover frequencies */
            for (int j = 0; j < 512; j++) {
                int k = j + 256;
                if (k > 511)
                    k = k - 512;
                ps[j][i] = fftout[k][0] * fftout[k][0] + fftout[k][1] * fftout[k][1];
            }
        }

        // Compute average spectrum
        float psavg[512] = {0};
        for (int i = 0; i < blocks; i++) {
            for (int j = 0; j < 512; j++) {
                psavg[j] += ps[j][i];
            }
        }

        // Already restricted by previous FIR
        // Smooth with 7-point window and limit spectrum to +/-150 Hz
        int32_t window[7] = {1, 1, 1, 1, 1, 1, 1};
        float smspec[411];
        for (int i = 0; i < 411; i++) {
            smspec[i] = 0.0;
            for (int j = -3; j <= 3; j++) {
                int k = 256 - 205 + i + j;
                smspec[i] += window[j + 3] * psavg[k];
            }
        }

        // Sort spectrum values, then pick off noise level as a percentile
        float tmpsort[411];
        for (int j = 0; j < 411; j++) {
            tmpsort[j] = smspec[j];
        }
        qsort(tmpsort, 411, sizeof(float), floatcomp);

        // Noise level of spectrum is estimated as 123/411= 30'th percentile
        float noise_level = tmpsort[122];

        /* Renormalize spectrum so that (large) peaks represent an estimate of snr.
         * We know from experience that threshold snr is near -7dB in wspr bandwidth,
         * corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
         * The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15. */

        float min_snr = powf(10.0, -8.0 / 10.0);  // this is min snr in wspr bw
        float snr_scaling_factor = 26.3;

        for (int j = 0; j < 411; j++) {
            smspec[j] = smspec[j] / noise_level - 1.0;
            if (smspec[j] < min_snr) smspec[j] = 0.1 * min_snr;
            continue;
        }

        // Find all local maxima in smoothed spectrum.
        for (int i = 0; i < 200; i++) {
            candidates[i].freq  = 0.0;
            candidates[i].snr   = 0.0;
            candidates[i].drift = 0.0;
            candidates[i].shift = 0;
            candidates[i].sync  = 0.0;
        }

        int npk = 0;
        unsigned char candidate;
        for (int j = 1; j < 410; j++) {
            candidate = (smspec[j] > smspec[j - 1]) &&
                        (smspec[j] > smspec[j + 1]) &&
                        (npk < 200);
            if (candidate) {
                candidates[npk].freq = (j - 205) * (DF / 2.0);
                candidates[npk].snr = 10.0 * log10f(smspec[j]) - snr_scaling_factor;
                npk++;
            }
        }

        // Don't waste time on signals outside of the range [fmin,fmax].
        int i = 0;
        for (int j = 0; j < npk; j++) {
            if (candidates[j].freq >= fmin && candidates[j].freq <= fmax) {
                candidates[i] = candidates[j];
                i++;
            }
        }
        npk = i;

        // bubble sort on snr, bringing freq along for the ride
        struct cand tmp;
        for (int pass = 1; pass <= npk - 1; pass++) {
            for (int k = 0; k < npk - pass; k++) {
                if (candidates[k].snr < candidates[k + 1].snr) {
                    tmp = candidates[k];
                    candidates[k] = candidates[k + 1];
                    candidates[k + 1] = tmp;
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
        for (int j = 0; j < npk; j++) {  // For each candidate...
            float sync, sync_max = -1e30;
            int if0 = candidates[j].freq / (DF / 2.0) + NSPERSYM;
            for (int ifr = if0 - 1; ifr <= if0 + 1; ifr++) {                      // Freq search
                for (int k0 = -10; k0 < 22; k0++) {                               // Time search
                    for (int idrift = -maxdrift; idrift <= maxdrift; idrift++) {  // Drift search
                        float ss = 0.0;
                        float pow = 0.0;
                        for (int k = 0; k < NSYM; k++) {  // Sum over symbols
                            int ifd = ifr + ((float)k - (float)NBITS) / (float)NBITS * ((float)idrift) / DF;
                            int kindex = k0 + 2 * k;
                            if (kindex < blocks) {
                                float p0 = sqrtf(ps[ifd - 3][kindex]);
                                float p1 = sqrtf(ps[ifd - 1][kindex]);
                                float p2 = sqrtf(ps[ifd + 1][kindex]);
                                float p3 = sqrtf(ps[ifd + 3][kindex]);

                                ss = ss + (2 * pr3vector[k] - 1) * ((p1 + p3) - (p0 + p2));
                                pow = pow + p0 + p1 + p2 + p3;
                                sync = ss / pow;
                            }
                        }
                        if (sync > sync_max) {  // Save coarse parameters
                            sync_max = sync;
                            candidates[j].shift = 128 * (k0 + 1);
                            candidates[j].drift = idrift;
                            candidates[j].freq = (ifr - NSPERSYM) * (DF / 2.0);
                            candidates[j].sync = sync;
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

        for (int j = 0; j < npk; j++) {
            memset(callsign, 0, sizeof(char) * 13);
            memset(call_loc_pow, 0, sizeof(char) * 23);
            memset(call, 0, sizeof(char) * 13);
            memset(loc, 0, sizeof(char) * 7);
            memset(pwr, 0, sizeof(char) * 3);

            float freq  = candidates[j].freq;
            float drift = candidates[j].drift;
            float sync  = candidates[j].sync;
            int   shift = candidates[j].shift;

            // Search for best sync lag (mode 0)
            fstep = 0.0;
            ifmin = 0;
            ifmax = 0;
            lagmin = shift - 128;
            lagmax = shift + 128;
            lagstep = 8;
            if (options.quickmode)
                lagstep = 16;
            sync_and_demodulate(idat, qdat, samples, symbols, &freq, ifmin, ifmax, fstep, &shift,
                                lagmin, lagmax, lagstep, &drift, symfac, &sync, 0);

            // Search for frequency peak (mode 1)
            fstep = 0.1;
            ifmin = -2;
            ifmax = 2;
            sync_and_demodulate(idat, qdat, samples, symbols, &freq, ifmin, ifmax, fstep, &shift,
                                lagmin, lagmax, lagstep, &drift, symfac, &sync, 1);

            candidates[j].freq  = freq;
            candidates[j].shift = shift;
            candidates[j].drift = drift;
            candidates[j].sync  = sync;

            if (sync > minsync1) {
                worth_a_try = 1;
            } else {
                worth_a_try = 0;
            }

            int idt = 0, ii = 0;
            int not_decoded = 1;
            while (worth_a_try && not_decoded && idt <= (128 / iifac)) {
                ii = (idt + 1) / 2;
                if (idt % 2 == 1) ii = -ii;
                ii = iifac * ii;
                int jiggered_shift = shift + ii;

                // Use mode 2 to get soft-decision symbols
                sync_and_demodulate(idat, qdat, samples, symbols, &freq, ifmin, ifmax, fstep,
                                    &jiggered_shift, lagmin, lagmax, lagstep, &drift, symfac,
                                    &sync, 2);
                float sq = 0.0;
                for (i = 0; i < NSYM; i++) {
                    float y = (float)symbols[i] - 128.0;
                    sq += y * y;
                }
                float rms = sqrtf(sq / (float)NSYM);

                if ((sync > minsync2) && (rms > minrms)) {
                    deinterleave(symbols);
                    not_decoded = fano(&metric, &cycles, &maxnp, decdata, symbols, NBITS,
                                       mettab, delta, maxcycles);
                }
                idt++;
                if (options.quickmode) 
                    break;
            }

            if (worth_a_try && !not_decoded) {
                for (i = 0; i < 11; i++) {
                    if (decdata[i] > 127) {
                        message[i] = decdata[i] - 256;
                    } else {
                        message[i] = decdata[i];
                    }
                }

                // Unpack the decoded message, update the hashtable, apply
                // sanity checks on grid and power, and return
                // call_loc_pow string and also callsign (for de-duping).
                int32_t noprint = unpk_(message, hashtab, loctab, call_loc_pow, call, loc, pwr, callsign);
                if (options.subtraction && (ipass == 0) && !noprint) {
                    unsigned char channel_symbols[NSYM];

                    if (get_wspr_channel_symbols(call_loc_pow, hashtab, loctab, channel_symbols)) {
                        subtract_signal2(idat, qdat, samples, freq, shift, drift, channel_symbols);
                    } else {
                        break;
                    }
                }

                // Avoid this incorrect pattern
                if (!strcmp(loc, "A000AA"))
                    break;

                // Remove dupes (same callsign and freq within 3 Hz)
                int32_t dupe = 0;
                for (i = 0; i < uniques; i++) {
                    if (!strcmp(callsign, allcalls[i]) && (fabs(freq - allfreqs[i]) < 3.0))
                        dupe = 1;
                }

                if (!dupe) {
                    strcpy(allcalls[uniques], callsign);
                    allfreqs[uniques] = freq;
                    uniques++;

                    double dialfreq = (double)options.freq / 1e6;
                    double freq_print = dialfreq + (1500.0 + freq) / 1e6;

                    decodes[uniques - 1].sync   = candidates[j].sync;
                    decodes[uniques - 1].snr    = candidates[j].snr;
                    decodes[uniques - 1].dt     = shift * DT - 2.0;
                    decodes[uniques - 1].freq   = freq_print;
                    decodes[uniques - 1].drift  = drift;
                    decodes[uniques - 1].cycles = cycles;
                    decodes[uniques - 1].jitter = ii;
                    strcpy(decodes[uniques - 1].message, call_loc_pow);
                    strcpy(decodes[uniques - 1].call, call);
                    strcpy(decodes[uniques - 1].loc, loc);
                    strcpy(decodes[uniques - 1].pwr, pwr);
                }
            }
        }
    }

    /* Sort the result */
    struct decoder_results temp;
    for (int j = 1; j <= uniques - 1; j++) {
        for (int k = 0; k < uniques - j; k++) {
            if (decodes[k].snr < decodes[k + 1].snr) {
                temp = decodes[k];
                decodes[k] = decodes[k + 1];
                decodes[k + 1] = temp;
            }
        }
    }

    /* Return number of spots to the calling fct */
    *n_results = uniques;

    fftwf_free(fftin);
    fftwf_free(fftout);

    if ((fp_fftw_wisdom_file = fopen("fftw_wisdom.dat", "w"))) {
        fftwf_export_wisdom_to_file(fp_fftw_wisdom_file);
        fclose(fp_fftw_wisdom_file);
    }

    fftwf_destroy_plan(PLAN);

    if (options.usehashtable) {
        fhash = fopen("hashtable.txt", "w");
        for (int i = 0; i < 32768; i++) {
            if (strncmp(hashtab + i * 13, "\0", 1) != 0) {
                fprintf(fhash, "%5d %s %s\n", i, hashtab + i * 13, loctab + i * 5);
            }
        }
        fclose(fhash);
    }

    return 0;
}
