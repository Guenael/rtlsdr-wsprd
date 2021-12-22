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

#define NFILT               256
#define NSIG                NSYM * NSPERSYM


/* Possible PATIENCE options: FFTW_ESTIMATE, FFTW_ESTIMATE_PATIENT, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE */
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


//***************************************************************************
void sync_and_demodulate(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
                         int *shift1, int lagmin, int lagmax, int lagstep,
                         float *drift1, int symfac, float *sync, int mode) {
    /***********************************************************************
     * mode = 0: no frequency or drift search. find best time lag.          *
     *        1: no time lag or drift search. find best frequency.          *
     *        2: no frequency or time lag search. calculate soft-decision   *
     *           symbols using passed frequency and shift.                  *
     ************************************************************************/

    static float fplast = -10000.0;
    float f0 = 0.0, fp, ss, 
          fbest = 0.0, 
          fsum = 0.0, 
          f2sum = 0.0, 
          fsymb[NSYM];
    int best_shift = 0;
    float i0[NSYM], q0[NSYM], 
          i1[NSYM], q1[NSYM], 
          i2[NSYM], q2[NSYM], 
          i3[NSYM], q3[NSYM];
    float p0, p1, p2, p3, cmet, totp, syncmax, fac;
    float c0[NSPERSYM], s0[NSPERSYM], 
          c1[NSPERSYM], s1[NSPERSYM], 
          c2[NSPERSYM], s2[NSPERSYM], 
          c3[NSPERSYM], s3[NSPERSYM];
    float dphi0, cdphi0, sdphi0, 
          dphi1, cdphi1, sdphi1, 
          dphi2, cdphi2, sdphi2,
          dphi3, cdphi3, sdphi3;

    syncmax = -1e30;
    if (mode == 0) {
        ifmin = 0;
        ifmax = 0;
        fstep = 0.0;
        f0 = *f1;
    }
    if (mode == 1) {
        lagmin = *shift1;
        lagmax = *shift1;
        f0 = *f1;
    }
    if (mode == 2) {
        lagmin = *shift1;
        lagmax = *shift1;
        ifmin = 0;
        ifmax = 0;
        f0 = *f1;
    }

    for (int ifreq = ifmin; ifreq <= ifmax; ifreq++) {
        f0 = *f1 + ifreq * fstep;
        for (int lag = lagmin; lag <= lagmax; lag = lag + lagstep) {
            ss = 0.0;
            totp = 0.0;
            for (int i = 0; i < NSYM; i++) {
                fp = f0 + (*drift1 / 2.0) * ((float)i - (float)NBITS) / (float)NBITS;
                if (i == 0 || (fp != fplast)) {  // only calculate sin/cos if necessary
                    dphi0 = TWOPIDT * (fp - DF15);
                    cdphi0 = cos(dphi0);
                    sdphi0 = sin(dphi0);

                    dphi1 = TWOPIDT * (fp - DF05);
                    cdphi1 = cos(dphi1);
                    sdphi1 = sin(dphi1);

                    dphi2 = TWOPIDT * (fp + DF05);
                    cdphi2 = cos(dphi2);
                    sdphi2 = sin(dphi2);

                    dphi3 = TWOPIDT * (fp + DF15);
                    cdphi3 = cos(dphi3);
                    sdphi3 = sin(dphi3);

                    c0[0] = 1;
                    s0[0] = 0;
                    c1[0] = 1;
                    s1[0] = 0;
                    c2[0] = 1;
                    s2[0] = 0;
                    c3[0] = 1;
                    s3[0] = 0;

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

                i0[i] = 0.0;
                q0[i] = 0.0;
                i1[i] = 0.0;
                q1[i] = 0.0;
                i2[i] = 0.0;
                q2[i] = 0.0;
                i3[i] = 0.0;
                q3[i] = 0.0;

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
                p0 = i0[i] * i0[i] + q0[i] * q0[i];
                p1 = i1[i] * i1[i] + q1[i] * q1[i];
                p2 = i2[i] * i2[i] + q2[i] * q2[i];
                p3 = i3[i] * i3[i] + q3[i] * q3[i];

                p0 = sqrt(p0);
                p1 = sqrt(p1);
                p2 = sqrt(p2);
                p3 = sqrt(p3);

                totp = totp + p0 + p1 + p2 + p3;
                cmet = (p1 + p3) - (p0 + p2);
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
    }      // freq loop

    if (mode <= 1) {  // Send best params back to caller
        *sync = syncmax;
        *shift1 = best_shift;
        *f1 = fbest;
        return;
    }

    if (mode == 2) {
        *sync = syncmax;
        for (int i = 0; i < NSYM; i++) {  // Normalize the soft symbols
            fsum = fsum + fsymb[i] / NSYM;
            f2sum = f2sum + fsymb[i] * fsymb[i] / NSYM;
        }
        fac = sqrt(f2sum - fsum * fsum);
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

void noncoherent_sequence_detection(float *id, float *qd, long np,
                                    unsigned char *symbols, float *f1, int *shift1,
                                    float *drift1, int symfac, int *nblocksize, int *bitmetric) {
    /************************************************************************
     *  Noncoherent sequence detection for wspr.                            *
     *  Allowed block lengths are nblock=1,2,3,6, or 9 symbols.             *
     *  Longer block lengths require longer channel coherence time.         *
     *  The whole block is estimated at once.                               *
     *  nblock=1 corresponds to noncoherent detection of individual symbols *
     *     like the original wsprd symbol demodulator.                      *
     ************************************************************************/
    static float fplast = -10000.0;

    int i, j, k, lag, itone, ib, b, nblock, nseq, imask;
    float xi[512], xq[512];
    float is[4][NSYM], qs[4][NSYM], 
          cf[4][NSYM], sf[4][NSYM], 
          cm, sm, cmp, smp;
    float p[512], fac, xm1, xm0;
    float c0[NSPERSYM+1], s0[NSPERSYM+1], 
          c1[NSPERSYM+1], s1[NSPERSYM+1], 
          c2[NSPERSYM+1], s2[NSPERSYM+1], 
          c3[NSPERSYM+1], s3[NSPERSYM+1];
    float dphi0, cdphi0, sdphi0, 
          dphi1, cdphi1, sdphi1, 
          dphi2, cdphi2, sdphi2,
          dphi3, cdphi3, sdphi3;
    float f0, fp, fsum = 0.0, f2sum = 0.0, fsymb[NSYM];

    f0 = *f1;
    lag = *shift1;
    nblock = *nblocksize;
    nseq = 1 << nblock;
    int bitbybit = *bitmetric;

    for (i = 0; i < NSYM; i++) {
        fp = f0 + (*drift1 / 2.0) * ((float)i - (float)NBITS) / (float)NBITS;
        if (i == 0 || (fp != fplast)) {  // only calculate sin/cos if necessary
            dphi0 = TWOPIDT * (fp - DF15);
            cdphi0 = cos(dphi0);
            sdphi0 = sin(dphi0);

            dphi1 = TWOPIDT * (fp - DF05);
            cdphi1 = cos(dphi1);
            sdphi1 = sin(dphi1);

            dphi2 = TWOPIDT * (fp + DF05);
            cdphi2 = cos(dphi2);
            sdphi2 = sin(dphi2);

            dphi3 = TWOPIDT * (fp + DF15);
            cdphi3 = cos(dphi3);
            sdphi3 = sin(dphi3);

            c0[0] = 1;
            s0[0] = 0;
            c1[0] = 1;
            s1[0] = 0;
            c2[0] = 1;
            s2[0] = 0;
            c3[0] = 1;
            s3[0] = 0;

            for (j = 1; j < (NSPERSYM+1); j++) {
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

        cf[0][i] = c0[NSPERSYM];
        sf[0][i] = s0[NSPERSYM];
        cf[1][i] = c1[NSPERSYM];
        sf[1][i] = s1[NSPERSYM];
        cf[2][i] = c2[NSPERSYM];
        sf[2][i] = s2[NSPERSYM];
        cf[3][i] = c3[NSPERSYM];
        sf[3][i] = s3[NSPERSYM];

        is[0][i] = 0.0;
        qs[0][i] = 0.0;
        is[1][i] = 0.0;
        qs[1][i] = 0.0;
        is[2][i] = 0.0;
        qs[2][i] = 0.0;
        is[3][i] = 0.0;
        qs[3][i] = 0.0;

        for (j = 0; j < NSPERSYM; j++) {
            k = lag + i * NSPERSYM + j;
            if ((k > 0) && (k < np)) {
                is[0][i] = is[0][i] + id[k] * c0[j] + qd[k] * s0[j];
                qs[0][i] = qs[0][i] - id[k] * s0[j] + qd[k] * c0[j];
                is[1][i] = is[1][i] + id[k] * c1[j] + qd[k] * s1[j];
                qs[1][i] = qs[1][i] - id[k] * s1[j] + qd[k] * c1[j];
                is[2][i] = is[2][i] + id[k] * c2[j] + qd[k] * s2[j];
                qs[2][i] = qs[2][i] - id[k] * s2[j] + qd[k] * c2[j];
                is[3][i] = is[3][i] + id[k] * c3[j] + qd[k] * s3[j];
                qs[3][i] = qs[3][i] - id[k] * s3[j] + qd[k] * c3[j];
            }
        }
    }

    for (i = 0; i < NSYM; i = i + nblock) {
        for (j = 0; j < nseq; j++) {
            xi[j] = 0.0;
            xq[j] = 0.0;
            cm = 1;
            sm = 0;
            for (ib = 0; ib < nblock; ib++) {
                b = (j & (1 << (nblock - 1 - ib))) >> (nblock - 1 - ib);
                itone = pr3vector[i + ib] + 2 * b;
                xi[j] = xi[j] + is[itone][i + ib] * cm + qs[itone][i + ib] * sm;
                xq[j] = xq[j] + qs[itone][i + ib] * cm - is[itone][i + ib] * sm;
                cmp = cf[itone][i + ib] * cm - sf[itone][i + ib] * sm;
                smp = sf[itone][i + ib] * cm + cf[itone][i + ib] * sm;
                cm = cmp;
                sm = smp;
            }
            p[j] = xi[j] * xi[j] + xq[j] * xq[j];
            p[j] = sqrt(p[j]);
        }
        for (ib = 0; ib < nblock; ib++) {
            imask = 1 << (nblock - 1 - ib);
            xm1 = 0.0;
            xm0 = 0.0;
            for (j = 0; j < nseq; j++) {
                if ((j & imask) != 0) {
                    if (p[j] > xm1) xm1 = p[j];
                }
                if ((j & imask) == 0) {
                    if (p[j] > xm0) xm0 = p[j];
                }
            }
            fsymb[i + ib] = xm1 - xm0;
            if (bitbybit == 1) {
                fsymb[i + ib] = fsymb[i + ib] / (xm1 > xm0 ? xm1 : xm0);
            }
        }
    }
    for (i = 0; i < NSYM; i++) {  // Normalize the soft symbols
        fsum = fsum + fsymb[i] / NSYM;
        f2sum = f2sum + fsymb[i] * fsymb[i] / NSYM;
    }
    fac = sqrt(f2sum - fsum * fsum);
    for (i = 0; i < NSYM; i++) {
        fsymb[i] = symfac * fsymb[i] / fac;
        if (fsymb[i] > 127) fsymb[i] = 127.0;
        if (fsymb[i] < -128) fsymb[i] = -128.0;
        symbols[i] = fsymb[i] + 128;
    }
    return;
}

/***************************************************************************
 symbol-by-symbol signal subtraction
 ****************************************************************************/
void subtract_signal(float *id, float *qd, long np,
                     float f0, int shift0, float drift0, unsigned char *channel_symbols) {

    float c0[NSPERSYM], s0[NSPERSYM];
    float dphi, cdphi, sdphi;

    for (int i = 0; i < NSYM; i++) {
        float fp = f0 + ((float)drift0 / 2.0) * ((float)i - (float)NBITS) / (float)NBITS;

        dphi  = TWOPIDT * (fp + ((float)channel_symbols[i] - 1.5) * DF);
        cdphi = cos(dphi);
        sdphi = sin(dphi);

        c0[0] = 1;
        s0[0] = 0;

        for (int j = 1; j < NSPERSYM; j++) {
            c0[j] = c0[j - 1] * cdphi - s0[j - 1] * sdphi;
            s0[j] = c0[j - 1] * sdphi + s0[j - 1] * cdphi;
        }

        float i0 = 0.0;
        float q0 = 0.0;

        for (int j = 0; j < NSPERSYM; j++) {
            int k = shift0 + i * NSPERSYM + j;
            if ((k > 0) & (k < np)) {
                i0 = i0 + id[k] * c0[j] + qd[k] * s0[j];
                q0 = q0 - id[k] * s0[j] + qd[k] * c0[j];
            }
        }

        // subtract the signal here.

        i0 = i0 / (float)NSPERSYM;  // will be wrong for partial symbols at the edges...
        q0 = q0 / (float)NSPERSYM;

        for (int j = 0; j < NSPERSYM; j++) {
            int k = shift0 + i * NSPERSYM + j;
            if ((k > 0) & (k < np)) {
                id[k] = id[k] - (i0 * c0[j] - q0 * s0[j]);
                qd[k] = qd[k] - (q0 * c0[j] + i0 * s0[j]);
            }
        }
    }
    return;
}


/******************************************************************************
  Subtract the coherent component of a signal
 *******************************************************************************/
void subtract_signal2(float *id,
                      float *qd,
                      long np,
                      float f0,
                      int shift0,
                      float drift0,
                      unsigned char *channel_symbols) {

    float phi = 0;
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

    // create reference wspr signal vector, centered on f0.
    //
    for (int i = 0; i < NSYM; i++) {
        float cs = (float)channel_symbols[i];

        float dphi = TWOPIDT * (f0 + (drift0 / 2.0) * ((float)i - (float)NSYM / 2.0) / ((float)NSYM / 2.0) + (cs - 1.5) * DF);

        for (int j = 0; j < NSPERSYM; j++) {
            int ii = NSPERSYM * i + j;
            refi[ii] = cos(phi);  // cannot precompute sin/cos because dphi is changing
            refq[ii] = sin(phi);
            phi = phi + dphi;
        }
    }

    float w[nfilt], norm = 0, partialsum[nfilt];

    // lowpass filter and remove startup transient
    for (int i = 0; i < nfilt; i++) {
        partialsum[i] = 0.0;
    }
    for (int i = 0; i < nfilt; i++) {
        w[i] = sin(M_PI * (float)i / (float)(nfilt - 1));
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
    // beginning of first symbol in received data is at shift0.
    // filter transient lasts nfilt samples
    // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    for (int i = 0; i < NSYM * NSPERSYM; i++) {
        int k = shift0 + i;
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
    // beginning of first symbol in received data is at shift0.
    for (int i = 0; i < NSYM * NSPERSYM; i++) {
        if (i < nfilt / 2) {  // take care of the end effect (LPF step response) here
            norm = partialsum[nfilt / 2 + i];
        } else if (i > (NSYM * NSPERSYM - 1 - nfilt / 2)) {
            norm = partialsum[nfilt / 2 + NSYM * NSPERSYM - 1 - i];
        } else {
            norm = 1.0;
        }
        int k = shift0 + i;
        int j = i + nfilt;
        if ((k > 0) && (k < np)) {
            id[k] = id[k] - (cfi[j] * refi[i] - cfq[j] * refq[i]) / norm;
            qd[k] = qd[k] - (cfi[j] * refq[i] + cfq[j] * refi[i]) / norm;
        }
    }
    return;
}


unsigned int count_hard_errors(unsigned char *symbols, unsigned char *channel_symbols) {
    int i, is;
    unsigned char cw[162];
    unsigned int nerrors;
    for (i = 0; i < 162; i++) {
        cw[i] = channel_symbols[i] >= 2 ? 1 : 0;
    }
    deinterleave(cw);
    nerrors = 0;
    for (i = 0; i < 162; i++) {
        is = symbols[i] > 127 ? 1 : 0;
        nerrors = nerrors + (is == cw[i] ? 0 : 1);
    }
    return nerrors;
}


int32_t wspr_decode(float *idat, float *qdat, uint32_t npoints,
                    struct decoder_options options, struct decoder_results *decodes,
                    int32_t *n_results) {
    int32_t i, j, k;
    uint8_t symbols[NBITS * 2] = {0};
    uint8_t decdata[(NBITS + 7) / 8] = {0};
    int8_t message[12] = {0};

    float freq0[200], freq1 = 0.0;
    float drift0[200], drift1 = 0.0;
    float sync0[200], sync1 = 0.0;
    float snr0[200];
    int32_t shift0[200], shift1 = 0;
    int32_t ifmin, ifmax;

    char callsign[13] = {0};
    char call_loc_pow[23] = {0};
    char call[13] = {0};
    char loc[7] = {0};
    char pwr[3] = {0};

    uint32_t metric, cycles, maxnp;
    int32_t worth_a_try;
    int32_t uniques = 0;
    float fmin = -110.0;
    float fmax = 110.0;

    // Hash table
    char hashtab[32768 * 13] = {0};
    char loctab[32768 * 5] = {0};
    int32_t nh;

    // Parameters used for performance-tuning:
    uint32_t maxcycles = 10000;              // Fano timeout limit
    double minsync1 = 0.10;                  // First sync limit
    double minsync2 = 0.12;                  // Second sync limit
    int32_t iifac = 3;                       // Step size in final DT peakup
    int32_t symfac = 50;                     // Soft-symbol normalizing factor
    int32_t maxdrift = 4;                    // Maximum (+/-) drift
    double minrms = 52.0 * (symfac / 64.0);  // Final test for plausible decoding
    int32_t delta = 60;                      // Fano threshold step

    // Results
    float allfreqs[100];
    char allcalls[100][13];
    memset(allfreqs, 0, sizeof(float) * 100);
    memset(allcalls, 0, sizeof(char) * 100 * 13);

    // Setup metric table
    int32_t mettab[2][256];
    float bias = 0.42;
    for (int i = 0; i < 256; i++) {
        mettab[0][i] = round(10 * (metric_tables[2][i] - bias));
        mettab[1][i] = round(10 * (metric_tables[2][255 - i] - bias));
    }

    // FFT buffer
    fftwf_complex *fftin, *fftout;

    FILE *fp_fftw_wisdom_file, *fhash;
    if ((fp_fftw_wisdom_file = fopen("fftw_wisdom.dat", "r"))) {  // Open FFTW wisdom
        fftwf_import_wisdom_from_file(fp_fftw_wisdom_file);
        fclose(fp_fftw_wisdom_file);
    }

    // Do windowed ffts over 2 symbols, stepped by half symbols
    int32_t nffts = 4 * floor(npoints / 512) - 1;
    fftin = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 512);
    fftout = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * 512);
    PLAN = fftwf_plan_dft_1d(512, fftin, fftout, FFTW_FORWARD, PATIENCE);

    float ps[512][nffts];
    float w[512];
    for (int i = 0; i < 512; i++) {
        w[i] = sinf(0.006147931 * i);
    }

    if (options.usehashtable) {
        char line[80], hcall[12];
        if ((fhash = fopen("hashtable.txt", "r+"))) {
            while (fgets(line, sizeof(line), fhash) != NULL) {
                sscanf(line, "%d %s", &nh, hcall);
                strcpy(hashtab + nh * 13, hcall);
            }
            fclose(fhash);
        }
    }

    // Main loop starts here
    for (int ipass = 0; ipass < options.npasses; ipass++) {
        if (ipass == 1 && uniques == 0)
            break;
        if (ipass == 1)  // otherwise we bog down on the second pass
            options.quickmode = 1;

        memset(ps, 0.0, sizeof(float) * 512 * nffts);
        for (int i = 0; i < nffts; i++) {
            for (j = 0; j < 512; j++) {
                k = i * 128 + j;
                fftin[j][0] = idat[k] * w[j];
                fftin[j][1] = qdat[k] * w[j];
            }

            fftwf_execute(PLAN);

            for (int j = 0; j < 512; j++) {
                k = j + 256;
                if (k > 511)
                    k = k - 512;
                ps[j][i] = fftout[k][0] * fftout[k][0] + fftout[k][1] * fftout[k][1];
            }
        }

        // Compute average spectrum
        float psavg[512] = {0};
        for (int i = 0; i < nffts; i++) {
            for (j = 0; j < 512; j++) {
                psavg[j] = psavg[j] + ps[j][i];
            }
        }

        // Smooth with 7-point window and limit spectrum to +/-150 Hz
        int32_t window[7] = {1, 1, 1, 1, 1, 1, 1};
        float smspec[411];
        for (int i = 0; i < 411; i++) {
            smspec[i] = 0.0;
            for (j = -3; j <= 3; j++) {
                k = 256 - 205 + i + j;
                smspec[i] = smspec[i] + window[j + 3] * psavg[k];
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

        float min_snr = powf(10.0, -7.0 / 10.0);  // this is min snr in wspr bw
        float snr_scaling_factor = 26.3;

        for (int j = 0; j < 411; j++) {
            smspec[j] = smspec[j] / noise_level - 1.0;
            if (smspec[j] < min_snr) smspec[j] = 0.1;
            continue;
        }

        // Find all local maxima in smoothed spectrum.
        for (int i = 0; i < 200; i++) {
            freq0[i] = 0.0;
            snr0[i] = 0.0;
            drift0[i] = 0.0;
            shift0[i] = 0;
            sync0[i] = 0.0;
        }

        int32_t npk = 0;
        for (int j = 1; j < 410; j++) {
            if ((smspec[j] > smspec[j - 1]) && (smspec[j] > smspec[j + 1]) && (npk < 200)) {
                freq0[npk] = (j - 205) * (DF / 2.0);
                snr0[npk] = 10.0 * log10f(smspec[j]) - snr_scaling_factor;
                npk++;
            }
        }

        /* Compute corrected fmin, fmax, accounting for dial frequency error
        float   dialfreq_error = 0.0; // dialfreq_error is in units of Hz
        fmin += dialfreq_error;
        fmax += dialfreq_error;
        */

        // Don't waste time on signals outside of the range [fmin,fmax].
        i = 0;
        for (int j = 0; j < npk; j++) {
            if (freq0[j] >= fmin && freq0[j] <= fmax) {
                freq0[i] = freq0[j];
                snr0[i] = snr0[j];
                i++;
            }
        }
        npk = i;

        // bubble sort on snr, bringing freq along for the ride
        float tmp;
        for (int pass = 1; pass <= npk - 1; pass++) {
            for (k = 0; k < npk - pass; k++) {
                if (snr0[k] < snr0[k + 1]) {
                    tmp = snr0[k];
                    snr0[k] = snr0[k + 1];
                    snr0[k + 1] = tmp;
                    tmp = freq0[k];
                    freq0[k] = freq0[k + 1];
                    freq0[k + 1] = tmp;
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
            float smax = -1e30;
            int if0 = freq0[j] / (DF / 2.0) + NSPERSYM;
            for (int ifr = if0 - 1; ifr <= if0 + 1; ifr++) {                      // Freq search
                for (int k0 = -10; k0 < 22; k0++) {                               // Time search
                    for (int idrift = -maxdrift; idrift <= maxdrift; idrift++) {  // Drift search
                        float ss = 0.0;
                        float pow = 0.0;
                        for (k = 0; k < NSYM; k++) {  // Sum over symbols
                            int ifd = ifr + ((float)k - (float)NBITS) / (float)NBITS * ((float)idrift) / DF;
                            int kindex = k0 + 2 * k;
                            if (kindex < nffts) {
                                float p0 = sqrtf(ps[ifd - 3][kindex]);
                                float p1 = sqrtf(ps[ifd - 1][kindex]);
                                float p2 = sqrtf(ps[ifd + 1][kindex]);
                                float p3 = sqrtf(ps[ifd + 3][kindex]);

                                ss = ss + (2 * pr3vector[k] - 1) * ((p1 + p3) - (p0 + p2));
                                pow = pow + p0 + p1 + p2 + p3;
                                sync1 = ss / pow;
                            }
                        }
                        if (sync1 > smax) {  // Save coarse parameters
                            smax = sync1;
                            shift0[j] = 128 * (k0 + 1);
                            drift0[j] = idrift;
                            freq0[j] = (ifr - NSPERSYM) * (DF / 2.0);
                            sync0[j] = sync1;
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
            memset(symbols, 0, sizeof(char) * NBITS * 2);
            memset(callsign, 0, sizeof(char) * 13);
            memset(call_loc_pow, 0, sizeof(char) * 23);
            memset(call, 0, sizeof(char) * 13);
            memset(loc, 0, sizeof(char) * 7);
            memset(pwr, 0, sizeof(char) * 3);

            freq1  = freq0[j];
            drift1 = drift0[j];
            shift1 = shift0[j];
            sync1  = sync0[j];

            // Fine search for best sync lag (mode 0)
            float fstep = 0.0;
            int32_t lagmin = shift1 - 144;
            int32_t lagmax = shift1 + 144;
            int32_t lagstep = 8;
            ifmin = 0;
            ifmax = 0;

            if (options.quickmode)
                lagstep = 16;

            sync_and_demodulate(idat, qdat, npoints, symbols, &freq1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);

            // Fine search for frequency peak (mode 1)
            fstep = 0.1;
            ifmin = 0;
            ifmin = -2;
            ifmax = 2;
            sync_and_demodulate(idat, qdat, npoints, symbols, &freq1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);

            if (sync1 > minsync1) {
                worth_a_try = 1;
            } else {
                worth_a_try = 0;
            }

            int32_t idt = 0, ii = 0, jiggered_shift;
            float y, sq, rms;
            int32_t not_decoded = 1;

            while (worth_a_try && not_decoded && idt <= (128 / iifac)) {
                ii = (idt + 1) / 2;
                if (idt % 2 == 1) ii = -ii;
                ii = iifac * ii;
                jiggered_shift = shift1 + ii;

                // Use mode 2 to get soft-decision symbols
                sync_and_demodulate(idat, qdat, npoints, symbols, &freq1, ifmin, ifmax, fstep,
                                    &jiggered_shift, lagmin, lagmax, lagstep, &drift1, symfac,
                                    &sync1, 2);
                sq = 0.0;
                for (i = 0; i < NSYM; i++) {
                    y = (float)symbols[i] - 128.0;
                    sq += y * y;
                }
                rms = sqrtf(sq / (float)NSYM);

                if ((sync1 > minsync2) && (rms > minrms)) {
                    deinterleave(symbols);
                    not_decoded = fano(&metric, &cycles, &maxnp, decdata, symbols, NBITS,
                                       mettab, delta, maxcycles);
                }
                idt++;
                if (options.quickmode) break;
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

                    //if (get_wspr_channel_symbols(call_loc_pow, hashtab, channel_symbols)) {
                    if (get_wspr_channel_symbols(call_loc_pow, hashtab, loctab, channel_symbols)) {
                        subtract_signal2(idat, qdat, npoints, freq1, shift1, drift1, channel_symbols);
                    } else {
                        break;
                    }
                }

                // Remove dupes (same callsign and freq within 3 Hz)
                int32_t dupe = 0;
                for (i = 0; i < uniques; i++) {
                    if (!strcmp(callsign, allcalls[i]) && (fabs(freq1 - allfreqs[i]) < 3.0))
                        dupe = 1;
                }

                if (!dupe) {
                    strcpy(allcalls[uniques], callsign);
                    allfreqs[uniques] = freq1;
                    uniques++;

                    double dialfreq = (double)options.freq / 1e6;
                    double freq_print = dialfreq + (1500 + freq1) / 1e6;
                    float dt_print = shift1 * DT - 2.0;

                    decodes[uniques - 1].sync = sync1;
                    decodes[uniques - 1].snr = snr0[j];
                    decodes[uniques - 1].dt = dt_print;
                    decodes[uniques - 1].freq = freq_print;
                    decodes[uniques - 1].drift = drift1;
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

    // sort the result
    struct decoder_results temp;
    for (j = 1; j <= uniques - 1; j++) {
        for (k = 0; k < uniques - j; k++) {
            if (decodes[k].snr < decodes[k + 1].snr) {
                temp = decodes[k];
                decodes[k] = decodes[k + 1];
                ;
                decodes[k + 1] = temp;
            }
        }
    }

    // Return number of spots to the calling fct
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
        for (i = 0; i < 32768; i++) {
            if (strncmp(hashtab + i * 13, "\0", 1) != 0) {
                fprintf(fhash, "%5d %s\n", i, hashtab + i * 13);
            }
        }
        fclose(fhash);
    }
    return 0;
}
