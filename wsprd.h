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


#pragma once


/* Option & config of decoder (Shared with the wsprd code) */
struct decoder_options {
    uint32_t freq;         // Dial frequency
    char     rcall[13];    // Callsign of the RX station
    char     rloc[7];      // Locator of the RX station
    char     date[7];      // Date & time of the processes samples
    char     uttime[5];    //  ''
    uint32_t quickmode;    // Decoder option & tweak
    uint32_t usehashtable; //  ''
    uint32_t npasses;      //  ''
    uint32_t subtraction;  //  ''
};


struct decoder_results {
    float    freq;
    float    sync;
    float    snr;
    float    dt;
    float    drift;
    int32_t  jitter;
    char     message[23];
    char     call[13];
    char     loc[7];
    char     pwr[3];
    uint32_t cycles;
};


void sync_and_demodulate(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1, float fstep,
                         int *shift1, int lagmin, int lagmax, int lagstep,
                         float *drift1, int symfac, float *sync, int mode);
void subtract_signal(float *id, float *qd, long np,
                     float f0, int shift0, float drift0, unsigned char* channel_symbols);
void subtract_signal2(float *id, float *qd, long np,
                      float f0, int shift0, float drift0, unsigned char* channel_symbols);
int wspr_decode(float *idat, float *qdat, unsigned int npoints,
                struct decoder_options options, struct decoder_results *decodes, int *n_results);
