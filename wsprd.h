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
    double   freq;
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
                         uint8_t *symbols, float *freq1, float fstep,
                         int32_t *shift1, int32_t lagmin, int32_t lagmax, int32_t lagstep,
                         float *drift1, int32_t symfac, float *sync, int32_t mode);
void subtract_signal(float *id, float *qd, long np,
                     float f0, int32_t shift0, float drift0, uint8_t* channel_symbols);
void subtract_signal2(float *id, float *qd, long np,
                      float f0, int32_t shift0, float drift0, uint8_t* channel_symbols);
int32_t wspr_decode(float *idat, float *qdat, uint32_t npoints,
                    struct decoder_options options, struct decoder_results *decodes,
                    int32_t *n_results);
