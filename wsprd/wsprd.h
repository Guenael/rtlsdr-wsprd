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
    uint32_t  freq;          // Dial frequency
    char      rcall[13];     // Callsign of the RX station
    char      rloc[7];       // Locator of the RX station
    uint32_t  quickmode;     // Decoder option & tweak
    uint32_t  usehashtable;  //  ''
    uint32_t  npasses;       //  ''
    uint32_t  subtraction;   //  ''
};

struct cand {
    float    freq;
    float    snr;
    int      shift;
    float    drift;
    float    sync;
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

// struct result {
//     char date[7];
//     char time[5];
//     float sync;
//     float snr;
//     float dt;
//     double freq;
//     char message[23];
//     float drift;
//     unsigned int cycles;
//     int jitter;
//     int blocksize;
//     unsigned int metric;
//     int nhardmin;
//     int ipass;
//     int decodetype;
// };

void sync_and_demodulate(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
                         int *shift1, int lagmin, int lagmax, int lagstep,
                         float *drift1, int symfac, float *sync, int mode);
void noncoherent_sequence_detection(float *id, float *qd, long np,
                                    unsigned char *symbols, float *f1, int *shift1,
                                    float *drift1, int symfac, int *nblocksize, int *bitmetric);
void subtract_signal(float *id, float *qd, long np,
                     float f0, int shift0, float drift0, unsigned char *channel_symbols);
void subtract_signal2(float *id, float *qd, long np,
                      float f0, int shift0, float drift0, unsigned char *channel_symbols);
unsigned int count_hard_errors(unsigned char *symbols, unsigned char *channel_symbols);
int32_t wspr_decode(float *idat, float *qdat, uint32_t npoints,
                    struct decoder_options options, struct decoder_results *decodes,
                    int32_t *n_results);








