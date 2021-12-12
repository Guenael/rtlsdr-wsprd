/*
 * rtlsrd-wsprd, WSPR daemon for RTL receivers, Guenael Jouchet (VA2GKA)
 * Copyright (C) 2016-2021, Guenael Jouchet (VA2GKA)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <unistd.h>

#ifndef bool
    typedef uint32_t bool;
    #define true 1
    #define false 0
#endif


struct receiver_state {
    /* Variables used for stop conditions */
    bool exit_flag;
    bool decode_flag;

    /* Buffer used for sampling */
    float *iSamples;
    float *qSamples;

    /* Simple index */
    uint32_t iqIndex;
};


/* Option & config of the receiver */
struct receiver_options {
    uint32_t dialfreq;
    uint32_t realfreq;
    int32_t  gain;
    int32_t  autogain;
    int32_t  ppm;
    int32_t  shift;
    int32_t  upconverter;
    int32_t  directsampling;
    int32_t  maxloop;
    int32_t  device;
    bool     selftest;
    bool     writefile;
    bool     readfile;
    char     filename[33];
    char     date[7];
    char     uttime[5];
};


static void rtlsdr_callback(unsigned char *samples, uint32_t sigLenght, void *ctx);
static void *rtlsdr_rx(void *arg);
void postSpots(uint32_t n_results);
void printSpots(uint32_t n_results);
void saveSample(float *iSamples, float *qSamples);
static void *decoder(void *arg);
double atofs(char *s);
int32_t parse_u64(char *s, uint64_t *const value);
void initSampleStorage();
void initDecoder_options();
void initrx_options();
void sigint_callback_handler(int signum);
int32_t readRawIQfile(float *iSamples, float *qSamples, char *filename);
int32_t writeRawIQfile(float *iSamples, float *qSamples, char *filename);
void decodeRecordedFile(char *filename);
float whiteGaussianNoise(float factor);
int32_t decoderSelfTest();
void usage(void);
int32_t readRawIQfile(float *iSamples, float *qSamples, char *filename);
int32_t writeRawIQfile(float *iSamples, float *qSamples, char *filename);
