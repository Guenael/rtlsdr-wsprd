/*
 * FreeBSD License
 * Copyright (c) 2016, Guenael
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
static void *wsprDecoder(void *arg);
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
