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


static void rtlsdr_callback(unsigned char *samples, uint32_t samples_count, void *ctx);
static void sigint_callback_handler(int signum);
static void *rtlsdr_rx(void *arg);
static void *decoder(void *arg);
void initSampleStorage();
void initrx_options();
void initDecoder_options();
void postSpots(uint32_t n_results);
void printSpots(uint32_t n_results);
void saveSample(float *iSamples, float *qSamples);
double atofs(char *s);
int32_t parse_u64(char *s, uint64_t *const value);
int32_t readRawIQfile(float *iSamples, float *qSamples, char *filename);
int32_t writeRawIQfile(float *iSamples, float *qSamples, char *filename);
int32_t readC2file(float *iSamples, float *qSamples, char *filename);
void decodeRecordedFile(char *filename);
float whiteGaussianNoise(float factor);
int32_t decoderSelfTest();
void usage(FILE *stream, int32_t status);
