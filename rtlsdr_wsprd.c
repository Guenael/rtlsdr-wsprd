/*
 * FreeBSD License
 * Copyright (c) 2016, Guenael
 * All rights reserved.
 *
 * This file is based on rtl-sdr project, contribution :
 *   Copyright (C) 2012 by Steve Markgraf <steve@steve-m.de>
 *   Copyright (C) 2012 by Hoernchen <la@tfc-server.de>
 *   Copyright (C) 2012 by Kyle Keen <keenerd@gmail.com>
 *   Copyright (C) 2013 by Elias Oenal <EliasOenal@gmail.com>
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


#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <curl/curl.h>
#include <pthread.h>
#include <rtl-sdr.h>

#include "rtlsdr_wsprd.h"
#include "wsprd.h"

/* TODO
 - multi device selection option
 - multispot report in one post
 - verbose option
*/


#define SIGNAL_LENGHT       116
#define SIGNAL_SAMPLE_RATE  375
#define SAMPLING_RATE       2400000
#define FS4_RATE            SAMPLING_RATE / 4                   // = 600 kHz
#define DOWNSAMPLING        SAMPLING_RATE / SIGNAL_SAMPLE_RATE  // = 6400
#define DEFAULT_BUF_LENGTH  (4 * 16384)                         // = 65536


/* Global declaration for these structs */
struct receiver_state   rx_state;
struct receiver_options rx_options;
struct decoder_options  dec_options;
struct decoder_results  dec_results[50];
static rtlsdr_dev_t *rtl_device = NULL;


/* Thread stuff for separate decoding */
struct decoder_state {
    pthread_t        thread;
    pthread_attr_t   tattr;

    pthread_rwlock_t rw;
    pthread_cond_t   ready_cond;
    pthread_mutex_t  ready_mutex;
};
struct decoder_state dec;


/* Thread stuff for separate RX (blocking function) */
struct dongle_state {
    pthread_t        thread;
};
struct dongle_state dongle;


/* Callback for each buffer received */
static void rtlsdr_callback(unsigned char *samples, uint32_t samples_count, void *ctx) {
    int8_t *sigIn = (int8_t*) samples;
    uint32_t sigLenght = samples_count;

    static uint32_t decimationIndex=0;

    /* CIC buffers */
    static int32_t  Ix1,Ix2,Qx1,Qx2;
    static int32_t  Iy1,It1y,It1z,Qy1,Qt1y,Qt1z;
    static int32_t  Iy2,It2y,It2z,Qy2,Qt2y,Qt2z;

    /* FIR compensation filter buffers */
    static float    firI[32], firQ[32];

    /* FIR compensation filter coefs
       Using : Octave/MATLAB code for generating compensation FIR coefficients
       URL : https://github.com/WestCoastDSP/CIC_Octave_Matlab
     */
    const static float zCoef[33] = {
        -0.0027772683, -0.0005058826,  0.0049745750, -0.0034059318,
        -0.0077557814,  0.0139375423,  0.0039896935, -0.0299394142,
         0.0162250643,  0.0405130860, -0.0580746013, -0.0272104968,
         0.1183705475, -0.0306029022, -0.2011241667,  0.1615898423,
         0.5000000000,
         0.1615898423, -0.2011241667, -0.0306029022,  0.1183705475,
        -0.0272104968, -0.0580746013,  0.0405130860,  0.0162250643,
        -0.0299394142,  0.0039896935,  0.0139375423, -0.0077557814,
        -0.0034059318,  0.0049745750, -0.0005058826, -0.0027772683
    };
    float Isum,Qsum;

    /* Convert unsigned to signed */
    for(uint32_t i=0; i<sigLenght; i++)
        sigIn[i] ^= 0x80;  // XOR with a binary mask to flip the first bit (sign)
        //sigIn[i] = (int8_t)((int32_t)samples[i] - 127);

    /* Economic mixer @ fs/4 (upper band)
       At fs/4, sin and cosin calculation are no longueur necessary.

               0   | pi/2 |  pi  | 3pi/2
             ----------------------------
       sin =   0   |  1   |  0   |  -1  |
       cos =   1   |  0   | -1   |   0  |

       out_I = in_I * cos(x) - in_Q * sin(x)
       out_Q = in_Q * cos(x) + in_I * sin(x)
       (Weaver technique, keep the upper band, IQ inverted on RTL devices)
    */
    int8_t tmp;
    for (uint32_t i=0; i<sigLenght; i+=8) {
        tmp = -sigIn[i+3];
        sigIn[i+3] = sigIn[i+2];
        sigIn[i+2] = tmp;

        sigIn[i+4] = -sigIn[i+4];
        sigIn[i+5] = -sigIn[i+5];

        tmp = -sigIn[i+6];
        sigIn[i+6] = sigIn[i+7];
        sigIn[i+7] = tmp;
    }

    /* CIC decimator (N=2)
       (could be not perfect in time for some sampling rate.
       Ex: AirSpy vs AirSpy Mini, but works fine in practice)
       Info: * Understanding CIC Compensation Filters
               https://www.altera.com/en_US/pdfs/literature/an/an455.pdf
             * Understanding cascaded integrator-comb filters
               http://www.embedded.com/design/configurable-systems/4006446/Understanding-cascaded-integrator-comb-filters
    */
    for(int32_t i=0; i<sigLenght/2; i++) {
        /* Integrator stages (N=2) */
        Ix1 += (int32_t)sigIn[i*2];
        Qx1 += (int32_t)sigIn[i*2+1];
        Ix2 += Ix1;
        Qx2 += Qx1;

        /* Decimation R=6400 */
        decimationIndex++;
        if (decimationIndex < DOWNSAMPLING) {
            continue;
        }

        // FIXME/TODO : some optimisition here
        /* 1st Comb */
        Iy1  = Ix2 - It1z;
        It1z = It1y;
        It1y = Ix2;
        Qy1  = Qx2 - Qt1z;
        Qt1z = Qt1y;
        Qt1y = Qx2;

        /* 2nd Comd */
        Iy2  = Iy1 - It2z;
        It2z = It2y;
        It2y = Iy1;
        Qy2  = Qy1 - Qt2z;
        Qt2z = Qt2y;
        Qt2y = Qy1;

        // FIXME/TODO : could be made with int32_t (8 bits, 20 bits)
        /* FIR compensation filter */
        Isum=0.0, Qsum=0.0;
        for (uint32_t j=0; j<32; j++) {
            Isum += firI[j]*zCoef[j];
            Qsum += firQ[j]*zCoef[j];
            if (j<31) {
                firI[j] = firI[j+1];
                firQ[j] = firQ[j+1];
            }
        }
        firI[31] = (float)Iy2;
        firQ[31] = (float)Qy2;
        Isum += firI[31]*zCoef[32];
        Qsum += firQ[31]*zCoef[32];

        /* Save the result in the buffer */
        if (rx_state.iqIndex < (SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE)) {
            /* Lock the buffer during writing */     // Overkill ?!
            pthread_rwlock_wrlock(&dec.rw);
            rx_state.iSamples[rx_state.iqIndex] = Isum;
            rx_state.qSamples[rx_state.iqIndex] = Qsum;
            pthread_rwlock_unlock(&dec.rw);
            rx_state.iqIndex++;
        } else {
            if (rx_state.decode_flag == false) {
                /* Send a signal to the other thread to start the decoding */
                pthread_mutex_lock(&dec.ready_mutex);
                pthread_cond_signal(&dec.ready_cond);
                pthread_mutex_unlock(&dec.ready_mutex);
                rx_state.decode_flag = true;
                //printf("RX done! [Buffer size: %d]\n", rx_state.iqIndex);
            }
        }
        decimationIndex = 0;
    }
}


/* Thread for RX blocking function */
static void *rtlsdr_rx(void *arg) {
    /* Read & blocking call */
    rtlsdr_read_async(rtl_device, rtlsdr_callback, NULL, 0, DEFAULT_BUF_LENGTH);
    exit(0);
    return 0;
}


void postSpots(uint32_t n_results) {
    CURL *curl;
    CURLcode res;
    char url[256]; // FIXME, possible buffer overflow

    for (uint32_t i=0; i<n_results; i++) {
        sprintf(url,"http://wsprnet.org/post?function=wspr&rcall=%s&rgrid=%s&rqrg=%.6f&date=%s&time=%s&sig=%.0f&dt=%.1f&tqrg=%.6f&tcall=%s&tgrid=%s&dbm=%s&version=0.2r_wsprd&mode=2",
                dec_options.rcall, dec_options.rloc, dec_results[i].freq, dec_options.date, dec_options.uttime,
                dec_results[i].snr, dec_results[i].dt, dec_results[i].freq,
                dec_results[i].call, dec_results[i].loc, dec_results[i].pwr);

        printf("Spot : %3.2f %4.2f %10.6f %2d  %-s\n",
               dec_results[i].snr, dec_results[i].dt, dec_results[i].freq,
               (int)dec_results[i].drift, dec_results[i].message);

        curl = curl_easy_init();
        if(curl) {
            curl_easy_setopt(curl, CURLOPT_URL, url);
            curl_easy_setopt(curl, CURLOPT_NOBODY, 1);
            res = curl_easy_perform(curl);

            if(res != CURLE_OK)
                fprintf(stderr, "curl_easy_perform() failed: %s\n",curl_easy_strerror(res));

            curl_easy_cleanup(curl);
        }
    }
    if (n_results == 0)
        printf("No spot\n");
}


static void *wsprDecoder(void *arg) {
    /* WSPR decoder use buffers of 45000 samples (hardcoded)
       (120 sec max @ 375sps = 45000 samples)
    */
    static float iSamples[45000]={0};
    static float qSamples[45000]={0};
    static uint32_t samples_len;
    int32_t n_results=0;

    while (!rx_state.exit_flag) {
        pthread_mutex_lock(&dec.ready_mutex);
        pthread_cond_wait(&dec.ready_cond, &dec.ready_mutex);
        pthread_mutex_unlock(&dec.ready_mutex);

        if(rx_state.exit_flag)  // Abord case, final sig
            break;

        /* Lock the buffer access and make a local copy */
        pthread_rwlock_wrlock(&dec.rw);
        memcpy(iSamples, rx_state.iSamples, rx_state.iqIndex * sizeof(float));
        memcpy(qSamples, rx_state.qSamples, rx_state.iqIndex * sizeof(float));
        samples_len = rx_state.iqIndex;  // Overkill ?
        pthread_rwlock_unlock(&dec.rw);

        /* Date and time will be updated/overload during the search & decoding process
           Make a simple copy
        */
        memcpy(dec_options.date, rx_options.date, sizeof(rx_options.date));
        memcpy(dec_options.uttime, rx_options.uttime, sizeof(rx_options.uttime));

        /* DEBUG -- Save samples
        printf("Writing file\n");
        FILE* fd = NULL;
        fd = fopen("samples.bin", "wb");
        int r=fwrite(rx_state.iSamples, sizeof(float), samples_len, fd);
        printf("%d samples written file\n", r);
        fclose(fd);
        */

        /* Search & decode the signal */
        wspr_decode(iSamples, qSamples, samples_len, dec_options, dec_results, &n_results);
        postSpots(n_results);

    }
    pthread_exit(NULL);
}


double atofs(char *s) {
    /* standard suffixes */
    char last;
    uint32_t len;
    double suff = 1.0;
    len = strlen(s);
    last = s[len-1];
    s[len-1] = '\0';
    switch (last) {
    case 'g':
    case 'G':
        suff *= 1e3;
    case 'm':
    case 'M':
        suff *= 1e3;
    case 'k':
    case 'K':
        suff *= 1e3;
        suff *= atof(s);
        s[len-1] = last;
        return suff;
    }
    s[len-1] = last;
    return atof(s);
}


int32_t parse_u64(char* s, uint64_t* const value) {
    uint_fast8_t base = 10;
    char* s_end;
    uint64_t u64_value;

    if( strlen(s) > 2 ) {
        if( s[0] == '0' ) {
            if( (s[1] == 'x') || (s[1] == 'X') ) {
                base = 16;
                s += 2;
            } else if( (s[1] == 'b') || (s[1] == 'B') ) {
                base = 2;
                s += 2;
            }
        }
    }

    s_end = s;
    u64_value = strtoull(s, &s_end, base);
    if( (s != s_end) && (*s_end == 0) ) {
        *value = u64_value;
        return 1;
    } else {
        return 0;
    }
}


/* Reset flow control variable & decimation variables */
void initSampleStorage() {
    rx_state.decode_flag = false;
    rx_state.iqIndex=0;
}


/* Default options for the decoder */
void initDecoder_options() {
    dec_options.usehashtable = 1;
    dec_options.npasses = 2;
    dec_options.subtraction = 1;
    dec_options.quickmode = 0;
}


/* Default options for the receiver */
void initrx_options() {
    rx_options.gain = 29;
    rx_options.autogain = 0;
    rx_options.ppm = 0;
    rx_options.shift = 0;
}


void sigint_callback_handler(int signum) {
    fprintf(stdout, "Caught signal %d\n", signum);
    rx_state.exit_flag = true;
}


void usage(void) {
    fprintf(stderr,
            "rtlsdr_wsprd, a simple WSPR daemon for RTL receivers\n\n"
            "Use:\trtlsdr_wsprd -f frequency -c callsign -l locator [options]\n"
            "\t-f dial frequency [(,k,M) Hz], check http://wsprnet.org/ for freq.\n"
            "\t-c your callsign (12 chars max)\n"
            "\t-l your locator grid (6 chars max)\n"
            "Receiver extra options:\n"
            "\t-g gain [0-49] (default: 29)\n"
            "\t-a auto gain (default: off)\n"
            "\t-o frequency offset (default: 0)\n"
            "\t-p crystal correction factor (ppm) (default: 0)\n"
            "\t-u upconverter (default: 0, example: 125M)\n"
            "\t-d direct dampling [0,1,2] (default: 0, 1 for I input, 2 for Q input)\n"
            "Decoder extra options:\n"
            "\t-H do not use (or update) the hash table\n"
            "\t-Q quick mode, doesn't dig deep for weak signals\n"
            "\t-S single pass mode, no subtraction (same as original wsprd)\n"
            "Example:\n"
            "\trtlsdr_wsprd -f 144.489M -c A1XYZ -l AB12cd -g 29 -o -4200\n");
    exit(1);
}


int main(int argc, char** argv) {
    uint32_t opt;

    int32_t  rtl_result;
    uint32_t rtl_index = 0; // By default, use the first RTLSDR
    int32_t  rtl_count;
    char     rtl_vendor[256], rtl_product[256], rtl_serial[256];

    initrx_options();
    initDecoder_options();

    /* RX buffer allocation */
    rx_state.iSamples=malloc(sizeof(float)*SIGNAL_LENGHT*SIGNAL_SAMPLE_RATE);
    rx_state.qSamples=malloc(sizeof(float)*SIGNAL_LENGHT*SIGNAL_SAMPLE_RATE);

    /* Stop condition setup */
    rx_state.exit_flag   = false;
    rx_state.decode_flag = false;

    if (argc <= 1)
        usage();

    while ((opt = getopt(argc, argv, "f:c:l:g:a:o:p:u:d:H:Q:S")) != -1) {
        switch (opt) {
        case 'f': // Frequency
            rx_options.dialfreq = (uint32_t)atofs(optarg);
            break;
        case 'c': // Callsign
            sprintf(dec_options.rcall, "%.12s", optarg);
            break;
        case 'l': // Locator / Grid
            sprintf(dec_options.rloc, "%.6s", optarg);
            break;
        case 'g': // Small signal amplifier gain
            rx_options.gain = atoi(optarg);
            if (rx_options.gain < 0) rx_options.gain = 0;
            if (rx_options.gain > 49 ) rx_options.gain = 49;
            rx_options.gain *= 10;
            break;
        case 'a': // Auto gain
            rx_options.autogain = atoi(optarg);
            if (rx_options.autogain < 0) rx_options.autogain = 0;
            if (rx_options.autogain > 1 ) rx_options.autogain = 1;
            break;
        case 'o': // Fine frequency correction
            rx_options.shift = atoi(optarg);
            break;
        case 'p':
            rx_options.ppm = atoi(optarg);
            break;
        case 'u': // Upconverter frequency
            rx_options.upconverter = (uint32_t)atofs(optarg);
            break;
        case 'd': // Direct Sampling
            rx_options.directsampling = (uint32_t)atofs(optarg);
            break;
        case 'H': // Decoder option, use a hastable
            dec_options.usehashtable = 0;
            break;
        case 'Q': // Decoder option, faster
            dec_options.quickmode = 1;
            break;
        case 'S': // Decoder option, single pass mode (same as original wsprd)
            dec_options.subtraction = 0;
            dec_options.npasses = 1;
            break;
        default:
            usage();
            break;
        }
    }

    if (rx_options.dialfreq == 0) {
        fprintf(stderr, "Please specify a dial frequency.\n");
        fprintf(stderr, " --help for usage...\n");
        exit(1);
    }

    if (dec_options.rcall[0] == 0) {
        fprintf(stderr, "Please specify your callsign.\n");
        fprintf(stderr, " --help for usage...\n");
        exit(1);
    }

    if (dec_options.rloc[0] == 0) {
        fprintf(stderr, "Please specify your locator.\n");
        fprintf(stderr, " --help for usage...\n");
        exit(1);
    }

    /* Calcule shift offset */
    rx_options.realfreq = rx_options.dialfreq + rx_options.shift + rx_options.upconverter;

    /* Store the frequency used for the decoder */
    dec_options.freq = rx_options.dialfreq;

    /* If something goes wrong... */
    signal(SIGINT, &sigint_callback_handler);
    signal(SIGTERM, &sigint_callback_handler);
    signal(SIGILL, &sigint_callback_handler);
    signal(SIGFPE, &sigint_callback_handler);
    signal(SIGSEGV, &sigint_callback_handler);
    signal(SIGABRT, &sigint_callback_handler);

    /* Init & parameter the device */
    rtl_count = rtlsdr_get_device_count();
    if (!rtl_count) {
        fprintf(stderr, "No supported devices found\n");
        return EXIT_FAILURE;
    }


    fprintf(stderr, "Found %d device(s):\n", rtl_count);
    for (uint32_t i=0; i<rtl_count; i++) {
        rtlsdr_get_device_usb_strings(i, rtl_vendor, rtl_product, rtl_serial);
        fprintf(stderr, "  %d:  %s, %s, SN: %s\n", i, rtl_vendor, rtl_product, rtl_serial);
    }
    fprintf(stderr, "\nUsing device %d: %s\n", rtl_index, rtlsdr_get_device_name(rtl_index));


    rtl_result = rtlsdr_open(&rtl_device, rtl_index);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to open rtlsdr device #%d.\n", rtl_index);
        return EXIT_FAILURE;
    }

    rtl_result = rtlsdr_set_direct_sampling(rtl_device, rx_options.directsampling);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to set sample rate\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }

    rtl_result = rtlsdr_set_sample_rate(rtl_device, SAMPLING_RATE);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to set sample rate\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }


    rtl_result = rtlsdr_set_tuner_gain_mode(rtl_device, 1);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to enable manual gain\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }


    if (rx_options.autogain) {
        rtl_result = rtlsdr_set_tuner_gain_mode(rtl_device, 0);
        if (rtl_result != 0) {
            fprintf(stderr, "ERROR: Failed to set tuner gain\n");
            rtlsdr_close(rtl_device);
            return EXIT_FAILURE;
        }
    } else {
        rtl_result = rtlsdr_set_tuner_gain(rtl_device, rx_options.gain);
        if (rtl_result != 0) {
            fprintf(stderr, "ERROR: Failed to set tuner gain\n");
            rtlsdr_close(rtl_device);
            return EXIT_FAILURE;
        }
    }


    if (rx_options.ppm != 0) {
        rtl_result = rtlsdr_set_freq_correction(rtl_device, rx_options.ppm);
        if (rtl_result < 0) {
            fprintf(stderr, "ERROR: Failed to set ppm error\n");
            rtlsdr_close(rtl_device);
            return EXIT_FAILURE;
        }
    }


    rtl_result = rtlsdr_set_center_freq(rtl_device, rx_options.realfreq + FS4_RATE + 1500);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to set frequency\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }


    rtl_result = rtlsdr_reset_buffer(rtl_device);
    if (rtl_result < 0) {
        fprintf(stderr, "ERROR: Failed to reset buffers.\n");
        rtlsdr_close(rtl_device);
        return EXIT_FAILURE;
    }

    /* Print used parameter */
    time_t rawtime;
    time ( &rawtime );
    struct tm *gtm = gmtime(&rawtime);
    printf("\nStarting rtlsdr-wsprd (%04d-%02d-%02d, %02d:%02dz) -- Version 0.2\n",
           gtm->tm_year + 1900, gtm->tm_mon + 1, gtm->tm_mday, gtm->tm_hour, gtm->tm_min);
    printf("  Callsign     : %s\n", dec_options.rcall);
    printf("  Locator      : %s\n", dec_options.rloc);
    printf("  Dial freq.   : %d Hz\n", rx_options.dialfreq);
    printf("  Real freq.   : %d Hz\n", rx_options.realfreq);
    printf("  PPM factor   : %d\n", rx_options.ppm);
    if(rx_options.autogain)
        printf("  Auto gain    : enable\n");
    else
        printf("  Gain         : %d dB\n", rx_options.gain/10);


    /* Time alignment stuff */
    struct timeval lTime;
    gettimeofday(&lTime, NULL);
    uint32_t sec   = lTime.tv_sec % 120;
    uint32_t usec  = sec * 1000000 + lTime.tv_usec;
    uint32_t uwait = 120000000 - usec;
    printf("Wait for time sync (start in %d sec)\n\n", uwait/1000000);

    /* Prepare a low priority param for the decoder thread */
    struct sched_param param;
    pthread_attr_init(&dec.tattr);
    pthread_attr_setschedpolicy(&dec.tattr, SCHED_RR);
    pthread_attr_getschedparam(&dec.tattr, &param);
    param.sched_priority = 90;  // = sched_get_priority_min();
    pthread_attr_setschedparam(&dec.tattr, &param);
    //int res=0;
    //printf("get: %d\n", res)

    /* Create a thread and stuff for separate decoding
       Info : https://computing.llnl.gov/tutorials/pthreads/
    */
    pthread_rwlock_init(&dec.rw, NULL);
    pthread_cond_init(&dec.ready_cond, NULL);
    pthread_mutex_init(&dec.ready_mutex, NULL);
    pthread_create(&dongle.thread, NULL, rtlsdr_rx, NULL);
    pthread_create(&dec.thread, &dec.tattr, wsprDecoder, NULL);

    /* Main loop : Wait, read, decode */
    while (!rx_state.exit_flag) {
        /* Wait for time Sync on 2 mins */
        gettimeofday(&lTime, NULL);
        sec   = lTime.tv_sec % 120;
        usec  = sec * 1000000 + lTime.tv_usec;
        uwait = 120000000 - usec + 10000;  // Adding 10ms, to be sure to reach this next minute
        usleep(uwait);
        //printf("SYNC! RX started\n");

        /* Use the Store the date at the begin of the frame */
        time ( &rawtime );
        gtm = gmtime(&rawtime);
        sprintf(rx_options.date,"%02d%02d%02d", gtm->tm_year - 100, gtm->tm_mon + 1, gtm->tm_mday);
        sprintf(rx_options.uttime,"%02d%02d", gtm->tm_hour, gtm->tm_min);

        /* Start to store the samples */
        initSampleStorage();

        while( (rx_state.exit_flag == false) &&
               (rx_state.iqIndex < (SIGNAL_LENGHT * SIGNAL_SAMPLE_RATE) ) ) {
            usleep(250000);
        }
    }

    /* Stop the RX and free the blocking function */
    rtlsdr_cancel_async(rtl_device);

    /* Close the RTL device */
    rtlsdr_close(rtl_device);

    printf("Bye!\n");

    /* Wait the thread join (send a signal before to terminate the job) */
    pthread_mutex_lock(&dec.ready_mutex);
    pthread_cond_signal(&dec.ready_cond);
    pthread_mutex_unlock(&dec.ready_mutex);
    pthread_join(dec.thread, NULL);
    pthread_join(dongle.thread, NULL);

    /* Destroy the lock/cond/thread */
    pthread_rwlock_destroy(&dec.rw);
    pthread_cond_destroy(&dec.ready_cond);
    pthread_mutex_destroy(&dec.ready_mutex);
    pthread_exit(NULL);

    return EXIT_SUCCESS;
}
