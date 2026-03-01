/*
 This file is part of program wsprd, a detector/demodulator/decoder
 for the Weak Signal Propagation Reporter (WSPR) mode.

 File name: wsprd_utils.c

 Copyright 2001-2015, Joe Taylor, K1JT

 Most of the code is based on work by Steven Franke, K9AN, which
 in turn was based on earlier work by K1JT.

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
#include <ctype.h>
#include <string.h>
#include <stdint.h>

#include "./nhash.h"
#include "./wsprd.h"
#include "./wsprd_utils.h"


void unpack50(signed char *dat, int32_t *n1, int32_t *n2) {
    int32_t i, i4;

    i = dat[0];
    i4 = i & 255;
    *n1 = i4 << 20;

    i = dat[1];
    i4 = i & 255;
    *n1 = *n1 + (i4 << 12);

    i = dat[2];
    i4 = i & 255;
    *n1 = *n1 + (i4 << 4);

    i = dat[3];
    i4 = i & 255;
    *n1 = *n1 + ((i4 >> 4) & 15);
    *n2 = (i4 & 15) << 18;

    i = dat[4];
    i4 = i & 255;
    *n2 = *n2 + (i4 << 10);

    i = dat[5];
    i4 = i & 255;
    *n2 = *n2 + (i4 << 2);

    i = dat[6];
    i4 = i & 255;
    *n2 = *n2 + ((i4 >> 6) & 3);
}

int unpackcall(int32_t ncall, char *call) {
    char c[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E',
                'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                'U', 'V', 'W', 'X', 'Y', 'Z', ' '};
    int32_t n;
    int i;
    char tmp[7];

    n = ncall;
    snprintf(call, 13, "......");
    if (n < 262177560) {
        i = n % 27 + 10;
        tmp[5] = c[i];
        n = n / 27;
        i = n % 27 + 10;
        tmp[4] = c[i];
        n = n / 27;
        i = n % 27 + 10;
        tmp[3] = c[i];
        n = n / 27;
        i = n % 10;
        tmp[2] = c[i];
        n = n / 10;
        i = n % 36;
        tmp[1] = c[i];
        n = n / 36;
        i = n;
        tmp[0] = c[i];
        tmp[6] = '\0';
        // remove leading whitespace
        for (i = 0; i < 5; i++) {
            if (tmp[i] != c[36])
                break;
        }
        snprintf(call, 13, "%-6s", &tmp[i]);
        // remove trailing whitespace
        for (i = 0; i < 6; i++) {
            if (call[i] == c[36]) {
                call[i] = '\0';
            }
        }
    } else {
        return 0;
    }
    return 1;
}

int unpackgrid(int32_t ngrid, char *grid) {
    char c[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E',
                'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                'U', 'V', 'W', 'X', 'Y', 'Z', ' '};
    int dlat, dlong;

    ngrid = ngrid >> 7;
    if (ngrid < 32400) {
        dlat = (ngrid % 180) - 90;
        dlong = (ngrid / 180) * 2 - 180 + 2;
        if (dlong < -180)
            dlong = dlong + 360;
        if (dlong > 180)
            dlong = dlong + 360;
        int nlong = 60.0 * (180.0 - dlong) / 5.0;
        int n1 = nlong / 240;
        int n2 = (nlong - 240 * n1) / 24;
        grid[0] = c[10 + n1];
        grid[2] = c[n2];

        int nlat = 60.0 * (dlat + 90) / 2.5;
        n1 = nlat / 240;
        n2 = (nlat - 240 * n1) / 24;
        grid[1] = c[10 + n1];
        grid[3] = c[n2];
    } else {
        snprintf(grid, 5, "XXXX");
        return 0;
    }
    return 1;
}

int unpackpfx(int32_t nprefix, char *call) {
    char nc, pfx[4] = {'\0'}, tmpcall[13];
    int i;
    int32_t n;

    snprintf(tmpcall, sizeof(tmpcall), "%s", call);
    if (nprefix < 60000) {
        // add a prefix of 1 to 3 characters
        n = nprefix;
        for (i = 2; i >= 0; i--) {
            nc = n % 37;
            if ((nc >= 0) && (nc <= 9)) {
                pfx[i] = nc + 48;
            } else if ((nc >= 10) && (nc <= 35)) {
                pfx[i] = nc + 55;
            } else {
                pfx[i] = ' ';
            }
            n = n / 37;
        }

        char *p = strrchr(pfx, ' ');
        snprintf(call, 13, "%s/%s", p ? p + 1 : pfx, tmpcall);

    } else {
        // add a suffix of 1 or 2 characters
        nc = nprefix - 60000;
        if ((nc >= 0) && (nc <= 9)) {
            pfx[0] = nc + 48;
            snprintf(call, 13, "%s/%c", tmpcall, pfx[0]);
        } else if ((nc >= 10) && (nc <= 35)) {
            pfx[0] = nc + 55;
            snprintf(call, 13, "%s/%c", tmpcall, pfx[0]);
        } else if ((nc >= 36) && (nc <= 125)) {
            pfx[0] = (nc - 26) / 10 + 48;
            pfx[1] = (nc - 26) % 10 + 48;
            snprintf(call, 13, "%s/%c%c", tmpcall, pfx[0], pfx[1]);
        } else {
            return 0;
        }
    }
    return 1;
}

void deinterleave(unsigned char *sym) {
    unsigned char tmp[162];
    unsigned char p, i, j;

    p = 0;
    i = 0;
    while (p < 162) {
        j = ((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
        if (j < 162) {
            tmp[p] = sym[j];
            p = p + 1;
        }
        i = i + 1;
    }
    for (i = 0; i < 162; i++) {
        sym[i] = tmp[i];
    }
}

// used by qsort
int doublecomp(const void *elem1, const void *elem2) {
    if (*(const double *)elem1 < *(const double *)elem2)
        return -1;
    return *(const double *)elem1 > *(const double *)elem2;
}

int floatcomp(const void *elem1, const void *elem2) {
    if (*(const float *)elem1 < *(const float *)elem2)
        return -1;
    return *(const float *)elem1 > *(const float *)elem2;
}

int unpk_(signed char *message, char *hashtab, char *loctab, char *call_loc_pow, char *call, char *loc, char *pwr, char *callsign) {
    int n1, n2, n3, ndbm, ihash, nadd, noprint = 0;
    char grid[5], grid6[7], cdbm[4];

    unpack50(message, &n1, &n2);
    if (!unpackcall(n1, callsign))
        return 1;
    if (!unpackgrid(n2, grid))
        return 1;
    int ntype = (n2 & 127) - 64;
    callsign[12] = 0;
    grid[4] = 0;

    /*
     Based on the value of ntype, decide whether this is a Type 1, 2, or
     3 message.
     * Type 1: 6 digit call, grid, power - ntype is positive and is a member
     of the set {0,3,7,10,13,17,20...60}
     * Type 2: extended callsign, power - ntype is positive but not
     a member of the set of allowed powers
     * Type 3: hash, 6 digit grid, power - ntype is negative.
     */

    if ((ntype >= 0) && (ntype <= 62)) {
        int nu = ntype % 10;
        if (nu == 0 || nu == 3 || nu == 7) {
            ndbm = ntype;
            snprintf(cdbm, sizeof(cdbm), "%02d", ndbm);
            snprintf(call_loc_pow, 23, "%s %s %s", callsign, grid, cdbm);
            ihash = nhash(callsign, strlen(callsign), (uint32_t)146);
            snprintf(hashtab + ihash * HASHTAB_ENTRY_LEN, HASHTAB_ENTRY_LEN, "%s", callsign);
            snprintf(loctab + ihash * LOCTAB_ENTRY_LEN, LOCTAB_ENTRY_LEN, "%s", grid);

            snprintf(call, HASHTAB_ENTRY_LEN, "%s", callsign);
            snprintf(loc, 7, "%s", grid);
            snprintf(pwr, 3, "%s", cdbm);
        } else {
            nadd = nu;
            if (nu > 3) nadd = nu - 3;
            if (nu > 7) nadd = nu - 7;
            n3 = n2 / 128 + HASHTAB_SIZE * (nadd - 1);
            if (!unpackpfx(n3, callsign)) return 1;
            ndbm = ntype - nadd;
            snprintf(cdbm, sizeof(cdbm), "%2d", ndbm);
            snprintf(call_loc_pow, 23, "%s %s", callsign, cdbm);
            int nu = ndbm % 10;
            if (nu == 0 || nu == 3 || nu == 7) {
                ihash = nhash(callsign, strlen(callsign), (uint32_t)146);
                snprintf(hashtab + ihash * HASHTAB_ENTRY_LEN, HASHTAB_ENTRY_LEN, "%s", callsign);
            } else
                noprint = 1;
        }
    } else if (ntype < 0) {
        ndbm = -(ntype + 1);
        memset(grid6, 0, sizeof(char) * 7);
        size_t len = 6;
        snprintf(grid6, sizeof(grid6), "%c%.*s", callsign[len - 1], (int)(len - 1), callsign);
        int nu = ndbm % 10;
        if ((nu != 0 && nu != 3 && nu != 7) ||
            !isalpha(grid6[0]) || !isalpha(grid6[1]) ||
            !isdigit(grid6[2]) || !isdigit(grid6[3])) {
            // not testing 4'th and 5'th chars because of this case: <PA0SKT/2> JO33 40
            // grid is only 4 chars even though this is a hashed callsign...
            noprint = 1;
        }

        ihash = (n2 - ntype - 64) / 128;
        if (strncmp(hashtab + ihash * HASHTAB_ENTRY_LEN, "\0", 1) != 0) {
            snprintf(callsign, HASHTAB_ENTRY_LEN, "<%s>", hashtab + ihash * HASHTAB_ENTRY_LEN);
        } else {
            snprintf(callsign, HASHTAB_ENTRY_LEN, "<...>");
        }

        snprintf(cdbm, sizeof(cdbm), "%2d", ndbm);
        snprintf(call_loc_pow, 23, "%s %s %s", callsign, grid6, cdbm);

        snprintf(call, HASHTAB_ENTRY_LEN, "%s", callsign);
        snprintf(loc, 7, "%s", grid6);
        snprintf(pwr, 3, "%s", cdbm);

        // I don't know what to do with these... They show up as "A000AA" grids.
        if (ntype == -64) 
            noprint = 1;
    }
    return noprint;
}
