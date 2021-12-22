/*
 This file is part of program wsprd, a detector/demodulator/decoder
 for the Weak Signal Propagation Reporter (WSPR) mode.

 File name: wsprd_utils.h

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

#pragma once

void unpack50(signed char *dat, int32_t *n1, int32_t *n2);
int unpackcall(int32_t ncall, char *call);
int unpackgrid(int32_t ngrid, char *grid);
int unpackpfx(int32_t nprefix, char *call);
void deinterleave(unsigned char *sym);

// used by qsort
int doublecomp(const void *elem1, const void *elem2);
int floatcomp(const void *elem1, const void *elem2);

int unpk_(signed char *message, char *hashtab, char *loctab, char *call_loc_pow, char *call, char *loc, char *pwr, char *callsign);
