#pragma once

void unpack50( signed char *dat, int32_t *n1, int32_t *n2 );
int unpackcall( int32_t ncall, char *call );
int unpackgrid( int32_t ngrid, char *grid);
int unpackpfx( int32_t nprefix, char *call);
void deinterleave(unsigned char *sym);

// used by qsort
int doublecomp(const void* elem1, const void* elem2);
int floatcomp(const void* elem1, const void* elem2);

int unpk_( signed char *message, char *hashtab, char *call_loc_pow, char *call, char *loc, char *pwr, char *callsign);
