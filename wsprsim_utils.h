#ifndef WSPRSIM_UTILS_H
#define WSPRSIM_UTILS_H

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

extern int printdata;

char get_locator_character_code(char ch);

char get_callsign_character_code(char ch);

long unsigned int pack_grid4_power(char *grid4, int power);

long unsigned int pack_call(char *callsign);

void pack_prefix(char *callsign, int32_t *n, int32_t *m, int32_t *nadd );

void interleave(unsigned char *sym);

int get_wspr_channel_symbols(char* message, char* hashtab, unsigned char* symbols);

#endif
