#pragma once

char get_locator_character_code(char ch);
char get_callsign_character_code(char ch);
long unsigned int pack_grid4_power(char const *grid4, int power);
long unsigned int pack_call(char const *callsign);
void pack_prefix(char *callsign, int32_t *n, int32_t *m, int32_t *nadd);
void interleave(unsigned char *sym);
int get_wspr_channel_symbols(char *rawmessage, char *hashtab, char *loctab, unsigned char *symbols);
