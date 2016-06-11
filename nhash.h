#ifndef NHASH_H_
#define NHASH_H_

#ifdef Win32
#include "win_stdint.h"	/* defines uint32_t etc */
#else
#include <stddef.h>
#include <stdint.h>	/* defines uint32_t etc */
#endif

uint32_t nhash( const void * key, size_t length, uint32_t initval);
uint32_t nhash_( void const * key, size_t const * length, uint32_t const * initval);

#endif
