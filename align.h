#ifndef _BIOLIBC_ALIGN_H_
#define _BIOLIBC_ALIGN_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

typedef struct
{
    size_t min_match;
    unsigned max_mismatch_percent;
}   bl_align_t;

#include "align-rvs.h"
#include "align-accessors.h"
#include "align-mutators.h"

/* align.c */
size_t bl_align_map_seq_sub(const bl_align_t *params, const char *big, size_t big_len, const char *little, size_t little_len);
size_t bl_align_map_seq_exact(const bl_align_t *params, const char *big, size_t big_len, const char *little, size_t little_len);

#ifdef __cplusplus
}
#endif

#endif  // _BIOLIBC_ALIGN_H_
