#ifndef _BIOLIBC_OVERLAP_H_
#define _BIOLIBC_OVERLAP_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _INTTYPES_H_
#include <inttypes.h>
#endif

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

// 1-based, inclusive at both ends
typedef struct
{
    uint64_t    feature1_len;
    uint64_t    feature2_len;
    uint64_t    overlap_start;
    uint64_t    overlap_end;
    uint64_t    overlap_len;
}   bl_overlap_t;

#include "overlap-rvs.h"
#include "overlap-accessors.h"
#include "overlap-mutators.h"

/* overlap.c */
int bl_overlap_set_all(bl_overlap_t *overlap, uint64_t feature1_len, uint64_t feature2_len, uint64_t overlap_start, uint64_t overlap_end);
int bl_overlap_print(bl_overlap_t *overlap, FILE *stream, char *feature1_name, char *feature2_name);

#endif // _BIOLIBC_OVERLAP_H_
