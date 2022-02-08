#ifndef _BIOLIBC_GFF_INDEX_H_
#define _BIOLIBC_GFF_INDEX_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _INTTYPES_H_
#include <inttypes.h>
#endif

#ifndef _BIOLIBC_GFF_H_
#include "gff.h"
#endif

#define BL_GFF_INDEX_INIT { 0, 0, NULL, NULL }

#define BL_GFF_INDEX_OK             0
#define BL_GFF_INDEX_MALLOC_FAILED  -1
#define BL_GFF_INDEX_BAD_ARG        -2

typedef struct
{
    size_t      array_size;
    size_t      count;
    long        *file_pos;  // Return type of ftell()
    char        **seqid;
    uint64_t    *start;
    uint64_t    *end;
}   bl_gff_index_t;

#include "gff-index-rvs.h"
#include "gff-index-accessors.h"
#include "gff-index-mutators.h"

int bl_gff_index_add(bl_gff_index_t *gi, bl_gff_t *feature);
int bl_gff_index_seek_reverse(bl_gff_index_t *gi, FILE *stream, bl_gff_t *feature, uint64_t gene_count, uint64_t max_nt);
uint64_t str2u64(const char *str);

#endif  // _BIOLIBC_GFF_INDEX_H_
