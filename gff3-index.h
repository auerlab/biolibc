#ifndef _BIOLIBC_GFF3_INDEX_H_
#define _BIOLIBC_GFF3_INDEX_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _INTTYPES_H_
#include <inttypes.h>
#endif

#ifndef _BIOLIBC_GFF3_H_
#include "gff3.h"
#endif

#define BL_GFF3_INDEX_INIT { 0, 0, NULL, NULL }

#define BL_GFF3_INDEX_OK             0
#define BL_GFF3_INDEX_MALLOC_FAILED  -1
#define BL_GFF3_INDEX_BAD_ARG        -2

typedef struct
{
    size_t      array_size;
    size_t      count;
    long        *file_pos;  // Return type of ftell()
    char        **seqid;
    int64_t     *start;
    int64_t     *end;
}   bl_gff3_index_t;

#include "gff3-index-rvs.h"
#include "gff3-index-accessors.h"
#include "gff3-index-mutators.h"

int bl_gff3_index_add(bl_gff3_index_t *gi, bl_gff3_t *feature);
int bl_gff3_index_seek_reverse(bl_gff3_index_t *gi, FILE *stream, bl_gff3_t *feature, int64_t gene_count, int64_t max_nt);

#ifdef __cplusplus
}
#endif

#endif  // _BIOLIBC_GFF3_INDEX_H_
