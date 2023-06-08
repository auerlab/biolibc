#ifndef _BIOLIBC_FASTA_H_
#define _BIOLIBC_FASTA_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

typedef struct
{
    char    *desc;
    char    *seq;
    size_t  desc_array_size,
	    seq_array_size,
	    desc_len,
	    seq_len;
}   bl_fasta_t;

#define BL_FASTA_INIT           { NULL, NULL, 0, 0, 0, 0 }
#define BL_FASTA_LINE_UNLIMITED 0

#include "fasta-rvs.h"
#include "fasta-accessors.h"
#include "fasta-mutators.h"

/* fasta.c */
int bl_fasta_read(bl_fasta_t *record, FILE *fasta_stream);
int bl_fasta_write(bl_fasta_t *record, FILE *fasta_stream, size_t chars_per_line);
void bl_fasta_free(bl_fasta_t *record);
void bl_fasta_init(bl_fasta_t *record);

#ifdef __cplusplus
}
#endif

#endif // _BIOLIBC_FASTA_H_
