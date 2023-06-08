#ifndef _BIOLIBC_FASTQ_H_
#define _BIOLIBC_FASTQ_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

typedef struct
{
    char    *desc,
	    *seq,
	    *plus,
	    *qual;
    size_t  desc_array_size,
	    seq_array_size,
	    plus_array_size,
	    qual_array_size,
	    desc_len,
	    seq_len,
	    plus_len,
	    qual_len;
}   bl_fastq_t;

#define BL_FASTQ_INIT           { NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0, 0 }

// Must be same as BL_FASTA_LINE_UNLIMITED
#define BL_FASTQ_LINE_UNLIMITED 0

#include "fastq-rvs.h"
#include "fastq-accessors.h"
#include "fastq-mutators.h"

/* fastq.c */
int bl_fastq_read(bl_fastq_t *record, FILE *fastq_stream);
int bl_fastq_write(bl_fastq_t *record, FILE *fastq_stream, size_t max_line_len);
void bl_fastq_free(bl_fastq_t *record);
void bl_fastq_init(bl_fastq_t *record);
size_t bl_fastq_find_adapter_smart(const bl_fastq_t *read, const char *adapter, size_t min_match, unsigned max_mismatch_percent);
size_t bl_fastq_find_adapter_exact(const bl_fastq_t *read, const char *adapter, size_t min_match, unsigned max_mismatch_percent);
size_t bl_fastq_3p_trim(bl_fastq_t *read, size_t new_len);
size_t bl_fastq_find_3p_low_qual(const bl_fastq_t *read, unsigned min_qual, unsigned phred_base);
size_t bl_fastq_name_cmp(bl_fastq_t *read1, bl_fastq_t *read2);

#ifdef __cplusplus
}
#endif

#endif  // _BIOLIBC_FASTQ_H_
