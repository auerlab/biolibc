
#ifndef _BIOLIBC_FASTA_H_
#include "fasta.h"
#endif

#ifndef _BIOLIBC_FASTQ_H_
#include "fastq.h"
#endif

typedef struct
{
    int     format;
    union
    {
	bl_fasta_t  fasta;
	bl_fastq_t  fastq;
    };
}   bl_fastx_t;

#define BL_FASTX_FORMAT_UNKNOWN 0
#define BL_FASTX_FORMAT_FASTA   1
#define BL_FASTX_FORMAT_FASTQ   2

#define BL_FASTX_INIT           { BL_FASTX_FORMAT_UNKNOWN }

#define BL_FASTX_FORMAT(fx) ((fx)->format)

#define BL_FASTX_LINE_UNLIMITED BL_FASTA_LINE_UNLIMITED

/* fastx.c */
int bl_fastx_read(bl_fastx_t *record, FILE *fastx_stream);
int bl_fastx_write(bl_fastx_t *record, FILE *fastx_stream, size_t max_line_len);
void bl_fastx_free(bl_fastx_t *record);
void bl_fastx_init(bl_fastx_t *record, FILE *fastx_stream);
char *bl_fastx_desc(bl_fastx_t *record);
size_t bl_fastx_desc_len(bl_fastx_t *record);
char *bl_fastx_seq(bl_fastx_t *record);
size_t bl_fastx_seq_len(bl_fastx_t *record);
char *bl_fastx_plus(bl_fastx_t *record);
size_t bl_fastx_plus_len(bl_fastx_t *record);
char *bl_fastx_qual(bl_fastx_t *record);
size_t bl_fastx_qual_len(bl_fastx_t *record);
