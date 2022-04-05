#ifndef _BIOLIBC_SAM_H_
#define _BIOLIBC_SAM_H_

#ifndef _INTTYPES_H
#include <inttypes.h>
#endif

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

#define BL_SAM_MAPQ_MAX_CHARS  5
#define BL_SAM_QNAME_MAX_CHARS 4096
#define BL_SAM_RNAME_MAX_CHARS 4096
#define BL_SAM_FLAG_MAX_DIGITS 4096    // What should this really be?
#define BL_SAM_CIGAR_MAX_CHARS 4096

// Keep this for initializing static objects, where we don't want to
// call bl_sam_init() every time.
#define BL_SAM_INIT { "", 0, "", 0, 0, "", "", 0, 0, NULL, NULL, 0, 0, 0, 0 }

typedef struct
{
    /* SAM specification fields.  Meet or exceed published ranges. */
    char            qname[BL_SAM_QNAME_MAX_CHARS + 1];
    unsigned        flag;
    char            rname[BL_SAM_RNAME_MAX_CHARS + 1];
    int64_t         pos;
    unsigned char   mapq;
    char            cigar[BL_SAM_CIGAR_MAX_CHARS + 1];
    char            rnext[BL_SAM_RNAME_MAX_CHARS + 1];
    int64_t         pnext;
    long            tlen;   // Max size?
    char            *seq;   // This can be large, so malloc() it
    char            *qual;  // PHRED scores, same length as seq if present
    
    /* Additional data */
    size_t          seq_array_size,
		    seq_len,
		    qual_array_size,
		    qual_len;
}   bl_sam_t;

typedef unsigned int        sam_field_mask_t;

#define BL_SAM_FIELD_ALL    0xfff
#define BL_SAM_FIELD_QNAME  0x001
#define BL_SAM_FIELD_FLAG   0x002
#define BL_SAM_FIELD_RNAME  0x004
#define BL_SAM_FIELD_POS    0x008
#define BL_SAM_FIELD_MAPQ   0x010
#define BL_SAM_FIELD_CIGAR  0x020
#define BL_SAM_FIELD_RNEXT  0x040
#define BL_SAM_FIELD_PNEXT  0x080
#define BL_SAM_FIELD_TLEN   0x100
#define BL_SAM_FIELD_SEQ    0x200
#define BL_SAM_FIELD_QUAL   0x400

#include "sam-rvs.h"
#include "sam-accessors.h"
#include "sam-mutators.h"

/* sam.c */
int bl_sam_read(bl_sam_t *sam_alignment, FILE *sam_stream, sam_field_mask_t field_mask);
void bl_sam_copy(bl_sam_t *dest, bl_sam_t *src);
void bl_sam_free(bl_sam_t *sam_alignment);
void bl_sam_init(bl_sam_t *sam_alignment);
int bl_sam_write(bl_sam_t *sam_alignment, FILE *sam_stream, sam_field_mask_t field_mask);
FILE *bl_sam_fopen(const char *filename, const char *mode);
int bl_sam_fclose(FILE *stream);

#endif // _BIOLIBC_SAM_H_
