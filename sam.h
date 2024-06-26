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
#define BL_SAM_FLAG_MAX_DIGITS 32    // What should this really be?

// Keep this for initializing static objects, where we don't want to
// call bl_sam_init() every time.
#define BL_SAM_INIT { "", 0, "", 0, 0, NULL, "", 0, 0, NULL, NULL, 0, 0, 0, 0, 0, 0 }

typedef struct
{
    /* SAM specification fields.  Meet or exceed published ranges. */
    // Query template (DNA fragment):
    // Same name means likely from the same template
    char            qname[BL_SAM_QNAME_MAX_CHARS + 1];
    
    // Bit flags indicating mapping results
    unsigned        flag;
    
    // Sequence to which mapped (e.g. chromosome)
    char            rname[BL_SAM_RNAME_MAX_CHARS + 1];
    int64_t         pos;
    
    // Mapping quality
    unsigned char   mapq;

    // Be sure to update bl_sam_free if changed to dynamic allocation
    // Alignment report (more detailed info than flag)
    char            *cigar;
    
    // Be sure to update bl_sam_free if changed to dynamic allocation
    // Seq and pos of next read in the template
    char            rnext[BL_SAM_RNAME_MAX_CHARS + 1];
    int64_t         pnext;
    
    long            tlen;   // Template length.  FIXME: Max size?
    
    char            *seq;   // This can be large, so malloc() it
    char            *qual;  // PHRED scores, same length as seq if present
    
    /* Additional data */
    size_t          cigar_array_size,
		    cigar_len,
		    seq_array_size,
		    seq_len,
		    qual_array_size,
		    qual_len;
}   bl_sam_t;

typedef unsigned int        sam_field_mask_t;

/* Bit flags to select fields in bl_sam_read() or bl_sam_write() */
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

/* Bit flags from SAM specification */
#define BL_SAM_FLAG_PAIRED          0x001
#define BL_SAM_FLAG_PROPER_PAIR     0x002
#define BL_SAM_FLAG_UNMAP           0x004
#define BL_SAM_FLAG_MUNMAP          0x008
#define BL_SAM_FLAG_REVERSE         0x010
#define BL_SAM_FLAG_MREVERSE        0x020
#define BL_SAM_FLAG_READ1           0x040
#define BL_SAM_FLAG_READ2           0x080
#define BL_SAM_FLAG_SECONDARY       0x100
#define BL_SAM_FLAG_QCFAIL          0x200
#define BL_SAM_FLAG_DUP             0x400
#define BL_SAM_FLAG_SUPPLEMENTARY   0x800

#include "sam-rvs.h"
#include "sam-accessors.h"
#include "sam-mutators.h"

// After bl_sam_t def due to mutual recursion
#ifndef _BIOLIBC_GFF3_H_
#include "gff3.h"
#endif

/* sam.c */
FILE *bl_sam_skip_header(FILE *sam_stream);
int bl_sam_copy_header(FILE *header_stream, FILE *sam_stream);
int bl_sam_read(bl_sam_t *sam_alignment, FILE *sam_stream, sam_field_mask_t field_mask);
void bl_sam_copy(bl_sam_t *dest, bl_sam_t *src);
void bl_sam_free(bl_sam_t *sam_alignment);
void bl_sam_init(bl_sam_t *sam_alignment);
int bl_sam_write(bl_sam_t *sam_alignment, FILE *sam_stream, sam_field_mask_t field_mask);
FILE *bl_sam_fopen(const char *filename, const char *mode, char *samtools_flags);
int bl_sam_fclose(FILE *stream);
int64_t bl_sam_gff3_overlap(bl_sam_t *alignment, bl_gff3_t *feature);
int bl_sam_gff3_cmp(bl_sam_t *alignment, bl_gff3_t *feature);

#endif // _BIOLIBC_SAM_H_
