#ifndef _sam_h_
#define _sam_h_

#ifndef _dsv_h_
#include "dsv.h"
#endif

#define SAM_MAPQ_MAX_CHARS  3
#define SAM_QNAME_MAX_CHARS 4096
#define SAM_RNAME_MAX_CHARS 4096
#define SAM_FLAG_MAX_DIGITS 4096    // What should this really be?
#define SAM_CIGAR_MAX_CHARS 4096
// Usually < 200 for Illumina data, but a few oddballs in SRA CRAMs
#define SAM_SEQ_MAX_CHARS   1024*1024

#define SAM_QNAME(s)        ((s)->qname)
#define SAM_RNAME(s)        ((s)->rname)
#define SAM_POS(s)          ((s)->pos)
#define SAM_MAPQ(s)         ((s)->mapq)
#define SAM_SEQ(s)          ((s)->seq)
#define SAM_QUAL(s)         ((s)->qual)
#define SAM_SEQ_LEN(s)      ((s)->seq_len)
#define SAM_QUAL_LEN(s)     ((s)->qual_len)

// Use this or the function for every new object
#define SAM_ALIGNMENT_INIT  { "", 0, "", 0, 0, "", "", 0, 0, NULL, NULL, 0 }

typedef unsigned int    sam_field_mask_t;

#define SAM_FIELD_ALL   0xfff
#define SAM_FIELD_QNAME 0x001
#define SAM_FIELD_FLAG  0x002
#define SAM_FIELD_RNAME 0x004
#define SAM_FIELD_POS   0x008
#define SAM_FIELD_MAPQ  0x010
#define SAM_FIELD_CIGAR 0x020
#define SAM_FIELD_RNEXT 0x040
#define SAM_FIELD_PNEXT 0x080
#define SAM_FIELD_TLEN  0x100
#define SAM_FIELD_SEQ   0x200
#define SAM_FIELD_QUAL  0x400

typedef struct
{
    /* SAM fields */
    char            qname[SAM_QNAME_MAX_CHARS + 1];
    unsigned        flag;
    char            rname[SAM_RNAME_MAX_CHARS + 1];
    size_t          pos;
    unsigned char   mapq;
    char            cigar[SAM_CIGAR_MAX_CHARS + 1];
    char            rnext[SAM_RNAME_MAX_CHARS + 1];
    size_t          pnext;
    size_t          tlen;   // Max size?
    char            *seq;   // This can be large, so malloc() it
    char            *qual;  // PHRED scores, same length as seq if present
    
    /* Additional data */
    size_t  seq_len;        // Qual len should be the same
    size_t  qual_len;
}   bl_sam_t;

/* sam.c */
int     sam_read_alignment(FILE *sam_stream, bl_sam_t *sam_alignment, sam_field_mask_t field_mask);
void    sam_copy_alignment(bl_sam_t *dest, bl_sam_t *src);
void    sam_free_alignment(bl_sam_t *sam_alignment);
void    sam_init_alignment(bl_sam_t *sam_alignment, size_t seq_len, sam_field_mask_t field_mask);

#endif // _sam_h_
