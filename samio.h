#ifndef __samio_h__
#define __samio_h__

#ifndef __tsvio_h__
#include "tsvio.h"
#endif

#ifndef __biostring_h__
#include "biostring.h"
#endif

#define SAM_MAPQ_MAX_CHARS  3
#define SAM_QNAME_MAX_CHARS 4096
#define SAM_FLAG_MAX_CHARS  4096
#define SAM_RNAME_MAX_CHARS 4096
#define SAM_POS_MAX_DIGITS  4096
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

#define SAM_READ_OK                 0
#define SAM_READ_EOF                -1
#define SAM_READ_OVERFLOW           -2
#define SAM_READ_TRUNCATED          -3

// Use this or the function for every new object
#define SAM_ALIGNMENT_INIT  { "", 0, "", 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0 }

typedef struct
{
    /* SAM fields */
    char            qname[SAM_QNAME_MAX_CHARS + 1];
    unsigned        flag;
    char            rname[SAM_RNAME_MAX_CHARS + 1];
    size_t          pos;
    unsigned char   mapq;
    char            *cigar; // FIXME: Should cigar and rnext be static size?
    char            *rnext; // [SAM_RNAME_MAX_CHARS + 1];
    size_t          pnext;
    size_t          tlen;   // Max size?
    char            *seq;   // This can be large, so malloc() it
    char            *qual;  // PHRED scores, same length as seq if present
    
    /* Additional data */
    size_t  seq_len;        // Qual len should be the same
    size_t  qual_len;
}   sam_alignment_t;

/* samio.c */
int     sam_alignment_read(FILE *sam_stream, sam_alignment_t *sam_alignment);
void    sam_alignment_copy(sam_alignment_t *dest, sam_alignment_t *src);
void    sam_alignment_free(sam_alignment_t *sam_alignment);
void    sam_alignment_init(sam_alignment_t *sam_alignment, size_t seq_len);

#endif // __samio_h__
