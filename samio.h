#ifndef __samio_h__
#define __samio_h__

#ifndef __tsvio_h__
#include "tsvio.h"
#endif

#define SAM_QNAME_MAX_CHARS 4096
#define SAM_FLAG_MAX_CHARS  4096
#define SAM_RNAME_MAX_CHARS 4096
#define SAM_POS_MAX_DIGITS  4096
#define SAM_CIGAR_MAX_CHARS 4096
#define SAM_SEQ_MAX_CHARS   1024*1024

#define SAM_QNAME(s)        ((s)->qname)
#define SAM_RNAME(s)        ((s)->rname)
#define SAM_SEQ(s)          ((s)->seq)
#define SAM_POS(s)          ((s)->pos)
#define SAM_SEQ_LEN(s)      ((s)->seq_len)

#define SAM_READ_OK                 0
#define SAM_READ_EOF                -1
#define SAM_READ_OVERFLOW           -2
#define SAM_READ_TRUNCATED          -3

typedef struct
{
    char    qname[SAM_QNAME_MAX_CHARS + 1],
	    rname[SAM_RNAME_MAX_CHARS + 1],
	    *seq;
    size_t  pos;
    size_t  seq_len;
}   sam_alignment_t;

/* samio.c */
int     sam_alignment_read(FILE *sam_stream, sam_alignment_t *sam_alignment);
void    sam_alignment_copy(sam_alignment_t *dest, sam_alignment_t *src);
void    sam_alignment_free(sam_alignment_t *sam_alignment);

#endif // __samio_h__
