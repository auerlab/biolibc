#ifndef __samio_h__
#define __samio_h__

#ifndef __tsvio_h__
#include "tsvio.h"
#endif

#define SAM_QNAME_MAX       4096
#define SAM_FLAG_MAX        4096
#define SAM_RNAME_MAX       4096
#define SAM_POS_MAX_DIGITS  4096
#define SAM_CIGAR_MAX       4096
#define SAM_SEQ_MAX         1024*1024

#define SAM_QNAME(s)        ((s)->qname)
#define SAM_RNAME(s)        ((s)->rname)
#define SAM_SEQ(s)          ((s)->seq)
#define SAM_POS(s)          ((s)->pos)
#define SAM_SEQ_LEN(s)      ((s)->seq_len)

typedef struct
{
    char    qname[SAM_QNAME_MAX + 1],
	    rname[SAM_RNAME_MAX + 1],
	    seq[SAM_SEQ_MAX + 1];
    size_t  pos;
    size_t  seq_len;
}   sam_alignment_t;

/* samio.c */
int sam_read_alignment(const char *argv[], FILE *sam_stream, sam_alignment_t *sam_alignment);

#endif // __samio_h__
