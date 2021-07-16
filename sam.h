#ifndef _sam_h_
#define _sam_h_

#ifndef __xtend_h__
#include <xtend.h>
#endif

#define BL_SAM_MAPQ_MAX_CHARS  3
#define BL_SAM_QNAME_MAX_CHARS 4096
#define BL_SAM_RNAME_MAX_CHARS 4096
#define BL_SAM_FLAG_MAX_DIGITS 4096    // What should this really be?
#define BL_SAM_CIGAR_MAX_CHARS 4096
// Usually < 200 for Illumina data, but a few oddballs in SRA CRAMs
#define BL_SAM_SEQ_MAX_CHARS   1024*1024

#define BL_SAM_QNAME(ptr)   ((ptr)->qname)
#define BL_SAM_FLAG(ptr)    ((ptr)->flag)
#define BL_SAM_RNAME(ptr)   ((ptr)->rname)
#define BL_SAM_POS(ptr)     ((ptr)->pos)
#define BL_SAM_MAPQ(ptr)    ((ptr)->mapq)
#define BL_SAM_CIGAR(ptr)   ((ptr)->cigar)
#define BL_SAM_RNEXT(ptr)   ((ptr)->rnext)
#define BL_SAM_PNEXT(ptr)   ((ptr)->pnext)
#define BL_SAM_TLEN(ptr)    ((ptr)->tlen)
#define BL_SAM_SEQ(ptr)     ((ptr)->seq)
#define BL_SAM_QUAL(ptr)    ((ptr)->qual)
#define BL_SAM_SEQ_LEN(ptr) ((ptr)->seq_len)
#define BL_SAM_QUAL_LEN(ptr) ((ptr)->qual_len)

#define BL_SAM_SET_QNAME(ptr,qname) strlcpy((ptr)->qname,qname,BL_SAM_QNAME_MAX_CHARS+1)
#define BL_SAM_SET_FLAG(ptr,flag)   ((ptr)->flag = (flag))
#define BL_SAM_SET_RNAME(ptr,rname) strlcpy((ptr)->rname,rname,BL_SAM_RNAME_MAX_CHARS+1)
#define BL_SAM_SET_POS(ptr,pos)     ((ptr)->pos = (pos))
#define BL_SAM_SET_MAPQ(ptr,mapq)   ((ptr)->mapq = (mapq))
#define BL_SAM_SET_CIGAR(ptr,cigar) strlcpy((ptr)->cigar,cigar,BL_SAM_CIGAR_MAX_CHARS+1)
#define BL_SAM_SET_RNEXT(ptr,rnext) strlcpy((ptr)->rnext,rnext,BL_SAM_RNAME_MAX_CHARS+1)
#define BL_SAM_SET_PNEXT(ptr,pnext) ((ptr)->pnext = (pnext))
#define BL_SAM_SET_TLEN(ptr,tlen)   ((ptr)->tlen = (tlen))
#define BL_SAM_SET_SEQ(ptr,seq)     ((ptr)->seq = (seq))
#define BL_SAM_SET_QUAL(ptr,qual)   ((ptr)->qual = (qual))
#define BL_SAM_SET_SEQ_LEN(ptr,seq_len)     ((ptr)->seq_len = (seq_len))
#define BL_SAM_SET_QUAL_LEN(ptr,qual_len)   ((ptr)->qual_len = (qual_len))

// Use this or the function for every new object
#define BL_SAM_ALIGNMENT_INIT  { "", 0, "", 0, 0, "", "", 0, 0, NULL, NULL, 0 }

typedef unsigned int    sam_field_mask_t;

#define BL_SAM_FIELD_ALL   0xfff
#define BL_SAM_FIELD_QNAME 0x001
#define BL_SAM_FIELD_FLAG  0x002
#define BL_SAM_FIELD_RNAME 0x004
#define BL_SAM_FIELD_POS   0x008
#define BL_SAM_FIELD_MAPQ  0x010
#define BL_SAM_FIELD_CIGAR 0x020
#define BL_SAM_FIELD_RNEXT 0x040
#define BL_SAM_FIELD_PNEXT 0x080
#define BL_SAM_FIELD_TLEN  0x100
#define BL_SAM_FIELD_SEQ   0x200
#define BL_SAM_FIELD_QUAL  0x400

typedef struct
{
    /* SAM fields */
    char            qname[BL_SAM_QNAME_MAX_CHARS + 1];
    unsigned        flag;
    char            rname[BL_SAM_RNAME_MAX_CHARS + 1];
    size_t          pos;
    unsigned char   mapq;
    char            cigar[BL_SAM_CIGAR_MAX_CHARS + 1];
    char            rnext[BL_SAM_RNAME_MAX_CHARS + 1];
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
