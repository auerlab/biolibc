
#ifndef _sam_buff_h_
#define _sam_buff_h_

#ifndef sam_h_
#include "sam.h"
#endif

/*
    256k was not enough for a few of the SRA CRAMs.
    NWD976804 needed more than 512k.  Bad data?
    Of 55k samples, only about 11% reached a sam buffer of > 8k.
    Set an upper limit to prevent runaway memory use and error out
    with EX_DATAERR to tell script not to retry.

    From sam-buff-stats script:
    
    Total      54982
    >4096      13723
    >8192       5824
    >16384      2795
    >32768      1612
    >65536       138
    >131072       33
    >262144       12
    >524288        0
*/
#define     SAM_BUFF_START_SIZE     4096
#define     SAM_BUFF_MAX_SIZE       524288

/*
 *  Copied from htslib/sam.h to avoid an htslib dependency.  It should be
 *  safe to assume this will never change, since changing it would break
 *  all existing SAM/BAM/CRAM files.
 */
#define     BAM_FUNMAP  4

// FIXME: Move this to samio when complete?
typedef struct
{
    size_t          buff_size;;
    bl_sam_t **alignments;
    size_t          buffered_count;
    size_t          max_count;
    size_t          previous_pos;
    char            previous_rname[SAM_RNAME_MAX_CHARS + 1];
    
    // Use 64 bits to accommodate large sums
    size_t          mapq_min,
		    mapq_low,
		    mapq_high,
		    mapq_sum,
		    reads_used,
		    total_alignments,
		    trailing_alignments,
		    discarded_alignments,
		    discarded_score_sum,
		    discarded_trailing,
		    min_discarded_score,
		    max_discarded_score,
		    unmapped_alignments;
}   bl_sam_buff_t;

// FIXME: Make sure all fields have accessors and mutators
#define SAM_BUFF_MAPQ_MIN(b)    ((b)->mapq_min)
#define SAM_BUFF_MAPQ_SUM(b)    ((b)->mapq_sum)
#define SAM_BUFF_MAPQ_LOW(b)    ((b)->mapq_low)
#define SAM_BUFF_MAPQ_HIGH(b)   ((b)->mapq_high)
#define SAM_BUFF_READS_USED(b)  ((b)->reads_used)

#define SAM_BUFF_TOTAL_ALIGNMENTS(b)        ((b)->total_alignments)
#define SAM_BUFF_UNMAPPED_ALIGNMENTS(b)     ((b)->unmapped_alignments)
#define SAM_BUFF_DISCARDED_ALIGNMENTS(b)    ((b)->discarded_alignments)
#define SAM_BUFF_DISCARDED_TRAILING(b)      ((b)->discarded_trailing)
#define SAM_BUFF_MIN_DISCARDED_SCORE(b)     ((b)->min_discarded_score)
#define SAM_BUFF_MAX_DISCARDED_SCORE(b)     ((b)->max_discarded_score)
#define SAM_BUFF_DISCARDED_SCORE_SUM(b)     ((b)->discarded_score_sum)
#define SAM_BUFF_BUFFERED_COUNT(b)          ((b)->buffered_count)
#define SAM_BUFF_MAX_COUNT(b)               ((b)->max_count)
#define SAM_BUFF_ALIGNMENTS(b,c)            ((b)->alignments[c])

#define SAM_BUFF_INC_TOTAL_ALIGNMENTS(b)    (++(b)->total_alignments)
#define SAM_BUFF_INC_TRAILING_ALIGNMENTS(b) (++(b)->trailing_alignments)
#define SAM_BUFF_INC_DISCARDED_TRAILING(b)  (++(b)->discarded_trailing)

/* sam-buff.c */
void sam_buff_check_order(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void sam_buff_init(bl_sam_buff_t *sam_buff, unsigned int mapq_min);
void sam_buff_add_alignment(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void sam_buff_out_of_order(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void sam_buff_free_alignment(bl_sam_buff_t *sam_buff, size_t c);
void sam_buff_shift(bl_sam_buff_t *sam_buff, size_t c);
_Bool sam_buff_alignment_ok(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);

#endif  // _sam_buff_h_
