
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
#define     BL_SAM_BUFF_START_SIZE     4096
#define     BL_SAM_BUFF_MAX_SIZE       524288

/*
 *  Copied from htslib/sam.h to avoid an htslib dependency.  It should be
 *  safe to assume this will never change, since changing it would break
 *  all existing SAM/BAM/CRAM files.
 */
#define     BAM_FUNMAP  4

// FIXME: Move this to samio when complete?
typedef struct
{
    size_t      buff_size;;
    bl_sam_t    **alignments;
    size_t      buffered_count;
    size_t      max_count;
    size_t      previous_pos;
    char        previous_rname[BL_SAM_RNAME_MAX_CHARS + 1];
    
    // Use 64 bits to accommodate large sums
    size_t      mapq_min,
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

#define BL_SAM_BUFF_BUFF_SIZE(ptr)              ((ptr)->buff_size)
#define BL_SAM_BUFF_ALIGNMENTS(ptr)             ((ptr)->alignments)
#define BL_SAM_BUFF_BUFFERED_COUNT(ptr)         ((ptr)->buffered_count)
#define BL_SAM_BUFF_MAX_COUNT(ptr)              ((ptr)->max_count)
#define BL_SAM_BUFF_PREVIOUS_POS(ptr)           ((ptr)->previous_pos)
#define BL_SAM_BUFF_PREVIOUS_RNAME(ptr)         ((ptr)->previous_rname)
#define BL_SAM_BUFF_MAPQ_MIN(ptr)               ((ptr)->mapq_min)
#define BL_SAM_BUFF_MAPQ_LOW(ptr)               ((ptr)->mapq_low)
#define BL_SAM_BUFF_MAPQ_HIGH(ptr)              ((ptr)->mapq_high)
#define BL_SAM_BUFF_MAPQ_SUM(ptr)               ((ptr)->mapq_sum)
#define BL_SAM_BUFF_READS_USED(ptr)             ((ptr)->reads_used)
#define BL_SAM_BUFF_TOTAL_ALIGNMENTS(ptr)       ((ptr)->total_alignments)
#define BL_SAM_BUFF_TRAILING_ALIGNMENTS(ptr)    ((ptr)->trailing_alignments)
#define BL_SAM_BUFF_DISCARDED_ALIGNMENTS(ptr)   ((ptr)->discarded_alignments)
#define BL_SAM_BUFF_DISCARDED_SCORE_SUM(ptr)    ((ptr)->discarded_score_sum)
#define BL_SAM_BUFF_DISCARDED_TRAILING(ptr)     ((ptr)->discarded_trailing)
#define BL_SAM_BUFF_MIN_DISCARDED_SCORE(ptr)    ((ptr)->min_discarded_score)
#define BL_SAM_BUFF_MAX_DISCARDED_SCORE(ptr)    ((ptr)->max_discarded_score)
#define BL_SAM_BUFF_UNMAPPED_ALIGNMENTS(ptr)    ((ptr)->unmapped_alignments)

#define BL_SAM_BUFF_SET_BUFF_SIZE(ptr,buff_size)        ((ptr)->buff_size = (buff_size))
#define BL_SAM_BUFF_SET_ALIGNMENTS(ptr,alignments)      ((ptr)->alignments = (alignments))
#define BL_SAM_BUFF_SET_BUFFERED_COUNT(ptr,buffered_count)  ((ptr)->buffered_count = (buffered_count))
#define BL_SAM_BUFF_SET_MAX_COUNT(ptr,max_count)        ((ptr)->max_count = (max_count))
#define BL_SAM_BUFF_SET_PREVIOUS_POS(ptr,previous_pos)  ((ptr)->previous_pos = (previous_pos))
#define BL_SAM_BUFF_SET_PREVIOUS_RNAME(ptr,previous_rname)  strlcpy(ptr->previous_rname,previous_rname,BL_SAM_RNAME_MAX_CHARS+1)
#define BL_SAM_BUFF_SET_MAPQ_MIN(ptr,mapq_min)          ((ptr)->mapq_min = (mapq_min))
#define BL_SAM_BUFF_SET_MAPQ_LOW(ptr,mapq_low)          ((ptr)->mapq_low = (mapq_low))
#define BL_SAM_BUFF_SET_MAPQ_HIGH(ptr,mapq_high)        ((ptr)->mapq_high = (mapq_high))
#define BL_SAM_BUFF_SET_MAPQ_SUM(ptr,mapq_sum)          ((ptr)->mapq_sum = (mapq_sum))
#define BL_SAM_BUFF_SET_READS_USED(ptr,reads_used)      ((ptr)->reads_used = (reads_used))
#define BL_SAM_BUFF_SET_TOTAL_ALIGNMENTS(ptr,total_alignments)  ((ptr)->total_alignments = (total_alignments))
#define BL_SAM_BUFF_SET_TRAILING_ALIGNMENTS(ptr,trailing_alignments)    ((ptr)->trailing_alignments = (trailing_alignments))
#define BL_SAM_BUFF_SET_DISCARDED_ALIGNMENTS(ptr,discarded_alignments)  ((ptr)->discarded_alignments = (discarded_alignments))
#define BL_SAM_BUFF_SET_DISCARDED_SCORE_SUM(ptr,discarded_score_sum)    ((ptr)->discarded_score_sum = (discarded_score_sum))
#define BL_SAM_BUFF_SET_DISCARDED_TRAILING(ptr,discarded_trailing)  ((ptr)->discarded_trailing = (discarded_trailing))
#define BL_SAM_BUFF_SET_MIN_DISCARDED_SCORE(ptr,min_discarded_score)    ((ptr)->min_discarded_score = (min_discarded_score))
#define BL_SAM_BUFF_SET_MAX_DISCARDED_SCORE(ptr,max_discarded_score)    ((ptr)->max_discarded_score = (max_discarded_score))
#define BL_SAM_BUFF_SET_UNMAPPED_ALIGNMENTS(ptr,unmapped_alignments)    ((ptr)->unmapped_alignments = (unmapped_alignments))

#define BL_SAM_BUFF_INC_TOTAL_ALIGNMENTS(b)    (++(b)->total_alignments)
#define BL_SAM_BUFF_INC_TRAILING_ALIGNMENTS(b) (++(b)->trailing_alignments)
#define BL_SAM_BUFF_INC_DISCARDED_TRAILING(b)  (++(b)->discarded_trailing)

/* sam-buff.c */
void sam_buff_check_order(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void sam_buff_init(bl_sam_buff_t *sam_buff, unsigned int mapq_min);
void sam_buff_add_alignment(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void sam_buff_out_of_order(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void sam_buff_free_alignment(bl_sam_buff_t *sam_buff, size_t c);
void sam_buff_shift(bl_sam_buff_t *sam_buff, size_t c);
bool sam_buff_alignment_ok(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);

#endif  // _sam_buff_h_
