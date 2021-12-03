
#ifndef _sam_buff_h_
#define _sam_buff_h_

#ifndef __bool_true_false_are_defined
#include <stdbool.h>
#endif

#ifndef sam_h_
#include "sam.h"
#endif

#ifndef _biolibc_h_
#include "biolibc.h"
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
#define BL_SAM_BUFF_START_SIZE  4096

#define BL_SAM_BUFF_OK          0
#define BL_SAM_BUFF_ADD_FAILED  1

/*
 *  Copied from htslib/sam.h to avoid an htslib dependency.  It should be
 *  safe to assume this will never change, since changing it would break
 *  all existing SAM/BAM/CRAM files.
 */
#define     BAM_FUNMAP  4

// FIXME: Move this to samio when complete?
typedef struct
{
    size_t      buff_size;
    size_t      max_alignments;
    bl_sam_t    **alignments;
    size_t      buffered_count;
    size_t      max_count;
    uint64_t    previous_pos;
    char        previous_rname[BL_SAM_RNAME_MAX_CHARS + 1];
    
    // Use 64 bits to accommodate large sums
    uint64_t    mapq_min,
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

/*
 *  Generated by /home/bacon/scripts/gen-get-set
 *
 *  Accessor macros.  Use these to access structure members from functions
 *  outside the bl_sam_buff_t class.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
  */

#define BL_SAM_BUFF_BUFF_SIZE(ptr)      ((ptr)->buff_size)
#define BL_SAM_BUFF_ALIGNMENTS(ptr)     ((ptr)->alignments)
#define BL_SAM_BUFF_ALIGNMENTS_AE(ptr,c) ((ptr)->alignments[c])
#define BL_SAM_BUFF_BUFFERED_COUNT(ptr) ((ptr)->buffered_count)
#define BL_SAM_BUFF_MAX_COUNT(ptr)      ((ptr)->max_count)
#define BL_SAM_BUFF_PREVIOUS_POS(ptr)   ((ptr)->previous_pos)
#define BL_SAM_BUFF_PREVIOUS_RNAME(ptr) ((ptr)->previous_rname)
#define BL_SAM_BUFF_PREVIOUS_RNAME_AE(ptr,c) ((ptr)->previous_rname[c])
#define BL_SAM_BUFF_MAPQ_MIN(ptr)       ((ptr)->mapq_min)
#define BL_SAM_BUFF_MAPQ_LOW(ptr)       ((ptr)->mapq_low)
#define BL_SAM_BUFF_MAPQ_HIGH(ptr)      ((ptr)->mapq_high)
#define BL_SAM_BUFF_MAPQ_SUM(ptr)       ((ptr)->mapq_sum)
#define BL_SAM_BUFF_READS_USED(ptr)     ((ptr)->reads_used)
#define BL_SAM_BUFF_TOTAL_ALIGNMENTS(ptr) ((ptr)->total_alignments)
#define BL_SAM_BUFF_TRAILING_ALIGNMENTS(ptr) ((ptr)->trailing_alignments)
#define BL_SAM_BUFF_DISCARDED_ALIGNMENTS(ptr) ((ptr)->discarded_alignments)
#define BL_SAM_BUFF_DISCARDED_SCORE_SUM(ptr) ((ptr)->discarded_score_sum)
#define BL_SAM_BUFF_DISCARDED_TRAILING(ptr) ((ptr)->discarded_trailing)
#define BL_SAM_BUFF_MIN_DISCARDED_SCORE(ptr) ((ptr)->min_discarded_score)
#define BL_SAM_BUFF_MAX_DISCARDED_SCORE(ptr) ((ptr)->max_discarded_score)
#define BL_SAM_BUFF_UNMAPPED_ALIGNMENTS(ptr) ((ptr)->unmapped_alignments)

/*
 *  Generated by /home/bacon/scripts/gen-get-set
 *
 *  Mutator macros for setting with no sanity checking.  Use these to
 *  set structure members from functions outside the bl_sam_buff_t
 *  class.  These macros perform no data validation.  Hence, they achieve
 *  maximum performance where data are guaranteed correct by other means.
 *  Use the mutator functions (same name as the macro, but lower case)
 *  for more robust code with a small performance penalty.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
  */

#define BL_SAM_BUFF_SET_BUFF_SIZE(ptr,val)      ((ptr)->buff_size = (val))
#define BL_SAM_BUFF_SET_ALIGNMENTS(ptr,val)     ((ptr)->alignments = (val))
// FIXME: Assuming all elements should be copied
#define BL_SAM_BUFF_SET_ALIGNMENTS_CPY(ptr,val,array_size) \
    for (size_t c = 0; c < (array_size); ++c) (ptr)->alignments[c] = val[c];
#define BL_SAM_BUFF_SET_ALIGNMENTS_AE(ptr,c,val) ((ptr)->alignments[c] = (val))
#define BL_SAM_BUFF_SET_BUFFERED_COUNT(ptr,val) ((ptr)->buffered_count = (val))
#define BL_SAM_BUFF_SET_MAX_COUNT(ptr,val)      ((ptr)->max_count = (val))
#define BL_SAM_BUFF_SET_PREVIOUS_POS(ptr,val)   ((ptr)->previous_pos = (val))
// FIXME: Assuming char array is a null-terminated string
#define BL_SAM_BUFF_SET_PREVIOUS_RNAME_CPY(ptr,val,array_size) strlcpy((ptr)->previous_rname,val,array_size)
#define BL_SAM_BUFF_SET_PREVIOUS_RNAME_AE(ptr,c,val) ((ptr)->previous_rname[c] = (val))
#define BL_SAM_BUFF_SET_MAPQ_MIN(ptr,val)       ((ptr)->mapq_min = (val))
#define BL_SAM_BUFF_SET_MAPQ_LOW(ptr,val)       ((ptr)->mapq_low = (val))
#define BL_SAM_BUFF_SET_MAPQ_HIGH(ptr,val)      ((ptr)->mapq_high = (val))
#define BL_SAM_BUFF_SET_MAPQ_SUM(ptr,val)       ((ptr)->mapq_sum = (val))
#define BL_SAM_BUFF_SET_READS_USED(ptr,val)     ((ptr)->reads_used = (val))
#define BL_SAM_BUFF_SET_TOTAL_ALIGNMENTS(ptr,val) ((ptr)->total_alignments = (val))
#define BL_SAM_BUFF_SET_TRAILING_ALIGNMENTS(ptr,val) ((ptr)->trailing_alignments = (val))
#define BL_SAM_BUFF_SET_DISCARDED_ALIGNMENTS(ptr,val) ((ptr)->discarded_alignments = (val))
#define BL_SAM_BUFF_SET_DISCARDED_SCORE_SUM(ptr,val) ((ptr)->discarded_score_sum = (val))
#define BL_SAM_BUFF_SET_DISCARDED_TRAILING(ptr,val) ((ptr)->discarded_trailing = (val))
#define BL_SAM_BUFF_SET_MIN_DISCARDED_SCORE(ptr,val) ((ptr)->min_discarded_score = (val))
#define BL_SAM_BUFF_SET_MAX_DISCARDED_SCORE(ptr,val) ((ptr)->max_discarded_score = (val))
#define BL_SAM_BUFF_SET_UNMAPPED_ALIGNMENTS(ptr,val) ((ptr)->unmapped_alignments = (val))

/* Not generated by gen-get-set */
#define BL_SAM_BUFF_INC_TOTAL_ALIGNMENTS(b)    (++(b)->total_alignments)
#define BL_SAM_BUFF_INC_TRAILING_ALIGNMENTS(b) (++(b)->trailing_alignments)
#define BL_SAM_BUFF_INC_DISCARDED_TRAILING(b)  (++(b)->discarded_trailing)

/* sam-buff.c */
void bl_sam_buff_check_order(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void bl_sam_buff_init(bl_sam_buff_t *sam_buff, unsigned int mapq_min, size_t max_alignments);
int bl_sam_buff_add_alignment(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void bl_sam_buff_out_of_order(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);
void bl_sam_buff_free_alignment(bl_sam_buff_t *sam_buff, size_t c);
void bl_sam_buff_shift(bl_sam_buff_t *sam_buff, size_t nelem);
bool bl_sam_buff_alignment_ok(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment);

/* sam-buff-mutators.c */
int bl_sam_buff_set_buff_size(bl_sam_buff_t *bl_sam_buff_ptr, size_t new_buff_size);
int bl_sam_buff_set_alignments(bl_sam_buff_t *bl_sam_buff_ptr, bl_sam_t **new_alignments);
int bl_sam_buff_set_alignments_ae(bl_sam_buff_t *bl_sam_buff_ptr, size_t c, bl_sam_t *new_alignments_element);
int bl_sam_buff_set_alignments_cpy(bl_sam_buff_t *bl_sam_buff_ptr, bl_sam_t **new_alignments, size_t array_size);
int bl_sam_buff_set_buffered_count(bl_sam_buff_t *bl_sam_buff_ptr, size_t new_buffered_count);
int bl_sam_buff_set_max_count(bl_sam_buff_t *bl_sam_buff_ptr, size_t new_max_count);
int bl_sam_buff_set_previous_pos(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_previous_pos);
int bl_sam_buff_set_previous_rname_ae(bl_sam_buff_t *bl_sam_buff_ptr, size_t c, char new_previous_rname_element);
int bl_sam_buff_set_previous_rname_cpy(bl_sam_buff_t *bl_sam_buff_ptr, char new_previous_rname[], size_t array_size);
int bl_sam_buff_set_mapq_min(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_mapq_min);
int bl_sam_buff_set_mapq_low(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_mapq_low);
int bl_sam_buff_set_mapq_high(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_mapq_high);
int bl_sam_buff_set_mapq_sum(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_mapq_sum);
int bl_sam_buff_set_reads_used(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_reads_used);
int bl_sam_buff_set_total_alignments(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_total_alignments);
int bl_sam_buff_set_trailing_alignments(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_trailing_alignments);
int bl_sam_buff_set_discarded_alignments(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_discarded_alignments);
int bl_sam_buff_set_discarded_score_sum(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_discarded_score_sum);
int bl_sam_buff_set_discarded_trailing(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_discarded_trailing);
int bl_sam_buff_set_min_discarded_score(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_min_discarded_score);
int bl_sam_buff_set_max_discarded_score(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_max_discarded_score);
int bl_sam_buff_set_unmapped_alignments(bl_sam_buff_t *bl_sam_buff_ptr, uint64_t new_unmapped_alignments);

#endif  // _sam_buff_h_
