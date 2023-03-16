#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <xtend/string.h>   // strlcpy() on Linux
#include <xtend/mem.h>
#include <xtend/math.h>     // XT_MIN()
#include "sam-buff.h"
#include "biostring.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Check that the newly read SAM alignment comes after the previous
 *      one, assuming the input is sorted first by chrom and then
 *      position.  The previous chrom and position are stored in
 *      sam_buff (and initialized so that the first SAM alignment read is
 *      always OK).
 *  
 *  Arguments:
 *      sam_buff        Pointer to a SAM buffer with recent alignments
 *      sam_alignment   Pointer to the most recently read alignment
 *
 *  See also:
 *      bl_sam_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_check_order(bl_sam_buff_t *sam_buff,
			     bl_sam_t *sam_alignment)

{
    /*fprintf(stderr, "Previous SAM: %s %zu, Current SAM: %s %zu\n",
	    sam_buff->previous_rname, sam_buff->previous_pos,
	    sam_alignment->rname, sam_alignment->pos);*/
    if ( strcmp(sam_alignment->rname, sam_buff->previous_rname) == 0 )
    {
	// Silly to assign when already ==, but sillier to add another check
	if (sam_alignment->pos < sam_buff->previous_pos )
	    bl_sam_buff_out_of_order(sam_buff, sam_alignment);
	else
	    sam_buff->previous_pos = sam_alignment->pos;
    }
    else if ( bl_chrom_name_cmp(sam_alignment->rname,
				  sam_buff->previous_rname) < 0 )
	bl_sam_buff_out_of_order(sam_buff, sam_alignment);
    else
    {
	strlcpy(sam_buff->previous_rname, sam_alignment->rname, BL_SAM_RNAME_MAX_CHARS);
	sam_buff->previous_pos = sam_alignment->pos;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a SAM alignment buffer for holding recently read SAM
 *      alignments.  This is useful, for example, when scanning a SAM
 *      stream for alignments overlapping a certain region or position.
 *      The buffer array is set to a
 *      reasonable initial size and extended as far as BL_SAM_BUFF_MAX_SIZE
 *      by bl_sam_buff_add_alignment(3) if needed.  A minimum MAPQ value
 *      is stored in the bl_sam_buff_t structure for filtering with
 *      bl_sam_buff_alignment_ok(3).
 *  
 *  Arguments:
 *      sam_buff    Pointer to a the bl_sam_buff_t structure to initialize
 *      mapq_min    User-selected minimum MAPQ value
 *
 *  See also:
 *      bl_sam_buff_check_order(3), bl_sam_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_init(bl_sam_buff_t *sam_buff, unsigned int mapq_min,
			 size_t max_alignments)

{
    size_t  c;
    
    sam_buff->buff_size = BL_SAM_BUFF_START_SIZE;
    sam_buff->max_alignments = max_alignments;
    sam_buff->buffered_count = 0;
    sam_buff->max_count = 0;
    sam_buff->previous_pos = 0;
    *sam_buff->previous_rname = '\0';
    
    sam_buff->mapq_min = mapq_min;
    sam_buff->mapq_low = UINT64_MAX;
    sam_buff->mapq_high = 0;
    sam_buff->mapq_sum = 0;
    sam_buff->reads_used = 0;

    sam_buff->total_alignments = 0;
    sam_buff->trailing_alignments = 0;
    sam_buff->discarded_alignments = 0;
    sam_buff->discarded_score_sum = 0;
    sam_buff->min_discarded_score = SIZE_MAX;
    sam_buff->max_discarded_score = 0;
    sam_buff->discarded_trailing = 0;
    sam_buff->unmapped_alignments = 0;
    
    /*
     *  Dynamically allocating the pointers is probably senseless since they
     *  take very little space compared to the alignment data.  By the time
     *  the pointer array takes a significant amount of RAM, you're probably
     *  already thrashing to accomodate the sequence data.  The pointer array
     *  size is capped by BL_SAM_BUFF_MAX_SIZE to prevent memory exhaustion.
     *  We may save a few megabytes with this, though.
     */
    sam_buff->alignments =
	(bl_sam_t **)xt_malloc(sam_buff->buff_size,
				   sizeof(bl_sam_t **));
    for (c = 0; c < sam_buff->buff_size; ++c)
	sam_buff->alignments[c] = NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Add a new alignment to the buffer, expanding the array as needed
 *      up to BL_SAM_BUFF_MAX_SIZE.
 *  
 *  Arguments:
 *      sam_buff    Pointer to bl_sam_buff_t structure where alignments are buffered
 *      sam_alignment   New SAM alignment to add to buffer
 *
 *  See also:
 *      bl_sam_buff_init(3), bl_sam_buff_check_order(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_buff_add_alignment(bl_sam_buff_t *sam_buff,
			       bl_sam_t *sam_alignment)

{
    size_t  old_buff_size,
	    c;

    bl_sam_buff_check_order(sam_buff, sam_alignment);
    
    sam_buff->mapq_low = XT_MIN(sam_buff->mapq_low, BL_SAM_MAPQ(sam_alignment));
    sam_buff->mapq_high = XT_MAX(sam_buff->mapq_high, BL_SAM_MAPQ(sam_alignment));
    sam_buff->mapq_sum += BL_SAM_MAPQ(sam_alignment);
    ++sam_buff->reads_used;

    // Just allocate the static fields, bl_sam_copy() does the rest
    if ( sam_buff->alignments[sam_buff->buffered_count] == NULL )
    {
	//fprintf(stderr, "Allocating alignment #%zu\n", sam_buff->buffered_count);
	sam_buff->alignments[sam_buff->buffered_count] = 
	    xt_malloc(1, sizeof(bl_sam_t));
	if ( sam_buff->alignments[sam_buff->buffered_count] == NULL )
	{
	    fprintf(stderr, "bl_sam_buff_add_alignment(): Could not allocate alignments.\n");
	    exit(EX_UNAVAILABLE);
	}
	// Redundant to bl_sam_copy()
	// bl_sam_init(sam_buff->alignments[sam_buff->buffered_count], 0);
    }
    else
	bl_sam_free(sam_buff->alignments[sam_buff->buffered_count]);

    //fprintf(stderr, "Adding alignment #%zu...\n", sam_buff->buffered_count);
    //fprintf(stderr, "buff_size = %zu\n", sam_buff->buff_size);
    bl_sam_copy(sam_buff->alignments[sam_buff->buffered_count], sam_alignment);
    
    ++sam_buff->buffered_count;

    if ( sam_buff->buffered_count > sam_buff->max_count )
    {
	sam_buff->max_count = sam_buff->buffered_count;
	// fprintf(stderr, "sam_buff->max_count = %zu\n", sam_buff->max_count);
    }
    
    if ( sam_buff->buffered_count == sam_buff->max_alignments )
    {
	fprintf(stderr,
		"bl_sam_buff_add_alignment(): Hit maximum alignments=%zu.\n",
		sam_buff->max_alignments);
	fprintf(stderr, "Aborting add to prevent runaway memory use.\n");
	fprintf(stderr, "Check your SAM input.\n");
	return BL_SAM_BUFF_ADD_FAILED;
    }
    
    if ( sam_buff->buffered_count == sam_buff->buff_size )
    {
	fprintf(stderr,
		"bl_sam_buff_add_alignment(): Hit buff_size=%zu, doubling buffer size.\n",
		sam_buff->buff_size);
	fprintf(stderr, "RNAME: %s  POS: %" PRId64 " LEN: %zu\n",
		BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
		BL_SAM_SEQ_LEN(sam_alignment));
	old_buff_size = sam_buff->buff_size;
	sam_buff->buff_size *= 2;
	sam_buff->alignments =
	    (bl_sam_t **)xt_realloc(sam_buff->alignments,
					sam_buff->buff_size,
					sizeof(bl_sam_t **));
	for (c = old_buff_size; c < sam_buff->buff_size; ++c)
	    sam_buff->alignments[c] = NULL;
    }
    return BL_SAM_BUFF_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Report SAM input sort error and terminate process.
 *  
 *  Arguments:
 *      sam_buff        Pointer to bl_sam_buff_t structure
 *      sam_alignment   Offending alignment out of order with previous
 *
 *  Returns:
 *      Does not return.
 *
 *  See also:
 *      bl_sam_buff_alignment_ok(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_out_of_order(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment)

{
    fprintf(stderr, "Error: SAM input must be sorted by chrom and then position.\n");
    fprintf(stderr, "Found %s,%" PRId64 " after %s,%" PRId64 ".\n",
	    BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
	    sam_buff->previous_rname, sam_buff->previous_pos);
    exit(EX_DATAERR);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free an element of the SAM alignment array by first freeing all
 *      memory allocated by the bl_sam_t structure and then freeing
 *      memory allocated for the structure itself.
 *  
 *  Arguments:
 *      sam_buff    Pointer to the bl_sam_buff_t structure holding alignments
 *      c           Index of the alignment to be freed (0-based)
 *
 *  See also:
 *      bl_sam_buff_init(3), bl_sam_buff_add_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_free_alignment(bl_sam_buff_t *sam_buff, size_t c)

{
    bl_sam_free(sam_buff->alignments[c]);
    bl_sam_init(sam_buff->alignments[c]);
    if ( sam_buff->alignments[c] != NULL )
    {
	free(sam_buff->alignments[c]);
	sam_buff->alignments[c] = NULL;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free nelem SAM alignments at the head of the queue and shift
 *      remaining elements forward nelem positions.
 *  
 *  Arguments:
 *      sam_buff    Pointer to bl_sam_buff_t structure holding alignments
 *      nelem       Number of alignments to free
 *
 *  See also:
 *      bl_sam_buff_free_alignment(3)
 *
 *  FIXME: Use circular queuing for better efficiency
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_shift(bl_sam_buff_t *sam_buff, size_t nelem)

{
    size_t  c;

    /*
     *  FIXME: A circular queue would be more efficient, but won't matter
     *  much to the bottom line to avoid shifting a few pointers on top
     *  of everything else going on here
     */
    
    /* Make susre elements to be removed are freed */
    for (c = 0; c < nelem; ++c)
	bl_sam_buff_free_alignment(sam_buff, c);

    /* Shift elements */
    for (c = 0; c < sam_buff->buffered_count - nelem; ++c)
	sam_buff->alignments[c] = sam_buff->alignments[c + nelem];
    
    /* Clear vacated elements */
    while ( c < sam_buff->buffered_count )
	sam_buff->alignments[c++] = NULL;
    
    sam_buff->buffered_count -= nelem;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Verify that an alignment meets the quality requirements for the
 *      given SAM buffer.  Each bl_sam_buff_t structure contains
 *      specifications for minimum quality scores as well as statistics
 *      on alignments checked, including the number of discarded alignments
 *      and a sum of the alignment scores (for calculating the average
 *      score of discarded alignments).
 *  
 *  Arguments:
 *      sam_buff    Pointer to bl_sam_buff_t structure with quality specs
 *      sam_alignment   Pointer to new alignment
 *
 *  Returns:
 *      true if the alignment meets the quality requirements
 *      false otherwise
 *
 *  See also:
 *      bl_sam_buff_check_order(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/
 
bool    bl_sam_buff_alignment_ok(bl_sam_buff_t *sam_buff,
			      bl_sam_t *sam_alignment)

{
    if ( sam_alignment->flag & BAM_FUNMAP )
    {
	++sam_buff->unmapped_alignments;
#ifdef DEBUG
	fprintf(stderr, "Discarding unmapped read: %s,%" PRId64 ",0x%04x\n",
		BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
		BL_SAM_FLAG(sam_alignment));
#endif
	return false;
    }
    else if ( BL_SAM_MAPQ(sam_alignment) < BL_SAM_BUFF_MAPQ_MIN(sam_buff) )
    {
	++sam_buff->discarded_alignments;
	sam_buff->discarded_score_sum += BL_SAM_MAPQ(sam_alignment);
	if ( BL_SAM_MAPQ(sam_alignment) < sam_buff->min_discarded_score )
	    sam_buff->min_discarded_score = BL_SAM_MAPQ(sam_alignment);
	if ( BL_SAM_MAPQ(sam_alignment) > sam_buff->max_discarded_score )
	    sam_buff->max_discarded_score = BL_SAM_MAPQ(sam_alignment);

#ifdef DEBUG
	fprintf(stderr, "bl_sam_buff_alignment_ok(): Discarding low quality alignment: %s,%" PRId64 " MAPQ=%u\n",
		BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
		BL_SAM_MAPQ(sam_alignment));
#endif
	return false;
    }
    else
	return true;
}
