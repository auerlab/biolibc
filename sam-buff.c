#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>  // MIN()
#include <xtend.h>
#include "sam-buff.h"
#include "biostring.h"

/***************************************************************************
 *  Description:
 *      Check order of sam alignments being read into buffer.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_check_order(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment)

{
    /*fprintf(stderr, "Previous SAM: %s %zu, Current SAM: %s %zu\n",
	    sam_buff->previous_rname, sam_buff->previous_pos,
	    sam_alignment->rname, sam_alignment->pos);*/
    if ( strcmp(sam_alignment->rname, sam_buff->previous_rname) == 0 )
    {
	// Silly to assign when ==, but sillier to add another check
	if (sam_alignment->pos < sam_buff->previous_pos )
	    sam_buff_out_of_order(sam_buff, sam_alignment);
	else
	    sam_buff->previous_pos = sam_alignment->pos;
    }
    else if ( chromosome_name_cmp(sam_alignment->rname, sam_buff->previous_rname) < 0 )
	sam_buff_out_of_order(sam_buff, sam_alignment);
    else
    {
	strlcpy(sam_buff->previous_rname, sam_alignment->rname, SAM_RNAME_MAX_CHARS);
	sam_buff->previous_pos = sam_alignment->pos;
    }
}


/***************************************************************************
 *  Description:
 *      Initialize a SAM alignment buffer
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_init(sam_buff_t *sam_buff, unsigned int mapq_min)

{
    size_t  c;
    
    sam_buff->buff_size = SAM_BUFF_START_SIZE;
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
     *  size is capped by SAM_BUFF_MAX_SIZE to prevent memory exhaustion.
     *  We may save a few megabytes with this, though.
     */
    sam_buff->alignments =
	(sam_alignment_t **)xt_malloc(sam_buff->buff_size,
				   sizeof(sam_alignment_t **));
    for (c = 0; c < sam_buff->buff_size; ++c)
	sam_buff->alignments[c] = NULL;
}


/***************************************************************************
 *  Description:
 *      Add a new alignment to the buffer
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_add_alignment(sam_buff_t *sam_buff,
			       sam_alignment_t *sam_alignment)

{
    size_t  old_buff_size,
	    c;

    sam_buff_check_order(sam_buff, sam_alignment);
    
    sam_buff->mapq_low = MIN(sam_buff->mapq_low, SAM_MAPQ(sam_alignment));
    sam_buff->mapq_high = MAX(sam_buff->mapq_high, SAM_MAPQ(sam_alignment));
    sam_buff->mapq_sum += SAM_MAPQ(sam_alignment);
    ++sam_buff->reads_used;

    // Just allocate the static fields, sam_alignment_copy() does the rest
    if ( sam_buff->alignments[sam_buff->buffered_count] == NULL )
    {
	//fprintf(stderr, "Allocating alignment #%zu\n", sam_buff->buffered_count);
	sam_buff->alignments[sam_buff->buffered_count] = 
	    xt_malloc(1, sizeof(sam_alignment_t));
	if ( sam_buff->alignments[sam_buff->buffered_count] == NULL )
	{
	    fprintf(stderr, "sam_buff_add_alignment(): malloc() failed.\n");
	    exit(EX_UNAVAILABLE);
	}
	// Redundant to sam_alignment_copy()
	// sam_alignment_init(sam_buff->alignments[sam_buff->buffered_count], 0);
    }
    else
	sam_alignment_free(sam_buff->alignments[sam_buff->buffered_count]);
    
    sam_alignment_copy(sam_buff->alignments[sam_buff->buffered_count], sam_alignment);
    
    ++sam_buff->buffered_count;

    if ( sam_buff->buffered_count > sam_buff->max_count )
    {
	sam_buff->max_count = sam_buff->buffered_count;
	// fprintf(stderr, "sam_buff->max_count = %zu\n", sam_buff->max_count);
    }
    
    if ( sam_buff->buffered_count == SAM_BUFF_MAX_SIZE )
    {
	fprintf(stderr,
		"sam_buff_add_alignment(): Hit SAM_BUFF_MAX_SIZE=%u.\n",
		SAM_BUFF_MAX_SIZE);
	fprintf(stderr, "Terminating to prevent runaway memory use.\n");
	fprintf(stderr, "Check your SAM input.\n");
	exit(EX_DATAERR);
    }
    
    if ( sam_buff->buffered_count == sam_buff->buff_size )
    {
	fprintf(stderr,
		"sam_buff_add_alignment(): Hit buff_size=%zu, doubling buffer size.\n",
		sam_buff->buff_size);
	fprintf(stderr, "RNAME: %s  POS: %zu  LEN: %zu\n",
		SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
		SAM_SEQ_LEN(sam_alignment));
	old_buff_size = sam_buff->buff_size;
	sam_buff->buff_size *= 2;
	sam_buff->alignments =
	    (sam_alignment_t **)realloc(sam_buff->alignments,
					sam_buff->buff_size *
					sizeof(sam_alignment_t **));
	for (c = old_buff_size; c < sam_buff->buff_size; ++c)
	    sam_buff->alignments[c] = NULL;
    }
}


/***************************************************************************
 *  Description:
 *      Explain SAM input sort error and exit.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_out_of_order(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment)

{
    fprintf(stderr, "ad2vcf: Error: SAM input must be sorted by chromosome and then position.\n");
    fprintf(stderr, "Found %s,%zu after %s,%zu.\n",
	    SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
	    sam_buff->previous_rname, sam_buff->previous_pos);
    exit(EX_DATAERR);
}


/***************************************************************************
 *  Description:
 *      Free an element of the SAM alignment array
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_free_alignment(sam_buff_t *sam_buff, size_t c)

{
    sam_alignment_free(sam_buff->alignments[c]);
    sam_alignment_init(sam_buff->alignments[c], 0);
    if ( sam_buff->alignments[c] != NULL )
    {
	free(sam_buff->alignments[c]);
	sam_buff->alignments[c] = NULL;
    }
}


/***************************************************************************
 *  Description:
 *      Shift SAM array elements up c positions.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_shift(sam_buff_t *sam_buff, size_t c)

{
    size_t  c2;

    /*
     *  FIXME: A circular queue would be more efficient, but won't matter
     *  much to the bottom line to avoid shifting a few pointers on top
     *  of everything else going on here
     */
    
    /* Make sure elements to be removed are freed */
    for (c2 = 0; c2 < c; ++c2)
	sam_buff_free_alignment(sam_buff, c2);

    /* Shift elements */
    for (c2 = 0; c2 < sam_buff->buffered_count - c; ++c2)
	sam_buff->alignments[c2] = sam_buff->alignments[c2 + c];
    
    /* Clear vacated elements */
    while ( c2 < sam_buff->buffered_count )
	sam_buff->alignments[c2++] = NULL;
    
    sam_buff->buffered_count -= c;
}


bool    sam_buff_alignment_ok(sam_buff_t *sam_buff,
			      sam_alignment_t *sam_alignment)

{
    if ( sam_alignment->flag & BAM_FUNMAP )
    {
	++sam_buff->unmapped_alignments;
#ifdef DEBUG
	fprintf(stderr, "Discarding unmapped read: %s,%zu\n",
		SAM_RNAME(sam_alignment), SAM_POS(sam_alignment));
#endif
	return false;
    }
    else if ( SAM_MAPQ(sam_alignment) < SAM_BUFF_MAPQ_MIN(sam_buff) )
    {
	++sam_buff->discarded_alignments;
	sam_buff->discarded_score_sum += SAM_MAPQ(sam_alignment);
	if ( SAM_MAPQ(sam_alignment) < sam_buff->min_discarded_score )
	    sam_buff->min_discarded_score = SAM_MAPQ(sam_alignment);
	if ( SAM_MAPQ(sam_alignment) > sam_buff->max_discarded_score )
	    sam_buff->max_discarded_score = SAM_MAPQ(sam_alignment);

#ifdef DEBUG
	fprintf(stderr, "Discarding low quality read: %s,%zu MAPQ=%u\n",
		SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
		SAM_MAPQ(sam_alignment));
#endif
	return false;
    }
    else
	return true;
}

