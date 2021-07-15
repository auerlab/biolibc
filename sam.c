#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include <xtend.h>      // strlcpy() on Linux
#include "sam.h"
#include "biolibc.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read next alignment (line) from a SAM stream.
 *
 *      If field_mask is not SAM_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in sam_alignment.
 *      That field in the structure is then populated with an appropriate
 *      marker, such as '.'.  Possible mask values are:
 *
 *      SAM_FIELD_ALL
 *      SAM_FIELD_QNAME
 *      SAM_FIELD_FLAG
 *      SAM_FIELD_RNAME
 *      SAM_FIELD_POS
 *      SAM_FIELD_MAPQ
 *      SAM_FIELD_CIGAR
 *      SAM_FIELD_RNEXT
 *      SAM_FIELD_PNEXT
 *      SAM_FIELD_TLEN
 *      SAM_FIELD_SEQ
 *      SAM_FIELD_QUAL
 *
 *  Arguments:
 *      sam_stream:     A FILE stream from which to read the line
 *      sam_alignment:  Pointer to a bl_sam_t structure
 *      field_mask:     Bit mask indicating which fields to store in sam_alignment
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered after a complete feature
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      sam_read_alignment(stdin, &sam_alignment, SAM_FIELD_ALL);
 *      sam_read_alignment(sam_stream, &sam_alignment,
 *                         SAM_FIELD_QNAME|SAM_FIELD_POS|SAM_FIELD_TLEN);
 *
 *  See also:
 *      sam_alignment_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     sam_read_alignment(FILE *sam_stream, bl_sam_t *sam_alignment,
			   sam_field_mask_t field_mask)

{
    char    mapq_str[SAM_MAPQ_MAX_CHARS + 1],
	    temp_seq_or_qual[SAM_SEQ_MAX_CHARS + 1],
	    pos_str[BL_POSITION_MAX_DIGITS + 1],
	    flag_str[SAM_FLAG_MAX_DIGITS + 1],
	    *end;
    size_t  len;
    static size_t   previous_pos = 0;
    int     delim;
    
    if ( field_mask & SAM_FIELD_QNAME )
	delim = tsv_read_field(sam_stream, sam_alignment->qname,
			SAM_QNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	*sam_alignment->qname = '\0';
    }
    if ( delim == EOF )
	return BL_READ_EOF;

    // 2 Flag
    if ( field_mask & SAM_FIELD_FLAG )
	delim = tsv_read_field(sam_stream, flag_str, SAM_FLAG_MAX_DIGITS, &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading flag: %s.\n",
		flag_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & SAM_FIELD_FLAG )
    {
	sam_alignment->flag = strtoul(flag_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "sam_read_alignment(): Invalid position: %s\n",
		    flag_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %zu\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	sam_alignment->flag = 0;    // FIXME: Is there a better choice?
    
    // 3 RNAME
    if ( field_mask & SAM_FIELD_RNAME )
	delim = tsv_read_field(sam_stream, sam_alignment->rname,
			       SAM_RNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	*sam_alignment->rname = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading rname: %s.\n",
		sam_alignment->rname);
	return BL_READ_TRUNCATED;
    }
    
    // 4 POS
    if ( field_mask & SAM_FIELD_POS )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading pos: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & SAM_FIELD_POS )
    {
	sam_alignment->pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "sam_read_alignment(): Invalid position: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %zu\n", previous_pos);
	    exit(EX_DATAERR);
	}
	previous_pos = sam_alignment->pos;
    }
    else
	sam_alignment->pos = 0;
    
    // 5 MAPQ
    if ( field_mask & SAM_FIELD_MAPQ )
	delim = tsv_read_field(sam_stream, mapq_str, SAM_MAPQ_MAX_CHARS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading mapq: %s.\n",
		mapq_str);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & SAM_FIELD_MAPQ )
    {
	sam_alignment->mapq = strtoul(mapq_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "sam_read_alignment(): Invalid mapq: %s\n",
		    mapq_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %zu\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	sam_alignment->mapq = 0;
    
    // 6 CIGAR
    if ( field_mask & SAM_FIELD_CIGAR )
	delim = tsv_read_field(sam_stream, sam_alignment->cigar,
			       SAM_CIGAR_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	*sam_alignment->cigar = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading cigar: %s.\n",
		sam_alignment->cigar);
	return BL_READ_TRUNCATED;
    }
    
    // 7 RNEXT
    if ( field_mask & SAM_FIELD_RNEXT )
	delim = tsv_read_field(sam_stream, sam_alignment->rnext,
			       SAM_RNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	*sam_alignment->rnext = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading rnext: %s.\n",
		sam_alignment->rnext);
	return BL_READ_TRUNCATED;
    }
    
    // 8 PNEXT
    if ( field_mask & SAM_FIELD_PNEXT )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading pnext: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & SAM_FIELD_PNEXT )
    {
	sam_alignment->pnext = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "sam_read_alignment(): Invalid pnext: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %zu\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	sam_alignment->pnext = 0;
    
    // 9 TLEN
    if ( field_mask & SAM_FIELD_TLEN )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading tlen: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & SAM_FIELD_TLEN )
    {
	sam_alignment->tlen = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "sam_read_alignment(): Invalid tlen: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %zu\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	sam_alignment->tlen = 0;
    
    // 10 SEQ
    if ( field_mask & SAM_FIELD_SEQ )
	delim = tsv_read_field(sam_stream, temp_seq_or_qual, SAM_SEQ_MAX_CHARS,
			       &sam_alignment->seq_len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	// sam_alignment->seq = NULL;   Mem leak if allocated elsewhere
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading seq: %s.\n",
		temp_seq_or_qual);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & SAM_FIELD_SEQ )
    {
	// May be allocated by sam_init_alignment() or sam_copy_alignment()
	if ( sam_alignment->seq == NULL )
	{
	    if ( (sam_alignment->seq = xt_malloc(sam_alignment->seq_len + 1,
		    sizeof(*sam_alignment->seq))) == NULL )
	    {
		fprintf(stderr, "sam_read_alignment(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	memcpy(sam_alignment->seq, temp_seq_or_qual, sam_alignment->seq_len + 1);
    }
    
    // 11 QUAL, should be last field
    if ( field_mask & SAM_FIELD_QUAL )
	delim = tsv_read_field(sam_stream, temp_seq_or_qual, SAM_SEQ_MAX_CHARS,
			       &sam_alignment->qual_len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	// sam_alignment->qual = NULL; Mem leak if allocated elsewhere
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "sam_read_alignment(): Got EOF reading qual: %s.\n",
		temp_seq_or_qual);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & SAM_FIELD_QUAL )
    {
	// May be allocated by sam_init_alignment() or sam_copy_alignment()
	if ( sam_alignment->qual == NULL )
	{
	    if ( (sam_alignment->qual = xt_malloc(sam_alignment->qual_len + 1,
		    sizeof(*sam_alignment->qual))) == NULL )
	    {
		fprintf(stderr, "sam_read_alignment(): Could not allocate qual.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	memcpy(sam_alignment->qual, temp_seq_or_qual,
	       sam_alignment->qual_len + 1);
    
	if ( (sam_alignment->qual_len != 1) &&
	     (sam_alignment->seq_len != sam_alignment->qual_len) )
	    fprintf(stderr, "sam_read_alignment(): Warning: qual_len != seq_len for %s,%zu\n",
		    sam_alignment->rname, sam_alignment->pos);
    }
    
    // Some SRA CRAMs have 11 fields, most have 12
    // Discard everything after the 11th
    if ( delim == '\t' )
	while ( getc(sam_stream) != '\n' )
	    ;

    /*fprintf(stderr,"sam_read_alignment(): %s,%zu,%zu\n",
	    SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
	    SAM_SEQ_LEN(sam_alignment));*/
    
    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc
 *
 *  Description:
 *      Copy a SAM alignment as efficiently as possible, allocating memory
 *      as needed.
 *
 *  Arguments:
 *      dest:   Pointer to bl_sam_t structure to receive copy
 *      src:    Pointer to bl_sam_t structure to be copied
 *
 *  See also:
 *      sam_read_alignment(3), sam_init_alignment(3), sam_free_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_copy_alignment(bl_sam_t *dest, bl_sam_t *src)

{
    strlcpy(dest->qname, src->qname, SAM_QNAME_MAX_CHARS + 1);
    dest->flag = src->flag;
    strlcpy(dest->rname, src->rname, SAM_RNAME_MAX_CHARS + 1);
    dest->pos = src->pos;
    dest->mapq = src->mapq;
    // FIXME: Add cigar and RNEXT
    strlcpy(dest->cigar, src->cigar, SAM_CIGAR_MAX_CHARS + 1);
    strlcpy(dest->rnext, src->rnext, SAM_RNAME_MAX_CHARS + 1);
    dest->pnext = src->pnext;
    dest->tlen = src->tlen;
    
    if ( (dest->seq = xt_malloc(src->seq_len + 1,
	    sizeof(*dest->seq))) == NULL )
    {
	fprintf(stderr, "sam_copy_alignment(): Could not allocate seq.\n");
	exit(EX_UNAVAILABLE);
    }
    memcpy(dest->seq, src->seq, src->seq_len + 1);
    
    if ( (dest->qual = xt_malloc(src->seq_len + 1,
	    sizeof(*dest->qual))) == NULL )
    {
	fprintf(stderr, "sam_copy_alignment(): Could not allocate qual.\n");
	exit(EX_UNAVAILABLE);
    }
    memcpy(dest->qual, src->qual, src->qual_len + 1);
    
    dest->seq_len = src->seq_len;
    dest->qual_len = src->qual_len;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc
 *
 *  Description:
 *      Free memory allocated by sam_read_alignment() or
 *      sam_init_alignment().
 *
 *  Arguments:
 *      sam_alignment:  Pointer to bl_sam_t structure to be freed.
 *
 *  See also:
 *      sam_read_alignment(3), sam_init_alignment(3), sam_copy_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_free_alignment(bl_sam_t *sam_alignment)

{
    if ( sam_alignment->seq != NULL )
	free(sam_alignment->seq);
    if ( sam_alignment->qual != NULL )
	free(sam_alignment->qual);
    // FIXME: Cigar and rnext?
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc
 *
 *  Description:
 *      Initialize a bl_sam_t structure, allocating memory for
 *      sequence and quality strings according to seq_len.  Passing a
 *      seq_len of 0 prevents memory allocation from occurring.
 *
 *      Only SAM_FIELD_SEQ and SAM_FIELD_QUAL are meaningful bits in
 *      field_mask, as they determine whether memory is allocated.  All
 *      other fields are unconditionally initialized to 0, NULL, or blank.
 *
 *  Arguments:
 *      sam_alignment:  Pointer to bl_sam_t structure to initialize
 *      seq_len:        Length of sequence and quality strings
 *      field_mask:     Bit mask indicating which fields will be used
 *
 *  See also:
 *      sam_read_alignment(3), sam_free_alignment(3), sam_copy_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_init_alignment(bl_sam_t *sam_alignment, size_t seq_len,
			   sam_field_mask_t field_mask)

{
    *sam_alignment->qname = '\0';
    sam_alignment->flag = 0;
    *sam_alignment->rname = '\0';
    sam_alignment->pos = 0;
    sam_alignment->mapq = 0;
    *sam_alignment->cigar = '\0';
    *sam_alignment->rnext = '\0';
    sam_alignment->pnext = 0;
    sam_alignment->tlen = 0;
    if ( seq_len == 0 )
    {
	sam_alignment->seq = NULL;
	sam_alignment->qual = NULL;
    }
    else
    {
	if ( seq_len != 0 )
	{
	    if ( (field_mask & SAM_FIELD_SEQ) && 
		 ((sam_alignment->seq = xt_malloc(seq_len + 1,
			sizeof(*sam_alignment->seq))) == NULL) )
	    {
		fprintf(stderr, "sam_init_alignment(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	    if ( (field_mask & SAM_FIELD_QUAL) &&
		 ((sam_alignment->qual = xt_malloc(seq_len + 1,
			sizeof(*sam_alignment->qual))) == NULL) )
	    {
		fprintf(stderr, "sam_init_alignment(): Could not allocate qual.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
    }
    sam_alignment->seq_len = seq_len;
}
