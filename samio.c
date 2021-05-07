#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include "samio.h"

/***************************************************************************
 *  Description:
 *      Read next alignment from a SAM stream
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     sam_alignment_read(FILE *sam_stream, sam_alignment_t *sam_alignment)

{
    char    mapq_str[SAM_MAPQ_MAX_CHARS + 1],
	    temp_seq_or_qual[SAM_SEQ_MAX_CHARS + 1],
	    pos_str[BIO_POSITION_MAX_DIGITS + 1],
	    flag_str[SAM_FLAG_MAX_DIGITS + 1],
	    *end;
    size_t  len;
    static size_t   previous_pos = 0;
    int     last_ch;
    
    if ( tsv_read_field(sam_stream, sam_alignment->qname, SAM_QNAME_MAX_CHARS, &len) != EOF )
    {
	// 2 Flag
	tsv_read_field(sam_stream, flag_str, SAM_FLAG_MAX_DIGITS, &len);
	sam_alignment->flag = strtoul(flag_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "sam_read_alignment(): Invalid alignment position: %s\n",
		    flag_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %zu\n", previous_pos);
	    exit(EX_DATAERR);
	}
	
	// 3 RNAME
	tsv_read_field(sam_stream, sam_alignment->rname, SAM_RNAME_MAX_CHARS, &len);
	
	// 4 POS
	tsv_read_field(sam_stream, pos_str, BIO_POSITION_MAX_DIGITS, &len);
	sam_alignment->pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "sam_read_alignment(): Invalid alignment position: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %zu\n", previous_pos);
	    exit(EX_DATAERR);
	}
	previous_pos = sam_alignment->pos;
	
	// 5 MAPQ
	tsv_read_field(sam_stream, mapq_str, SAM_MAPQ_MAX_CHARS, &len);
	sam_alignment->mapq = strtoul(mapq_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "sam_read_alignment(): Invalid alignment mapq: %s\n",
		    mapq_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %zu\n", previous_pos);
	    exit(EX_DATAERR);
	}
	
	// 6 CIGAR
	tsv_skip_field(sam_stream);
	
	// 7 RNEXT
	tsv_skip_field(sam_stream);
	
	// 8 PNEXT
	tsv_skip_field(sam_stream);
	
	// 9 TLEN
	tsv_skip_field(sam_stream);
	
	// 10 SEQ
	tsv_read_field(sam_stream, temp_seq_or_qual, SAM_SEQ_MAX_CHARS,
	    &sam_alignment->seq_len);
	if ( sam_alignment->seq == NULL )
	{
	    //fprintf(stderr, "sam_alignment_read() allocating seq...\n");
	    if ( (sam_alignment->seq = malloc(sam_alignment->seq_len + 1)) == NULL )
	    {
		fprintf(stderr, "sam_alignment_read(): malloc() failed.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	memcpy(sam_alignment->seq, temp_seq_or_qual, sam_alignment->seq_len + 1);
	
	// 11 QUAL, should be last field
	last_ch = tsv_read_field(sam_stream, temp_seq_or_qual, SAM_SEQ_MAX_CHARS,
	    &sam_alignment->qual_len);
	if ( sam_alignment->qual == NULL )
	{
	    //fprintf(stderr, "sam_alignment_read() allocating qual...\n");
	    if ( (sam_alignment->qual = malloc(sam_alignment->qual_len + 1)) == NULL )
	    {
		fprintf(stderr, "sam_alignment_read(): malloc() failed.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	memcpy(sam_alignment->qual, temp_seq_or_qual,
	       sam_alignment->qual_len + 1);

	if ( (sam_alignment->qual_len != 1) &&
	     (sam_alignment->seq_len != sam_alignment->qual_len) )
	    fprintf(stderr, "sam_alignment_read(): Warning: qual_len != seq_len for %s,%zu\n",
		    sam_alignment->rname, sam_alignment->pos);
	
	// Some SRA CRAMs have 11 fields, most have 12
	// Discard everything after the 11th
	if ( last_ch == '\t' )
	    while ( getc(sam_stream) != '\n' )
		;

	/*fprintf(stderr,"sam_alignment_read(): %s,%zu,%zu\n",
		SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
		SAM_SEQ_LEN(sam_alignment));*/
	
	return 1;
    }
    else
	return 0;
}


/***************************************************************************
 *  Description:
 *      Copy a SAM alignment as efficiently as possible
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_alignment_copy(sam_alignment_t *dest, sam_alignment_t *src)

{
    strlcpy(dest->qname, src->qname, SAM_QNAME_MAX_CHARS);
    dest->flag = src->flag;
    strlcpy(dest->rname, src->rname, SAM_RNAME_MAX_CHARS);
    dest->pos = src->pos;
    dest->mapq = src->mapq;
    // FIXME: Add cigar and RNEXT
    dest->cigar = NULL;
    dest->rnext = NULL;
    dest->pnext = src->pnext;
    dest->tlen = src->tlen;
    
    if ( (dest->seq = malloc(src->seq_len + 1)) == NULL )
    {
	fprintf(stderr, "sam_alignment_copy(): malloc() failed.\n");
	exit(EX_UNAVAILABLE);
    }
    memcpy(dest->seq, src->seq, src->seq_len + 1);
    
    if ( (dest->qual = malloc(src->seq_len + 1)) == NULL )
    {
	fprintf(stderr, "sam_alignment_copy(): malloc() failed.\n");
	exit(EX_UNAVAILABLE);
    }
    memcpy(dest->qual, src->qual, src->qual_len + 1);
    
    dest->seq_len = src->seq_len;
    dest->qual_len = src->qual_len;
}


/***************************************************************************
 *  Description:
 *      Free memory allocated by sam_alignment_read()
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_alignment_free(sam_alignment_t *sam_alignment)

{
    if ( sam_alignment->seq != NULL )
	free(sam_alignment->seq);
    if ( sam_alignment->qual != NULL )
	free(sam_alignment->qual);
    // FIXME: Cigar and rnext?
}


/***************************************************************************
 *  Description:
 *      Initialize a sam_alignment_t structure
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_alignment_init(sam_alignment_t *sam_alignment, size_t seq_len)

{
    
    *sam_alignment->qname = '\0';
    sam_alignment->flag = 0;
    *sam_alignment->rname = '\0';
    sam_alignment->pos = 0;
    sam_alignment->mapq = 0;
    sam_alignment->cigar = NULL;
    sam_alignment->rnext = NULL;
    sam_alignment->pnext = 0;
    sam_alignment->tlen = 0;
    if ( seq_len == 0 )
    {
	sam_alignment->seq = NULL;
	sam_alignment->qual = NULL;
    }
    else
    {
	if ( (sam_alignment->seq = malloc(seq_len + 1)) == NULL )
	{
	    fprintf(stderr, "sam_alignment_init(): malloc() failed.\n");
	    exit(EX_UNAVAILABLE);
	}
	if ( (sam_alignment->qual = malloc(seq_len + 1)) == NULL )
	{
	    fprintf(stderr, "sam_alignment_init(): malloc() failed.\n");
	    exit(EX_UNAVAILABLE);
	}
    }
    sam_alignment->seq_len = seq_len;
}
