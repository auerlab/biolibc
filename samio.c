#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include "samio.h"

/***************************************************************************
 *  Description:
 *      Read next alignment from a SAM stream
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     sam_alignment_read(FILE *sam_stream, sam_alignment_t *sam_alignment)

{
    char    temp_seq[SAM_SEQ_MAX_CHARS + 1],
	    pos_str[SAM_POS_MAX_DIGITS + 1],
	    *end;
    size_t  len;
    static size_t   previous_pos = 0;
    
    if ( tsv_read_field(sam_stream, sam_alignment->qname, SAM_QNAME_MAX_CHARS, &len) != EOF )
    {
	// Flag
	tsv_skip_field(sam_stream);
	
	// RNAME
	tsv_read_field(sam_stream, sam_alignment->rname, SAM_RNAME_MAX_CHARS, &len);
	
	// POS
	tsv_read_field(sam_stream, pos_str, SAM_POS_MAX_DIGITS, &len);
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
	
	// MAPQ
	tsv_skip_field(sam_stream);
	
	// CIGAR
	tsv_skip_field(sam_stream);
	
	// RNEXT
	tsv_skip_field(sam_stream);
	
	// PNEXT
	tsv_skip_field(sam_stream);
	
	// TLEN
	tsv_skip_field(sam_stream);
	
	// SEQ
	tsv_read_field(sam_stream, temp_seq, SAM_SEQ_MAX_CHARS,
	    &sam_alignment->seq_len);
	if ( (sam_alignment->seq = strdup(temp_seq)) == NULL )
	{
	    fprintf(stderr, "sam_alignment_read(): malloc() failed.\n");
	    exit(EX_UNAVAILABLE);
	}
	
	// QUAL
	// Some SRA CRAMs have 11 fields, most have 12
	if ( tsv_skip_field(sam_stream) == '\t' )
	    while ( getc(sam_stream) != '\n' )
		;
	return 1;
    }
    else
	return 0;
}


/***************************************************************************
 *  Description:
 *      Copy a SAM alignment as efficiently as possible
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_alignment_copy(sam_alignment_t *dest, sam_alignment_t *src)

{
    strlcpy(dest->qname, src->qname, SAM_QNAME_MAX_CHARS);
    strlcpy(dest->rname, src->rname, SAM_RNAME_MAX_CHARS);
    dest->seq = src->seq;
    dest->pos = src->pos;
    dest->seq_len = src->seq_len;
}


/***************************************************************************
 *  Description:
 *      Free memory allocated by sam_alignment_read()
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_alignment_free(sam_alignment_t *sam_alignment)

{
    free(sam_alignment->seq);
}
