
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <xtend.h>
#include "fasta.h"
#include "biolibc.h"


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read a FASTA record from a FILE stream.  Each record must begin
 *      with a description line (beginning with '>'), which is then
 *      followed by one or more lines of sequence data.  The end of the
 *      sequence is marked either by the next description line or EOF.
 *  
 *  Arguments:
 *      fasta_stream:   FILE stream from which FASTA data are read
 *      record:         Pointer to a bl_fast_t structure to receive data
 *
 *  Returns:
 *      BL_READ_OK upon successful read of description and sequence
 *      BL_READ_BAD_DATA if something is amiss with input format
 *      BL_READ_EOF if no more data are available
 *
 *  Examples:
 *      bl_fasta_t  rec;
 *
 *      while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fasta_write(stdout, &rec, BL_FASTA_LINE_UNLIMITED);
 *
 *  See also:
 *      bl_fasta_write(3), bl_fastq_read(3), bl_fastq_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fasta_read(FILE *fasta_stream, bl_fasta_t *record)

{
    int     ch;
    char    *p;
    
    /* Skip comment lines */
    while ( (ch = getc(fasta_stream)) == ';' )
	while ( getc(fasta_stream) != '\n' )
	    ;
    
    if ( ch == EOF )
	return BL_READ_EOF;
    
    /* Every record should begin with a '>' */
    if ( ch == '>' )    // Desc
    {
	record->desc_array_size = 1024;
	record->desc = xt_malloc(record->desc_array_size, sizeof(*record->desc));
	if ( record->desc == NULL )
	{
	    fprintf(stderr, "bl_fasta_read(): Could not allocate desc.\n");
	    exit(EX_UNAVAILABLE);
	}

	p = record->desc;
	while ( ((ch = getc(fasta_stream)) != '\n') && (ch != EOF) )
	{
	    *p++ = ch;
	    if ( p - record->desc == record->desc_array_size )
	    {
		record->desc_array_size *= 2;
		record->desc = xt_realloc(record->desc, record->desc_array_size,
		    sizeof(*record->desc));
		if ( record->desc == NULL )
		{
		    fprintf(stderr, "bl_fasta_read(): Could not reallocate desc.\n");
		    exit(EX_UNAVAILABLE);
		}
	    }
	}
	*p = '\0';
	record->desc_len = p - record->desc;

	/* Trim array */
	record->desc = xt_realloc(record->desc, record->desc_len+1,
	    sizeof(*record->desc));
	//printf("Trimmed desc to %zu\n", record->desc_len+1);
	
	/* Should not encounter EOF while reading description line */
	/* Every description should be followed by at least one seq line */
	if ( ch == EOF )
	    return BL_READ_TRUNCATED;

	/*
	 *  Read sequence lines
	 */
	
	record->seq_array_size = 1024;
	record->seq = xt_malloc(record->seq_array_size, sizeof(*record->seq));
	if ( record->seq == NULL )
	{
	    fprintf(stderr, "bl_fasta_read(): Could not allocate seq.\n");
	    exit(EX_UNAVAILABLE);
	}
	p = record->seq;
	do
	{
	    if ( ch != '\n' )
		*p++ = ch;
	    if ( p - record->seq == record->seq_array_size )
	    {
		record->seq_array_size *= 2;
		record->seq = xt_realloc(record->seq, record->seq_array_size,
		    sizeof(*record->seq));
		if ( record->seq == NULL )
		{
		    fprintf(stderr, "bl_fasta_read(): Could not reallocate seq.\n");
		    exit(EX_UNAVAILABLE);
		}
	    }
	}   while ( ((ch = getc(fasta_stream)) != '>') && (ch != EOF) );
	*p = '\0';
	record->seq_len = p - record->seq;

	/* Trim array */
	record->seq = xt_realloc(record->seq, record->seq_len+1,
	    sizeof(*record->desc));
	//printf("Trimmed sequence to %zu\n", record->seq_len+1);

	if ( ch == '>' )
	    ungetc(ch, fasta_stream);
	return BL_READ_OK;
    }
    else
	return BL_READ_BAD_DATA;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc
 *
 *  Description:
 *      Write a FASTA record to the specified FILE stream, writing at most
 *      max_line_len sequence characters per line.  The special value
 *      BL_FASTA_LINE_UNLIMITED indicates no line length limit.
 *  
 *  Arguments:
 *      fasta_stream:   FILE stream to which data are written
 *      record:         Pointer to a bl_fasta_t structure to be written
 *      max_line_len:   Maximum length of a sequence line in output
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE if a write error occurs.
 *
 *  Examples:
 *      while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fasta_write(stdout, &rec, BL_FASTA_LINE_UNLIMITED);
 *
 *  See also:
 *      bl_fasta_read(3), bl_fastq_read(3), bl_fastq_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fasta_write(FILE *fasta_stream, bl_fasta_t *record,
		       size_t max_line_len)

{
    size_t  c;
    int     save_ch;
    
    if ( fprintf(fasta_stream, ">%s\n", record->desc) < 0 )
	return BL_WRITE_FAILURE;
    
    for (c = 0; c < record->seq_len; c += max_line_len)
    {
	if ( record->seq_len - c > max_line_len )
	{
	    save_ch = record->seq[c + max_line_len];
	    record->seq[c + max_line_len] = '\0';
	}
	if ( fprintf(fasta_stream, "%s\n", record->seq + c) < 0 )
	    return BL_WRITE_FAILURE;
	if ( record->seq_len - c > max_line_len )
	    record->seq[c + max_line_len] = save_ch;
    }
    return BL_WRITE_OK;
}
