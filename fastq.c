
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <xtend.h>
#include "fastq.h"
#include "biolibc.h"


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a FASTQ record from a FILE stream.  Each record must begin
 *      with a description line (beginning with '@'), which is then
 *      followed by one or more lines of sequence data, a separator line
 *      beginning with '+', and a line of quality scores.  The end of the
 *      sequence is marked either by the next description line or EOF.
 *      If desc_len and seq_len are 0 (e.g. the structure is initialized
 *      with BL_FASTQ_INIT or has been freed with bl_fastq_free(3), then
 *      memory is allocated for each line.
 *
 *      Otherwise, the existing allocated buffers are reused.  Hence, when
 *      reading many FASTQ records of the same length, only one allocation
 *      is needed.  In any case, the buffers are automatically enlarged if
 *      they become full and automatically trimmed to the actual data size
 *      after reading is complete.
 *
 *      Buffer memory should be freed as soon as possible by calling
 *      bl_fastq_free(3).
 *  
 *  Arguments:
 *      fastq_stream    FILE stream from which FASTQ data are read
 *      record          Pointer to a bl_fastq_t structure to receive data
 *
 *  Returns:
 *      BL_READ_OK upon successful read of description and sequence
 *      BL_READ_BAD_DATA if something is amiss with input format
 *      BL_READ_EOF if no more data are available
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_write(3), bl_fastq_read(3), bl_fastq_write(3),
 *      bl_fastq_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastq_read(FILE *fastq_stream, bl_fastq_t *record)

{
    int     ch;
    char    *p;
    
    /* Skip comment lines */
    while ( (ch = getc(fastq_stream)) == ';' )
	while ( getc(fastq_stream) != '\n' )
	    ;
    
    if ( ch == EOF )
	return BL_READ_EOF;
    
    /* Every record should begin with a '@' */
    if ( ch == '@' )    // Desc
    {
	/*
	 *  Read description
	 */
	
	fprintf(stderr, "New record...\n");
	if ( record->desc_array_size == 0 )
	{
	    record->desc_array_size = 1024;
	    record->desc = xt_malloc(record->desc_array_size, sizeof(*record->desc));
	    if ( record->desc == NULL )
	    {
		fprintf(stderr, "bl_fastq_read(): Could not allocate desc.\n");
		exit(EX_UNAVAILABLE);
	    }
	}

	p = record->desc;
	while ( ((ch = getc(fastq_stream)) != '\n') && (ch != EOF) )
	{
	    *p++ = ch;
	    if ( p - record->desc == record->desc_array_size )
	    {
		record->desc_array_size *= 2;
		record->desc = xt_realloc(record->desc, record->desc_array_size,
		    sizeof(*record->desc));
		if ( record->desc == NULL )
		{
		    fprintf(stderr, "bl_fastq_read(): Could not reallocate desc.\n");
		    exit(EX_UNAVAILABLE);
		}
	    }
	}
	*p = '\0';
	record->desc_len = p - record->desc;

	/* Trim array */
	record->desc_array_size = record->desc_len + 1;
	record->desc = xt_realloc(record->desc, record->desc_array_size,
	    sizeof(*record->desc));
	fprintf(stderr, "desc = %s\n", record->desc);
	
	/* Should not encounter EOF while reading description line */
	/* Every description should be followed by at least one seq line */
	if ( ch == EOF )
	    return BL_READ_TRUNCATED;

	/*
	 *  Read sequence lines
	 */
	
	if ( record->seq_array_size == 0 )
	{
	    record->seq_array_size = 1024;
	    record->seq = xt_malloc(record->seq_array_size, sizeof(*record->seq));
	    if ( record->seq == NULL )
	    {
		fprintf(stderr, "bl_fastq_read(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
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
		    fprintf(stderr, "bl_fastq_read(): Could not reallocate seq.\n");
		    exit(EX_UNAVAILABLE);
		}
	    }
	}   while ( ((ch = getc(fastq_stream)) != '+') && (ch != EOF) );
	*p = '\0';
	record->seq_len = p - record->seq;

	/* Trim array */
	record->seq_array_size = record->seq_len + 1;
	record->seq = xt_realloc(record->seq, record->seq_array_size,
	    sizeof(*record->desc));
	fprintf(stderr, "seq = %s\n", record->seq);

	/* Should not encounter EOF while reading sequence lines */
	/* Every sequence should be followed by a + separator line */
	if ( ch == EOF )
	    return BL_READ_TRUNCATED;
	    
	/*
	 *  Read + separator
	 */
	
	if ( ch != '+' )
	    return BL_READ_BAD_DATA;
	
	if ( record->plus_array_size == 0 )
	{
	    record->plus_array_size = 1024;
	    record->plus = xt_malloc(record->plus_array_size, sizeof(*record->plus));
	    if ( record->plus == NULL )
	    {
		fprintf(stderr, "bl_fastq_read(): Could not allocate plus.\n");
		exit(EX_UNAVAILABLE);
	    }
	}

	p = record->plus;
	while ( ((ch = getc(fastq_stream)) != '\n') && (ch != EOF) )
	{
	    *p++ = ch;
	    if ( p - record->plus == record->plus_array_size )
	    {
		record->plus_array_size *= 2;
		record->plus = xt_realloc(record->plus, record->plus_array_size,
		    sizeof(*record->plus));
		if ( record->plus == NULL )
		{
		    fprintf(stderr, "bl_fastq_read(): Could not reallocate plus.\n");
		    exit(EX_UNAVAILABLE);
		}
	    }
	}
	*p = '\0';
	record->plus_len = p - record->plus;

	/* Trim array */
	record->plus_array_size = record->plus_len + 1;
	record->plus = xt_realloc(record->plus, record->plus_array_size,
	    sizeof(*record->plus));
	fprintf(stderr, "plus = %s\n", record->plus);
	
	/* Should not encounter EOF while reading plus line */
	/* Every plus should be followed by at least one qual line */
	if ( ch == EOF )
	    return BL_READ_TRUNCATED;

	/*
	 *  Read quality string
	 */
	
	if ( record->qual_array_size == 0 )
	{
	    // Must match sequence len and no need to trim
	    record->qual_array_size = record->seq_array_size;
	    record->qual = xt_malloc(record->qual_array_size, sizeof(*record->qual));
	    if ( record->qual == NULL )
	    {
		fprintf(stderr, "bl_fastq_read(): Could not allocate qual.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	p = record->qual;
	do
	{
	    /* Read at least one full line, since '@' can be a quality score */
	    while ( ((ch = getc(fastq_stream)) != '\n') && (ch != EOF) )
	    {
		*p++ = ch;
		if ( p - record->qual == record->qual_array_size )
		{
		    record->qual_array_size *= 2;
		    record->qual = xt_realloc(record->qual, record->qual_array_size,
			sizeof(*record->qual));
		    if ( record->qual == NULL )
		    {
			fprintf(stderr, "bl_fastq_read(): Could not reallocate qual.\n");
			exit(EX_UNAVAILABLE);
		    }
		}
	    }
	    if ( ch == EOF )
		return BL_READ_BAD_DATA;
	    
	}   while ( ((ch = getc(fastq_stream)) != '@') && (ch != EOF) );
	*p = '\0';
	record->qual_len = p - record->qual;
	fprintf(stderr, "qual = %s\n", record->qual);
	// No need to trim since qual must be the same size as seq

	if ( ch == '@' )
	    ungetc(ch, fastq_stream);
	return BL_READ_OK;
    }
    else
	return BL_READ_BAD_DATA;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a FASTQ record to the specified FILE stream, writing at most
 *      max_line_len sequence characters per line.  The special value
 *      BL_FASTQ_LINE_UNLIMITED indicates no line length limit.
 *  
 *  Arguments:
 *      fastq_stream    FILE stream to which data are written
 *      record          Pointer to a bl_fastq_t structure to be written
 *      max_line_len    Maximum length of a sequence line in output
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE if a write error occurs.
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_read(3), bl_fastq_read(3), bl_fastq_write(3),
 *
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastq_write(FILE *fastq_stream, bl_fastq_t *record,
		       size_t max_line_len)

{
    size_t  c;
    int     save_ch;
    
    if ( fprintf(fastq_stream, "@%s\n", record->desc) < 0 )
	return BL_WRITE_FAILURE;
    
    for (c = 0; c < record->seq_len; c += max_line_len)
    {
	if ( record->seq_len - c > max_line_len )
	{
	    save_ch = record->seq[c + max_line_len];
	    record->seq[c + max_line_len] = '\0';
	}
	if ( fprintf(fastq_stream, "%s\n", record->seq + c) < 0 )
	    return BL_WRITE_FAILURE;
	if ( record->seq_len - c > max_line_len )
	    record->seq[c + max_line_len] = save_ch;
    }
    
    if ( fprintf(fastq_stream, "+%s\n", record->plus) < 0 )
	return BL_WRITE_FAILURE;
    
    for (c = 0; c < record->qual_len; c += max_line_len)
    {
	if ( record->qual_len - c > max_line_len )
	{
	    save_ch = record->qual[c + max_line_len];
	    record->qual[c + max_line_len] = '\0';
	}
	if ( fprintf(fastq_stream, "%s\n", record->qual + c) < 0 )
	    return BL_WRITE_FAILURE;
	if ( record->qual_len - c > max_line_len )
	    record->qual[c + max_line_len] = save_ch;
    }
    
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fast.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_fastq_read()
 *  
 *  Arguments:
 *      record  Pointer to a previously populated bl_fastq_t structure
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_read(3), bl_fastq_write(3)
 *      bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastq_free(bl_fastq_t *record)

{
    free(record->seq);
    free(record->desc);
    free(record->plus);
    free(record->qual);
    record->desc = record->seq = record->plus = record->qual = NULL;
    record->desc_array_size = record->seq_array_size = 
	record->plus_array_size = record->qual_array_size = 0;
    record->desc_len = record->seq_len = 
	record->plus_len = record->qual_len = 0;
}
