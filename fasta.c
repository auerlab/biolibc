#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <xtend/dsv.h>
#include <xtend/mem.h>
#include "fasta.h"
#include "biolibc.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a FASTA record from a FILE stream.  Each record must begin
 *      with a description line (beginning with '>'), which is then
 *      followed by one or more lines of sequence data.  The end of the
 *      sequence is marked either by the next description line or EOF.
 *      If desc_len and seq_len are 0 (e.g. the structure is initialized
 *      with BL_FASTA_INIT or bl_fasta_init(3), or has been freed with
 *      bl_fasta_free(3), then
 *      memory is allocated for the description and sequence.
 *
 *      Otherwise, the existing allocated buffers are reused.  Hence, when
 *      reading many FASTA records of the same length, only one allocation
 *      is needed.  In any case, the buffers are automatically enlarged if
 *      they become full and automatically trimmed to the actual data size
 *      after reading is complete.
 *
 *      Buffer memory should be freed as soon as possible by calling
 *      bl_fasta_free(3).
 *  
 *  Arguments:
 *      fasta_stream    FILE stream from which FASTA data are read
 *      record          Pointer to a bl_fasta_t structure to receive data
 *
 *  Returns:
 *      BL_READ_OK upon successful read of description and sequence
 *      BL_READ_BAD_DATA if something is amiss with input format
 *      BL_READ_EOF if no more data are available
 *
 *  Examples:
 *      bl_fasta_t  rec = BL_FASTA_INIT;
 *
 *      while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fasta_write(stdout, &rec, BL_FASTA_LINE_UNLIMITED);
 *      bl_fasta_free(&rec);
 *
 *  See also:
 *      bl_fasta_write(3), bl_fastq_read(3), bl_fastq_write(3),
 *      bl_fasta_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fasta_read(FILE *fasta_stream, bl_fasta_t *record)

{
    int     ch;
    size_t  len;
    
    /* Skip comment lines */
    while ( (ch = getc(fasta_stream)) == ';' )
	while ( getc(fasta_stream) != '\n' )
	    ;
    
    if ( ch == EOF )
	return BL_READ_EOF;
    
    /* Every record should begin with a '>' */
    if ( ch == '>' )    // Desc
    {
	ungetc(ch, fasta_stream);
	ch = dsv_read_field_malloc(fasta_stream, &record->desc,
			    &record->desc_array_size, "", &record->desc_len);
	if ( record->desc == NULL )
	{
	    fprintf(stderr, "bl_fasta_read(): Could not allocate desc.\n");
	    exit(EX_UNAVAILABLE);
	}
	
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
		fprintf(stderr, "bl_fasta_read(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	
	len = 0;
	do
	{
	    if ( ch != '\n' )
		record->seq[len++] = ch;
	    if ( len == record->seq_array_size )
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
	record->seq[len] = '\0';
	record->seq_len = len;

	/* Trim array */
	record->seq_array_size = record->seq_len + 1;
	record->seq = xt_realloc(record->seq, record->seq_array_size,
	    sizeof(*record->desc));

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
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a FASTA record to the specified FILE stream, writing at most
 *      max_line_len sequence characters per line.  The special value
 *      BL_FASTA_LINE_UNLIMITED indicates no line length limit.
 *  
 *  Arguments:
 *      fasta_stream    FILE stream to which data are written
 *      record          Pointer to a bl_fasta_t structure to be written
 *      max_line_len    Maximum length of a sequence line in output
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE if a write error occurs.
 *
 *  Examples:
 *      bl_fasta_t  rec = BL_FASTA_INIT;
 *
 *      while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fasta_write(stdout, &rec, BL_FASTA_LINE_UNLIMITED);
 *      bl_fasta_free(&rec);
 *
 *  See also:
 *      bl_fasta_read(3), bl_fastq_read(3), bl_fastq_write(3),
 *
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
    
    if ( fprintf(fasta_stream, "%s\n", record->desc) < 0 )
	return BL_WRITE_FAILURE;
    
    for (c = 0; c < record->seq_len; c += max_line_len)
    {
	// Temporarily null-terminate segment of string to be printed
	if ( record->seq_len - c > max_line_len )
	{
	    save_ch = record->seq[c + max_line_len];
	    record->seq[c + max_line_len] = '\0';
	}
	
	// Print segment
	if ( fprintf(fasta_stream, "%s\n", record->seq + c) < 0 )
	    return BL_WRITE_FAILURE;
	
	// Remove temporary null-termination
	if ( record->seq_len - c > max_line_len )
	    record->seq[c + max_line_len] = save_ch;
    }
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fast.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_fasta_read()
 *  
 *  Arguments:
 *      record  Pointer to a previously populated bl_fasta_t structure
 *
 *  Examples:
 *      bl_fasta_t  rec = BL_FASTA_INIT;
 *
 *      while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fasta_write(stdout, &rec, BL_FASTA_LINE_UNLIMITED);
 *      bl_fasta_free(&rec);
 *
 *  See also:
 *      bl_fasta_read(3), bl_fasta_write(3)
 *      bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_fasta_free(bl_fasta_t *record)

{
    free(record->seq);
    free(record->desc);
    record->desc = record->seq = NULL;
    record->desc_array_size = record->seq_array_size = 0;
    record->desc_len = record->seq_len = 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_fasta_t structure.  This must be done before
 *      passing it to bl_fasta_read() for the first time, so that
 *      bl_fasta_read() will know to allocate memory for the fields.
 *  
 *  Arguments:
 *      record  Pointer to the bl_fasta_t structure to initialize.
 *
 *  Examples:
 *      bl_fasta_t  rec;
 *
 *      bl_fasta_init(&rec);
 *      bl_fasta_read(stdin, &rec);
 *
 *  See also:
 *      bl_fasta_read(3), bl_fasta_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_fasta_init(bl_fasta_t *record)

{
    record->desc = record->seq = NULL;
    record->desc_array_size = record->seq_array_size = 0;
    record->desc_len = record->seq_len = 0;
}
