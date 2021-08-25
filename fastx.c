#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include "fastx.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a FASTA or FASTQ record from a FILE stream by calling
 *      bl_read_fasta(3) or bl_read_fastq(3).  The bl_fastx_t structure
 *      must first be initialized by assigning it BL_FASTX_INIT and
 *      calling bl_fastx_init(3).
 *      See bl_fasta_read(3) and bl_fastq_read(3) for further details.
 *
 *  Arguments:
 *      fastx_stream    FILE stream from which FASTA data are read
 *      record          Pointer to a bl_fastx_t structure to receive data
 *
 *  Returns:
 *      BL_READ_OK upon successful read of description and sequence
 *      BL_READ_BAD_DATA if something is amiss with input format
 *      BL_READ_EOF if no more data are available
 *
 *  Examples:
 *      bl_fastx_t  rec = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &rec);
 *      while ( bl_fastx_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastx_write(stdout, &rec, BL_FASTX_LINE_UNLIMITED);
 *      bl_fastx_free(&rec);
 *
 *  See also:
 *      bl_fastx_write(3), bl_fastq_read(3), bl_fastq_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastx_read(FILE *fastx_stream, bl_fastx_t *record)

{
    switch(BL_FASTX_FORMAT(record))
    {
	case BL_FASTX_FORMAT_FASTA:
	    record->format = BL_FASTX_FORMAT_FASTA;
	    return bl_fasta_read(fastx_stream, &record->fasta);
	case BL_FASTX_FORMAT_FASTQ:
	    record->format = BL_FASTX_FORMAT_FASTQ;
	    return bl_fastq_read(fastx_stream, &record->fastq);
    }
    fprintf(stderr, "bl_fastx_read(): Input format is unknown.  Call bl_fastx_init() first.\n");
    return BL_READ_UNKNOWN_FORMAT;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a FASTA or FASTQ record from a FILE stream by calling
 *      bl_fasta_write(3) or bl_fastq_write(3).  The bl_fastx_t structure
 *      must first be initialized by assigning it BL_FASTX_INIT and
 *      calling bl_fastx_init(3), and then populated by bl_fastx_read(3)
 *      or other means.
 *      See bl_fasta_write(3) and bl_fastq_write(3) for further details.
 *  
 *  Arguments:
 *      fastx_stream    FILE stream to which data are written
 *      record          Pointer to a bl_fastx_t structure to be written
 *      max_line_len    Maximum length of a sequence line in output
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE if a write error occurs.
 *
 *  Examples:
 *      bl_fastx_t  rec;
 *
 *      bl_fastx_init(stdin, &rec);
 *      while ( bl_fastx_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastx_write(stdout, &rec, BL_FASTX_LINE_UNLIMITED);
 *      bl_fastx_free(&rec);
 *
 *  See also:
 *      bl_fastx_read(3), bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastx_write(FILE *fastx_stream, bl_fastx_t *record,
		       size_t max_line_len)

{
    switch(BL_FASTX_FORMAT(record))
    {
	case BL_FASTX_FORMAT_FASTA:
	    return bl_fasta_write(fastx_stream, &record->fasta, max_line_len);
	    break;
	case BL_FASTX_FORMAT_FASTQ:
	    return bl_fastq_write(fastx_stream, &record->fastq, max_line_len);
	    break;
    }
    fprintf(stderr, "bl_fasta_write(): File format is unknown.\n");
    return BL_WRITE_FAILURE;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fast.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_fastx_read()
 *  
 *  Arguments:
 *      record  Pointer to a previously populated bl_fastx_t structure
 *
 *  Examples:
 *      bl_fastx_t  rec = BL_FASTX_INIT;
 *
 *      while ( bl_fastx_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastx_write(stdout, &rec, BL_FASTX_LINE_UNLIMITED);
 *      bl_fastx_free(&rec);
 *
 *  See also:
 *      bl_fastx_read(3), bl_fastx_write(3)
 *      bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastx_free(bl_fastx_t *record)

{
    switch(BL_FASTX_FORMAT(record))
    {
	case BL_FASTX_FORMAT_FASTA:
	    bl_fasta_free(&record->fasta);
	    break;
	case BL_FASTX_FORMAT_FASTQ:
	    bl_fastq_free(&record->fastq);
	    break;
    }
    record->format = BL_FASTX_FORMAT_UNKNOWN;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_fastx_t structure.  This must be done before
 *      passing it to bl_fastx_read() for the first time, so that
 *      bl_fastx_read() will know to allocate memory for the fields.
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure to initialize.
 *
 *  Examples:
 *      bl_fastx_t  rec;
 *
 *      bl_fastx_init(&rec);
 *      bl_fastx_read(stdin, &rec);
 *
 *  See also:
 *      bl_fastx_read(3), bl_fastx_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastx_init(FILE *fastx_stream, bl_fastx_t *record)

{
    int     ch;
    
    ch = getc(fastx_stream);
    ungetc(ch, fastx_stream);
    switch(ch)
    {
	case '>':
	    record->format = BL_FASTX_FORMAT_FASTA;
	    bl_fasta_init(&record->fasta);
	    break;
	case '@':
	    record->format = BL_FASTX_FORMAT_FASTQ;
	    bl_fastq_init(&record->fastq);
	    break;
    }
}
