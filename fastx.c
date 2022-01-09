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

int     bl_fastx_read(bl_fastx_t *record, FILE *fastx_stream)

{
    switch(BL_FASTX_FORMAT(record))
    {
	case BL_FASTX_FORMAT_FASTA:
	    record->format = BL_FASTX_FORMAT_FASTA;
	    return bl_fasta_read(&record->fasta, fastx_stream);
	case BL_FASTX_FORMAT_FASTQ:
	    record->format = BL_FASTX_FORMAT_FASTQ;
	    return bl_fastq_read(&record->fastq, fastx_stream);
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
 *      or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
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
 *      bl_fastx_t  rec = BL_FASTX_INIT;
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

int     bl_fastx_write(bl_fastx_t *record, FILE *fastx_stream,
	    size_t max_line_len)

{
    switch(BL_FASTX_FORMAT(record))
    {
	case BL_FASTX_FORMAT_FASTA:
	    return bl_fasta_write(&record->fasta, fastx_stream, max_line_len);
	    break;
	case BL_FASTX_FORMAT_FASTQ:
	    return bl_fastq_write(&record->fastq, fastx_stream, max_line_len);
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
 *      Initialize a bl_fastx_t structure by peaking at the first character
 *      of the description string to determine whether the stream is FASTA
 *      or FASTQ, and then initializing the appropriate structure within
 *      the bl_fastx_t structure.  This must be done before
 *      passing it to bl_fastx_read() for the first time, so that
 *      bl_fastx_read() will know to allocate memory for the fields.
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure to initialize.
 *
 *  Examples:
 *      bl_fastx_t  rec = BL_FASTX_INIT;
 *
 *      bl_fastx_init(&rec);
 *      bl_fastx_read(stdin, &rec);
 *      ...
 *      bl_fastx_free(&rec);
 *
 *  See also:
 *      bl_fastx_read(3), bl_fastx_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastx_init(bl_fastx_t *record, FILE *fastx_stream)

{
    int     ch;
    
    if ( record->format != BL_FASTX_FORMAT_UNKNOWN )
    {
	fputs("bl_fastx_init(): Warning: format should be unknown.\n"
	      "bl_fastx_t variables should be initialized with BL_FASTX_INIT.\n"
	      "This could also indicate a previously used structure that has\n"
	      "not been freed, which is a memory leak.\n", stderr);
    }
    
    /* Skip comment lines */
    while ( ((ch = getc(fastx_stream)) == ';') && (ch != EOF) )
	while ( ((ch = getc(fastx_stream)) != '\n') && (ch != EOF) )
	    ;
    if ( ch == EOF )
    {
	fputs("bl_fastq_init(): EOF encountered.\n", stderr);
	exit(EX_DATAERR);
    }
    ungetc(ch, fastx_stream);
    switch(ch)
    {
	case '>':
	    fputs("File format is FASTA.\n", stderr);
	    record->format = BL_FASTX_FORMAT_FASTA;
	    bl_fasta_init(&record->fasta);
	    break;
	case '@':
	    fputs("File format is FASTQ.\n", stderr);
	    record->format = BL_FASTX_FORMAT_FASTQ;
	    bl_fastq_init(&record->fastq);
	    break;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return a pointer to the description string of a FASTA or FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Pointer to the description string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Desc string = %s\n", bl_fastx_desc(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

char    *bl_fastx_desc(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    return BL_FASTA_DESC(&record->fasta);
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_DESC(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_desc(): File format is unknown.\n");
    return NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the length of the description string of a FASTA or FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Length of the description string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Desc string length = %zu\n", bl_fastx_desc_len(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastx_desc_len(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    return BL_FASTA_DESC_LEN(&record->fasta);
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_DESC_LEN(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_desc_len(): File format is unknown.\n");
    return 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return a pointer to the seq string of a FASTA or FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Pointer to the seq string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Seq string = %s\n", bl_fastx_seq(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

char    *bl_fastx_seq(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    return BL_FASTA_SEQ(&record->fasta);
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_SEQ(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_seq(): File format is unknown.\n");
    return NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the length of the seq string of a FASTA or FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Length of the seq string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Seq string length = %zu\n", bl_fastx_seq_len(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastx_seq_len(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    return BL_FASTA_SEQ_LEN(&record->fasta);
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_SEQ_LEN(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_seq_len(): File format is unknown.\n");
    return 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return a pointer to the + string of a FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  If the format
 *      if the fastx stream is FASTA, a warning is generated and NULL
 *      is returned.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Pointer to the + string, or NULL if FASTA
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("+ string = %s\n", bl_fastx_plus(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

char    *bl_fastx_plus(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    fputs("bl_fastx_plus(): Warning: Attempt to access + field in a FASTA stream.\n", stderr);
	    return NULL;
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_PLUS(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_plus(): File format is unknown.\n");
    return NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the length of the + string of a FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  If the format
 *      if the fastx stream is FASTA, a warning is generated and 0
 *      is returned.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Length of the + string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("+ string length = %zu\n", bl_fastx_plus_len(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastx_plus_len(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    fputs("bl_fastx_plus_len(): Warning: Attempt to access + length field in a FASTA stream.\n", stderr);
	    return 0;
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_PLUS_LEN(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_plus_len(): File format is unknown.\n");
    return 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return a pointer to the qual string of a FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  If the format
 *      if the fastx stream is FASTA, a warning is generated and NULL
 *      is returned.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Pointer to the qual string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Qual string = %s\n", bl_fastx_qual(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

char    *bl_fastx_qual(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    fputs("bl_fastx_qual(): Warning: Attempt to access + field in a FASTA stream.\n", stderr);
	    return NULL;
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_QUAL(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_qual(): File format is unknown.\n");
    return NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the length of the qual string of a FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  If the format
 *      if the fastx stream is FASTA, a warning is generated and 0
 *      is returned.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Length of the qual string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Qual string length = %zu\n", bl_fastx_qual_len(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastx_qual_len(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    fputs("bl_fastx_qual_len(): Warning: Attempt to access + length field in a FASTA stream.\n", stderr);
	    return 0;
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_QUAL_LEN(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_qual_len(): File format is unknown.\n");
    return 0;
}
