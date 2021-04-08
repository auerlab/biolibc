#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include "bedio.h"
#include "dsvio.h"

/***************************************************************************
 *  Description:
 *      Skip over header lines in bed input stream.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-05  Jason Bacon Begin
 ***************************************************************************/

FILE    *bed_skip_header(FILE *bed_stream)

{
    char    start[7] = "xxxxxx";
    size_t  count;
    int     ch;
    FILE    *header_stream = tmpfile();

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like peak-classifier to replicate the
     *  header in output files.
     */
    
    while ( ((count=fread(start, 6, 1, bed_stream)) == 1) && 
	    ((memcmp(start, "browser", 6) == 0) ||
	    (memcmp(start, "track", 5) == 0) ||
	    (*start == '#')) )
    {
	fwrite(start, 6, 1, header_stream);
	do
	{
	    ch = getc(bed_stream);
	    putc(ch, header_stream);
	}   while ( ch != '\n' );
    }
    return header_stream;
}


/***************************************************************************
 *  Description:
 *      Read static fields from one line of a single-entry BED file.
 *      Does not read sample data.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bed_read_feature(FILE *bed_stream, bed_feature_t *bed_feature)

{
    char    *end;
    size_t  len;
    
    // Chromosome
    if ( tsv_read_field(bed_stream, bed_feature->chromosome,
			BED_CHROMOSOME_MAX_CHARS, &len) == EOF )
    {
	fputs("bed_read_feature(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return BIO_READ_EOF;
    }
    
    // Feature position
    if ( tsv_read_field(bed_stream, bed_feature->start_pos_str,
			BED_POSITION_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "bed_read_feature(): Got EOF reading start POS: %s.\n",
		bed_feature->start_pos_str);
	return BIO_READ_TRUNCATED;
    }
    else
    {
	bed_feature->start_pos = strtoul(bed_feature->start_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bed_read_feature(): Invalid feature position: %s\n",
		    bed_feature->start_pos_str);
	    return BIO_READ_TRUNCATED;
	}
    }
    
    // Feature position
    if ( tsv_read_field(bed_stream, bed_feature->end_pos_str,
			BED_POSITION_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "bed_read_feature(): Got EOF reading end POS: %s.\n",
		bed_feature->end_pos_str);
	return BIO_READ_TRUNCATED;
    }
    else
    {
	bed_feature->end_pos = strtoul(bed_feature->end_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bed_read_feature(): Invalid feature position: %s\n",
		    bed_feature->end_pos_str);
	    return BIO_READ_TRUNCATED;
	}
    }
    
    return BIO_READ_OK;
}


/***************************************************************************
 *  Description:
 *      Write static fields from one line of a single-entry bed file.
 *      Does not write sample data.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bed_write_FEATURE(FILE *bed_stream, bed_feature_t *bed_feature,
				bed_field_mask_t field_mask)

{
    return 0;
}
