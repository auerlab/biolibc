#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <inttypes.h>   // PRIu64
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
    int     ch, c;
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
	}   while ( (ch != '\n') && (ch != EOF) );
    }
    
    // Rewind to start of first non-header line
    if ( count == 1 )
	for (c = 5; c >= 0; --c)
	    ungetc(start[c], bed_stream);
    return header_stream;
}


/***************************************************************************
 *  Description:
 *      Read fields from one line of a BED file.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bed_read_feature(FILE *bed_stream, bed_feature_t *bed_feature)

{
    char    *end;
    size_t  len;
    int     delim;
    
    // Chromosome
    if ( tsv_read_field(bed_stream, bed_feature->chromosome,
			BIO_CHROMOSOME_MAX_CHARS, &len) == EOF )
    {
	fputs("bed_read_feature(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return BIO_READ_EOF;
    }
    
    // Feature start position
    if ( tsv_read_field(bed_stream, bed_feature->start_pos_str,
			BIO_POSITION_MAX_DIGITS, &len) == EOF )
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
    
    // Feature end position
    if ( (delim = tsv_read_field(bed_stream, bed_feature->end_pos_str,
			BIO_POSITION_MAX_DIGITS, &len)) == EOF )
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

    bed_feature->fields = 3;
    
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, bed_feature->name,
			    BED_NAME_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading end NAME: %s.\n",
		    bed_feature->name);
	    return BIO_READ_TRUNCATED;
	}
	++bed_feature->fields;
    }
    
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, bed_feature->score_str,
			    BIO_POSITION_MAX_DIGITS, &len)) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading end SCORE: %s.\n",
		    bed_feature->score_str);
	    return BIO_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->score = strtoul(bed_feature->score_str, &end, 10);
	    if ( (*end != '\0') || (bed_feature->score > 1000) )
	    {
		fprintf(stderr,
			"bed_read_feature(): Invalid feature score: %s\n",
			bed_feature->score_str);
		return BIO_READ_TRUNCATED;
	    }
	}
	++bed_feature->fields;
    }
    
    if ( delim != '\n' )
	dsv_skip_rest_of_line(bed_stream);
    return BIO_READ_OK;
}


/***************************************************************************
 *  Description:
 *      Write fields from one line of a bed file.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bed_write_feature(FILE *bed_stream, bed_feature_t *bed_feature,
				bed_field_mask_t field_mask)

{
    fprintf(bed_stream, "%s\t%" PRIu64 "\t%" PRIu64,
	    bed_feature->chromosome,
	    bed_feature->start_pos, bed_feature->end_pos);
    if ( bed_feature->fields > 3 )
	fprintf(bed_stream, "\t%s", bed_feature->name);
    if ( bed_feature->fields > 4 )
	fprintf(bed_stream, "\t%u", bed_feature->score);
    putc('\n', bed_stream);
    return 0;
}
