#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <inttypes.h>   // PRIu64
#include "gffio.h"
#include "dsvio.h"

/***************************************************************************
 *  Description:
 *      Skip over header lines in gff input stream.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-08  Jason Bacon Begin
 ***************************************************************************/

FILE    *gff_skip_header(FILE *gff_stream)

{
    int     ch;
    FILE    *header_stream = tmpfile();

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like peak-classifier to replicate the
     *  header in output files.
     */
    
    while ( (ch = getc(gff_stream)) == '#' )
    {
	putc(ch, header_stream);
	do
	{
	    ch = getc(gff_stream);
	    putc(ch, header_stream);
	}   while ( (ch != '\n') && (ch != EOF) );
    }
    
    // Rewind to start of first non-header line
    if ( ch != EOF )
	ungetc(ch, gff_stream);
    return header_stream;
}


/***************************************************************************
 *  Description:
 *      Read fields from one line of a GFF file.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-08  Jason Bacon Begin
 ***************************************************************************/

int     gff_read_feature(FILE *gff_stream, gff_feature_t *gff_feature)

{
    char    *end;
    size_t  len;
    int     delim;
    
    // Chromosome
    if ( tsv_read_field(gff_stream, gff_feature->sequence,
			BIO_CHROMOSOME_MAX_CHARS, &len) == EOF )
    {
	fputs("gff_read_feature(): Info: Got EOF reading SEQUENCE, as expected.\n", stderr);
	return BIO_READ_EOF;
    }
    
    // Source
    if ( (delim = tsv_read_field(gff_stream, gff_feature->source,
			GFF_NAME_MAX_CHARS, &len)) == EOF )
    {
	fprintf(stderr, "gff_read_feature(): Got EOF reading SOURCE: %s.\n",
		gff_feature->source);
	return BIO_READ_TRUNCATED;
    }

    // Feature
    if ( (delim = tsv_read_field(gff_stream, gff_feature->feature,
			GFF_NAME_MAX_CHARS, &len)) == EOF )
    {
	fprintf(stderr, "gff_read_feature(): Got EOF reading feature: %s.\n",
		gff_feature->feature);
	return BIO_READ_TRUNCATED;
    }
    
    // Feature start position
    if ( tsv_read_field(gff_stream, gff_feature->start_pos_str,
			BIO_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "gff_read_feature(): Got EOF reading start POS: %s.\n",
		gff_feature->start_pos_str);
	return BIO_READ_TRUNCATED;
    }
    else
    {
	gff_feature->start_pos = strtoul(gff_feature->start_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "gff_read_feature(): Invalid feature position: %s\n",
		    gff_feature->start_pos_str);
	    return BIO_READ_TRUNCATED;
	}
    }
    
    // Feature end position
    if ( (delim = tsv_read_field(gff_stream, gff_feature->end_pos_str,
			BIO_POSITION_MAX_DIGITS, &len)) == EOF )
    {
	fprintf(stderr, "gff_read_feature(): Got EOF reading end POS: %s.\n",
		gff_feature->end_pos_str);
	return BIO_READ_TRUNCATED;
    }
    else
    {
	gff_feature->end_pos = strtoul(gff_feature->end_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "gff_read_feature(): Invalid feature position: %s\n",
		    gff_feature->end_pos_str);
	    return BIO_READ_TRUNCATED;
	}
    }

    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(gff_stream, gff_feature->score_str,
			    BIO_POSITION_MAX_DIGITS, &len)) == EOF )
	{
	    fprintf(stderr, "gff_read_feature(): Got EOF reading end SCORE: %s.\n",
		    gff_feature->score_str);
	    return BIO_READ_TRUNCATED;
	}
	else
	{
	    gff_feature->score = strtoul(gff_feature->score_str, &end, 10);
	    if ( (*end != '\0') || (gff_feature->score > 1000) )
	    {
		fprintf(stderr,
			"gff_read_feature(): Invalid feature score: %s\n",
			gff_feature->score_str);
		return BIO_READ_TRUNCATED;
	    }
	}
    }
    
    if ( delim != '\n' )
	dsv_skip_rest_of_line(gff_stream);
    return BIO_READ_OK;
}


/***************************************************************************
 *  Description:
 *      Write fields from one line of a gff file.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-08  Jason Bacon Begin
 ***************************************************************************/

int     gff_write_feature(FILE *gff_stream, gff_feature_t *gff_feature,
				gff_field_mask_t field_mask)

{
    return fprintf(gff_stream, "%s\t%s\t%s\t%" PRIu64 "\t%" PRIu64 "\t\n%u",
	    gff_feature->sequence, gff_feature->source, gff_feature->feature,
	    gff_feature->start_pos, gff_feature->end_pos, gff_feature->score);
}
