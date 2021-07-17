#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <inttypes.h>   // PRIu64
#include <xtend.h>
#include "gff.h"
#include "bed.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc
 *
 *  Description:
 *      Skip over header lines in gff input stream.  The FILE pointer
 *      gff_stream is advanced to the first character of the first line
 *      after the header.  The header is copied to a temporary file and and
 *      the function returns a FILE pointer to the header stream.
 *
 *  Arguments:
 *      gff_stream: FILE pointer to the open GFF file
 *
 *  Returns:
 *      A FILE pointer to a temporary file containing a copy of the header
 *
 *  See also:
 *      bl_gff_read_feature(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-08  Jason Bacon Begin
 ***************************************************************************/

FILE    *bl_gff_skip_header(FILE *gff_stream)

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
    rewind(header_stream);
    return header_stream;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read next feature (line) from a GFF file.
 *
 *      If field_mask is not BL_GFF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in gff_feature.
 *      That field in the structure is then populated with an appropriate
 *      marker, such as '.'.  Possible mask values are:
 *
 *      BL_GFF_FIELD_ALL
 *      BL_GFF_FIELD_SEQUENCE
 *      BL_GFF_FIELD_SOURCE
 *      BL_GFF_FIELD_NAME
 *      BL_GFF_FIELD_START_POS
 *      BL_GFF_FIELD_END_POS
 *      BL_GFF_FIELD_SCORE
 *      BL_GFF_FIELD_STRAND
 *      BL_GFF_FIELD_PHASE
 *      BL_GFF_FIELD_ATTRIBUTES
 *
 *  Arguments:
 *      gff_stream:     A FILE stream from which to read the line
 *      gff_feature:    Pointer to a bl_gff_t structure
 *      field_mask:     Bit mask indicating which fields to store in gff_feature
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered after a complete feature
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      bl_gff_read_feature(stdin, &gff_feature, BL_GFF_FIELD_ALL);
 *      bl_gff_read_feature(gff_stream, &gff_feature,
 *                       BL_GFF_FIELD_SEQUENCE|BL_GFF_FIELD_START_POS|BL_GFF_FIELD_END_POS);
 *
 *  See also:
 *      bl_gff_write_feature(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_read_feature(FILE *gff_stream, bl_gff_t *gff_feature,
			 gff_field_mask_t field_mask)

{
    char    *end,
	    line[BL_GFF_LINE_MAX_CHARS + 1],
	    temp_attributes[BL_GFF_ATTRIBUTES_MAX_CHARS + 1],
	    strand_str[BL_GFF_STRAND_MAX_CHARS + 1],
	    phase_str[BL_GFF_PHASE_MAX_DIGITS + 1],
	    *id_start,
	    *id_end,
	    start_pos_str[BL_POSITION_MAX_DIGITS + 1],
	    end_pos_str[BL_POSITION_MAX_DIGITS + 1],
	    score_str[BL_GFF_SCORE_MAX_DIGITS + 1];
    size_t  len;
    int     delim,
	    ch;
    
    // Check for group terminators (Line with just ###)
    // FIXME: Rely on parent ID instead of ###?
    if ( (ch = getc(gff_stream)) == '#' )
    {
	fgets(line, BL_GFF_LINE_MAX_CHARS, gff_stream);
	if ( strcmp(line, "##\n") == 0 )
	{
	    strlcpy(gff_feature->name, "###", BL_GFF_NAME_MAX_CHARS);
	    return BL_READ_OK;
	}
    }
    else if ( ch != EOF )
	ungetc(ch, gff_stream);
    
    // 1 Chromosome
    if ( tsv_read_field(gff_stream, gff_feature->sequence,
			BL_CHROMOSOME_MAX_CHARS, &len) == EOF )
    {
	//fputs("bl_gff_read_feature(): Info: Got EOF reading SEQUENCE, as expected.\n", stderr);
	return BL_READ_EOF;
    }
    
    // 2 Source
    if ( (delim = tsv_read_field(gff_stream, gff_feature->source,
			BL_GFF_NAME_MAX_CHARS, &len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read_feature(): Got EOF reading SOURCE: %s.\n",
		gff_feature->source);
	return BL_READ_TRUNCATED;
    }

    // 3 Feature
    if ( (delim = tsv_read_field(gff_stream, gff_feature->name,
			BL_GFF_NAME_MAX_CHARS, &len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read_feature(): Got EOF reading feature: %s.\n",
		gff_feature->name);
	return BL_READ_TRUNCATED;
    }
    
    // 4 Feature start position
    if ( tsv_read_field(gff_stream, start_pos_str,
			BL_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read_feature(): Got EOF reading start POS: %s.\n",
		start_pos_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	gff_feature->start_pos = strtoul(start_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_gff_read_feature(): Invalid feature position: %s\n",
		    start_pos_str);
	    return BL_READ_TRUNCATED;
	}
    }
    
    // 5 Feature end position
    if ( (delim = tsv_read_field(gff_stream, end_pos_str,
			BL_POSITION_MAX_DIGITS, &len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read_feature(): Got EOF reading end POS: %s.\n",
		end_pos_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	gff_feature->end_pos = strtoul(end_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_gff_read_feature(): Invalid feature position: %s\n",
		    end_pos_str);
	    return BL_READ_TRUNCATED;
	}
    }

    // 6 Score
    if ( (delim = tsv_read_field(gff_stream, score_str,
			BL_GFF_SCORE_MAX_DIGITS, &len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read_feature(): Got EOF reading SCORE: %s.\n",
		score_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	gff_feature->score = strtod(score_str, &end);
	if ( *end != '\0' )
	    gff_feature->score = BL_GFF_SCORE_UNAVAILABLE;
    }
    
    
    // 7 Strand
    if ( (delim = tsv_read_field(gff_stream, strand_str,
			BL_GFF_STRAND_MAX_CHARS, &len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read_feature(): Got EOF reading STRAND: %s.\n",
		strand_str);
	return BL_READ_TRUNCATED;
    }
    else
	gff_feature->strand = *strand_str;
    
    // 8 Phase (bases to start of next codon: 0, 1, or 2. "." if unavailable)
    if ( (delim = tsv_read_field(gff_stream, phase_str,
			BL_GFF_PHASE_MAX_DIGITS, &len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read_feature(): Got EOF reading PHASE: %s.\n",
		phase_str);
	return BL_READ_TRUNCATED;
    }
    else
	gff_feature->phase = *phase_str;

    // 9 Attributes
    if ( (delim = tsv_read_field(gff_stream, temp_attributes,
			BL_GFF_ATTRIBUTES_MAX_CHARS, &len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read_feature(): Got EOF reading ATTRIBUTES: %s.\n",
		temp_attributes);
	return BL_READ_TRUNCATED;
    }
    else if ( (gff_feature->attributes = strdup(temp_attributes)) == NULL )
    {
	fprintf(stderr, "bl_gff_read_feature(): Could not strdup attributes: %s.\n",
		temp_attributes);
	return BL_READ_TRUNCATED;
    }
    
    // printf("delim = %u\n", delim);
    if ( delim != '\n' )
	dsv_skip_rest_of_line(gff_stream);

    if ( (id_start = strchr(gff_feature->attributes, ':')) != NULL )
    {
	++id_start; // Start is first char after ':'
	if ( (id_end = strchr(id_start, ';')) != NULL )
	{
	    *id_end = '\0';
	    gff_feature->feature_id = strdup(id_start);
	    // FIXME: Report malloc() failure somehow
	    *id_end = ';';
	}
    }
    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc
 *
 *  Description:
 *      Write fields from a GFF feature to the specified FILE
 *      stream.
 *
 *      If field_mask is not BL_GFF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate marker for that field,
 *      such as a '.', rather than writing the real data.
 *      Possible mask values are:
 *
 *      BL_GFF_FIELD_ALL
 *      BL_GFF_FIELD_SEQUENCE
 *      BL_GFF_FIELD_SOURCE
 *      BL_GFF_FIELD_NAME
 *      BL_GFF_FIELD_START_POS
 *      BL_GFF_FIELD_END_POS
 *      BL_GFF_FIELD_SCORE
 *      BL_GFF_FIELD_STRAND
 *      BL_GFF_FIELD_PHASE
 *      BL_GFF_FIELD_ATTRIBUTES
 *
 *  Arguments:
 *      gff_stream:     FILE stream to which TSV gff line is written
 *      gff_feature:    Pointer to the bl_gff_t structure to output
 *      field_mask:     Bit mask indicating which fields to output
 *
 *  Returns:
 *      BL_WRITE_OK on success
 *      BL_WRITE_ERROR on failure (errno may provide more information)
 *
 *  Examples:
 *      bl_gff_write_feature(stdout, &gff_feature, BL_GFF_FIELD_ALL);
 *      bl_gff_write_feature(gff_stream, &gff_feature,
 *                        BL_GFF_FIELD_SEQUENCE|BL_GFF_FIELD_START_POS|BL_GFF_FIELD_END_POS);
 *
 *  See also:
 *      bl_gff_read_feature(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_write_feature(FILE *gff_stream, bl_gff_t *gff_feature,
				gff_field_mask_t field_mask)

{
    return fprintf(gff_stream,
	"%s\t%s\t%s\t%" PRIu64 "\t%" PRIu64 "\t%f\t%c\t%c\t%s\n",
	gff_feature->sequence, gff_feature->source, gff_feature->name,
	gff_feature->start_pos, gff_feature->end_pos, gff_feature->score,
	gff_feature->strand, gff_feature->phase, gff_feature->attributes);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc
 *
 *  Description:
 *      Copy GFF fields to a BED structure to the extent possible.  Since
 *      GFF and BED files do not necessarily contain the same information,
 *      some information may be lost or filled in with appropriate markers.
 *
 *  Arguments:
 *      bed_feature: Pointer to the bl_bed_t structure to receive data
 *      gff_feature: Pointer to the bl_gff_t structure to copy
 *
 *  See also:
 *      bl_bed_read_feature(3), bl_gff_read_feature(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-19  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_to_bed(bl_bed_t *bed_feature, bl_gff_t *gff_feature)

{
    char    name[BL_BED_NAME_MAX_CHARS + 1],
	    strand = BL_GFF_STRAND(gff_feature);
    
    BL_BED_SET_CHROMOSOME(bed_feature, BL_GFF_SEQUENCE(gff_feature));
    /*
     *  BED start is 0-based and inclusive
     *  GFF is 1-based and inclusive
     */
    BL_BED_SET_START_POS(bed_feature, BL_GFF_START_POS(gff_feature) - 1);
    /*
     *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
     *  GFF is the same
     */
    BL_BED_SET_END_POS(bed_feature, BL_GFF_END_POS(gff_feature));
    snprintf(name, BL_BED_NAME_MAX_CHARS, "%s", BL_GFF_NAME(gff_feature));
    BL_BED_SET_NAME(bed_feature, name);
    BL_BED_SET_SCORE(bed_feature, 0);  // FIXME: Take as arg?
    if ( BL_BED_SET_STRAND(bed_feature, strand) != BL_DATA_OK )
    {
	fputs("bl_gff_to_bed().\n", stderr);
	exit(EX_DATAERR);
    }
}
