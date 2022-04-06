#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <inttypes.h>       // PRId64
#include <xtend/string.h>   // strlcpy() on Linux
#include <xtend/dsv.h>
#include <xtend/mem.h>
#include "gff.h"
#include "bed.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over header lines in gff input stream.  The FILE pointer
 *      gff_stream is advanced to the first character of the first line
 *      after the header.  The header is copied to a temporary file and and
 *      the function returns a FILE pointer to the header stream.
 *
 *  Arguments:
 *      gff_stream  FILE pointer to the open GFF file
 *
 *  Returns:
 *      A FILE pointer to a temporary file containing a copy of the header
 *
 *  See also:
 *      bl_gff_read(3), bl_gff_copy_header(3)
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
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Copy GFF header from one FILE stream to another.  This is meant to
 *      be used in conjunction with bl_gff_skip_header(), which stores the
 *      header in a temporary file.
 *
 *  Arguments:
 *      header_stream   Open FILE stream of GFF header
 *      gff_stream      FILE stream to which header is copied
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE or BL_READ_* on failure
 *
 *  See also:
 *      bl_gff_skip_header(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-03  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_copy_header(FILE *header_stream, FILE *gff_stream)

{
    int     ch;
    
    rewind(header_stream);
    while ( (ch = getc(header_stream)) != EOF )
	if ( putc(ch, gff_stream) == EOF )
	    return BL_WRITE_FAILURE;
    rewind(header_stream);
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read next feature (line) from a GFF file.
 *
 *      gff_feature must be initialized using BL_GFF_INIT or bl_gff_init()
 *      before being passed to this function.
 *
 *      bl_gff_read() will allocate memory for string fields as needed.
 *      The object should be passed to bl_gff_free() as soon as possible
 *      after the data are no longer needed.
 *
 *      If passed an object that is not in an initialized state,
 *      bl_gff_read() will free and initialize it before repopulating it
 *      with a new feature.
 *
 *      If field_mask is not BL_GFF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in gff_feature.
 *      That field in the structure is then populated with an appropriate
 *      marker, such as '.'.  Possible mask values are:
 *
 *      BL_GFF_FIELD_ALL
 *      BL_GFF_FIELD_SEQID
 *      BL_GFF_FIELD_SOURCE
 *      BL_GFF_FIELD_TYPE
 *      BL_GFF_FIELD_START
 *      BL_GFF_FIELD_END
 *      BL_GFF_FIELD_SCORE
 *      BL_GFF_FIELD_STRAND
 *      BL_GFF_FIELD_PHASE
 *      BL_GFF_FIELD_ATTRIBUTES
 *
 *  Arguments:
 *      gff_feature     Pointer to a bl_gff_t structure
 *      gff_stream      A FILE stream from which to read the line
 *      field_mask      Bit mask indicating which fields to store in gff_feature
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered after a complete feature
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      bl_gff_skip_header(stdin);
 *      bl_gff_init(&gff_feature);
 *
 *      bl_gff_read(&gff_feature, stdin, BL_GFF_FIELD_ALL);
 *      bl_gff_read(&gff_feature, gff_stream,
 *          BL_GFF_FIELD_SEQID|BL_GFF_FIELD_START|BL_GFF_FIELD_END);
 *
 *  See also:
 *      bl_gff_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_read(bl_gff_t *gff_feature, FILE *gff_stream,
	    gff_field_mask_t field_mask)

{
    char    *end,
	    line[BL_GFF_LINE_MAX_CHARS + 1],
	    temp_attributes[BL_GFF_ATTRIBUTES_MAX_CHARS + 1],
	    strand_str[BL_GFF_STRAND_MAX_CHARS + 1],
	    phase_str[BL_GFF_PHASE_MAX_DIGITS + 1],
	    start_str[BL_POSITION_MAX_DIGITS + 1],
	    end_str[BL_POSITION_MAX_DIGITS + 1],
	    score_str[BL_GFF_SCORE_MAX_DIGITS + 1];
    size_t  len;
    int     delim,
	    ch;
    
    // Use this as a model for other _read() functions?
    // Makes reusing a structure easy without risk of memory leaks
    if ( gff_feature->attributes != NULL )
	bl_gff_free(gff_feature);
    
    // Check for group terminators (Line with just ###)
    // FIXME: Rely on parent ID instead of ###?
    if ( (ch = getc(gff_stream)) == '#' )
    {
	fgets(line, BL_GFF_LINE_MAX_CHARS, gff_stream);
	if ( strcmp(line, "##\n") == 0 )
	{
	    strlcpy(gff_feature->type, "###", BL_GFF_TYPE_MAX_CHARS);
	    return BL_READ_OK;
	}
    }
    else if ( ch != EOF )
	ungetc(ch, gff_stream);

    gff_feature->file_pos = ftell(gff_stream);
    
    // FIXME: Respect field_mask
    
    // 1 Chromosome
    if ( tsv_read_field(gff_stream, gff_feature->seqid,
			BL_CHROM_MAX_CHARS, &len) == EOF )
    {
	return BL_READ_EOF;
    }
    
    // 2 Source
    if ( tsv_read_field(gff_stream, gff_feature->source,
			BL_GFF_SOURCE_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading SOURCE: %s.\n",
		gff_feature->source);
	return BL_READ_TRUNCATED;
    }

    // 3 Feature
    if ( tsv_read_field(gff_stream, gff_feature->type,
			BL_GFF_TYPE_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading feature: %s.\n",
		gff_feature->type);
	return BL_READ_TRUNCATED;
    }
    
    // 4 Feature start position
    if ( tsv_read_field(gff_stream, start_str,
			BL_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading start POS: %s.\n",
		start_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	gff_feature->start = strtoul(start_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_gff_read(): Invalid feature position: %s\n",
		    start_str);
	    return BL_READ_TRUNCATED;
	}
    }
    
    // 5 Feature end position
    if ( tsv_read_field(gff_stream, end_str,
			BL_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading end POS: %s.\n",
		end_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	gff_feature->end = strtoul(end_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_gff_read(): Invalid feature position: %s\n",
		    end_str);
	    return BL_READ_TRUNCATED;
	}
    }

    // 6 Score
    if ( tsv_read_field(gff_stream, score_str,
			BL_GFF_SCORE_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading SCORE: %s.\n",
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
    if ( tsv_read_field(gff_stream, strand_str,
			BL_GFF_STRAND_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading STRAND: %s.\n",
		strand_str);
	return BL_READ_TRUNCATED;
    }
    else
	gff_feature->strand = *strand_str;
    
    // 8 Phase (bases to start of next codon: 0, 1, or 2. "." if unavailable)
    if ( tsv_read_field(gff_stream, phase_str,
			BL_GFF_PHASE_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading PHASE: %s.\n",
		phase_str);
	return BL_READ_TRUNCATED;
    }
    else
	gff_feature->phase = *phase_str;

    // 9 Attributes
    if ( (delim = tsv_read_field(gff_stream, temp_attributes,
			BL_GFF_ATTRIBUTES_MAX_CHARS, &len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading ATTRIBUTES: %s.\n",
		temp_attributes);
	return BL_READ_TRUNCATED;
    }
    if ( (gff_feature->attributes = strdup(temp_attributes)) == NULL )
    {
	fprintf(stderr, "bl_gff_read(): Could not strdup attributes: %s.\n",
		temp_attributes);
	return BL_READ_TRUNCATED;
    }
    //fprintf(stderr, "%s %zu\n", gff_feature->attributes,
    //        strlen(gff_feature->attributes));
    
    // printf("delim = %u\n", delim);
    if ( delim != '\n' )
	dsv_skip_rest_of_line(gff_stream);

    // Extract feature ID from attributes
    gff_feature->feature_id = bl_gff_extract_attribute(gff_feature, "ID");

    // Extract feature ID from attributes
    gff_feature->feature_name = bl_gff_extract_attribute(gff_feature, "Name");
    if ( gff_feature->feature_name == NULL )
    {
	if ( (gff_feature->feature_name = strdup("unnamed")) == NULL )
	    fprintf(stderr, "bl_gff_read(): Could not strdup() feature_name.\n");
    }
    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
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
 *      BL_GFF_FIELD_SEQID
 *      BL_GFF_FIELD_SOURCE
 *      BL_GFF_FIELD_TYPE
 *      BL_GFF_FIELD_START
 *      BL_GFF_FIELD_END
 *      BL_GFF_FIELD_SCORE
 *      BL_GFF_FIELD_STRAND
 *      BL_GFF_FIELD_PHASE
 *      BL_GFF_FIELD_ATTRIBUTES
 *
 *  Arguments:
 *      gff_feature     Pointer to the bl_gff_t structure to output
 *      gff_stream      FILE stream to which TSV gff line is written
 *      field_mask      Bit mask indicating which fields to output
 *
 *  Returns:
 *      BL_WRITE_OK on success
 *      BL_WRITE_ERROR on failure (errno may provide more information)
 *
 *  Examples:
 *      bl_gff_write(&gff_feature, stdout, BL_GFF_FIELD_ALL);
 *      bl_gff_write(&gff_feature, gff_stream,
 *          BL_GFF_FIELD_SEQID|BL_GFF_FIELD_START|BL_GFF_FIELD_END);
 *
 *  See also:
 *      bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_write(bl_gff_t *gff_feature, FILE *gff_stream,
	    gff_field_mask_t field_mask)

{
    int     printed = 0;
    
    /* FIXME: Fully test and enable this
    if ( field_mask & BL_GFF_FIELD_SEQID )
	printed += fprintf(gff_stream, "%s", gff_feature->seqid);
    else
	printed += putc('.', gff_stream);
	
    if ( field_mask & BL_GFF_FIELD_SOURCE )
	printed += fprintf(gff_stream, "\t%s", gff_feature->source);
    else
	printed += fprintf(gff_stream, "\t.");

    if ( field_mask & BL_GFF_FIELD_TYPE )
	printed += fprintf(gff_stream, "\t%s", gff_feature->type);
    else
	printed += fprintf(gff_stream, "\t.");

    if ( field_mask & BL_GFF_FIELD_START )
	printed += fprintf(gff_stream, "\t%" PRId64, gff_feature->start);
    else
	printed += fprintf(gff_stream, "\t-1");
    
    if ( field_mask & BL_GFF_FIELD_END )
	printed += fprintf(gff_stream, "\t%" PRId64, gff_feature->end);
    else
	printed += fprintf(gff_stream, "\t-1");
    
    if ( field_mask & BL_GFF_FIELD_SCORE )
	printed += fprintf(gff_stream, "\t%f", gff_feature->score);
    else
	printed += fprintf(gff_stream, "\t.");
    
    if ( field_mask & BL_GFF_FIELD_STRAND )
	printed += fprintf(gff_stream, "\t%c", gff_feature->strand);
    else
	printed += fprintf(gff_stream, "\t.");
    
    if ( field_mask & BL_GFF_FIELD_PHASE )
	printed += fprintf(gff_stream, "\t%c", gff_feature->phase);
    else
	printed += fprintf(gff_stream, "\t.");
    
    if ( field_mask & BL_GFF_FIELD_ATTRIBUTES )
	printed += fprintf(gff_stream, "\t%s", gff_feature->attributes);
    else
	printed += fprintf(gff_stream, "\t.");
    putc('\n', gff_stream);
    */
    fprintf(gff_stream,
	"%s\t%s\t%s\t%" PRId64 "\t%" PRId64 "\t%f\t%c\t%c\t%s\n",
	gff_feature->seqid, gff_feature->source, gff_feature->type,
	gff_feature->start, gff_feature->end, gff_feature->score,
	gff_feature->strand, gff_feature->phase, gff_feature->attributes);
    return printed;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Copy GFF fields to a BED structure to the extent possible.  Since
 *      GFF and BED files do not necessarily contain the same information,
 *      some information may be lost or filled in with appropriate markers.
 *
 *  Arguments:
 *      gff_feature  Pointer to the bl_gff_t structure to copy
 *      bed_feature  Pointer to the bl_bed_t structure to receive data
 *
 *  See also:
 *      bl_bed_read(3), bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-19  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_to_bed(bl_gff_t *gff_feature, bl_bed_t *bed_feature)

{
    char    name[BL_BED_NAME_MAX_CHARS + 1],
	    strand = BL_GFF_STRAND(gff_feature);
    
    // Update this if/when more fields are converted
    bl_bed_set_fields(bed_feature, 6);
    bl_bed_set_score(bed_feature, 0);
    
    bl_bed_set_chrom_cpy(bed_feature, BL_GFF_SEQID(gff_feature), BL_CHROM_MAX_CHARS + 1);
    /*
     *  BED start is 0-based and inclusive
     *  GFF is 1-based and inclusive
     */
    bl_bed_set_chrom_start(bed_feature, BL_GFF_START(gff_feature) - 1);
    /*
     *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
     *  GFF is the same
     */
    bl_bed_set_chrom_end(bed_feature, BL_GFF_END(gff_feature));
    snprintf(name, BL_BED_NAME_MAX_CHARS + 1, "%s", BL_GFF_TYPE(gff_feature));
    bl_bed_set_name_cpy(bed_feature, name, BL_BED_NAME_MAX_CHARS + 1);
    bl_bed_set_score(bed_feature, 0);  // FIXME: Take as arg?
    if ( bl_bed_set_strand(bed_feature, strand) != BL_BED_DATA_OK )
    {
	fputs("bl_gff_to_bed(): bl_bed_set_strand() failed.\n", stderr);
	exit(EX_DATAERR);
    }
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated for a bl_gff_t object
 *  
 *  Arguments:
 *      gff_feature     Pointer to the bl_gff_t object
 *
 *  Examples:
 *      bl_gff_t    feature;
 *
 *      bl_gff_free(&feature);
 *
 *  See also:
 *      bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_free(bl_gff_t *gff_feature)

{
    if ( gff_feature->attributes != NULL )
    {
	/*fprintf(stderr, "Freeing %s %p %zu %s\n",
		gff_feature->type, gff_feature->attributes,
		strlen(gff_feature->attributes), gff_feature->attributes);
	fflush(stderr);
	*/
	free(gff_feature->attributes);
    }
    if ( gff_feature->feature_id != NULL )
	free(gff_feature->feature_id);
    if ( gff_feature->feature_name != NULL )
	free(gff_feature->feature_name);
    bl_gff_init(gff_feature);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Extract an attribute value of a feature given the attribute name.
 *      Common attribute names include "ID" and "Name".  Attributes are
 *      embedded in the GFF attributes field in the form name=value;, e.g.
 *      ID=gene:ENSDARG00000029944;Name=parpbp.
 *  
 *  Arguments:
 *      Attribute name, such as "ID" or "Name"
 *
 *  Returns:
 *      Attribute value (text after '='), or NULL if name is not found
 *
 *  Examples:
 *      bl_gff_t    gff_feature;
 *
 *      if ( bl_gff_extract_attribute(&gff_feature, "Name") != NULL )
 *      {
 *      }
 *
 *  See also:
 *      bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-05  Jason Bacon Begin
 ***************************************************************************/

char    *bl_gff_extract_attribute(bl_gff_t *gff_feature, const char *attr_name)

{
    char    *attribute = NULL,
	    *start,
	    *val_start,
	    *end;
    size_t  len = strlen(attr_name);

    // Find attribute beginning with "attr_name="
    for (start = gff_feature->attributes; (*start != '\0'); ++start)
    {
	if ( (memcmp(start, attr_name, len) == 0) && (start[len] == '=') )
	{
	    val_start = start + len + 1;
	    // ; separates attributes, last one terminated by null byte
	    if ( (end = strchr(val_start, ';')) != NULL )
		*end = '\0';    // Not thread safe
	    if ( (attribute = strdup(val_start)) == NULL )
		fprintf(stderr, "%s: strdup() failed.\n", __FUNCTION__);
	    if ( end != NULL )
		*end = ';';
	}
    }
    return attribute;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_gff_t object, setting all fields to sentinel
 *      values such as 0, NULL, or '.' as appropriate.  Note that bl_gff_t
 *      objects defined as structures, not pointers to structures, can
 *      also be initialized with the BL_GFF_INIT macro.
 *  
 *  Arguments:
 *      feature     Address of a bl_gff_t structure
 *
 *  Examples:
 *      bl_gff_t    feature1 = BL_GFF_INIT,
 *                  *feature2;
 *
 *      if ( (feature2 = xt_malloc(1, sizeof(*feature2))) != NULL )
 *          bl_gff_init(feature2);
 *
 *  See also:
 *      bl_gff_read(3), bl_gff_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-16  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_init(bl_gff_t *feature)

{
    feature->seqid[0] = feature->source[0] = feature->type[0] = '.';
    feature->seqid[1] = feature->source[1] = feature->type[1] = '\0';
    feature->start = feature->end = 0;
    feature->score = 0.0;
    feature->strand = feature->phase = '.';
    feature->attributes = feature->feature_id = feature->feature_name = NULL;
    feature->file_pos = 0;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-23  Jason Bacon Begin
 ***************************************************************************/

bl_gff_t    *bl_gff_dup(bl_gff_t *feature)

{
    bl_gff_t    *copy;
    
    if ( (copy = xt_malloc(1, sizeof(bl_gff_t))) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate new bl_gff_t object.\n",
		__FUNCTION__);
	return NULL;
    }
    bl_gff_init(copy);
    return bl_gff_copy(copy, feature);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-23  Jason Bacon Begin
 ***************************************************************************/

bl_gff_t    *bl_gff_copy(bl_gff_t *copy, bl_gff_t *feature)

{
    strlcpy(copy->seqid, feature->seqid, BL_CHROM_MAX_CHARS + 1);
    strlcpy(copy->source, feature->source, BL_GFF_SOURCE_MAX_CHARS + 1);
    strlcpy(copy->type, feature->type, BL_GFF_TYPE_MAX_CHARS + 1);
    copy->start = feature->start;
    copy->end = feature->end;
    copy->score = feature->score;
    copy->strand = feature->strand;
    copy->phase = feature->phase = '.';
    
    if ( (copy->attributes = strdup(feature->attributes)) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate attributes.\n", __FUNCTION__);
	free(copy);
	return NULL;
    }
    
    if ( feature->feature_id == NULL )
	copy->feature_id = NULL;
    else if ( (copy->feature_id = strdup(feature->feature_id)) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate attributes.\n", __FUNCTION__);
	free(copy->attributes);
	free(copy);
	return NULL;
    }

    if ( feature->feature_name == NULL )
	copy->feature_name = NULL;
    else if ( (copy->feature_name = strdup(feature->feature_name)) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate attributes.\n", __FUNCTION__);
	free(copy->attributes);
	free(copy->feature_id);
	free(copy);
	return NULL;
    }
    
    copy->file_pos = feature->file_pos;
    
    return copy;
}
