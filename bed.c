#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <inttypes.h>   // PRId64
#include <xtend/math.h> // XT_MAX(), XT_MIN()
#include <xtend/dsv.h>
#include <xtend/mem.h>
#include "bed.h"
#include "gff3.h"
#include "biostring.h"

/***************************************************************************
 *  Name:
 *      bl_bed_skip_header() - Read past BED header
 *
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over header lines in bed input stream, leaving the FILE
 *      structure pointing to the first character in the first line of data.
 *      The header is copied to a temporary file whose FILE pointer 
 *      is returned.
 *
 *  Arguments:
 *      stream  Pointer to the FILE structure for reading the BED stream
 *
 *  Returns:
 *      Pointer to the FILE structure of the temporary file.
 *
 *  Examples:
 *      FILE    *header, *bed_stream;
 *      ...
 *      header = bl_bed_skip_header(bed_stream);
 *
 *  See also:
 *      bl_bed_read(3), xt_fopen(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

FILE    *bl_bed_skip_header(FILE *bed_stream)

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
    
    while ( ((count=fread(start, 1, 7, bed_stream)) == 7) && 
	    ((memcmp(start, "browser", 7) == 0) ||
	    (memcmp(start, "track", 5) == 0) ||
	    (*start == '#')) )
    {
	fwrite(start, count, 1, header_stream);
	do
	{
	    ch = getc(bed_stream);
	    putc(ch, header_stream);
	}   while ( (ch != '\n') && (ch != EOF) );
    }
    
    // Rewind to start of first non-header line
    for (c = count - 1; c >= 0; --c)
	ungetc(start[c], bed_stream);
    rewind(header_stream);
    return header_stream;
}


/***************************************************************************
 *  Name:
 *      bl_bed_read() - Read a BED record
 *
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read next entry (line) from a BED file.  The line must have at
 *      least the first 3 fields (chrom, start, and end).  It may
 *      have up to 12 fields, all of which must be in the correct order
 *      according to the BED specification.
 *
 *      If field_mask is not BL_BED_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in bed_feature.
 *      Possible mask values are:
 *
 *      BL_BED_FIELD_ALL
 *      BL_BED_FIELD_NAME
 *      BL_BED_FIELD_SCORE
 *      BL_BED_FIELD_STRAND
 *      BL_BED_FIELD_THICK
 *      BL_BED_FIELD_RGB
 *      BL_BED_FIELD_BLOCK
 *
 *      The chrom, start, and end fields are required and therefore have
 *      no corresponding mask bits. The thickStart and thickEnd fields must
 *      occur together or not at all, so only a single bit BL_BED_FIELD_THICK
 *      selects both of them.  Likewise, blockCount, blockSizes and
 *      blockStarts must all be present or omitted, so BL_BED_FIELD_BLOCK
 *      masks all three.
 *
 *  Arguments:
 *      bed_feature     Pointer to a bl_bed_t structure
 *      bed_stream      A FILE stream from which to read the line
 *      field_mask      Bit mask indicating which fields to store in bed_feature
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered at the start of a line
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      bl_bed_read(stdin, &bed_feature, BL_BED_FIELD_ALL);
 *      bl_bed_read(bed_stream, &bed_feature,
 *                       BL_BED_FIELD_NAME|BL_BED_FIELD_SCORE);
 *
 *  See also:
 *      bl_bed_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_bed_read(bl_bed_t *bed_feature, FILE *bed_stream,
	    bed_field_mask_t field_mask)

{
    char    *end,
	    strand[BL_BED_STRAND_MAX_CHARS + 1],
	    block_count_str[BL_BED_BLOCK_COUNT_MAX_DIGITS + 1],
	    block_size_str[BL_BED_BLOCK_SIZE_MAX_DIGITS + 1],
	    block_start_str[BL_BED_BLOCK_START_MAX_DIGITS + 1],
	    chrom_start_str[BL_POSITION_MAX_DIGITS + 1],
	    chrom_end_str[BL_POSITION_MAX_DIGITS + 1],
	    score_str[BL_BED_SCORE_MAX_DIGITS + 1],
	    thick_start_str[BL_POSITION_MAX_DIGITS + 1],
	    thick_end_str[BL_POSITION_MAX_DIGITS + 1];
    size_t  len;
    int     delim;
    unsigned long   block_count;
    unsigned    c;
    
    // FIXME: Respect field_mask
    
    // Chromosome
    if ( xt_tsv_read_field(bed_stream, bed_feature->chrom,
			BL_CHROM_MAX_CHARS, &len) == EOF )
    {
	// fputs("bl_bed_read(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return BL_READ_EOF;
    }
    
    // Feature start position
    if ( xt_tsv_read_field(bed_stream, chrom_start_str,
			BL_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_bed_read(): Got EOF reading start position: %s.\n",
		chrom_start_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	bed_feature->chrom_start = strtoul(chrom_start_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_bed_read(): Invalid start position: %s\n",
		    chrom_start_str);
	    return BL_READ_TRUNCATED;
	}
    }
    
    // Feature end position
    // FIXME: Check for > or < start if strand + or -
    if ( (delim = xt_tsv_read_field(bed_stream, chrom_end_str,
			BL_POSITION_MAX_DIGITS, &len)) == EOF )
    {
	fprintf(stderr, "bl_bed_read(): Got EOF reading end position: %s.\n",
		chrom_end_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	bed_feature->chrom_end = strtoul(chrom_end_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_bed_read(): Invalid end position: %s\n",
		    chrom_end_str);
	    return BL_READ_TRUNCATED;
	}
    }

    bed_feature->fields = 3;
    
    // Read NAME field if present
    if ( delim != '\n' )
    {
	if ( (delim = xt_tsv_read_field(bed_stream, bed_feature->name,
			    BL_BED_NAME_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading name: %s.\n",
		    bed_feature->name);
	    return BL_READ_TRUNCATED;
	}
	++bed_feature->fields;
    }
    
    // Read SCORE if present
    if ( delim != '\n' )
    {
	if ( (delim = xt_tsv_read_field(bed_stream, score_str,
			    BL_POSITION_MAX_DIGITS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading score: %s.\n",
		    score_str);
	    return BL_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->score = strtoul(score_str, &end, 10);
	    if ( (*end != '\0') || (bed_feature->score > 1000) )
	    {
		fprintf(stderr,
			"bl_bed_read(): Invalid feature score: %s\n",
			score_str);
		return BL_READ_TRUNCATED;
	    }
	}
	++bed_feature->fields;
    }
    
    // Read strand if present
    if ( delim != '\n' )
    {
	if ( (delim = xt_tsv_read_field(bed_stream, strand,
			    BL_BED_STRAND_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading strand: %s.\n",
		    bed_feature->name);
	    return BL_READ_TRUNCATED;
	}
	if ( (len != 1) || ((*strand != '+') && (*strand != '-') && (*strand != '.')) )
	{
	    fprintf(stderr, "bl_bed_read(): Strand must be + or - or .: %s\n",
		    strand);
	    return BL_READ_TRUNCATED;
	}
	bed_feature->strand = *strand;
	++bed_feature->fields;
    }
    
    // Read thick start position if present
    // Must be followed by thick end position, > or < for + or - strand
    // Feature start position
    if ( delim != '\n' )
    {
	if ( xt_tsv_read_field(bed_stream, thick_start_str,
			    BL_POSITION_MAX_DIGITS, &len) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading thick start "
		    "POS: %s.\n", thick_start_str);
	    return BL_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->thick_start =
		strtoul(thick_start_str, &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bl_bed_read(): Invalid thick start "
				"position: %s\n",
				thick_start_str);
		return BL_READ_TRUNCATED;
	    }
	}
	
	if ( delim == '\n' )
	{
	    fprintf(stderr, "bl_bed_read(): Found thick start, but no thick end.\n");
	    return BL_READ_TRUNCATED;
	}
    
	if ( xt_tsv_read_field(bed_stream, thick_end_str,
			    BL_POSITION_MAX_DIGITS, &len) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading thick end "
		    "POS: %s.\n", thick_end_str);
	    return BL_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->thick_end =
		strtoul(thick_end_str, &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bl_bed_read(): Invalid thick end "
				"position: %s\n",
				thick_end_str);
		return BL_READ_TRUNCATED;
	    }
	}
	bed_feature->fields += 2;
    }

    // Read RGB string field if present
    if ( delim != '\n' )
    {
	if ( (delim = xt_tsv_read_field(bed_stream, bed_feature->item_rgb,
			    BL_BED_ITEM_RGB_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading RGB: %s.\n",
		    bed_feature->name);
	    return BL_READ_TRUNCATED;
	}
	++bed_feature->fields;
    }

    /*
     *  Read block count if present
     *  Must be followed by comma-separated list of sizes
     *  and comma-separated list of start positions
     */
    if ( delim != '\n' )
    {
	if ( (delim = xt_tsv_read_field(bed_stream, block_count_str,
			    BL_BED_BLOCK_COUNT_MAX_DIGITS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading block count: %s.\n",
		    score_str);
	    return BL_READ_TRUNCATED;
	}
	else
	{
	    block_count = strtoul(block_count_str, &end, 10);
	    if ( (*end != '\0') || (block_count > 65535) )
	    {
		fprintf(stderr,
			"bl_bed_read(): Invalid block count: %s\n",
			score_str);
		return BL_READ_TRUNCATED;
	    }
	    bed_feature->block_count = block_count;
	}
	bed_feature->block_sizes = xt_malloc(bed_feature->block_count,
					sizeof(*bed_feature->block_sizes));
	if ( bed_feature->block_sizes == NULL )
	{
	    fputs("bl_bed_read(): Cannot allocate block_sizes.\n", stderr);
	    exit(EX_UNAVAILABLE);
	}
	bed_feature->block_starts = xt_malloc(bed_feature->block_count,
					sizeof(*bed_feature->block_starts));
	if ( bed_feature->block_starts == NULL )
	{
	    fputs("bl_bed_read(): Cannot allocate block_starts.\n", stderr);
	    exit(EX_UNAVAILABLE);
	}
	if ( delim == '\n' )
	{
	    fputs("bl_bed_read(): Found block count, but no sizes.\n", stderr);
	    return BL_READ_TRUNCATED;
	}
	
	// Read comma-separated sizes
	c = 0;
	do
	{
	    delim = xt_dsv_read_field(bed_stream, block_size_str,
			    BL_BED_BLOCK_SIZE_MAX_DIGITS, ",\t", &len);
	    bed_feature->block_sizes[c++] = strtoul(block_size_str, &end, 10);
	    //fprintf(stderr, "Block size[%u] = %s\n", c-1, block_size_str);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bl_bed_read(): Invalid block size: %s\n",
			block_size_str);
		return BL_READ_TRUNCATED;
	    }
	}   while ( delim == ',' );
	if ( c != bed_feature->block_count )
	{
	    fprintf(stderr, "bl_bed_read(): Block count = %u  Sizes = %u\n",
		    bed_feature->block_count, c);
	    return BL_READ_MISMATCH;
	}
	if ( delim == '\n' )
	{
	    fputs("bl_bed_read(): Found block sizes, but no starts.\n", stderr);
	    return BL_READ_TRUNCATED;
	}
	
	// Read comma-separated starts
	c = 0;
	do
	{
	    delim = xt_dsv_read_field(bed_stream, block_start_str,
			    BL_BED_BLOCK_START_MAX_DIGITS, ",\t", &len);
	    bed_feature->block_starts[c++] = strtoul(block_start_str, &end, 10);
	    //fprintf(stderr, "Block start[%u] = %s\n", c-1, block_start_str);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bl_bed_read(): Invalid block start: %s\n",
			block_start_str);
		return BL_READ_TRUNCATED;
	    }
	}   while ( delim == ',' );
	if ( c != bed_feature->block_count )
	{
	    fprintf(stderr, "bl_bed_read(): Block count = %u  Sizes = %u\n",
		    bed_feature->block_count, c);
	    return BL_READ_MISMATCH;
	}
	bed_feature->fields += 3;
    }

    //fprintf(stderr, "Bed fields = %u\n", bed_feature->fields);
    /*
     *  There shouldn't be anything left at this point.  Once block reads
     *  are implemented, we should error out of delim != '\n'
     */
    
    if ( delim != '\n' )
    {
	fputs("bl_bed_read(): Extra columns found.\n", stderr);
	return BL_READ_EXTRA_COLS;
    }
    return BL_READ_OK;
}


/***************************************************************************
 *  Name:
 *      bl_bed_write() - Write a BED record
 *
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write fields from one line of a bed file to the specified FILE
 *      stream.  If field_mask is not BL_BED_FIELD_ALL, only selected fields
 *      are written.
 *
 *      If field_mask is not BL_BED_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate marker for that field,
 *      such as a '.', rather than writing the real data.
 *      Possible mask values are:
 *
 *      BL_BED_FIELD_NAME
 *      BL_BED_FIELD_SCORE
 *      BL_BED_FIELD_STRAND
 *      BL_BED_FIELD_THICK
 *      BL_BED_FIELD_RGB
 *      BL_BED_FIELD_BLOCK
 *
 *      The chrom, start, and end fields are required and therefore have
 *      no corresponding mask bits. The thickStart and thickEnd fields must
 *      occur together or not at all, so only a single bit BL_BED_FIELD_THICK
 *      selects both of them.  Likewise, blockCount, blockSizes and
 *      blockStarts must all be present or omitted, so BL_BED_FIELD_BLOCK
 *      masks all three.
 *
 *  Arguments:
 *      bed_feature     Pointer to the bl_bed_t structure to output
 *      bed_stream      FILE stream to which TSV bed line is written
 *      field_mask      Bit mask indicating which fields to output
 *
 *  Returns:
 *      BL_WRITE_OK on success
 *      BL_WRITE_ERROR on failure (errno may provide more information)
 *
 *  Examples:
 *      bl_bed_write(stdout, &bed_feature, BL_BED_FIELD_ALL);
 *      bl_bed_write(bed_stream, &bed_feature,
 *                        BL_BED_FIELD_NAME|BL_BED_FIELD_SCORE);
 *
 *  See also:
 *      bl_bed_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_bed_write(bl_bed_t *bed_feature, FILE *bed_stream,
	    bed_field_mask_t field_mask)

{
    unsigned    c;
    
    // FIXME: Respect field_mask
    // FIXME: Check fprintf() return codes
    fprintf(bed_stream, "%s\t%" PRId64 "\t%" PRId64,
	    bed_feature->chrom,
	    bed_feature->chrom_start, bed_feature->chrom_end);
    if ( bed_feature->fields > 3 )
	fprintf(bed_stream, "\t%s", bed_feature->name);
    if ( bed_feature->fields > 4 )
	fprintf(bed_stream, "\t%u", bed_feature->score);
    if ( bed_feature->fields > 5 )
	fprintf(bed_stream, "\t%c", bed_feature->strand);
    if ( bed_feature->fields > 6 )
	fprintf(bed_stream, "\t%" PRId64 "\t%" PRId64,
		bed_feature->thick_start, bed_feature->thick_end);
    if ( bed_feature->fields > 8 )
	fprintf(bed_stream, "\t%s", bed_feature->item_rgb);
    if ( bed_feature->fields > 9 )
    {
	fprintf(bed_stream, "\t%u\t", bed_feature->block_count);
	for (c = 0; c < bed_feature->block_count - 1; ++c)
	    fprintf(bed_stream, "%" PRId64 ",", bed_feature->block_sizes[c]);
	fprintf(bed_stream, "%" PRId64 "\t", bed_feature->block_sizes[c]);
	for (c = 0; c < bed_feature->block_count - 1; ++c)
	    fprintf(bed_stream, "%" PRId64 ",", bed_feature->block_starts[c]);
	fprintf(bed_stream, "%" PRId64, bed_feature->block_starts[c]);
    }
    putc('\n', bed_stream);
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Name:
 *      bl_bed_check_order() - Compare positions of two bed records
 *
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Make sure the BED input is sorted by chrom and start position.
 *
 *  Arguments:
 *      bed_feature     Pointer to BED structure containing current entry
 *      last_chrom      Chromosome of the previous BED entry
 *      last_start      Start position of the previous BED entry
 *
 *  Returns:
 *      Nothing: Terminates process if input is out of order
 *
 *  See also:
 *      bl_bed_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-08  Jason Bacon Begin
 ***************************************************************************/

void    bl_bed_check_order(bl_bed_t *bed_feature, char last_chrom[],
			   int64_t last_start)

{
    if ( bl_chrom_name_cmp(bed_feature->chrom, last_chrom) == 0 )
    {
	if ( bed_feature->chrom_start < last_start )
	{
	    fprintf(stderr, "peak-classifier: BED file not sorted by start position.\n");
	    exit(EX_DATAERR);
	}
    }
    else if ( bl_chrom_name_cmp(bed_feature->chrom, last_chrom) < 0 )
    {
	fprintf(stderr, "peak-classifier: BED file not sorted by chrom.\n");
	fprintf(stderr, "%s, %s\n", bed_feature->chrom, last_chrom);
	exit(EX_DATAERR);
    }
}

/***************************************************************************
 *  Name:
 *      bl_bed_gff3_cmp() - Compare positions of BED and GFF3 objects
 *
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Compare the position of a BED feature to that of a GFF feature.
 *      Return 0 if the features overlap, < 0 if the BED feature is upstream
 *      of the GFF feature, > 0 if the BED feature is downstream of the GFF
 *      feature.
 *
 *      If the features overlap, populate the bl_overlap_t structure
 *      pointed to by overlap.  The structure contains the lengths of the
 *      two features, the start and end positions of the overlapping region,
 *      and the length of the overlap.  Positions in overlap are 1-based and
 *      inclusive at both ends (like most bioinformatics formats and unlike
 *      BED).
 *
 *  Arguments:
 *      bed_feature     Pointer to the bl_bed_t structure to compare
 *      gff3_feature     Pointer to the bl_gff3_t structure to compare
 *      overlap         Pointer to the bl_overlap_t structure to receive
 *                      comparison results
 *
 *  Returns:
 *      A value < 0 if the BED feature is upstream of the GFF feature
 *      A value > 0 if the BED feature is downstream of the GFF feature
 *      0 if the BED feature overlaps the GFF feature
 *
 *  Examples:
 *      if ( bl_bed_gff3_cmp(&bed_feature, &gff3_feature, *overlap) == 0 )
 *
 *  See also:
 *      bl_bed_read(3), bl_gff3_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_bed_gff3_cmp(bl_bed_t *bed_feature, bl_gff3_t *gff3_feature,
		       bl_overlap_t *overlap)

{
    int         chrom_cmp;
    int64_t    bed_start, bed_end, bed_len,
		gff3_start, gff3_end, gff3_len;
    
    chrom_cmp = bl_chrom_name_cmp(BL_BED_CHROM(bed_feature),
					 BL_GFF3_SEQID(gff3_feature));
    if ( chrom_cmp == 0 )
    {
	/*
	 *  BED positions are 0-based, with end non-inclusive, which can
	 *  also be viewed as an inclusive 1-based coordinate
	 *  GFF is 1-based, both ends inclusive
	 */
	
	if ( BL_BED_CHROM_END(bed_feature) < BL_GFF3_START(gff3_feature) )
	{
	    bl_overlap_set_all(overlap, 0, 0, 0, 0);
	    return -1;
	}
	else if ( BL_BED_CHROM_START(bed_feature) + 1 > BL_GFF3_END(gff3_feature) )
	{
	    bl_overlap_set_all(overlap, 0, 0, 0, 0);
	    return 1;
	}
	else
	{
	    bed_start = BL_BED_CHROM_START(bed_feature);
	    bed_end = BL_BED_CHROM_END(bed_feature);
	    gff3_start = BL_GFF3_START(gff3_feature);
	    gff3_end = BL_GFF3_END(gff3_feature);
	    bed_len = bed_end - bed_start;
	    gff3_len = gff3_end - gff3_start + 1;
	    bl_overlap_set_all(overlap, bed_len, gff3_len,
			    XT_MAX(bed_start+1, gff3_start),
			    XT_MIN(bed_end, gff3_end));
	    return 0;
	}
    }
    return chrom_cmp;
}
