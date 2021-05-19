#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <inttypes.h>   // PRIu64
#include <sys/param.h>  // MAX(), MIN()
#include "bedio.h"
#include "dsvio.h"
#include "gffio.h"
#include "biostring.h"

/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Skip over header lines in bed input stream.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
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
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read fields from one line of a BED file.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bed_read_feature(FILE *bed_stream, bed_feature_t *bed_feature)

{
    char    *end,
	    strand[BED_STRAND_MAX_CHARS + 1],
	    block_count_str[BED_BLOCK_COUNT_MAX_DIGITS + 1],
	    block_size_str[BED_BLOCK_SIZE_MAX_DIGITS + 1],
	    block_start_str[BED_BLOCK_START_MAX_DIGITS + 1];
    size_t  len;
    int     delim;
    unsigned long   block_count;
    unsigned    c;
    
    // Chromosome
    if ( tsv_read_field(bed_stream, bed_feature->chromosome,
			BIO_CHROMOSOME_MAX_CHARS, &len) == EOF )
    {
	// fputs("bed_read_feature(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return BIO_READ_EOF;
    }
    
    // Feature start position
    if ( tsv_read_field(bed_stream, bed_feature->start_pos_str,
			BIO_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bed_read_feature(): Got EOF reading start position: %s.\n",
		bed_feature->start_pos_str);
	return BIO_READ_TRUNCATED;
    }
    else
    {
	bed_feature->start_pos = strtoul(bed_feature->start_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bed_read_feature(): Invalid start position: %s\n",
		    bed_feature->start_pos_str);
	    return BIO_READ_TRUNCATED;
	}
    }
    
    // Feature end position
    // FIXME: Check for > or < start if strand + or -
    if ( (delim = tsv_read_field(bed_stream, bed_feature->end_pos_str,
			BIO_POSITION_MAX_DIGITS, &len)) == EOF )
    {
	fprintf(stderr, "bed_read_feature(): Got EOF reading end position: %s.\n",
		bed_feature->end_pos_str);
	return BIO_READ_TRUNCATED;
    }
    else
    {
	bed_feature->end_pos = strtoul(bed_feature->end_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bed_read_feature(): Invalid end position: %s\n",
		    bed_feature->end_pos_str);
	    return BIO_READ_TRUNCATED;
	}
    }

    bed_feature->fields = 3;
    
    // Read NAME field if present
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, bed_feature->name,
			    BED_NAME_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading name: %s.\n",
		    bed_feature->name);
	    return BIO_READ_TRUNCATED;
	}
	++bed_feature->fields;
    }
    
    // Read SCORE if present
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, bed_feature->score_str,
			    BIO_POSITION_MAX_DIGITS, &len)) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading score: %s.\n",
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
    
    // Read strand if present
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, strand,
			    BED_STRAND_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading strand: %s.\n",
		    bed_feature->name);
	    return BIO_READ_TRUNCATED;
	}
	if ( (len != 1) || ((*strand != '+') && (*strand != '-')) )
	{
	    fprintf(stderr, "bed_read_feature(): Strand must be + or -: %s\n",
		    strand);
	    return BIO_READ_TRUNCATED;
	}
	bed_feature->strand = *strand;
	++bed_feature->fields;
    }
    
    // Read thick start position if present
    // Must be followed by thick end position, > or < for + or - strand
    // Feature start position
    if ( delim != '\n' )
    {
	if ( tsv_read_field(bed_stream, bed_feature->thick_start_pos_str,
			    BIO_POSITION_MAX_DIGITS, &len) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading thick start "
		    "POS: %s.\n", bed_feature->thick_start_pos_str);
	    return BIO_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->thick_start_pos =
		strtoul(bed_feature->thick_start_pos_str, &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bed_read_feature(): Invalid thick start "
				"position: %s\n",
				bed_feature->thick_start_pos_str);
		return BIO_READ_TRUNCATED;
	    }
	}
	
	if ( delim == '\n' )
	{
	    fprintf(stderr, "bed_read_feature(): Found thick start, but no thick end.\n");
	    return BIO_READ_TRUNCATED;
	}
    
	if ( tsv_read_field(bed_stream, bed_feature->thick_end_pos_str,
			    BIO_POSITION_MAX_DIGITS, &len) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading thick end "
		    "POS: %s.\n", bed_feature->thick_end_pos_str);
	    return BIO_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->thick_end_pos =
		strtoul(bed_feature->thick_end_pos_str, &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bed_read_feature(): Invalid thick end "
				"position: %s\n",
				bed_feature->thick_end_pos_str);
		return BIO_READ_TRUNCATED;
	    }
	}
    }

    // Read RGB string field if present
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, bed_feature->name,
			    BED_RGB_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading RGB: %s.\n",
		    bed_feature->name);
	    return BIO_READ_TRUNCATED;
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
	if ( (delim = tsv_read_field(bed_stream, block_count_str,
			    BED_BLOCK_COUNT_MAX_DIGITS, &len)) == EOF )
	{
	    fprintf(stderr, "bed_read_feature(): Got EOF reading block count: %s.\n",
		    bed_feature->score_str);
	    return BIO_READ_TRUNCATED;
	}
	else
	{
	    block_count = strtoul(block_count_str, &end, 10);
	    if ( (*end != '\0') || (block_count > 65535) )
	    {
		fprintf(stderr,
			"bed_read_feature(): Invalid block count: %s\n",
			bed_feature->score_str);
		return BIO_READ_TRUNCATED;
	    }
	    bed_feature->block_count = block_count;
	}
	++bed_feature->fields;
	bed_feature->block_sizes = xt_malloc(bed_feature->block_count,
					sizeof(*bed_feature->block_sizes));
	if ( bed_feature->block_sizes == NULL )
	{
	    fputs("bed_read_feature(): Cannot allocate block_sizes.\n", stderr);
	    exit(EX_UNAVAILABLE);
	}
	bed_feature->block_starts = xt_malloc(bed_feature->block_count,
					sizeof(*bed_feature->block_starts));
	if ( bed_feature->block_starts == NULL )
	{
	    fputs("bed_read_feature(): Cannot allocate block_starts.\n", stderr);
	    exit(EX_UNAVAILABLE);
	}
	if ( delim == '\n' )
	{
	    fputs("bed_read_feature(): Found block count, but no sizes.\n", stderr);
	    return BIO_READ_TRUNCATED;
	}
	
	// Read comma-separated sizes
	c = 0;
	while ( (delim = dsv_read_field(bed_stream, block_size_str,
			    BED_BLOCK_SIZE_MAX_DIGITS, ",\t", &len)) == ',' )
	{
	    bed_feature->block_sizes[c++] = strtoul(block_size_str, &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bed_read_feature(): Invalid block size: %s\n",
			block_size_str);
		return BIO_READ_TRUNCATED;
	    }
	}
	if ( c != bed_feature->block_count )
	{
	    fprintf(stderr, "bed_read_feature(): Block count = %u  Sizes = %u\n",
		    bed_feature->block_count, c);
	    return BIO_READ_MISMATCH;
	}
	if ( delim == '\n' )
	{
	    fputs("bed_read_feature(): Found block sizes, but no starts.\n", stderr);
	    return BIO_READ_TRUNCATED;
	}
	
	// Read comma-separated starts
	c = 0;
	while ( (delim = dsv_read_field(bed_stream, block_start_str,
			    BED_BLOCK_START_MAX_DIGITS, ",\t", &len)) == ',' )
	{
	    bed_feature->block_starts[c++] = strtoul(block_start_str, &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bed_read_feature(): Invalid block start: %s\n",
			block_start_str);
		return BIO_READ_TRUNCATED;
	    }
	}
	if ( c != bed_feature->block_count )
	{
	    fprintf(stderr, "bed_read_feature(): Block count = %u  Sizes = %u\n",
		    bed_feature->block_count, c);
	    return BIO_READ_MISMATCH;
	}
    }

    /*
     *  There shouldn't be anything left at this point.  Once block reads
     *  are implemented, we should error out of delim != '\n'
     */
    
    if ( delim != '\n' )
    {
	fputs("bed_read_feature(): Extra columns found.\n", stderr);
	return BIO_READ_EXTRA_COLS;
    }
    return BIO_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Write fields from one line of a bed file.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
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
    if ( bed_feature->fields > 5 )
	fprintf(bed_stream, "\t%c", bed_feature->strand);
    if ( bed_feature->fields > 6 )
	fprintf(bed_stream, "\t%s\t%s", bed_feature->thick_start_pos_str,
		bed_feature->thick_end_pos_str);
    if ( bed_feature->fields > 8 )
	fprintf(bed_stream, "\t%s", bed_feature->rgb_str);
    putc('\n', bed_stream);
    return 0;
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Make sure the BED input is sorted
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-08  Jason Bacon Begin
 ***************************************************************************/

void    bed_check_order(bed_feature_t *bed_feature, char last_chrom[],
			uint64_t last_pos)

{
    if ( chromosome_name_cmp(BED_CHROMOSOME(bed_feature), last_chrom) == 0 )
    {
	if ( BED_START_POS(bed_feature) < last_pos )
	{
	    fprintf(stderr, "peak-classifier: BED file not sorted by start position.\n");
	    exit(EX_DATAERR);
	}
    }
    else if ( chromosome_name_cmp(BED_CHROMOSOME(bed_feature), last_chrom) < 0 )
    {
	fprintf(stderr, "peak-classifier: BED file not sorted by chromosome.\n");
	fprintf(stderr, "%s, %s\n", BED_CHROMOSOME(bed_feature), last_chrom);
	exit(EX_DATAERR);
    }
}

/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Compare the position of a BED feature to that of a GFF feature.
 *      Return < 0 if the BED feature is "earlier" (lower chromosome or
 *      same chromosome and lower position), etc. much like strcmp().
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bed_gff_cmp(bed_feature_t *bed_feature, gff_feature_t *gff_feature,
		    bio_overlap_t *overlap)

{
    int         chromosome_cmp;
    uint64_t    bed_start, bed_end, bed_len,
		gff_start, gff_end, gff_len;
    
    chromosome_cmp = chromosome_name_cmp(BED_CHROMOSOME(bed_feature),
					 GFF_SEQUENCE(gff_feature));
    if ( chromosome_cmp == 0 )
    {
	/*
	 *  BED positions are 0-based, with end non-inclusive, which can
	 *  also be viewed as an inclusive 1-based coordinate
	 *  GFF is 1-based, both ends inclusive
	 */
	
	if ( BED_END_POS(bed_feature) < GFF_START_POS(gff_feature) )
	{
	    bio_set_overlap(overlap, 0, 0, 0, 0);
	    return -1;
	}
	else if ( BED_START_POS(bed_feature) + 1 > GFF_END_POS(gff_feature) )
	{
	    bio_set_overlap(overlap, 0, 0, 0, 0);
	    return 1;
	}
	else
	{
	    bed_start = BED_START_POS(bed_feature);
	    bed_end = BED_END_POS(bed_feature);
	    gff_start = GFF_START_POS(gff_feature);
	    gff_end = GFF_END_POS(gff_feature);
	    bed_len = bed_end - bed_start;
	    gff_len = gff_end - gff_start + 1;
	    bio_set_overlap(overlap, bed_len, gff_len,
			    MAX(bed_start+1, gff_start),
			    MIN(bed_end, gff_end));
	    return 0;
	}
    }
    return chromosome_cmp;
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for fields
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-15  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_fields(bed_feature_t *bed_feature, unsigned fields)

{
    if ( (fields < 3) || (fields > 9) )
	return BIO_OUT_OF_RANGE;
    else
    {
	bed_feature->fields = fields;
	return BIO_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for chromosome
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-15  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_chromosome(bed_feature_t *bed_feature, char *chromosome)

{
    if ( chromosome == NULL )
	return BIO_INVALID_DATA;
    else
    {
	strlcpy(bed_feature->chromosome, chromosome, BIO_CHROMOSOME_MAX_CHARS);
	return BIO_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for start_pos_str
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-15  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_start_pos_str(bed_feature_t *bed_feature, char *start_pos_str)

{
    char    *end;
    long long   n;
    
    if ( start_pos_str == NULL )
	return BIO_INVALID_DATA;
    else
    {
	n = strtoull(start_pos_str, &end, 10);
	if ( *end == '\0' )
	{
	    strlcpy(bed_feature->start_pos_str, start_pos_str,
		    BIO_POSITION_MAX_DIGITS);
	    bed_feature->start_pos = n;
	    return BIO_DATA_OK;
	}
	else
	    return BIO_INVALID_DATA;
    }
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for start_pos
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-15  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_start_pos(bed_feature_t *bed_feature, uint64_t start_pos)

{
    bed_feature->start_pos = start_pos;
    snprintf(bed_feature->start_pos_str, BIO_POSITION_MAX_DIGITS, "%" PRIu64,
	    start_pos);
    return BIO_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for end_pos_str
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-15  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_end_pos_str(bed_feature_t *bed_feature, char *end_pos_str)

{
    char    *end;
    long long   n;
    
    if ( end_pos_str == NULL )
	return BIO_INVALID_DATA;
    else
    {
	n = strtoull(end_pos_str, &end, 10);
	if ( *end == '\0' )
	{
	    strlcpy(end_pos_str, end_pos_str,
		    BIO_POSITION_MAX_DIGITS);
	    bed_feature->end_pos = n;
	    return BIO_DATA_OK;
	}
	else
	    return BIO_INVALID_DATA;
    }
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for end_pos
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-15  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_end_pos(bed_feature_t *bed_feature, uint64_t end_pos)

{
    bed_feature->end_pos = end_pos;
    snprintf(bed_feature->end_pos_str, BIO_POSITION_MAX_DIGITS, "%" PRIu64,
	    end_pos);
    return BIO_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for name
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-15  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_name(bed_feature_t *bed_feature, char *name)

{
    if ( name == NULL )
	return BIO_INVALID_DATA;
    else
    {
	strlcpy(bed_feature->name, name, BED_NAME_MAX_CHARS);
	return BIO_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for score
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-28  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_score(bed_feature_t *feature, unsigned score)

{
    if ( score > 1000 )
    {
	fprintf(stderr, "bed_set_score(): Score must be 0 to 1000: %u\n",
		score);
	return BIO_INVALID_DATA;
    }
    feature->score = score;
    return BIO_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <bedio.h>
 *      -lbiolibc
 *
 *  Description:
 *      Mutator for strand
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-28  Jason Bacon Begin
 ***************************************************************************/

int     bed_set_strand(bed_feature_t *feature, int strand)

{
    if ( (strand != '+') && (strand != '-') && (strand != '.') )
    {
	fprintf(stderr, "bed_set_strand(): Strand must be '+' or '-': %c\n",
		strand);
	return BIO_INVALID_DATA;
    }
    feature->strand = strand;
    return BIO_DATA_OK;
}
