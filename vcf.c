#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <xtend/dsv.h>
#include <xtend/string.h>   // ltostrn()
#include <xtend/mem.h>
#include "vcf.h"
#include "biostring.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over header lines in VCF input stream, leaving the FILE
 *      structure pointing to the first character in the first line of data
 *      or the first character of the header line starting with #CHROM if
 *      one is present.  The header line is typically read using
 *      bl_vcf_get_sample_ids(3). The skipped header is copied to a temporary
 *      file whose FILE pointer is returned.
 *
 *  Arguments:
 *      vcf_stream  FILE pointer of VCF stream to be read
 *
 *  Returns:
 *      A FILE pointer to the temporary file with a copy of the header
 *
 *  See also:
 *      bl_vcf_get_sample_ids(3), bl_vcf_read_static_fields(3), bl_vcf_read_ss_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_skip_meta_data(FILE *vcf_stream, FILE **meta_stream)

{
    int     ch,
	    c,
	    count;
    char    start[6];

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like vcf-split to replicate the
     *  header in output files.
     */

    *meta_stream = tmpfile();
    
    while ( (ch = getc(vcf_stream)) == '#' )
    {
	count = fread(start, 1, 5, vcf_stream);
	
	// Put back "CHROM" or whatever was read regardless
	for (c = count - 1; c >= 0; --c)
	    ungetc(start[c], vcf_stream);
	// Something is seriously wrong if we don't find at least 5 chars
	if ( count != 5 )
	    return BL_READ_TRUNCATED;
	
	if ( memcmp(start, "CHROM", 5) == 0 )
	{
	    // After return, read should start with #CHROM
	    ungetc(ch, vcf_stream);
	    rewind(*meta_stream);
	    return BL_READ_OK;
	}
	else
	{
	    // No #CHROM, transfer entire line to temp file
	    putc('#', *meta_stream);
	    do
	    {
		ch = getc(vcf_stream);
		putc(ch, *meta_stream);
	    }   while ( (ch != '\n') && (ch != EOF) );
	    if ( ch == EOF )
	    {
		fprintf(stderr,
		    "bl_vcf_skip_header(): EOF reached reading meta-data.\n");
		rewind(*meta_stream);
		return BL_READ_TRUNCATED;
	    }
	}
    }
    
    fprintf(stderr, "bl_vcf_skip_header(): Warning: No #CHROM found in header.\n");
    rewind(*meta_stream);
    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Extract sample IDs from a VCF input header line.  This is typically
 *      done following bl_vcf_skip_header(3), which will leave the FILE
 *      pointer pointing to the beginning of the header line, if one is
 *      present.
 *
 *      The arguments first_col and last_col represent the first and
 *      last sample columns, both inclusive, from which sample IDs should
 *      be extracted.  A value of 1 represents the first sample column.
 *      This feature allows a VCF file with many columns to be processed
 *      in multiple stages.  For example, the vcf-split tool, based on
 *      biolibc, cannot efficiently process more than abou1 10,000 samples
 *      at once, since each sample requires an open output file.  A VCF
 *      with 150,000 samples can be processed in 15 separate passes.
 *
 *  Arguments:
 *      vcf_stream  FILE pointer to the VCF input stream
 *      sample_ids  Array if character pointers to receive sample IDs
 *      first_col   First column from which a sample ID should be saved
 *      last_col    Last column from which a sample ID should be saved
 *
 *  See also:
 *      bl_vcf_skip_header(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    bl_vcf_get_sample_ids(FILE *vcf_stream, char *sample_ids[],
			   size_t first_col, size_t last_col)

{
    size_t  c,
	    len;
    char    temp_sample_id[BL_VCF_ID_MAX_CHARS + 1];
    int     delimiter = 0;
    
    // Skip standard header tags to get to sample IDs
    for (c = 0; c < 9; ++c)
	tsv_skip_field(vcf_stream);
    
    // Skip sample IDs before first_col
    for (c = 1; c < first_col; ++c)
	tsv_skip_field(vcf_stream);
    
    for (; (c <= last_col) &&
	   (delimiter = tsv_read_field(vcf_stream, temp_sample_id,
				     BL_VCF_ID_MAX_CHARS, &len)) != EOF; ++c)
    {
	sample_ids[c - first_col] = strdup(temp_sample_id);
	// fprintf(stderr, "'%s'\n", temp_sample_id);
    }
    
    if ( delimiter == 0 )
    {
	fprintf(stderr, "Reached last_col before reading any sample IDs.\n");
	fprintf(stderr, "Check your first_col and last_col values.\n");
	exit(EX_DATAERR);
    }
    
    // Skip any remaining fields after last_col
    if ( delimiter != '\n' )
	tsv_skip_rest_of_line(vcf_stream);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read static fields (columns 1 to 9) from one line of a VCF file.
 *      This function does not read any of the sample data in columns 10
 *      and on.  Samples can be read using a loop with tsv_read_field(3).
 *
 *      If field_mask is not BL_VCF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in bed_feature.
 *      Possible mask values are:
 *
 *      BL_VCF_FIELD_ALL
 *      BL_VCF_FIELD_CHROM
 *      BL_VCF_FIELD_POS
 *      BL_VCF_FIELD_ID
 *      BL_VCF_FIELD_REF
 *      BL_VCF_FIELD_ALT
 *      BL_VCF_FIELD_QUAL
 *      BL_VCF_FIELD_FILTER
 *      BL_VCF_FIELD_INFO
 *      BL_VCF_FIELD_FORMAT
 *
 *  Arguments:
 *      vcf_stream  FILE stream for VCF input
 *      vcf_call    Pointer to bl_vcf_t structure to receive fields
 *      field_mask  Bit mask indicating which fields should be stored
 *
 *  Returns:
 *      BL_READ_OK upon success
 *      BL_READ_TRUNCATED if EOF is encountered while reading a call
 *      BL_READ_EOF if EOF is encountered between calls as it should be
 *
 *  Examples:
 *      FILE        *stream;
 *      bl_vcf_t  vcf_call;
 *      char        sample_data[MAX_CHARS + 1];
 *      size_t      len;
 *
 *      bl_vcf_read_static_fields(stream, &vcf_call, BL_VCF_FIELD_ALL);
 *      while ( tsv_read_field(stream, sample_data, MAX_CHARS, &len) != '\n' )
 *      {
 *          ...
 *      }
 *
 *  See also:
 *      bl_vcf_write_static_fields(3), bl_vcf_read_ss_call(3), bl_vcf_write_ss_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_read_static_fields(bl_vcf_t *vcf_call, FILE *vcf_stream, 
	    vcf_field_mask_t field_mask)

{
    char    *end,
	    pos_str[BL_POSITION_MAX_DIGITS + 1];
    size_t  len;
    int     status;
    
    vcf_call->ref_count = vcf_call->alt_count = vcf_call->other_count = 0;
    
    // Chromosome
    if ( field_mask & BL_VCF_FIELD_CHROM )
	status = tsv_read_field(vcf_stream, vcf_call->chrom,
			BL_CHROM_MAX_CHARS, &len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	strlcpy(vcf_call->chrom, ".", 2);
    }
    if ( status == EOF )
    {
	// fputs("bl_vcf_read_static_fields(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return BL_READ_EOF;
    }
    
    // Call position
    if ( field_mask & BL_VCF_FIELD_POS )
	status = tsv_read_field(vcf_stream, pos_str,
			BL_POSITION_MAX_DIGITS, &len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	strlcpy(pos_str, "0", 2);
    }
    if ( status == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading POS: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	vcf_call->pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_vcf_read_static_fields(): Invalid call position: %s\n",
		    pos_str);
	    return BL_READ_TRUNCATED;
	}
    }
    
    // ID
    if ( field_mask & BL_VCF_FIELD_ID )
	status = tsv_read_field(vcf_stream, vcf_call->id,
			BL_VCF_ID_MAX_CHARS, &len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	strlcpy(vcf_call->id, ".", 2);
    }
    if ( status == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading ID.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Ref
    if ( field_mask & BL_VCF_FIELD_REF )
	status = tsv_read_field(vcf_stream, vcf_call->ref,
			BL_VCF_REF_MAX_CHARS, &len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	strlcpy(vcf_call->ref, ".", 2);
    }
    if ( status == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading REF.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Alt
    if ( field_mask & BL_VCF_FIELD_ALT )
	status = tsv_read_field(vcf_stream, vcf_call->alt,
		   BL_VCF_ALT_MAX_CHARS, &len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	strlcpy(vcf_call->alt, ".", 2);
    }
    if ( status == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading ALT.\n");
	return BL_READ_TRUNCATED;
    }

    // Qual
    if ( field_mask & BL_VCF_FIELD_QUAL )
	status = tsv_read_field(vcf_stream, vcf_call->qual,
		   BL_VCF_QUAL_MAX_CHARS, &len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	strlcpy(vcf_call->qual, ".", 2);
    }
    if ( status == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading QUAL.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Filter
    if ( field_mask & BL_VCF_FIELD_FILTER )
	status = tsv_read_field(vcf_stream, vcf_call->filter,
		   BL_VCF_FILTER_MAX_CHARS, &len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	strlcpy(vcf_call->filter, ".", 2);
    }
    if ( status == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading FILTER.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Info
    if ( field_mask & BL_VCF_FIELD_INFO )
	status = tsv_read_field(vcf_stream, vcf_call->info,
		   vcf_call->info_max, &vcf_call->info_len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	vcf_call->info = strdup(".");
	vcf_call->info_len = 1;
    }
    if ( status == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading INFO.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Format
    if ( field_mask & BL_VCF_FIELD_FORMAT )
	status = tsv_read_field(vcf_stream, vcf_call->format,
		   vcf_call->format_max, &vcf_call->format_len);
    else
    {
	status = tsv_skip_field(vcf_stream);
	strlcpy(vcf_call->format, ".", 2);
	vcf_call->format_len = 1;
    }
    if ( status == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading FORMAT.\n");
	return BL_READ_TRUNCATED;
    }

    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a single-sample VCF call (static fields and one sample column).
 *      This should only be used with VCF inputs that have exactly one
 *      sample column.  For multisample VCFs, use bl_vcf_read_static_fields()
 *      followed by a loop to read the sample data.
 *
 *      If field_mask is not BL_VCF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in bed_feature.
 *      Possible mask values are:
 *
 *      BL_VCF_FIELD_ALL
 *      BL_VCF_FIELD_CHROM
 *      BL_VCF_FIELD_POS
 *      BL_VCF_FIELD_ID
 *      BL_VCF_FIELD_REF
 *      BL_VCF_FIELD_ALT
 *      BL_VCF_FIELD_QUAL
 *      BL_VCF_FIELD_FILTER
 *      BL_VCF_FIELD_INFO
 *      BL_VCF_FIELD_FORMAT
 *
 *  Arguments:
 *      vcf_stream  FILE pointer to VCF input stream
 *      vcf_call    bl_vcf_t structure to receive VCF data
 *      field_mask  Bit mask to indicate which fields to store
 *
 *  Returns:
 *      BL_READ_OK upon success
 *      BL_READ_TRUNCATED if EOF is encountered while reading a call
 *      BL_READ_EOF if EOF is encountered between calls as it should be
 *
 *  See also:
 *      bl_vcf_read_static_fields(3), tsv_read_field(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-11  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_read_ss_call(bl_vcf_t *vcf_call, FILE *vcf_stream,
	    vcf_field_mask_t field_mask)

{
    size_t  len;
    int     status;
    
    status = bl_vcf_read_static_fields(vcf_call, vcf_stream, field_mask);
    if ( status == BL_READ_OK )
    {
	if ( tsv_read_field(vcf_stream, vcf_call->single_sample,
			vcf_call->sample_max, &len) != EOF )
	    return BL_READ_OK;
	else
	{
	    fprintf(stderr, "bl_vcf_read_ss_call(): Got EOF reading sample.\n");
	    return BL_READ_TRUNCATED;
	}
    }
    else
	return status;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write static fields from one line of a single-entry VCF file.
 *      Does not write sample data.
 *
 *      If field_mask is not BL_VCF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate placeholder such as '.'
 *      rather than the actual data.  Possible mask values are:
 *
 *      BL_VCF_FIELD_ALL
 *      BL_VCF_FIELD_CHROM
 *      BL_VCF_FIELD_POS
 *      BL_VCF_FIELD_ID
 *      BL_VCF_FIELD_REF
 *      BL_VCF_FIELD_ALT
 *      BL_VCF_FIELD_QUAL
 *      BL_VCF_FIELD_FILTER
 *      BL_VCF_FIELD_INFO
 *      BL_VCF_FIELD_FORMAT
 *
 *  Arguments:
 *      vcf_stream  FILE pointer to the VCF output stream
 *      vcf_call    Pointer to the bl_vcf_t structure to output
 *      field_mask  Bit mask indicating which fields to output
 *
 *  Returns:
 *      The number of items output (as returned by fprintf())
 *
 *  See also:
 *      bl_vcf_read_static_fields(3), bl_vcf_write_ss_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_write_static_fields(bl_vcf_t *vcf_call, FILE *vcf_stream,
	    vcf_field_mask_t field_mask)

{
    char    *chrom = ".",
	    pos_str[BL_POSITION_MAX_DIGITS+1] = ".",
	    *id = ".",
	    *ref = ".",
	    *alt = ".",
	    *qual = ".",
	    *filter = ".",
	    *info = ".",
	    *format = ".";
    
    if ( field_mask & BL_VCF_FIELD_CHROM )
	chrom = vcf_call->chrom;
    if ( field_mask & BL_VCF_FIELD_POS )
	ltostrn(pos_str, vcf_call->pos, 10, BL_POSITION_MAX_DIGITS);
    if ( field_mask & BL_VCF_FIELD_ID )
	id = vcf_call->id;
    if ( field_mask & BL_VCF_FIELD_REF )
	ref = vcf_call->ref;
    if ( field_mask & BL_VCF_FIELD_ALT )
	alt = vcf_call->alt;
    if ( field_mask & BL_VCF_FIELD_QUAL )
	qual = vcf_call->qual;
    if ( field_mask & BL_VCF_FIELD_FILTER )
	filter = vcf_call->filter;
    if ( field_mask & BL_VCF_FIELD_INFO )
	info = vcf_call->info;
    if ( field_mask & BL_VCF_FIELD_FORMAT )
	format = vcf_call->format;
    
    return fprintf(vcf_stream,
	    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
	    chrom, pos_str,
	    id, ref, alt, 
	    qual, filter, info,
	    format);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a single-sample VCF call to vcf_stream.
 *      This should only be used with VCF calls that have exactly one
 *      sample column.  For multisample VCFs, use bl_vcf_write_static_fields()
 *      followed by a loop to write the sample data.
 *
 *      If field_mask is not BL_VCF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate placeholder such as '.'
 *      rather than the actual data.  Possible mask values are:
 *
 *      BL_VCF_FIELD_ALL
 *      BL_VCF_FIELD_CHROM
 *      BL_VCF_FIELD_POS
 *      BL_VCF_FIELD_ID
 *      BL_VCF_FIELD_REF
 *      BL_VCF_FIELD_ALT
 *      BL_VCF_FIELD_QUAL
 *      BL_VCF_FIELD_FILTER
 *      BL_VCF_FIELD_INFO
 *      BL_VCF_FIELD_FORMAT
 *
 *  Arguments:
 *      vcf_stream  FILE pointer to the VCF output stream
 *      vcf_call    Pointer to the bl_vcf_t structure to output
 *      field_mask  Bit mask indicating which fields to output
 *
 *  Returns:
 *      The number of items output (as returned by fprintf())
 *
 *  See also:
 *      bl_vcf_read_ss_call(3), bl_vcf_write_static_fields(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_write_ss_call(bl_vcf_t *vcf_call, FILE *vcf_stream,
	    vcf_field_mask_t field_mask)

{
    bl_vcf_write_static_fields(vcf_call, vcf_stream, field_mask);
    return fprintf(vcf_stream, "%s\n", vcf_call->single_sample);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Allocate an array for multiple samples in a VCF call.
 *
 *  Arguments:
 *      vcf_call    Pointer to a bl_vcf_t structure with multiple samples
 *      samples     The number of samples that must be accommodated
 *
 *  Returns:
 *      Address of the allocated array (NULL if malloc failed)
 *
 *  See also:
 *      bl_vcf_init(3), bl_vcf_read_static_fields(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

char    **bl_vcf_sample_alloc(bl_vcf_t *vcf_call, size_t samples)

{
    size_t  c;
    
    if ( (vcf_call->multi_samples =
	 (char **)xt_malloc(samples,
		    sizeof(*vcf_call->multi_samples))) != NULL )
    {
	for (c = 0; c < samples; ++c)
	{
	    if ( (vcf_call->multi_samples[c] =
		 (char *)xt_malloc(vcf_call->sample_max + 1,
				sizeof(*vcf_call->multi_samples[c]))) == NULL )
		return NULL;
	}
    }
    return vcf_call->multi_samples;
}


/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

#if 0
int     vcf_phred_add(bl_vcf_t *vcf_call, unsigned char score)

{
    if ( vcf_call->phreds == NULL )
    {
	// fprintf(stderr, "vcf_phred_add(): Allocating initial buffer.\n");
	if ( (vcf_call->phreds = xt_malloc(vcf_call->phred_buff_size,
				    sizeof(*vcf_call->phreds))) == NULL )
	{
	    fprintf(stderr, "vcf_phred_add(): Could not allocate phreds.\n");
	    exit(EX_UNAVAILABLE);
	}
    }
    
    // fprintf(stderr, "vcf_phred_add(): Adding '%c' at %zu\n", score, vcf_call->phred_count);
    vcf_call->phreds[vcf_call->phred_count++] = score;
    vcf_call->phreds[vcf_call->phred_count] = '\0';
    
    if ( vcf_call->phred_count == vcf_call->phred_buff_size )
    {
	vcf_call->phred_buff_size *= 2;
	if ( (vcf_call->phreds = xt_realloc(vcf_call->phreds,
		    vcf_call->phred_buff_size,
		    sizeof(*vcf_call->phreds))) == NULL )
	{
	    fprintf(stderr, "vcf_phred_add(): Could not reallocate phreds.\n");
	    exit(EX_UNAVAILABLE);
	}
    }
    return BL_READ_OK;
}


/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    vcf_phred_blank(bl_vcf_t *vcf_call)

{
    memcpy(vcf_call->phreds, "z", 2);
    vcf_call->phred_count = 0;
}

    
/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    vcf_phred_free(bl_vcf_t *vcf_call)

{
    if ( vcf_call->phreds != NULL )
    {
	free(vcf_call->phreds);
	vcf_call->phreds = NULL;
	vcf_call->phred_buff_size = BL_VCF_PHRED_BUFF_SIZE;
    }
    vcf_phred_blank(vcf_call);
}
#endif


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free all memory associated with a VCF call.
 *
 *  Arguments:
 *      vcf_call    Pointer to the bl_vcf_t structure to free.
 *
 *  See also:
 *      bl_vcf_init(3), bl_vcf_sample_alloc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    bl_vcf_free(bl_vcf_t *vcf_call)

{
    free(vcf_call->info);
    free(vcf_call->format);
    free(vcf_call->single_sample);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_vcf_t structure, allocating default buffer
 *      sizes for some fields.
 *
 *  Arguments:
 *      vcf_call    Pointer to the bl_vcf_t structure to initialize
 *      info_max    Maximum size of INFO field in bytes
 *      format_max  Maximum size of FORMAT field in bytes
 *      sample_max  Maxixum size of SAMPLE field in bytes
 *
 *  See also:
 *      bl_vcf_free(3), vcf_read_call(3), bl_vcf_sample_alloc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    bl_vcf_init(bl_vcf_t *vcf_call,
	    size_t info_max, size_t format_max, size_t sample_max)

{
    vcf_call->chrom[0] = '\0';
    vcf_call->id[0] = '\0';
    vcf_call->ref[0] = '\0';
    vcf_call->alt[0] = '\0';
    vcf_call->qual[0] = '\0';
    vcf_call->filter[0] = '\0';
    vcf_call->pos = 0;
    vcf_call->info_len = 0;
    vcf_call->ref_count = 0;
    vcf_call->alt_count = 0;
    vcf_call->other_count = 0;

    if ( (vcf_call->info = xt_malloc(info_max + 1,
	    sizeof(*vcf_call->info))) == NULL )
    {
	fprintf(stderr, "bl_vcf_init(): Could not allocate info field.\n");
	exit(EX_UNAVAILABLE);
    }
    if ( (vcf_call->format = xt_malloc(format_max + 1,
	    sizeof(*vcf_call->format))) == NULL )
    {
	fprintf(stderr, "bl_vcf_init(): Could not allocate format field.\n");
	exit(EX_UNAVAILABLE);
    }
    if ( (vcf_call->single_sample = xt_malloc(sample_max + 1,
	    sizeof(*vcf_call->single_sample))) == NULL )
    {
	fprintf(stderr, "bl_vcf_init(): Could not allocate single_sample field.\n");
	exit(EX_UNAVAILABLE);
    }
    
    vcf_call->info_max = info_max;
    vcf_call->format_max = format_max;
    vcf_call->sample_max = sample_max;
    
    vcf_call->info[0] = '\0';
    vcf_call->format[0] = '\0';
    vcf_call->single_sample[0] = '\0';
    vcf_call->multi_samples = NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Convert a comma-separated list of VCF fields to a field bit mask
 *
 *  Arguments:
 *      spec    Character string containing comma-separated field list
 *
 *  Returns:
 *      A vcf_field_mask_t value with bits set for specified fields
 *
 *  See also:
 *      vcf_read_call(3), vcf_write_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

vcf_field_mask_t    bl_vcf_parse_field_spec(char *spec)

{
    vcf_field_mask_t    field_mask;
    char            *field_name;
    
    if ( strcmp(spec, "all") == 0 )
    {
	field_mask = BL_VCF_FIELD_ALL;
    }
    else
    {
	field_mask = 0x0;
	while ((field_name = strsep(&spec, ",")) != NULL)
	{
	    if ( strcmp(field_name, "chrom") == 0 )
		field_mask |= BL_VCF_FIELD_CHROM;
	    else if ( strcmp(field_name, "pos") == 0 )
		field_mask |= BL_VCF_FIELD_POS;
	    else if ( strcmp(field_name, "id") == 0 )
		field_mask |= BL_VCF_FIELD_ID;
	    else if ( strcmp(field_name, "ref") == 0 )
		field_mask |= BL_VCF_FIELD_REF;
	    else if ( strcmp(field_name, "alt") == 0 )
		field_mask |= BL_VCF_FIELD_ALT;
	    else if ( strcmp(field_name, "qual") == 0 )
		field_mask |= BL_VCF_FIELD_QUAL;
	    else if ( strcmp(field_name, "filter") == 0 )
		field_mask |= BL_VCF_FIELD_FILTER;
	    else if ( strcmp(field_name, "info") == 0 )
		field_mask |= BL_VCF_FIELD_INFO;
	    else if ( strcmp(field_name, "format") == 0 )
		field_mask |= BL_VCF_FIELD_FORMAT;
	    else
		field_mask = BL_VCF_FIELD_ERROR;
	}
    }
    return field_mask;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Determine if a VCF call is within a SAM alignment, i.e. on the
 *      same chrom and between the start and end positions of the
 *      alignment.
 *
 *  Arguments:
 *      vcf_call    Pointer to bl_vcf_t structure containing VCF call
 *      sam_alignment   Pointer to bl_sam_t structure containing alignment
 *
 *  Returns:
 *      true if the call is within the alignment
 *      false otherwise
 *
 *  See also:
 *      bl_vcf_call_downstream_of_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/


bool    bl_vcf_call_in_alignment(bl_vcf_t *vcf_call, bl_sam_t *sam_alignment)

{
    if ( (strcmp(BL_VCF_CHROM(vcf_call), BL_SAM_RNAME(sam_alignment)) == 0) &&
	 (BL_VCF_POS(vcf_call) >= BL_SAM_POS(sam_alignment)) &&
	 (BL_VCF_POS(vcf_call) <
	    BL_SAM_POS(sam_alignment) + BL_SAM_SEQ_LEN(sam_alignment)) )
	return true;
    else
	return false;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Determine if a VCF call is downstream of a SAM alignment.
 *      For the purpose of this function, this could mean on the same
 *      chrom and higher position, or on a later chrom.
 *
 *  Arguments:
 *      vcf_call    Pointer to bl_vcf_t structure containing VCF call
 *      sam_alignment   Pointer to bl_sam_t structure containing alignment
 *
 *  Returns:
 *      true if the call is downstream of the alignment
 *      false otherwise
 *
 *  See also:
 *      bl_vcf_call_in_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    bl_vcf_call_downstream_of_alignment(bl_vcf_t *vcf_call,
	    bl_sam_t *alignment)

{
    /*fprintf(stderr, "bl_vcf_call_downstream_of_alignment(): %s,%zu,%zu %s,%zu\n",
	    BL_SAM_RNAME(sam_alignment),BL_SAM_POS(sam_alignment),
	    BL_SAM_SEQ_LEN(sam_alignment),
	    BL_VCF_CHROM(vcf_call),BL_VCF_POS(vcf_call));*/
    if ( (BL_SAM_POS(alignment) + BL_SAM_SEQ_LEN(alignment) <= BL_VCF_POS(vcf_call)) &&
	  (strcmp(BL_SAM_RNAME(alignment), BL_VCF_CHROM(vcf_call)) == 0) )
	return true;
    else if ( bl_chrom_name_cmp(BL_SAM_RNAME(alignment), BL_VCF_CHROM(vcf_call)) < 0 )
	return true;
    else
	return false;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Report VCF input sort error and terminate the process.
 *
 *  Arguments:
 *      vcf_call        Pointer to bl_vcf_t structure with latest call
 *      previous_chrom  Chromosome of previous VCF call
 *      previous_pos    Position of previous VCF call
 *
 *  Returns:
 *      Does not return
 *
 *  See also:
 *      vcf_read_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_vcf_call_out_of_order(bl_vcf_t *vcf_call,
	    char *previous_chrom, uint64_t previous_pos)

{
    fprintf(stderr, "ad2vcf: Error: VCF input must be sorted by chrom and then position.\n");
    fprintf(stderr, "Found %s,%" PRIu64 " after %s,%" PRIu64 ".\n",
	    BL_VCF_CHROM(vcf_call), BL_VCF_POS(vcf_call),
	    previous_chrom, previous_pos);
    exit(EX_DATAERR);
}
