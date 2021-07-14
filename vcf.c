#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <xtend.h>
#include "vcf.h"
#include "biostring.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *      Skip over header lines in VCF input stream, leaving the FILE
 *      structure pointing to the first character in the first line of data
 *      or the first character of the header line starting with #CHROM if
 *      one is present.  The header line is typically read using
 *      vcf_get_sample_ids(3). The skipped header is copied to a temporary
 *      file whose FILE pointer is returned.
 *
 *  Arguments:
 *      vcf_stream: FILE pointer of VCF stream to be read
 *
 *  Returns:
 *      A FILE pointer to the temporary file with a copy of the header
 *
 *  See also:
 *      vcf_get_sample_ids(3), vcf_read_static_fields(3), vcf_read_ss_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

FILE    *vcf_skip_header(FILE *vcf_stream)

{
    char    start[7] = "xxxxxx";
    size_t  count;
    int     ch;
    FILE    *header_stream = tmpfile();

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like vcf-split to replicate the
     *  header in output files.
     */
    
    while ( ((count=fread(start, 6, 1, vcf_stream)) == 1) && 
	    (memcmp(start, "#CHROM", 6) != 0) )
    {
	fwrite(start, 6, 1, header_stream);
	do
	{
	    ch = getc(vcf_stream);
	    putc(ch, header_stream);
	}   while ( ch != '\n' );
    }
    
    // puts(start);
    if ( count == 0 )
    {
	fprintf(stderr, "vcf_skip_header(): No #CHROM header found.\n");
	exit(EX_DATAERR);
    }
    rewind(header_stream);
    return header_stream;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *      Extract sample IDs from a VCF input header line.  This is typically
 *      done following vcf_skip_header(3), which will leave the FILE
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
 *      vcf_stream: FILE pointer to the VCF input stream
 *      sample_ids: Array if character pointers to receive sample IDs
 *      first_col:  First column from which a sample ID should be saved
 *      last_col:   Last column from which a sample ID should be saved
 *
 *  See also:
 *      vcf_skip_header(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    vcf_get_sample_ids(FILE *vcf_stream,
			   char *sample_ids[],
			   size_t first_col, size_t last_col)

{
    size_t  c,
	    len;
    char    temp_sample_id[VCF_ID_MAX_CHARS + 1];
    int     delimiter = 0;
    
    // Skip standard header tags to get to sample IDs
    for (c = 0; c < 9; ++c)
	tsv_skip_field(vcf_stream);
    
    // Skip sample IDs before first_col
    for (c = 1; c < first_col; ++c)
	tsv_skip_field(vcf_stream);
    
    for (; (c <= last_col) &&
	   (delimiter = tsv_read_field(vcf_stream, temp_sample_id,
				     VCF_ID_MAX_CHARS, &len)) != EOF; ++c)
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
 *      -lbiolibc
 *
 *  Description:
 *      Read static fields from one line of a single-entry VCF file.
 *      Does not read sample data.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     vcf_read_static_fields(FILE *vcf_stream, vcf_call_t *vcf_call,
			       vcf_field_mask_t field_mask)

{
    char    *end;
    size_t  len;
    
    vcf_call->ref_count = vcf_call->alt_count = vcf_call->other_count = 0;
    
    // Chromosome
    if ( tsv_read_field(vcf_stream, vcf_call->chromosome,
			BIO_CHROMOSOME_MAX_CHARS, &len) == EOF )
    {
	// fputs("vcf_read_static_fields(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return BIO_READ_EOF;
    }
    
    // Call position
    if ( tsv_read_field(vcf_stream, vcf_call->pos_str,
			BIO_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading POS: %s.\n",
		vcf_call->pos_str);
	return BIO_READ_TRUNCATED;
    }
    else
    {
	vcf_call->pos = strtoul(vcf_call->pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "vcf_read_static_fields(): Invalid call position: %s\n",
		    vcf_call->pos_str);
	    return BIO_READ_TRUNCATED;
	}
    }
    
    // ID
    if ( tsv_read_field(vcf_stream, vcf_call->id,
			VCF_ID_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading ID.\n");
	return BIO_READ_TRUNCATED;
    }
    
    // Ref
    if ( tsv_read_field(vcf_stream, vcf_call->ref,
			VCF_REF_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading REF.\n");
	return BIO_READ_TRUNCATED;
    }
    
    // Alt
    if ( tsv_read_field(vcf_stream, vcf_call->alt,
		   VCF_ALT_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading ALT.\n");
	return BIO_READ_TRUNCATED;
    }

    // Qual
    if ( tsv_read_field(vcf_stream, vcf_call->quality,
		   VCF_QUALITY_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading QUAL.\n");
	return BIO_READ_TRUNCATED;
    }
    
    // Filter
    if ( tsv_read_field(vcf_stream, vcf_call->filter,
		   VCF_FILTER_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading FILTER.\n");
	return BIO_READ_TRUNCATED;
    }
    
    // Info
    if ( tsv_read_field(vcf_stream, vcf_call->info,
		   vcf_call->info_max, &vcf_call->info_len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading INFO.\n");
	return BIO_READ_TRUNCATED;
    }
    
    // Format
    if ( tsv_read_field(vcf_stream, vcf_call->format,
		   vcf_call->format_max, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading FORMAT.\n");
	return BIO_READ_TRUNCATED;
    }

    return BIO_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read a single-sample VCF call.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-11  Jason Bacon Begin
 ***************************************************************************/

int     vcf_read_ss_call(FILE *vcf_stream, vcf_call_t *vcf_call,
			 vcf_field_mask_t field_mask)

{
    size_t  len;
    int     status;
    
    status = vcf_read_static_fields(vcf_stream, vcf_call, field_mask);
    if ( status == BIO_READ_OK )
    {
	if ( tsv_read_field(vcf_stream, vcf_call->single_sample,
			vcf_call->sample_max, &len) != EOF )
	    return BIO_READ_OK;
	else
	{
	    fprintf(stderr, "vcf_read_ss_call(): Got EOF reading sample.\n");
	    return BIO_READ_TRUNCATED;
	}
    }
    else
	return status;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *      Write static fields from one line of a single-entry VCF file.
 *      Does not write sample data.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     vcf_write_static_fields(FILE *vcf_stream, vcf_call_t *vcf_call,
				vcf_field_mask_t field_mask)

{
    char    *chromosome = ".",
	    *pos_str = ".",
	    *id = ".",
	    *ref = ".",
	    *alt = ".",
	    *quality = ".",
	    *filter = ".",
	    *info = ".",
	    *format = ".";
    
    if ( field_mask & VCF_FIELD_CHROM )
	chromosome = vcf_call->chromosome;
    if ( field_mask & VCF_FIELD_POS )
	pos_str = vcf_call->pos_str;
    if ( field_mask & VCF_FIELD_ID )
	id = vcf_call->id;
    if ( field_mask & VCF_FIELD_REF )
	ref = vcf_call->ref;
    if ( field_mask & VCF_FIELD_ALT )
	alt = vcf_call->alt;
    if ( field_mask & VCF_FIELD_QUAL )
	quality = vcf_call->quality;
    if ( field_mask & VCF_FIELD_FILTER )
	filter = vcf_call->filter;
    if ( field_mask & VCF_FIELD_INFO )
	info = vcf_call->info;
    if ( field_mask & VCF_FIELD_FORMAT )
	format = vcf_call->format;
    
    return fprintf(vcf_stream,
	    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
	    chromosome, pos_str,
	    id, ref, alt, 
	    quality, filter, info,
	    format);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *      Write a single-sample VCF call to vcf_stream.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     vcf_write_ss_call(FILE *vcf_stream, vcf_call_t *vcf_call,
			  vcf_field_mask_t field_mask)

{
    vcf_write_static_fields(vcf_stream, vcf_call, field_mask);
    fprintf(vcf_stream, "%s\n", vcf_call->single_sample);
    return 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

char    **vcf_sample_alloc(vcf_call_t *vcf_call, size_t samples)

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
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

#if 0
int     vcf_phred_add(vcf_call_t *vcf_call, unsigned char score)

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
    return BIO_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    vcf_phred_blank(vcf_call_t *vcf_call)

{
    memcpy(vcf_call->phreds, "z", 2);
    vcf_call->phred_count = 0;
}

    
/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    vcf_phred_free(vcf_call_t *vcf_call)

{
    if ( vcf_call->phreds != NULL )
    {
	free(vcf_call->phreds);
	vcf_call->phreds = NULL;
	vcf_call->phred_buff_size = VCF_PHRED_BUFF_SIZE;
    }
    vcf_phred_blank(vcf_call);
}
#endif


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    vcf_call_free(vcf_call_t *vcf_call)

{
    free(vcf_call->info);
    free(vcf_call->format);
    free(vcf_call->single_sample);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    vcf_call_init(vcf_call_t *vcf_call,
		      size_t info_max, size_t format_max, size_t sample_max)

{
    vcf_call->chromosome[0] = '\0';
    vcf_call->pos_str[0] = '\0';
    vcf_call->id[0] = '\0';
    vcf_call->ref[0] = '\0';
    vcf_call->alt[0] = '\0';
    vcf_call->quality[0] = '\0';
    vcf_call->filter[0] = '\0';
    vcf_call->pos = 0;
    vcf_call->info_len = 0;
    vcf_call->ref_count = 0;
    vcf_call->alt_count = 0;
    vcf_call->other_count = 0;

    if ( (vcf_call->info = xt_malloc(info_max + 1,
	    sizeof(*vcf_call->info))) == NULL )
    {
	fprintf(stderr, "vcf_call_init(): Could not allocate info field.\n");
	exit(EX_UNAVAILABLE);
    }
    if ( (vcf_call->format = xt_malloc(format_max + 1,
	    sizeof(*vcf_call->format))) == NULL )
    {
	fprintf(stderr, "vcf_call_init(): Could not allocate format field.\n");
	exit(EX_UNAVAILABLE);
    }
    if ( (vcf_call->single_sample = xt_malloc(sample_max + 1,
	    sizeof(*vcf_call->single_sample))) == NULL )
    {
	fprintf(stderr, "vcf_call_init(): Could not allocate single_sample field.\n");
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
 *      -lbiolibc
 *
 *  Description:
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

vcf_field_mask_t    vcf_parse_field_spec(char *spec)

{
    vcf_field_mask_t    field_mask = VCF_FIELD_ALL;
    char            *field_name;
    
    if ( strcmp("spec", "all") != 0 )
    {
	while ((field_name = strsep(&spec, ",")) != NULL)
	{
	    if ( strcmp(field_name, "chrom") == 0 )
		field_mask |= VCF_FIELD_CHROM;
	    else if ( strcmp(field_name, "pos") == 0 )
		field_mask |= VCF_FIELD_POS;
	    else if ( strcmp(field_name, "id") == 0 )
		field_mask |= VCF_FIELD_ID;
	    else if ( strcmp(field_name, "ref") == 0 )
		field_mask |= VCF_FIELD_REF;
	    else if ( strcmp(field_name, "alt") == 0 )
		field_mask |= VCF_FIELD_ALT;
	    else if ( strcmp(field_name, "qual") == 0 )
		field_mask |= VCF_FIELD_QUAL;
	    else if ( strcmp(field_name, "filter") == 0 )
		field_mask |= VCF_FIELD_FILTER;
	    else if ( strcmp(field_name, "info") == 0 )
		field_mask |= VCF_FIELD_INFO;
	    else if ( strcmp(field_name, "format") == 0 )
		field_mask |= VCF_FIELD_FORMAT;
	    else
		field_mask = VCF_FIELD_ERROR;
	}
    }
    return field_mask;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *      Determine whether a VCF call is within a SAM alignment.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/


bool    vcf_call_in_alignment(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment)

{
    if ( (strcmp(VCF_CHROMOSOME(vcf_call), SAM_RNAME(sam_alignment)) == 0) &&
	 (VCF_POS(vcf_call) >= SAM_POS(sam_alignment)) &&
	 (VCF_POS(vcf_call) <
	    SAM_POS(sam_alignment) + SAM_SEQ_LEN(sam_alignment)) )
	return true;
    else
	return false;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *      Determine SAM alignment is completely upstream of a VCF call position,
 *      i.e. not overlapping and at a lower position or chromosome.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    vcf_call_downstream_of_alignment(vcf_call_t *vcf_call, sam_alignment_t *alignment)

{
    /*fprintf(stderr, "vcf_call_downstream_of_alignment(): %s,%zu,%zu %s,%zu\n",
	    SAM_RNAME(sam_alignment),SAM_POS(sam_alignment),
	    SAM_SEQ_LEN(sam_alignment),
	    VCF_CHROMOSOME(vcf_call),VCF_POS(vcf_call));*/
    if ( (SAM_POS(alignment) + SAM_SEQ_LEN(alignment) <= VCF_POS(vcf_call)) &&
	  (strcmp(SAM_RNAME(alignment), VCF_CHROMOSOME(vcf_call)) == 0) )
	return true;
    else if ( chromosome_name_cmp(SAM_RNAME(alignment), VCF_CHROMOSOME(vcf_call)) < 0 )
	return true;
    else
	return false;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc
 *
 *  Description:
 *      Explain VCF input sort error and exit.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    vcf_out_of_order(vcf_call_t *vcf_call,
			 char *previous_chromosome, size_t previous_pos)

{
    fprintf(stderr, "ad2vcf: Error: VCF input must be sorted by chromosome and then position.\n");
    fprintf(stderr, "Found %s,%zu after %s,%zu.\n",
	    VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call),
	    previous_chromosome, previous_pos);
    exit(EX_DATAERR);
}

