#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include "vcfio.h"

/***************************************************************************
 *  Description:
 *      Skip over header lines in VCF input stream.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    vcf_skip_header(FILE *vcf_stream)

{
    char    start[7] = "xxxxxx";
    size_t  count;

    while ( ((count=fread(start, 6, 1, vcf_stream)) == 1) && 
	    (memcmp(start, "#CHROM", 6) != 0) )
	tsv_skip_rest_of_line(vcf_stream);
    
    // puts(start);
    if ( count == 0 )
    {
	fprintf(stderr, "vcf_skip_header(): No #CHROM header found.\n");
	exit(EX_DATAERR);
    }
}


/***************************************************************************
 *  Description:
 *      Extract sample IDs from input header line.
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
 *  Description:
 *      Read static fields from one line of a single-entry VCF file.
 *      Does not read sample data.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     vcf_read_static_fields(FILE *vcf_stream, vcf_call_t *vcf_call)

{
    char    *end;
    size_t  len;
    
    vcf_call->ref_count = vcf_call->alt_count = vcf_call->other_count = 0;
    
    // Chromosome
    if ( tsv_read_field(vcf_stream, vcf_call->chromosome,
			VCF_CHROMOSOME_MAX_CHARS, &len) == EOF )
    {
	fputs("vcf_read_static_fields(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return VCF_READ_EOF;
    }
    
    // Call position
    if ( tsv_read_field(vcf_stream, vcf_call->pos_str,
			VCF_POSITION_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading POS: %s.\n",
		vcf_call->pos_str);
	return VCF_READ_TRUNCATED;
    }
    else
    {
	vcf_call->pos = strtoul(vcf_call->pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "vcf_read_static_fields(): Invalid call position: %s\n",
		    vcf_call->pos_str);
	    return VCF_READ_TRUNCATED;
	}
    }
    
    // ID
    if ( tsv_read_field(vcf_stream, vcf_call->id,
			VCF_ID_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading ID.\n");
	return VCF_READ_TRUNCATED;
    }
    
    // Ref
    if ( tsv_read_field(vcf_stream, vcf_call->ref,
			VCF_REF_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading REF.\n");
	return VCF_READ_TRUNCATED;
    }
    
    // Alt
    if ( tsv_read_field(vcf_stream, vcf_call->alt,
		   VCF_ALT_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading ALT.\n");
	return VCF_READ_TRUNCATED;
    }

    // Qual
    if ( tsv_read_field(vcf_stream, vcf_call->quality,
		   VCF_QUALITY_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading QUAL.\n");
	return VCF_READ_TRUNCATED;
    }
    
    // Filter
    if ( tsv_read_field(vcf_stream, vcf_call->filter,
		   VCF_FILTER_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading FILTER.\n");
	return VCF_READ_TRUNCATED;
    }
    
    // Info
    if ( tsv_read_field(vcf_stream, vcf_call->info,
		   vcf_call->info_max, &vcf_call->info_len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading INFO.\n");
	return VCF_READ_TRUNCATED;
    }
    
    // Format
    if ( tsv_read_field(vcf_stream, vcf_call->format,
		   vcf_call->format_max, &len) == EOF )
    {
	fprintf(stderr, "vcf_read_static_fields(): Got EOF reading FORMAT.\n");
	return VCF_READ_TRUNCATED;
    }

    return VCF_OK;
}


/***************************************************************************
 *  Description:
 *      Read a single-sample VCF call.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-11  Jason Bacon Begin
 ***************************************************************************/

int     vcf_read_ss_call(FILE *vcf_stream, vcf_call_t *vcf_call)

{
    size_t  len;
    int     status;
    
    status = vcf_read_static_fields(vcf_stream, vcf_call);
    if ( status == VCF_OK )
    {
	if ( tsv_read_field(vcf_stream, vcf_call->single_sample,
			vcf_call->sample_max, &len) != EOF )
	    return VCF_OK;
	else
	{
	    fprintf(stderr, "vcf_read_ss_call(): Got EOF reading sample.\n");
	    return VCF_READ_TRUNCATED;
	}
    }
    else
	return status;
}


/***************************************************************************
 *  Description:
 *      Write static fields from one line of a single-entry VCF file.
 *      Does not write sample data.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     vcf_write_static_fields(FILE *vcf_stream, vcf_call_t *vcf_call)

{
    return fprintf(vcf_stream,
	    "%s\t%s/%zu\t%s\t%s\t%s\t%s\t%s\t%s/%zu\t%s\n",
	    vcf_call->chromosome, vcf_call->pos_str, vcf_call->pos,
	    vcf_call->id, vcf_call->ref, vcf_call->alt, 
	    vcf_call->quality, vcf_call->filter, vcf_call->info,
	    vcf_call->info_len, vcf_call->format);
}


/***************************************************************************
 *  Description:
 *      Write a single-sample VCF call to vcf_stream.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     vcf_write_ss_call(FILE *vcf_stream, vcf_call_t *vcf_call)

{
    vcf_write_static_fields(vcf_stream, vcf_call);
    fprintf(vcf_stream, "%s\n", vcf_call->single_sample);
    return 0;
}


char    **vcf_sample_alloc(vcf_call_t *vcf_call, size_t samples)

{
    size_t  c;
    
    if ( (vcf_call->multi_samples =
	 (char **)malloc(samples * sizeof(char *))) != NULL )
    {
	for (c = 0; c < samples; ++c)
	{
	    if ( (vcf_call->multi_samples[c] =
		 (char *)malloc(vcf_call->sample_max + 1)) == NULL )
		return NULL;
	}
    }
    return vcf_call->multi_samples;
}


#if 0
int     vcf_phred_add(vcf_call_t *vcf_call, unsigned char score)

{
    if ( vcf_call->phreds == NULL )
    {
	// fprintf(stderr, "vcf_phred_add(): Allocating initial buffer.\n");
	if ( (vcf_call->phreds = malloc(vcf_call->phred_buff_size)) == NULL )
	{
	    fprintf(stderr, "vcf_phred_add(): malloc() failure.\n");
	    exit(EX_UNAVAILABLE);
	}
    }
    
    // fprintf(stderr, "vcf_phred_add(): Adding '%c' at %zu\n", score, vcf_call->phred_count);
    vcf_call->phreds[vcf_call->phred_count++] = score;
    vcf_call->phreds[vcf_call->phred_count] = '\0';
    
    if ( vcf_call->phred_count == vcf_call->phred_buff_size )
    {
	vcf_call->phred_buff_size *= 2;
	if ( (vcf_call->phreds = realloc(vcf_call->phreds, vcf_call->phred_buff_size)) == NULL )
	{
	    fprintf(stderr, "vcf_phred_add(): realloc() failure.\n");
	    exit(EX_UNAVAILABLE);
	}
    }
    return VCF_OK;
}


void    vcf_phred_blank(vcf_call_t *vcf_call)

{
    memcpy(vcf_call->phreds, "z", 2);
    vcf_call->phred_count = 0;
}

    
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


void    vcf_call_free(vcf_call_t *vcf_call)

{
    free(vcf_call->info);
    free(vcf_call->format);
    free(vcf_call->single_sample);
}


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

    if ( (vcf_call->info = malloc(info_max)) == NULL )
    {
	fprintf(stderr, "vcf_call_init(): Could not allocate info field.\n");
	exit(EX_UNAVAILABLE);
    }
    if ( (vcf_call->format = malloc(format_max)) == NULL )
    {
	fprintf(stderr, "vcf_call_init(): Could not allocate format field.\n");
	exit(EX_UNAVAILABLE);
    }
    if ( (vcf_call->single_sample = malloc(sample_max)) == NULL )
    {
	fprintf(stderr, "vcf_call_init(): Could not allocate sample field.\n");
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


#ifdef __linux__
size_t  strlcpy(char *dest, const char *src, size_t len)

{
    char   *save_dest, *end;

    save_dest = dest;
    end = (char *)src + len - 1;
    while ((*src != '\0') && (src < end))
	*dest++ = *src++;
    *dest = '\0';
    return dest - save_dest;
}
#endif
