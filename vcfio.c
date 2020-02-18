#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include "tsvio.h"
#include "vcfio.h"

/***************************************************************************
 *  Description:
 *      Skip over header lines in VCF input stream.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    vcf_skip_header(const char *argv[], FILE *vcf_stream)

{
    char    start[7] = "xxxxxx";
    size_t  count;

    while ( ((count=fread(start, 6, 1, vcf_stream)) == 1) && 
	    (memcmp(start, "#CHROM", 6) != 0) )
	tsv_skip_rest_of_line(argv, vcf_stream);
    
    // puts(start);
    if ( count == 0 )
    {
	fprintf(stderr, "%s: vcf_skip_header(): No #CHROM header found.\n", argv[0]);
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

void    vcf_get_sample_ids(const char *argv[], FILE *vcf_stream,
			   char *sample_ids[],
			   size_t first_col, size_t last_col)

{
    size_t  c,
	    len;
    char    temp_sample_id[VCF_ID_MAX_CHARS + 1];
    int     end_of_field;
    
    // Skip standard header tags to get to sample IDs
    for (c = 0; c < 9; ++c)
	tsv_skip_field(argv, vcf_stream);
    
    // Skip sample IDs before first_col
    for (c = 1; c < first_col; ++c)
	tsv_skip_field(argv, vcf_stream);
    
    for (; (c <= last_col) &&
	   (end_of_field = tsv_read_field(argv, vcf_stream, temp_sample_id,
				     VCF_ID_MAX_CHARS, &len)) != EOF; ++c)
    {
	sample_ids[c - first_col] = strdup(temp_sample_id);
	// fprintf(stderr, "'%s'\n", temp_sample_id);
    }
    
    // Skip any remaining fields after last_col
    if ( end_of_field != '\n' )
	tsv_skip_rest_of_line(argv, vcf_stream);
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

int     vcf_read_static_fields(const char *argv[],
		      FILE *vcf_stream, vcf_call_t *vcf_call)

{
    char    *end;
    size_t  len;
    
    vcf_call->ref_count = vcf_call->alt_count = vcf_call->other_count = 0;
    
    // Chromosome
    if ( tsv_read_field(argv, vcf_stream, vcf_call->chromosome,
			VCF_CHROMOSOME_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading CHROM: %s.\n",
		argv[0], vcf_call->chromosome);
	fputs("This is normal.\n", stderr);
	return 0;
    }
    
    // Call position
    if ( tsv_read_field(argv, vcf_stream, vcf_call->pos_str,
			VCF_POSITION_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading POS: %s.\n",
		argv[0], vcf_call->pos_str);
	return 0;
    }
    else
    {
	vcf_call->pos = strtoul(vcf_call->pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "%s: vcf_read_static_fields(): Invalid call position: %s\n",
		    argv[0], vcf_call->pos_str);
	    return 0;
	}
    }
    
    // ID
    if ( tsv_read_field(argv, vcf_stream, vcf_call->id,
			VCF_ID_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading ID.\n",
		argv[0]);
	return 0;
    }
    
    // Ref
    if ( tsv_read_field(argv, vcf_stream, vcf_call->ref,
			VCF_REF_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading REF.\n",
		argv[0]);
	return 0;
    }
    
    // Alt
    if ( tsv_read_field(argv, vcf_stream, vcf_call->alt,
		   VCF_ALT_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading ALT.\n",
		argv[0]);
	return 0;
    }

    // Qual
    if ( tsv_read_field(argv, vcf_stream, vcf_call->quality,
		   VCF_QUALITY_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading QUAL.\n",
		argv[0]);
	return 0;
    }
    
    // Filter
    if ( tsv_read_field(argv, vcf_stream, vcf_call->filter,
		   VCF_FILTER_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading FILTER.\n",
		argv[0]);
	return 0;
    }
    
    // Info
    if ( tsv_read_field(argv, vcf_stream, vcf_call->info,
		   VCF_INFO_MAX_CHARS, &vcf_call->info_len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading INFO.\n",
		argv[0]);
	return 0;
    }
    
    // Format
    if ( tsv_read_field(argv, vcf_stream, vcf_call->format,
		   VCF_FORMAT_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr,
		"%s: vcf_read_static_fields(): Hit EOF reading FORMAT.\n",
		argv[0]);
	return 0;
    }

#if 0
    fprintf(stderr, "%s %s %s %s %s %s\n",
	vcf_call->chromosome,
	vcf_call->pos_str,
	vcf_call->ref,
	vcf_call->alt,
	vcf_call->format);
#endif
    return 1;
}


/***************************************************************************
 *  Description:
 *      Read a single-sample VCF call.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-11  Jason Bacon Begin
 ***************************************************************************/

int     vcf_read_ss_call(const char *argv[],
		      FILE *vcf_stream, vcf_call_t *vcf_call)

{
    size_t  len;
    
    if ( vcf_read_static_fields(argv, vcf_stream, vcf_call) )
    {
	if ( vcf_sample_alloc(vcf_call, 1) != NULL )
	{
	    if ( tsv_read_field(argv, vcf_stream, vcf_call->samples[0],
			    VCF_SAMPLE_MAX_CHARS, &len) != EOF )
		return 1;
	    else
	    {
		fprintf(stderr,
			"%s: vcf_read_ss_call(): Hit EOF reading sample.\n",
			argv[0]);
		return 0;
	    }
	}
	else
	{
	    fprintf(stderr, "%s: vcf_read_ss_call(): malloc() failed.\n",
		    argv[0]);
	    return 0;
	}
    }
    else
	return 0;
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

int     vcf_write_static_fields(const char *argv[],
		      FILE *vcf_stream, vcf_call_t *vcf_call)

{
    return 0;
}


/***************************************************************************
 *  Description:
 *      Write a single-sample VCF call to vcf_stream.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     vcf_write_ss_call(const char *argv[],
		      FILE *vcf_stream, vcf_call_t *vcf_call)

{
    return 0;
}


/***************************************************************************
 *  Description:
 *      Read in all consecutive calls with the same position.
 *
 *  Returns:
 *      Array of calls and number of calls.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-11  Jason Bacon Begin
 ***************************************************************************/

size_t  vcf_read_duplicate_calls(const char *argv[], FILE *vcf_stream,
				 vcf_duplicate_call_t *vcf_duplicate_calls)

{
    // Cache the next VCF call after the last one returned
    static vcf_call_t   vcf_call = VCF_CALL_INIT;
    static size_t       buffered_calls = 0; // 0 or 1
    size_t              c;
    
    // Prime cache with the first call
    if ( buffered_calls == 0 )
	if ( (buffered_calls = vcf_read_ss_call(argv, vcf_stream, &vcf_call)) == 0 )
	    return 0;
    
    /*
     *  Read all VCF calls with the same position.  The first one with a
     *  different position is not added to vcf_duplicate_calls and is left in
     *  the static vcf_call for the next invocation of this function.
     */
    c = 0;
    do
    {
	vcf_duplicate_calls->call[c++] = vcf_call;
	buffered_calls = vcf_read_ss_call(argv, vcf_stream, &vcf_call);
    }   while ( (buffered_calls == 1) &&
		(vcf_call.pos == vcf_duplicate_calls->call[c-1].pos) &&
		(strcmp(vcf_call.chromosome,
			vcf_duplicate_calls->call[c-1].chromosome) == 0) );

    // Return the number of calls with the same position
    return vcf_duplicate_calls->count = c;
}


char    **vcf_sample_alloc(vcf_call_t *vcf_call, size_t samples)

{
    size_t  c;
    
    if ( (vcf_call->samples =
	 (char **)malloc(samples * sizeof(char *))) != NULL )
    {
	for (c = 0; c < samples; ++c)
	{
	    if ( (vcf_call->samples[c] =
		 (char *)malloc(VCF_SAMPLE_MAX_CHARS + 1)) == NULL )
		return NULL;
	}
    }
    return vcf_call->samples;
}


#ifdef __linux__
size_t  strlcpy(char * _RESTRICT dest, const char * _RESTRICT src, size_t len)

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
