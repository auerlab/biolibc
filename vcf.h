#ifndef _BIOLIBC_VCF_H_
#define _BIOLIBC_VCF_H_

#ifndef _BIOLIBC_SAM_H_
#include "sam.h"
#endif

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

#ifndef _bool_true_false_are_defined
#include <stdbool.h>
#endif

/*
 *  vcfio is meant to provide a very simple and fast method for processing
 *  VCF streams one call at a time.  As such, there should generally be only
 *  one or a few bl_vcf_t structures substantiated at a given moment, and
 *  we can afford to be generous with the max sizes.
 *  If you're writing programs that inhale many VCF calls into memory, vcfio
 *  is not for you.
 */

// FIXME: Arbitrary guess
// Only used for temp variables.  Replace with read_field_malloc().
#define BL_VCF_SAMPLE_ID_MAX_CHARS    4096

// Hack:
// Use different sizes for each so xt_dsv_read_field() buffer overflow errors
// will point to a specific field.  Eventually should have xt_dsv_read_field()
// return an error code rather than exit with an error message
//#define BL_VCF_REF_MAX_CHARS        32
//#define BL_VCF_ALT_MAX_CHARS        33
//#define BL_VCF_QUAL_MAX_CHARS       34
//#define BL_VCF_FILTER_MAX_CHARS     64

// We actually saw INFO fields over 512k in some dbGap BCFs
typedef struct
{
    char        *chrom,
		*id,
		*ref,
		*alt,
		*qual,
		*filter,
		*info,
		*format,
		*single_sample,     // Simpler than using multi_samples
		**multi_samples;
    int64_t     pos;
    size_t      chrom_array_size,
		chrom_len,
		id_array_size,
		id_len,
		ref_array_size,
		ref_len,
		alt_array_size,
		alt_len,
		qual_array_size,
		qual_len,
		filter_array_size,
		filter_len,
		info_array_size,
		info_len,
		format_array_size,
		format_len,
		single_sample_array_size,
		single_sample_len,
		multi_sample_pointer_array_size,
		multi_sample_count,
		*multi_sample_array_sizes,
		*multi_sample_lens;
    unsigned    ref_count,
		alt_count,
		other_count;
    
    // Apps can buffer phred scores from reads to collect stats
    unsigned char   *phreds;
    size_t          phred_count;
    size_t          phred_buff_size;
}   bl_vcf_t;

typedef unsigned int vcf_field_mask_t;

#define BL_VCF_FIELD_ALL        0xfff
#define BL_VCF_FIELD_CHROM      0x001
#define BL_VCF_FIELD_POS        0x002
#define BL_VCF_FIELD_ID         0x004
#define BL_VCF_FIELD_REF        0x008
#define BL_VCF_FIELD_ALT        0x010
#define BL_VCF_FIELD_QUAL       0x020
#define BL_VCF_FIELD_FILTER     0x040
#define BL_VCF_FIELD_INFO       0x080
#define BL_VCF_FIELD_FORMAT     0x100
#define BL_VCF_FIELD_ERROR      0x000

#include "vcf-rvs.h"
#include "vcf-accessors.h"
#include "vcf-mutators.h"

/* vcf.c */
FILE *bl_vcf_skip_meta_data(FILE *vcf_stream);
FILE *bl_vcf_skip_header(FILE *vcf_stream);
void bl_vcf_get_sample_ids(FILE *vcf_stream, char *sample_ids[], size_t first_col, size_t last_col);
int bl_vcf_read_static_fields(bl_vcf_t *vcf_call, FILE *vcf_stream, vcf_field_mask_t field_mask);
int bl_vcf_read_ss_call(bl_vcf_t *vcf_call, FILE *vcf_stream, vcf_field_mask_t field_mask);
int bl_vcf_write_static_fields(bl_vcf_t *vcf_call, FILE *vcf_stream, vcf_field_mask_t field_mask);
int bl_vcf_write_ss_call(bl_vcf_t *vcf_call, FILE *vcf_stream, vcf_field_mask_t field_mask);
char **bl_vcf_sample_alloc(bl_vcf_t *vcf_call, size_t samples);
void bl_vcf_free(bl_vcf_t *vcf_call);
void bl_vcf_init(bl_vcf_t *vcf_call);
vcf_field_mask_t bl_vcf_parse_field_spec(char *spec);
_Bool bl_vcf_call_in_alignment(bl_vcf_t *vcf_call, bl_sam_t *sam_alignment);
_Bool bl_vcf_call_downstream_of_alignment(bl_vcf_t *vcf_call, bl_sam_t *alignment);
void bl_vcf_call_out_of_order(bl_vcf_t *vcf_call, char *previous_chrom, int64_t previous_pos);

#endif // _BIOLIBC_VCF_H_
