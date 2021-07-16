#ifndef _vcf_h_
#define _vcf_h_

#ifndef __xtend_h__
#include <xtend.h>
#endif

#ifndef _sam_h_
#include "sam.h"
#endif

#ifndef _biolibc_h_
#include "biolibc.h"
#endif

#ifndef _bool_true_false_are_defined
#include <stdbool.h>
#endif

// FIXME: Are there limits defined by the VCF format?
#define BL_VCF_ID_MAX_CHARS            256
// FIXME: What's the real maximum?  Maybe 3 since there are only 3 alternate
// alleles possible with standard bases?
#define BL_VCF_DUP_CALL_MAX            10

/*
 *  vcfio is meant to provide a very simple and fast method for processing
 *  VCF streams one call at a time.  As such, there should generally be only
 *  one or a few bl_vcf_t structures substantiated at a given moment, and
 *  we can afford to be generous with the max sizes.
 *  If you're writing programs that inhale many VCF calls into memory, vcfio
 *  is not for you.
 */

// Hack:
// Use different sizes for each so dsv_read_field() buffer overflow errors
// will point to a specific field.  Eventually should have dsv_read_field()
// return an error code rather than exit with an error message
#define BL_VCF_REF_MAX_CHARS           32
#define BL_VCF_ALT_MAX_CHARS           33
#define BL_VCF_QUALITY_MAX_CHARS       34
#define BL_VCF_FILTER_MAX_CHARS        64

#define BL_VCF_CALL_INIT   { "", "", "", "", "", "", "", NULL, NULL, NULL, \
			    0, 0, 0, 0, 0, 0, 0, 0, \
			    NULL, NULL, 0, 0 \
			}
typedef struct
{
    char    chromosome[BL_CHROMOSOME_MAX_CHARS + 1],
	    id[BL_VCF_ID_MAX_CHARS + 1],
	    ref[BL_VCF_REF_MAX_CHARS + 1],
	    alt[BL_VCF_ALT_MAX_CHARS + 1],
	    quality[BL_VCF_QUALITY_MAX_CHARS + 1],
	    filter[BL_VCF_FILTER_MAX_CHARS + 1],
	    // We actually saw INFO fields over 512k in some dbGap BCFs
	    *info,
	    *format,
	    *single_sample; // Avoid using multi_samples
    size_t  pos,
	    info_len,
	    info_max,
	    format_max,
	    sample_max;
    int     ref_count,
	    alt_count,
	    other_count;
    // Use vcf_sample_alloc() to initialize this for multi-sample VCFs
    char    **multi_samples;
    // Apps can buffer phred scores from reads to collect stats
    unsigned char   *phreds;
    size_t  phred_count;
    size_t  phred_buff_size;
}   bl_vcf_t;

#define BL_VCF_CHROMOSOME(ptr)  ((ptr)->chromosome)
#define BL_VCF_ID(ptr)  ((ptr)->id)
#define BL_VCF_REF(ptr) ((ptr)->ref)
#define BL_VCF_ALT(ptr) ((ptr)->alt)
#define BL_VCF_QUALITY(ptr) ((ptr)->quality)
#define BL_VCF_FILTER(ptr)  ((ptr)->filter)
#define BL_VCF_POS(ptr) ((ptr)->pos)
#define BL_VCF_INFO_LEN(ptr)    ((ptr)->info_len)
#define BL_VCF_INFO_MAX(ptr)    ((ptr)->info_max)
#define BL_VCF_FORMAT_MAX(ptr)  ((ptr)->format_max)
#define BL_VCF_SAMPLE_MAX(ptr)  ((ptr)->sample_max)
#define BL_VCF_REF_COUNT(ptr)   ((ptr)->ref_count)
#define BL_VCF_ALT_COUNT(ptr)   ((ptr)->alt_count)
#define BL_VCF_OTHER_COUNT(ptr) ((ptr)->other_count)
#define BL_VCF_MULTI_SAMPLES(ptr,c)   ((ptr)->multi_samples[c])
#define BL_VCF_PHREDS(ptr)  ((ptr)->phreds)
#define BL_VCF_PHRED_COUNT(ptr) ((ptr)->phred_count)
#define BL_VCF_PHRED_BUFF_SIZE(ptr) ((ptr)->phred_buff_size)

#define BL_VCF_SET_CHROMOSOME(ptr,chromosome)   strlcpy(ptr->chromosome,chromosome,BL_CHROMOSOME_MAX_CHARS+1)
#define BL_VCF_SET_ID(ptr,id)                   strlcpy(ptr->id,id,BL_VCF_ID_MAX_CHARS+1)
#define BL_VCF_SET_REF(ptr,ref)                 strlcpy(ptr->ref,ref,BL_VCF_REF_MAX_CHARS+1)
#define BL_VCF_SET_ALT(ptr,alt)                 strlcpy(ptr->alt,alt,BL_VCF_ALT_MAX_CHARS+1)
#define BL_VCF_SET_QUALITY(ptr,quality)         strlcpy(ptr->quality,quality,BL_VCF_QUALITY_MAX_CHARS+1)
#define BL_VCF_SET_FILTER(ptr,filter)           strlcpy(ptr->filter,filter,BL_VCF_FILTER_MAX_CHARS+1)
#define BL_VCF_SET_POS(ptr,pos)                 ((ptr)->pos = (pos))
#define BL_VCF_SET_INFO_LEN(ptr,info_len)       ((ptr)->info_len = (info_len))
#define BL_VCF_SET_INFO_MAX(ptr,info_max)       ((ptr)->info_max = (info_max))
#define BL_VCF_SET_FORMAT_MAX(ptr,format_max)   ((ptr)->format_max = (format_max))
#define BL_VCF_SET_SAMPLE_MAX(ptr,sample_max)   ((ptr)->sample_max = (sample_max))
#define BL_VCF_SET_REF_COUNT(ptr,ref_count)     ((ptr)->ref_count = (ref_count))
#define BL_VCF_SET_ALT_COUNT(ptr,alt_count)     ((ptr)->alt_count = (alt_count))
#define BL_VCF_SET_OTHER_COUNT(ptr,other_count) ((ptr)->other_count = (other_count))
#define BL_VCF_SET_MULTI_SAMPLES(ptr, c, samples)   ((ptr)->multi_samples[c] = (sample))
#define BL_VCF_SET_PHREDS(ptr,phreds)           ((ptr)->phreds = (phreds))
#define BL_VCF_SET_PHRED_COUNT(ptr,phred_count) ((ptr)->phred_count = (phred_count))
#define BL_VCF_SET_PHRED_BUFF_SIZE(ptr,phred_buff_size) ((ptr)->phred_buff_size = (phred_buff_size))

typedef unsigned int        vcf_field_mask_t;

#define BL_VCF_FIELD_ALL       0xfff
#define BL_VCF_FIELD_CHROM     0x001
#define BL_VCF_FIELD_POS       0x002
#define BL_VCF_FIELD_ID        0x004
#define BL_VCF_FIELD_REF       0x008
#define BL_VCF_FIELD_ALT       0x010
#define BL_VCF_FIELD_QUAL      0x020
#define BL_VCF_FIELD_FILTER    0x040
#define BL_VCF_FIELD_INFO      0x080
#define BL_VCF_FIELD_FORMAT    0x100
#define BL_VCF_FIELD_ERROR     0x000

// Future expansion: Copy all or part of header
typedef unsigned int        vcf_header_t;

#define BL_VCF_HEADER_NONE     0x0
#define BL_VCF_HEADER_FORMAT   0x1
#define BL_VCF_HEADER_ALL      0x1

/* vcf.c */
FILE *vcf_skip_header(FILE *vcf_stream);
void vcf_get_sample_ids(FILE *vcf_stream, char *sample_ids[], size_t first_col, size_t last_col);
int vcf_read_static_fields(FILE *vcf_stream, bl_vcf_t *vcf_call, vcf_field_mask_t field_mask);
int vcf_read_ss_call(FILE *vcf_stream, bl_vcf_t *vcf_call, vcf_field_mask_t field_mask);
int vcf_write_static_fields(FILE *vcf_stream, bl_vcf_t *vcf_call, vcf_field_mask_t field_mask);
int vcf_write_ss_call(FILE *vcf_stream, bl_vcf_t *vcf_call, vcf_field_mask_t field_mask);
char **vcf_sample_alloc(bl_vcf_t *vcf_call, size_t samples);
void vcf_call_free(bl_vcf_t *vcf_call);
void vcf_call_init(bl_vcf_t *vcf_call,
		   size_t info_max, size_t format_max, size_t sample_max);
vcf_field_mask_t vcf_parse_field_spec(char *spec);
bool vcf_call_in_alignment(bl_vcf_t *vcf_call, bl_sam_t *sam_alignment);
bool vcf_call_downstream_of_alignment(bl_vcf_t *vcf_call, bl_sam_t *alignment);
void vcf_out_of_order(bl_vcf_t *vcf_call,
			 char *previous_chromosome, size_t previous_pos);

#endif // _vcf_h_
