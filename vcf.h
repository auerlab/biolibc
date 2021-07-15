#ifndef _vcf_h_
#define _vcf_h_

#ifndef _dsv_h_
#include "dsv.h"
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
#define VCF_ID_MAX_CHARS            256
// FIXME: What's the real maximum?  Maybe 3 since there are only 3 alternate
// alleles possible with standard bases?
#define VCF_DUP_CALL_MAX            10

/*
 *  vcfio is meant to provide a very simple and fast method for processing
 *  VCF streams one call at a time.  As such, there should generally be only
 *  one or a few vcf_call_t structures substantiated at a given moment, and
 *  we can afford to be generous with the max sizes.
 *  If you're writing programs that inhale many VCF calls into memory, vcfio
 *  is not for you.
 */

// Hack:
// Use different sizes for each so dsv_read_field() buffer overflow errors
// will point to a specific field.  Eventually should have dsv_read_field()
// return an error code rather than exit with an error message
#define VCF_REF_MAX_CHARS           32
#define VCF_ALT_MAX_CHARS           33
#define VCF_QUALITY_MAX_CHARS       34
#define VCF_FILTER_MAX_CHARS        64
// Yes, we actually saw INFO fields over 512k in some dbGap BCFs
//#define VCF_INFO_MAX_CHARS          1048576
//#define VCF_FORMAT_MAX_CHARS        4096
//#define VCF_SAMPLE_MAX_CHARS        2048

// Access macros.  Separate interface from implementation, so client programs
// don't reference structure members explicitly.
#define VCF_CHROMOSOME(vcf_call)    ((vcf_call)->chromosome)
#define VCF_POS(vcf_call)           ((vcf_call)->pos)
#define VCF_POS_STR(vcf_call)       ((vcf_call)->pos_str)
#define VCF_ID(vcf_call)            ((vcf_call)->id)
#define VCF_REF(vcf_call)           ((vcf_call)->ref)
#define VCF_ALT(vcf_call)           ((vcf_call)->alt)
#define VCF_QUAL(vcf_call)          ((vcf_call)->quality)
#define VCF_FILTER(vcf_call)        ((vcf_call)->filter)
#define VCF_INFO(vcf_call)          ((vcf_call)->info)
#define VCF_FORMAT(vcf_call)        ((vcf_call)->format)
#define VCF_INFO_LEN(vcf_call)      ((vcf_call)->info_len)
#define VCF_REF_COUNT(vcf_call)     ((vcf_call)->ref_count)
#define VCF_ALT_COUNT(vcf_call)     ((vcf_call)->alt_count)
#define VCF_OTHER_COUNT(vcf_call)   ((vcf_call)->other_count)
#define VCF_SINGLE_SAMPLE(vcf_call) ((vcf_call)->single_sample)
#define VCF_MULTI_SAMPLE(vcf_call,i) ((vcf_call)->multi_samples, (i))
#define VCF_PHREDS(vcf_call)        ((vcf_call)->phreds)
#define VCF_PHRED_VAL(vcf_call,i)   ((vcf_call)->phreds[i])
#define VCF_PHRED_COUNT(vcf_call)   ((vcf_call)->phred_count)

#define VCF_PHRED_BUFF_SIZE         256

#define VCF_CALL_INIT   { "", "", "", "", "", "", "", NULL, NULL, NULL, \
			    0, 0, 0, 0, 0, 0, 0, 0, \
			    NULL, NULL, 0, 0 \
			}
typedef struct
{
    char    chromosome[BIO_CHROMOSOME_MAX_CHARS + 1],
	    pos_str[BIO_POSITION_MAX_DIGITS + 1],
	    id[VCF_ID_MAX_CHARS + 1],
	    ref[VCF_REF_MAX_CHARS + 1],
	    alt[VCF_ALT_MAX_CHARS + 1],
	    quality[VCF_QUALITY_MAX_CHARS + 1],
	    filter[VCF_FILTER_MAX_CHARS + 1],
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
}   vcf_call_t;

typedef unsigned int        vcf_field_mask_t;

#define VCF_FIELD_ALL       0xfff
#define VCF_FIELD_CHROM     0x001
#define VCF_FIELD_POS       0x002
#define VCF_FIELD_ID        0x004
#define VCF_FIELD_REF       0x008
#define VCF_FIELD_ALT       0x010
#define VCF_FIELD_QUAL      0x020
#define VCF_FIELD_FILTER    0x040
#define VCF_FIELD_INFO      0x080
#define VCF_FIELD_FORMAT    0x100
#define VCF_FIELD_ERROR     0x000

// Future expansion: Copy all or part of header
typedef unsigned int        vcf_header_t;

#define VCF_HEADER_NONE     0x0
#define VCF_HEADER_FORMAT   0x1
#define VCF_HEADER_ALL      0x1

/* vcf.c */
FILE *vcf_skip_header(FILE *vcf_stream);
void vcf_get_sample_ids(FILE *vcf_stream, char *sample_ids[], size_t first_col, size_t last_col);
int vcf_read_static_fields(FILE *vcf_stream, vcf_call_t *vcf_call, vcf_field_mask_t field_mask);
int vcf_read_ss_call(FILE *vcf_stream, vcf_call_t *vcf_call, vcf_field_mask_t field_mask);
int vcf_write_static_fields(FILE *vcf_stream, vcf_call_t *vcf_call, vcf_field_mask_t field_mask);
int vcf_write_ss_call(FILE *vcf_stream, vcf_call_t *vcf_call, vcf_field_mask_t field_mask);
char **vcf_sample_alloc(vcf_call_t *vcf_call, size_t samples);
void vcf_call_free(vcf_call_t *vcf_call);
void vcf_call_init(vcf_call_t *vcf_call,
		   size_t info_max, size_t format_max, size_t sample_max);
vcf_field_mask_t vcf_parse_field_spec(char *spec);
bool vcf_call_in_alignment(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment);
bool vcf_call_downstream_of_alignment(vcf_call_t *vcf_call, sam_alignment_t *alignment);
void vcf_out_of_order(vcf_call_t *vcf_call,
			 char *previous_chromosome, size_t previous_pos);

#endif // _vcf_h_
