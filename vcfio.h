// FIXME: Are there limits defined by the VCF format?
#define VCF_ID_MAX_CHARS            256
#define VCF_DUP_CALL_MAX            10  // FIXME: What's the real maximum?

/*
 *  vcfio is meant to provide a very simple and fast method for processing
 *  VCF streams one call at a time.  As such, there should generally be only
 *  one or a few vcf_call_t structures substantiated at a given moment, and
 *  we can afford to be generous with the max sizes.
 *  If you're writing programs that inhale many VCF calls into memory, vcfio
 *  is not for you.
 */

#define VCF_CHROMOSOME_MAX_CHARS    256
#define VCF_POSITION_MAX_CHARS      32
#define VCF_REF_MAX_CHARS           32
#define VCF_ALT_MAX_CHARS           32
#define VCF_QUALITY_MAX_CHARS       32
#define VCF_FILTER_MAX_CHARS        64
#define VCF_INFO_MAX_CHARS          65536
#define VCF_FORMAT_MAX_CHARS        4096
#define VCF_SAMPLE_MAX_CHARS        1024

#define VCF_CALL_INIT \
	{ \
	    "", "", "", "", "", "", "", "", "", NULL, 0, 0, 0, 0 \
	}

// Access macros.  Separate interface from implementation, so client programs
// don't reference structure members explicitly.
#define VCF_GET_CHROMOSOME(call)    ((call).chromosome)
#define VCF_GET_POS_STR(call)       ((call).pos_str)
#define VCF_GET_ID(call)            ((call).id)
#define VCF_GET_REF(call)           ((call).ref)
#define VCF_GET_ALT(call)           ((call).alt)
#define VCF_GET_QUALITY_STR(call)   ((call).quality_str)
#define VCF_GET_FILTER_STR(call)    ((call).filter_str)
#define VCF_GET_INFO(call)          ((call).info)
#define VCF_GET_FORMAT(call)        ((call).format)
#define VCF_GET_SAMPLE(call, index) ((call).samples)[index])

typedef struct
{
    char    chromosome[VCF_CHROMOSOME_MAX_CHARS + 1],
	    pos_str[VCF_POSITION_MAX_CHARS + 1],
	    id[VCF_ID_MAX_CHARS + 1],
	    ref[VCF_REF_MAX_CHARS + 1],
	    alt[VCF_ALT_MAX_CHARS + 1],
	    quality[VCF_QUALITY_MAX_CHARS + 1],
	    filter[VCF_FILTER_MAX_CHARS + 1],
	    info[VCF_INFO_MAX_CHARS + 1],
	    format[VCF_FORMAT_MAX_CHARS + 1],
	    **samples;
    size_t  pos,
	    info_len;
    int     ref_count,
	    alt_count,
	    other_count;
}   vcf_call_t;

typedef struct
{
    size_t      count;
    vcf_call_t  call[VCF_DUP_CALL_MAX];
}   vcf_duplicate_call_t;

// CentOS 7 gcc does not support restrict, which helps the optimizer produce
// faster code.  Keep _RESTRICT def separate from strlcpy() prototype in case
// other platforms are missing one but not the other.
#ifdef __linux__
#define _RESTRICT
#else
#define _RESTRICT   restrict
#endif

#ifdef __linux__
size_t strlcpy(char * _RESTRICT dest, const char * _RESTRICT src, size_t len);
#endif

#include "vcfio-protos.h"
