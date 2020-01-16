// FIXME: Are there limits defined by the VCF format?
#define VCF_ID_MAX_LEN          256
#define VCF_CHROMOSOME_NAME_MAX 256
#define VCF_POSITION_MAX_DIGITS 32
#define VCF_REF_NAME_MAX        32
#define VCF_ALT_NAME_MAX        32
#define VCF_GENOTYPE_NAME_MAX   256
#define VCF_FORMAT_MAX          4096
#define VCF_DUP_CALL_MAX        10  // FIXME: What's the real maximum?

typedef struct
{
    char    chromosome[VCF_CHROMOSOME_NAME_MAX + 1],
	    pos_str[VCF_POSITION_MAX_DIGITS + 1],
	    ref[VCF_REF_NAME_MAX + 1],
	    alt[VCF_ALT_NAME_MAX + 1],
	    format[VCF_FORMAT_MAX + 1],
	    genotype[VCF_GENOTYPE_NAME_MAX + 1];    // Only used for SS calls
    size_t  pos;
    int     ref_count,
	    alt_count,
	    other_count;
}   vcf_call_t;

typedef struct
{
    size_t      count;
    vcf_call_t  call[VCF_DUP_CALL_MAX];
}   vcf_duplicate_call_t;

// CentOS 7 gcc does not support restrict
#ifdef __linux__
#define _RESTRICT
#else
#define _RESTRICT   restrict
#endif

#ifdef __linux__
size_t strlcpy(char * _RESTRICT dest, const char * _RESTRICT src, size_t len);
#endif

#include "vcfio-protos.h"
