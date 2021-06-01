#ifndef _biolibc_h_
#define _biolibc_h_

#ifndef _STDIO_H_
#include <stdio.h>          // FILE
#endif

#ifndef _SYS_STDINT_H_
#include <stdint.h>         // uint64_t
#endif

#ifndef _INTTYPES_H_
#include <inttypes.h>       // PRIu64
#endif

#define BIO_READ_OK                 0
#define BIO_READ_EOF                -1
#define BIO_READ_OVERFLOW           -2
#define BIO_READ_TRUNCATED          -3
#define BIO_READ_GFF_TERMINATOR     -4
#define BIO_READ_EXTRA_COLS         -5
#define BIO_READ_MISMATCH           -6

#define BIO_WRITE_OK                0
#define BIO_WRITE_FAILURE           -1

#define BIO_DATA_OK                 0
#define BIO_DATA_INVALID            -1      // Catch-all for non-specific error
#define BIO_DATA_OUT_OF_RANGE       -2

#define BIO_CHROMOSOME_MAX_CHARS    256
#define BIO_POSITION_MAX_DIGITS     32

#define BIO_CMD_MAX                 4096    // Arbitrary

// 1-based, inclusive at both ends
typedef struct
{
    uint64_t    f1_len,
		f2_len,
		ov_start,
		ov_end,
		ov_len;
}   bio_overlap_t;

/* chromosome-name-cmp.c */
int chromosome_name_cmp(const char *n1, const char *n2);
void bio_set_overlap(bio_overlap_t *overlap,
			uint64_t f1_len, uint64_t f2_len,
			uint64_t ov_start, uint64_t ov_end);
void bio_print_overlap(bio_overlap_t *overlap, char *f1_name, char *f2_name);

#endif  // _biolibc_h_
