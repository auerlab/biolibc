#ifndef __biolibc_h__
#define __biolibc_h__

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

#define BIO_DATA_OK                 0
#define BIO_OUT_OF_RANGE            -1
#define BIO_INVALID_DATA            -2

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

void    bio_set_overlap(bio_overlap_t *overlap,
			uint64_t f1_len, uint64_t f2_len,
			uint64_t ov_start, uint64_t ov_end);
void    bio_print_overlap(bio_overlap_t *overlap, char *f1_name, char *f2_name);
FILE    *bio_fopen(char *filename, char *mode);
int     bio_fclose(FILE *stream);

#endif  // __biolibc_h__
