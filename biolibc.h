#ifndef _BIOLIBC_H_
#define _BIOLIBC_H_

#ifndef _STDIO_H_
#include <stdio.h>          // FILE
#endif

#ifndef _SYS_STDINT_H_
#include <stdint.h>         // int64_t
#endif

#ifndef _INTTYPES_H_
#include <inttypes.h>       // PRId64
#endif

#define BL_READ_OK              0
#define BL_READ_EOF             -1
#define BL_READ_OVERFLOW        -2
#define BL_READ_TRUNCATED       -3
#define BL_READ_GFF3_TERMINATOR  -4
#define BL_READ_EXTRA_COLS      -5
#define BL_READ_MISMATCH        -6
#define BL_READ_BAD_DATA        -7
#define BL_READ_UNKNOWN_FORMAT  -8

#define BL_WRITE_OK             0
#define BL_WRITE_FAILURE        -1

#define BL_CHROM_MAX_CHARS      256
#define BL_POSITION_MAX_DIGITS  32

#define BL_CMD_MAX              4096    // Arbitrary

#endif  // _BIOLIBC_H_
