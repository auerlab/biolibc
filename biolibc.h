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

#define BL_READ_OK              0
#define BL_READ_EOF             -1
#define BL_READ_OVERFLOW        -2
#define BL_READ_TRUNCATED       -3
#define BL_READ_GFF_TERMINATOR  -4
#define BL_READ_EXTRA_COLS      -5
#define BL_READ_MISMATCH        -6
#define BL_READ_BAD_DATA        -7

#define BL_WRITE_OK             0
#define BL_WRITE_FAILURE        -1

#define BL_DATA_OK              0
#define BL_DATA_INVALID         -1      // Catch-all for non-specific error
#define BL_DATA_OUT_OF_RANGE    -2

#define BL_CHROM_MAX_CHARS      256
#define BL_POSITION_MAX_DIGITS  32

#define BL_CMD_MAX              4096    // Arbitrary

#endif  // _biolibc_h_
