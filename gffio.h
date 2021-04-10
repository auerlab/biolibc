#ifndef __gffio_h__
#define __gffio_h__

#ifndef __biolibc_h__
#include "biolibc.h"
#endif

#define GFF_NAME_MAX_CHARS          256
#define GFF_SCORE_MAX_DIGITS        64      // Floating point
#define GFF_SOURCE_MAX_CHARS        1024    // Guess
#define GFF_FEATURE_MAX_CHARS       1024    // Guess
#define GFF_SCORE_UNAVAILABLE       -1.0

typedef unsigned int        gff_field_mask_t;

#define GFF_FIELD_ALL       0x0
#define GFF_FIELD_CHROM     0x1
#define GFF_FIELD_START_POS 0x2
#define GFF_FIELD_END_POS   0x4

#define GFF_SEQUENCE(gff_feature)   ((gff_feature)->sequence)
#define GFF_START_POS(gff_feature)  ((gff_feature)->start_pos)
#define GFF_END_POS(gff_feature)    ((gff_feature)->end_pos)

typedef struct
{
    char            sequence[BIO_CHROMOSOME_MAX_CHARS + 1];
    char            source[GFF_SOURCE_MAX_CHARS + 1];
    char            feature[GFF_FEATURE_MAX_CHARS + 1];
    char            start_pos_str[BIO_POSITION_MAX_DIGITS + 1],  // 0-based
		    end_pos_str[BIO_POSITION_MAX_DIGITS + 1];   
    uint64_t        start_pos,
		    end_pos;
    char            score_str[GFF_SCORE_MAX_DIGITS + 1];
    double          score;
    
    // char strand; '+' or '-' or '.'
    // phase 0, 1, 2 (for CDS features) or "." (for everything else)
    // char attributes[]
}   gff_feature_t;

FILE *gff_skip_header(FILE *gff_stream);
int gff_read_feature(FILE *gff_stream, gff_feature_t *gff_feature);
int gff_write_feature(FILE *gff_stream, gff_feature_t *gff_feature, gff_field_mask_t field_mask);

#endif  // __gffio_h__
