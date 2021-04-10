#ifndef __bedio_h__
#define __bedio_h__

#ifndef __biolibc_h__
#include "biolibc.h"
#endif

#ifndef __gffio_h__
#include "gffio.h"
#endif

#define BED_NAME_MAX_CHARS          256
#define BED_SCORE_MAX_DIGITS        4   // 0 to 1000

typedef unsigned int        bed_field_mask_t;

#define BED_FIELD_ALL       0x0
#define BED_FIELD_CHROM     0x1
#define BED_FIELD_START_POS 0x2
#define BED_FIELD_END_POS   0x4

#define BED_CHROMOSOME(bed_feature) ((bed_feature)->chromosome)
#define BED_START_POS(bed_feature)  ((bed_feature)->start_pos)
#define BED_END_POS(bed_feature)    ((bed_feature)->end_pos)

typedef struct
{
    unsigned short  fields;
    char            chromosome[BIO_CHROMOSOME_MAX_CHARS + 1];
    /*
     *      12345
     *      ACCGT
     *      01234
     *
     *      chr1 0 5
     */
    char            start_pos_str[BIO_POSITION_MAX_DIGITS + 1],  // 0-based
		    end_pos_str[BIO_POSITION_MAX_DIGITS + 1];   
    uint64_t        start_pos,
		    end_pos;
    char            name[BED_NAME_MAX_CHARS + 1];
    char            score_str[BED_SCORE_MAX_DIGITS + 1];
    unsigned short  score; // 0 to 1000
    
    // char strand; '+' or '-'
    // uint64_t     thick_start,
    //              thick_end;
    // uint32_t     rgb;
    // unsigned     block_count;
    // uint64_t     *block_sizes;   // array
    // uint64_t     *block_starts;  // array
}   bed_feature_t;

FILE *bed_skip_header(FILE *bed_stream);
int bed_read_feature(FILE *bed_stream, bed_feature_t *bed_feature);
int bed_write_feature(FILE *bed_stream, bed_feature_t *bed_feature, bed_field_mask_t field_mask);
void    bed_check_order(bed_feature_t *bed_feature, char last_chrom[],
			uint64_t last_pos);
int     bed_gff_cmp(bed_feature_t *bed_feature, gff_feature_t *gff_feature,
		    bio_overlap_t *overlap);

#endif  // __bedio_h__
