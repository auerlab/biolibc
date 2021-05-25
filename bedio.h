#ifndef __bedio_h__
#define __bedio_h__

#ifndef __biolibc_h__
#include "biolibc.h"
#endif

#define BED_NAME_MAX_CHARS          256
#define BED_SCORE_MAX_DIGITS        4   // 0 to 1000
#define BED_STRAND_MAX_CHARS        2
#define BED_RGB_STR_MAX_CHARS       11  // 255,255,255
#define BED_BLOCK_COUNT_MAX_DIGITS  5
#define BED_BLOCK_SIZE_MAX_DIGITS   20  // 2^64
#define BED_BLOCK_START_MAX_DIGITS  20  // 2^64

typedef unsigned int        bed_field_mask_t;

/*
 *  Chromosome, start, and end are required so no mask bits are defined
 */
#define BED_FIELD_ALL       0x00
#define BED_FIELD_NAME      0x01
#define BED_FIELD_SCORE     0x02
#define BED_FIELD_STRAND    0x04
#define BED_FIELD_THICK     0X08
#define BED_FIELD_RGB       0x10
#define BED_FIELD_BLOCK     0x20

// FIXME: Write man pages for accessors
#define BED_CHROMOSOME(bf)      ((bf)->chromosome)
#define BED_START_POS(bf)       ((bf)->start_pos)
#define BED_END_POS(bf)         ((bf)->end_pos)
#define BED_SCORE(bf)           ((bf)->score)
#define BED_STRAND(bf)          ((bf)->strand)
#define BED_THICK_START_POS(bf) ((bf)->thick_start_pos)
#define BED_THICK_END_POS(bf)   ((bf)->thick_end_pos)
#define BED_RGB_STR(bf)         ((bf)->rgb_str)
#define BED_BLOCK_COUNT(bf)     ((bf)->block_count)
#define BED_BLOCK_SIZES(bf,i)   ((bf)->block_sizes[i])
#define BED_BLOCK_STARTS(bf,i)  ((bf)->block_starts[i])

#define BED_INIT \
	{ 0, "", 0, 0, "", 0, '.', 0, 0, "", 0, NULL, NULL }

typedef struct
{
    unsigned short  fields;     // 3 to 9
    char            chromosome[BIO_CHROMOSOME_MAX_CHARS + 1];
    /*
     *      12345
     *      ACCGT
     *      01234
     *
     *      chr1 0 5
     */
    uint64_t        start_pos,
		    end_pos;
    char            name[BED_NAME_MAX_CHARS + 1];
    unsigned short  score; // 0 to 1000
    char            strand;
    uint64_t        thick_start_pos,
		    thick_end_pos;
    // FIXME: Store RGB in a more compact format
    char            rgb_str[BED_RGB_STR_MAX_CHARS+1];
    unsigned short  block_count;
    uint64_t        *block_sizes;
    uint64_t        *block_starts;
}   bed_feature_t;

// After bed_feature_t for prototypes
#ifndef __gffio_h__
#include "gffio.h"
#endif

/* bedio.c */
FILE *bed_skip_header(FILE *bed_stream);
int bed_read_feature(FILE *bed_stream, bed_feature_t *bed_feature, bed_field_mask_t field_mask);
int bed_write_feature(FILE *bed_stream, bed_feature_t *bed_feature, bed_field_mask_t field_mask);
void bed_check_order(bed_feature_t *bed_feature, char last_chrom[], uint64_t last_start);
int bed_gff_cmp(bed_feature_t *bed_feature, gff_feature_t *gff_feature, bio_overlap_t *overlap);
int bed_set_fields(bed_feature_t *bed_feature, unsigned fields);
int bed_set_chromosome(bed_feature_t *bed_feature, char *chromosome);
int bed_set_start_pos(bed_feature_t *bed_feature, uint64_t start_pos);
int bed_set_end_pos(bed_feature_t *bed_feature, uint64_t end_pos);
int bed_set_name(bed_feature_t *bed_feature, char *name);
int bed_set_score(bed_feature_t *feature, unsigned score);
int bed_set_strand(bed_feature_t *feature, int strand);
int bed_set_thick_start_pos(bed_feature_t *bed_feature, uint64_t thick_start_pos);
int bed_set_thick_end_pos(bed_feature_t *bed_feature, uint64_t thick_end_pos);
int bed_set_rgb_str(bed_feature_t *bed_feature, char *rgb_str);
#endif  // __bedio_h__
