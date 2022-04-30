#ifndef _BIOLIBC_BED_H_
#define _BIOLIBC_BED_H_

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

#define BL_BED_NAME_MAX_CHARS          256
#define BL_BED_SCORE_MAX_DIGITS        4   // 0 to 1000
#define BL_BED_STRAND_MAX_CHARS        2
#define BL_BED_ITEM_RGB_MAX_CHARS      11  // 255,255,255
#define BL_BED_BLOCK_COUNT_MAX_DIGITS  5
#define BL_BED_BLOCK_SIZE_MAX_DIGITS   20  // 2^64
#define BL_BED_BLOCK_START_MAX_DIGITS  20  // 2^64

#define BL_BED_INIT \
	{ "", 0, 0, "", 0, '.', 0, 0, "", 0, NULL, NULL, 0 }

typedef struct
{
    char            chrom[BL_CHROM_MAX_CHARS + 1];
    /*
     *      12345
     *      ACCGT
     *      01234
     *
     *      chr1 0 5
     */
    int64_t         chrom_start,
		    chrom_end;
    char            name[BL_BED_NAME_MAX_CHARS + 1];
    unsigned short  score;      // aggs:0:1000
    char            strand;
    int64_t         thick_start,
		    thick_end;
    // FIXME: Store RGB in a more compact format
    char            item_rgb[BL_BED_ITEM_RGB_MAX_CHARS+1];
    unsigned short  block_count;
    int64_t         *block_sizes;
    int64_t         *block_starts;

    // Not part of BED spec
    unsigned short  fields;     // aggs:3:9
}   bl_bed_t;

typedef unsigned int            bed_field_mask_t;

/*
 *  Chromosome, start, and end are required so no mask bits are defined
 */
#define BL_BED_FIELD_NAME      0x01
#define BL_BED_FIELD_SCORE     0x02
#define BL_BED_FIELD_STRAND    0x04
#define BL_BED_FIELD_THICK     0X08
#define BL_BED_FIELD_RGB       0x10
#define BL_BED_FIELD_BLOCK     0x20
#define BL_BED_FIELD_ALL       0xff

// After bl_bed_t def for prototypes
#ifndef _BIOLIBC_GFF_H_
#include "gff.h"
#endif

#ifndef _bl_overlap_h
#include "overlap.h"
#endif

#include "bed-rvs.h"
#include "bed-accessors.h"
#include "bed-mutators.h"

/* bed.c */
FILE *bl_bed_skip_header(FILE *bed_stream);
int bl_bed_read(bl_bed_t *bed_feature, FILE *bed_stream, bed_field_mask_t field_mask);
int bl_bed_write(bl_bed_t *bed_feature, FILE *bed_stream, bed_field_mask_t field_mask);
void bl_bed_check_order(bl_bed_t *bed_feature, char last_chrom[], int64_t last_start);
int bl_bed_gff_cmp(bl_bed_t *bed_feature, bl_gff_t *gff_feature, bl_overlap_t *overlap);

#endif  // _BIOLIBC_BED_H_
