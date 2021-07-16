#ifndef _bed_h_
#define _bed_h_

#ifndef _biolibc_h_
#include "biolibc.h"
#endif

#define BL_BED_NAME_MAX_CHARS          256
#define BL_BED_SCORE_MAX_DIGITS        4   // 0 to 1000
#define BL_BED_STRAND_MAX_CHARS        2
#define BL_BED_RGB_STR_MAX_CHARS       11  // 255,255,255
#define BL_BED_BLOCK_COUNT_MAX_DIGITS  5
#define BL_BED_BLOCK_SIZE_MAX_DIGITS   20  // 2^64
#define BL_BED_BLOCK_START_MAX_DIGITS  20  // 2^64

typedef unsigned int        bed_field_mask_t;

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

#define BL_BED_INIT \
	{ 0, "", 0, 0, "", 0, '.', 0, 0, "", 0, NULL, NULL }

typedef struct
{
    unsigned short  fields;     // 3 to 9
    char            chromosome[BL_CHROMOSOME_MAX_CHARS + 1];
    /*
     *      12345
     *      ACCGT
     *      01234
     *
     *      chr1 0 5
     */
    uint64_t        start_pos,
		    end_pos;
    char            name[BL_BED_NAME_MAX_CHARS + 1];
    unsigned short  score; // 0 to 1000
    char            strand;
    uint64_t        thick_start_pos,
		    thick_end_pos;
    // FIXME: Store RGB in a more compact format
    char            rgb_str[BL_BED_RGB_STR_MAX_CHARS+1];
    unsigned short  block_count;
    uint64_t        *block_sizes;
    uint64_t        *block_starts;
}   bl_bed_t;

#define BL_BED_FIELDS(ptr)          ((ptr)->fields)
#define BL_BED_CHROMOSOME(ptr)      ((ptr)->chromosome)
#define BL_BED_START_POS(ptr)       ((ptr)->start_pos)
#define BL_BED_END_POS(ptr)         ((ptr)->end_pos)
#define BL_BED_NAME(ptr)            ((ptr)->name)
#define BL_BED_SCORE(ptr)           ((ptr)->score)
#define BL_BED_STRAND(ptr)          ((ptr)->strand)
#define BL_BED_THICK_START_POS(ptr) ((ptr)->thick_start_pos)
#define BL_BED_THICK_END_POS(ptr)   ((ptr)->thick_end_pos)
#define BL_BED_RGB_STR(ptr)         ((ptr)->rgb_str)
#define BL_BED_BLOCK_COUNT(ptr)     ((ptr)->block_count)
#define BL_BED_BLOCK_SIZES(ptr)     ((ptr)->block_sizes)
#define BL_BED_BLOCK_STARTS(ptr)    ((ptr)->block_starts)

#define BL_BED_SET_FIELDS(ptr,v)        ((ptr)->fields = (v))
#define BL_BED_SET_CHROMOSOME(ptr,v)    strlcpy((ptr)->chromosome,v,BL_CHROMOSOME_MAX_CHARS+1)
#define BL_BED_SET_START_POS(ptr,v)     ((ptr)->start_pos = (v))
#define BL_BED_SET_END_POS(ptr,v)       ((ptr)->end_pos = (v))
#define BL_BED_SET_NAME(ptr,v)          strlcpy((ptr)->name,v,BL_BED_NAME_MAX_CHARS+1)
#define BL_BED_SET_SCORE(ptr,v)         ((ptr)->score = (v))
#define BL_BED_SET_STRAND(ptr,v)        ((ptr)->strand = (v))
#define BL_BED_SET_THICK_START_POS(ptr,v)   ((ptr)->thick_start_pos = (v))
#define BL_BED_SET_THICK_END_POS(ptr,v) ((ptr)->thick_end_pos = (v))
#define BL_BED_SET_RGB_STR(ptr,v)       strlcpy((ptr)->rgb_str,v,BL_BED_RGB_STR_MAX_CHARS+1)
#define BL_BED_SET_BLOCK_COUNT(ptr,v)   ((ptr)->block_count = (v))
#define BL_BED_SET_BLOCK_SIZES(ptr,v)   ((ptr)->block_sizes = (v))
#define BL_BED_SET_BLOCK_STARTS(ptr,v)  ((ptr)->block_starts = (v))

// After bl_bed_t for prototypes
#ifndef _gff_h_
#include "gff.h"
#endif

#ifndef _bl_overlap_h
#include "overlap.h"
#endif

/* bed.c */
FILE *bed_skip_header(FILE *bed_stream);
int bed_read_feature(FILE *bed_stream, bl_bed_t *bed_feature, bed_field_mask_t field_mask);
int bed_write_feature(FILE *bed_stream, bl_bed_t *bed_feature, bed_field_mask_t field_mask);
void bed_check_order(bl_bed_t *bed_feature, char last_chrom[], uint64_t last_start);
int bed_gff_cmp(bl_bed_t *bed_feature, bl_gff_t *gff_feature, bl_overlap_t *overlap);
#endif  // _bed_h_
