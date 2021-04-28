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

/*
 *  Chromosome, start, and end are required so no mask bits are defined
 */
#define BED_FIELD_ALL       0x0
#define BED_FIELD_NAME      0x1
#define BED_FIELD_SCORE     0x2

#define BED_CHROMOSOME(bf) ((bf)->chromosome)
#define BED_START_POS(bf)  ((bf)->start_pos)
#define BED_END_POS(bf)    ((bf)->end_pos)

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
    char            start_pos_str[BIO_POSITION_MAX_DIGITS + 1],  // 0-based
		    end_pos_str[BIO_POSITION_MAX_DIGITS + 1];   
    uint64_t        start_pos,
		    end_pos;
    char            name[BED_NAME_MAX_CHARS + 1];
    char            score_str[BED_SCORE_MAX_DIGITS + 1];
    unsigned short  score; // 0 to 1000
    char            strand;
    uint64_t        thick_start_pos,
		    thick_end_pos;
    
    // char strand; '+' or '-' or '.'
    // uint64_t     thick_start,
    //              thick_end;
    // uint32_t     rgb;
    // unsigned     block_count;
    // uint64_t     *block_sizes;   // array
    // uint64_t     *block_starts;  // array
}   bed_feature_t;

/* bedio.c */
FILE *bed_skip_header(FILE *bed_stream);
int bed_read_feature(FILE *bed_stream, bed_feature_t *bed_feature);
int bed_write_feature(FILE *bed_stream, bed_feature_t *bed_feature, bed_field_mask_t field_mask);
void bed_check_order(bed_feature_t *bed_feature, char last_chrom[], uint64_t last_pos);
int bed_gff_cmp(bed_feature_t *bed_feature, gff_feature_t *gff_feature, bio_overlap_t *overlap);
int bed_set_fields(bed_feature_t *bed_feature, unsigned fields);
int bed_set_chromosome(bed_feature_t *bed_feature, char *chromosome);
int bed_set_start_pos_str(bed_feature_t *bed_feature, char *start_pos_str);
int bed_set_start_pos(bed_feature_t *bed_feature, uint64_t start_pos);
int bed_set_end_pos_str(bed_feature_t *bed_feature, char *end_pos_str);
int bed_set_end_pos(bed_feature_t *bed_feature, uint64_t end_pos);
int bed_set_name(bed_feature_t *bed_feature, char *name);
int bed_set_score(bed_feature_t *bed_feature, unsigned score);
int bed_set_strand(bed_feature_t *bed_feature, int strand);

#endif  // __bedio_h__
