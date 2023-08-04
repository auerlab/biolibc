#ifndef _BIOLIBC_GFF3_H_
#define _BIOLIBC_GFF3_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

#define BL_GFF3_SOURCE_MAX_CHARS     1024     // Guess
#define BL_GFF3_TYPE_MAX_CHARS       256      // Guess
#define BL_GFF3_SCORE_MAX_DIGITS     64       // Floating point
#define BL_GFF3_STRAND_MAX_CHARS     2
#define BL_GFF3_LINE_MAX_CHARS       32768
#define BL_GFF3_PHASE_MAX_DIGITS     2

#define BL_GFF3_SCORE_UNAVAILABLE    -1.0
#define BL_GFF3_PHASE_UNAVAILABLE    '.'

typedef struct
{
    char            seqid[BL_CHROM_MAX_CHARS + 1];
    char            source[BL_GFF3_SOURCE_MAX_CHARS + 1];
    char            type[BL_GFF3_TYPE_MAX_CHARS + 1];
    int64_t         start,
		    end;
    double          score;
    char            strand;         // '+' or '-' or '.'
    char            phase;          // 0, 1, 2, or '.' (bases to condon start)
    char            *attributes;    // Ensembl ID, Name, etc.
    size_t          attributes_array_size;
    size_t          attributes_len;
    
    /*
     *  Fields below are not part of GFF3.  They are extracted from attributes
     *  because they are frequently useful.
     */
    char            *feature_id;    // In every feature of Ensemble GFF3s
    char            *feature_name;  // Extract from gene features and look
				    // up using Ensemble ID for others
    char            *feature_parent;    // Transcripts, exons, etc.
    
    // Offset of the feature in the GFF3 file for indexing
    long            file_pos;
}   bl_gff3_t;

typedef unsigned int            gff3_field_mask_t;

#define BL_GFF3_FIELD_SEQID      0x001
#define BL_GFF3_FIELD_SOURCE     0x002
#define BL_GFF3_FIELD_TYPE       0x004
#define BL_GFF3_FIELD_START      0x008
#define BL_GFF3_FIELD_END        0x010
#define BL_GFF3_FIELD_SCORE      0x020
#define BL_GFF3_FIELD_STRAND     0x040
#define BL_GFF3_FIELD_PHASE      0x080
#define BL_GFF3_FIELD_ATTRIBUTES 0x100
#define BL_GFF3_FIELD_ALL        0xfff

// After bl_gff3_t for prototypes
#ifndef _BIOLIBC_BED_H_
#include "bed.h"
#endif

#include "gff3-rvs.h"
#include "gff3-accessors.h"
#include "gff3-mutators.h"

// After bl_sam_t def due to mutual recursion
#ifndef _BIOLIBC_SAM_H_
#include "sam.h"
#endif

/* gff3.c */
FILE *bl_gff3_skip_header(FILE *gff3_stream);
int bl_gff3_copy_header(FILE *header_stream, FILE *gff3_stream);
int bl_gff3_read(bl_gff3_t *gff3_feature, FILE *gff3_stream, gff3_field_mask_t field_mask);
int bl_gff3_write(bl_gff3_t *gff3_feature, FILE *gff3_stream, gff3_field_mask_t field_mask);
void bl_gff3_to_bed(bl_gff3_t *gff3_feature, bl_bed_t *bed_feature);
void bl_gff3_free(bl_gff3_t *gff3_feature);
char *bl_gff3_extract_attribute(bl_gff3_t *feature, const char *attr_name);
void bl_gff3_init(bl_gff3_t *feature);
bl_gff3_t *bl_gff3_dup(bl_gff3_t *feature);
bl_gff3_t *bl_gff3_copy(bl_gff3_t *copy, bl_gff3_t *feature);
int bl_gff3_sam_cmp(bl_gff3_t *feature, bl_sam_t *alignment);
int64_t bl_gff3_sam_overlap(bl_gff3_t *feature, bl_sam_t *alignment);

#ifdef __cplusplus
}
#endif

#endif  // _BIOLIBC_GFF3_H_
