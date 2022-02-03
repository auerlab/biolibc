#ifndef _BIOLIBC_GFF_H__
#define _BIOLIBC_GFF_H__

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

#define BL_GFF_SCORE_MAX_DIGITS        64       // Floating point
#define BL_GFF_SOURCE_MAX_CHARS        1024     // Guess
#define BL_GFF_FEATURE_MAX_CHARS       256      // Guess
#define BL_GFF_STRAND_MAX_CHARS        2
#define BL_GFF_LINE_MAX_CHARS          4096
#define BL_GFF_PHASE_MAX_DIGITS        2
#define BL_GFF_ATTRIBUTES_MAX_CHARS    8192     // For temp vars only.
						// Structure uses malloc()

#define BL_GFF_SCORE_UNAVAILABLE       -1.0
#define BL_GFF_PHASE_UNAVAILABLE       '.'

#define BL_GFF_INIT \
	{ "", "", "", 0, 0, 0.0, '.', '.', NULL, NULL, NULL }

typedef struct
{
    char            sequence[BL_CHROM_MAX_CHARS + 1];
    char            source[BL_GFF_SOURCE_MAX_CHARS + 1];
    char            feature[BL_GFF_FEATURE_MAX_CHARS + 1];
    uint64_t        start,
		    end;
    double          score;
    char            strand;         // '+' or '-' or '.'
    char            phase;          // 0, 1, 2, or '.' (bases to condon start)
    char            *attributes;    // Ensembl ID, Name, etc.
    
    /*
     *  Fields below are not part of GFF.  They are extracted from attributes
     *  and may be useful.
     */
    char            *feature_id;    // In every feature of Ensemble GFFs
    char            *gene_name;     // Extract from gene features and look
				    // up using Ensemble ID for others

    long            file_pos;       // Offset in the GFF file
}   bl_gff_t;

typedef unsigned int            gff_field_mask_t;

#define BL_GFF_FIELD_SEQUENCE   0x001
#define BL_GFF_FIELD_SOURCE     0x002
#define BL_GFF_FIELD_FEATURE    0x004
#define BL_GFF_FIELD_START      0x008
#define BL_GFF_FIELD_END        0x010
#define BL_GFF_FIELD_SCORE      0x020
#define BL_GFF_FIELD_STRAND     0x040
#define BL_GFF_FIELD_PHASE      0x080
#define BL_GFF_FIELD_ATTRIBUTES 0x100
#define BL_GFF_FIELD_ALL        0xfff

// After bl_gff_t for prototypes
#ifndef _BIOLIBC_BED_H_
#include "bed.h"
#endif

#include "gff-rvs.h"
#include "gff-accessors.h"
#include "gff-mutators.h"

/* gff.c */
FILE *bl_gff_skip_header(FILE *gff_stream);
int bl_gff_copy_header(FILE *header_stream, FILE *gff_stream);
int bl_gff_read(bl_gff_t *gff_feature, FILE *gff_stream, gff_field_mask_t field_mask);
int bl_gff_write(bl_gff_t *gff_feature, FILE *gff_stream, gff_field_mask_t field_mask);
void bl_gff_to_bed(bl_gff_t *gff_feature, bl_bed_t *bed_feature);
void bl_gff_free(bl_gff_t *gff_feature);

#endif  // _BIOLIBC_GFF_H__
