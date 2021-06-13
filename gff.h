#ifndef _gff_h_
#define _gff_h_

#ifndef _biolibc_h_
#include "biolibc.h"
#endif

#ifndef _xtend_h_
#include <xtend.h>  // strlcpy() on Linux
#endif

#define GFF_SCORE_MAX_DIGITS        64      // Floating point
#define GFF_SOURCE_MAX_CHARS        1024    // Guess
#define GFF_NAME_MAX_CHARS          1024    // Guess
#define GFF_STRAND_MAX_CHARS        2
#define GFF_LINE_MAX_CHARS          4096
#define GFF_PHASE_MAX_DIGITS        2
#define GFF_ATTRIBUTES_MAX_CHARS    8192    // For temp vars only.
					    // Structure uses malloc()

#define GFF_SCORE_UNAVAILABLE       -1.0
#define GFF_PHASE_UNAVAILABLE       '.'

typedef unsigned int        gff_field_mask_t;

#define GFF_FIELD_ALL       0x0
#define GFF_FIELD_CHROM     0x1
#define GFF_FIELD_START_POS 0x2
#define GFF_FIELD_END_POS   0x4

#define GFF_SEQUENCE(gf)        ((gf)->sequence)
#define GFF_START_POS(gf)       ((gf)->start_pos)
#define GFF_END_POS(gf)         ((gf)->end_pos)
#define GFF_NAME(gf)            ((gf)->name)
#define GFF_SCORE(gf)           ((gf)->score)
#define GFF_STRAND(gf)          ((gf)->strand)

#define GFF_SET_START_POS(gf, p)    ((gf)->start_pos = (p))
#define GFF_SET_END_POS(gf, p)      ((gf)->end_pos = (p))
#define GFF_SET_NAME(gf, n) \
	(strlcpy((gf)->name, (n), GFF_NAME_MAX_CHARS))

#define GFF_INIT \
	{ "", "", "", "", "", 0, 0, "", 0.0, '.', '.', NULL, NULL, NULL }

typedef struct
{
    char            sequence[BIO_CHROMOSOME_MAX_CHARS + 1];
    char            source[GFF_SOURCE_MAX_CHARS + 1];
    char            name[GFF_NAME_MAX_CHARS + 1];
    char            start_pos_str[BIO_POSITION_MAX_DIGITS + 1],  // 0-based
		    end_pos_str[BIO_POSITION_MAX_DIGITS + 1];   
    uint64_t        start_pos,
		    end_pos;
    char            score_str[GFF_SCORE_MAX_DIGITS + 1];
    double          score;
    char            strand;         // '+' or '-' or '.'
    char            phase;          // 0, 1, 2, or '.' (bases to condon start)
    char            *attributes;    // Ensembl ID, Name, etc.
    
    /*
     *  Fields below are not part of GFF.  They are extracted from attributes
     *  and may be useful.
     */
    char            *feature_id;    // In every feature of Ensemble GFFs
    char            *gane_name;     // Extract from gene features and look
				    // up using Ensemble ID for others
}   gff_feature_t;

// After gff_feature_t for prototypes
#ifndef _bed_h_
#include "bed.h"
#endif

FILE *gff_skip_header(FILE *gff_stream);
int gff_read_feature(FILE *gff_stream, gff_feature_t *gff_feature);
int gff_write_feature(FILE *gff_stream, gff_feature_t *gff_feature, gff_field_mask_t field_mask);
void gff_to_bed(bed_feature_t *bed_feature, gff_feature_t *gff_feature);

#endif  // _gff_h_
