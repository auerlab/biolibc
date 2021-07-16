#ifndef _gff_h_
#define _gff_h_

#ifndef _biolibc_h_
#include "biolibc.h"
#endif

#ifndef _xtend_h_
#include <xtend.h>  // strlcpy() on Linux
#endif

#define BL_GFF_SCORE_MAX_DIGITS        64      // Floating point
#define BL_GFF_SOURCE_MAX_CHARS        1024    // Guess
#define BL_GFF_NAME_MAX_CHARS          1024    // Guess
#define BL_GFF_STRAND_MAX_CHARS        2
#define BL_GFF_LINE_MAX_CHARS          4096
#define BL_GFF_PHASE_MAX_DIGITS        2
#define BL_GFF_ATTRIBUTES_MAX_CHARS    8192    // For temp vars only.
					    // Structure uses malloc()

#define BL_GFF_SCORE_UNAVAILABLE       -1.0
#define BL_GFF_PHASE_UNAVAILABLE       '.'

typedef unsigned int        gff_field_mask_t;

#define BL_GFF_FIELD_SEQUENCE      0x001
#define BL_GFF_FIELD_SOURCE        0x002
#define BL_GFF_FIELD_NAME          0x004
#define BL_GFF_FIELD_START_POS     0x008
#define BL_GFF_FIELD_END_POS       0x010
#define BL_GFF_FIELD_SCORE         0x020
#define BL_GFF_FIELD_STRAND        0x040
#define BL_GFF_FIELD_PHASE         0x080
#define BL_GFF_FIELD_ATTRIBUTES    0x100
#define BL_GFF_FIELD_ALL           0xfff

#define BL_GFF_SEQUENCE(ptr)    ((ptr)->sequence)
#define BL_GFF_SOURCE(ptr)      ((ptr)->source)
#define BL_GFF_NAME(ptr)        ((ptr)->name)
#define BL_GFF_START_POS(ptr)   ((ptr)->start_pos)
#define BL_GFF_END_POS(ptr)     ((ptr)->end_pos)
#define BL_GFF_SCORE(ptr)       ((ptr)->score)
#define BL_GFF_STRAND(ptr)      ((ptr)->strand)
#define BL_GFF_PHASE(ptr)       ((ptr)->phase)
#define BL_GFF_ATTRIBUTES(ptr)  ((ptr)->attributes)
#define BL_GFF_FEATURE_ID(ptr)  ((ptr)->feature_id)
#define BL_GFF_GENE_NAME(ptr)   ((ptr)->gene_name)

#define BL_GFF_SET_SEQUENCE(ptr,v)  strlcpy((ptr)->sequence,v,BL_CHROMOSOME_MAX_CHARS+1)
#define BL_GFF_SET_SOURCE(ptr,v)    strlcpy((ptr)->source,v,BL_GFF_SOURCE_MAX_CHARS+1)
#define BL_GFF_SET_NAME(ptr,v)      strlcpy((ptr)->name,v,BL_GFF_NAME_MAX_CHARS+1)
#define BL_GFF_SET_START_POS(ptr,v) ((ptr)->start_pos = (v))
#define BL_GFF_SET_END_POS(ptr,v)   ((ptr)->end_pos = (v))
#define BL_GFF_SET_SCORE(ptr,v)     ((ptr)->score = (v))
#define BL_GFF_SET_STRAND(ptr,v)    ((ptr)->strand = (v))
#define BL_GFF_SET_PHASE(ptr,v)     ((ptr)->phase = (v))
#define BL_GFF_SET_ATTRIBUTES(ptr,v)    ((ptr)->attributes = (v))
#define BL_GFF_SET_FEATURE_ID(ptr,v)    ((ptr)->feature_id = (v))
#define BL_GFF_SET_GENE_NAME(ptr,v) ((ptr)->gene_name = (v))

#define BL_GFF_INIT \
	{ "", "", "", 0, 0, 0.0, '.', '.', NULL, NULL, NULL }

typedef struct
{
    char            sequence[BL_CHROMOSOME_MAX_CHARS + 1];
    char            source[BL_GFF_SOURCE_MAX_CHARS + 1];
    char            name[BL_GFF_NAME_MAX_CHARS + 1];
    uint64_t        start_pos,
		    end_pos;
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
}   bl_gff_t;

// After bl_gff_t for prototypes
#ifndef _bed_h_
#include "bed.h"
#endif

FILE *gff_skip_header(FILE *gff_stream);
int gff_read_feature(FILE *gff_stream, bl_gff_t *gff_feature, gff_field_mask_t field_mask);
int gff_write_feature(FILE *gff_stream, bl_gff_t *gff_feature, gff_field_mask_t field_mask);
void gff_to_bed(bl_bed_t *bed_feature, bl_gff_t *gff_feature);

#endif  // _gff_h_
