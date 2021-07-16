#include <string.h>
#include <sys/stat.h>
#include <xtend.h>      // strlcpy() on Linux
#include "biolibc.h"
#include "overlap.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc
 *
 *  Description:
 *      Set all fields in a bl_overlap_t structure.  Start and end
 *      positions are 1-based regardless of the feature type.  (BED file
 *      positions are 0-based and must be adjusted before passed to this
 *      function.)
 *
 *  Arguments:
 *      feature1_len:     Length of feature 1
 *      feature2_len:     Length of feature 2
 *      overlap_start:   Start position of overlap relative to start of feature 1
 *      overlap_end:     End position of overlap relative to start of feature 1
 *
 *  Returns:
 *      BL_DATA_OK upon success.
 *      BL_DATA_INVALID is arguments don't make sense.
 *
 *  Examples:
 *          bed_start = BL_BED_START_POS(bed_feature);
 *          bed_end = BL_BED_END_POS(bed_feature);
 *          gff_start = BL_GFF_START_POS(gff_feature);
 *          gff_end = BL_GFF_END_POS(gff_feature);
 *          bed_len = bed_end - bed_start;
 *          gff_len = gff_end - gff_start + 1;
 *          bl_overlap_set_all(overlap, bed_len, gff_len,
 *                          MAX(bed_start+1, gff_start),
 *                          MIN(bed_end, gff_end));
 *
 *  See also:
 *      bl_overlap_print(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_overlap_set_all(bl_overlap_t *overlap,
			uint64_t feature1_len, uint64_t feature2_len,
			uint64_t overlap_start, uint64_t overlap_end)

{
    overlap->feature1_len = feature1_len;
    overlap->feature2_len = feature2_len;
    overlap->overlap_start = overlap_start;
    overlap->overlap_end = overlap_end;
    overlap->overlap_len = overlap_end - overlap_start + 1;
    
    // FIXME: Return BL_DATA_INVALID if sanity checks fail
    return BL_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc
 *
 *  Description:
 *      Print all fields in a bl_overlap_t structure for debugging.
 *
 *  Arguments:
 *      stream:     FILE stream to which data are printed (e.g. stderr)
 *      overlap:    Address of a bl_overlap_t structure
 *      feature1_name:    Name of field 1 to print with data
 *      feature2_name:    Name of field 2 to print with data
 *
 *  Returns:
 *      Return status from fprintf(3)
 *
 *  See also:
 *      bl_overlap_set_all(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_overlap_print(FILE *stream, bl_overlap_t *overlap,
			  char *feature1_name, char *feature2_name)

{
    char    feature1_len[16], feature2_len[16];
    
    strlcpy(feature1_len, feature1_name, 12);
    strlcat(feature1_len, " len", 16);
    strlcpy(feature2_len, feature2_name, 12);
    strlcat(feature2_len, " len", 16);
    return fprintf(stream, "%-16s: %" PRIu64 "\n"
	   "%-16s: %" PRIu64 "\n"
	   "Overlap start   : %" PRIu64 "\n"
	   "Overlap end     : %" PRIu64 "\n"
	   "Overlap length  : %" PRIu64 "\n",
	   feature1_len, overlap->feature1_len,
	   feature2_len, overlap->feature2_len,
	   overlap->overlap_start, overlap->overlap_end, overlap->overlap_len);
}
