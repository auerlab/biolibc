#include <string.h>
#include <sys/stat.h>
#include <xtend/string.h>      // strlcpy() on Linux
#include "biolibc.h"
#include "overlap.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Set all fields in a bl_overlap_t structure.  Start and end
 *      positions are 1-based regardless of the feature type.  (BED file
 *      positions are 0-based and must be adjusted before passed to this
 *      function.)
 *
 *  Arguments:
 *      feature1_len      Length of feature 1
 *      feature2_len      Length of feature 2
 *      overlap_start    Start position of overlap relative to start of feature 1
 *      overlap_end      End position of overlap relative to start of feature 1
 *
 *  Returns:
 *      BL_DATA_OK upon success.
 *      BL_DATA_INVALID is arguments don't make sense.
 *
 *  Examples:
 *          bed_start = BL_BED_CHROM_START(bed_feature);
 *          bed_end = BL_BED_CHROM_END(bed_feature);
 *          gff3_start = BL_GFF3_CHROM_START(gff3_feature);
 *          gff3_end = BL_GFF3_CHROM_END(gff3_feature);
 *          bed_len = bed_end - bed_start;
 *          gff3_len = gff3_end - gff3_start + 1;
 *          bl_overlap_set_all(overlap, bed_len, gff3_len,
 *                          XT_MAX(bed_start+1, gff3_start),
 *                          XT_MIN(bed_end, gff3_end));
 *
 *  See also:
 *      bl_overlap_print(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_overlap_set_all(bl_overlap_t *overlap,
			int64_t feature1_len, int64_t feature2_len,
			int64_t overlap_start, int64_t overlap_end)

{
    overlap->feature1_len = feature1_len;
    overlap->feature2_len = feature2_len;
    overlap->overlap_start = overlap_start;
    overlap->overlap_end = overlap_end;
    overlap->overlap_len = overlap_end - overlap_start + 1;
    
    // FIXME: Return BL_DATA_INVALID if sanity checks fail
    return BL_OVERLAP_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Print all fields in a bl_overlap_t structure for debugging.
 *
 *  Arguments:
 *      stream      FILE stream to which data are printed (e.g. stderr)
 *      overlap     Address of a bl_overlap_t structure
 *      feature1_name     Name of field 1 to print with data
 *      feature2_name     Name of field 2 to print with data
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

int     bl_overlap_print(bl_overlap_t *overlap, FILE *stream,
			  char *feature1_name, char *feature2_name)

{
    char    feature1_len[16], feature2_len[16];
    
    strlcpy(feature1_len, feature1_name, 12);
    strlcat(feature1_len, " len", 16);
    strlcpy(feature2_len, feature2_name, 12);
    strlcat(feature2_len, " len", 16);
    return fprintf(stream, "%-16s: %" PRId64 "\n"
	   "%-16s: %" PRId64 "\n"
	   "Overlap start   : %" PRId64 "\n"
	   "Overlap end     : %" PRId64 "\n"
	   "Overlap length  : %" PRId64 "\n",
	   feature1_len, overlap->feature1_len,
	   feature2_len, overlap->feature2_len,
	   overlap->overlap_start, overlap->overlap_end, overlap->overlap_len);
}
