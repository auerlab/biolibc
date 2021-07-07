#include <string.h>
#include <sys/stat.h>
#include <xtend.h>      // strlcpy() on Linux
#include "biolibc.h"
#include "bio-overlap.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/bio-overlap.h>
 *      -lbiolibc
 *
 *  Description:
 *      Set all fields in a bio_overlap_t structure.  Start and end
 *      positions are 1-based regardless of the feature type.  (BED file
 *      positions are 0-based and must be adjusted before passed to this
 *      function.)
 *
 *  Arguments:
 *      f1_len:     Length of feature 1
 *      f2_len:     Length of feature 2
 *      ov_start:   Start position of overlap relative to start of feature 1
 *      ov_end:     End position of overlap relative to start of feature 1
 *
 *  Returns:
 *      BIO_DATA_OK upon success.
 *      BIO_DATA_INVALID is arguments don't make sense.
 *
 *  Examples:
 *          bed_start = BED_START_POS(bed_feature);
 *          bed_end = BED_END_POS(bed_feature);
 *          gff_start = GFF_START_POS(gff_feature);
 *          gff_end = GFF_END_POS(gff_feature);
 *          bed_len = bed_end - bed_start;
 *          gff_len = gff_end - gff_start + 1;
 *          bio_overlap_set_all(overlap, bed_len, gff_len,
 *                          MAX(bed_start+1, gff_start),
 *                          MIN(bed_end, gff_end));
 *
 *  See also:
 *      bio_overlap_print(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bio_overlap_set_all(bio_overlap_t *overlap,
			uint64_t f1_len, uint64_t f2_len,
			uint64_t ov_start, uint64_t ov_end)

{
    overlap->f1_len = f1_len;
    overlap->f2_len = f2_len;
    overlap->ov_start = ov_start;
    overlap->ov_end = ov_end;
    overlap->ov_len = ov_end - ov_start + 1;
    
    // FIXME: Return BIO_DATA_INVALID if sanity checks fail
    return BIO_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bio-overlap.h>
 *      -lbiolibc
 *
 *  Description:
 *      Print all fields in a bio_overlap_t structure for debugging.
 *
 *  Arguments:
 *      stream:     FILE stream to which data are printed (e.g. stderr)
 *      overlap:    Address of a bio_overlap_t structure
 *      f1_name:    Name of field 1 to print with data
 *      f2_name:    Name of field 2 to print with data
 *
 *  Returns:
 *      Return status from fprintf(3)
 *
 *  See also:
 *      bio_overlap_set_all(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bio_overlap_print(FILE *stream, bio_overlap_t *overlap,
			  char *f1_name, char *f2_name)

{
    char    f1_len[16], f2_len[16];
    
    strlcpy(f1_len, f1_name, 12);
    strlcat(f1_len, " len", 16);
    strlcpy(f2_len, f2_name, 12);
    strlcat(f2_len, " len", 16);
    return fprintf(stream, "%-16s: %" PRIu64 "\n"
	   "%-16s: %" PRIu64 "\n"
	   "Overlap start   : %" PRIu64 "\n"
	   "Overlap end     : %" PRIu64 "\n"
	   "Overlap length  : %" PRIu64 "\n",
	   f1_len, overlap->f1_len,
	   f2_len, overlap->f2_len,
	   overlap->ov_start, overlap->ov_end, overlap->ov_len);
}
