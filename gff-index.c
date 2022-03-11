#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>  // MIN()
#include <xtend/mem.h>
#include <xtend/math.h>
#include "gff-index.h"

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      The gff_index_t class maintains an in-memory index of GFF
 *      features, containing the GFF fields SEQ_ID, START, and END,
 *      and the offset into the file as reported by ftell(3), or by
 *      bl_gff_read(3), which records the file position of each GFF
 *      feature it reads.
 *
 *      .B bl_gff_index_add_pos(3)
 *      adds a GFF feature with file position file_pos to the index.
 *      Features of interest, perhaps only genes or only exons, can
 *      be added to the index on-the fly while reading through a GFF
 *      file with bl_gff_read(3).
 *
 *      The index can later be searched or traversed forward or backward
 *      to quickly find
 *      the location of a feature and reposition a FILE pointer to it
 *      using fseek(3).  This system eliminates the need to inhale
 *      large numbers of GFF features into memory.
 *  
 *  Arguments:
 *      gi      Pointer to gff_index_t object to which a record will be added
 *      feature Pointer to GFF feature to be indexed
 *
 *  Returns:
 *      BL_GFF_INDEX_OK on success, BL_GFF_MALLOC_FAILED if memory could
 *      not be allocated
 *
 *  Examples:
 *      bl_gff_index_t  gi;
 *      bl_gff_t        feature;
 *
 *      if ( bl_gff_read(&feature, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK )
 *      {
 *          if ( bl_gff_index_add(&gi, &feature) != BL_GFF_INDEX_OK )
 *              fprintf(stderr, "Error addind to GFF index.\n");
 *      }
 *
 *  See also:
 *      bl_gff_index_seek_reverse(3), bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_index_add(bl_gff_index_t *gi, bl_gff_t *feature)

{
    if ( gi->count == gi->array_size )
    {
	gi->array_size += 65536;
	gi->file_pos = xt_realloc(gi->file_pos, gi->array_size, sizeof(*gi->file_pos));
	if ( gi->file_pos == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->start = xt_realloc(gi->start, gi->array_size, sizeof(*gi->start));
	if ( gi->start == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->end = xt_realloc(gi->end, gi->array_size, sizeof(*gi->end));
	if ( gi->end == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->seqid = xt_realloc(gi->seqid, gi->array_size, sizeof(*gi->seqid));
	if ( gi->seqid == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
    }
    gi->file_pos[gi->count] = BL_GFF_FILE_POS(feature);
    gi->start[gi->count] = BL_GFF_START(feature);
    gi->end[gi->count] = BL_GFF_END(feature);
    
    // Hash chr to a 64-bit integer, padding with 0s beyond the end
    if ( (gi->seqid[gi->count] = strdup(BL_GFF_SEQID(feature))) == NULL )
	return BL_GFF_INDEX_MALLOC_FAILED;
    ++gi->count;
    return BL_GFF_INDEX_OK;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      The gff_index_t class maintains an in-memory index of GFF
 *      features, containing the GFF fields SEQ_ID, START, and END,
 *      and the offset into the file as reported by ftell(3), or by
 *      bl_gff_read(3), which records the file position of each GFF
 *      feature it reads.
 *
 *      .B bl_gff_index_seek_reverse(3)
 *      moved the FILE pointer stream to feature_count indexed features
 *      upstream of feature, or to the most upstream feature within
 *      max_nt of feature.  A max_nt of 0 indicates no maximum distance,
 *      i.e. the search will proceed to feature_count features behind
 *      feature or to the beginning of the file, whichever is encontered
 *      first.
 *
 *      The max_nt parameter refers to the END of a feature, i.e.
 *      .B bl_gff_index_seek_reverse()
 *      will back up to a feature that overlaps the position of feature
 *      minus max_nt.  The START position of the feature moved to could
 *      be more than max_nt nucleotides behind the START of feature.
 *
 *      Note that this function counts only *indexed* features, i.e. those
 *      added to gi by bl_gff_index_add(3), not all features in the GFF
 *      file.  An application may only add genes to the index, for example,
 *      ignoring exons, etc.
 *  
 *  Arguments:
 *      gi              Pointer to the gff_index_t object used to search
 *      feature         Feature from which search starts
 *      feature_count   Number of indexed features to back up from feature
 *      max_nt          Maximum number of nucleotides to back up
 *
 *  Returns:
 *      The offset into stream of the feature moved to.  This is only
 *      informative, since the fseek() operation has already been done.
 *
 *  Examples:
 *      bl_gff_index_t  gi;
 *      bl_gff_t        feature;
 *      long            new_pos;
 *
 *      new_pos = bl_gff_index_seek_reverse(&gi, &feature, 4, 200000);
 *
 *  See also:
 *      bl_gff_index_add(3), fseek(3), bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_index_seek_reverse(bl_gff_index_t *gi, FILE *stream,
	    bl_gff_t *feature, int64_t feature_count, int64_t max_nt)

{
    ssize_t     c;
    char        *ref_seqid = BL_GFF_SEQID(feature);
    int64_t     ref_start = BL_GFF_START(feature),
		end = MAX(ref_start - max_nt, 0),
		f;

    // First find the reference feature where the search begins
    for (c = gi->count - 1; (c >= 0) &&
			    (gi->start[c] != ref_start) &&
			    (strcmp(gi->seqid[c], ref_seqid) != 0); --c)
	;
    
    // Now back up feature_count features or to the leftmost feature
    // overlapping with the ref feature start - max_nt
    for (f = feature_count; (f > 0) &&
			    (c > 0) &&
			    (strcmp(gi->seqid[c], ref_seqid) == 0) &&
			    (gi->end[c] > end); --f, --c)
	;
    if ( gi->end[c] < end )
	++c;

    return fseek(stream, gi->file_pos[c], SEEK_SET);
}

