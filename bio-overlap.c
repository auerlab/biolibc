#include <string.h>
#include <sys/stat.h>
#include <xtend.h>      // strlcpy() on Linux
#include "bio-overlap.h"

/***************************************************************************
 *  Description:
 *      Set overlap parameters
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

void    bio_overlap_set_all(bio_overlap_t *overlap,
			uint64_t f1_len, uint64_t f2_len,
			uint64_t ov_start, uint64_t ov_end)

{
    overlap->f1_len = f1_len;
    overlap->f2_len = f2_len;
    overlap->ov_start = ov_start;
    overlap->ov_end = ov_end;
    overlap->ov_len = ov_end - ov_start + 1;
}


/***************************************************************************
 *  Description:
 *      Print overlap parameters
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  Files:
 *
 *  Environment:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

void    bio_overlap_print(bio_overlap_t *overlap, char *f1_name, char *f2_name)

{
    char    f1_len[16], f2_len[16];
    
    strlcpy(f1_len, f1_name, 12);
    strlcat(f1_len, " len", 16);
    strlcpy(f2_len, f2_name, 12);
    strlcat(f2_len, " len", 16);
    printf("%-16s: %" PRIu64 "\n"
	   "%-16s: %" PRIu64 "\n"
	   "Overlap start   : %" PRIu64 "\n"
	   "Overlap end     : %" PRIu64 "\n"
	   "Overlap length  : %" PRIu64 "\n",
	   f1_len, overlap->f1_len,
	   f2_len, overlap->f2_len,
	   overlap->ov_start, overlap->ov_end, overlap->ov_len);
}
