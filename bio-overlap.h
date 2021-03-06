#ifndef _bio_overlap_h_

#ifndef _INTTYPES_H_    // Save needless rereading
#include <inttypes.h>
#endif

// 1-based, inclusive at both ends
typedef struct
{
    uint64_t    f1_len,
		f2_len,
		ov_start,
		ov_end,
		ov_len;
}   bio_overlap_t;

/* chromosome-name-cmp.c */
int bio_overlap_set_all(bio_overlap_t *overlap,
			uint64_t f1_len, uint64_t f2_len,
			uint64_t ov_start, uint64_t ov_end);
int bio_overlap_print(FILE *stream, bio_overlap_t *overlap,
		      char *f1_name, char *f2_name);

#endif // _bio_overlap_h_
