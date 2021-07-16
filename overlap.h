#ifndef _bl_overlap_h_

#ifndef _INTTYPES_H_    // Save needless rereading
#include <inttypes.h>
#endif

#define BL_OVERLAP_FEATURE1_LEN(ptr)    ((ptr)->feature1_len)
#define BL_OVERLAP_FEATURE2_LEN(ptr)    ((ptr)->feature2_len)
#define BL_OVERLAP_OVERLAP_START(ptr)   ((ptr)->overlap_start)
#define BL_OVERLAP_OVERLAP_END(ptr)     ((ptr)->overlap_end)
#define BL_OVERLAP_OVERLAP_LEN(ptr)     ((ptr)->overlap_len)

#define BL_OVERLAP_SET_FEATURE1_LEN(ptr,feature1_len)   ((ptr)->feature1_len = (feature1_len))
#define BL_OVERLAP_SET_FEATURE2_LEN(ptr,feature2_len)   ((ptr)->feature2_len = (feature2_len))
#define BL_OVERLAP_SET_OVERLAP_START(ptr,overlap_start) ((ptr)->overlap_start = (overlap_start))
#define BL_OVERLAP_SET_OVERLAP_END(ptr,overlap_end)     ((ptr)->overlap_end = (overlap_end))
#define BL_OVERLAP_SET_OVERLAP_LEN(ptr,overlap_len)     ((ptr)->overlap_len = (overlap_len))

// 1-based, inclusive at both ends
typedef struct
{
    uint64_t    feature1_len;
    uint64_t    feature2_len;
    uint64_t    overlap_start;
    uint64_t    overlap_end;
    uint64_t    overlap_len;
}   bl_overlap_t;

/* chromosome-name-cmp.c */
int bl_overlap_set_all(bl_overlap_t *overlap,
			uint64_t f1_len, uint64_t f2_len,
			uint64_t ov_start, uint64_t ov_end);
int bl_overlap_print(FILE *stream, bl_overlap_t *overlap,
		      char *f1_name, char *f2_name);

#endif // _bl_overlap_h_
