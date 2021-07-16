
#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDINT_H_
#include <stdint.h>
#endif

#define BL_POS_LIST_INIT  { 0, 0, NULL };

#define BL_POS_LIST_ARRAY_SIZE(ptr)     ((ptr)->array_size)
#define BL_POS_LIST_COUNT(ptr)          ((ptr)->count)
#define BL_POS_LIST_POSITIONS(ptr,c)    ((ptr)->positions[c])

#define BL_POS_LIST_SET_ARRAY_SIZE(ptr,v)   ((ptr)->array_size = (v))
#define BL_POS_LIST_SET_COUNT(ptr,v)        ((ptr)->count = (v))
#define BL_POS_LIST_SET_POSITIONS(ptr,v)    ((ptr)->positions = (v))

#define BL_POS_LIST_ASCENDING  0
#define BL_POS_LIST_DESCENDING 1

typedef int bl_pos_list_sort_order_t;

typedef struct
{
    size_t      array_size;
    size_t      count;
    uint64_t    *positions;
}   bl_pos_list_t;

/* pos-list.c */
void bl_pos_list_allocate(bl_pos_list_t *pos_list, size_t max_positions);
void bl_pos_list_free(bl_pos_list_t *pos_list);
int bl_pos_list_add_position(bl_pos_list_t *pos_list, uint64_t position);
int bl_pos_list_from_csv(bl_pos_list_t *pos_list, const char *bounds_str, size_t max_positions);
void bl_pos_list_sort(bl_pos_list_t *pos_list, bl_pos_list_sort_order_t order);
