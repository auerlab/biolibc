
#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDINT_H_
#include <stdint.h>
#endif

#define POS_LIST_INIT  { 0, 0, NULL };

#define POS_LIST_ARRAY_SIZE(p)      ((p)->array_size)
#define POS_LIST_COUNT(p)           ((p)->count)
#define POS_LIST_POSITIONS(p, c)    ((p)->positions[c])

typedef enum
{
    POS_LIST_ASCENDING,
    POS_LIST_DESCENDING
}   pos_list_sort_order_t;

typedef struct
{
    size_t      array_size;
    size_t      count;
    uint64_t    *positions;
}   pos_list_t;

/* pos-list.c */
void pos_list_allocate(pos_list_t *pos_list, size_t max_positions);
void pos_list_free(pos_list_t *pos_list);
int pos_list_add_position(pos_list_t *pos_list, uint64_t position);
int pos_list_from_csv(pos_list_t *pos_list, const char *bounds_str, size_t max_positions);
void pos_list_sort(pos_list_t *pos_list, pos_list_sort_order_t order);
