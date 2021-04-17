
#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDINT_H_
#include <stdint.h>
#endif

#define PLIST_INIT  { 0, 0, NULL };

#define PLIST_MAX_POSITIONS(p)  ((p)->max_positions)
#define PLIST_COUNT(p)          ((p)->count)
#define PLIST_POSITIONS(p, c)   ((p)->positions[c])

typedef enum
{
    PLIST_ASCENDING,
    PLIST_DESCENDING
}   plist_sort_order_t;

typedef struct
{
    size_t      max_positions;
    size_t      count;
    uint64_t    *positions;
}   plist_t;

/* plist.c */
void plist_allocate(plist_t *plist, size_t max_positions);
void plist_free(plist_t *plist);
int plist_add_position(plist_t *plist, uint64_t position);
int plist_from_csv(plist_t *plist, const char *bounds_str, size_t max_positions);
void    plist_sort(plist_t *plist, plist_sort_order_t order);
