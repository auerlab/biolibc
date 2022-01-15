
#ifndef _BIOLIBC_POS_LIST_H_
#define _BIOLIBC_POS_LIST_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDINT_H_
#include <stdint.h>
#endif

#ifndef _BIOLIBC_H_
#include "biolibc.h"
#endif

#define BL_POS_LIST_INIT  { 0, 0, NULL };

typedef struct
{
    size_t      array_size;
    size_t      count;
    uint64_t    *positions;
}   bl_pos_list_t;

typedef int bl_pos_list_sort_order_t;

#define BL_POS_LIST_ASCENDING  0
#define BL_POS_LIST_DESCENDING 1

#include "pos-list-rvs.h"
#include "pos-list-accessors.h"
#include "pos-list-mutators.h"

/* pos-list.c */
void bl_pos_list_allocate(bl_pos_list_t *pos_list, size_t array_size);
void bl_pos_list_free(bl_pos_list_t *pos_list);
int bl_pos_list_add_position(bl_pos_list_t *pos_list, uint64_t position);
int bl_pos_list_from_csv(bl_pos_list_t *pos_list, const char *bounds_str, size_t array_size);
int position_cmp_ascending(const uint64_t *pos1, const uint64_t *pos2);
int position_cmp_descending(const uint64_t *pos1, const uint64_t *pos2);
void bl_pos_list_sort(bl_pos_list_t *pos_list, bl_pos_list_sort_order_t order);

#endif  // _BIOLIBC_POS_LIST_H_
