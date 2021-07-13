#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <string.h>
#include <xtend.h>
#include "pos-list.h"
#include "biolibc.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc
 *
 *  Description:
 *      Initialize a position list with an initial array size of
 *      array_size.  The size will be increased by pos_list_add_position()
 *      if necessary.
 *
 *  Arguments:
 *      pos_list:   Pointer to the pos_list_t structure to initialize
 *      array_size: Initial size of the array of positions
 *
 *  See also:
 *      pos_list_add_position(3), pos_list_free(3), pos_list_from_csv(3),
 *      pos_list_sort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    pos_list_allocate(pos_list_t *pos_list, size_t array_size)

{
    if ( (pos_list->count != 0) || (pos_list->array_size != 0) ||
	 (pos_list->positions != NULL) )
    {
	fputs("pos_list_allocate(): List is not blank.\n", stderr);
	fputs("Was it previously allocated?\n", stderr);
	fputs("Did you forget to initialize it with POS_LIST_INIT?\n", stderr);
	exit(EX_SOFTWARE);
    }
    pos_list->positions = xt_malloc(array_size, sizeof(uint64_t));
    if ( pos_list->positions == NULL )
    {
	fputs("pos_list_allocate(): Could not allocate positions.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    pos_list->array_size = array_size;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc
 *
 *  Description:
 *      Free the array of positions in position list and reset array_size
 *      and position count to 0.
 *
 *  Arguments:
 *      pos_list:   Pointer to the position_list_t structure to reset
 *
 *  See also:
 *      pos_list_allocate(3), pos_list_add_position(3), pos_list_from_csv(3),
 *      pos_list_sort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    pos_list_free(pos_list_t *pos_list)

{
    if ( pos_list->positions == NULL )
    {
	fputs("pos_list_free(): List pointer is NULL.\n", stderr);
	fputs("Was it previously allocated?\n", stderr);
	exit(EX_SOFTWARE);
    }
    pos_list->count = 0;
    pos_list->array_size = 0;
    free(pos_list->positions);
    pos_list->positions = NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc
 *
 *  Description:
 *      Add another position to a pos_list, expanding the array if needed.
 *
 *  Arguments:
 *      pos_list:   Pointer to the pos_list_t structure to add to
 *      position:   New position to be added to the list
 *
 *  See also:
 *      pos_list_allocate(3), pos_list_free(3), pos_list_from_csv(3),
 *      pos_list_sort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

int     pos_list_add_position(pos_list_t *pos_list, uint64_t position)

{
    if ( pos_list->count == pos_list->array_size )
    {
	pos_list->array_size *= 2;
	pos_list->positions = xt_realloc(pos_list->positions,
	    pos_list->array_size, sizeof(*pos_list->positions));
	if ( pos_list == NULL )
	{
	    fputs("pos_list_add_position(): Could not reallocate positions.\n", stderr);
	    exit(EX_UNAVAILABLE);
	}
    }
    pos_list->positions[pos_list->count++] = position;
    return BIO_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc
 *
 *  Description:
 *      Convert a comma-separated list of positions to a pos_list_t list.
 *      The array_size argument should be your best guess at the final size
 *      of the list.  Choosing a large enough value is not critical since
 *      it will be extended by pos_list_add_position() if needed.
 *
 *  Arguments:
 *      pos_list:   Pointer to the pos_list_t to receive the list
 *      list_str:   Character string containing comma-separated list of positions
 *      array_size: Initial array size
 *
 *  Returns:
 *      The number of positions added to the list
 *
 *  See also:
 *      pos_list_allocate(3), pos_list_add_position(3), pos_list_free(3),
 *      pos_list_sort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

int     pos_list_from_csv(pos_list_t *pos_list, const char *bounds_str,
		       size_t array_size)

{
    char        *copy, *p, *token, *end;
    size_t      c;
    uint64_t    position;
    
    if ( (copy = strdup(bounds_str)) == NULL )
    {
	fputs("peak-classifier: Cannot allocate temporary bounds string.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    pos_list_allocate(pos_list, array_size);
    for (p = copy, c = 0; (c < POS_LIST_ARRAY_SIZE(pos_list)) &&
			  ((token = strsep(&p, ",")) != NULL); ++c)
    {
	position = strtoull(token, &end, 10);
	if ( *end != '\0' )
	    return BIO_DATA_INVALID;
	else
	    pos_list_add_position(pos_list, position);
    }
    return c;
}


/***************************************************************************
 *  Compare two uint64_t values for sort functions.  Difference may
 *  exceed the range of an int, so don't just subtract.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

int     position_cmp_ascending(const uint64_t *pos1, const uint64_t *pos2)

{
    if ( *pos1 == *pos2 )
	return 0;
    else if ( *pos1 > *pos2 )
	return 1;
    else
	return -1;
}


int     position_cmp_descending(const uint64_t *pos1, const uint64_t *pos2)

{
    return -position_cmp_ascending(pos1, pos2);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc
 *
 *  Description:
 *      Sort a position list in either ascending or descending order.
 *
 *  Arguments:
 *      pos_list:   Pointer to the position_list_t structure to sort
 *      order:      POS_LIST_ASCENDING or POS_LIST_DESCENDING
 *
 *  See also:
 *      pos_list_allocate(3), pos_list_add_position(3), pos_list_free(3),
 *      pos_list_from_csv(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    pos_list_sort(pos_list_t *pos_list, pos_list_sort_order_t order)

{
    if ( order == POS_LIST_ASCENDING )
	qsort(pos_list->positions, pos_list->count, sizeof(pos_list->positions[0]),
	     (int (*)(const void *,const void *))position_cmp_ascending);
    else
	qsort(pos_list->positions, pos_list->count, sizeof(pos_list->positions[0]),
	     (int (*)(const void *,const void *))position_cmp_descending);
}
