#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <string.h>
#include "plist.h"
#include "biolibc.h"

/***************************************************************************
 *  Description:
 *      Initialize a position list
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
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    plist_allocate(plist_t *plist, size_t max_positions)

{
    if ( (plist->count != 0) || (plist->max_positions != 0) ||
	 (plist->positions != NULL) )
    {
	fputs("plist_allocate(): List is not blank.\n", stderr);
	fputs("Was it previously allocated?\n", stderr);
	fputs("Did you forget to initialize it with PLIST_INIT?\n", stderr);
	exit(EX_SOFTWARE);
    }
    plist->positions = malloc(max_positions * sizeof(uint64_t));
    if ( plist->positions == NULL )
    {
	fputs("plist_allocate(): Cannot allocate list.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    plist->max_positions = max_positions;
}


/***************************************************************************
 *  Description:
 *      Free a position list
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
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    plist_free(plist_t *plist)

{
    if ( plist->positions == NULL )
    {
	fputs("plist_free(): List pointer is NULL.\n", stderr);
	fputs("Was it previously allocated?\n", stderr);
	exit(EX_SOFTWARE);
    }
    plist->count = 0;
    plist->max_positions = 0;
    free(plist->positions);
}


/***************************************************************************
 *  Description:
 *      Add another position to a plist
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
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

int     plist_add_position(plist_t *plist, uint64_t position)

{
    if ( plist->count == plist->max_positions )
    {
	fputs("plist_add_position(): List is full.\n", stderr);
	exit(EX_SOFTWARE);
    }
    plist->positions[plist->count++] = position;
    return BIO_DATA_OK;
}


/***************************************************************************
 *  Description:
 *      Convert command-line syntax for boundaries to an array of positions.
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
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

int     plist_from_csv(plist_t *plist, const char *bounds_str,
		       size_t max_positions)

{
    char        *copy, *p, *token, *end;
    size_t      c;
    uint64_t    position;
    
    if ( (copy = strdup(bounds_str)) == NULL )
    {
	fputs("peak-classifier: Cannot allocate temporary bounds string.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    plist_allocate(plist, max_positions);
    for (p = copy, c = 0; (c < PLIST_MAX_POSITIONS(plist)) &&
			  ((token = strsep(&p, ",")) != NULL); ++c)
    {
	position = strtoull(token, &end, 10);
	if ( *end != '\0' )
	    return BIO_INVALID_DATA;
	else
	    plist_add_position(plist, position);
    }
    return c;
}


/***************************************************************************
 *  Description:
 *      Compare two uint64_t values for sort functions.  Difference may
 *      exceed the range of an int, so don't just subtract.
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
 *  Description:
 *      Sort a position list
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
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    plist_sort(plist_t *plist, plist_sort_order_t order)

{
    if ( order == PLIST_ASCENDING )
	qsort(plist->positions, plist->count, sizeof(plist->positions[0]),
	     (int (*)(const void *,const void *))position_cmp_ascending);
    else
	qsort(plist->positions, plist->count, sizeof(plist->positions[0]),
	     (int (*)(const void *,const void *))position_cmp_descending);
}

