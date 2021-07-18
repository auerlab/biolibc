#include <stdio.h>
#include <sysexits.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include "biolibc.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/biolibc.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Perform a numeric comparison of two chromosome names.
 *
 *      The names may contain a prefix of non-digits, such as "chr".
 *      Characters that follow must be a chromosome number or letter.
 *      Numbers are considered less than letters (e.g. 22 < X).  As such,
 *      if either is a letter, they are compared lexically.  If both are
 *      numbers, they are converted to integers and compared numerically.
 *
 *      Use bl_chromosome_name_cmp() only if you need to know which string is
 *      < or >.  If only checking for equality/inequality, strcmp() will be
 *      faster.
 *
 *  Arguments:
 *      name1, name2:   Names of two chromosomes
 *
 *  Returns:
 *      A value < 1 if name1 is numerically < name2
 *      A value > 1 if name1 is numerically > name2
 *      0 if name1 == name2
 *
 *  See also:
 *      strcmp(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-07  Jason Bacon Begin
 ***************************************************************************/

int     bl_chromosome_name_cmp(const char *name1, const char *name2)

{
    const char      *p1 = name1, *p2 = name2;
    char            *end;
    unsigned long   c1, c2;

    /* Skip identical portions of strings, e.g. "chr" prefix */
    while ( (*p1 == *p2) && (*p1 != '\0') )
	++p1, ++p2;
    
    /*
     *  Next should be a number or letter ID such as X or Y.  If either
     *  ID is not a number, simply compare lexically.
     *  ISO character order will take care of it since letters come after
     *  digits (chrX > chr22) and everything comes after null
     *  (chr22 > chr2).  This also handles the case where the names are
     *  the same (both *p1 and *p2 are '\0') or we reached the end of one
     *  of them (chr2, chr22).
     */
    if ( !isdigit(*p1) || !isdigit(*p2) )
	return *p1 - *p2;

    /* Both IDs are numeric, so perform an integer compare */
    c1 = strtoul(p1, &end, 10);
    if ( *end != '\0' )
    {
	fprintf(stderr,
		"bl_chromosome_name_cmp(): Invalid chromosome ID: %s\n", name1);
	exit(EX_DATAERR);
    }
    c2 = strtoul(p2, &end, 10);
    if ( *end != '\0' )
    {
	fprintf(stderr,
		"bl_chromosome_name_cmp(): Invalid chromosome ID: %s\n", name2);
	exit(EX_DATAERR);
    }
    return c1 - c2;
}
