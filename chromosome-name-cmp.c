#include <stdio.h>
#include <sysexits.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include "biostring.h"

/*
 *  Perform a numeric comparison of two chromosome names.
 *
 *  The names may contain a prefix of non-digits, such as "chr".
 *  Characters that follow must be a chromosome number or letter.
 *  Numbers come before letters (e.g. 22 before X).  As such, if either
 *  is a letter, they are compared lexically.  If both are numbers, they
 *  are converted to integers and compared numerically.
 *
 *  Use this only if you need to know which string is < or >.
 *  If only checking for equality/inequality, strcmp() will be faster.
 *
 *  FIXME: Do some sort of input validation and designate an error code
 *         such as MAXINT
 */

int     chromosome_name_cmp(const char *n1, const char *n2)

{
    const char      *p1 = n1, *p2 = n2;
    char            *end;
    unsigned long   c1, c2;
    
    /* Skip identical portions of strings */
    while ( (*p1 == *p2) && (*p1 != '\0') )
	++p1, ++p2;
    
    /* 
     *  Compare letters, letters to numbers, or letters to null lexically.
     *  ISO character order will take care of it since letters come after
     *  digits (chrX > chr22) and everything comes after null
     *  (chr22 > chr2).
     */
    if ( !isdigit(*p1) || !isdigit(*p2) )
	return *p1 - *p2;
    
    if ( (*p1 == '\0') || (*p2 == '\0') )
    {
	fprintf(stderr, "Invalid argument: chromosome_name_cmp(%s, %s).\n", n1, n2);
	exit(EX_DATAERR);
    }
    
    c1 = strtoul(p1, &end, 10);
    c2 = strtoul(p2, &end, 10);
    return c1 - c2;
}
