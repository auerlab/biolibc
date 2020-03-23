#include <stdio.h>
#include <sysexits.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>

/*
 *  Perform a numeric comparison of two chromosome names.
 *  The names must contain the chromosome number in the first digits present.
 *  Use this only if you need to know which string is < or >.
 *  If only checking for equality/inequality, use strcmp().
 */

int     chromosome_name_cmp(const char *n1, const char *n2)

{
    const char      *p1, *p2;
    char            *end;
    unsigned long   c1, c2;
    
    for (p1 = n1; !isdigit(*p1) && (*p1 != '\0'); ++p1)
	;
    for (p2 = n2; !isdigit(*p2) && (*p2 != '\0'); ++p2)
	;
    
    if ( (*p1 == '\0') || (*p2 == '\0') )
    {
	fprintf(stderr, "Invalid argument: chromosome_name_cmp(%s, %s).\n", n1, n2);
	exit(EX_DATAERR);
    }
    
    c1 = strtoul(p1, &end, 10);
    c2 = strtoul(p2, &end, 10);
    return c1 - c2;
}
