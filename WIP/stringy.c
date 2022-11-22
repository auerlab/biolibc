#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
    int     low;
    int     high;
}   range_t;

int     main(int argc,char *argv[])

{
    range_t ranges[] = { {5,10}, {14,18}, {22,27} };
    char    string[] = "1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19-20-21-22-23-24-25-26-27-28-29-30",
	    expect[] = "A-A-A-A-X-X-X-X-X-X-B-B-B-X-X-X-X-X-A-A-A-X-X-X-X-X-X-B-B-B",
	    *p, *token,
	    alt[] = "AB";
    int     alt_index = 0, n, c, ch,
	    range_count = sizeof(ranges) / sizeof(range_t);
    
    for (p = string; (token = strsep(&p, "-")) != NULL;)
    {
	n = atoi(token);
	for (c = 0; (c < range_count); ++c)
	{
	    if ( (n >= ranges[c].low) && (n <= ranges[c].high) )
	    {
		ch = 'X';
		break;
	    }
	}
	if ( c == range_count )
	{
	    if ( ch == 'X' )
		alt_index = (alt_index + 1) % 2;
	    ch = alt[alt_index];
	}
	if ( token != string ) putchar('-');
	printf("%c", ch);
    }
    printf("\nExpected:\n%s\n", expect);
    return EX_OK;
}
