#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include "tsvio.h"

/***************************************************************************
 *  Description:
 *      Read next tab-separated field
 *      Return delimiter ending the field (should be tab or newline)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     tsv_read_field(const char *argv[], FILE *infile,
		       char buff[], size_t buff_size, size_t *len)

{
    size_t  c;
    char    *p;
    int     ch;
    
    for (c = 0, p = buff; (c < buff_size) && ((ch = getc(infile)) != '\t') &&
			  (ch != '\n') && (ch != EOF); ++c, ++p )
	*p = ch;
    *p = '\0';
    
    if ( c == buff_size )
    {
	fprintf(stderr, "%s: tsv_read_field(): Buffer overflow reading field.\n", argv[0]);
	fprintf(stderr, "Buffer size = %zu\n", buff_size);
	fputs(buff, stderr);
	// FIXME: Replace this with another sentinal value?
	// Would require all callers to handle both EOF and overflow
	return EOF;
    }
    
    *len = c;
    return ch;
}


/***************************************************************************
 *  Description:
 *      Discard next field
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     tsv_skip_field(const char *argv[], FILE *infile)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != '\t') && (ch != '\n') && (ch != EOF) )
	;
    
    return ch;
}


/***************************************************************************
 *  Description:
 *      Discard the rest of the current input line.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     tsv_skip_rest_of_line(const char *argv[], FILE *infile)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != EOF) && (ch != '\n') )
	;
    return ch;
}
