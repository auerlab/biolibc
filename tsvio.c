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
    
    if ( c == buff_size )
    {
	fprintf(stderr, "%s: tsv_read_field(): Buffer overflow reading field.\n", argv[0]);
	fputs(buff, stderr);
	exit(EX_DATAERR);
    }
    
    if ( ch == EOF )
    {
	fprintf(stderr, "%s: tsv_read_field(): Encountered EOF.  This should not happen.\n", argv[0]);
	fprintf(stderr, "Either there's a program bug or your VCF input is corrupt.\n");
    }
    
    *p = '\0';
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

void    tsv_skip_field(const char *argv[], FILE *infile)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != '\t') && (ch != '\n') && (ch != EOF) )
	;
    
    if ( ch == EOF )
    {
	fprintf(stderr, "%s: tsv_skip_field(): EOF encounterd unexpectedly.\n", argv[0]);
	exit(EX_DATAERR);
    }
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
