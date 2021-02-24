#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include "dsvio.h"

/***************************************************************************
 *  Description:
 *      Read next delim-separated field
 *      Return delimiter ending the field, possibly newline or EOF
 *      This function is separate from mdsv_read_field() for performance
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     dsv_read_field(FILE *infile, char buff[], size_t buff_size,
		       int delim, size_t *len)

{
    size_t  c;
    char    *p;
    int     ch;
    
    for (c = 0, p = buff; (c < buff_size) && ((ch = getc(infile)) != delim) &&
			  (ch != '\n') && (ch != EOF); ++c, ++p )
	*p = ch;
    *p = '\0';
    
    if ( c == buff_size )
    {
	fprintf(stderr, "dsv_read_field(): Buffer overflow reading field.\n");
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
 *      Discard next field separated by the character delim
 *      Return the delimiter encountered, possibly newline or EOF
 *      This function is separate from mdsv_skip_field() for performance
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     dsv_skip_field(FILE *infile, int delim)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != delim) && (ch != '\n') && (ch != EOF) )
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

int     dsv_skip_rest_of_line(FILE *infile)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != EOF) && (ch != '\n') )
	;
    return ch;
}


/***************************************************************************
 *  Description:
 *      Read next multiple-delim-separated field
 *      Fields may be separated by any character in the string delim
 *      Return delimiter ending the field (member of delim or newline)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-24  Jason Bacon Begin
 ***************************************************************************/

int     mdsv_read_field(FILE *infile, char buff[], size_t buff_size,
		       char *delim, size_t *len)

{
    size_t  c;
    char    *p;
    int     ch;
    
    for (c = 0, p = buff; (c < buff_size) && 
			  ( strchr(delim, ch = getc(infile)) == NULL) &&
			  (ch != '\n') && (ch != EOF); ++c, ++p )
	*p = ch;
    *p = '\0';
    
    if ( c == buff_size )
    {
	fprintf(stderr, "dsv_read_field(): Buffer overflow reading field.\n");
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
 *      Discard next field separated by any character in the string delim
 *      Return the delimiter encountered, possibly newline or EOF
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-24  Jason Bacon Begin
 ***************************************************************************/

int     mdsv_skip_field(FILE *infile, char *delim)

{
    int     ch;
    
    while ( (strchr(delim, ch = getc(infile)) == NULL) &&
	    (ch != '\n') && (ch != EOF) )
	;
    
    return ch;
}


int     tsv_read_field(FILE *infile, char buff[], size_t buff_size,
		       size_t *len)

{
    return dsv_read_field(infile, buff, buff_size, '\t', len);
}


int     tsv_skip_field(FILE *infile)

{
    return dsv_skip_field(infile, '\t');
}


int     tsv_skip_rest_of_line(FILE *infile)

{
    return dsv_skip_rest_of_line(infile);
}


int     csv_read_field(FILE *infile, char buff[], size_t buff_size,
		       size_t *len)

{
    return dsv_read_field(infile, buff, buff_size, ',', len);
}


int     csv_skip_field(FILE *infile)

{
    return dsv_skip_field(infile, ',');
}


int     csv_skip_rest_of_line(FILE *infile)

{
    return dsv_skip_rest_of_line(infile);
}
