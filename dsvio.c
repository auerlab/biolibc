#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include "dsvio.h"

/***************************************************************************
 *  Description:
 *      Read next multiple-delim-separated field
 *      Fields may be separated by any character in the string delim
 *      Return delimiter ending the field (member of delim or newline)
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
 *  2021-02-24  Jason Bacon Begin
 ***************************************************************************/

int     dsv_read_field(FILE *infile, char buff[], size_t buff_size,
		       const char *delims, size_t *len)

{
    size_t  c;
    char    *p;
    int     ch;
    
    for (c = 0, p = buff; (c < buff_size) && 
			  ( strchr(delims, ch = getc(infile)) == NULL) &&
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
 *      Discard next field separated by any character in the string delims
 *      Return the delimiter encountered, possibly newline or EOF
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
 *  2021-02-24  Jason Bacon Begin
 ***************************************************************************/

int     dsv_skip_field(FILE *infile, const char *delims)

{
    int     ch;
    
    while ( (strchr(delims, ch = getc(infile)) == NULL) &&
	    (ch != '\n') && (ch != EOF) )
	;
    
    return ch;
}


/***************************************************************************
 *  Description:
 *      Discard the rest of the current input line.
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
 *      Read a line of an arbitrary DSV file
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
 *  2021-04-30  Jason Bacon Begin
 ***************************************************************************/

int     dsv_read_line(FILE *infile, dsv_line_t *dsv_line, const char *delims)

{
    int     actual_delim;
    char    field[DSV_FIELD_MAX_CHARS + 1];
    size_t  actual_len;
    
    dsv_line->array_size = 32;
    dsv_line->num_fields = 0;
    
    if ( (dsv_line->fields = malloc(dsv_line->array_size * sizeof(char *))) == NULL )
    {
	fputs("dsv_read_line(): Cannot allocate fields.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    
    if ( (dsv_line->delims = malloc(dsv_line->array_size * sizeof(char))) == NULL )
    {
	fputs("dsv_read_line(): Cannot allocate delims.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    
    // FIXME: Check actual_delim and/or actual_len to detect truncation
    while ( ((actual_delim = dsv_read_field(infile,
		field, DSV_FIELD_MAX_CHARS, delims, &actual_len)) != EOF) )
    {
	if ( (dsv_line->fields[dsv_line->num_fields] = strdup(field)) == NULL )
	{
	    fprintf(stderr, "dsv_read_line(): Cannot strdup() field %zu.\n",
		    dsv_line->num_fields - 1);
	    exit(EX_UNAVAILABLE);
	}
	dsv_line->delims[dsv_line->num_fields++] = actual_delim;
	if ( dsv_line->num_fields == dsv_line->array_size )
	{
	    dsv_line->array_size *= 2;
	    if ( (dsv_line->fields = realloc(dsv_line->fields,
		    dsv_line->array_size * sizeof(char *))) == NULL )
	    {
		fputs("dsv_read_line(): Cannot reallocate fields.\n", stderr);
		exit(EX_UNAVAILABLE);
	    }
	    
	    if ( (dsv_line->delims = realloc(dsv_line->delims,
		    dsv_line->array_size * sizeof(char))) == NULL )
	    {
		fputs("dsv_read_line(): Cannot reallocate delims.\n", stderr);
		exit(EX_UNAVAILABLE);
	    }
	}
	if ( actual_delim == '\n' )
	    break;
    }
    return actual_delim;
}


/***************************************************************************
 *  Description:
 *      Print an arbitrary DSV line
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
 *  2021-05-01  Jason Bacon Begin
 ***************************************************************************/

void    dsv_write_line(FILE *outfile, dsv_line_t *dsv_line)

{
    int     c;
    
    for (c = 0; c < dsv_line->num_fields; ++c)
	fprintf(outfile, "%s%c", dsv_line->fields[c], dsv_line->delims[c]);
}


/***************************************************************************
 *  Description:
 *      Duplicate an arbitrary DSV line
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
 *  2021-05-01  Jason Bacon Begin
 ***************************************************************************/

void    dsv_copy_line(dsv_line_t *dest, dsv_line_t *src)

{
    size_t  c;
    
    // Prune unused pointers in src
    dest->array_size = dest->num_fields = src->num_fields;
    
    dest->fields = malloc(dest->array_size * sizeof(char *));
    dest->delims = malloc(dest->array_size * sizeof(char));
    
    for (c = 0; c < src->num_fields; ++c)
    {
	dest->fields[c] = strdup(src->fields[c]);
	dest->delims[c] = src->delims[c];
    }
}


/***************************************************************************
 *  Description:
 *      Free allocated memory for a DSV object
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
 *  2021-05-01  Jason Bacon Begin
 ***************************************************************************/

void    dsv_free_line(dsv_line_t *dsv_line)

{
    int     c;
    
    if ( dsv_line->fields != NULL )
    {
	for (c = 0; c < dsv_line->num_fields; ++c)
	    free(dsv_line->fields[c]);
	free(dsv_line->fields);
    }
    dsv_line->num_fields = 0;
}


int     tsv_read_field(FILE *infile, char buff[], size_t buff_size,
		       size_t *len)

{
    return dsv_read_field(infile, buff, buff_size, "\t", len);
}


int     tsv_skip_field(FILE *infile)

{
    return dsv_skip_field(infile, "\t");
}


int     tsv_skip_rest_of_line(FILE *infile)

{
    return dsv_skip_rest_of_line(infile);
}


int     csv_read_field(FILE *infile, char buff[], size_t buff_size,
		       size_t *len)

{
    return dsv_read_field(infile, buff, buff_size, ",", len);
}


int     csv_skip_field(FILE *infile)

{
    return dsv_skip_field(infile, ",");
}


int     csv_skip_rest_of_line(FILE *infile)

{
    return dsv_skip_rest_of_line(infile);
}
