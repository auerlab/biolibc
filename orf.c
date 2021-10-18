#include <stdio.h>
#include <sysexits.h>
#include <ctype.h>
#include "translate.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/translate.h>
 *      -lbiolibc
 *
 *  Description:
 *      Locate the next start codon in stream and report its position in
 *      the input.  Reported positions are 1-based.
 *  
 *  Arguments:
 *      rna_stream  FILE stream containing RNA sequence data
 *
 *  Returns:
 *      1-based position of the next start codon within the stream
 *      or EOF if no codon is found.
 *
 *  Examples:
 *      unsigned long   pos;
 *      char            codon[4];
 *      
 *      if ( (pos = next_start_codon(stdin)) != EOF )
 *      {
 *          printf("Start codon at %lu.\n", pos);
 *          if ( (pos = next_stop_codon(stdin, codon)) != EOF )
 *              printf("Stop codon %s at %lu.\n", codon, pos);
 *          else
 *              puts("No stop codon found.");
 *      }
 *      else
 *          puts("No start codon found.");
 *
 *  See also:
 *      next_stop_codon(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-17  Jason Bacon Begin
 ***************************************************************************/

unsigned long   next_start_codon(FILE *rna_stream)

{
    int     ch;
    
    /*
     *  Not actually scanning to EOF.  Returns from inside loop if
     *  start codon is found.
     *  Start codon is AUG.
     */
    
    while ( !feof(rna_stream) )
    {
	while ( ((ch = getc(rna_stream)) != EOF) && (tolower(ch) != 'a') )
	    ;
	if ( ch == EOF )
	    return EOF;
	if ( tolower(ch = getc(rna_stream)) == 'u' )
	{
	    if ( tolower(ch = getc(rna_stream)) == 'g' )
		return ftell(rna_stream) - 2;
	    else if ( ch == EOF )
		return EOF;
	    else
		ungetc(ch, rna_stream);
	}
	else if ( ch == EOF )
	    return EOF;
	else
	    ungetc(ch, rna_stream);
    }
    return EOF;
}

/***************************************************************************
 *  Library:
 *      #include <biolibc/translate.h>
 *      -lbiolibc
 *
 *  Description:
 *      Locate the next stop codon in stream and report its position in
 *      the input.  Reported positions are 1-based.
 *  
 *  Arguments:
 *      rna_stream  FILE stream containing RNA sequence data
 *
 *  Returns:
 *      1-based position of the next stop codon within the stream
 *      or EOF if no codon is found.
 *
 *  Examples:
 *      unsigned long   pos;
 *      char            codon[4];
 *      
 *      if ( (pos = next_start_codon(stdin)) != EOF )
 *      {
 *          printf("Start codon at %lu.\n", pos);
 *          if ( (pos = next_stop_codon(stdin, codon)) != EOF )
 *              printf("Stop codon %s at %lu.\n", codon, pos);
 *          else
 *              puts("No stop codon found.");
 *      }
 *      else
 *          puts("No start codon found.");
 *
 *  See also:
 *      next_stop_codon(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-17  Jason Bacon Begin
 ***************************************************************************/

unsigned long   next_stop_codon(FILE *rna_stream, char codon[4])

{
    int     ch1, ch2, ch3;
    
    /*
     *  Not actually scanning to EOF.  Returns from inside loop if
     *  stop codon is found.
     *  Valid codons are UAG, UAA, UGA.
     */
    
    codon[0] = codon[3] = '\0';
    
    while ( !feof(rna_stream) )
    {
	while ( ((ch1 = tolower(getc(rna_stream))) != EOF) && (ch1 != 'u') )
	    ;
	if ( ch1 == EOF )
	    return EOF;
	if ( (ch2 = tolower(getc(rna_stream))) == 'a' )
	{
	    if ( (ch3 = tolower(getc(rna_stream))) == 'g' || (ch3 == 'a') )
	    {
		codon[0] = ch1; codon[1] = ch2; codon[2] = ch3;
		return ftell(rna_stream) - 2;
	    }
	    else if ( ch3 == EOF )
		return EOF;
	    else
		ungetc(ch3, rna_stream);
	}
	else if ( ch2 == 'g' )
	{
	    if ( (ch3 = tolower(getc(rna_stream))) == 'a' )
	    {
		codon[0] = ch1; codon[1] = ch2; codon[2] = ch3;
		return ftell(rna_stream) - 2;
	    }
	    else if ( ch3 == EOF )
		return EOF;
	    else
		ungetc(ch3, rna_stream);
	}
	else if ( ch2 == EOF )
	    return EOF;
	else
	    ungetc(ch2, rna_stream);
    }
    return EOF;
}


int     main(int argc,char *argv[])

{
    unsigned long   pos;
    char            codon[4];
    
    if ( (pos = next_start_codon(stdin)) != EOF )
    {
	printf("Start codon at %lu.\n", pos);
	if ( (pos = next_stop_codon(stdin, codon)) != EOF )
	    printf("Stop codon %s at %lu.\n", codon, pos);
	else
	    puts("No stop codon found.");
    }
    else
	puts("No start codon found.");
    return EX_OK;
}
