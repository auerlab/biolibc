#include <stdio.h>
#include <ctype.h>
#include "translate.h"

/***************************************************************************
 *  Name:
 *      bl_next_start_codon() - Find next start codon
 *
 *  Library:
 *      #include <biolibc/translate.h>
 *      -lbiolibc
 *
 *  Description:
 *      Locate the next start codon in stream and report its position in
 *      the input.  Reported positions are 0-based offsets from the file
 *      position at the time of the call.  Hence, to find the absolute
 *      positions of codons within a file stream across multiple calls to
 *      next_start_codon() or next_stop_codon(), their return values should
 *      added to the previous return value + 3 (the size of the previous
 *      codon).  For example, givent the following input:
 *
 *      acaucauguucguggugacc
 *
 *      A call to next_start_codon() will return 5 (off set of AUG from the
 *      beginning of the sequence and a subsequent call to next_stop_codon()
 *      will return 7 (offset of UGA from the first base after the AUG).
 *  
 *  Arguments:
 *      rna_stream  FILE stream containing RNA sequence data
 *      codon       4-character buffer to receive codon sequence
 *                  Set to "" if no codon found
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

long    bl_next_start_codon(FILE *rna_stream, char codon[4])

{
    int     ch1,
	    ch2,
	    ch3;
    long    pos = 0;
    
    // Blank and null-terminate
    codon[0] = codon[3] = '\0';
    
    /*
     *  Not actually scanning to EOF.  Returns from inside loop if
     *  start codon is found.
     *  Start codon is AUG.
     */
    
    while ( !feof(rna_stream) )
    {
	while ( ((ch1 = toupper(getc(rna_stream))) != EOF) && (ch1 != 'A') )
	    ++pos;
	if ( ch1 != EOF )
	{
	    ++pos;  // Count the 'A'
	    if ( ((ch2 = toupper(getc(rna_stream))) == 'U') || (ch2 == 'T') )
	    {
		if ( (ch3 = toupper(getc(rna_stream))) == 'G' )
		{
		    codon[0] = ch1; codon[1] = ch2; codon[2] = ch3;
		    return pos - 1;
		}
		else if ( ch3 != EOF )
		{
		    ungetc(ch3, rna_stream);
		    ungetc(ch2, rna_stream);
		}
	    }
	    else if ( ch2 != EOF )
		ungetc(ch2, rna_stream);
	}
    }
    return EOF;
}


/***************************************************************************
 *  Name:
 *      bl_next_stop_codon() - Find next stop codon
 *
 *  Library:
 *      #include <biolibc/translate.h>
 *      -lbiolibc
 *
 *  Description:
 *      Locate the next stop codon in stream and report its position in
 *      the input.  Reported positions are 0-based offsets from the file
 *      position at the time of the call.  Hence, to find the absolute
 *      positions of codons within a file stream across multiple calls to
 *      next_start_codon() or next_stop_codon(), their return values should
 *      added to the previous return value + 3 (the size of the previous
 *      codon).  For example, givent the following input:
 *
 *      acaucauguucguggugacc
 *
 *      A call to next_start_codon() will return 5 (off set of AUG from the
 *      beginning of the sequence and a subsequent call to next_stop_codon()
 *      will return 7 (offset of UGA from the first base after the AUG).
 *  
 *  Arguments:
 *      rna_stream  FILE stream containing RNA sequence data
 *      codon       4-character buffer to receive codon sequence
 *                  Set to "" if no codon found
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

long    bl_next_stop_codon(FILE *rna_stream, char codon[4])

{
    int     ch1,
	    ch2,
	    ch3;
    long    pos = 0;
    
    // Blank and null-terminate
    codon[0] = codon[3] = '\0';
    
    /*
     *  Not actually scanning to EOF.  Returns from inside loop if
     *  stop codon is found.
     *  Valid codons are UAG, UAA, UGA.
     */
    
    while ( !feof(rna_stream) )
    {
	while ( ((ch1 = toupper(getc(rna_stream))) != EOF) &&
		(ch1 != 'U') && (ch1 != 'T') )
	    ++pos;
	if ( ch1 != EOF )
	{
	    ++pos;  // Count the 'U' or 'T'
	    if ( (ch2 = toupper(getc(rna_stream))) == 'A' )
	    {
		if ( (ch3 = toupper(getc(rna_stream))) == 'G' || (ch3 == 'A') )
		{
		    codon[0] = ch1; codon[1] = ch2; codon[2] = ch3;
		    return pos - 1;
		}
		else if ( ch3 != EOF )
		{
		    ungetc(ch3, rna_stream);
		    ungetc(ch2, rna_stream);
		}
	    }
	    else if ( ch2 == 'G' )
	    {
		if ( (ch3 = toupper(getc(rna_stream))) == 'A' )
		{
		    codon[0] = ch1; codon[1] = ch2; codon[2] = ch3;
		    return pos - 1;
		}
		else if ( ch3 != EOF )
		{
		    ungetc(ch3, rna_stream);
		    ungetc(ch2, rna_stream);
		}
	    }
	    else if ( ch2 != EOF )
		ungetc(ch2, rna_stream);
	}
    }
    return EOF;
}
