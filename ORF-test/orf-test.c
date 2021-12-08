/***************************************************************************
 *  Description:
 *      Find ORF within a stream of RNA nulceotides
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-14  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <ctype.h>
#include <biolibc/translate.h>

int     main(int argc,char *argv[])

{
    unsigned long   pos;
    char            codon[4];
    
    if ( (pos = bl_next_start_codon(stdin, codon)) != EOF )
    {
	printf("Start codon %s at %lu.\n", codon, pos);
	if ( (pos = bl_next_stop_codon(stdin, codon)) != EOF )
	    printf("Stop codon %s at %lu.\n", codon, pos);
	else
	    puts("No stop codon found.");
    }
    else
	puts("No start codon found.");
    return EX_OK;
}
