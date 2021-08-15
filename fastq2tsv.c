/***************************************************************************
 *  Description:
 *      Remove duplicates from a fastq stream
 *  
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <string.h>
#include <xtend/string.h>
#include <biolibc/fastq.h>
#include <biolibc/biolibc.h>

int     main(int argc, char *argv[])

{
    bl_fastq_t  rec = BL_FASTQ_INIT;
    unsigned long   records = 0;
    
    while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
    {
	// Replace TABs in description to avoid interpretation as separators
	strtr(BL_FASTQ_DESC(&rec), "\t", " ", 0);
	
	printf("%s\t%s\t%s\t%s\n", BL_FASTQ_DESC(&rec), BL_FASTQ_SEQ(&rec),
		"+", BL_FASTQ_QUAL(&rec));
	
	++records;
    }
    bl_fastq_free(&rec);
    fprintf(stderr, "%lu records processed.\n", records);
    return EX_OK;
}
