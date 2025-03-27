/***************************************************************************
 *  Description:
 *      Test basic fastq functionality
 *  
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>

#include <biolibc/fastq.h>
#include <biolibc/biolibc.h>

int     main(int argc,char *argv[])

{
    bl_fastq_t  rec = BL_FASTQ_INIT;
    int         min_qual;
    
    if ( argc != 2 )
    {
	fprintf(stderr, "Usage: %s min-qual\n", argv[0]);
	return 1;
    }
    min_qual = atoi(argv[1]);
    
    // Copy with trimming: Run multiple times with different min qual
    while ( bl_fastq_read(&rec, stdin) != BL_READ_EOF )
    {
	if ( min_qual >0 )
	    puts("Raw read:");
	bl_fastq_write(&rec, stdout, 100);
	
	ssize_t  pos = bl_fastq_find_5p_low_qual(&rec, min_qual, 33);
	if ( min_qual >0 )
	    printf("5' trim pos = %zd\n", pos);
	bl_fastq_5p_trim(&rec, pos);
	if ( min_qual >0 )
	{
	    puts("5' trimmed:");
	    bl_fastq_write(&rec, stdout, 100);
	}
	
	pos = bl_fastq_find_3p_low_qual(&rec, min_qual, 33);
	if ( min_qual >0 )
	    printf("3' trim pos = %zd\n", pos);
	bl_fastq_3p_trim(&rec, pos);
	if ( min_qual >0 )
	{
	    puts("3' trimmed:");
	    bl_fastq_write(&rec, stdout, 100);
	}
    }
    
    // bl_fastq_name_cmp()

    // bl_fastq_free(&rec);
    
    return EX_OK;
}

