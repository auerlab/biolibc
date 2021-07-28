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
#include <biolibc/fastq.h>
#include <biolibc/biolibc.h>

int     main(int argc,char *argv[])

{
    bl_fastq_t  rec = BL_FASTQ_INIT;
    
    while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
	bl_fastq_write(stdout, &rec, 100);
    bl_fastq_free(&rec);
    return EX_OK;
}

