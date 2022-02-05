/***************************************************************************
 *  Description:
 *      Test basic fasta functionality
 *  
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <biolibc/fasta.h>
#include <biolibc/biolibc.h>

int     main(int argc,char *argv[])

{
    bl_fasta_t  rec = BL_FASTA_INIT;
    
    while ( bl_fasta_read(&rec, stdin) != BL_READ_EOF )
	bl_fasta_write(&rec, stdout, 100);
    bl_fasta_free(&rec);
    return EX_OK;
}

