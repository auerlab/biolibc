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
    
    while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
    {
	bl_fasta_write(stdout, &rec, 100);
	bl_fasta_free(&rec);
    }
    return EX_OK;
}

