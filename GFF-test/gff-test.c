/***************************************************************************
 *  Description:
 *      Test gff.c functions
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-19  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include "gff.h"

int     main(int argc,char *argv[])

{
    bl_gff_t    gff_feature;
    FILE        *tmpfile;
    int         ch;
    
    tmpfile = bl_gff_skip_header(stdin);
    while ( (ch = getc(tmpfile)) != EOF )
	putchar(ch);
    while ( bl_gff_read(&gff_feature, stdin, BL_GFF_FIELD_ALL) != EOF )
    {
	fprintf(stderr, "%s %s\n", BL_GFF_FEATURE_ID(&gff_feature),
		BL_GFF_GENE_NAME(&gff_feature));
	bl_gff_write(&gff_feature, stdout, BL_GFF_FIELD_ALL);
    }
    return EX_OK;
}
