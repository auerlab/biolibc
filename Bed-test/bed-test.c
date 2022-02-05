/***************************************************************************
 *  Description:
 *      Test bedio.c functions
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-19  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include "bed.h"

int     main(int argc,char *argv[])

{
    bl_bed_t   bed_feature;
    
    bl_bed_skip_header(stdin);
    while ( bl_bed_read(&bed_feature, stdin, BL_BED_FIELD_ALL) != EOF )
	bl_bed_write(&bed_feature, stdout, BL_BED_FIELD_ALL);
    return EX_OK;
}

