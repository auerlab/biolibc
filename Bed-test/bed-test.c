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
#include "bedio.h"

int     main(int argc,char *argv[])

{
    bed_feature_t   bed_feature;
    
    bed_skip_header(stdin);
    while ( bed_read_feature(stdin, &bed_feature) != EOF )
	bed_write_feature(stdout, &bed_feature, BED_FIELD_ALL);
    return EX_OK;
}

