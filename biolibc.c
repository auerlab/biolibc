#include <string.h>
#include <sys/stat.h>
#include "biolibc.h"

/***************************************************************************
 *  Description:
 *      Set overlap parameters
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

void    bio_set_overlap(bio_overlap_t *overlap,
			uint64_t f1_len, uint64_t f2_len,
			uint64_t ov_start, uint64_t ov_end)

{
    overlap->f1_len = f1_len;
    overlap->f2_len = f2_len;
    overlap->ov_start = ov_start;
    overlap->ov_end = ov_end;
    overlap->ov_len = ov_end - ov_start + 1;
}


/***************************************************************************
 *  Description:
 *      Print overlap parameters
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

void    bio_print_overlap(bio_overlap_t *overlap, char *f1_name, char *f2_name)

{
    char    f1_len[16], f2_len[16];
    
    strlcpy(f1_len, f1_name, 12);
    strlcat(f1_len, " len", 16);
    strlcpy(f2_len, f2_name, 12);
    strlcat(f2_len, " len", 16);
    printf("%-16s: %" PRIu64 "\n"
	   "%-16s: %" PRIu64 "\n"
	   "Overlap start   : %" PRIu64 "\n"
	   "Overlap end     : %" PRIu64 "\n"
	   "Overlap length  : %" PRIu64 "\n",
	   f1_len, overlap->f1_len,
	   f2_len, overlap->f2_len,
	   overlap->ov_start, overlap->ov_end, overlap->ov_len);
}


/***************************************************************************
 *  Description:
 *      Open a raw data file using foipen() or a gzipped, bzipped, or
 *      xzipped file using popen().  Must be used in conjunction with
 *      bio_fclose() to ensure that fclose() or pclose() is called where
 *      appropriate.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

FILE    *bio_fopen(const char *filename, const char *mode)

{
    char    *ext = strrchr(filename, '.'),
	    cmd[BIO_CMD_MAX + 1];
    
    if ( (strcmp(mode, "r") != 0 ) && (strcmp(mode, "w") != 0) )
    {
	fprintf(stderr, "bio_open(): Only \"r\" and \"w\" modes supported.\n");
	return NULL;
    }
    
    if ( ext == NULL )
    {
	fprintf(stderr, "bio_open(): No filename extension on %s.\n", filename);
	return NULL;
    }

    if ( *mode == 'r' )
    {
	if ( strcmp(ext, ".gz") == 0 )
	{
	    snprintf(cmd, BIO_CMD_MAX, "zcat %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, BIO_CMD_MAX, "bzcat %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, BIO_CMD_MAX, "xzcat %s", filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
    else    // "w"
    {
	if ( strcmp(ext, ".gz") == 0 )
	{
	    snprintf(cmd, BIO_CMD_MAX, "gzip -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, BIO_CMD_MAX, "bzip2 -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, BIO_CMD_MAX, "xz -c > %s", filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
}


/***************************************************************************
 *  Description:
 *      Close a FILE stream with fclose() or pclose() as appropriate
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-10  Jason Bacon Begin
 ***************************************************************************/

int     bio_fclose(FILE *stream)

{
    struct stat stat;
    
    fstat(fileno(stream), &stat);
    if ( S_ISFIFO(stat.st_mode) )
	return pclose(stream);
    else
	return fclose(stream);
}

