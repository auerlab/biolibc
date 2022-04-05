#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include <xtend/string.h>      // strlcpy() on Linux
#include <xtend/dsv.h>
#include <xtend/mem.h>
#include <xtend/file.h>
#include "sam.h"
#include "biolibc.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read next alignment (line) from a SAM stream.
 *
 *      If field_mask is not BL_SAM_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in sam_alignment.
 *      That field in the structure is then populated with an appropriate
 *      marker, such as '.'.  Possible mask values are:
 *
 *      BL_SAM_FIELD_ALL
 *      BL_SAM_FIELD_QNAME
 *      BL_SAM_FIELD_FLAG
 *      BL_SAM_FIELD_RNAME
 *      BL_SAM_FIELD_POS
 *      BL_SAM_FIELD_MAPQ
 *      BL_SAM_FIELD_CIGAR
 *      BL_SAM_FIELD_RNEXT
 *      BL_SAM_FIELD_PNEXT
 *      BL_SAM_FIELD_TLEN
 *      BL_SAM_FIELD_SEQ
 *      BL_SAM_FIELD_QUAL
 *
 *  Arguments:
 *      sam_stream      A FILE stream from which to read the line
 *      sam_alignment   Pointer to a bl_sam_t structure
 *      field_mask      Bit mask indicating which fields to store in sam_alignment
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered after a complete feature
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      bl_sam_read(stdin, &sam_alignment, BL_SAM_FIELD_ALL);
 *      bl_sam_read(sam_stream, &sam_alignment,
 *                         BL_SAM_FIELD_QNAME|BL_SAM_FIELD_POS|BL_SAM_FIELD_TLEN);
 *
 *  See also:
 *      bl_sam_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_read(bl_sam_t *sam_alignment, FILE *sam_stream,
			   sam_field_mask_t field_mask)

{
    char    mapq_str[BL_SAM_MAPQ_MAX_CHARS + 1],
	    pos_str[BL_POSITION_MAX_DIGITS + 1],
	    flag_str[BL_SAM_FLAG_MAX_DIGITS + 1],
	    *end;
    size_t  len;
    static int64_t   previous_pos = 0;
    int     delim;
    
    if ( field_mask & BL_SAM_FIELD_QNAME )
	delim = tsv_read_field(sam_stream, sam_alignment->qname,
			BL_SAM_QNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	*sam_alignment->qname = '\0';
    }
    if ( delim == EOF )
	return BL_READ_EOF;

    // 2 Flag
    if ( field_mask & BL_SAM_FIELD_FLAG )
	delim = tsv_read_field(sam_stream, flag_str, BL_SAM_FLAG_MAX_DIGITS, &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading flag: %s.\n",
		flag_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_FLAG )
    {
	sam_alignment->flag = strtoul(flag_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid position: %s\n",
		    flag_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	sam_alignment->flag = 0;    // FIXME: Is there a better choice?
    
    // 3 RNAME
    if ( field_mask & BL_SAM_FIELD_RNAME )
	delim = tsv_read_field(sam_stream, sam_alignment->rname,
			       BL_SAM_RNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	*sam_alignment->rname = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading rname: %s.\n",
		sam_alignment->rname);
	return BL_READ_TRUNCATED;
    }
    
    // 4 POS
    if ( field_mask & BL_SAM_FIELD_POS )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading pos: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_POS )
    {
	sam_alignment->pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid position: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
	previous_pos = sam_alignment->pos;
    }
    else
	sam_alignment->pos = 0;
    
    // 5 MAPQ
    if ( field_mask & BL_SAM_FIELD_MAPQ )
	delim = tsv_read_field(sam_stream, mapq_str, BL_SAM_MAPQ_MAX_CHARS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading mapq: %s.\n",
		mapq_str);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_MAPQ )
    {
	sam_alignment->mapq = strtoul(mapq_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid mapq: %s\n",
		    mapq_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	sam_alignment->mapq = 0;
    
    // 6 CIGAR
    if ( field_mask & BL_SAM_FIELD_CIGAR )
	delim = tsv_read_field(sam_stream, sam_alignment->cigar,
			       BL_SAM_CIGAR_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	*sam_alignment->cigar = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading cigar: %s.\n",
		sam_alignment->cigar);
	return BL_READ_TRUNCATED;
    }
    
    // 7 RNEXT
    if ( field_mask & BL_SAM_FIELD_RNEXT )
	delim = tsv_read_field(sam_stream, sam_alignment->rnext,
			       BL_SAM_RNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	*sam_alignment->rnext = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading rnext: %s.\n",
		sam_alignment->rnext);
	return BL_READ_TRUNCATED;
    }
    
    // 8 PNEXT
    if ( field_mask & BL_SAM_FIELD_PNEXT )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading pnext: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_PNEXT )
    {
	sam_alignment->pnext = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid pnext: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	sam_alignment->pnext = 0;
    
    // 9 TLEN
    if ( field_mask & BL_SAM_FIELD_TLEN )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading tlen: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_TLEN )
    {
	sam_alignment->tlen = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid tlen: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    sam_alignment->qname, sam_alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	sam_alignment->tlen = 0;
    
    // 10 SEQ
    if ( field_mask & BL_SAM_FIELD_SEQ )
	delim = tsv_read_field_malloc(sam_stream, &sam_alignment->seq,
		    &sam_alignment->seq_array_size, &sam_alignment->seq_len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	// sam_alignment->seq = NULL;   Mem leak if allocated elsewhere
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading seq: %s.\n",
		sam_alignment->seq);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_SEQ )
    {
	// May be allocated by bl_sam_init() or bl_sam_copy()
	if ( sam_alignment->seq == NULL )
	{
	    if ( (sam_alignment->seq = xt_malloc(sam_alignment->seq_len + 1,
		    sizeof(*sam_alignment->seq))) == NULL )
	    {
		fprintf(stderr, "bl_sam_read(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
    }
    
    // 11 QUAL, should be last field
    if ( field_mask & BL_SAM_FIELD_QUAL )
	delim = tsv_read_field_malloc(sam_stream, &sam_alignment->qual,
		    &sam_alignment->qual_array_size,
		    &sam_alignment->qual_len);
    else
    {
	delim = tsv_skip_field(sam_stream);
	// sam_alignment->qual = NULL; Mem leak if allocated elsewhere
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading qual: %s.\n",
		sam_alignment->qual);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_QUAL )
    {
	// May be allocated by bl_sam_init() or bl_sam_copy()
	if ( sam_alignment->qual == NULL )
	{
	    if ( (sam_alignment->qual = xt_malloc(sam_alignment->qual_len + 1,
		    sizeof(*sam_alignment->qual))) == NULL )
	    {
		fprintf(stderr, "bl_sam_read(): Could not allocate qual.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
    
	if ( (sam_alignment->qual_len != 1) &&
	     (sam_alignment->seq_len != sam_alignment->qual_len) )
	    fprintf(stderr, "bl_sam_read(): Warning: qual_len != seq_len for %s,%" PRId64 "\n",
		    sam_alignment->rname, sam_alignment->pos);
    }
    
    // Some SRA CRAMs have 11 fields, most have 12
    // Discard everything after the 11th
    if ( delim == '\t' )
	while ( getc(sam_stream) != '\n' )
	    ;

    /*fprintf(stderr,"bl_sam_read(): %s,%" PRId64 ",%zu\n",
	    BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
	    BL_SAM_SEQ_LEN(sam_alignment));*/
    
    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Copy a SAM alignment as efficiently as possible, allocating memory
 *      as needed.
 *
 *  Arguments:
 *      dest    Pointer to bl_sam_t structure to receive copy
 *      src     Pointer to bl_sam_t structure to be copied
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_init(3), bl_sam_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_copy(bl_sam_t *dest, bl_sam_t *src)

{
    strlcpy(dest->qname, src->qname, BL_SAM_QNAME_MAX_CHARS + 1);
    dest->flag = src->flag;
    strlcpy(dest->rname, src->rname, BL_SAM_RNAME_MAX_CHARS + 1);
    dest->pos = src->pos;
    dest->mapq = src->mapq;
    // FIXME: Add cigar and RNEXT
    strlcpy(dest->cigar, src->cigar, BL_SAM_CIGAR_MAX_CHARS + 1);
    strlcpy(dest->rnext, src->rnext, BL_SAM_RNAME_MAX_CHARS + 1);
    dest->pnext = src->pnext;
    dest->tlen = src->tlen;
    
    if ( (dest->seq = xt_malloc(src->seq_len + 1,
	    sizeof(*dest->seq))) == NULL )
    {
	fprintf(stderr, "bl_sam_copy(): Could not allocate seq.\n");
	exit(EX_UNAVAILABLE);
    }
    memcpy(dest->seq, src->seq, src->seq_len + 1);
    
    if ( (dest->qual = xt_malloc(src->seq_len + 1,
	    sizeof(*dest->qual))) == NULL )
    {
	fprintf(stderr, "bl_sam_copy(): Could not allocate qual.\n");
	exit(EX_UNAVAILABLE);
    }
    
    /* qual is an optional field */
    if ( src->qual_len > 0 )
	memcpy(dest->qual, src->qual, src->qual_len + 1);
    
    dest->seq_len = src->seq_len;
    dest->qual_len = src->qual_len;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_sam_read() or
 *      bl_sam_init().
 *
 *  Arguments:
 *      sam_alignment   Pointer to bl_sam_t structure to be freed.
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_init(3), bl_sam_copy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_free(bl_sam_t *sam_alignment)

{
    if ( sam_alignment->seq != NULL )
	free(sam_alignment->seq);
    if ( sam_alignment->qual != NULL )
	free(sam_alignment->qual);
    // FIXME: Cigar and rnext?
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_sam_t structure, allocating memory for
 *      sequence and quality strings according to seq_len.  Passing a
 *      seq_len of 0 prevents memory allocation from occurring.
 *
 *      Only BL_SAM_FIELD_SEQ and BL_SAM_FIELD_QUAL are meaningful bits in
 *      field_mask, as they determine whether memory is allocated.  All
 *      other fields are unconditionally initialized to 0, NULL, or blank.
 *
 *  Arguments:
 *      sam_alignment   Pointer to bl_sam_t structure to initialize
 *      seq_len         Length of sequence and quality strings
 *      field_mask      Bit mask indicating which fields will be used
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_free(3), bl_sam_copy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_init(bl_sam_t *sam_alignment)

{
    *sam_alignment->qname = '\0';
    sam_alignment->flag = 0;
    *sam_alignment->rname = '\0';
    sam_alignment->pos = 0;
    sam_alignment->mapq = 0;
    *sam_alignment->cigar = '\0';
    *sam_alignment->rnext = '\0';
    sam_alignment->pnext = 0;
    sam_alignment->tlen = 0;
    sam_alignment->seq = NULL;
    sam_alignment->qual = NULL;
    sam_alignment->seq_array_size = 0;
    sam_alignment->seq_len = 0;
    sam_alignment->qual_array_size = 0;
    sam_alignment->qual_len = 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write an alignment (line) to a SAM stream.
 *
 *      If field_mask is not BL_SAM_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate placeholder such as '.'
 *      rather than stored in sam_alignment.  Possible mask values are:
 *
 *      BL_SAM_FIELD_ALL
 *      BL_SAM_FIELD_QNAME
 *      BL_SAM_FIELD_FLAG
 *      BL_SAM_FIELD_RNAME
 *      BL_SAM_FIELD_POS
 *      BL_SAM_FIELD_MAPQ
 *      BL_SAM_FIELD_CIGAR
 *      BL_SAM_FIELD_RNEXT
 *      BL_SAM_FIELD_PNEXT
 *      BL_SAM_FIELD_TLEN
 *      BL_SAM_FIELD_SEQ
 *      BL_SAM_FIELD_QUAL
 *
 *  Arguments:
 *      sam_stream      A FILE stream to which to write the line
 *      sam_alignment   Pointer to a bl_sam_t structure
 *      field_mask      Bit mask indicating which fields to store in sam_alignment
 *
 *  Returns:
 *      Number of items written (per fprintf() output)
 *
 *  See also:
 *      bl_sam_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_write(bl_sam_t *sam_alignment, FILE *sam_stream,
			   sam_field_mask_t field_mask)

{
    int     count;
    
    count = fprintf(sam_stream, "%s\t%u\t%s\t%" PRId64
		    "\t%u\t%s\t%s\t%" PRId64 "zu\t%zu\t%s\t%s\t%zu\t%zu\n",
		    sam_alignment->qname,
		    sam_alignment->flag,
		    sam_alignment->rname,
		    sam_alignment->pos,
		    sam_alignment->mapq,
		    sam_alignment->cigar,
		    sam_alignment->rnext,
		    sam_alignment->pnext,
		    sam_alignment->tlen,
		    sam_alignment->seq,
		    sam_alignment->qual,
		    sam_alignment->seq_len,
		    sam_alignment->qual_len);
    return count;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      Open a raw SAM file using fopen() or a gzipped, bzipped, or
 *      xzipped SAM file or BAM or CRAM file using popen().
 *      Must be used in conjunction with
 *      bl_sam_fclose() to ensure that fclose() or pclose() is called where
 *      appropriate.
 *
 *  Arguments:
 *      filename:   Name of the file to be opened
 *      mode:       "r" or "w", passed to fopen() or popen()
 *
 *  Returns:
 *      A pointer to the FILE structure or NULL if open failed
 *
 *  See also:
 *      fopen(3), popen(3), gzip(1), bzip2(1), xz(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-05  Jason Bacon Derived from xt_fclose()
 ***************************************************************************/

FILE    *bl_sam_fopen(const char *filename, const char *mode)

{
    char    *ext = strrchr(filename, '.'),
	    cmd[XT_CMD_MAX_CHARS + 1];
    
    if ( (strcmp(mode, "r") != 0 ) && (strcmp(mode, "w") != 0) )
    {
	fprintf(stderr, "bl_sam_fopen(): Only \"r\" and \"w\" modes supported.\n");
	return NULL;
    }
    
    if ( ext == NULL )
    {
	fprintf(stderr, "bl_sam_fopen(): No filename extension on %s.\n", filename);
	return NULL;
    }

    if ( *mode == 'r' )
    {
	if ( strcmp(ext, ".gz") == 0 )
	{
// Big Sur zcat requires a .Z extension and CentOS 7 lacks gzcat
#ifdef __APPLE__
	    snprintf(cmd, XT_CMD_MAX_CHARS, "gzcat %s", filename);
#else
	    snprintf(cmd, XT_CMD_MAX_CHARS, "zcat %s", filename);
#endif
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "bzcat %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "xzcat %s", filename);
	    return popen(cmd, mode);
	}
	else if ( (strcmp(ext, ".bam") == 0) || (strcmp(ext, ".cram") == 0) )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --with-header %s", filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
    else    // "w"
    {
	if ( strcmp(ext, ".gz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "gzip -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "bzip2 -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "xz -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bam") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --bam %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".cram") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --cram %s", filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      Close a FILE stream with fclose() or pclose() as appropriate.
 *      Automatically determines the proper close function to call using
 *      S_ISFIFO on the stream stat structure.
 *
 *  Arguments:
 *      stream: The FILE structure to be closed
 *
 *  Returns:
 *      The value returned by fclose() or pclose()
 *
 *  See also:
 *      fopen(3), popen(3), gzip(1), bzip2(1), xz(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-05  Jason Bacon Derived from xt_fclose()
 ***************************************************************************/

int     bl_sam_fclose(FILE *stream)

{
    return xt_fclose(stream);
}

