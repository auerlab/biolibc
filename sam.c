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
#include "biostring.h"

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Name:
 *      bl_sam_skip_header() - Read past SAM header
 *
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over header lines in SAM input stream.  The FILE pointer
 *      sam_stream is advanced to the first character of the first line
 *      after the header.  The header is copied to a temporary file and and
 *      the function returns a FILE pointer to the header stream.
 *
 *  Arguments:
 *      sam_stream  FILE pointer to the open sam file
 *
 *  Returns:
 *      A FILE pointer to a temporary file containing a copy of the header
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_copy_header(3)
 *  
 *  History: 
 *  Date        Name        Modification
 *  2022-04-08  Jason Bacon Begin
 ***************************************************************************/

FILE    *bl_sam_skip_header(FILE *sam_stream)

{
    int     ch;
    FILE    *header_stream = tmpfile();

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like peak-classifier to replicate the
     *  header in output files.
     */
    
    while ( (ch = getc(sam_stream)) == '@' )
    {
	putc(ch, header_stream);
	do
	{
	    ch = getc(sam_stream);
	    putc(ch, header_stream);
	}   while ( (ch != '\n') && (ch != EOF) );
    }
    
    // Rewind to start of first non-header line
    if ( ch != EOF )
	ungetc(ch, sam_stream);
    rewind(header_stream);
    return header_stream;
}


/***************************************************************************
 *  Name:
 *      bl_sam_copy_header() - Copy SAM header to another stream
 *
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Copy SAM header from one FILE stream to another.  This is meant to
 *      be used in conjunction with bl_sam_skip_header(), which stores the
 *      header in a temporary file.
 *
 *  Arguments:
 *      header_stream   Open FILE stream of SAM header
 *      sam_stream      FILE stream to which header is copied
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE or BL_READ_* on failure
 *
 *  See also:
 *      bl_sam_skip_header(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-03  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_copy_header(FILE *header_stream, FILE *sam_stream)

{
    int     ch;
    
    rewind(header_stream);
    while ( (ch = getc(header_stream)) != EOF )
	if ( putc(ch, sam_stream) == EOF )
	    return BL_WRITE_FAILURE;
    rewind(header_stream);
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Name:
 *      bl_sam_read() - Read one SAM record
 *
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read next alignment (line) from a SAM stream.
 *
 *      If field_mask is not BL_SAM_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in alignment.
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
 *      sam_stream  A FILE stream from which to read the line
 *      alignment   Pointer to a bl_sam_t structure
 *      field_mask  Bit mask indicating which fields to store in alignment
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered after a complete feature
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      bl_sam_read(stdin, &alignment, BL_SAM_FIELD_ALL);
 *      bl_sam_read(sam_stream, &alignment,
 *                         BL_SAM_FIELD_QNAME|BL_SAM_FIELD_POS|BL_SAM_FIELD_TLEN);
 *
 *  See also:
 *      bl_sam_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_read(bl_sam_t *alignment, FILE *sam_stream,
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
	delim = tsv_read_field(sam_stream, alignment->qname,
			BL_SAM_QNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream, &len);
	*alignment->qname = '\0';
    }
    if ( delim == EOF )
	return BL_READ_EOF;

    // 2 Flag
    if ( field_mask & BL_SAM_FIELD_FLAG )
	delim = tsv_read_field(sam_stream, flag_str, BL_SAM_FLAG_MAX_DIGITS, &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading flag: %s.\n",
		flag_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_FLAG )
    {
	alignment->flag = strtoul(flag_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid position: %s\n",
		    flag_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	alignment->flag = 0;    // FIXME: Is there a better choice?
    
    // 3 RNAME
    if ( field_mask & BL_SAM_FIELD_RNAME )
	delim = tsv_read_field(sam_stream, alignment->rname,
			       BL_SAM_RNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream, &len);
	*alignment->rname = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading rname: %s.\n",
		alignment->rname);
	return BL_READ_TRUNCATED;
    }
    
    // 4 POS
    if ( field_mask & BL_SAM_FIELD_POS )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading pos: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_POS )
    {
	alignment->pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid position: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
	previous_pos = alignment->pos;
    }
    else
	alignment->pos = 0;
    
    // 5 MAPQ
    if ( field_mask & BL_SAM_FIELD_MAPQ )
	delim = tsv_read_field(sam_stream, mapq_str, BL_SAM_MAPQ_MAX_CHARS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading mapq: %s.\n",
		mapq_str);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_MAPQ )
    {
	alignment->mapq = strtoul(mapq_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid mapq: %s\n",
		    mapq_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	alignment->mapq = 0;
    
    // 6 CIGAR
    if ( field_mask & BL_SAM_FIELD_CIGAR )
	delim = tsv_read_field_malloc(sam_stream, &alignment->cigar,
			       &alignment->cigar_array_size,
			       &alignment->cigar_len);
    else
    {
	delim = tsv_skip_field(sam_stream, &len);
	alignment->cigar_len = 0;
	// Do not set to NULL or set array_size to 0.  Leave buffer
	// allocated for reuse.
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading cigar: %s.\n",
		alignment->cigar);
	return BL_READ_TRUNCATED;
    }
    
    // 7 RNEXT
    if ( field_mask & BL_SAM_FIELD_RNEXT )
	delim = tsv_read_field(sam_stream, alignment->rnext,
			       BL_SAM_RNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream, &len);
	*alignment->rnext = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading rnext: %s.\n",
		alignment->rnext);
	return BL_READ_TRUNCATED;
    }
    
    // 8 PNEXT
    if ( field_mask & BL_SAM_FIELD_PNEXT )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading pnext: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_PNEXT )
    {
	alignment->pnext = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid pnext: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	alignment->pnext = 0;
    
    // 9 TLEN
    if ( field_mask & BL_SAM_FIELD_TLEN )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading tlen: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_TLEN )
    {
	alignment->tlen = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid tlen: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	alignment->tlen = 0;
    
    // 10 SEQ
    if ( field_mask & BL_SAM_FIELD_SEQ )
	delim = tsv_read_field_malloc(sam_stream, &alignment->seq,
		    &alignment->seq_array_size, &alignment->seq_len);
    else
    {
	delim = tsv_skip_field(sam_stream, &alignment->seq_len);
	alignment->seq_len = 0;
	// Do not set to NULL or set array_size to 0.  Leave buffer
	// allocated for reuse.
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading seq: %s.\n",
		alignment->seq);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_SEQ )
    {
	// May be allocated by bl_sam_init() or bl_sam_copy()
	if ( alignment->seq == NULL )
	{
	    if ( (alignment->seq = xt_malloc(alignment->seq_len + 1,
		    sizeof(*alignment->seq))) == NULL )
	    {
		fprintf(stderr, "bl_sam_read(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
    }
    
    // 11 QUAL, should be last field
    if ( field_mask & BL_SAM_FIELD_QUAL )
    {
	delim = tsv_read_field_malloc(sam_stream, &alignment->qual,
		    &alignment->qual_array_size,
		    &alignment->qual_len);
    }
    else
    {
	delim = tsv_skip_field(sam_stream, &alignment->qual_len);
	alignment->qual_len = 0;
	// Do not set to NULL or set array_size to 0.  Leave buffer
	// allocated for reuse.
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading qual: %s.\n",
		alignment->qual);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_QUAL )
    {
	// May be allocated by bl_sam_init() or bl_sam_copy()
	if ( alignment->qual == NULL )
	{
	    if ( (alignment->qual = xt_malloc(alignment->qual_len + 1,
		    sizeof(*alignment->qual))) == NULL )
	    {
		fprintf(stderr, "bl_sam_read(): Could not allocate qual.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
    
	if ( (alignment->qual_len != 1) &&
	     (alignment->seq_len != alignment->qual_len) )
	    fprintf(stderr, "bl_sam_read(): Warning: qual_len != seq_len for %s,%" PRId64 "\n",
		    alignment->rname, alignment->pos);
    }
    
    // Some SRA CRAMs have 11 fields, most have 12
    // Discard everything after the 11th
    if ( delim == '\t' )
	while ( getc(sam_stream) != '\n' )
	    ;

    /*fprintf(stderr,"bl_sam_read(): %s,%" PRId64 ",%zu\n",
	    BL_SAM_RNAME(alignment), BL_SAM_POS(alignment),
	    BL_SAM_SEQ_LEN(alignment));*/
    
    return BL_READ_OK;
}


/***************************************************************************
 *  Name:
 *      bl_sam_copy() - Copy a SAM object
 *
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

    if ( src->cigar != NULL )
    {
	dest->cigar = strdup(src->cigar);
	if ( dest->cigar == NULL )
	{
	    fprintf(stderr, "bl_sam_copy(): Could not allocate cigar.\n");
	    exit(EX_UNAVAILABLE);
	}
	dest->cigar_array_size = src->cigar_len + 1;
	dest->cigar_len = src->cigar_len;
    }
    dest->cigar_array_size = src->cigar_array_size;
    dest->cigar_len = src->cigar_len;
    
    strlcpy(dest->rnext, src->rnext, BL_SAM_RNAME_MAX_CHARS + 1);
    dest->pnext = src->pnext;
    dest->tlen = src->tlen;

    if ( src->seq != NULL )
    {
	if ( (dest->seq = strdup(src->seq)) == NULL )
	{
	    fprintf(stderr, "bl_sam_copy(): Could not allocate seq.\n");
	    exit(EX_UNAVAILABLE);
	}
	dest->seq_array_size = src->seq_len + 1;
	dest->seq_len = src->seq_len;
    }
    dest->seq_array_size = src->seq_array_size;
    dest->seq_len = src->seq_len;
    
    //fprintf(stderr, "src->seq = %s %zu\n", src->seq, src->seq_len);
    //fprintf(stderr, "src->qual = %s %zu\n", src->qual, src->qual_len);
    
    /*
     *  qual is an optional field.  If skipped using tsv_skip_field()
     *  qual_len will be non-zero, but qual will be NULL.
     */
    if ( src->qual != NULL )
    {
	if ( (dest->qual = strdup(src->qual)) == NULL )
	{
	    fprintf(stderr, "bl_sam_copy(): Could not allocate qual.\n");
	    exit(EX_UNAVAILABLE);
	}
	dest->qual_array_size = src->qual_len + 1;
	dest->qual_len = src->qual_len;
    }
    dest->qual_array_size = src->qual_array_size;
    dest->qual_len = src->qual_len;
}


/***************************************************************************
 *  Name:
 *      bl_sam_free() - Destroy a SAM object
 *
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_sam_read() or
 *      bl_sam_init().
 *
 *  Arguments:
 *      alignment   Pointer to bl_sam_t structure to be freed.
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_init(3), bl_sam_copy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_free(bl_sam_t *alignment)

{
    if ( alignment->cigar != NULL )
	free(alignment->cigar);
    if ( alignment->seq != NULL )
	free(alignment->seq);
    if ( alignment->qual != NULL )
	free(alignment->qual);
}


/***************************************************************************
 *  Name:
 *      bl_sam_init() - Initialize all fields of a SAM object
 *
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
 *      alignment   Pointer to bl_sam_t structure to initialize
 *      seq_len     Length of sequence and quality strings
 *      field_mask  Bit mask indicating which fields will be used
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_free(3), bl_sam_copy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_init(bl_sam_t *alignment)

{
    *alignment->qname = '\0';
    alignment->flag = 0;
    *alignment->rname = '\0';
    alignment->pos = 0;
    alignment->mapq = 0;
    alignment->cigar = NULL;
    *alignment->rnext = '\0';
    alignment->pnext = 0;
    alignment->tlen = 0;
    alignment->seq = NULL;
    alignment->qual = NULL;
    alignment->seq_array_size = 0;
    alignment->seq_len = 0;
    alignment->qual_array_size = 0;
    alignment->qual_len = 0;
}


/***************************************************************************
 *  Name:
 *      bl_sam_write() - Write a SAM object to a file stream
 *
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write an alignment (line) to a SAM stream.
 *
 *      If field_mask is not BL_SAM_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate placeholder such as '.'
 *      rather than stored in alignment.  Possible mask values are:
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
 *      sam_stream  A FILE stream to which to write the line
 *      alignment   Pointer to a bl_sam_t structure
 *      field_mask  Bit mask indicating which fields to store in alignment
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

int     bl_sam_write(bl_sam_t *alignment, FILE *sam_stream,
			   sam_field_mask_t field_mask)

{
    int     count;
    
    // FIXME: Respect field_mask
    count = fprintf(sam_stream, "%s\t%u\t%s\t%" PRId64
		    "\t%u\t%s\t%s\t%" PRId64 "zu\t%lu\t%s\t%s\t%zu\t%zu\n",
		    alignment->qname,
		    alignment->flag,
		    alignment->rname,
		    alignment->pos,
		    alignment->mapq,
		    alignment->cigar,
		    alignment->rnext,
		    alignment->pnext,
		    alignment->tlen,
		    alignment->seq,
		    alignment->qual,
		    alignment->seq_len,
		    alignment->qual_len);
    return count;
}


/***************************************************************************
 *  Name:
 *      bl_sam_fopen() - Open a SAM/BAM/CRAM file
 *
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      Open a raw SAM file using fopen() or a compressed
 *      SAM file, or BAM or CRAM file using popen().  If the
 *      file extension is .bam or .cram, or samtools_args is not
 *      NULL or "", data will be piped through "samtools view" with
 *      the given samtools_args as arguments.  The flag --with-header
 *      is always added for consistency with the case of reading a
 *      raw SAM file without piping through samtools.  Programs that
 *      don't want the header can filter it out by other means, such
 *      as bl_sam_skip_header(3).
 *
 *      bl_sam_fopen() must be used in conjunction with
 *      bl_sam_fclose() to ensure that fclose() or pclose() is called where
 *      appropriate.
 *
 *  Arguments:
 *      filename:       Name of the file to be opened
 *      mode:           "r" or "w", passed to fopen() or popen()
 *      samtools_args   Flags to pass to samtools view
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

FILE    *bl_sam_fopen(const char *filename, const char *mode,
		      char *samtools_args)

{
    char    *ext = strrchr(filename, '.'),
	    cmd[XT_CMD_MAX_CHARS + 1];
    struct stat sb;
    
    if ( samtools_args == NULL )
	samtools_args = "";
    
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

    // popen() does not return NULL when the file does not exist
    if ( stat(filename, &sb) != 0 )
	return NULL;
    
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
	else if ( (strcmp(ext, ".bam") == 0) || (strcmp(ext, ".cram") == 0) 
		  || ! strblank(samtools_args) )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --with-header %s %s",
		    samtools_args, filename);
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
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --bam --with-header %s %s",
		    samtools_args, filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".cram") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --cram --with-header %s %s",
		    samtools_args, filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
}


/***************************************************************************
 *  Name:
 *      bl_sam_fclose() - Close a stream opened by bl_sam_fopen(3)
 *
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



/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Name:
 *      bl_sam_gff3_overlap() - Compute SAM/GFF3 overlap
 *
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the amount of overlap between a GFF3 feature and a SAM
 *      alignment.
 *  
 *  Arguments:
 *      alignment   Pointer to a bl_sam_t object
 *      feature     Pointer to a bl_gff3_t object
 *
 *  Returns:
 *      The number of bases of overlap between the feature and alignment.
 *      A zero or negative return value indicates no overlap.
 *
 *  See also:
 *      bl_gff3_sam_overlap(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-07  Jason Bacon Begin
 ***************************************************************************/

int64_t bl_sam_gff3_overlap(bl_sam_t *alignment, bl_gff3_t *feature)

{
    return bl_gff3_sam_overlap(feature, alignment);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Name:
 *      bl_sam_gff3_cmp() - Compare positions of SAM and GFF3 records
 *
 *  Library:
 *      #include <sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Compare the positions of a SAM alignment and a GFF feature and
 *      return a status value much like strcmp().  0 is returned if the
 *      alignment and feature overlap.  A value < 0 is returned if the
 *      alignment is entirely "before" the feature, i.e. it is on an
 *      earlier chromosome according to bl_chrom_name_cmp(3), or on the
 *      same chromosome at a lower position.  A value > 0 is returned
 *      if the alignment is entirely "after" the feature, i.e. on a later
 *      chromosome or same chromosome and higher position.
 *
 *      This function is mainly intended for programs that sweep properly
 *      sorted GFF and SAM files locating overlaps in a single pass.
 *  
 *      A converse function, bl_gff3_sam_cmp(3) is also provided so that
 *      the programmer can choose the more intuitive interface.
 *  
 *  Arguments:
 *      alignment   Pointer to a bl_sam_t object
 *      feature     Pointer to a bl_gff3_t object
 *
 *  Returns:
 *      A value < 0 if the the alignment is entirely before the feature
 *      A value > 0 if the the alignment is entirely after the feature
 *      0 if the alignment and the feature overlap
 *
 *  See also:
 *      bl_gff3_sam_cmp(3), bl_chrom_name_cmp(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-06  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_gff3_cmp(bl_sam_t *alignment, bl_gff3_t *feature)

{
    int     status = bl_chrom_name_cmp(BL_SAM_RNAME(alignment),
					    BL_GFF3_SEQID(feature));
    
    if ( status != 0 )
	// Different chromosomes
	return status;
    else if ( BL_SAM_POS(alignment) + BL_SAM_SEQ_LEN(alignment) - 1
		< BL_GFF3_START(feature) )
	// Alignment ends before the start of feature
	return -1;
    else if ( BL_SAM_POS(alignment) > BL_GFF3_END(feature) )
	// Alignment starts after the end of feature
	return 1;
    else
	// Overlap
	return 0;
}

