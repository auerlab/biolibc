#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <sysexits.h>
#include <sys/param.h>  // MIN()
#include <xtend/dsv.h>
#include <xtend/mem.h>
#include "fastq.h"
#include "biolibc.h"

/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a FASTQ record from a FILE stream.  Each record must begin
 *      with a description line (beginning with '@'), which is then
 *      followed by one or more lines of sequence data, a separator line
 *      beginning with '+', and a line of quality scores.  The end of the
 *      sequence is marked either by the next description line or EOF.
 *      If desc_len and seq_len are 0 (e.g. the structure is initialized
 *      with BL_FASTQ_INIT or bl_fastq_init(3), has been freed with
 *      bl_fastq_free(3), then memory is allocated for each line.
 *
 *      Otherwise, the existing allocated buffers are reused.  Hence, when
 *      reading many FASTQ records of the same length, only one allocation
 *      is needed.  In any case, the buffers are automatically enlarged if
 *      they become full and automatically trimmed to the actual data size
 *      after reading is complete.
 *
 *      Buffer memory should be freed as soon as possible by calling
 *      bl_fastq_free(3).
 *  
 *  Arguments:
 *      fastq_stream    FILE stream from which FASTQ data are read
 *      record          Pointer to a bl_fastq_t structure to receive data
 *
 *  Returns:
 *      BL_READ_OK upon successful read of description and sequence
 *      BL_READ_BAD_DATA if something is amiss with input format
 *      BL_READ_EOF if no more data are available
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_write(3), bl_fastq_read(3), bl_fastq_write(3),
 *      bl_fastq_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastq_read(bl_fastq_t *record, FILE *fastq_stream)

{
    int     ch,
	    last_ch;
    size_t  len;
    
    /* Skip comment lines */
    while ( ((ch = getc(fastq_stream)) == ';') && (ch != EOF) )
	while ( ((ch = getc(fastq_stream)) != '\n') && (ch != EOF) )
	    ;
    
    if ( ch == EOF )
	return BL_READ_EOF;
    
    /* Every record should begin with a '@' */
    if ( ch == '@' )    // Desc
    {
	/*
	 *  Read description
	 */

	ungetc(ch, fastq_stream);
	ch = dsv_read_field_malloc(fastq_stream, &record->desc,
			    &record->desc_array_size, "", &record->desc_len);
	if ( record->desc == NULL )
	{
	    fprintf(stderr, "bl_fastq_read(): Could not allocate desc.\n");
	    exit(EX_UNAVAILABLE);
	}
	
	/* Should not encounter EOF while reading description line */
	/* Every description should be followed by at least one seq line */
	if ( ch == EOF )
	{
	    fprintf(stderr, "bl_fastq_read(): Record truncated in desc %s.\n",
		    record->desc);
	    return BL_READ_TRUNCATED;
	}
	else if ( ch != '\n' )
	{
	    fprintf(stderr, "bl_fastq_read(): Bad data after desc %s\n", record->desc);
	    return BL_READ_BAD_DATA;
	}

	/*
	 *  Read sequence lines.  May span multiple lines so can't use
	 *  dsv_read_field_malloc().
	 */
	
	if ( record->seq_array_size == 0 )
	{
	    // Easily hold a short read sequence and minimize reallocs for long
	    record->seq_array_size = 1024;
	    record->seq = xt_malloc(record->seq_array_size,
				    sizeof(*record->seq));
	    if ( record->seq == NULL )
	    {
		fprintf(stderr, "bl_fastq_read(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	
	len = 0;
	do
	{
	    if ( ch != '\n' )
		record->seq[len++] = ch;
	    if ( len == record->seq_array_size - 1 )
	    {
		record->seq_array_size *= 2;
		record->seq = xt_realloc(record->seq, record->seq_array_size,
		    sizeof(*record->seq));
		if ( record->seq == NULL )
		{
		    fprintf(stderr, "bl_fastq_read(): Could not reallocate seq.\n");
		    exit(EX_UNAVAILABLE);
		}
	    }
	    last_ch = ch;
	}   while ( ((ch = getc(fastq_stream)) != '+') && (ch != EOF) );
	record->seq[len] = '\0';
	record->seq_len = len;

	if ( last_ch != '\n' )
	    fprintf(stderr, "bl_fasta_read(): Missing newline at end of qual %s.\n",
		    record->qual);

	/* 
	 * Trim array.  realloc() can carry a significant cost, but it does
	 * not affect overall performance here, probably because I/O is
	 * the major bottleneck.
	 */
	if ( record->seq_array_size != record->seq_len + 1 )
	{
	    record->seq_array_size = record->seq_len + 1;
	    record->seq = xt_realloc(record->seq, record->seq_array_size,
		sizeof(*record->seq));
	}

	/* Should not encounter EOF while reading sequence lines */
	/* Every sequence should be followed by a + separator line */
	if ( ch == EOF )
	{
	    fprintf(stderr, "bl_fastq_read(): Record truncated in seq %s.\n",
		    record->seq);
	    return BL_READ_TRUNCATED;
	}
	else if (ch != '+')
	{
	    fprintf(stderr, "bl_fasq_read(): Bad data after seq %s\n", record->seq);
	    return BL_READ_BAD_DATA;
	}
	// Put '+' back so it's read into plus field
	ungetc(ch, fastq_stream);
	    
	/*
	 *  Read + separator
	 */
	
	ch = dsv_read_field_malloc(fastq_stream, &record->plus,
			    &record->plus_array_size, "", &record->plus_len);
	if ( record->plus == NULL )
	{
	    fprintf(stderr, "bl_fastq_read(): Could not allocate plus.\n");
	    exit(EX_UNAVAILABLE);
	}
	
	/* Should not encounter EOF while reading plus line */
	/* Every plus should be followed by at least one qual line */
	if ( ch == EOF )
	{
	    fprintf(stderr, "bl_fastq_read(): Record truncated in plus %s.\n",
		    record->plus);
	    return BL_READ_TRUNCATED;
	}
	else if ( ch != '\n' )
	{
	    fprintf(stderr, "bl_fasq_read(): Bad data after plus %s\n",
		    record->plus);
	    return BL_READ_BAD_DATA;
	}

	/*
	 *  Read quality string.  May span multiple lines so can't use
	 *  dsv_read_field_malloc().
	 */

	// FIXME: This could be problematic with bad data where qual len
	// doesn't match seq len
	if ( record->qual_array_size == 0 )
	{
	    // Must match sequence len and no need to trim
	    record->qual_array_size = record->seq_array_size;
	    record->qual = xt_malloc(record->qual_array_size, sizeof(*record->qual));
	    if ( record->qual == NULL )
	    {
		fprintf(stderr, "bl_fastq_read(): Could not allocate qual.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	
	len = 0;
	do
	{
	    /* Read at least one full line, since '@' can be a quality score */
	    while ( ((ch = getc(fastq_stream)) != '\n') && (ch != EOF) )
	    {
		record->qual[len++] = ch;
		if ( len == record->qual_array_size - 1 )
		{
		    record->qual_array_size *= 2;
		    record->qual = xt_realloc(record->qual, record->qual_array_size,
			sizeof(*record->qual));
		    if ( record->qual == NULL )
		    {
			fprintf(stderr, "bl_fastq_read(): Could not reallocate qual.\n");
			exit(EX_UNAVAILABLE);
		    }
		}
	    }
	    last_ch = ch;
	}   while ( ((ch = getc(fastq_stream)) != '@') && (ch != EOF) );
	record->qual[len] = '\0';
	record->qual_len = len;
	
	if ( last_ch != '\n' )
	    fprintf(stderr, "bl_fasta_read(): Missing newline at end of qual %s.\n",
		    record->qual);

	/*
	 *  This is where EOF should occur since we read past newlines
	 *  No need to trim since qual must be the same size as seq
	 */

	// Put '@' back so it's read into next desc
	if ( ch == '@' )
	    ungetc(ch, fastq_stream);

	return BL_READ_OK;
    }
    else
	return BL_READ_BAD_DATA;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a FASTQ record to the specified FILE stream, writing at most
 *      max_line_len sequence characters per line.  The special value
 *      BL_FASTQ_LINE_UNLIMITED indicates no line length limit.
 *  
 *  Arguments:
 *      fastq_stream    FILE stream to which data are written
 *      record          Pointer to a bl_fastq_t structure to be written
 *      max_line_len    Maximum length of a sequence line in output
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE if a write error occurs.
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_read(3), bl_fastq_read(3), bl_fastq_write(3),
 *
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastq_write(bl_fastq_t *record, FILE *fastq_stream,
		       size_t max_line_len)

{
    size_t  c;
    int     save_ch = ' ';  // Silence false uninit warning on CentOS
    
    if ( fprintf(fastq_stream, "%s\n", record->desc) < 0 )
	return BL_WRITE_FAILURE;
    
    if ( max_line_len == BL_FASTQ_LINE_UNLIMITED )
    {
	if ( fprintf(fastq_stream, "%s\n", record->seq) < 0 )
	    return BL_WRITE_FAILURE;
    }
    else
    {
	for (c = 0; c < record->seq_len; c += max_line_len)
	{
	    // Temporarily null-terminate segment of string to be printed
	    if ( record->seq_len - c > max_line_len )
	    {
		save_ch = record->seq[c + max_line_len];
		record->seq[c + max_line_len] = '\0';
	    }
	    
	    // Print segment
	    if ( fprintf(fastq_stream, "%s\n", record->seq + c) < 0 )
		return BL_WRITE_FAILURE;
    
	    // Remove temporary null-termination
	    if ( record->seq_len - c > max_line_len )
		record->seq[c + max_line_len] = save_ch;
	}
    }
    
    if ( fprintf(fastq_stream, "%s\n", record->plus) < 0 )
	return BL_WRITE_FAILURE;
    
    if ( max_line_len == BL_FASTQ_LINE_UNLIMITED )
    {
	if ( fprintf(fastq_stream, "%s\n", record->qual) < 0 )
	    return BL_WRITE_FAILURE;
    }
    else
    {
	for (c = 0; c < record->qual_len; c += max_line_len)
	{
	    // Temporarily null-terminate segment of string to be printed
	    if ( record->qual_len - c > max_line_len )
	    {
		save_ch = record->qual[c + max_line_len];
		record->qual[c + max_line_len] = '\0';
	    }
    
	    // Print segment
	    if ( fprintf(fastq_stream, "%s\n", record->qual + c) < 0 )
		return BL_WRITE_FAILURE;
    
	    // Remove temporary null-termination
	    if ( record->qual_len - c > max_line_len )
		record->qual[c + max_line_len] = save_ch;
	}
    }
    
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_fastq_read()
 *  
 *  Arguments:
 *      record  Pointer to a previously populated bl_fastq_t structure
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_read(3), bl_fastq_write(3)
 *      bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastq_free(bl_fastq_t *record)

{
    free(record->seq);
    free(record->desc);
    free(record->plus);
    free(record->qual);
    record->desc = record->seq = record->plus = record->qual = NULL;
    record->desc_array_size = record->seq_array_size = 
	record->plus_array_size = record->qual_array_size = 0;
    record->desc_len = record->seq_len = 
	record->plus_len = record->qual_len = 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_fastq_t structure.  This must be done before
 *      passing it to bl_fastq_read() for the first time, so that
 *      bl_fastq_read() will know to allocate memory for the fields.
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastq_t structure to initialize.
 *
 *  Examples:
 *      bl_fastq_t  rec;
 *
 *      bl_fastq_init(&rec);
 *      bl_fastq_read(stdin, &rec);
 *
 *  See also:
 *      bl_fastq_read(3), bl_fastq_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastq_init(bl_fastq_t *record)

{
    record->desc = record->seq = record->plus = record->qual = NULL;
    record->desc_array_size = record->seq_array_size = 
	record->plus_array_size = record->qual_array_size = 0;
    record->desc_len = record->seq_len = 
	record->plus_len = record->qual_len = 0;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate adapter (or a portion thereof if the end of the read is
 *      reached) in a FASTQ sequence.
 *
 *      The content of adapter is assumed to be all upper case.  This
 *      improves speed by avoiding millions of redundant toupper()
 *      conversions on the same string.  Use strlupper(3) or strupper(3)
 *      before calling this function if necessary.
 *
 *      A minimum of min_match bases must match between the adapter and
 *      the read.  This mainly matters near the end of a read, where
 *      remaining bases are fewer than the length of the adapter.
 *
 *      A maximum of max_mismatch_percent mismatched bases are tolerated
 *      to allow for read errors in the adapter, which do happen, albeit
 *      rarely.  This is taken as a percent of the adapter length, or the
 *      remaining bases in the sequence, whichever is smaller.
 *      Note that the NUMBER of mismatched bases tolerated is
 *      truncated from the percent calculation.  E.g. using 10% tolerance,
 *      0 mismatched base is tolerated among 9 total bases, 1 among 10.
 *
 *      Higher values of max_mismatch_percent will results in slightly
 *      longer run times, more adapters detected, and a higher risk of
 *      false-positives (falsely identifying natural sequences as adapters).
 *
 *      Indels (insertions and deletions) are not currently handled.
 *
 *      Note that adapter matching is not an exact science.  We cannot remove
 *      every adapter without also removing many natural sequences, since
 *      it is impossible to know whether any given sequence is really an
 *      adapter or naturally occurring.  The best we can do is guestimate
 *      what will remove the most adapters and fewest natural sequences.
 *      It is also not usually important to remove every adapter, but only to
 *      minimize adapter contamination.  Failing to align a small percentage
 *      of sequences due to adapter contamination will not change the story
 *      told by the downstream analysis.  Nor will erroneously trimming off
 *      the 3' end of a small percentage of reads containing natural
 *      sequences resembling adapters.  Just trimming exact matches of
 *      the adapter sequence will generally remove 99% or more of the
 *      adapter contamination and minimize false-positives.  Tolerating
 *      1 or 2 differences has been shown to do slightly better overall.
 *
 *  Arguments:
 *      read        FASTQ read to be searched
 *      adapter     Adapter sequence to be located
 *      min_match   Minimum number of characters to match in adapter
 *      max_mismatch_percent    Maximum percent of unequal bases in overlap
 *
 *  Returns:
 *      Index of adapter sequence if found, index of null terminator otherwise
 *
 *  Examples:
 *      bl_fastq_t  read;
 *      char        *adapter;
 *      size_t      index;
 *
 *      index = bl_fastq_find_adapter_smart(&read, adapter, 3, 10);
 *      if ( BL_FASTQ_SEQ_AE(&read, index) != '\0' )
 *          bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_fastq_find_adapter_exact(3), bl_fastq_3p_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_find_adapter_smart(const bl_fastq_t *read,
	    const char *adapter, size_t min_match,
	    unsigned max_mismatch_percent)

{
    // The strlen() looks expensive, but tests show that eliminating it
    // doesn't reduce run time measurably
    size_t      match, mismatch, max_mismatch,
		adapter_len = strlen(adapter),
		start, rc, ac,
		md, adapter_mm;
    
    // Start at 5' end assuming 5' adapters already removed
    // Cutadapt uses a semiglobal alignment algorithm to find adapters.
    // Not sure what the benefit of this is over exact matching. I would
    // assume that errors in adapter sequences are extremely rare.
    // https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm

    // Convert max mismatch percentage to a divisor for the string len
    md = 100 / max_mismatch_percent;
    adapter_mm = adapter_len / md;  // Max mismatch based on adapter len
    // Could stop at read->seq_len - min_match, but the extra math
    // outweights the few iterations saved
    for (start = 0; start < read->seq_len; ++start)
    {
	// Terminate loop as soon as max_mismatch is reached, before
	// checking other conditions
	max_mismatch = MIN((read->seq_len - start) / md, adapter_mm);
	for (rc = start, ac = 0, match = mismatch = 0;
	     (mismatch <= max_mismatch) &&
	     (adapter[ac] != '\0') && (read->seq[rc] != '\0'); ++rc, ++ac)
	{
	    if ( toupper(read->seq[rc]) != adapter[ac] )
		++mismatch;
	}
	if ( mismatch <= max_mismatch )
	{
	    match = ac - mismatch;
	    if ( match >= min_match )
		return start;
	}
    }
    return read->seq_len;   // Location of '\0' terminator
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate adapter (or a portion thereof if the end of the read is
 *      reached) in a FASTQ sequence.
 *
 *      Functionally equivalent to but faster than
 *      bl_fastq_find_adapter_exact(&read, adapter, min_match, 0);
 *
 *      The content of adapter is assumed to be all upper case.  This
 *      improves speed by avoiding millions of redundant toupper()
 *      conversions on the same string.  Use strlupper(3) or strupper(3)
 *      before calling this function if necessary.
 *
 *      A minimum of min_match bases must match between the adapter and
 *      the read.  This mainly matters near the end of a read, where
 *      remaining bases are fewer than the length of the adapter.
 *
 *      Note that adapter matching is not an exact science.  We cannot remove
 *      every adapter without also removing many natural sequences, since
 *      it is impossible to know whether any given sequence is really an
 *      adapter or naturally occurring.  The best we can do is guestimate
 *      what will remove the most adapters and fewest natural sequences.
 *      It is also not usually important to remove every adapter, but only to
 *      minimize adapter contamination.  Failing to align a small percentage
 *      of sequences due to adapter contamination will not change the story
 *      told by the downstream analysis.  Nor will erroneously trimming off
 *      the 3' end of a small percentage of reads containing natural
 *      sequences resembling adapters.  Just trimming exact matches of
 *      the adapter sequence will generally remove 99% or more of the
 *      adapter contamination and minimize false-positives.  Tolerating
 *      1 or 2 differences has been shown to do slightly better overall.
 *
 *  Arguments:
 *      read        FASTQ read to be searched
 *      adapter     Adapter sequence to be located
 *      min_match   Minimum number of characters to match in adapter
 *      max_mismatch_percent    Ignored: Dummy argument for compatibility
 *                              with bl_fastq_fund_adapter_smart(3)
 *
 *  Returns:
 *      Index of adapter sequence if found, index of null terminator otherwise
 *
 *  Examples:
 *      bl_fastq_t  read;
 *      char        *adapter;
 *      size_t      index;
 *
 *      index = bl_fastq_find_adapter_exact(&read, adapter, 3, 10);
 *      if ( BL_FASTQ_SEQ_AE(&read, index) != '\0' )
 *          bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_fastq_find_adapter_smart(3), bl_fastq_3p_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_find_adapter_exact(const bl_fastq_t *read,
	    const char *adapter, size_t min_match,
	    unsigned max_mismatch_percent)

{
    size_t  start, rc, ac;
    
    // Start at 5' end assuming 5' adapters already removed
    // Cutadapt uses a semiglobal alignment algorithm to find adapters.
    // Not sure what the benefit of this is over exact matching. I would
    // assume that errors in adapter sequences are extremely rare.
    // https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm

    // Could stop at read->seq_len - min_match, but the extra math
    // outweights the few iterations saved
    for (start = 0; start < read->seq_len; ++start)
    {
	for (rc = start, ac = 0; (toupper(read->seq[rc]) == adapter[ac]) &&
	     (adapter[ac] != '\0'); ++rc, ++ac)
	    ;
	if ( (adapter[ac] == '\0') || ((read->seq[rc] == '\0') && (ac >= min_match)) )
	    return start;
    }
    return read->seq_len;   // Location of '\0' terminator
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Trim the 3' end of a FASTQ sequence and qualit string at location
 *      new_len.
 *  
 *  Arguments:
 *      read        FASTQ read to be trimmed
 *      new_len     New length and location of the null terminators
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if new_len is between 0 and original length,
 *      BL_FASTQ_DATA_INVALID otherwise.
 *
 *  Examples:
 *      bl_fastq_t  read;
 *      char        *adapter;
 *      size_t      index;
 *
 *      index = bl_fastq_find_adapter_smart(&read, adapter, 3, 10);
 *      if ( BL_FASTQ_SEQ_AE(&read, index) != '\0' )
 *          bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_fastq_find_adapter_smart(3), bl_fastq_find_adapter_exact(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_3p_trim(bl_fastq_t *read, size_t new_len)

{
    if ( (new_len >= 0) && (new_len <= read->seq_len) )
    {
	read->seq_len = read->qual_len = new_len;
	read->seq[new_len] = read->qual[new_len] = '\0';
	// FIXME: realloc?
	return BL_FASTQ_DATA_OK;
    }
    else
	return BL_FASTQ_DATA_INVALID;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate start of a low-quality 3' end in a FASTQ read.  This
 *      function uses the same algorithm as fastq and cutadapt as of the
 *      time of writing.  Namely, it starts at the 3' end of the quality
 *      string and sums (base quality - minimum quality) while moving in
 *      the 5' direction.  This sum will be < 0 as long as the average
 *      base quality is < minimum quality.  It also keeps track of where
 *      the minimum of this sum occurs.  When the sum become > 0, we have
 *      reached a point where the average quality of the 3' end is
 *      satisfactory, and it is assumed it will remain that way if we
 *      continue in the 5' direction.  ( Illumina reads tend to drop in
 *      quality near the 3' end. )  The location of the minimum sum is
 *      then returned, since the average quality of everything in the 5'
 *      direction must be satisfactory.
 *  
 *  Arguments:
 *      read        FASTQ read to be searched
 *      min_qual    Minimum quality of bases to keep
 *      phred_base  Offset into the ISO character set used by PHRED scores
 *                  (usually 33 for modern data)
 *
 *  Returns:
 *      Index of first low-quality base at the 3' end if found,
 *      index of NULL terminator otherwise
 *
 *  Examples:
 *      bl_fastq_t  read;
 *      
 *      ...
 *      index = bl_fastq_find_3p_low_qual(&read, 20, 33);
 *      bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_fastq_find_adapter_smart(3), bl_fastq_find_adapter_exact(3),
 *      bl_fastq_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_find_3p_low_qual(const bl_fastq_t *read, unsigned min_qual,
			unsigned phred_base)

{
    ssize_t     c,
		cut_pos;
    long        sum,
		min_sum;
    
    /*
     *  Use same algorithm as BWA/cutadapt
     *  https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
     *  score-minimum will be negative for bases we want to remove
     *  Sum score-minimum for each base starting at end until sum > 0
     *  Use the position of the minimum sum as the trim point
     *  Verified using 42, 40, 26, 27, 8, 7, 11, 4, 2, 3 example from link
     */

    if ( read->seq_len != read->qual_len )
    {
	fprintf(stderr, "bl_fastq_find_3p_low_qual(): qual_len != seq_len.\n");
	exit(EX_DATAERR);
    }
    
    sum = min_sum = 0;
    c = read->qual_len - 1;
    cut_pos = read->seq_len;
    while ( (c >= 0) && (sum <= 0) )
    {
	// Revert promotions to unsigned
	sum = (long)sum + read->qual[c] - phred_base - min_qual;
	if ( sum < min_sum )
	{
	    // fprintf(stderr, "%zu %c %c %ld\n", c, read->seq[c], read->qual[c], sum);
	    min_sum = sum;
	    cut_pos = c;
	}
	--c;
    }
    // fprintf(stderr, "Returning %zd\n", cut_pos);
    return cut_pos;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Compare the read names of two FASTQ reads.  This is useful when
 *      processing paired-end data, which must be kept in-sync.  I.e.
 *      if a sequence if removed from a 5' file, the same sequence should
 *      be removed from the 3' file whether or not it meets quality
 *      minimums.
 *  
 *  Arguments:
 *      read1, read2    FASTQ reads to compare   
 *
 *  Returns:
 *      0 if read1 and read2 have the same name
 *      < 0 if read1 name is lexically less than read2 name
 *      > 0 if read1 name is lexically greater than read2 name
 *
 *  Examples:
 *      s1 = bl_fastq_read(&fastq_rec[0], tp->instream1);
 *      s2 = bl_fastq_read(&fastq_rec[1], tp->instream2);
 *      if ( (s1 == BL_READ_OK) && (s2 == BL_READ_OK) )
 *      {
 *          if ( bl_fastq_name_cmp(&fastq_rec[0], &fastq_rec[1]) != 0 )
 *          {
 *              fprintf(stderr, "fastq-trim: Paired files out of sync.\n");
 *              trim_close_files(tp);
 *              exit(EX_DATAERR);
 *          }
 *          ...
 *      }
 *
 *  See also:
 *      bl_fastq_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_name_cmp(bl_fastq_t *read1, bl_fastq_t *read2)

{
    // FIXME: This is a hack based on test data.  Find out how to 
    // properly compare names in arbitrary FASTQ files
    // Description up to first space char is the same for R1 and R2 files
    char    *p1 = strchr(BL_FASTQ_DESC(read1), ' ');
    char    *p2 = strchr(BL_FASTQ_DESC(read2), ' ');
    int     save_p1, save_p2, status;
    
    // Temporarily null-terminate descriptions at first space char
    // Not thread-safe
    save_p1 = *p1;
    save_p2 = *p2;
    *p1 = *p2 = '\0';
    status = strcmp(BL_FASTQ_DESC(read1), BL_FASTQ_DESC(read2));
    *p1 = save_p1;
    *p2 = save_p2;
    return status;
}

