

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Name:
 *      bl_fastq_5p_trim() - Trim 5' end of a FASTQ object
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Trim the 5' end of a FASTQ sequence and quality string at location
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
 *          bl_fastq_5p_trim(&read, index);
 *
 *  See also:
 *      bl_fastq_find_adapter_smart(3), bl_fastq_find_adapter_exact(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_5p_trim(bl_fastq_t *read, size_t new_len)

{
    if ( (new_len >= 0) && (new_len <= read->seq_len) )
    {
	size_t  trimmed = read->seq_len - new_len;
	read->seq_len = read->qual_len = new_len;
	// Include null byte in move
	memmove(read->seq, read->seq + trimmed, new_len + 1);
	// FIXME: realloc?
	return BL_FASTQ_DATA_OK;
    }
    else
	return BL_FASTQ_DATA_INVALID;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Name:
 *      bl_fastq_find_5p_low_qual() - Find start of low-quality 5' end
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate start of a low-quality 5' end in a FASTQ read.  This
 *      function uses the same algorithm as fastq and cutadapt as of the
 *      time of writing.  Namely, it starts at the 5' end of the quality
 *      string and sums (base quality - minimum quality) while moving in
 *      the 5' direction.  This sum will be < 0 as long as the average
 *      base quality is < minimum quality.  It also keeps track of where
 *      the minimum of this sum occurs.  When the sum become > 0, we have
 *      reached a point where the average quality of the 5' end is
 *      satisfactory, and it is assumed it will remain that way if we
 *      continue in the 5' direction.  ( Illumina reads tend to drop in
 *      quality near the 5' end. )  The location of the minimum sum is
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
 *      Index of first low-quality base at the 5' end if found,
 *      index of NULL terminator otherwise
 *
 *  Examples:
 *      bl_fastq_t  read;
 *      
 *      ...
 *      index = bl_fastq_find_5p_low_qual(&read, 20, 33);
 *      bl_fastq_5p_trim(&read, index);
 *
 *  See also:
 *      bl_fastq_find_adapter_smart(3), bl_fastq_find_adapter_exact(3),
 *      bl_fastq_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_find_5p_low_qual(const bl_fastq_t *read, unsigned min_qual,
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
	fprintf(stderr, "bl_fastq_find_5p_low_qual(): qual_len != seq_len.\n");
	exit(EX_DATAERR);
    }
    
    sum = min_sum = 0;
    c = 0; // read->qual_len - 1;
    cut_pos = -1; // read->seq_len;
    while ( (c < read->seq_len) && (sum <= 0) )
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

