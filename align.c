#include <stdio.h>
#include <ctype.h>
#include <xtend/math.h> // XT_MIN()
#include "align.h"

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/align.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate the leftmost (farthest 5') match for sequence little within
 *      sequence big, tolerating the given percentage of mismatched bases.
 *
 *      The content of little is assumed to be all upper case.  This
 *      improves speed by avoiding numerous redundant toupper()
 *      conversions on the same string, assuming multiple big strings will
 *      be searched for little, as in adapter removal and read mapping.
 *      Use strlupper(3) or strupper(3) before calling this function if
 *      necessary.
 *
 *      A minimum of min_match bases must match between little and
 *      big.  This mainly matters near the end of big, where
 *      remaining bases are fewer than the length of little.
 *
 *      A maximum of max_mismatch_percent mismatched bases are tolerated
 *      to allow for read errors. This is taken as a percent of little, or
 *      the same percent of remaining bases in big, whichever is smaller.
 *      Note that the NUMBER of allowed mismatched bases tolerated is
 *      truncated from the percent calculation.  E.g. using 10% tolerance,
 *      0 mismatched bases are tolerated among 9 total bases, or 1 mismatch
 *      among 10 total.
 *
 *      Higher values of max_mismatch_percent will results in slightly
 *      longer run times, more alignments detected, and a higher risk of
 *      false-positives (falsely identifying other big sequences as matching
 *      little.
 *
 *      Indels (insertions and deletions) are not currently handled.
 *
 *      Note that alignment is not an exact science.  We cannot detect every
 *      true little sequence without falsely detecting other sequences, since
 *      it is impossible to know whether any given sequence is really from
 *      the source of interest (e.g. an adapter) or naturally
 *      occurring from another source.  The best we can do is guestimate
 *      what will provide the most true positives (best statistical power)
 *      and fewest false positives.
 *
 *      In the case of adapter removal,
 *      it is also not usually important to remove every adapter, but only to
 *      minimize adapter contamination.  Failing to align a small percentage
 *      of sequences due to adapter contamination will not change the story
 *      told by the downstream analysis.  Nor will erroneously trimming off
 *      the 3' end of a small percentage of reads containing natural
 *      sequences resembling adapters.  Just trimming exact matches of
 *      the adapter sequence will generally remove 99% or more of the
 *      adapter contamination and minimize false-positives.  Tolerating
 *      1 or 2 differences has been shown to do slightly better overall.
 *      Modern read mapping software is also tolerant of adapter
 *      contamination and can clip adapters as needed.
 *
 *  Arguments:
 *      params      bl_align_t parameters.  Only min_match and
 *                  max_mismatch_percent are used.
 *      big         Sequence to be searched for matches to little
 *      little      Sequence to be located within big
 *
 *  Returns:
 *      Index of little sequence within big if found, index of null
 *      terminator of big otherwise
 *
 *  Examples:
 *      bl_param_t  params;
 *      bl_fastq_t  read;
 *      char        *adapter;
 *      size_t      index;
 *
 *      bl_align_set_min_match(&params, 3);
 *      bl_align_set_max_mismatch_percent(&params, 10);
 *      index = bl_align_map_seq_sub(&params,
 *          BL_FASTQ_SEQ(&read), BL_FASTQ_SEQ_LEN(&read),
 *          little, strlen(adapter)3, 10);
 *      if ( index != BL_FASTQ_SEQ_LEN(&read) )
 *          bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_align_map_seq_exact(3), bl_align_set_min_match(3),
 *      bl_align_set_max_mismatch_percent(3), bl_fastq_3p_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_align_map_seq_sub(const bl_align_t *params,
	    const char *big, size_t big_len,
	    const char *little, size_t little_len)

{
    // The strlen() looks expensive, but tests show that eliminating it
    // doesn't reduce run time measurably
    size_t      mismatch, max_mismatch,
		start, bc, lc,
		md, little_mm,
		min_match = params->min_match;
    
    // Start at 5' end assuming 5' littles already removed
    // Cutadapt uses a semiglobal alignment algorithm to find littles.
    // Not sure what the benefit of this is over exact matching. I would
    // assume that errors in little sequences are extremely rare.
    // https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm

    // Convert max mismatch percentage to a divisor for the string len
    md = 100 / params->max_mismatch_percent;
    little_mm = little_len / md;  // Max mismatch based on little len
    // Could stop at big_len - min_match, but the extra math
    // outweights the few iterations saved
    for (start = 0; start < big_len; ++start)
    {
	// Terminate loop as soon as max_mismatch is reached, before
	// checking other conditions
	max_mismatch = XT_MIN((big_len - start) / md, little_mm);
	for (bc = start, lc = 0, mismatch = 0;
	     (mismatch <= max_mismatch) &&
	     (lc < little_len) && (bc < big_len); ++bc, ++lc)
	{
	    if ( toupper(big[bc]) != little[lc] )
		++mismatch;
	}
	if ( mismatch <= max_mismatch )
	{
	    if ( lc - mismatch >= min_match )
		return start;
	}
    }
    return big_len;   // Location of '\0' terminator
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/align.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate the leftmost (farthest 5') match for sequence little within
 *      sequence big, using exact matching only.
 *
 *      The content of little is assumed to be all upper case.  This
 *      improves speed by avoiding numerous redundant toupper()
 *      conversions on the same string, assuming multiple big strings will
 *      be searched for little, as in adapter removal and read mapping.
 *      Use strlupper(3) or strupper(3) before calling this function if
 *      necessary.
 *
 *      A minimum of min_match bases must match between little and
 *      big.  This mainly matters near the end of big, where
 *      remaining bases are fewer than the length of little.
 *
 *      Note that alignment is not an exact science.  We cannot detect every
 *      true little sequence without falsely detecting other sequences, since
 *      it is impossible to know whether any given sequence is really from
 *      the source of interest (e.g. an adapter) or naturally
 *      occurring from another source.  The best we can do is guestimate
 *      what will provide the most true positives (best statistical power)
 *      and fewest false positives.
 *
 *      In the case of adapter removal,
 *      it is also not usually important to remove every adapter, but only to
 *      minimize adapter contamination.  Failing to align a small percentage
 *      of sequences due to adapter contamination will not change the story
 *      told by the downstream analysis.  Nor will erroneously trimming off
 *      the 3' end of a small percentage of reads containing natural
 *      sequences resembling adapters.  Just trimming exact matches of
 *      the adapter sequence will generally remove 99% or more of the
 *      adapter contamination and minimize false-positives.  Tolerating
 *      1 or 2 differences has been shown to do slightly better overall.
 *      Modern read mapping software is also tolerant of adapter
 *      contamination and can clip adapters as needed.
 *
 *  Arguments:
 *      params      bl_align_t parameters.  Only min_match is used.
 *      big         Sequence to be searched for matches to little
 *      little      Sequence to be located within big
 *
 *  Returns:
 *      Index of little sequence within big if found, index of null
 *      terminator of big otherwise
 *
 *  Examples:
 *      bl_param_t  params;
 *      bl_fastq_t  read;
 *      char        *adapter;
 *      size_t      index;
 *
 *      bl_align_set_min_match(&params, 3);
 *      index = bl_align_map_seq_exact(&params,
 *          BL_FASTQ_SEQ(&read), BL_FASTQ_SEQ_LEN(&read),
 *          little, strlen(adapter)3, 10);
 *      if ( index != BL_FASTQ_SEQ_LEN(&read) )
 *          bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_align_map_seq_sub(3), bl_align_set_min_match(3),
 *      bl_fastq_3p_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_align_map_seq_exact(const bl_align_t *params,
	    const char *big, size_t big_len,
	    const char *little, size_t little_len)

{
    size_t  start, bc, lc;
    
    // Start at 5' end assuming 5' adapters already removed
    // Cutadapt uses a semiglobal alignment algorithm to find adapters.
    // Not sure what the benefit of this is over exact matching. I would
    // assume that errors in adapter sequences are extremely rare.
    // https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm

    // Could stop at big_len - min_match, but the extra math
    // outweights the few iterations saved
    for (start = 0; start < big_len; ++start)
    {
	for (bc = start, lc = 0; (toupper(big[bc]) == little[lc]) &&
	     (lc < little_len); ++bc, ++lc)
	    ;
	if ( (lc == little_len) || ((bc == big_len) &&
	     (lc >= params->min_match)) )
	    return start;
    }
    return big_len;   // Location of '\0' terminator
}

