
/*
 *  Generated by /usr/local/bin/auto-gen-get-set
 *
 *  Mutator functions for setting with no sanity checking.  Use these to
 *  set structure members from functions outside the bl_fastx_t
 *  class.  These macros perform no data validation.  Hence, they achieve
 *  maximum performance where data are guaranteed correct by other means.
 *  Use the mutator functions (same name as the macro, but lower case)
 *  for more robust code with a small performance penalty.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
 */

/* temp-fastx-mutators.c */
int bl_fastx_set_format(bl_fastx_t *bl_fastx_ptr, int new_format);
int bl_fastx_set_fasta(bl_fastx_t *bl_fastx_ptr, bl_fasta_t new_fasta);
int bl_fastx_set_fastq(bl_fastx_t *bl_fastx_ptr, bl_fastq_t new_fastq);
