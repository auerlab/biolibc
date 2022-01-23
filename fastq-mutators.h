
/*
 *  Generated by /usr/local/bin/auto-gen-get-set
 *
 *  Mutator macros for setting with no sanity checking.  Use these to
 *  set structure members from functions outside the bl_fastq_t
 *  class.  These macros perform no data validation.  Hence, they achieve
 *  maximum performance where data are guaranteed correct by other means.
 *  Use the mutator functions (same name as the macro, but lower case)
 *  for more robust code with a small performance penalty.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
 */

/* temp-fastq-mutators.c */
int bl_fastq_set_desc(bl_fastq_t *bl_fastq_ptr, char *new_desc);
int bl_fastq_set_desc_ae(bl_fastq_t *bl_fastq_ptr, size_t c, char new_desc_element);
int bl_fastq_set_desc_cpy(bl_fastq_t *bl_fastq_ptr, char *new_desc, size_t array_size);
int bl_fastq_set_seq(bl_fastq_t *bl_fastq_ptr, char *new_seq);
int bl_fastq_set_seq_ae(bl_fastq_t *bl_fastq_ptr, size_t c, char new_seq_element);
int bl_fastq_set_seq_cpy(bl_fastq_t *bl_fastq_ptr, char *new_seq, size_t array_size);
int bl_fastq_set_plus(bl_fastq_t *bl_fastq_ptr, char *new_plus);
int bl_fastq_set_plus_ae(bl_fastq_t *bl_fastq_ptr, size_t c, char new_plus_element);
int bl_fastq_set_plus_cpy(bl_fastq_t *bl_fastq_ptr, char *new_plus, size_t array_size);
int bl_fastq_set_qual(bl_fastq_t *bl_fastq_ptr, char *new_qual);
int bl_fastq_set_qual_ae(bl_fastq_t *bl_fastq_ptr, size_t c, char new_qual_element);
int bl_fastq_set_qual_cpy(bl_fastq_t *bl_fastq_ptr, char *new_qual, size_t array_size);
int bl_fastq_set_desc_array_size(bl_fastq_t *bl_fastq_ptr, size_t new_desc_array_size);
int bl_fastq_set_seq_array_size(bl_fastq_t *bl_fastq_ptr, size_t new_seq_array_size);
int bl_fastq_set_plus_array_size(bl_fastq_t *bl_fastq_ptr, size_t new_plus_array_size);
int bl_fastq_set_qual_array_size(bl_fastq_t *bl_fastq_ptr, size_t new_qual_array_size);
int bl_fastq_set_desc_len(bl_fastq_t *bl_fastq_ptr, size_t new_desc_len);
int bl_fastq_set_seq_len(bl_fastq_t *bl_fastq_ptr, size_t new_seq_len);
int bl_fastq_set_plus_len(bl_fastq_t *bl_fastq_ptr, size_t new_plus_len);
int bl_fastq_set_qual_len(bl_fastq_t *bl_fastq_ptr, size_t new_qual_len);