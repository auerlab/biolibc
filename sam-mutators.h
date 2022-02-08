
/*
 *  Generated by /usr/local/bin/auto-gen-get-set
 *
 *  Mutator functions for setting with no sanity checking.  Use these to
 *  set structure members from functions outside the bl_sam_t
 *  class.  These macros perform no data validation.  Hence, they achieve
 *  maximum performance where data are guaranteed correct by other means.
 *  Use the mutator functions (same name as the macro, but lower case)
 *  for more robust code with a small performance penalty.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
 */

/* temp-sam-mutators.c */
int bl_sam_set_qname_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_qname_element);
int bl_sam_set_qname_cpy(bl_sam_t *bl_sam_ptr, char new_qname[], size_t array_size);
int bl_sam_set_flag(bl_sam_t *bl_sam_ptr, unsigned new_flag);
int bl_sam_set_rname_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_rname_element);
int bl_sam_set_rname_cpy(bl_sam_t *bl_sam_ptr, char new_rname[], size_t array_size);
int bl_sam_set_pos(bl_sam_t *bl_sam_ptr, uint64_t new_pos);
int bl_sam_set_mapq(bl_sam_t *bl_sam_ptr, unsigned char new_mapq);
int bl_sam_set_cigar_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_cigar_element);
int bl_sam_set_cigar_cpy(bl_sam_t *bl_sam_ptr, char new_cigar[], size_t array_size);
int bl_sam_set_rnext_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_rnext_element);
int bl_sam_set_rnext_cpy(bl_sam_t *bl_sam_ptr, char new_rnext[], size_t array_size);
int bl_sam_set_pnext(bl_sam_t *bl_sam_ptr, uint64_t new_pnext);
int bl_sam_set_tlen(bl_sam_t *bl_sam_ptr, long new_tlen);
int bl_sam_set_seq(bl_sam_t *bl_sam_ptr, char *new_seq);
int bl_sam_set_seq_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_seq_element);
int bl_sam_set_seq_cpy(bl_sam_t *bl_sam_ptr, char *new_seq, size_t array_size);
int bl_sam_set_qual(bl_sam_t *bl_sam_ptr, char *new_qual);
int bl_sam_set_qual_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_qual_element);
int bl_sam_set_qual_cpy(bl_sam_t *bl_sam_ptr, char *new_qual, size_t array_size);
int bl_sam_set_seq_len(bl_sam_t *bl_sam_ptr, size_t new_seq_len);
int bl_sam_set_qual_len(bl_sam_t *bl_sam_ptr, size_t new_qual_len);
