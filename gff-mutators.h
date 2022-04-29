
/*
 *  Generated by /usr/local/bin/auto-gen-get-set
 *
 *  Mutator functions for setting with no sanity checking.  Use these to
 *  set structure members from functions outside the bl_gff_t
 *  class.  These macros perform no data validation.  Hence, they achieve
 *  maximum performance where data are guaranteed correct by other means.
 *  Use the mutator functions (same name as the macro, but lower case)
 *  for more robust code with a small performance penalty.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
 */

/* temp-gff-mutators.c */
int bl_gff_set_seqid_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_seqid_element);
int bl_gff_set_seqid_cpy(bl_gff_t *bl_gff_ptr, char new_seqid[], size_t array_size);
int bl_gff_set_source_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_source_element);
int bl_gff_set_source_cpy(bl_gff_t *bl_gff_ptr, char new_source[], size_t array_size);
int bl_gff_set_type_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_type_element);
int bl_gff_set_type_cpy(bl_gff_t *bl_gff_ptr, char new_type[], size_t array_size);
int bl_gff_set_start(bl_gff_t *bl_gff_ptr, int64_t new_start);
int bl_gff_set_end(bl_gff_t *bl_gff_ptr, int64_t new_end);
int bl_gff_set_score(bl_gff_t *bl_gff_ptr, double new_score);
int bl_gff_set_strand(bl_gff_t *bl_gff_ptr, char new_strand);
int bl_gff_set_phase(bl_gff_t *bl_gff_ptr, char new_phase);
int bl_gff_set_attributes(bl_gff_t *bl_gff_ptr, char *new_attributes);
int bl_gff_set_attributes_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_attributes_element);
int bl_gff_set_attributes_cpy(bl_gff_t *bl_gff_ptr, char *new_attributes, size_t array_size);
int bl_gff_set_attributes_array_size(bl_gff_t *bl_gff_ptr, size_t new_attributes_array_size);
int bl_gff_set_attributes_len(bl_gff_t *bl_gff_ptr, size_t new_attributes_len);
int bl_gff_set_feature_id(bl_gff_t *bl_gff_ptr, char *new_feature_id);
int bl_gff_set_feature_id_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_feature_id_element);
int bl_gff_set_feature_id_cpy(bl_gff_t *bl_gff_ptr, char *new_feature_id, size_t array_size);
int bl_gff_set_feature_name(bl_gff_t *bl_gff_ptr, char *new_feature_name);
int bl_gff_set_feature_name_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_feature_name_element);
int bl_gff_set_feature_name_cpy(bl_gff_t *bl_gff_ptr, char *new_feature_name, size_t array_size);
int bl_gff_set_feature_parent(bl_gff_t *bl_gff_ptr, char *new_feature_parent);
int bl_gff_set_feature_parent_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_feature_parent_element);
int bl_gff_set_feature_parent_cpy(bl_gff_t *bl_gff_ptr, char *new_feature_parent, size_t array_size);
int bl_gff_set_file_pos(bl_gff_t *bl_gff_ptr, long new_file_pos);
