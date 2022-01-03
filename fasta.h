
#ifndef _biolibc_h_
#include "biolibc.h"
#endif

typedef struct
{
    char    *desc;
    char    *seq;
    size_t  desc_array_size,
	    seq_array_size,
	    desc_len,
	    seq_len;
}   bl_fasta_t;

#define BL_FASTA_INIT           { NULL, NULL, 0, 0, 0, 0 }
#define BL_FASTA_LINE_UNLIMITED 0

/*
 *  Generated by /home/bacon/scripts/gen-get-set
 *
 *  Accessor macros.  Use these to access structure members from functions
 *  outside the bl_fasta_t class.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
  */

#define BL_FASTA_DESC(ptr)              ((ptr)->desc)
#define BL_FASTA_DESC_AE(ptr,c)         ((ptr)->desc[c])
#define BL_FASTA_SEQ(ptr)               ((ptr)->seq)
#define BL_FASTA_SEQ_AE(ptr,c)          ((ptr)->seq[c])
#define BL_FASTA_DESC_ARRAY_SIZE(ptr)   ((ptr)->desc_array_size)
#define BL_FASTA_SEQ_ARRAY_SIZE(ptr)    ((ptr)->seq_array_size)
#define BL_FASTA_DESC_LEN(ptr)          ((ptr)->desc_len)
#define BL_FASTA_SEQ_LEN(ptr)           ((ptr)->seq_len)

/* fasta.c */
int bl_fasta_read(FILE *fasta_stream, bl_fasta_t *record);
int bl_fasta_write(FILE *fasta_stream, bl_fasta_t *record, size_t chars_per_line);
void bl_fasta_free(bl_fasta_t *record);
void bl_fasta_init(bl_fasta_t *record);

/* fasta-mutators.c */
int bl_fasta_set_desc(bl_fasta_t *bl_fasta_ptr, char *new_desc);
int bl_fasta_set_desc_ae(bl_fasta_t *bl_fasta_ptr, size_t c, char new_desc_element);
int bl_fasta_set_desc_cpy(bl_fasta_t *bl_fasta_ptr, char *new_desc, size_t array_size);
int bl_fasta_set_seq(bl_fasta_t *bl_fasta_ptr, char *new_seq);
int bl_fasta_set_seq_ae(bl_fasta_t *bl_fasta_ptr, size_t c, char new_seq_element);
int bl_fasta_set_seq_cpy(bl_fasta_t *bl_fasta_ptr, char *new_seq, size_t array_size);
int bl_fasta_set_desc_array_size(bl_fasta_t *bl_fasta_ptr, size_t new_desc_array_size);
int bl_fasta_set_seq_array_size(bl_fasta_t *bl_fasta_ptr, size_t new_seq_array_size);
int bl_fasta_set_desc_len(bl_fasta_t *bl_fasta_ptr, size_t new_desc_len);
int bl_fasta_set_seq_len(bl_fasta_t *bl_fasta_ptr, size_t new_seq_len);
