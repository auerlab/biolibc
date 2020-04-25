#ifndef __tsvio_h__
#define __tsvio_h__

/* tsvio.c */
int tsv_read_field(FILE *infile, char buff[], size_t buff_size, size_t *len);
int tsv_skip_field(FILE *infile);
int tsv_skip_rest_of_line(FILE *infile);

#endif  // __tsvio_h__
