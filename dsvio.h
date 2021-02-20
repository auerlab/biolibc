#ifndef __dsvio_h__
#define __dsvio_h__

/* dsvio.c */
int dsv_read_field(FILE *infile, char buff[], size_t buff_size,
		   int delim, size_t *len);
int dsv_skip_field(FILE *infile, int delim);
int dsv_skip_rest_of_line(FILE *infile);
int tsv_read_field(FILE *infile, char buff[], size_t buff_size,
		   size_t *len);
int tsv_skip_field(FILE *infile);
int tsv_skip_rest_of_line(FILE *infile);
int csv_read_field(FILE *infile, char buff[], size_t buff_size,
		   size_t *len);
int csv_skip_field(FILE *infile);
int csv_skip_rest_of_line(FILE *infile);

#endif  // __dsvio_h__
