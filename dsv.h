#ifndef _dsv_h_
#define _dsv_h_

#define DSV_INIT                { 0, 0, NULL, NULL }
#define DSV_FIELD_MAX_CHARS     32767

#define DSV_FIELD(line,col)     ((line)->fields[col-1]) // 1-based column

typedef struct
{
    size_t      array_size,
		num_fields;
    char        **fields,
		*delims;
}   dsv_line_t;

/* dsv.c */
int dsv_read_line(FILE *infile, dsv_line_t *dsv_line, const char *delims);
void dsv_write_line(FILE *outfile, dsv_line_t *dsv_line);
void dsv_free_line(dsv_line_t *dsv_line);
void dsv_copy_line(dsv_line_t *dest, dsv_line_t *src);
int dsv_read_field(FILE *infile, char buff[], size_t buff_size,
		   const char *delims, size_t *len);
int dsv_skip_field(FILE *infile, const char *delims);
int dsv_skip_rest_of_line(FILE *infile);
int tsv_read_field(FILE *infile, char buff[], size_t buff_size,
		   size_t *len);
int tsv_skip_field(FILE *infile);
int tsv_skip_rest_of_line(FILE *infile);
int csv_read_field(FILE *infile, char buff[], size_t buff_size,
		   size_t *len);
int csv_skip_field(FILE *infile);
int csv_skip_rest_of_line(FILE *infile);

#endif  // _dsv_h_
