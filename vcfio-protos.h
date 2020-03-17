/* vcfio.c */
void vcf_skip_header(const char *argv[], FILE *vcf_stream);
void vcf_get_sample_ids(const char *argv[], FILE *vcf_stream, char *sample_ids[], size_t first_col, size_t last_col);
int vcf_read_static_fields(const char *argv[], FILE *vcf_stream, vcf_call_t *vcf_call);
int vcf_read_ss_call(const char *argv[], FILE *vcf_stream, vcf_call_t *vcf_call);
int vcf_write_static_fields(const char *argv[], FILE *vcf_stream, vcf_call_t *vcf_call);
int vcf_write_ss_call(const char *argv[], FILE *vcf_stream, vcf_call_t *vcf_call);
size_t vcf_read_calls_for_position(const char *argv[], FILE *vcf_stream, vcf_calls_for_position_t *vcf_calls_for_position);
char **vcf_sample_alloc(vcf_call_t *vcf_call, size_t samples);
