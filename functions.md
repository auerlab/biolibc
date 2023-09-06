# Biolibc function list

| Function | Purpose |
|----------|---------|
| BIOLIBC  |  Library of high |
| bl_align_map_seq_exact(3)  |  Locate little sequence in big sequence |
| bl_align_map_seq_sub(3)  |  Locate little sequence in big sequence |
| bl_bed_check_order(3)  |  Compare positions of two bed records |
| bl_bed_gff3_cmp(3)  |  Compare positions of BED and GFF3 objects |
| bl_bed_read(3)  |  Read a BED record |
| bl_bed_skip_header(3)  |  Read past BED header |
| bl_bed_write(3)  |  Write a BED record |
| bl_chrom_name_cmp(3)  |  Compare chromosome names numerically or lexically |
| bl_fasta_free(3)  |  Free memory for a FASTA object |
| bl_fasta_init(3)  |  Initialize all fields of a FASTA object |
| bl_fasta_read(3)  |  Read a FASTA record |
| bl_fasta_write(3)  |  Write a FASTA object |
| bl_fastq_3p_trim(3)  |  Trim 3' end of a FASTQ object |
| bl_fastq_find_3p_low_qual(3)  |  Find start of low |
| bl_fastq_free(3)  |  Free memory for a FASTQ object |
| bl_fastq_init(3)  |  Initialize all fields in a FASTQ object |
| bl_fastq_name_cmp(3)  |  Compare read names of two FASTQ objects |
| bl_fastq_read(3)  |  Read a FASTQ record |
| bl_fastq_write(3)  |  Write a FASTQ record |
| bl_fastx_desc(3)  |  Return  description of a FASTX object |
| bl_fastx_desc_len(3)  |  Return length of FASTX description |
| bl_fastx_free(3)  |  Free memory for a FASTX object |
| bl_fastx_init(3)  |  Initialize a FASTX object |
| bl_fastx_plus(3)  |  Return '+' line of a FASTQ object, NULL if FASTA |
| bl_fastx_plus_len(3)  |  Return length of FASTQ '+' line, 0 if FASTA |
| bl_fastx_qual(3)  |  Return FASTQ quality line, NULL if FASTA |
| bl_fastx_qual_len(3)  |  Return length of FASTQ quality line, 0 if FASTA |
| bl_fastx_read(3)  |  Read a FASTX record |
| bl_fastx_seq(3)  |  Return sequence of a FASTX object |
| bl_fastx_seq_len(3)  |  Return length of a FASTX sequence object |
| bl_fastx_write(3)  |  Write a FASTX record |
| bl_gff3_copy(3)  |  Copy a GFF3 object |
| bl_gff3_copy_header(3)  |  Read and copy a GFF3 header |
| bl_gff3_dup(3)  |  Duplicate a GFF3 object |
| bl_gff3_extract_attribute(3)  |  Extract GFF3 attribute by name |
| bl_gff3_free(3)  |  Free memory for a GFF3 object |
| bl_gff3_index_add(3)  |  Add a GFF3 feature to an in |
| bl_gff3_index_seek_reverse(3)  |  Search backward through GFF3 index |
| bl_gff3_init(3)  |  Initialize all fields in a GFF3 object |
| bl_gff3_read(3)  |  Read a GFF3 feature |
| bl_gff3_sam_cmp(3)  |  Compare SAM/GFF3 positions |
| bl_gff3_sam_overlap(3)  |  Compute SAM/GFF3 overlap |
| bl_gff3_skip_header(3)  |  Read past header in a GFF3 file |
| bl_gff3_to_bed(3)  |  Convert a GFF3 featuer to a BED object |
| bl_gff3_write(3)  |  Write a GFF3 feature |
| bl_next_start_codon(3)  |  Find next start codon |
| bl_next_stop_codon(3)  |  Find next stop codon |
| bl_overlap_print(3)  |  Print overlap summary for two features |
| bl_overlap_set_all(3)  |  Set overlap fields for two features |
| bl_pos_list_add_position(3)  |  Add a position to a list |
| bl_pos_list_allocate(3)  |  Initialize position list object |
| bl_pos_list_free(3)  |  Free a position list object |
| bl_pos_list_from_csv(3)  |  Convert comma |
| bl_pos_list_sort(3)  |  Sort a position list |
| bl_sam_buff_add_alignment(3)  |  Add alignment to SAM buffer |
| bl_sam_buff_alignment_ok(3)  |  Verify alignment quality |
| bl_sam_buff_check_order(3)  |  Check sort order of SAM records |
| bl_sam_buff_free_alignment(3)  |  Free an alignment in a SAM buffer |
| bl_sam_buff_init(3)  |  Initialize a SAM buffer object |
| bl_sam_buff_out_of_order(3)  |  Print sort order message and exit |
| bl_sam_buff_shift(3)  |  Close gap after removing a SAM alignment |
| bl_sam_copy(3)  |  Copy a SAM object |
| bl_sam_copy_header(3)  |  Copy SAM header to another stream |
| bl_sam_fclose(3)  |  Close a stream opened by bl_sam_fopen(3) |
| bl_sam_fopen(3)  |  Open a SAM/BAM/CRAM file |
| bl_sam_free(3)  |  Destroy a SAM object |
| bl_sam_gff3_cmp(3)  |  Compare positions of SAM and GFF3 records |
| bl_sam_gff3_overlap(3)  |  Compute SAM/GFF3 overlap |
| bl_sam_init(3)  |  Initialize all fields of a SAM object |
| bl_sam_read(3)  |  Read one SAM record |
| bl_sam_skip_header(3)  |  Read past SAM header |
| bl_sam_write(3)  |  Write a SAM object to a file stream |
| bl_vcf_call_downstream_of_alignment(3)  |  Return true if VCF call is downstream of alignment |
| bl_vcf_call_in_alignment(3)  |  Return true if VCF call is within alignment |
| bl_vcf_call_out_of_order(3)  |  Terminate with VCF sort error message |
| bl_vcf_free(3)  |  Destroy a VCF object |
| bl_vcf_get_sample_ids(3)  |  Extract sample IDs from a VCF header |
| bl_vcf_init(3)  |  Initialize fields in a VCF object |
| bl_vcf_parse_field_spec(3)  |  Convert comma |
| bl_vcf_read_ss_call(3)  |  Read a single |
| bl_vcf_read_static_fields(3)  |  Read static VCF fields |
| bl_vcf_skip_header(3)  |  Read past VCF header |
| bl_vcf_skip_meta_data(3)  |  Read past VCF metadata |
| bl_vcf_write_ss_call(3)  |  Write a single |
| bl_vcf_write_static_fields(3)  |  Write VCF static fields |
