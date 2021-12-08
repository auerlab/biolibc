#ifndef _BIOLIBC_TRANSLATE_H
#define _BIOLIBC_TRANSLATE_H

long   bl_next_start_codon(FILE *rna_stream, char codon[4]);
long   bl_next_stop_codon(FILE *rna_stream, char codon[4]);

#endif  // _BIOLIBC_TRANSLATE_H
