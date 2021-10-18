#ifndef _BIOLIBC_TRANSLATE_H
#define _BIOLIBC_TRANSLATE_H

unsigned long   next_start_codon(FILE *rna_stream);
unsigned long   next_stop_codon(FILE *rna_stream, char codon[4]);

#endif  // _BIOLIBC_TRANSLATE_H
