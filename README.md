# biolibc

## Status

The biolibc library is new and evolving rapidly.  It will likely undergo
occasional API changes as we discover suboptimal design
decisions that have made it into the implementation.

API changes will become less frequent as the software matures.
For now, be prepared to update dependent software in order to stay with
the latest biolibc.

Keep in mind that older releases are always available on Github, so API changes
won't impact you by surprise.  You only need to deal with them when
updating to a newer version of biolibc.

## Description

Many bioinformatics pipelines are essentially made of mud and straw;
quick-and-dirty single-use scripts written in various interpreted languages,
with little or no meaningful error reporting or documentation, and generally
poor performance.

Biolibc is a collection of high-quality bricks that can be used to build
efficient, robust software applications to replace disposable scripts.
Using biolibc, you can easily develop permanent solutions that are easy to
use and install, with near-optimal performance, so that no one ever need
reinvent that particular wheel.  Biolibc also facilitates development of
more complex applications by providing many commonly used building blocks,
thus releasing you from low-level coding.

For some examples, see
[biolibc-tools](https://github.com/auerlab/biolibc-tools),
[vcf-split](https://github.com/auerlab/vcf-split),
[ad2vcf](https://github.com/auerlab/ad2vcf),
[vcf2hap](https://github.com/auerlab/vcf2hap),
[haploh-vcf-depths](https://github.com/auerlab/haploh-vcf-depths),
[peak-classifier](https://github.com/auerlab/peak-classifier),
[generand](https://github.com/auerlab/generand),
[ad-matrix](https://github.com/auerlab/ad-matrix), and
[fastq-trim](https://github.com/auerlab/fastq-trim).

In more technical terms, biolibc is a library of fast, memory-efficient C
functions for processing biological data.  Like libc, it consists of numerous
disparate, general-purpose functions that can be used by a wide variety of
applications.

These include functions for reading and writing common file formats such as
BED, GFF, FASTA, FASTQ, SAM and VCF, string functions specific to
bioinformatics such as chromosome_name_cmp(), detecting feature overlaps,
etc.

Biolibc, and the more generic [libxtend](https://github.com/outpaddling/libxtend),
move your C coding to a higher level by providing the most basic
building blocks and additional software layers commonly needed in
bioinformatics programming.  Using C and biolibc, you can write simple,
near-optimal C programs that may be orders of magnitude faster than scripting
languages such as Perl, Python, and R, while requiring comparable coding effort.
[fastx-derep.c](https://github.com/auerlab/biolibc-tools/blob/main/fastx-derep.c),
a 1-page C program for dereplicating FASTA/FASTQ files, offers a great example
of how the right libraries eliminate the need for low-level programming in C.

See the
[Research Computing User's Guide](https://acadix.biz/publications.php)
[Compiled vs Interpreted Languages section](https://acadix.biz/RCUG/HTML/lang-selection.html#compiled-interpreted)
for a comparison of language performance.

## Summary list of functions

This summary lists currently available functions.
It is also available via "man biolibc" when biolibc is properly installed
via a package manager.  Each function listed has its own man page with
a more detailed description.

| Function | Purpose |
|--|--|
| bl_align_map_seq_exact(3)  |  Locate little sequence within big sequence |
| bl_align_map_seq_sub(3)  |  Locate little sequence within big sequence |
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
| bl_fastx_desc(3)  |  Return  description of a FASTX (FASTA or FASTQ) object |
| bl_fastx_desc_len(3)  |  Return length of a FASTX (FASTA or FASTQ) |
| bl_fastx_free(3)  |  Free memory for a FASTX (FASTA or FASTQ) object |
| bl_fastx_init(3)  |  Initialize a FASTX (FASTA or FASTQ) object |
| bl_fastx_plus(3)  |  Return '+' line of a FASTQ object, NULL if FASTA |
| bl_fastx_plus_len(3)  |  Return length of FASTQ '+' line, 0 if FASTA |
| bl_fastx_qual(3)  |  Return FASTQ quality line, NULL if FASTA |
| bl_fastx_qual_len(3)  |  Return length of FASTQ quality line, 0 if FASTA |
| bl_fastx_read(3)  |  Read FASTA or FASTQ record |
| bl_fastx_seq(3)  |  Return sequence of a FASTX (FASTA or FASTQ) object |
| bl_fastx_seq_len(3)  |  Return sequence length of a FASTX (FASTA or FASTQ) |
| bl_fastx_write(3)  |  Write FASTA or FASTQ record |
| bl_gff3_copy(3)  |  Copy a GFF3 object |
| bl_gff3_copy_header(3)  |  Read and copy a GFF3 header |
| bl_gff3_dup(3)  |  Duplicate a GFF3 object |
| bl_gff3_extract_attribute(3)  |  Extract GFF3 attribute by name |
| bl_gff3_free(3)  |  Free memory for a GFF3 object |
| bl_gff3_index_add(3)  |  Add a GFF3 feature to an in |
| bl_gff3_index_seek_reverse(3)  |  Search backward through GFF3 index |
| bl_gff3_init(3)  |  Initialize all fields in a GFF3 object |
| bl_gff3_read(3)  |  Read a GFF3 feature |
| bl_gff3_sam_cmp(3)  |  Compare positions of a SAM alignment and GFF3 |
| bl_gff3_sam_overlap(3)  |  Compute overlap between a SAM alignment and a GFF3 feature |
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
| bl_pos_list_from_csv(3)  |  Convert from comma |
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
| bl_sam_gff3_overlap(3)  |  Return the amount of overlap between a SAM alignment and a GFF3 feature |
| bl_sam_init(3)  |  Initialize all fields of a SAM object |
| bl_sam_read(3)  |  Read one SAM record |
| bl_sam_skip_header(3)  |  Read past SAM header |
| bl_sam_write(3)  |  Write a SAM object to a file stream |
| bl_vcf_call_downstream_of_alignment(3)  |  Return true if the location of a VCF call is downstream of an alignment |
| bl_vcf_call_in_alignment(3)  |  Return true if location of VCF call is |
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
## Design and Implementation

The code is organized following basic object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with accessor macros and mutator functions
provided, so dependent applications and libraries need not access
structure members directly.  Since the C language cannot enforce this, it's
up to application programmers to exercise self-discipline.

For detailed coding standards, see
https://github.com/outpaddling/Coding-Standards/.

## Building and installing

biolibc is intended to build cleanly in any POSIX environment on
any CPU architecture.  Please
don't hesitate to open an issue if you encounter problems on any
Unix-like system.

Primary development is done on FreeBSD with clang, but the code is frequently
tested on Linux, MacOS, NetBSD, and OpenIndiana as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [pkgsrc](https://pkgsrc.org/), etc.
End users should install via one of these if at all possible.

I maintain a FreeBSD port and a pkgsrc package, which is sufficient to install
cleanly on virtually any POSIX platform.  If you would like to see a
biolibc package in another package manager, please consider creating a package
yourself.  This will be one of the easiest packages in the collection and
hence a good vehicle to learn how to create packages.

For an overview of available package managers, see the
[Repology website](https://repology.org/).

### Installing biolibc on FreeBSD:

FreeBSD is a highly underrated platform for scientific computing, with over
2,000 scientific libraries and applications in the FreeBSD ports collection
(of more than 30,000 total), modern clang compiler, fully-integrated ZFS
file system, and renowned security, performance, and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like [Ubuntu](https://ubuntu.com/).
However, if you're a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).
[GhostBSD](https://ghostbsd.org/) offers and experience very similar
to Ubuntu, but is build on FreeBSD rather than Debian Linux.  GhostBSD
packages lag behind FreeBSD by a few months, but this is not generally
an issue and there are workarounds.

To install the binary package on FreeBSD:

```
pkg install biolibc
```

You can just as easily build and install from source.  This is useful for
FreeBSD ports with special build options, for building with non-portable
optimizations such as -march=native, and for 
[work-in-progress ports](https://github.com/outpaddling/freebsd-ports-wip),
for which binary packages are not yet maintained.

```
cd /usr/ports/biology/biolibc && env CFLAGS='-march=native -O2' make install
cd /usr/ports/wip/biolibc && make install
```

### Installing via pkgsrc

pkgsrc is a cross-platform package manager that works on any Unix-like
platform. It is native to [NetBSD](https://www.netbsd.org/) and well-supported
on [Illumos](https://illumos.org/), [MacOS](https://www.apple.com/macos/),
[RHEL](https://www.redhat.com)/[CentOS](https://www.centos.org/), and
many other Linux distributions.
Using pkgsrc does not require admin privileges.  You can install a pkgsrc
tree in any directory to which you have write access and easily install any
of the nearly 20,000 packages in the collection.

The
[auto-pkgsrc-setup](https://github.com/outpaddling/auto-admin/blob/master/User-scripts/auto-pkgsrc-setup)
script will help you install pkgsrc in about 10 minutes.  Just download it
and run

```
sh auto-pkgsrc-setup
```

Then, assuming you selected current packages and the default prefix

```
source ~/Pkgsrc/pkg/etc/pkgsrc.sh   # Or pkgsrc.csh for csh or tcsh
cd ~/Pkgsrc/biology/biolibc
sbmake install clean clean-depends
```

See the pkgsrc documentation for more information.

Community support for pkgsrc is available through the
[pkgsrc-users](http://netbsd.org/mailinglists) mailing list.

### Other package managers

Packages for libxtend are known to exist in the following package managers.
These are maintained by third parties and not directly supported here.

[BioArchLinux](https://github.com/BioArchLinux/Packages)

### Building biolibc locally

Below are caveman install instructions for development purposes, not
recommended for regular use.

biolibc depends on [libxtend](https://github.com/outpaddling/libxtend).
Install libxtend before attempting to build biolibc.

1. Clone the repository
2. Run "make depend" to update Makefile.depend
3. Run "make install"

The default install prefix is ../local.  Clone biolibc, libxtend and dependent
apps into sibling directories so that ../local represents a common path to all
of them.

To facilitate incorporation into package managers, the Makefile respects
standard make/environment variables such as CC, CFLAGS, PREFIX, etc.  

The library, headers, and man pages are installed under
`${DESTDIR}${PREFIX}`.  DESTDIR is empty by default and is primarily used by
package managers to stage installations.  PREFIX defaults to ../local.

Add-on libraries required for the build, such as biolibc, should be found
under either `${PREFIX}` or `${LOCALBASE}`, which defaults to `${PREFIX}`.
LOCALBASE can be set independently if you want to use libraries installed
by FreeBSD ports (/usr/local), MacPorts (/opt/local), pkgsrc (/usr/pkg), etc.

To install directly to /myprefix, assuming biolibc is installed there as well,
using a make variable:

```
make PREFIX=/myprefix clean depend install
```

Using an environment variable:

```
# C-shell and derivatives
setenv PREFIX /myprefix
make clean depend install

# Bourne shell and derivatives
PREFIX=/myprefix
export PREFIX
make clean depend install
```

View the Makefile for full details.
