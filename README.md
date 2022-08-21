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

### Building biolibc locally

Below are cave man install instructions for development purposes, not
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
