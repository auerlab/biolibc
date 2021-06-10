# biolibc

Purpose
=======

Biolibc is a library of fast, memory-efficient, low-level functions for
processing biological data.

Like libc, it consists of numerous disparate, general-purpose functions which
could be used by a wide variety of applications.

These include functions for streaming common file formats such as SAM and
VCF, string functions specific to bioinformatics, etc.

Design and Implementation
=========================

The code is organized following basic object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with accessor and mutator functions
(or macros) provided, so dependent applications and libraries need not access
structure members directly.  Since the C language cannot enforce this, it's
up to application programmers to exercise self-discipline.

Building and installing
=======================

biolibc depends on [libxtend](https://github.com/outpaddling/libxtend).
Install libxtend before attempting to build biolibc.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [pkgsrc](https://pkgsrc.org/), etc.
End users should install via one of these if at all possible.

I maintain a FreeBSD port and a pkgsrc package.

FreeBSD is a highly underrated platform for scientific computing, with over
1,800 scientific libraries and applications in the FreeBSD ports collection,
(of more than 30,000 total),
fully-integrated ZFS, and renowned performance and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like Ubuntu.
However, if you are a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).

To install biolibc on FreeBSD:

```
pkg install biolibc
```

pkgsrc is a cross-platform package manager that works on any Unix-like
platform. It is native to NetBSD and well-supported on Illumos, Linux, and
MacOS.  Using pkgsrc does not require admin privileges.  The
[auto-pkgsrc-setup](http://netbsd.org/~bacon/) script
can assist you with basic setup.

To install via pkgsrc, first bootstrap pkgsrc using auto-pkgsrc-setup or any
other means.  Then run the following commands:

```
cd pkgsrc-dir/sysutils/auto-admin
bmake install clean
cd pkgsrc-dir/wip/biolibc
bmake install clean
```

To build biolibc locally (for development purposes, not recommended for
regular use):

1. Clone the repository
2. Run "make depend" to update Makefile.depend
3. Run "make install"

The default install prefix is ../local.  Clone biolibc, libxtend and dependent
apps into sibling directories so that ../local represents a common path to all
of them.

To facilitate incorporation into package managers, the Makefile respects
standard make/environment variables such as CC, CFLAGS, PREFIX, etc.  

Add-on libraries required for the build, such as libxtend, should be found
under ${LOCABASE}, which defaults to ../local.
The library, headers, and man pages are installed under
${DESTDIR}${PREFIX}.  DESTDIR is empty by default and is primarily used by
package managers to stage installations.  PREFIX defaults to ${LOCALBASE}.

To install directly to /myprefix, assuming libxtend is installed there as well,
using a make variable:

```
make LOCALBASE=/myprefix clean depend install
```

Using an environment variable:

```
# C-shell and derivatives
setenv LOCALBASE /myprefix
make clean depend install

# Bourne shell and derivatives
LOCALBASE=/myprefix
export LOCALBASE
make clean depend install
```

View the Makefile for full details.
