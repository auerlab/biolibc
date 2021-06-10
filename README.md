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

The code is organized following object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with access functions or macros provided,
so applications should not access data members directly.  Since the C language
cannot enforce this, it's up to application programmers to exercise
self-discipline.

Building and installing
=======================

biolibc depends on [libxtend](https://github.com/outpaddling/libxtend).
Install libxtend before attempting to build biolibc.

The Makefile is designed to be friendly to package managers, such as
Debian packages, [FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [pkgsrc](https://pkgsrc.org/), etc.
End users should install via one of these if at all possible.

To build locally for development purposes:

1. Clone the repository
2. Run "make depend" to update Makefile.depend
3. Run "make install"

The default install prefix is ../local.  Clone libxtend and dependent apps
into sibling directories so that ../local represents a common path to all of
them.

To facilitate easy packaging, the Makefile respects standard make/environment
variables such as CC, CFLAGS, PREFIX, etc.  For example, to install to
/myprefix:

```
make PREFIX=/myprefix install
```

If libxtend is installed under /myprefix:

```
make LOCALBASE=/myprefix depend
```

View the Makefile for full details.
