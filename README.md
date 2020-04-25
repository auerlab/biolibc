# biolibc

Purpose
=======

Biolibc is a library of fast, memory-efficient, low-level functions for
processing biological data.

Like libc, it consists of numerous disparate, general-purpose functions which
could be used by a wide variety of applications.

Design and Implementation
=========================

The code is organized following object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with access functions or macros provided,
to applications should not access data members directly.  Since the C language
cannot enforce this, it's up to application programmers to exercise
self-discipline.

Building and installing
=======================

Set PREFIX to the prefix where you would like to install.  Default is ../local.

Then simply run

```sh
make install
```
