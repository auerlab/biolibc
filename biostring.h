#ifndef __biostring_h__
#define __biostring_h__

#ifndef __biolibc_h__
#include "biolibc.h"
#endif

/* chromosome-name-cmp.c */
int chromosome_name_cmp(const char *n1, const char *n2);

/* strptrcmp.c */
int strptrcmp(const char **p1, const char **p2);

// CentOS 7 gcc does not support restrict, which helps the optimizer produce
// faster code.  Keep _RESTRICT def separate from strlcpy() prototype in case
// other platforms are missing one but not the other.
#ifdef __linux__
#define _RESTRICT
#else
#define _RESTRICT   restrict
#endif

#ifdef __linux__
// size_t strlcpy(char *dest, const char *src, size_t len);
#define strlcpy(dest,src,len)   strcpy(dest,src)
#endif

#endif // __biostring_h__
