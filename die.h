#ifndef DIE_H
#define DIE_H

#include <errno.h>
#include <string.h>

void die(const char *fmt, ...);

#define die_errno(msg)          die("%s: %s", msg, strerror(errno))

#endif
