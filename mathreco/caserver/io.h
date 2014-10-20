#ifndef IO_H_
#define IO_H_

#include "expr.h"

int fullread(int fd, void *buf, size_t n);
int fullwrite(int fd, void *buf, size_t n);
int readexpr(int fd, exprtree **tree);

#endif
