#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>

void log_init();
void logmsg(const char *s, ...);
void logtime();
FILE *logfp();

#endif 
