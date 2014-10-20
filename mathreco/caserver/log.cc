#include <stdio.h>
#include <stdarg.h>
#include <time.h>

static FILE *fp = 0;

void
log_init() {
	fp = fopen("/tmp/caslog", "a");
	if (fp) fprintf(fp, "-------------------------------------------------\n");
}

FILE *
logfp() {
	return fp;
}

void
logtime() {
	if (fp) {
		time_t rawtime;
		struct tm *curtime;
		time(&rawtime);
		curtime = localtime(&rawtime);
		fprintf(fp, "Timestamp: %d:%.2d:%.2d - %.d/%.d/%d\n", curtime->tm_hour, curtime->tm_min, curtime->tm_sec, curtime->tm_mon+1, curtime->tm_mday, curtime->tm_year+1900);
	}
}

void
logmsg(const char *s, ...) {
	if (fp) {
		va_list va;
		va_start(va, s);
		vfprintf(fp, s, va);
		va_end(va);
		fflush(fp);
	}
}
