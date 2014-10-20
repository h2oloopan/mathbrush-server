#include "io.h"
#include "cmdcode.h"
#include "expr.h"
#include "log.h"
#ifdef WIN32
#define NOMINMAX
#include <winsock.h>
#include <io.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#endif
#include <errno.h>
#include "int32.h"

#define min(a,b) (((a) < (b)) ? (a) : (b))

int
fullread(int fd, void *buf, size_t n) {
	size_t nr = 0;
	int e;
	while (nr < n) {
		e = recv(fd, ((char *)buf + nr), n - nr, 0);
		if (e <= 0) {
			if (e == EINTR || e == EAGAIN 
#ifndef WIN32
			 || e == EWOULDBLOCK
#endif
			 || (e == -1 && (errno == EINTR || errno == EAGAIN
#ifndef WIN32
			 || errno == EWOULDBLOCK
#endif
			 ))) {
				continue;
			}
			else {
				return CASCMD_HUP;
			}
		}
		nr += e;
	}
	return nr;
}

int
fullwrite(int fd, void *buf, size_t n) {
	size_t nr = 0;
	int e;
	while (nr < n) {
		e = send(fd, ((char *)buf + nr), n - nr, 0);
		if (e <= 0) {
			if (e == EINTR || e == EAGAIN
#ifndef WIN32
			 || e == EWOULDBLOCK
#endif
			 || (e == -1 && (errno == EINTR || errno == EAGAIN
#ifndef WIN32
			 || errno == EWOULDBLOCK
#endif
			 ))) {
				continue;
			}
			else {
				return CASCMD_HUP;
			}
		}
		nr += e;
	}
	return nr;
}

int
readexpr(int fd, exprtree **tree) {
#define CMDBUFLEN 16384
	char buf[CMDBUFLEN];
	uint32_t mathlen;
	int e;
	int32_t rep = CASREP_OK;

	e = fullread(fd, (void *)&mathlen, sizeof(mathlen));
	if (e != sizeof(mathlen)) {
		return e;
	}
	
	if (mathlen > CMDBUFLEN) { // or just hup()?
		while (mathlen > 0) {
			size_t readamt = min(sizeof(buf), mathlen);
			e = fullread(fd, (void *)buf, readamt);
			if (e < 0) {
				rep = CASREP_SYSERR;
				break;
			}
			if (e != readamt) {
				rep = CASREP_MSGERR;
				break;
			}
			mathlen -= e;
		}
		if (mathlen == 0) {
			rep = CASREP_TOOLARGE;
		}
	}
	else {
		e = fullread(fd, (void *)buf, mathlen);
		if (e < 0) {
			rep = CASREP_SYSERR;
		}
		if (e != mathlen) {
			rep = CASREP_MSGERR;
		}
	}
	
	if (rep == CASREP_OK) {
		char *bufp = buf;
		size_t len = mathlen;
		*tree = mkexprtree((const char **)&bufp, &len);
		if (!*tree || len != 0) {
			logmsg("Could not parse tree ");
			rep = CASREP_MSGERR;
		}
	}

	return rep;
}
