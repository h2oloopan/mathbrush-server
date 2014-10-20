#ifndef TIMING_H_
#define TIMING_H_


#ifdef WIN32
#	include <windows.h>
#else
#	include <sys/time.h>
#	include <sys/resource.h>
#endif

namespace scg {


#ifdef WIN32

#define TIME_OPERATION(code) \
	double ELAPSED; \
	do { \
		FILETIME dummy, start, end; \
		GetProcessTimes(GetCurrentProcess(), &dummy, &dummy, &dummy, &start); \
		code; \
		GetProcessTimes(GetCurrentProcess(), &dummy, &dummy, &dummy, &end); \
		ULARGE_INTEGER istart, iend; \
		istart.LowPart = start.dwLowDateTime; istart.HighPart = start.dwHighDateTime; \
		iend.LowPart = end.dwLowDateTime; iend.HighPart = end.dwHighDateTime; \
		istart.QuadPart = iend.QuadPart - istart.QuadPart; \
		ELAPSED = 1e-7 * istart.QuadPart; \
	} while(0)

#else

#define TIME_OPERATION(code) \
	double ELAPSED; \
	struct timeval ELAPSED_; \
	do { \
		struct rusage start, end; \
		getrusage(RUSAGE_SELF, &start); \
		code; \
		getrusage(RUSAGE_SELF, &end); \
		timersub(&end.ru_utime, &start.ru_utime, &ELAPSED_); \
		ELAPSED = ELAPSED_.tv_sec + 1e-6 * ELAPSED_.tv_usec; \
	} while(0)


//int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

#endif

}


#endif
