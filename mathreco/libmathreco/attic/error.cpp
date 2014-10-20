#include "error.h"

#ifdef WIN32
#include <windows.h>
#endif

#include <iostream>


namespace scg
{


int errval;
std::ostream &errout = std::cerr;


std::string error_message(int code)
{
#ifdef WIN32
	if (IS_API_ERROR(code)) {
		DWORD api_error = GET_API_ERROR(code);
		
		LPTSTR buf;
		if (!FormatMessage(
			FORMAT_MESSAGE_ALLOCATE_BUFFER | 
			FORMAT_MESSAGE_FROM_SYSTEM | 
			FORMAT_MESSAGE_IGNORE_INSERTS,
			NULL,
			api_error,
			MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
			(LPTSTR)&buf,
			0,
			0))
		{
			return "API error";
		}
		
		char *cbuf;
#ifdef UNICODE
        size_t n = 1+wcstombs(0, buf, 0);
        cbuf = new char[n];
        wcstombs(cbuf, buf, n);
#else
        cbuf = buf;
#endif
		std::string ret(cbuf);
#ifdef UNICODE
        delete[] cbuf;
#endif
		LocalFree(buf);
		return ret;
	}
#endif
	
	switch (code) {
	case E_NONE:
		return "no error";
	case E_INVALID:
		return "invalid parameter";
	case E_IO:
		return "I/O error";
	case E_NOTFOUND:
	    return "the requested item was not found";
	case E_OUTOFMEM:
	    return "out of memory";
	case E_INTERNAL:
	    return "internal error";
	case E_NOTREADY:
	    return "object not ready";
	default:
		return "unknown error code";
	}
}


}

