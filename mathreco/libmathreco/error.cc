#include "error.h"

#ifdef WIN32
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>


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
	case E_OUTOFMEM:
	    return "out of memory";
	case E_NOTFOUND:
	    return "the requested item was not found";
	case E_INVALID:
		return "invalid parameter";
	case E_IO:
		return "I/O error";
	case E_INTERRUPT:
		return "interrupted";
	case E_INTERNAL:
	    return "internal error";
	case E_NOTREADY:
	    return "not ready";
	case E_IGNORE:
		return "ignoring previous error";
	case E_EOF:
		return "end of file";
	case E_ALREADY:
		return "operation already complete";
	case E_SYSTEM:
		return "hidden system error";
	case E_UNSUPPORTED:
		return "unsupported operation";
	case E_OUTOFDATE:
		return "out of date";
	default:
		return "unknown error code";
	}
}


static error curr_error(E_INFO, E_NONE, std::string());
static std::ostream *error_stream = 0;
static int log_level = E_WARNING;


const error &get_error() { return curr_error; }

void
set_error(int level, int code, const std::string &msg)
{
	curr_error.level = level;
	curr_error.code = code;
	curr_error.msg = msg;

	if (error_stream && level >= log_level) {
		*error_stream << msg << std::endl;
		if (!*error_stream && level != E_FATAL) {
			set_error(E_ERROR, E_IO, "while writing to error stream");
			throw get_error().code;
		}
	}

	if (level == E_FATAL) {
		std::abort();
	}
}


void
set_log_level(int level)
{
	log_level = level;
}

void
set_error_log(const std::string &filename)
{
	std::ofstream *ofs = new std::ofstream(filename.c_str());
	if (!ofs->is_open() || !*ofs) {
		ERR(E_IO, "cannot open new error stream " << filename);
		delete ofs;
	}
	else {
		error_stream = ofs;
	}
}


void
set_error_log(std::ostream &os)
{
	error_stream = &os;
}

std::ostream &
operator<<(std::ostream &os, const error &e)
{
	return os << e.msg;
}


}

