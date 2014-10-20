#ifndef ERROR_H_
#define ERROR_H_


#include <string>
#include <ostream>
#include <sstream>


// NB: on Win32, bit 29 of all error codes is reserved for applications, so use it
// to flag an API-specific error.
#ifdef WIN32

#define MAKE_API_ERROR(code) ((code) | E_APIFAILURE)
#define IS_API_ERROR(code) ((code) & E_APIFAILURE)
#define GET_API_ERROR(code) ((code) & ~E_APIFAILURE)

#endif

#define FAILURE(code) ((code) != 0)
#define SUCCESS(code) !FAILURE(code)



#define E_NONE         0
#define E_OK           0
#define E_OUTOFMEM     1
#define E_NOTFOUND     2
#define E_INVALID      3
#define E_IO           4
#define E_INTERRUPT    5
#define E_INTERNAL     6
#define E_NOTREADY     7
#define E_IGNORE       8
#define E_EOF          9
#define E_ALREADY     10
#define E_SYSTEM      11
#define E_UNSUPPORTED 12
#define E_OUTOFDATE   13


#define E_APIFAILURE (1 << 29)



#define E_INFO    0
#define E_WARNING 1
#define E_ERROR   2
#define E_FATAL   3



#define GENERIC_ERROR(level, code, msgcons) do { \
	std::stringstream ss; \
	ss << ::scg::error_message(code) << " : " << msgcons; \
	::scg::set_error(level, code, ss.str()); \
	} while (0)


// Windows defines ERROR, so we can't use it; use ERR instead

#define INFO(code, msgcons) GENERIC_ERROR(E_INFO, code, msgcons)
#define WARNING(code, msgcons) GENERIC_ERROR(E_WARNING, code, msgcons)
#define ERR(code, msgcons) GENERIC_ERROR(E_ERROR, code, msgcons)
#define FATAL(code, msgcons) GENERIC_ERROR(E_FATAL, code, msgcons)

#define THROW_LAST() throw ::scg::get_error().code
#define ERROR_CODE() ::scg::get_error().code

#define THROW_INFO(code, msgcons) \
	GENERIC_ERROR(E_INFO, code, msgcons); \
	THROW_LAST()

#define THROW_WARNING(code, msgcons) \
	GENERIC_ERROR(E_WARNING, code, msgcons); \
	THROW_LAST()

#define THROW_ERROR(code, msgcons) \
	GENERIC_ERROR(E_ERROR, code, msgcons); \
	THROW_LAST()

#define THROW_FATAL(code, msgcons) \
	GENERIC_ERROR(E_FATAL, code, msgcons); \
	THROW_LAST()

#define ENSURE_ERROR(code) do {	if (ERROR_CODE() != e) ERR(e, ""); } while (0)


namespace scg
{


struct error {
	int level;
	int code;
	std::string msg;

	error(int level_, int code_, const std::string &msg_)
		: level(level_), code(code_), msg(msg_)
		{ }
};

std::ostream &operator<<(std::ostream &os, const error &e);


const error &get_error();
void set_error(int level, int code, const std::string &msg);

void set_log_level(int level);
void set_error_log(const std::string &filename);
void set_error_log(std::ostream &os);




extern int errval;

extern std::ostream &errout;

std::string error_message(int code);


}


#endif

