#ifndef VERB_H_
#define VERB_H_


#include <ostream>

namespace scg
{
extern std::ostream *verb_out;
extern int verb_level;
}

#ifndef STRIP_VERBOSE
# define VERBOSE_AT(level, expr) do { if (::scg::verb_out && level <= ::scg::verb_level) { expr; } } while (false)
#else
# define VERBOSE_AT(level, expr) 
#endif
//# define VERBOSE_AT(level, expr) 

#define VERBOSE0(expr) VERBOSE_AT(0, expr)
#define VERBOSE(expr) VERBOSE_AT(1, expr)
#define VERBOSE1(expr) VERBOSE_AT(1, expr)
#define VERBOSE2(expr) VERBOSE_AT(2, expr)
#define VERBOSE3(expr) VERBOSE_AT(3, expr)
#define VERBOSE4(expr) VERBOSE_AT(4, expr)
#define VERBOSE5(expr) VERBOSE_AT(5, expr)
#define VERBOSE6(expr) VERBOSE_AT(6, expr)
#define VERBOSE7(expr) VERBOSE_AT(7, expr)
#define VERBOSE8(expr) VERBOSE_AT(8, expr)
#define VERBOSE9(expr) VERBOSE_AT(9, expr)

#endif
