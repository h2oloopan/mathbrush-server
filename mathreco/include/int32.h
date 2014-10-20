#ifndef INT32_H_
#define INT32_H_

#ifdef WIN32
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int16 int16_t;
typedef unsigned __int16 uint16_t;
#else
#include <stdint.h>
#endif

#endif