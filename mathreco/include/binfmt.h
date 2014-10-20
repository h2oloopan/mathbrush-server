#ifndef BINFMT_H_
#define BINFMT_H_

#ifdef WIN32
#define NOMINMAX
//#include <winsock2.h>
typedef __int32 bin32_t;
/*extern "C" {
extern bin32_t ntohl(bin32_t);
extern bin32_t htonl(bin32_t);
}*/
#else
#include <arpa/inet.h>
#include <stdint.h>
typedef uint32_t bin32_t;
#endif

#endif