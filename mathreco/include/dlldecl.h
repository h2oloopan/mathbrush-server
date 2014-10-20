#ifndef DLLDECL_H_
#define DLLDECL_H_


#ifdef WIN32
# ifdef STATIC_LIBMATHRECO
#  define DLLDECL 
# else
#  ifdef BUILDING_LIBMATHRECO
#   define DLLDECL __declspec(dllexport)
#  else
#   define DLLDECL __declspec(dllimport)
#  endif
# endif
#else
# define DLLDECL 
#endif

#endif

