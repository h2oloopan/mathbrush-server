#ifndef WIN32_UTIL_H_
#define WIN32_UTIL_H_


#include <istream>
#include <string>

#include "group.h"


namespace scg
{


int GetProfilePath(std::string &path);
int GetProfilePathW(std::wstring &path);
int GetUserProfilePath(std::string &path);
int GetUserProfilePathW(std::wstring &path);

//void SetTrainingPath(const char *path);


int ReadMicrosoftStrokeGroupInk(std::istream &is, RawStrokeGroup &strokes, unsigned ink_sz);
int ReadMicrosoftStrokeGroupInk(std::istream &is, RawStrokeGroup &strokes);


}


#endif
