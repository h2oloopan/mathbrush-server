#include "strutils.h"
#include <cstdlib>
#include <iostream>

namespace scg {

std::wstring
str2wstr(const std::string &s) {
	size_t mbsz = mbstowcs(0, s.c_str(), 0);
	std::wstring ws;
	//ws.resize(mbsz);
	wchar_t *wcs = new wchar_t[mbsz+1];
	if (mbstowcs(wcs, s.c_str(), mbsz+1) != mbsz || wcs[mbsz] != 0) {
		abort();
	}
	ws = wcs;
	delete[] wcs;
	return ws;
}

std::string
wstr2str(const std::wstring &ws) {
	size_t sz = wcstombs(0, ws.c_str(), 0);
	std::string s;
	char *cs = new char[sz+1];
	if (wcstombs(cs, ws.c_str(), sz+1) != sz || cs[sz] != 0) {
		abort();
	}
	s = cs;
	delete[] cs;
	return s;
}

}
