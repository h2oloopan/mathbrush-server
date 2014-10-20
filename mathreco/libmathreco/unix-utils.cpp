#include "utils.h"
#include "error.h"
#include "memory.h"

#include <cstdlib>


namespace scg
{


static std::string profile_path;

int
GetProfilePath(std::string &path) {
	if (!profile_path.empty()) {
		path = profile_path;
	}
	else {
	 path = std::getenv("LIBMATHRECO_ROOT");
	 path.append("/data-files/");
	}
	return 0;
}

void
SetProfilePath(const char *path) {
	profile_path = path;
}

int
ReadMicrosoftStrokeGroupInk(std::istream &is, RawStrokeGroup &strokes, unsigned ink_sz)
{
	return E_INVALID;
}

int
ReadMicrosoftStrokeGroupInk(std::istream &is, RawStrokeGroup &strokes)
{
	return E_INVALID;
}

}

