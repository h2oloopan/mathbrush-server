#include "utils.h"
#include <string>

namespace scg {

static std::string user_profile_path;

void
SetUserProfilePath(const char *path) {
	user_profile_path = path;
}

int
GetUserProfilePath(std::string &path) {
	if (!user_profile_path.empty()) {
		path = user_profile_path;
		return 0;
	}
	return GetProfilePath(path);
}

}
