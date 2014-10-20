#include "cache.h"


namespace scg {


void
cachegrp::manage(cacheable &value) {
	values.push_back(&value);
}

void
cachegrp::invalidate() {
	for (std::list<cacheable *>::iterator i = values.begin(); i != values.end(); ++i) {
		(*i)->invalidate();
	}
}


}
