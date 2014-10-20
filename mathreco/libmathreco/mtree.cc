#include "mtree.h"
#include <cmath>

namespace scg {

double
splitweight(double x) {
	x -= 0.5;
	return 1 - 4*x*x;//2*std::abs(x);
	return 1;
	x -= 0.5;
	return std::exp(-x*x / 8.0);
}

}
