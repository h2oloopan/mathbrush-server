#include "norm.h"

double
NormalScale()
{
	return 1.0;
}

int
NormToScreen(double x)
{
	return (int)(x * 3.5);
}

double
NormalizeToScale(double x, double scale)
{
	return (x / scale) * 1.0;
}

double
NormalClamp(double v)
{
	if (v < 0.0) {
		return 0.0;
	}
	else if (v > 1.0) {
		return 1.0;
	}
	else {
		return v;
	}
}

