#ifndef NORMALIZE_H_
#define NORMALIZE_H_

#include "stroke.h"
#include "group.h"
#include "rect.h"
#include <algorithm>

namespace scg {

template <typename T>
NormalizedStroke
normalize(const Stroke<T> &s, const Rect<T> &bounds) {
	double *x = new double[s.npoints];
	double *y = new double[s.npoints];
	unsigned long *t = (s.time) ? new unsigned long[s.npoints] : 0;

	if (s.time) {
		std::copy(s.time, s.time + s.npoints, t);
	}
	
	double sc = std::max(bounds.width(), bounds.height());
	//double h = bounds.height();
	for (size_t i = 0; i < s.npoints; ++i) {
		x[i] = (s.x[i] - bounds.left) / sc;
		y[i] = (s.y[i] - bounds.top) / sc;
	}

	return NormalizedStroke(x, y, t, s.npoints);
}

template <typename T>
NormalizedStrokeGroup
normalize(const StrokeGroup<T> &g) {
	Rect<T> bounds = g.bounds();
	NormalizedStroke *s = new NormalizedStroke[g.nstrokes];
	for (size_t i = 0; i < g.nstrokes; ++i) {
		NormalizedStroke S = normalize(g.strokes[i], bounds);
		s[i] = S;
	}
	return NormalizedStrokeGroup(s, g.nstrokes);
}

}

#endif
