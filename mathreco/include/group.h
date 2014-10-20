#ifndef GROUP_H_
#define GROUP_H_

#include "stroke.h"
#include <numeric>
#include <cmath>
#include <cassert>

namespace scg {

template <typename T>
struct StrokeGroup {
	Stroke<T> *strokes;
	size_t nstrokes;
	
	StrokeGroup() : strokes(0), nstrokes(0) { }
	StrokeGroup(Stroke<T> *strokes_, size_t nstrokes_) : strokes(strokes_), nstrokes(nstrokes_) { }
	~StrokeGroup() { clear(); }
	
	bool empty() const { return nstrokes == 0; }

	Stroke<T> *release() {
		Stroke<T> *p = strokes;
		strokes = 0;
		nstrokes = 0;
		return p;
	}

	void clear()	{
		delete[] strokes;
		strokes = 0;
		nstrokes = 0;
	}

	void set_strokes(Stroke<T> *s, size_t n) {
		clear();
		strokes = s;
		nstrokes = n;
	}

	StrokeGroup &operator=(const StrokeGroup &g) {
		if (this != &g) {
			strokes = g.strokes;
			nstrokes = g.nstrokes;
			StrokeGroup &gg = (StrokeGroup &)g;
			gg.strokes = 0;
			gg.nstrokes = 0;
		}
		return *this;
	}

	size_t npoints() const {
		size_t np = 0;
		for (size_t i = 0; i < nstrokes; ++i) {
			np += strokes[i].npoints;
		}
		return np;
	}

	T x_center() const {
		T xc = 0;
		for (size_t i = 0; i < nstrokes; ++i) {
			xc += std::accumulate(strokes[i].x, strokes[i].x + strokes[i].npoints, T(0));
		}
		return xc / npoints();
	}
	T y_center() const {
		T yc = 0;
		for (size_t i = 0; i < nstrokes; ++i) {
			yc += std::accumulate(strokes[i].y, strokes[i].y + strokes[i].npoints, T(0));
		}
		return yc / npoints();
	}

	StrokeGroup<T> &translate(T dx, T dy) {
		for (size_t i = 0; i < nstrokes; ++i) {
			strokes[i].translate(dx, dy);
		}
		return *this;
	}

	Rect<T> bounds() const {
		if (nstrokes == 0) {
			return Rect<T>();
		}
		Rect<T> bbox = strokes[0].bounds();
		for (size_t i = 1; i < nstrokes; ++i) {
			bbox = merge(bbox, strokes[i].bounds());
		}
		return bbox;
	}

	StrokeGroup<T> copy() const {
		Stroke<T> *cps = new Stroke<T>[nstrokes];
		for (size_t i = 0; i < nstrokes; ++i) {
			Stroke<T> S = strokes[i].copy();
			cps[i] = S;
		}
		return StrokeGroup<T>(cps, nstrokes);
	}
};

typedef StrokeGroup<long> RawStrokeGroup;
typedef StrokeGroup<double> NormalizedStrokeGroup;

template <typename T>
double
lspdist(const StrokeGroup<T> &g1, const StrokeGroup<T> &g2, const std::vector<size_t> &inputorder) {
	assert(g1.nstrokes <= g2.nstrokes);
	double d = 0;
	for (size_t i = 0; i < g1.nstrokes; ++i) {
		d += lspdist(g1.strokes[i], g2.strokes[inputorder[i]]);
	}
	return d/g1.nstrokes;//std::pow(d, 1.0/g1.nstrokes);
}

template <typename T>
double
lspscore(const StrokeGroup<T> &g1, const StrokeGroup<T> &g2) {
	assert(g1.nstrokes <= g2.nstrokes);
	double d = 1;
	for (size_t i = 0; i < g1.nstrokes; ++i) {
		d *= lspscore(g1.strokes[i], g2.strokes[i]);
	}
	return std::pow(d, 1.0/g1.nstrokes);
}

}


#endif

