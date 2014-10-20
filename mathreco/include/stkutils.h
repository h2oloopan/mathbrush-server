#ifndef STKUTILS_H_
#define STKUTILS_H_

#include "rect.h"
#include "stroke.h"
#include "group.h"

namespace scg {

RawStroke triple(const RawStroke &stk);
RawStrokeGroup triple(const RawStrokeGroup &stks);

bool bbox_is_dot(const Rect<long> &bbox);
bool stroke_is_dot(const RawStroke &stroke);

template <typename T>
bool
stroke_is_dot(const Stroke<T> &stroke, const Rect<T> &group_bounds) {
	if (stroke.npoints == 0) {
		return false;
	}
	const Rect<T> &bounds = stroke.bounds();
	return (stroke.npoints <= 8 && ((double)bounds.width()/group_bounds.width() <= 0.12 || (double)bounds.height()/group_bounds.height() <= 0.12));
}

template <typename T>
bool
stroke_is_closed(const Stroke<T> &stroke, const Rect<T> &bounds) {
	if (stroke.npoints == 0 || stroke_is_dot(stroke, bounds)) {
		return false;
	}
	double dx = std::abs(stroke.x[stroke.npoints - 1] - stroke.x[0]);
	double dy = std::abs(stroke.y[stroke.npoints - 1] - stroke.y[0]);
	return dx <= 0.2 * bounds.width() && dy <= 0.2 * bounds.height();
}

template <typename T>
double *
tangents_along_stroke(const Stroke<T> &s, double *tangents) {
    if (s.npoints == 1) {
        tangents[0] = 0.0;
    }
    else {        
        T *x = s.x, *y = s.y;
        tangents[0] = std::atan2(*(y + 1) - *y, *(x + 1) - *x);
        x++, y++;
        size_t i;
        for (i = 1; i != s.npoints - 1; x++, y++, i++) {
            tangents[i] = std::atan2(*(y + 1) - *(y - 1), *(x + 1) - *(x - 1));
        }
        tangents[i] = std::atan2(*y - *(y - 1), *x - *(x - 1));
    }
    return tangents;
}

}

#endif
