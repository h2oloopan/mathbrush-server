#ifndef DIST_H_
#define DIST_H_


#include <limits>
#include <utility>
#include <vector>

#include "group.h"
#include "rect.h"
#include "stroke.h"
#include "verb.h"

namespace scg
{


// distance squared between two points.
template <typename T>
T dist_sq(T x0, T y0, T x1, T y1) {
	T dx = x1 - x0;
	T dy = y1 - y0;
	return dx * dx + dy * dy;
}


// distance squared between two rectangles.
template <typename T>
T
dist_sq(const Rect<T> &r0, const Rect<T> &r1) {
	/*
		There are four cases:
			Case 1 - the rectangles overlap in y-space.  The distance is then the length of a horizontal line
			         between the rectangles or 0 if they also overlap in x-space.
			Case 2 - the rectangles overlap in x-space.  Essentially the same as case #1
			Case 3 - "r1" is completely below "r0" in y-space (WOLOG). The distance is then the straight-line
			         distance from a top corner of "r1" to a bottom corner of "r0" (choose whichever corners
			         give the shortest distance).
			Case 4 - "r1" is completely above "r0" in y-space (WOLOG). Swap "top" and "bottom" in case #3 and
			         then this is the same thing.
	*/
	
	// CASE #1
	if ((r1.top >= r0.top && r1.top <= r0.bottom) || (r1.bottom >= r0.top && r1.bottom <= r0.bottom)) {
		if (r1.right < r0.left) {
			T d = r0.left - r1.right;
			return d * d;
		}
		else if (r1.left > r0.right) {
			T d = r1.left - r0.right;
			return d * d;
		}
		else {
			return T(0);
		}
	}
	
	// CASE #2
	if ((r1.left >= r0.left && r1.left <= r0.right) || (r1.right >= r0.left && r1.right <= r0.right)) {
		if (r1.bottom < r0.top) {
			T d = r0.top - r1.bottom;
			return d * d;
		}
		else if (r1.top > r0.bottom) {
			T d = r1.top - r0.bottom;
			return d * d;
		}
		else {
			return T(0);
		}
	}

	// CASE #3
	if (r1.top > r0.bottom) {
		if (r1.right < r0.left) {
			return dist_sq(r1.right, r1.top, r0.left, r0.bottom);
		}
	
		if (r1.left > r0.right) {
			return dist_sq(r1.left, r1.top, r0.right, r0.bottom);
		}
	}

	// CASE #4
	if (r1.bottom < r0.top) {
		if (r1.right < r0.left) {
			return dist_sq(r1.right, r1.bottom, r0.left, r0.top);
		}

		if (r1.left > r0.right) {
			return dist_sq(r1.left, r1.bottom, r0.right, r0.top);
		}
	}

	// the above code only tests each case for one of the two possibilities, so check the others
	// if nothing got picked up
	return dist_sq(r1, r0);
}


template <typename T>
T
dist_sq(T x, T y, T u1, T v1, T u2, T v2) {
	double du = u2 - u1;
	double dv = v2 - v1;
	T dx = x - u1;
	T dy = y - v1;
	double t = (du*dx + dv*dy)/(du*du + dv*dv);
	double u, v;
	if (t <= 0) {
		u = u1;
		v = v1;
	}
	else if (t >= 1) {
		u = u2;
		v = v2;
	}
	else {
		u = u1 + t*du;
		v = v1 + t*dv;
	}
	du = u-x;
	dv = v-y;
	return T(du*du+dv*dv);
}

template <typename T>
T
dist_sq(T x1, T y1, T x2, T y2,
        T u1, T v1, T u2, T v2) {
	T denom = u1*y1+v1*x2+v2*x1+u2*y2-u1*y2-v1*x1-u2*y1-v2*x2;
	if (denom == 0) {
		T dx = x1 - u1;
		T dy = y1 - v1;
		return dx*dx + dy*dy;
	}
	double t1 = (double)(v2*x1+u2*v1+u1*y1-v1*x1-u1*v2-u2*y1)/denom;
	double t2 = (double)(u1*y1+x1*y2+x2*v1-u1*y2-x2*y1-v1*x1)/denom;
	t1 = std::min(1.0, std::max(0.0, t1));
	t2 = std::min(1.0, std::max(0.0, t2));
	double x = x1 + t1*(x2-x1);
	double y = y1 + t1*(y2-y1);
	double u = u1 + t2*(u2-u1);
	double v = v1 + t2*(v2-v1);
	double dx = x-u;
	double dy = y-v;
	return T(dx*dx+dy*dy);
}


template <typename T>
T
dist_sq(const Stroke<T> &s1, const Stroke<T> &s2) {
	static const T zero = T(0);

	T min_dist = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
    
	if (s1.npoints == 1) {
		if (s2.npoints == 1) {
			return dist_sq(s1.x[0], s1.y[0], s2.x[0], s2.y[0]);
		}
		else {
			const T *x2, *y2;
			const T *x2end = s2.x + s2.npoints - 1;
			for (x2 = s2.x, y2 = s2.y; x2 != x2end; x2++, y2++) {
				T dist = dist_sq(s1.x[0], s1.y[0], x2[0], y2[0], x2[1], y2[1]);
				if (dist < min_dist) {
					if (dist == zero) {
						return zero;
					}
					min_dist = dist;
				}
			}
		}
	}
	else if (s2.npoints == 1) {
		const T *x1, *y1;
		const T *x1end = s1.x + s1.npoints - 1;
		for (x1 = s1.x, y1 = s1.y; x1 != x1end; x1++, y1++) {
			T dist = dist_sq(s2.x[0], s2.y[0], x1[0], y1[0], x1[1], y1[1]);
			if (dist < min_dist) {
				if (dist == zero) {
					return zero;
				}
				min_dist = dist;
			}
		}
	}
	else {
		const T *x1, *y1;
		const T *x1end = s1.x + s1.npoints - 1;
		const T *x2, *y2;
		const T *x2end = s2.x + s2.npoints - 1;
		for (x1 = s1.x, y1 = s1.y; x1 != x1end; x1++, y1++) {
			for (x2 = s2.x, y2 = s2.y; x2 != x2end; x2++, y2++) {
				T dist = dist_sq(x1[0], y1[0], x1[1], y1[1], x2[0], y2[0], x2[1], y2[1]);
				if (dist < min_dist) {
					if (dist == zero) {
						return zero;
					}
					min_dist = dist;
				}
			}
		}
	}
	return min_dist;
}


template <typename T>
T
dist_sq(T x, T y, const Stroke<T> &s){
    static const T zero = T(0);
	T min_dist = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();    
    
    const T *x1, *y1;
    const T *x1end = s.x + s.npoints - 1;
    
    for (x1 = s.x, y1 = s.y; x1 != x1end; x1++, y1++) {
        T dist = dist_sq(x, y, *x1, *y1);
        if (dist < min_dist) {
			if (dist == zero) {
				return zero;
			}
            min_dist = dist;
        }
    }

    return min_dist;
}

/*
template <typename T>
unsigned
num_intersections(const Stroke<T> &s1, const Stroke<T> &s2)
{
    unsigned intersections = 0;
        
    const T zero = T(0);
    
    const T *x1, *y1;
    const T *x1end = s1.x + num_points(s1) - 1;
    const T *x2, *y2;
    const T *x2end = s2.x + num_points(s2) - 1;
    
    for (x1 = s1.x, y1 = s1.y; x1 != x1end; x1++, y1++) {
        for (x2 = s2.x, y2 = s2.y; x2 != x2end; x2++, y2++) {
            T dist = dist_sq(*x1, *y1, *(x1 + 1), *(y1 + 1),
                             *x2, *y2, *(x2 + 1), *(y2 + 1));
            if (dist == zero) {
                intersections++;
            }
        }
    }

    return intersections;
}


template <typename T, template <typename U> class G>
unsigned
num_intersections(const G<T> &g1, const G<T> &g2)
{
    unsigned intersections = 0;
    
    for (typename G<T>::const_iterator i = g1.begin(); i != g1.end(); ++i) {
        for (typename G<T>::const_iterator j = g2.begin(); j != g2.end(); ++j) {
            intersections += num_intersections(*i, *j);
        }
    }
    
    return intersections;
}


template <typename T>
bool
strokes_intersect(const Stroke<T> &s1, const Stroke<T> &s2)
{
    return dist_sq(s1, s2) == T(0);
}

// This next method returns the point index greater or equal to "from" at which
// the stroke self-intersects, or s1.size() if no such point exists.
template <typename T>
unsigned
stroke_self_intersection(const Stroke<T> &stroke, unsigned from, T dist_threshold = T(0))
{
    unsigned i;
    T *x = stroke.x + from + 2;
    T *y = stroke.y + from + 2;
    for (i = from + 2; i < num_points(stroke); i++) {
        if (dist_sq(stroke.x[from], stroke.y[from], stroke.x[from + 1], stroke.y[from + 1],
                    *x, *y, *(x + 1), *(y + 1)) <= dist_threshold) {
            break;
        }
    }
    
    return i;
}

template <typename T>
std::vector<std::pair<unsigned, unsigned> >
find_self_intersections(const Stroke<T> &stroke, T dist_threshold = T(0))
{
    std::vector<std::pair<unsigned, unsigned> > intersections;
    
    for (unsigned i = 0; i < stroke.size() - 1; i++) {
        unsigned j = stroke_self_intersection(stroke, i, dist_threshold);
        if (j < stroke.size()) {
            intersections.push_back(std::make_pair(i, j));
        }
    }
    
    return intersections;
}
*/

}


#endif

