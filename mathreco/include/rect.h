/*
 * rect.h
 *
 * This file defines a generic structure for manipuating rectangles as well
 * as several helper functions.
 */


#ifndef RECT_H_
#define RECT_H_


#include "dlldecl.h"
#include "verb.h"

#include <algorithm>
#include <cmath>
#include <ostream>


namespace scg
{


template <typename T>
struct DLLDECL Rect
{
    T left;
    T top;
    T right;
    T bottom;
    
    Rect() {}
    Rect(T l, T t, T r, T b) : left(l), top(t), right(r), bottom(b) {}

    template <typename U>
    inline bool operator==(const Rect<U> &rhs) const {
        return left==rhs.left && top == rhs.top && right == rhs.right && bottom == rhs.bottom;
    }

	template <typename U>
	inline bool operator!=(const Rect<U> &rhs) const {
		return !(*this == rhs);
	}
    
    template <typename U>
	inline bool operator<(const Rect<U> &rhs) const	{
		return (left < rhs.left) || ((left == rhs.left) && ((top < rhs.top || ((top == rhs.top) && ((right < rhs.right) || ((right == rhs.right) && (bottom < rhs.bottom)))))));
	}

	T width() const { return right - left; }
	T height() const { return bottom - top; }
	T x_center() const { return (left + right) / 2; }
	T y_center() const { return (top + bottom) / 2; }

	T area() const { return std::max(T(0), width()) * std::max(T(0), height()); }

	double height_to_width_ratio() const { return (double)height() / width(); }
	double width_to_height_ratio() const { return (double)width() / height(); }
};

template <typename T>
std::ostream &
operator<<(std::ostream &os, const Rect<T> &r) {
    return os << "(" << r.left << "," << r.top << ")-(" << r.right << "," << r.bottom << ")";
}


template <typename T>
Rect<T>
mkrect(T l, T t, T r, T b) {
    return Rect<T>(l, t, r, b);
}


template <typename T>
bool
rectangles_overlap(const Rect<T> &r1, const Rect<T> &r2) {
	if (r1.right < r2.left || r1.left > r2.right || r1.bottom < r2.top || r1.top > r2.bottom) {
		return false;
	}
	return true;  
}


template <typename T>
bool
rectangles_overlap_completely(const Rect<T> &r1, const Rect<T> &r2) {
	if (r1.left >= r2.left && r1.right <= r2.right
	 && r1.top >= r2.top && r1.bottom <= r2.bottom) {
		return true;
	}

	if (r2.left >= r1.left && r2.right <= r1.right
	 && r2.top >= r1.top && r2.bottom <= r1.bottom) {
		return true;
	}

	return false;
}


template <typename T>
Rect<T>
intersection(const Rect<T> &r1, const Rect<T> &r2) {
    T left = std::max(r1.left, r2.left);
	T right = std::max(left, std::min(r1.right, r2.right));
    T top = std::max(r1.top, r2.top);
	T bottom = std::max(top, std::min(r1.bottom, r2.bottom));
    
    return Rect<T>(left, top, right, bottom);
}


template <typename T>
double
overlap_proportion(const Rect<T> &r1, const Rect<T> &r2) {
    if (!rectangles_overlap(r1, r2)) {
        return 0.0;
    }
    
		Rect<T> cp1 = r1, cp2 = r2;
		if (cp1.left == cp1.right) cp1.right++;
		if (cp1.top == cp1.bottom) cp1.bottom++;
		if (cp2.left == cp2.right) cp2.right++;
		if (cp2.top == cp2.bottom) cp2.bottom++;

    double intersect_area = static_cast<double>(intersection(cp1, cp2).area());
    return intersect_area / std::min(cp1.area(), cp2.area());
}


template <typename T>
Rect<T>
merge(const Rect<T> &r1, const Rect<T> &r2) {
    /*if (r1.right < r1.left) {
        return merge(Rect<T>(r1.right, r1.top, r1.left, r1.bottom), r2);
    }
    if (r2.right < r2.left) {
        return merge(r1, Rect<T>(r2.right, r2.top, r2.left, r2.bottom));
    }
    if (r1.bottom < r1.top) {
        return merge(Rect<T>(r1.left, r1.bottom, r1.right, r1.top), r2);
    }
    if (r2.bottom < r2.top) {
        return merge(r1, Rect<T>(r2.left, r2.bottom, r2.right, r2.top));
    }*/
    
    return Rect<T>(std::min(r1.left, r2.left),
                   std::min(r1.top, r2.top),
                   std::max(r1.right, r2.right),
                   std::max(r1.bottom, r2.bottom));
}


}


#endif
