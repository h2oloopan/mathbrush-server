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
    /*Rect &operator=(const Rect &r)
    {
        if (this != &r) {
            left = r.left;
            top = r.top;
            right = r.right;
            bottom = r.bottom;
        }
        return *this;
    }*/

    template <typename U>
    inline bool operator==(const Rect<U> &rhs) const
    {
        return left==rhs.left && top == rhs.top && right == rhs.right && bottom == rhs.bottom;
    }

		template <typename U>
		inline bool operator!=(const Rect<U> &rhs) const
		{
			return !(*this == rhs);
		}
    
    template <typename U>
	inline bool operator<(const Rect<U> &rhs) const
	{
		return (left < rhs.left) || ((left == rhs.left) && ((top < rhs.top || ((top == rhs.top) && ((right < rhs.right) || ((right == rhs.right) && (bottom < rhs.bottom)))))));
	}
};

template <typename T>
std::ostream &
operator<<(std::ostream &os, const Rect<T> &r)
{
    return os << "(" << r.left << "," << r.top << ")-(" << r.right << "," << r.bottom << ")";
}


template <typename T>
Rect<T>
make_rect(T l, T t, T r, T b)
{
    return Rect<T>(l, t, r, b);
}


template <typename T>
inline T
width(const Rect<T> &r)
{
    return std::abs(r.right - r.left);
}

template <typename T>
inline T
height(const Rect<T> &r)
{
    return std::abs(r.bottom - r.top);
}


template <typename T>
inline T
x_center(const Rect<T> &r)
{
    return (r.left + r.right) / T(2);
}

template <typename T>
inline T
y_center(const Rect<T> &r)
{
    return (r.top + r.bottom) / T(2);
}


template <typename T>
inline double
height_width_ratio(const Rect<T> &r)
{
    return static_cast<double>(height(r)) / width(r);
}

template <typename T>
inline double
width_height_ratio(const Rect<T> &r)
{
    return static_cast<double>(width(r)) / height(r);
}


template <typename T>
bool
rectangles_overlap(const Rect<T> &r1, const Rect<T> &r2)
{
	if (r1.right < r2.left || r1.left > r2.right || r1.bottom < r2.top || r1.top > r2.bottom) {
		return false;
	}
	return true;  
}


template <typename T>
bool
rectangles_overlap_completely(const Rect<T> &r1, const Rect<T> &r2)
{
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
T
area(const Rect<T> &rc)
{
    return (rc.right - rc.left) * (rc.bottom - rc.top);
}


template <typename T>
Rect<T>
intersection(const Rect<T> &r1, const Rect<T> &r2)
{
    T left = std::max(r1.left, r2.left);
    T right = std::min(r1.right, r2.right);
    T top = std::max(r1.top, r2.top);
    T bottom = std::min(r1.bottom, r2.bottom);
    
    return Rect<T>(left, top, right, bottom);
}


template <typename T>
double
overlap_proportion(const Rect<T> &r1, const Rect<T> &r2)
{
    if (!rectangles_overlap(r1, r2)) {
        return 0.0;
    }
    
    double intersect_area = static_cast<double>(area(intersection(r1, r2)));
    return intersect_area / std::min(area(r1), area(r2));
}


template <typename T>
Rect<T>
merge(const Rect<T> &r1, const Rect<T> &r2)
{
    if (r1.right < r1.left) {
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
    }
    
    return Rect<T>(std::min(r1.left, r2.left),
                   std::min(r1.top, r2.top),
                   std::max(r1.right, r2.right),
                   std::max(r1.bottom, r2.bottom));
}


}


#endif
