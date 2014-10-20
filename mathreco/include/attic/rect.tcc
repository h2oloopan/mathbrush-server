#include <algorithm>
#include <cmath>


namespace scg {


template <typename T>
rect<T>::rect()
	: left(0), top(0), right(0), bottom(0)
	{ }

template <typename T>
rect<T>::rect(T left_, T top_, T right_, T bottom_)
	: left(left_), top(top_), right(right_), bottom(bottom_)
	{ }

template <typename T>
T
rect<T>::width() const {
	return std::abs(right - left);
}

template <typename T>
T
rect<T>::height() const {
	return std::abs(bottom - top);
}

template <typename T>
T
rect<T>::area() const {
	return width() * height();
}


template <typename T>
bool
rect<T>::operator==(const rect<T> &rhs) const {
	return left == rhs.left && top == rhs.top &&
	       right == rhs.right && bottom == rhs.bottom;
}

template <typename T>
bool
rect<T>::operator!=(const rect<T> &rhs) const {
	return !(*this == rhs);
}

template <typename T>
bool
rect<T>::lexicographical_cmp::operator()(const rect<T> &lhs,
                                         const rect<T> &rhs) const {
	return (lhs.left < rhs.left) || (lhs.left == rhs.left &&
	        ((lhs.top < rhs.top) || (lhs.top == rhs.top &&
	         ((lhs.right < rhs.right) || (lhs.right == rhs.right &&
					  (lhs.bottom < rhs.bottom))))));
}


template <typename T>
rect<T> &
rect<T>::merge(const rect<T> &with) {
	left = std::min(left, with.left);
	right = std::max(right, with.right);
	top = std::min(top, with.top);
	bottom = std::max(bottom, with.bottom);
	return *this;
}


template <typename T>
rect<T>
mkrect(T left, T top, T right, T bottom) {
	return rect<T>(left, top, right, bottom);
}

template <typename From, typename To>
rect<To>
convert_rect(const rect<From> &rc) {
	return mkrect(To(rc.left), To(rc.top), To(rc.right), To(rc.bottom));
}

template <typename T>
rect<T>
merge(const rect<T> &first, const rect<T> &second) {
	rect<T> rect = first;
	return first.merge(second);
}

template <typename T>
rect<T>
intersection(const rect<T> &first, const rect<T> &second) {
	rect<T> n;
	n.left = std::max(first.left, second.left);
	n.right = std::max(n.left, std::min(first.right, second.right));
	n.top = std::max(first.top, second.top);
	n.bottom = std::max(n.top, std::min(first.bottom, second.bottom));
	return n;
}

template <typename T>
static T
distsq(T x1, T y1, T x2, T y2) {
	T dx = x2 - x1;
	T dy = y2 - y1;
	return dx * dx + dy * dy;
}

template <typename T>
T
distsq(const rect<T> &r0, const rect<T> &r1) {
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
	if ((r1.top >= r0.top && r1.top <= r0.bottom)
	 || (r1.bottom >= r0.top && r1.bottom <= r0.bottom)) {
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
			return distsq(r1.right, r1.top, r0.left, r0.bottom);
		}
	
		if (r1.left > r0.right) {
			return distsq(r1.left, r1.top, r0.right, r0.bottom);
		}
	}

	// CASE #4
	if (r1.bottom < r0.top) {
		if (r1.right < r0.left) {
			return distsq(r1.right, r1.bottom, r0.left, r0.top);
		}

		if (r1.left > r0.right) {
			return distsq(r1.left, r1.bottom, r0.right, r0.top);
		}
	}

	// the above code only tests each case for one of the two possibilities, so check the others
	// if nothing got picked up
	return distsq(r1, r0);
}

template <typename T>
double
overlap(const rect<T> &box1, const rect<T> &box2) {
	rect<T> n = intersection(box1, box2);
	return double(n.area()) / std::min(box1.area(), box2.area());
}

template <typename T>
std::ostream &
operator<<(std::ostream &os, const rect<T> &rc) {
	return os << rc.left << ' ' << rc.top << ' '
	          << rc.right << ' ' << rc.bottom;
}
template <typename T>
std::istream &
operator>>(std::istream &is, rect<T> &rc) {
	return is >> rc.left >> rc.top >> rc.right >> rc.bottom;
}


}
