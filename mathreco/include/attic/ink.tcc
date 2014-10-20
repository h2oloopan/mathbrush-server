#include "ink.h"
#include "func.h"
#include "lsp.h"
#include "error.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>

namespace scg {

extern size_t NEXTSTK;

template <typename T>
point<T>::point() {
}
template <typename T>
point<T>::point(T x_, T y_)
	: x(x_), y(y_) {
}

template <typename T>
template <typename U>
point<T>::point(const point<U> &pt)
	: x(pt.x), y(pt.y) {
}
	
template <typename T>
bool
point<T>::operator==(const point<T> &rhs) const {
	return x == rhs.x && y == rhs.y; 
}

template <typename T>
bool
point<T>::operator!=(const point<T> &rhs) const {
	return !(*this == rhs);
}

template <typename T>
bool
point<T>::lexicographic_cmp::operator()(const point<T> &lhs, const point<T> &rhs) const {
	return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
}

template <typename T>
point<T> &
point<T>::operator+=(const point<T> &rhs) {
	x += rhs.x;
	y += rhs.y;
	return *this;
}
template <typename T>
point<T> &
point<T>::operator-=(const point<T> &rhs) {
	x -= rhs.x;
	y -= rhs.y;
	return *this;
}
template <typename T>
template <typename U>
point<T> &
point<T>::operator*=(const U a) {
	x = T(x * a);
	y = T(y * a);
	return *this;
}
template <typename T>
template <typename U>
point<T> &
point<T>::operator/=(const U a) {
	x = T(x / a);
	y = T(y / a);
	return *this;
}
template <typename T>
point<T>
point<T>::operator+(const point<T> &rhs) const {
	return point<T>(x + rhs.x, y + rhs.y);
}
template <typename T>
point<T>
point<T>::operator-(const point<T> &rhs) const {
	return point<T>(x - rhs.x, y - rhs.y);
}
template <typename T>
template <typename U>
point<T>
point<T>::operator/(const U a) const {
	return point<T>(x / a, y / a);
}

template <typename T, typename U>
point<U>
operator*(const T lhs, const point<U> &rhs) {
	return point<T>(lhs * rhs.x, lhs * rhs.y);
}

template <typename T>
point<T> &
point<T>::apply(const affine_transform<T> &xf) {
	const T *m = xf.matrix;
	x = x * m[0] + y * m[1] + m[2];
	y = y * m[3] + y * m[4] + m[5];
	return *this;
}

template <typename T>
point<T>
mkpoint(T x, T y) {
	return point<T>(x, y);
}

template <typename From, typename To>
int
convert_point(point<To> &dst, const point<From> &pt) {
	dst.x = To(pt.x);
	dst.y = To(pt.y);
	return 0;
}

template <typename T>
std::ostream &
operator<<(std::ostream &os, const point<T> &pt) {
	return os << pt.x << ' ' << pt.y;
}
template <typename T>
std::istream &
operator>>(std::istream &is, point<T> &pt) {
	return is >> pt.x >> pt.y;
}


template <typename T>
affine_transform<T>::affine_transform() {
	matrix[0] = matrix[4] = T(1);
	matrix[1] = matrix[2] = matrix[3] = matrix[5] = T(0);
}

template <typename T>
affine_transform<T> &
affine_transform<T>::chain(const affine_transform<T> &next) {
	return *this = compose(*this, next);
}

template <typename T>
affine_transform<T> &
affine_transform<T>::scale(T x, T y) {
	matrix[0] *= x;
	matrix[1] *= x;
	matrix[2] *= x;
	matrix[3] *= y;
	matrix[4] *= y;
	matrix[5] *= y;
	return *this;
}
template <typename T>
affine_transform<T> &
affine_transform<T>::translate(T x, T y) {
	matrix[2] += x;
	matrix[5] += y;
	return *this;
}

template <typename T>
affine_transform<T>
compose(const affine_transform<T> &first,
        const affine_transform<T> &second) {
	affine_transform<T> ret();
	T *m1 = first.matrix;
	T *m2 = second.matrix;
	T *r = ret.matrix;
	r[0] = m2[0] * m1[0] + m2[1] * m1[3];
	r[1] = m2[0] * m1[1] + m2[1] * m1[4];
	r[2] = m2[0] * m1[2] + m2[1] * m1[5] + m2[2];
	r[3] = m2[3] * m1[0] + m2[4] * m1[3];
	r[4] = m2[3] * m1[1] + m2[4] * m1[4];
	r[5] = m2[3] * m1[2] + m2[4] * m1[5] + m2[5];
	return ret;
}


template <typename T>
void
stroke<T>::init_cache() {
	cache.manage(bbox);
}

template <typename T>
stroke<T>::stroke()
	: bbox(const_bound_mem_fun(&stroke<T>::compute_bbox, this)),
	  grp(0), constgrp(0), id_(NEXTSTK++) {
	init_cache();
}

template <typename T>
template <typename InIt>
stroke<T>::stroke(InIt first, InIt last)
	: points(first, last),
	  bbox(const_bound_mem_fun(&stroke<T>::compute_bbox, this)),
		grp(0), constgrp(0), id_(NEXTSTK++) {
	init_cache();
	mkcoeffs();
}

template <typename T>
stroke<T>::stroke(const stroke<T> &rhs)
	: bbox(const_bound_mem_fun(&stroke<T>::compute_bbox, this)),
	  points(rhs.points),
		lsp(rhs.lsp),
		grp(0), constgrp(0), id_(NEXTSTK++) {
	init_cache();
	bbox = rhs.bbox;
}

template <typename T>
stroke<T>::~stroke() {
	delete grp;
	delete constgrp;
}

template <typename T>
size_t
stroke<T>::npoints() const {
	return points.size();
}

template <typename T>
void
stroke<T>::clear() {
	cache.invalidate();
	points.clear();
}

template <typename T>
stroke<T> &
stroke<T>::apply(const affine_transform<T> &xf) {
	for (typename point_collection::iterator i = points.begin(); i != points.end(); ++i) {
		(*i).apply(xf);
	}
	for (std::vector<point<double> >::iterator i = lsp.begin(); i != lsp.end(); ++i) {
		(*i).apply(xf);
	}
	return *this;
}

template <typename T>
size_t
stroke<T>::id() const {
	return id_;
}

template <typename T>
stkgroup<stroke<T> > *
stroke<T>::selfgrp() {
	if (!grp) {
		grp = new stkgroup<stroke<T> >;
		grp->addstk(this);
	}
	return grp;
}
template <typename T>
stkgroup<const stroke<T> > *
stroke<T>::const_selfgrp() const {
	if (!constgrp) {
		constgrp = new stkgroup<const stroke<T> >;
		constgrp->addstk(this);
	}
	return constgrp;
}

template <typename T>
void
stroke<T>::update() {
	cache.invalidate();
	lsp.clear();
	mkcoeffs();
}

template <typename T>
rect<T>
stroke<T>::compute_bbox() const {
	rect<T> bbox;
	if (!points.empty()) {
		typename point_collection::const_iterator i = points.begin();
		bbox.left = bbox.right = i->x;
		bbox.top = bbox.bottom = i->y;
		for (++i; i != points.end(); ++i) {
			if (i->x < bbox.left) {
				bbox.left = i->x;
			}
			else if (i->x > bbox.right) {
				bbox.right = i->x;
			}
			if (i->y < bbox.top) {
				bbox.top = i->y;
			}
			else if (i->y > bbox.bottom) {
				bbox.bottom = i->y;
			}
		}
	}
	return bbox;
}

template <typename T>
void
stroke<T>::mkcoeffs() {
	rect<T> bounds = bbox();
	point<T> tr(bounds.left, bounds.top);
	double sc = 1.0 / std::max(bounds.right - bounds.left,
	                           bounds.bottom - bounds.top);
	lsp.insert(lsp.end(), LSPDEG + 1, point<double>(0.0, 0.0));
	std::vector<point<double> > acc;
	acc.insert(acc.end(), LSPDEG + 1, point<double>(0.0, 0.0));
	typename point_collection::const_iterator i;
	double L = 0.0;
	for (i = points.begin(); i != points.end() - 1; ++i) {
		const point<double> d = sc * point<double>(*(i+1) - *i);
		const point<double> s = sc * 0.5 * point<double>(*(i+1) + *i - 2*tr);
		double L0 = std::sqrt(d.x*d.x + d.y*d.y);
		double L1 = L;
		double L1k = L;
		L += L0;
		double L2k = L;
		for (size_t k = 0; k <= LSPDEG; ++k) {
			double al = (L2k - L1k) / (k+1);
			acc[k] += al * s;
			L1k *= L1;
			L2k *= L;
		}
	}

	double Lk = L;
	std::vector<point<double> >::iterator j;
	for (j = acc.begin(); j != acc.end(); ++j) {
		*j /= Lk;
		Lk *= L;
	}

	double n = 0.0;
	static double u = 0.125;
	for (size_t k = 1; k <= LSPDEG; ++k) {
		point<double> &coeff = lsp[k];
		for (size_t m = 0; m <= k; ++m) {
			coeff += B[k][m] * acc[m];
			if (m > 0) {
				double c = u * B[k][m] * m;
				coeff += c * sc * point<double>(*i - tr);
				if (m > 1) {
					c *= (m-1);
					coeff -= c * acc[m-2];
				}
			}
		}
		coeff -= u * B[k][1] * sc * point<double>(*points.begin() - tr);
		n += coeff.x*coeff.x + coeff.y*coeff.y;
	}

	n = std::sqrt(n);
	for (size_t k = 1; k <= LSPDEG; ++k) {
		lsp[k] /= n;
	}
}

template <typename T>
std::ostream &
operator<<(std::ostream &os, const stroke<T> &stk) {
	os << stk.npoints() << ' ';
	const typename stroke<T>::point_collection &pts = stk.points;
	for (typename stroke<T>::point_collection::const_iterator i = pts.begin(); i != pts.end(); ++i) {
		os << *i << ' ';
	}
	return os;
}

template <typename T>
std::istream &
operator>>(std::istream &is, stroke<T> &stk) {
	stk.clear();
	size_t n;
	is >> n;
	typename stroke<T>::point_collection &pts = stk.points;
	while (n--) {
		point<T> p;
		is >> p;
		// test for duplicate points in input stream
		if (pts.empty() || p != pts.back()) {
			pts.push_back(p);
		}
	}
	stk.update();
	return is;
}


template <typename From, typename To>
int
convert_stroke(stroke<To> &dst, const stroke<From> &stk) {
	dst.clear();
	for (typename stroke<From>::point_collection::const_iterator i = stk.points.begin(); i != stk.points.end(); ++i) {
		point<To> p;
		convert_point(p, *i);
		dst.points.push_back(p);
	}
	return 0;
}

template <typename T>
int
subdivide(stroke<T> &dst, const stroke<T> &stk, size_t npoints) {
	dst.clear();

	typename stroke<T>::point_collection &retpts = dst.points;

	double len = 0.0;
	typename stroke<T>::point_collection::const_iterator p;
	for (p = stk.points.begin(); p != stk.points.end() - 1; ++p) {
		T dx = (p+1)->x - p->x;
		T dy = (p+1)->y - p->y;
		len += std::sqrt(dx * dx + dy * dy);
	}

	double seglen = len / (npoints - 1);

	p = stk.points.begin();
	retpts.push_back(*p);

	double nextlen = seglen;
	double plen;
	len = 0.0;
	for (size_t i = 1; retpts.size() < npoints - 1; ) {
		while (len < nextlen) {
			plen = len;
			T dx = (p+1)->x - p->x;
			T dy = (p+1)->y - p->x;
			len += std::sqrt(dx * dx + dy * dy);
			++p;
		}
		
		while (nextlen <= len && retpts.size() < npoints - 1) {
			double t = (nextlen - plen) / (len - plen);
			retpts.push_back(mkpoint(T(t * p->x + (1.0 - t) * (p-1)->x),
			                         T(t * p->y + (1.0 - t) * (p-1)->y)));
			++i;
			nextlen += seglen;
		}
	}

	retpts.push_back(stk.points.back());
	return 0;
}


template <typename T>
affine_transform<float_t>
normalization_xform(const rect<T> &frombox, const float_rect &tobox) {
	float_t fromx = 0.5 * (frombox.left + frombox.right);
	float_t fromy = 0.5 * (frombox.top + frombox.bottom);
	float_t tox = 0.5 * (tobox.left + tobox.right);
	float_t toy = 0.5 * (tobox.top + tobox.bottom);

	float_t xscale = tobox.width() / frombox.width();
	float_t yscale = tobox.height() / frombox.height();

	affine_transform<float_t> xf;
	xf.translate(-fromx, -fromy).scale(xscale, yscale).translate(tox, toy);
	return xf;
}


template <typename T>
int
normalize(float_stroke &dst, const stroke<T> &stk, const float_rect &bounds) {
	rect<T> bbox = stk.bbox();
	affine_transform<float_t> xf = normalization_xform(bbox, bounds);
	int e;
	if (FAILURE(e = convert_stroke(dst, stk))) {
		ERR(geterr(), "Error normalizing ink stroke");
		return e;
	}
	dst.apply(xf);
	return 0;
}

template <typename T, typename U>
double
lspdist(const stroke<T> &lhs, const stroke<U> &rhs) {
	std::vector<point<double> >::const_iterator il = lhs.lsp.begin();
	std::vector<point<double> >::const_iterator ir = rhs.lsp.begin();
	double diff = 0.0;
	for ( ; il != lhs.lsp.end(); ++il, ++ir) {
		point<double> d = *il - *ir;
		diff += d.x*d.x + d.y*d.y;
	}
	return std::sqrt(diff);
}

template <typename T, typename U>
double
lspscore(const stroke<T> &lhs, const stroke<U> &rhs) {
	static const double norm = std::sqrt(2.0);
	return 1.0 - lspdist(lhs, rhs)/norm;
}

template <typename T>
void
stkgroup<T>::init_cache() {
	cache.manage(bbox);
	cache.manage(npoints);
}

template <typename T>
stkgroup<T>::stkgroup()
	: ownshp(own_weak),
	  bbox(const_bound_mem_fun(&stkgroup<T>::compute_bbox, this)),
	  npoints(const_bound_mem_fun(&stkgroup<T>::compute_npoints, this)) {
	init_cache();
}

template <typename T>
stkgroup<T>::stkgroup(const stkgroup<T> &rhs)
	: ownshp(own_weak),
	  strokes(rhs.strokes),
	  bbox(const_bound_mem_fun(&stkgroup<T>::compute_bbox, this)),
	  npoints(const_bound_mem_fun(&stkgroup<T>::compute_npoints, this)) {
	init_cache();
	bbox = rhs.bbox;
	npoints = rhs.npoints;
}

template <typename T>
stkgroup<T>::~stkgroup() {
	clear();
	for (typename std::map<const stkgroup<T> *, stkgroup<T> *>::iterator i = mergecache.begin(); i != mergecache.end(); ++i) {
		delete i->second;
	}
}


template <typename T>
typename stkgroup<T>::const_iterator
stkgroup<T>::begin() const {
	return strokes.begin();
}

template <typename T>
typename stkgroup<T>::const_iterator
stkgroup<T>::end() const {
	return strokes.end();
}

template <typename T>
T *
stkgroup<T>::stk(size_t i) {
	if (i >= strokes.size()) {
		THROW_ERROR(E_INVALID, "Stroke index " << i << " is too large for this stkgroup");
	}
	return strokes[i];
}

template <typename T>
const T *
stkgroup<T>::stk(size_t i) const {
	if (i >= strokes.size()) {
		THROW_ERROR(E_INVALID, "Stroke index " << i << " is too large for this stkgroup");
	}
	return strokes[i];
}

template <typename T>
bool
stkgroup<T>::hasstk(T *stk) const {
	return std::binary_search(strokes.begin(), strokes.end(), stk);
}

template <typename T>
int
stkgroup<T>::addstk(T *stk) {
	iterator i;
	i = std::upper_bound(strokes.begin(), strokes.end(), stk);
	strokes.insert(i, stk);
	cache.invalidate();
	return 0;
}

template <typename T>
int
stkgroup<T>::rmstk(T *stk) {
	iterator i;
	i = std::upper_bound(strokes.begin(), strokes.end(), stk);
	if (i != strokes.end() && *i == stk) {
		strokes.erase(i);
		cache.invalidate();
		return 0;
	}
	ERR(E_NOTFOUND, "The stroke " << stk << " was not found for removal in this stkgroup");
	return E_NOTFOUND;
}

template <typename T>
rect<typename stkgroup<T>::pt_t>
stkgroup<T>::compute_bbox() const {
	rect<pt_t> bbox;
	if (!strokes.empty()) {
		const_iterator i = strokes.begin();
		const stk_t *s = *i;
		bbox = s->bbox();
		for (++i; i != strokes.end(); ++i) {
			s = *i;
			bbox.merge(s->bbox());
		}
	}
	return bbox;
}

template <typename T>
size_t
stkgroup<T>::compute_npoints() const {
	size_t n = 0;
	for (const_iterator i = strokes.begin(); i != strokes.end(); ++i) {
		const stk_t *s = *i;
		n += s->npoints();
	}
	return n;
}


template <typename T>
size_t
stkgroup<T>::nstrokes() const {
	return strokes.size();
}

template <typename T>
void
stkgroup<T>::clear() {
	cache.invalidate();
	if (ownshp == own_ptrs) {
		for (iterator i = strokes.begin(); i != strokes.end(); ++i) {
			delete *i;
		}
	}
	else if (ownshp == own_array && !strokes.empty()) {
		delete[] strokes.front();
	}
	strokes.clear();
}

template <typename T>
stkgroup<T> &
stkgroup<T>::apply(const affine_transform<typename stkgroup<T>::pt_t> &xf) {
	for (iterator i = strokes.begin(); i != strokes.end(); ++i) {
		stk_t *s = *i;
		s->apply(xf);
	}
	return *this;
}

template <typename T>
void
stkgroup<T>::chown(typename stkgroup<T>::ownstat st) {
	ownshp = st;
}


extern const std::string INK_MAGIC_STR;

template <typename T>
std::ostream &
operator<<(std::ostream &os, const stkgroup<T> &grp) {
	typedef typename stkgroup<T>::stk_t stk_t;

	typename stkgroup<T>::const_iterator i;
	os << INK_MAGIC_STR << ' ' << grp.nstrokes() << std::endl;
	for (i = grp.begin(); i != grp.end(); ++i) {
		const stk_t *s = *i;
		os << *s << std::endl;
	}
	return os;
}

template <typename T>
std::istream &
operator>>(std::istream &is, stkgroup<T> &grp) {
	typedef typename stkgroup<T>::stk_t stk_t;

	std::string magic;
	is >> magic;
	if (magic != INK_MAGIC_STR) {
		THROW_ERROR(E_IO, "Magic ink string corrupt or missing");
	}
	size_t n;
	is >> n;
	stk_t *stks = new stk_t[n];
	stk_t *s = stks;
	try {
		while (n--) {
			is >> *s;
			grp.addstk(s++);
		}
	}
	catch (error &e) {
		delete[] stks;
		throw e;
	}

	grp.chown(stkgroup<T>::own_array);
	return is;
}


template <typename From, typename To>
int
convert_stkgroup(stkgroup<To> &dst, const stkgroup<From> &grp) {
	typedef typename stkgroup<To>::stk_t stk_t;
	dst.clear();
	stk_t *stks = new stk_t[grp.nstrokes()];
	stk_t *s = stks;
	for (typename stkgroup<From>::const_iterator i = grp.begin(); i != grp.end(); ++i) {
		int e;
		if (FAILURE(e = convert_stroke(*s, **i))) {
			ERR(geterr(), "Error converting stroke group");
			dst.clear();
			delete[] stks;
		}
		dst.addstk(s++);
	}
	dst.chown(stkgroup<To>::own_array);
	return 0;
}

template <typename T>
int
normalize(float_stkgroup &dst, const stkgroup<T> &grp, const float_rect &bounds) {
	rect<typename stkgroup<T>::pt_t> bbox = grp.bbox();
	affine_transform<float_t> xf = normalization_xform(bbox, bounds);
	int e;
	if (FAILURE(e = convert_stkgroup(dst, grp))) {
		ERR(geterr(), "Error normalizing stroke group");
		return e;
	}
	dst.apply(xf);
	return 0;
}


template <typename T>
int
subdivide(stkgroup<T> &dst, const stkgroup<T> &grp, size_t ptsperstk) {
	typedef typename stkgroup<T>::stk_t stk_t;

	stk_t *stks = new stk_t[grp.nstrokes()];
	stk_t *s = stks;
	for (typename stkgroup<T>::const_iterator i = grp.begin(); i != grp.end(); ++i) {
		int e;
		if (FAILURE(e = subdivide(*s, **i, ptsperstk))) {
			ERR(geterr(), "Error subdividing stroke group");
			delete[] stks;
			return e;
		}
		dst.addstk(s++);
	}
	dst.chown(stkgroup<T>::own_array);
	return 0;
}


template <typename T>
stkgroup<T> *
stkgroup<T>::merge(const stkgroup<T> *rhs) const {
	stkgroup<T> *&newgrp = mergecache[rhs];
	if (!newgrp) {
		newgrp = new stkgroup<T>;
		for (size_t i = 0; i < nstrokes(); ++i) {
			newgrp->addstk(stk(i));
		}
		for (size_t i = 0; i < rhs->nstrokes(); ++i) {
			if (!newgrp->hasstk(rhs->stk(i))) {
				newgrp->addstk(rhs->stk(i));
			}
		}
	}
	return newgrp;
}

template <typename T>
stkgroup<T> *
mergegrps(const stkgroup<T> *grp1, const stkgroup<T> *grp2) {
	if (!grp1) {
		return const_cast<stkgroup<T> *>(grp2);
	}
	if (!grp2) {
		return const_cast<stkgroup<T> *>(grp1);
	}
	return grp1->merge(grp2);
}

template <typename T, typename U>
double
lspscore(const stkgroup<T> &lhs, const stkgroup<U> &rhs) {
	if (lhs.nstrokes() != rhs.nstrokes()) {
		return 0.0;
	}

	double score = 1.0;
	for (size_t i = 0; i < lhs.nstrokes(); ++i) {
		score *= lspscore(*lhs.stk(i), *rhs.stk(i));
	}

	return score;
}

}
