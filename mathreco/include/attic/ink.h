#ifndef INK_H_
#define INK_H_

#include "inktype.h"
#include "cache.h"
#include "rect.h"

#include <vector>
#include <map>
#include <istream>
#include <ostream>


namespace scg {

template <typename T> struct affine_transform;

template <typename T>
struct point {
	T x;
	T y;

	point();
	point(T x_, T y_);
	template <typename U>
	point(const point<U> &pt);
	
	bool operator==(const point &rhs) const;
	bool operator!=(const point &rhs) const;

	struct lexicographic_cmp {
		bool operator()(const point &lhs, const point &rhs) const;
	};

	point &operator+=(const point &rhs);
	point &operator-=(const point &rhs);
	template <typename U>
	point &operator*=(const U a);
	template <typename U>
	point &operator/=(const U a);
	point operator+(const point &rhs) const;
	point operator-(const point &rhs) const;
	template <typename U>
	point operator/(const U a) const;

	point<T> &apply(const affine_transform<T> &xf);
};

template <typename T, typename U>
point<U> operator*(T lhs, const point<U> &rhs);

template <typename T>
point<T> mkpoint(T x, T y);

template <typename From, typename To>
int convert_point(point<To> &dst, const point<From> &pt);

template <typename T>
std::ostream &operator<<(std::ostream &os, const point<T> &pt);
template <typename T>
std::istream &operator>>(std::istream &is, point<T> &pt);


template <typename T>
struct affine_transform {
	T matrix[6];

	affine_transform();

	affine_transform &chain(const affine_transform<T> &next);

	affine_transform &scale(T x, T y);
	affine_transform &translate(T x, T y);
};

template <typename T>
affine_transform<T>
compose(const affine_transform<T> &first,
        const affine_transform<T> &second);


template <typename T> struct stkgroup;

template <typename T>
struct stroke {
public:
	typedef T pt_t;

public:
	cachegrp cache;
	void init_cache();

	rect<T> compute_bbox() const;

	stkgroup<stroke> *grp;
	mutable stkgroup<const stroke> *constgrp;

	size_t id_;

	std::vector<point<double> > lsp;

public:
	void mkcoeffs();

public:
	typedef std::vector<point<T> > point_collection;
	point_collection points;

	cached_value<rect<T>, const_bound_mem_fun_t<rect<T>, stroke<T> > > bbox;
	
public:
	stroke();
	template <typename InIt>
	stroke(InIt first, InIt last);
	stroke(const stroke &rhs);
	~stroke();

	size_t npoints() const;
	void clear();

	size_t id() const;

	stroke &apply(const affine_transform<T> &xf);

	stkgroup<stroke> *selfgrp();
	stkgroup<const stroke> *const_selfgrp() const;

	void update();

	const std::vector<point<double> > &lspvec() const { return lsp; }

private:
	template <typename S, typename U>
	friend double lspdist(const stroke<S> &, const stroke<U> &);
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const stroke<T> &stk);
template <typename T>
std::istream &operator>>(std::istream &is, stroke<T> &stk);

typedef stroke<int_t> int_stroke;
typedef stroke<float_t> float_stroke;

template <typename From, typename To>
int convert_stroke(stroke<To> &dst, const stroke<From> &stk);
template <typename T>
int subdivide(stroke<T> &dst, const stroke<T> &stk, size_t npoints);
template <typename T>
int normalize(float_stroke &dst, const stroke<T> &stk, const float_rect &bounds);

template <typename T, typename U>
double lspdist(const stroke<T> &lhs, const stroke<U> &rhs);
template <typename T, typename U>
double lspscore(const stroke<T> &lhs, const stroke<U> &rhs);

template <typename T>
struct stkgroup {
public:
	typedef T stk_t;
	typedef typename T::pt_t pt_t;

public:
	enum ownstat {
		own_weak,
		own_ptrs,
		own_array
	};

private:
	ownstat ownshp;

	cachegrp cache;
	void init_cache();

	rect<pt_t> compute_bbox() const;
	size_t compute_npoints() const;

private:
	typedef std::vector<T *> stroke_collection;
	stroke_collection strokes;

private:
	typedef typename stroke_collection::iterator iterator;
public:
	typedef typename stroke_collection::const_iterator const_iterator;

	const_iterator begin() const;
	const_iterator end() const;

	int addstk(T *stk);
	int insstk(T *stk, size_t pos);
	int rmstk(T *stk);
	bool hasstk(T *stk) const;
	T *stk(size_t i);
	const T *stk(size_t i) const;

public:
	cached_value<rect<pt_t>, const_bound_mem_fun_t<rect<pt_t>, stkgroup<T> > > bbox;
	cached_value<size_t, const_bound_mem_fun_t<size_t, stkgroup<T> > > npoints;

public:
	stkgroup();
	stkgroup(const stkgroup<T> &rhs);
	~stkgroup();

	size_t nstrokes() const;
	void clear();

	stkgroup &apply(const affine_transform<pt_t> &xf);

	void chown(ownstat st);

	stkgroup *merge(const stkgroup *rhs) const;

private:
	mutable std::map<const stkgroup *, stkgroup *> mergecache;
};


typedef stkgroup<stroke<int_t> > int_stkgroup;
typedef stkgroup<stroke<float_t> > float_stkgroup;
typedef stkgroup<const stroke<int_t> > const_int_stkgroup;
typedef stkgroup<const stroke<float_t> > const_float_stkgroup;


template <typename T>
std::ostream &operator<<(std::ostream &os, const stkgroup<T> &grp);
template <typename T>
std::istream &operator>>(std::istream &is, stkgroup<T> &grp);

template <typename From, typename To>
int convert_stkgroup(stkgroup<To> &dst, const stkgroup<From> &grp);

template <typename From>
int normalize(float_stkgroup &dst, const stkgroup<From> &grp,
              const float_rect &bounds);

template <typename From, typename To>
int subdivide(stkgroup<To> &dst, const stkgroup<From> &grp,
							size_t ptsperstk);

template <typename T>
stkgroup<T> *mergegrps(const stkgroup<T> *grp1, const stkgroup<T> *grp2);

template <typename T, typename U>
double lspscore(const stkgroup<T> &lhs, const stkgroup<U> &rhs);

}


#include "ink.tcc"

#endif
