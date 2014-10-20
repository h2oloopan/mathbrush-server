#ifndef STROKE_H_
#define STROKE_H_

#include "rect.h"
#include "lsp.h"

#include <limits>
#include <numeric>
#include <functional>
#include <algorithm>
#include <vector>

namespace scg {

template <typename T>
struct Stroke {
public:
	T *x;
	T *y;
	unsigned long *time;
	size_t npoints;
	mutable double *lspxc;
	mutable double *lspyc;

	Stroke() : x(0), y(0), time(0), npoints(0), lspxc(0), lspyc(0) { }
	Stroke(T *x_, T *y_, unsigned long *time_, size_t n_) : x(x_), y(y_), time(time_), npoints(n_), lspxc(0), lspyc(0) { }
	~Stroke() { clear(); }

	void clear() {
		delete[] x;
		delete[] y;
		delete[] time;
		delete[] lspxc;
		delete[] lspyc;
		npoints = 0;
		x = y = 0;
		time = 0;
		lspxc = lspyc = 0;
	}

   Stroke &operator=(const Stroke &s) {
		if (&s != this) {
		  	clear();
			x = s.x;
			y = s.y;
			time = s.time;
			npoints = s.npoints;
			lspxc = s.lspxc;
			lspyc = s.lspyc;
			
			Stroke &ss = (Stroke &)s;
			ss.x = ss.y = 0;
			ss.time = 0;
			ss.npoints = 0;
			ss.lspxc = ss.lspyc = 0;
		}
		return *this;
	}

    void set_points(T *x_, T *y_, size_t n) {
        clear();
        this->x = x_;
        this->y = y_;
        this->npoints = n;
    }
    
    void set_points(T *x_, T *y_, unsigned long *t, size_t n)
    {
        clear();
        this->x = x_;
        this->y = y_;
        this->time = t;
        this->npoints = n;
    }
	
	void ensurelsp() const {
		if (!lspxc) {
			lspxc = new double[LSPDEG+1];
			lspyc = new double[LSPDEG+1];
			mklsp();
		}
	}

	inline void reset() {
		x = 0;
		y = 0;
		time = 0;
		npoints = 0;
	}

	Stroke<T> &reverse() {
		T *xa, *xb;
		T *ya, *yb;
		xa = x;
		xb = x + npoints - 1;
		ya = y;
		yb = y + npoints - 1;
		while (xa < xb) {
			std::swap(*xa, *xb);
			std::swap(*ya, *yb);
			++xa;
			++ya;
			--xb;
			--yb;
		}
		clearlsp();
		return *this;
	}

   Stroke<T> substroke(size_t start, size_t length = (size_t)-1) const {
		if (length == static_cast<size_t>(-1)) {
			length = this->npoints - start;
		}
		T *newx = new T[length];
		T *newy = new T[length];
		unsigned long *newtime = 0;
		if (time) {
		 newtime = new unsigned long[length];
		} 
		std::copy(x + start, x + start + length, newx);
		std::copy(y + start, y + start + length, newy);
		if (time) {
		 std::copy(time + start, time + start + length, newtime);
		}        
		return Stroke<T>(newx, newy, newtime, length);
	}

	Rect<T> bounds() const {
		Rect<T> bbox;
		if (std::numeric_limits<T>::has_infinity) {
			bbox.left = bbox.top = std::numeric_limits<T>::infinity();
			bbox.right = bbox.bottom = -std::numeric_limits<T>::infinity();
		}
		else {
			bbox.left = bbox.top = std::numeric_limits<T>::max();
			bbox.right = bbox.bottom = std::numeric_limits<T>::min();
		}
		for (size_t i = 0; i < npoints; ++i) {
			bbox.left = std::min(bbox.left, x[i]);
			bbox.right = std::max(bbox.right, x[i]);
			bbox.top = std::min(bbox.top, y[i]);
			bbox.bottom = std::max(bbox.bottom, y[i]);
		}
		return bbox;
	}

	T x_center() const {
		return std::accumulate(x, x + npoints, T(0)) / npoints;
	}
	T y_center() const {
		return std::accumulate(y, y + npoints, T(0)) / npoints;
	}

	Stroke<T> &translate(T dx, T dy) {
		for (T *px = x; px != x + npoints; ++px) {
			*px += dx;
		}
		for (T *py = y; py != y + npoints; ++py) {
			*py += dy;
		}
		return *this;
	}

	Stroke<T> copy() const {
		T *cpx = new T[npoints];
		T *cpy = new T[npoints];
		unsigned long *cpt = time ? new unsigned long[npoints] : 0;
		std::copy(x, x + npoints, cpx);
		std::copy(y, y + npoints, cpy);
		if (time) {
			std::copy(time, time + npoints, cpt);
		}
		return Stroke<T>(cpx, cpy, cpt, npoints);
	}

	double length() const {
		if (npoints == 0) return 0;
		double L = 0;
		for (size_t i = 0; i < npoints-1; ++i) {
			T dx = x[i+1] - x[i];
			T dy = y[i+1] - y[i];
			L += std::sqrt((double)dx*dx + dy*dy);
		}
		return L;
	}

private:
	void mklsp() const {
		Rect<T> bbox = bounds();
		double sc = 1.0 / std::max(bbox.width(), bbox.height());
		std::fill(lspxc, lspxc + LSPDEG + 1, 0);
		std::fill(lspyc, lspyc + LSPDEG + 1, 0);
		std::vector<double> xacc(LSPDEG+1, 0);
		std::vector<double> yacc(LSPDEG+1, 0);
		double L = 0.0;
		for (size_t i = 0; i < npoints - 1; ++i) {
			double dx = sc*(x[i+1] - x[i]);
			double dy = sc*(y[i+1] - y[i]);
			double sx = sc*0.5*(x[i+1] + x[i] - 2*bbox.left);
			double sy = sc*0.5*(y[i+1] + y[i] - 2*bbox.top);
			double L0 = std::sqrt(dx*dx + dy*dy);
			double L1 = L;
			double L1k = L;
			L += L0;
			double L2k = L;
			for (size_t k = 0; k <= LSPDEG; ++k) {
				double al = (L2k - L1k) / (k+1);
				xacc[k] += al * sx;
				yacc[k] += al * sy;
				L1k *= L1;
				L2k *= L;
			}
		}

		double Lk = L;
		for (size_t j = 0; j <= LSPDEG; ++j) {
			xacc[j] /= Lk;
			yacc[j] /= Lk;
			Lk *= L;
		}

		double n = 0.0;
		const static double u = 0.125;
		for (size_t k = 1; k <= LSPDEG; ++k) {
			for (size_t m = 0; m <= k; ++m) {
				lspxc[k] += B[k][m] * xacc[m];
				lspyc[k] += B[k][m] * yacc[m];
				if (m > 0) {
					double c = u * B[k][m] * m;
					lspxc[k] += c * sc * (x[npoints-1] - bbox.left);
					lspyc[k] += c * sc * (y[npoints-1] - bbox.top);
					if (m > 1) {
						c *= (m-1);
						lspxc[k] -= c * xacc[m-2];
						lspyc[k] -= c * yacc[m-2];
					}
				}
			}
			lspxc[k] -= u * B[k][1] * sc * (x[0] - bbox.left);
			lspyc[k] -= u * B[k][1] * sc * (y[0] - bbox.top);
			n += lspxc[k]*lspxc[k] + lspyc[k]*lspyc[k];
		}

		n = std::sqrt(n);
		for (size_t k = 1; k <= LSPDEG; ++k) {
			lspxc[k] /= n;
			lspyc[k] /= n;
		}

		//lspxc[0] = x[0];
		//lspyc[0] = y[0];
	}

	void clearlsp() const {
		delete[] lspxc;
		delete[] lspyc;
		lspxc = lspyc = 0;
	}
};


typedef Stroke<long> RawStroke;
typedef Stroke<double> NormalizedStroke;

template <typename T>
Stroke<T>
concat(const Stroke<T> &s1, const Stroke<T> &s2) {
	T *x = new T[s1.npoints + s2.npoints];
	T *y = new T[s1.npoints + s2.npoints];
	unsigned long *t = (s1.time && s2.time) ? new unsigned long[s1.npoints + s2.npoints] : 0;
	std::copy(s1.x, s1.x + s1.npoints, x);
	std::copy(s1.y, s1.y + s1.npoints, y);
	if (t) {
		std::copy(s1.time, s1.time + s1.npoints, t);
	}
	std::copy(s2.x, s2.x + s2.npoints, x + s1.npoints);
	std::copy(s2.y, s2.y + s2.npoints, y + s1.npoints);
	if (t) {
		std::copy(s2.time, s2.time + s2.npoints, t + s1.npoints);
	}
	return Stroke<T>(x, y, t, s1.npoints + s2.npoints);
}

template <typename T>
double
lspdist(const Stroke<T> &s1, const Stroke<T> &s2) {
	s1.ensurelsp();
	s2.ensurelsp();
	double d = 0;
	for (size_t i = 0; i <= LSPDEG; ++i) {
		double dx = s1.lspxc[i] - s2.lspxc[i];
		double dy = s1.lspyc[i] - s2.lspyc[i];
		d += dx*dx+dy*dy;
	}
	return std::sqrt(d);
}

template <typename T>
double
lspscore(const Stroke<T> &s1, const Stroke<T> &s2) {
	return 1.0 - lspdist(s1, s2)/M_SQRT2;
}

}


#endif
