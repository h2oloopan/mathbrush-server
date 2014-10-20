#include "error.h"
#include <cassert>
#include <algorithm>
#include <iostream>

namespace scg {

template <typename T>
interpreter<T>::~interpreter() {
}

template <typename T>
void
interpreter<T>::reset() {
	localreset();
	for (std::set<interpreter<T> *>::iterator i = parents.begin(); i != parents.end(); ++i) {
		(*i)->reset();
	}
}

template <typename T>
void
interpreter<T>::addparent(interpreter<T> *p) {
	parents.insert(p);
}

template <typename T>
void
interpreter<T>::rmparent(interpreter<T> *p) {
	std::set<interpreter<T> *>::iterator i = parents.find(p);
	if (i != parents.end()) {
		parents.erase(i);
	}
}

template <typename T>
staticintrpr<T>::staticintrpr()
	: at(0) {
}

template <typename T>
staticintrpr<T>::staticintrpr(T *t, bool rm)
	: at(0) {
	addknown(t, rm);
}

template <typename T>
staticintrpr<T>::~staticintrpr() {
	typename std::vector<T *>::iterator i;
	for (size_t i = 0; i < ts.size(); ++i) {
		if (rmflags[i]) {
			delete ts[i];
		}
	}
}

template <typename T>
size_t
staticintrpr<T>::nknown() const {
	return at;
}

template <typename T>
T *
staticintrpr<T>::nth(size_t i) const {
	return (i < ts.size()) ? ts[i] : 0;
}

template <typename T>
size_t
staticintrpr<T>::nadded() const {
	return ts.size();
}

template <typename T>
T *
staticintrpr<T>::next() {
	if (at == ts.size()) {
		return 0;
	}
	else {
		return ts[at++];
	}
}

template <typename T>
void
staticintrpr<T>::localreset() {
	at = 0;
}

template <typename T>
bool
cmpscored(const T *lhs, const T *rhs) {
	return lhs->score() > rhs->score();
}

template <typename T>
void
staticintrpr<T>::addknown(T *t, bool rm) {
	typename std::vector<T *>::iterator i;
	i = std::upper_bound(ts.begin(), ts.end(), t, &cmpscored<T>);
	size_t n = i - ts.begin();
	ts.insert(i, t);
	rmflags.insert(rmflags.begin() + n, rm);
}


template <typename T>
void
parser<T>::freeknown() {
	typename std::vector<T *>::iterator i;
	for (i = known.begin(); i != known.end(); ++i) {
		delete *i;
	}
}

template <typename T>
size_t
parser<T>::nknown() const {
	return known.size();
}

template <typename T>
T *
parser<T>::nth(size_t i) const {
	if (i >= known.size()) {
		THROW_ERROR(E_INVALID, "Interpretation index " << i << " is too large for parser " << this);
	}
	return known[i];
}

template <typename T>
T *
parser<T>::next() {
	T *t = getnext();
	if (!t) {
		return 0;
	}

	typename std::vector<T *>::iterator i;
	known.push_back(t);
	return t;
}

template <typename T>
void
parser<T>::localreset() {
	known.clear();
}

template <typename T>
iterator<T>::iterator(interpreter<T> *src_, bool ownptr)
	: src(src_), curr(0), own(ownptr) {
	assert(src);
	src->addparent(this);
	//DEBUG_ONLY(std::cerr << "iterator: creating " << this << " for " << src << std::endl);
}

template <typename T>
iterator<T>::~iterator() {
	if (own) {
		delete src;
	}
	else {
		src->rmparent(this);
	}
}

template <typename T>
size_t
iterator<T>::nknown() const {
	return curr;
}

template <typename T>
T *
iterator<T>::nth(size_t i) const {
	return i < curr ? src->nth(i) : 0;
}

template <typename T>
T *
iterator<T>::getnext() {
	if (curr < src->nknown()) {
		return src->nth(curr++);
	}
	T *nx = src->next();
	if (nx) {
		++curr;
	}
	return nx;
}

template <typename T>
void
iterator<T>::localreset() {
	parser<T>::localreset();
	curr = 0;
	src->reset();
}

template <typename T, typename I>
indirect_iterator<T, I>::indirect_iterator(const I &src_, bool ownptr)
	: src(src_), curr(0), own(ownptr) {
	assert(src);
	(*src)->addparent(this);
	//DEBUG_ONLY(std::cerr << "iterator: creating " << this << " for " << src << std::endl);
}

template <typename T, typename I>
indirect_iterator<T, I>::~indirect_iterator() {
	if (own) {
		delete src;
	}
	else {
		(*src)->rmparent(this);
	}	
}

template <typename T, typename I>
size_t
indirect_iterator<T, I>::nknown() const {
	return curr;
}

template <typename T, typename I>
T *
indirect_iterator<T, I>::nth(size_t i) const {
	return i < curr ? (*src)->nth(i) : 0;
}

template <typename T, typename I>
T *
indirect_iterator<T, I>::getnext() {
	if (curr < (*src)->nknown()) {
		return (*src)->nth(curr++);
	}
	T *nx = (*src)->next();
	if (nx) {
		++curr;
	}
	return nx;
}

template <typename T, typename I>
void
indirect_iterator<T, I>::localreset() {
	parser<T>::localreset();
	curr = 0;
}

template <typename T>
T *
getnth(interpreter<T> *intrpr, size_t n) {
	if (intrpr->nknown() <= n) {
		T *t;
		do {
			t = intrpr->next();
			if (!t) {
				return 0;
			}
		} while (intrpr->nknown() <= n);
		return t;
	}
	return intrpr->nth(n);
}


template <typename T>
T *
getfirst(interpreter<T> *intrpr) {
	return (intrpr->nknown() > 0) ? intrpr->nth(0) : intrpr->next();
}

template <typename T>
multiplexor<T>::multiplexor()
	: parser<T>(), prevsrc(0) {
}

template <typename T>
multiplexor<T>::~multiplexor() {
	for (std::vector<intrpr_t>::iterator i = intrprs.begin(); i != intrprs.end(); ++i) {
		delete i->intrpr;
	}
}

template <typename T>
void
multiplexor<T>::addparser(interpreter<T> *intrpr) {
	intrpr->addparent(this);
	T *intrp = getfirst(intrpr);
	if (intrp) {
		intrprs.push_back(intrpr_t(intrpr));
		addintrp(intrp, intrprs.size());
	}
}

template <typename T>
size_t
multiplexor<T>::nparsers() const {
	return intrprs.size();
}

template <typename T>
T *
multiplexor<T>::getnext() {
	if (prevsrc) {
		intrpr_t &intrpr = intrprs[prevsrc - 1];
		++intrpr.i;
		T *intrp = getnth(intrpr.intrpr, intrpr.i);
		if (intrp) {
			addintrp(intrp, prevsrc);
		}
	}
	intrp_t best(0, 0);
	if (!known.empty()) {
		best = known.back();
		known.pop_back();
		prevsrc = best.src;
	}
	return best.intrp;
}

template <typename T>
void
multiplexor<T>::localreset() {
	parser<T>::localreset();
	prevsrc = 0;
	known.clear();
	for (size_t i = 0; i < intrprs.size(); ++i) {
		intrprs[i].i = 0;
		T *intrp = getfirst(intrprs[i].intrpr);
		if (intrp) {
			addintrp(intrp, i + 1);
		}
	}
}

template <typename T>
multiplexor<T>::intrpr_t::intrpr_t(interpreter<T> *intrpr_) : intrpr(intrpr_), i(0) { }

template <typename T>
multiplexor<T>::intrp_t::intrp_t(T *intrp_, size_t src_)
	: intrp(intrp_), src(src_) {
}

template <typename T>
bool
multiplexor<T>::intrp_t::operator<(const typename multiplexor::intrp_t &rhs) const {
	return intrp->score() < rhs.intrp->score();
}

template <typename T>
void
multiplexor<T>::addintrp(T *intrp, size_t src) {
	intrp_t val(intrp, src);
	typename std::vector<intrp_t>::iterator i;
	i = std::upper_bound(known.begin(), known.end(), val);
	known.insert(i, val);
}

}
