#include "cache.h"


namespace scg {


template <typename T, typename CalcFnT>
cached_value<T, CalcFnT>::cached_value(const CalcFnT &fn)
	: calc_fn(fn), dirty(true) {
}

template <typename T, typename CalcFnT>
cached_value<T, CalcFnT>::cached_value(const T &val_, const CalcFnT &fn)
	: val(val_), calc_fn(fn), dirty(false) {
}

template <typename T, typename CalcFnT>
cached_value<T, CalcFnT> &
cached_value<T, CalcFnT>::operator=(const cached_value<T, CalcFnT> &rhs) {
	val = rhs.val;
	dirty = rhs.dirty;
	return *this;
}

template <typename T, typename CalcFnT>
cached_value<T, CalcFnT> &
cached_value<T, CalcFnT>::operator=(const T &rhs) {
	val = rhs;
	dirty = false;
}

template <typename T, typename CalcFnT>
void
cached_value<T, CalcFnT>::invalidate() {
	dirty = true;
}

template <typename T, typename CalcFnT>
T
cached_value<T, CalcFnT>::operator()() const {
	if (dirty) {
		val = calc_fn();
		dirty = false;
	}
	return val;
}


}
