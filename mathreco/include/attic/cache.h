#ifndef CACHE_H_
#define CACHE_H_

#include "func.h"

#include <list>

namespace scg {


class cacheable {
public:
	virtual void invalidate() = 0;
};


template <typename T, typename CalcFnT>
class cached_value : public cacheable {
private:
	mutable T val;
	mutable bool dirty;
	const CalcFnT calc_fn;

public:
	cached_value(const CalcFnT &fn);
	cached_value(const T &val_, const CalcFnT &fn);
	void invalidate();
	T operator()() const;
	cached_value &operator=(const T &rhs);

	// This assignment operator copies only the state of the cached
	// value (i.e., the value and dirty specs). The calculation
	// function is immutable after construction.
	cached_value &operator=(const cached_value &rhs);
};

class cachegrp {
	std::list<cacheable *> values;

public:
	void manage(cacheable &value);
	void invalidate();
};


}


#include "cache.tcc"

#endif
