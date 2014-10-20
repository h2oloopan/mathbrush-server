#ifndef UTILS_H_
#define UTILS_H_


#include <iostream>
#include <utility>


namespace scg
{


template <typename T1, typename T2>
bool
subset(T1 startA, T1 endA, T2 startB, T2 endB)
{
	while (startA != endA) {
		if (std::find(startB, endB, *startA) == endB) {
			return false;
		}
		++startA;
	}

	return true;
}


template <typename T1, typename V, typename F, typename P>
std::pair<T1, V>
fn_first_by(T1 start, T1 end, V base, const F &f, const P &pred)
{
	T1 minelem = end;
	while (start != end) {
		V val = f(*start);
		if (pred(val, base)) {
			base = val;
			minelem = start;
		}
		++start;
	}
	return std::make_pair(minelem, base);
}


template <typename T1, typename V, typename F>
inline std::pair<T1, V>
minimizer(T1 start, T1 end, V base, const F &f)
{
	return fn_first_by<T1, V, F, std::less<V> >(start, end, base, f, std::less<V>());
}

template <typename T1, typename V, typename F>
inline std::pair<T1, V>
maximizer(T1 start, T1 end, V base, const F &f)
{
	return fn_first_by<T1, V, F, std::greater<V> >(start, end, base, f, std::greater<V>());
}

template <typename T1, typename P>
T1
first_by(T1 start, T1 end, const P &pred)
{
	if (start == end) {
		return end;
	}

	T1 max_elem = start++;
	while (start != end) {
		if (pred(*start, *max_elem)) {
			max_elem = start;
		}
		++start;
	}
	return max_elem;
}

template <typename V, typename T1, typename P>
V
first_by(T1 start, T1 end, V base, const P &pred)
{
	if (start == end) {
		return base;
	}

	V max_elem = *(start++);
	while (start != end) {
		if (pred(*start, max_elem)) {
			max_elem = *start;
		}
		++start;
	}
	return max_elem;
}

template <typename V, typename T1>
V
range_min(T1 start, T1 end, V base)
{ return first_by(start, end, base, std::less<V>()); }


template <typename V, typename T1>
V
range_max(T1 start, T1 end, V base)
{ return first_by(start, end, base, std::greater<V>()); }


}


#endif
