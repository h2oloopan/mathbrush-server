/*
   algorithm.h
   
 A supplement to the STL algorithms available in <algorithm>.  This file implements
   some simple and common operations generically.
*/

#ifndef ALGORITHM_H_
#define ALGORITHM_H_


namespace scg
{


template <class C, typename T, typename V, typename Pred>
T
sorted_insert(C &container, T first, T last, const V &value, const Pred &pred)
{
	while (first != last) {
		if (pred(value, *first)) {
			break;
		}
		++first;
	}
	return container.insert(first, value);
}

template <class C, typename T, typename V, typename Pred>
inline T
sorted_insert(C &container, T first, T last, const V &value)
{ return sorted_insert(container, first, last, value, std::less<V>()); }

template <typename T, class C, typename Pred>
inline typename C::iterator
sorted_insert(C &container, const T &value, const Pred &pred)
{ return sorted_insert(container, container.begin(), container.end(), value, pred); }

template <typename T, class C>
inline typename C::iterator
sorted_insert(C &container, const T &value)
{ return sorted_insert(container, value, std::less<T>()); }

template <typename T, typename U, typename P>
U
find_binary(U first, U last, const T &val, P pred)
{
    while (first != last) {
        if (pred(*first, val)) {
            break;
        }
        ++first;
    }
    
    return first;
}

template <typename T, typename U>
U
find(U first, U last, const T &val)
{
    while (first != last) {
        if (*first == val) {
            break;
        }
        ++first;
    }
    
    return first;
}

template <typename U, typename Pred>
U
find_if(U first, U last, Pred pred)
{
    while (first != last) {
        if (pred(*first)) {
            break;
        }
        ++first;
    }
    
    return first;
}


/*
 Find the minimum element over a range by using a user-provided binary function "pred".
   pred(a,b) should return true iff a is lexicographically "less than" b.
*/
template <typename T, typename U, typename F2>
T
min(U first, U last, F2 pred)
{
    T m = *first;
    while (first != last) {
        if (pred(*first, m)) {
            m = *first;
        }
        first++;
    }
    return m;
}

/*
 Find the minimum element over a range (STL provides only the minimum of two object)
   using the standard less-than operator.
*/
template <typename T, typename U>
T
min(U first, U last)
{
    return scg::min<T, U, std::less<T> >(first, last, std::less<T>());
}


/*
 Find the maximum element over a range by using a user-provided binary function "pred".
   pred(a,b) should return true iff a is lexicographically "less than" b.
*/
template <typename T, typename U, typename F2>
T
max(U first, U last, F2 pred)
{
    T m = *first;
    while (first != last) {
        if (pred(m, *first)) {
            m = *first;
        }
        first++;
    }
    return m;
}

/*
 Find the maximum element over a range (STL provides only the maximum of two object)
   using the standard less-than operator.
*/
template <typename T, typename U>
T
max(U first, U last)
{
    return scg::max<T, U, std::less<T> >(first, last, std::less<T>());
}



template <typename T, typename U, typename F>
U
copy_using(T first, T last, U out, F f)
{
	while (first != last) {
		*(out++) = f(*(first++));
	}
	return out;
}


}


#endif

