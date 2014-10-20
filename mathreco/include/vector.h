/*
 * vector.h
 *
 * This file defines helper methods to be used with the STL
 * 'vector' container.
 *
 * In particular, it implements generic sorted insertion assuming
 * a lexicographical ordering by the < operator.
 *
 */

#ifndef VECTOR_H_
#define VECTOR_H_


#include <vector>

namespace scg
{


template <typename T>
int
sorted_insert(const T& val, std::vector<T> &vector)
{
    typename std::vector<T>::iterator i;
    i = std::find(vector.begin(), vector.end(), val);
    if (i != vector.end()) {
        if (*i < val) {
            return 0;
        }
        vector.erase(i);
    }
    
    for (i = vector.begin(); i != vector.end() && *i < val; i++);
    
    vector.insert(i, val);
        
    return 0;
}


template <typename T>
int
sorted_index_insert(const T &val, std::vector<T> &vector)
{
	typename std::vector<T>::iterator i;
	for (i = vector.begin(); i != vector.end(); ++i) {
		if (val <= *i) {
			i = vector.insert(i, val);
			++i;
			break;
		}
	}

	if (i == vector.end()) {
		vector.insert(vector.end(), val);
		return 0;
	}

	for (; i != vector.end(); ++i) {
		++(*i);
	}

	return 0;
}


}


#endif

