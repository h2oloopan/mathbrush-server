#ifndef VECTOR_IO_H_
#define VECTOR_IO_H_


#include <istream>
#include <ostream>
#include <vector>

#include "error.h"


// need these in both the scg and root namespaces to be useful both
// within and without the scg namespace


template <typename T>
std::istream &
operator>>(std::istream &is, std::vector<T> &vec)
{
    unsigned count;
    
    is >> count;
    
    while (count--) {
        if (!is) {
            throw E_IO;
        }
        
        vec.push_back(T());
        is >> vec.back();
    }
    
    return is;
}


template <typename T>
std::ostream &
operator<<(std::ostream &os, const std::vector<T> &vec)
{
    os << static_cast<unsigned>(vec.size()) << std::endl;
    
    for (typename std::vector<T>::const_iterator i = vec.begin(); i != vec.end(); ++i) {
        os << *i << std::endl;
    }
    
    return os;
}


namespace scg
{

template <typename T>
std::istream &
operator>>(std::istream &is, std::vector<T> &vec)
{
    unsigned count;
    
    is >> count;
    
    while (count--) {
        if (!is) {
            throw E_IO;
        }
        
        vec.push_back(T());
        is >> vec.back();
    }
    
    return is;
}


template <typename T>
std::ostream &
operator<<(std::ostream &os, const std::vector<T> &vec)
{
    os << static_cast<unsigned>(vec.size()) << std::endl;
    
    for (typename std::vector<T>::const_iterator i = vec.begin(); i != vec.end(); ++i) {
        os << *i << std::endl;
    }
    
    return os;
}


}


#endif

