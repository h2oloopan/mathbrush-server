/*
 * gauss.h
 *
 * This file defines various statistical constructs and functions, especially
 * for using and manipulating gaussian distributions.
 */


#ifndef GAUSS_H_
#define GAUSS_H_


#include <cmath>


namespace scg
{


template <typename T>
struct gaussian_t
{
    T mean;
    T variance;
};


template <typename T>
inline double
gauss_pdf(T dist_sq, T variance)
{
	static const double SQRT_2PI = M_SQRT2 * std::sqrt(M_PI);
	
	return std::exp(static_cast<double>(dist_sq) / (-2 * variance)) / (SQRT_2PI * std::sqrt(static_cast<double>(variance)));
}


template <typename T>
inline double
gauss_pdf(T mean, T variance, T x)
{
	return gauss_pdf((x - mean) * (x - mean), variance);
}


template <typename T>
inline double
gauss_pdf(const gaussian_t<T> &gauss, T x)
{
    return gauss_pdf(gauss.mean, gauss.variance, x);
}



template <typename T, typename U>
T
mean(U first, U last, T bias)
{
    T total = T(0);
    unsigned count = 0;
       
    while (first != last) {
        total += *first;
        count++;
        first++;
    }
    
    if (count == 0) {
        throw E_INVALID;
    }
    
    return bias + total / count;
}


template <typename T, typename U>
T
variance(U first, U last, T mean)
{
    T total = T(0);
    unsigned count = 0;
    
    while (first != last) {
        total += (*first - mean) * (*first - mean);
        count++;
        first++;
    }
    
    if (count < 2) {
        throw E_INVALID;
    }
    
    return total / (count - 1);
}


}


#endif

