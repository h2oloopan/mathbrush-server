#ifndef SCG_MATH_FUNCTION_H_
#define SCG_MATH_FUNCTION_H_

#include <utility>
#include <vector>
#include <utility>
#include <math.h>
#include "error.h"

namespace scg
{


class BoundedUnivariateFunction {
protected:
	int lowerBound, upperBound; // defines the domain of the function

public:
	BoundedUnivariateFunction(int lowerBound, int upperBound) : lowerBound(lowerBound), upperBound(upperBound) {}

	// Get the function variable's domain as a pair of (lower, higher)
	std::pair<int,int> getDomain() const { return std::pair<int,int>(lowerBound, upperBound); }
	void setDomain(int lowerBound, int upperBound) { this->lowerBound = lowerBound; this->upperBound = upperBound; }

	// Output the value of f(x) to "result". Returns E_INVALID if x is not in the domain
	virtual int evaluate(double& result, int x) const = 0;

	virtual BoundedUnivariateFunction* copy() = 0;
};


// Define a bounded univariate function of the form f(x) = m*x + b
class BoundedLinearFunction : public BoundedUnivariateFunction {
	int m, b; // slope and y intercept in the equation f(x) = m*x + b

public:
	/* m - slope
	 * b - y intercept
	 * lowerBound - lower bound of the function
	 * upperBound - upper bound of the function
	 */
	BoundedLinearFunction(int m, int b, int lowerBound, int upperBound) : BoundedUnivariateFunction(lowerBound, upperBound), m(m), b(b) {}
	int evaluate(double& result, int x) const;

	BoundedUnivariateFunction* copy() { return new BoundedLinearFunction(m, b, lowerBound, upperBound); }
};

// A polynomial term defined by a coefficient and an exponent
struct PolyTerm {
	double coefficient;
	int exponent;

	PolyTerm(double coefficient = 0.0, int exponent = 0) : coefficient(coefficient), exponent(exponent) {}
	
	inline double operator() (double x) const { return coefficient * pow(x, (double) exponent); }
	inline bool operator < (const PolyTerm& other) const { return exponent < other.exponent; } 
	inline bool operator > (const PolyTerm& other) const { return exponent > other.exponent; } 
	inline bool operator == (const PolyTerm& other) const { return exponent == other.exponent; } 
};

// Defines a bounded polynomial function of arbitrary degree
class BoundedPolynomialFunction : public BoundedUnivariateFunction {
	std::vector<PolyTerm> terms;

public:
	explicit BoundedPolynomialFunction(int lowerBound = 0, int upperBound = 0) : BoundedUnivariateFunction(lowerBound, upperBound) {}

	size_t getNumTerms() const { return terms.size(); }
	unsigned getDegree() const { return (terms.size() > 0) ? terms.back().exponent : 0; }
	double getCoefficient(int exponent) const;
	std::string str() const;

	/* Add a term to the polynomial with a coefficient and exponent e.g. 2*x^3 has coefficient 2 and exponent 3
	 * If a term with the specified exponent already exists, their coefficients will be added
	 */
	void addTerm(double coefficient, int exponent);
	void addTerm(PolyTerm term);
	int evaluate(double& result, int x) const;

	BoundedUnivariateFunction* copy()
	{
		BoundedPolynomialFunction* copy = new BoundedPolynomialFunction(lowerBound, upperBound);
		copy->terms = terms;
		return copy;
	}

	BoundedPolynomialFunction operator* (const BoundedPolynomialFunction& other) const;
	BoundedPolynomialFunction operator* (double a) const;
	BoundedPolynomialFunction operator/ (double a) const;
	BoundedPolynomialFunction operator+ (const BoundedPolynomialFunction& other) const;
	BoundedPolynomialFunction operator+ (int a) const;
	BoundedPolynomialFunction operator- (const BoundedPolynomialFunction& other) const;
	BoundedPolynomialFunction operator- (int a) const;
	friend BoundedPolynomialFunction operator* (double a, const BoundedPolynomialFunction& other);
};



BoundedPolynomialFunction createLagrangePolynomial(std::vector<std::pair<int,int> > points);


}

#endif