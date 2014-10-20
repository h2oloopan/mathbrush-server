#include "functions.h"
#include <sstream>
#include <limits>

namespace scg
{

// ----- Linear Function -----
int BoundedLinearFunction::evaluate(double& result, int x) const
{
	if (x < this->lowerBound || x > this->upperBound) return E_INVALID;

	result = (double) m * x + b;
	return 0;
}

// ----- Polynomial Function -----
double BoundedPolynomialFunction::getCoefficient(int exponent) const
{
	// If the term with the specified exponent exists, return its coefficient
	for (std::vector<PolyTerm>::const_iterator it = terms.begin(); it != terms.end(); it++) {
		const PolyTerm& term = *it;
		if (term.exponent == exponent) {
			return term.coefficient;
		}
	}

	return 0;
}

std::string BoundedPolynomialFunction::str() const
{
	std::stringstream ss;
	for (size_t i = terms.size(); i > 0; i--) {
		const PolyTerm term = terms[i-1];
		if (term.coefficient == 0) continue;
		else if (term.coefficient != 1) {
			if (term.coefficient > 0 && i != terms.size()) {
				ss << "+";
			}

			ss << term.coefficient;
		}

		if (term.exponent == 1) ss << "x";
		else if (term.exponent > 0) ss << "x^" << term.exponent;
	}

	std::string out = ss.str();
	return out;
}

void BoundedPolynomialFunction::addTerm(double coefficient, int exponent)
{
	PolyTerm term(coefficient, exponent);
	addTerm(term);
}

void BoundedPolynomialFunction::addTerm(PolyTerm term)
{
	if (term.coefficient == 0) return;

	// Find where the term belongs
	unsigned i;
	for (i = 0; i < terms.size() && terms[i].exponent < term.exponent; i++) {}
	
	// Either insert the term or add its coefficient to the existing term
	if (i < terms.size() && terms[i].exponent == term.exponent) {
		PolyTerm& existingTerm = terms[i];
		existingTerm.coefficient += term.coefficient;
		if (existingTerm.coefficient == 0) terms.erase(terms.begin() + i);
	} else {
		terms.insert(terms.begin() + i, term);
	}
}

int BoundedPolynomialFunction::evaluate(double& result, int x) const
{
	if (x < this->lowerBound || x > this->upperBound) return E_INVALID;

	result = 0.0;
	// Add the evaluation of each term
	for (std::vector<PolyTerm>::const_iterator it = terms.begin(); it != terms.end(); it++) {
		result += (*it)(x);
	}

	return 0;
}

BoundedPolynomialFunction BoundedPolynomialFunction::operator* (const BoundedPolynomialFunction& other) const
{
	BoundedPolynomialFunction ret;

	// TODO O(n^2) for now -- do we need better? (presumably dealing with small polynomials)
	for (std::vector<PolyTerm>::const_iterator terms1 = terms.begin(); terms1 != terms.end(); terms1++) {
		for (std::vector<PolyTerm>::const_iterator terms2 = other.terms.begin(); terms2 != other.terms.end(); terms2++) {
			PolyTerm newTerm(terms1->coefficient * terms2->coefficient, terms1->exponent + terms2->exponent);
			ret.addTerm(newTerm);
		}
	}

	return ret;
}

BoundedPolynomialFunction operator* (double a, const BoundedPolynomialFunction& other)
{
	BoundedPolynomialFunction ret = other;
	for (std::vector<PolyTerm>::iterator it = ret.terms.begin(); it != ret.terms.end(); it++) {
		it->coefficient *= a;
	}

	return ret;
}

BoundedPolynomialFunction BoundedPolynomialFunction::operator* (double a) const
{
	return a * (*this);
}

BoundedPolynomialFunction BoundedPolynomialFunction::operator/ (double a) const
{
	return (*this) * (1.0 / a);
}

BoundedPolynomialFunction BoundedPolynomialFunction::operator+ (const BoundedPolynomialFunction& other) const
{
	BoundedPolynomialFunction ret;

	unsigned i = 0, j = 0;
	// Add the coefficients of the terms with the same exponent value
	while (i < terms.size() && j < other.terms.size()) {
		PolyTerm thisTerm = terms[i];
		PolyTerm otherTerm = other.terms[j];
		
		if (thisTerm.exponent == otherTerm.exponent) {
			PolyTerm newTerm(thisTerm.coefficient + otherTerm.coefficient, thisTerm.exponent);
			ret.addTerm(newTerm);
			i++; j++;
		} else if (thisTerm.exponent < otherTerm.exponent) {
			ret.addTerm(thisTerm);
			i++;
		} else if (thisTerm.exponent > otherTerm.exponent) {
			ret.addTerm(otherTerm);
			j++;
		}
	}

	// Finish off any unvisited terms
	for (; i < terms.size(); i++) {
		ret.addTerm(terms[i]);
	}
	for (; j < other.terms.size(); j++) {
		ret.addTerm(other.terms[j]);
	}

	return ret;
}

BoundedPolynomialFunction BoundedPolynomialFunction::operator+ (int a) const
{
	BoundedPolynomialFunction ret = *this;
	ret.addTerm(a, 0);
	return ret;
}

BoundedPolynomialFunction BoundedPolynomialFunction::operator- (const BoundedPolynomialFunction& other) const
{
	return this->operator+(-1*other);
}

BoundedPolynomialFunction BoundedPolynomialFunction::operator- (int a) const
{
	return (*this) + (-1 * a);
}


BoundedPolynomialFunction createLagrangePolynomial(std::vector<std::pair<int,int> > points)
{
	BoundedPolynomialFunction f(0, std::numeric_limits<int>::max());

	// Calculate the basis polynomials
	BoundedPolynomialFunction* basisPoly = new BoundedPolynomialFunction[points.size()];
	for (unsigned i = 0; i < points.size(); i++) {
		for (unsigned j = 0; j < points.size(); j++) {
			if (i == j) continue;

			BoundedPolynomialFunction polynomial;

			// Construct the polynomial of the form (x - x_j)/(x_i - x_j)
			polynomial.addTerm(1, 1);
			polynomial = polynomial - points[j].first;
			polynomial = polynomial / (points[i].first - points[j].first);

			if (basisPoly[i].getNumTerms() == 0) {
				basisPoly[i] = polynomial;
			} else {
				basisPoly[i] = basisPoly[i] * polynomial;
			}
		}
	}

	// Calculate the interpolating polynomial
	for (unsigned i = 0; i < points.size(); i++) {
		BoundedPolynomialFunction poly = points[i].second * basisPoly[i];
		f = f + poly;
	}

	delete [] basisPoly;

	return f;
}



}
