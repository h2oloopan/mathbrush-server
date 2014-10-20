#include "matrix.h"

#include <sstream>
#include "interval.h"
#include "MathRecoTypes.h"
#include "MathRecognizer.h"
#include "error.h"

namespace scg
{


// MatrixElement
MatrixElement::~MatrixElement()
{
	VERBOSE(*verb_out << "\tMatrixElement: dying " << this << std::endl);
}

void MatrixElement::addStroke(const RawStroke* stroke) {
	strokes.push_back(stroke);
	boundingRect = (strokes.size() == 1) ? stroke->bounds() : merge(stroke->bounds(), boundingRect);

	notifyObservers(this, MatrixElementEvent::ADD_STROKE);
}

void MatrixElement::registerObserver(MatrixElementObserver& obs) const { observers.push_back(&obs); }
void MatrixElement::removeObserver(MatrixElementObserver& obs) const
{
	std::vector<MatrixElementObserver*>::iterator it;
	for (it = observers.begin(); it != observers.end() && *it != &obs; it++) {}
	if (it != observers.end()) observers.erase(it);
}

void MatrixElement::notifyObservers(const MatrixElement* changedElem, MatrixElementEvent::ACTION action)
{
	MatrixElementEvent evt;
	evt.elem = changedElem;
	evt.action = action;

	for (std::vector<MatrixElementObserver*>::iterator it = observers.begin(); it != observers.end(); it++) {
		(*it)->notify(evt);
	}
}

Rect<long> MatrixElement::getBoundingRect() const { return boundingRect; }
std::vector<const RawStroke*> MatrixElement::getStrokes() const { return strokes; }


// Matrix
Matrix::Matrix(MathRecognizer* parser) : parser(parser)
{
}

//Matrix::Matrix(const Matrix& other)
//{
//	nrows = ncols = 0;
//	numPlaceholders = other.numPlaceholders;
//	parser = other.parser;
//
//	for (unsigned i = 0; i < other.nrows; i++) {
//		for (unsigned j = 0; j < other.ncols; j++) {
//			MatrixElement* elem = other(i, j);
//			if (elem != NULL) {
//				//insert(*elem, i, j, false, false);
//				foo(*elem, i, j, false, false);
//			}
//		}
//	}
//}

Matrix::~Matrix()
{
}

void Matrix::registerObserver(MatrixObserver& obs) const
{
	observers.push_back(&obs);
}

void Matrix::removeObserver(MatrixObserver& obs) const
{
	std::vector<MatrixObserver*>::iterator it;
	for (it = observers.begin(); it != observers.end() && *it != &obs; it++) {}
	if (it != observers.end()) observers.erase(it);
}

void Matrix::notifyObservers(const MatrixElement* changedElem, MatrixEvent::ACTION action)
{
	MatrixEvent evt;
	evt.elem = changedElem;
	evt.action = action;

	for (std::vector<MatrixObserver*>::iterator it = observers.begin(); it != observers.end(); it++) {
		(*it)->notify(evt);
	}
}


}