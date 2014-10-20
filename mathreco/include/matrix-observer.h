#ifndef SCG_MATRIX_OBSERVER_H
#define SCG_MATRIX_OBSERVER_H

#include "matrix.h"

namespace scg
{


class MatrixElement;

struct MatrixEvent
{
	enum ACTION {INSERT, MOVE, REMOVE} action;
	const MatrixElement* elem;

	MatrixEvent() : elem(NULL) {}
	MatrixEvent(MatrixEvent::ACTION action, MatrixElement* elem) : action(action), elem(elem) {}
};

// Defines the interface for an object who wants to observe changes to a Matrix instance
class MatrixObserver {
public:
	virtual void notify(MatrixEvent evt) = 0;
};


struct MatrixElementEvent
{
	enum ACTION {ADD_STROKE} action;
	const MatrixElement* elem;
};

// Defines the interface for an object who wants to observe changes to a MatrixElement
class MatrixElementObserver {
public:
	virtual void notify(MatrixElementEvent evt) = 0;
};


}

#endif