#ifndef SCG_MATRIX_H_
#define SCG_MATRIX_H_

#include <vector>
#include "rect.h"
#include "stroke.h"
#include "MathRecoTypes.h"
#include "matrix-observer.h"
	
namespace scg
{


// Defines the MatrixElement type that will make up a Matrix.
class MatrixElement
{
	Rect<long> boundingRect; //bounding rect of all the MatrixElement's strokes
	std::vector<const RawStroke*> strokes;
	mutable std::vector<MatrixElementObserver*> observers;

public:
	// TODO check for no memory leaks (counter on con/destructor)
	~MatrixElement();

	void addStroke(const RawStroke* stroke);

	Rect<long> getBoundingRect() const;
	std::vector<const RawStroke*> getStrokes() const;

	// Observer pattern methods
	void registerObserver(MatrixElementObserver& obs) const;
	void removeObserver(MatrixElementObserver& obs) const;

private:
	void notifyObservers(const MatrixElement* changedElem, MatrixElementEvent::ACTION action);
};


/* Defines a dynamic abstract matrix. At any point in time, each row must have the same number of
 * columns, where some of these cells may be empty placeholder cells (i.e. a NULL
 * MatrixElement instance).
 *
 * The Matrix destructor does not delete MatrixElement pointers - they may be being used externally.
 */
class Matrix
{
	mutable std::vector<MatrixObserver*> observers;

protected:
	MathRecognizer* parser;

public:
	/* A matrix will only deallocate the MatrixElement that are contained within the Matrix at the end of its
	 * life. If any were removed, they will not be handled
	 */
	Matrix(MathRecognizer* parser = NULL);
	//Matrix(const Matrix& other); // shallow copy (MatrixElement addresses are preserved i.e. not copied)
	virtual ~Matrix();

	// Get the dimensions that can be acquired by a call to getCell()
	virtual int getDimensions(size_t& rows, size_t& cols) const = 0;

	/* Get an ExpressionIterator of the contents of a cell. Note: "results" should be released by caller.
	 * Returns:
	 *  - E_IGNORE if there are no strokes (it is an empty/placeholder cell)
	 *  - E_NOTFOUND if no iterator for strokes was found
	 *  - E_OK if a valid result set was found
	 */
	virtual int getCell(size_t row, size_t col, ExpressionIterator*& results) const = 0;
	virtual void setRetrievalType(bool isShortForm) = 0;

	virtual void rebuild() = 0;

	/* Get the row and column of an element, returned in "row" and "col"
	 * Return E_NOTFOUND if the element was not found in the matrix
	 */
	virtual int getPosition(const MatrixElement& elem, size_t& row, size_t& col) const = 0;

	/* If there is a non-empty MatrixElement at (row, col), then this will insert
	 * elem at (row, col) and an empty (placeholder) element at (r2, col) for r2 != row, 
	 * pushing all following cells over one space in each row
	 */
	virtual void insert(MatrixElement& elem, size_t row, size_t col, bool newIntermediateRow=false, bool newIntermediateCol=false) = 0;

	// Remove an element from the matrix. Once removed, the matrix will not deallocate the element's memory
	virtual void remove(MatrixElement& elem) = 0;
	// Move an element in the matrix to a new row,col. If the move results in an empty row or column, that row/column will be removed.
	virtual void move(MatrixElement& elem, size_t newRow, size_t newCol) = 0;
	virtual void clear() = 0;

	virtual bool hasPlaceholders() const = 0;
	virtual bool isEmpty() const = 0;
	virtual const char* str() const = 0;

	// Allow access of the form matrix(row, col). If no element is stored, returns NULL
	virtual MatrixElement* operator() (const size_t row, const size_t col) const = 0;
	class InvalidBoundsException {}; // used in operator() overload

	// Observer pattern methods
	void registerObserver(MatrixObserver& obs) const;
	void removeObserver(MatrixObserver& obs) const;

protected:
	void notifyObservers(const MatrixElement* changedElem, MatrixEvent::ACTION action);
};

}

#endif
