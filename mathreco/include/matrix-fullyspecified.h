#ifndef SCG_MATRIX_FULLYSPECIFIED
#define SCG_MATRIX_FULLYSPECIFIED

#include "matrix.h"


namespace scg
{


class UnderspecifiedMatrix;

class FullyspecifiedMatrix : public Matrix
{
	friend class UnderspecifiedMatrix;

	size_t curr_row, curr_col;

	std::vector<std::vector<MatrixElement*> > cells; // MatrixElement cells, accessed by cells[row][col]
	size_t nrows, ncols; // number of rows and columns (each row will always have the same number of columns)
	size_t numPlaceholders; // total number of placeholder (NULL) cells

	mutable std::string matrixStr;

public:
	FullyspecifiedMatrix();
	FullyspecifiedMatrix(MathRecognizer* parser);
	virtual ~FullyspecifiedMatrix();

	void setRetrievalType(bool isShortForm);
	int getCell(size_t row, size_t col, ExpressionIterator*& results) const;
	int getPosition(const MatrixElement& elem, size_t& row, size_t& col) const;
	int getDimensions(size_t& rows, size_t& cols) const;

	void rebuild();

	void insert(MatrixElement& elem, size_t row, size_t col, bool newIntermediateRow=false, bool newIntermediateCol=false);

	void remove(MatrixElement& elem);
	void move(MatrixElement& elem, size_t newRow, size_t newCol); // if the move makes an empty row or column, it will be removed
	void clear();

	bool hasPlaceholders() const;
	bool isEmpty() const;
	const char* str() const;

	MatrixElement* operator() (const size_t row, const size_t col) const;

private:
	bool deleteRowIfEmpty(size_t row);
	bool deleteColIfEmpty(size_t col);
};


}

#endif