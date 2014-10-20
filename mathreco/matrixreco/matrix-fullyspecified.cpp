#include "matrix-fullyspecified.h"

#include <sstream>
#include "interval.h"
#include "MathRecoTypes.h"
#include "MathRecognizer.h"
#include "error.h"
#include "grammar.h"
#include "mathrecognizer-private.h"

namespace scg
{

static const size_t badpos = ~0;

// Matrix
FullyspecifiedMatrix::FullyspecifiedMatrix() : Matrix(), nrows(0), ncols(0), numPlaceholders(0), curr_row(badpos), curr_col(badpos)
{
}

FullyspecifiedMatrix::FullyspecifiedMatrix(MathRecognizer* parser) : Matrix(parser), nrows(0), ncols(0), numPlaceholders(0), curr_row(badpos), curr_col(badpos)
{
}

FullyspecifiedMatrix::~FullyspecifiedMatrix()
{
}

int FullyspecifiedMatrix::getDimensions(size_t& rows, size_t& cols) const
{
	rows = nrows;
	cols = ncols;
	return 0;
}

MatrixElement* FullyspecifiedMatrix::operator() (const size_t row, const size_t col) const
{
	if (row >= nrows || col >= ncols) {
		throw InvalidBoundsException();
	}
	
	return cells[row][col];
}

int FullyspecifiedMatrix::getCell(size_t row, size_t col, ExpressionIterator*& results) const
{
	results = NULL;

	MatrixElement* elem = NULL;
	try {
		elem = (*this)(row, col);
	} catch (InvalidBoundsException) {
		results = NULL;
		return E_NOTFOUND;
	}

	if (elem != NULL) {
		std::vector<const RawStroke*> strokeVect = elem->getStrokes();
		size_t numStrokes = strokeVect.size();
		// If there are strokes (it is not an empty/placeholder cell)
		if (numStrokes != 0) {
			// Get the strokes from the cell and create an ExpressionIterator for them
			const RawStroke** strokes = new const RawStroke*[numStrokes];
			for (size_t i = 0; i < numStrokes; i++) {
				strokes[i] = strokeVect[i];
			}
			math_recognizer_base *rec = (math_recognizer_base *)parser;
			const nonterminal *nt;
			if (row == curr_row && col == curr_col) {
				// allow partial expressions in the matrix cell currently being added to
				nt = GetMathGrammar()->single_expr_root;
			}
			else {
				// allow only complete expressions otherwise
				nt = GetMathGrammar()->matrix_cell_root;
				// XXX: this is a hack that relies on the grammar being in a certain form
				//compexpr_root->productions[0]->rhs[0];
			}
			assert(nt);
			results = rec->CreateIteratorForStrokes(nt, strokes, numStrokes, false);
			//results = parser->CreateIteratorForStrokes(strokes, numStrokes, false);
			delete [] strokes;
			return 0;
		}
	}

	if (results == NULL) {
		results = parser->CreateIteratorForExpression(CreatePlaceholderExpression(), true, false);
	}
	return 0;
}

void FullyspecifiedMatrix::setRetrievalType(bool isShortForm)
{
	// no-op
}

int FullyspecifiedMatrix::getPosition(const MatrixElement& elem, size_t& row, size_t& col) const
{
	// O(n) for now
	for (unsigned i = 0; i < nrows; i++) {
		for (unsigned j = 0; j < ncols; j++) {
			if (cells[i][j] == &elem) {
				row = i;
				col = j;
				return 0;
			}
		}
	}

	return E_NOTFOUND;
}

void FullyspecifiedMatrix::rebuild()
{
}

void FullyspecifiedMatrix::insert(MatrixElement& elem, size_t row, size_t col, bool newIntermediateRow, bool newIntermediateCol)
{
	// Insert new row
	if ((nrows < (row + 1)) || newIntermediateRow) {
		std::vector<MatrixElement*> emptyVect;
		for (size_t i = 0; i < ncols; i++) { emptyVect.push_back(NULL); }

		if (newIntermediateRow) {
			cells.insert(cells.begin() + row, emptyVect);
		} else {
			for (size_t i = cells.size(); i <= row; i++) {
				cells.push_back(emptyVect);
			}
		}
		numPlaceholders += (newIntermediateRow) ? ncols : ncols * ((row + 1) - nrows);
		nrows += (newIntermediateRow) ? 1 : (row + 1) - nrows;
	}
	// Insert placeholder column in all rows
	if (ncols < (col + 1) || newIntermediateCol) {
		MatrixElement* emptyElem = NULL;
		for (size_t row = 0; row < cells.size(); row++) {
			if (newIntermediateCol) {
				cells[row].insert(cells[row].begin() + col, emptyElem);
			} else {
				for (size_t j = cells[row].size(); j <= col; j++) {
					cells[row].push_back(emptyElem);
				}
			}
		}
		numPlaceholders += (newIntermediateCol) ? nrows : nrows * ((col + 1) - ncols);
		ncols += (newIntermediateCol) ? 1 : (col + 1) - ncols;
	}

	if (cells[row][col] == NULL) {
		VERBOSE(*verb_out << "\tMATRIX: replace placeholder with element " << &elem << std::endl);

		// If it's a placeholder cell, fill it
		cells[row][col] = &elem;
		numPlaceholders--;
	} else {
		VERBOSE(*verb_out << "\tMATRIX: insert element " << &elem << std::endl);

		// If not a placeholder, insert elem at (row, col) and an empty element at
		// (r2, col) for r2 != row, pushing all following cells over one space in each row
		for (size_t r = 0; r < nrows; r++) {
			MatrixElement* newElem = NULL;
			if (r == row) newElem = &elem;
			cells[r].insert(cells[r].begin() + col, newElem);
		}
		numPlaceholders += nrows - 1;
		ncols++;
	}

	curr_row = row;
	curr_col = col;

	// Notify observers that we have inserted a new element into the matrix
	notifyObservers(&elem, MatrixEvent::INSERT);
}

bool FullyspecifiedMatrix::deleteRowIfEmpty(size_t row)
{
	// If this row is now empty, remove it
	size_t c;
	for (c = 0; c < ncols; c++) {
		if (cells[row][c] != NULL) {
			break;
		}
	}
	if (c == ncols) {
		cells.erase(cells.begin() + row);
		nrows--;
		numPlaceholders -= ncols;
		return true;
	}

	return false;
}

bool FullyspecifiedMatrix::deleteColIfEmpty(size_t col)
{
	// If this column is now empty, remove it
	size_t r;
	for (r = 0; r < nrows; r++) {
		if (cells[r][col] != NULL) {
			break;
		}
	}
	if (r == nrows) {
		for (r = 0; r < nrows; r++) {
			cells[r].erase(cells[r].begin() + col);
		}
		ncols--;
		numPlaceholders -= nrows;
		return true;
	}

	return false;
}

void FullyspecifiedMatrix::remove(MatrixElement& elem)
{
	VERBOSE(*verb_out << "\tMATRIX: remove element " << &elem << std::endl);

	size_t rowDel, colDel;
	if (getPosition(elem, rowDel, colDel) == E_NOTFOUND) {
		return;
	}
	cells[rowDel][colDel] = NULL;
	numPlaceholders++;

	deleteRowIfEmpty(rowDel);
	deleteColIfEmpty(colDel);

	notifyObservers(&elem, MatrixEvent::REMOVE);
}

void FullyspecifiedMatrix::move(MatrixElement& elem, size_t newRow, size_t newCol)
{
	size_t oldRow, oldCol;
	if (getPosition(elem, oldRow, oldCol) == E_NOTFOUND) {
		return;
	}

	VERBOSE(*verb_out << "\tMATRIX: move element " << &elem << " from " << oldRow << "," << oldCol << " to " << newRow << "," << newCol << std::endl);

	cells[oldRow][oldCol] = NULL;
	numPlaceholders++;
	insert(elem, newRow, newCol);

	deleteRowIfEmpty(oldRow);
	deleteColIfEmpty(oldCol);

	notifyObservers(&elem, MatrixEvent::MOVE);
}

void FullyspecifiedMatrix::clear()
{
	for (size_t i = 0; i < cells.size(); i++) {
		cells[i].clear();
	}
	cells.clear();

	nrows = 0;
	ncols = 0;
	numPlaceholders = 0;
}

bool FullyspecifiedMatrix::hasPlaceholders() const
{
	return numPlaceholders > 0;
}

bool FullyspecifiedMatrix::isEmpty() const
{
	return cells.size() == 0;
}

const char* FullyspecifiedMatrix::str() const
{
	std::stringstream ss;

	for (size_t i = 0; i < nrows; i++) {
		for (size_t j = 0; j < ncols; j++) {
			MatrixElement* elem = NULL;
			try { elem = (*this)(i,j); }
			catch (InvalidBoundsException) { continue; }

			size_t numStrokes = 0;
			if (elem != NULL) {
				numStrokes = elem->getStrokes().size();
			}

			ss << i << "," << j << " (" << numStrokes << ")";

			ExpressionIterator* results;
			getCell(i, j, results);
			if (results != NULL) {
				const ExpressionTree* tree = results->next();
				if (tree != NULL) {
					ss << ": " << tree->str();
					// XXX tree->release();
					// XXX results->release();
				}
			}

			ss << "\t";
		}
		ss << std::endl;
	}

	matrixStr = ss.str();
	const char* output = matrixStr.c_str();
	return output;
}


}