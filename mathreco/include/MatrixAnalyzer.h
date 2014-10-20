#ifndef SCG_MATRIXRECOGNIZER_H_
#define SCG_MATRIXRECOGNIZER_H_

#include <vector>
#include <queue>
#include "MathRecoTypes.h"
#include "group.h"
#include "grouping-strat.h"
#include "matrix.h"
#include "matrix-underspecified.h"

namespace scg
{

class MatrixIterator
{
protected:
	size_t index;
	std::vector<Matrix*> matrices;

	MatrixIterator() {}

public:
	MatrixIterator(std::vector<Matrix*> matrices);
	virtual ~MatrixIterator();

	virtual Matrix* next();
	virtual bool hasNext();
};

class LongFormMatrixIterator : public MatrixIterator
{

public:
	LongFormMatrixIterator(std::vector<Matrix*> matrices);
	~LongFormMatrixIterator();

	Matrix* next();
	bool hasNext();
};

class MatrixAnalyzer
{
public:
	enum {
		RECOGNIZE_ROWS = 1,
		RECOGNIZE_COLUMNS = 2,
		RECOGNIZE_ELLIPSIS = 4,
		RECOGNIZE_FULL = 7
	};

private:
	enum MatrixMode { FULLY_SPECIFIED, UNDER_SPECIFIED } matrixMode;

	MathRecognizer* parser;
	ElemGroupingStrat* groupingStrat; // Strategy Pattern of matrix element grouping strategies
	std::vector<Matrix*> fullySpecifiedMatrices, underspecifiedMatrices; // vector of all the matrices being considered
	const RawStroke** lastThreeStrokes; // To check for ellipses
	unsigned lastThreeStrokesIndex; // Index of the next holder in lastThreeStrokes

	// Buffer to store entries in order if waiting for updateRecognition() call
	struct BufferEntry {
		enum CommandTypes { AddStrokeGroup, AddStroke, RemoveStroke } Command;

		const RawStrokeGroup* strokeGroup;
		const RawStroke* stroke;
	};
	std::queue<BufferEntry> buffer;

	int capabilities;

public:
	MatrixAnalyzer(MathRecognizer* parser, int caps = RECOGNIZE_FULL);
	~MatrixAnalyzer();

	/* Add a stroke for the matrix analyzer to consider
	 * @param update whether the matrix analyzer should perform logic update or not (blocks)
	 */
	int addStroke(const RawStroke* stroke, bool update=true);

	/* Add multiple strokes for the matrix analyzer to consider
	 * @param update whether the matrix analyzer should perform logic update or not (blocks)
	 */
	int addStrokes(const RawStrokeGroup& strokeGroup, bool update=true);

	/* Remove a stroke from the matrix analyzer. This does *not* delete the stroke
	 * (that should be done in MathRecognizer)
	 */
	int removeStroke(const RawStroke* stroke, bool update=true);

	// Allow the matrix analyzer to perform a logic update on the newly added/removed strokes
	int updateRecognition();

	// Get the confidence that a matrix "m" is the correct matrix
	int getConfidence(const Matrix* m, double& confidence) const;

	/* Returns if we are dealing with underspecified matrices or not. If we are, we can make use
	 * of their long forms through the iterators
	 */
	bool hasValidLongForm();

	// Create an iterator over all of the generated matrices
	MatrixIterator* createIterator();

	// Create an iterator over all of the underspecified matrices that will generate the long form of the matrix
	LongFormMatrixIterator* createLongFormIterator();

private:
	// Check the new stroke to see if it forms an ellipsis
	bool handleEllipsis(const RawStroke* stroke);

	void trackLastTwoStrokes();
};

MatrixAnalyzer* CreateMatrixAnalyzer(MathRecognizer* parser, int capabilities = MatrixAnalyzer::RECOGNIZE_FULL);
void DestroyMatrixAnalyzer(MatrixAnalyzer* ma);

}


#endif
