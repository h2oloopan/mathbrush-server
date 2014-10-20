#include "MatrixAnalyzer.h"
#include "grouping-strat.h"
#include "rect.h"
#include "verb.h"
#include "stkutils.h"
#include <cstring>

namespace scg
{

MatrixAnalyzer::MatrixAnalyzer(MathRecognizer* parser, int caps) : parser(parser), matrixMode(FULLY_SPECIFIED), capabilities(caps)
{
	lastThreeStrokes = new const RawStroke*[3];
	lastThreeStrokesIndex = 0;

	Matrix* m = new FullyspecifiedMatrix(parser);
	fullySpecifiedMatrices.push_back(m);

	MatrixAnalyzer::groupingStrat = 
		new LocalOptimizationGrouping(parser, fullySpecifiedMatrices, caps);
		//new InAirGrouping(parser, fullySpecifiedMatrices);
		//new GravitationalGrouping(parser, fullySpecifiedMatrices);
}

MatrixAnalyzer::~MatrixAnalyzer()
{
	while (underspecifiedMatrices.size() > 0) {
		Matrix* m = underspecifiedMatrices.back();
		underspecifiedMatrices.pop_back();

		size_t rows, cols;
		const UnderspecifiedMatrix* underspec = dynamic_cast<const UnderspecifiedMatrix*>(m);
		if (underspec == NULL) {
			m->getDimensions(rows, cols);
		} else {
			underspec->getShortFormDimensions(rows, cols);
		}

		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				try {
					MatrixElement* elem = (*m)(i,j);
					if (elem != NULL) {
						delete elem;
					}
				} catch (Matrix::InvalidBoundsException) {
				}
			}
		}

		delete m;
	}

	delete groupingStrat;
}

void MatrixAnalyzer::trackLastTwoStrokes()
{
	lastThreeStrokes[0] = lastThreeStrokes[1];
	lastThreeStrokes[1] = lastThreeStrokes[2];
	lastThreeStrokesIndex = 2;
}

bool MatrixAnalyzer::handleEllipsis(const RawStroke* stroke) {
	if (!(capabilities & RECOGNIZE_ELLIPSIS)) { return false; }

	if (!stroke_is_dot(*stroke)) {
		lastThreeStrokesIndex = 0;
		return false;
	}

	lastThreeStrokes[lastThreeStrokesIndex++] = stroke;

	if (lastThreeStrokesIndex != 3) {
		return false;
	}

	ExpressionIterator* eit = parser->CreateIteratorForStrokes(lastThreeStrokes, 3, false);
	if (!eit) {
		trackLastTwoStrokes();
		return false;
	}
		
	const ExpressionTree* tree = eit->next();
	if (!tree || tree->type() != DOTS_EXPR) {
		trackLastTwoStrokes();
		return false;
	}

	const char* name = tree->child(0)->str();
	VERBOSE(*verb_out << "\tMatrixAnalyzer: detected ellipsis" << std::endl);

	groupingStrat->forceGrouping(lastThreeStrokes, 3);

	UnderspecifiedMatrix::EllipsisType ellipsisType;
	if (std::strcmp(name, "&hellip;") == 0) {
		ellipsisType = UnderspecifiedMatrix::HORIZONTAL;
	}
	else if (std::strcmp(name, "&vellip;") == 0) {
		ellipsisType = UnderspecifiedMatrix::VERTICAL;
	}
	else if (std::strcmp(name, "&dtdots;") == 0) {
		ellipsisType = UnderspecifiedMatrix::DIAGONAL;
	}
	else if (std::strcmp(name, "&utdots;") == 0) {
		ellipsisType = UnderspecifiedMatrix::ANTIDIAGONAL;
	}
	else {
		THROW_ERROR(E_INVALID, "unknown DOTS_EXPR type " << name);
	}

	// Set underspecified matrix mode
	if (matrixMode != UNDER_SPECIFIED) {
		VERBOSE(*verb_out << "\tMatrixAnalyzer: convert to underspecified matrix" << std::endl);
		matrixMode = UNDER_SPECIFIED;

		// Convert each active matrix to an underspecified matrix
		for (size_t i = 0; i < fullySpecifiedMatrices.size(); i++) {
			Matrix* m = fullySpecifiedMatrices[i];
			UnderspecifiedMatrix* underspecifiedMatrix = new UnderspecifiedMatrix(parser, m);
			underspecifiedMatrix->addEllipsis(lastThreeStrokes, ellipsisType);

			underspecifiedMatrices.push_back(underspecifiedMatrix);
		}
	}
	else {
		// Convert each active matrix to an underspecified matrix
		for (size_t i = 0; i < underspecifiedMatrices.size(); i++) {
			UnderspecifiedMatrix* m = dynamic_cast<UnderspecifiedMatrix*>(underspecifiedMatrices[i]);
			if (m != NULL) {
				m->addEllipsis(lastThreeStrokes, ellipsisType);
			}
		}
	}

	lastThreeStrokesIndex = 0;

	groupingStrat->setMatrices(underspecifiedMatrices);
	return true;
}

int MatrixAnalyzer::addStroke(const RawStroke* stroke, bool update)
{
	// If we're not signalled to update, store stroke in buffer
	if (!update) {
		BufferEntry entry;
		entry.Command = BufferEntry::AddStroke;
		entry.stroke = stroke;
		buffer.push(entry);
		return 0;
	}

	VERBOSE(*verb_out << "\tMatrixAnalyzer: add stroke " << stroke << std::endl;);

	if (handleEllipsis(stroke)) return 0;

	return groupingStrat->addStroke(stroke);
}

int MatrixAnalyzer::addStrokes(const RawStrokeGroup& strokeGroup, bool update)
{
	// If we're not signalled to update, store stroke group in buffer
	if (!update) {
		BufferEntry entry;
		entry.Command = BufferEntry::AddStrokeGroup;
		entry.strokeGroup = &strokeGroup;
		buffer.push(entry);
		return 0;
	}

	VERBOSE(*verb_out << "\tMatrixAnalyzer: adding " << strokeGroup.nstrokes << " strokes" << std::endl;);

	// Detect an ellipsis
	for (size_t i = 0; i < strokeGroup.nstrokes; ++i) {
		const RawStroke* stroke = &strokeGroup.strokes[i];
		addStroke(stroke, update);
	}

	return 0;
}

int MatrixAnalyzer::removeStroke(const RawStroke* stroke, bool update)
{
	// If we're not signalled to update, store stroke in buffer
	if (!update) {
		BufferEntry entry;
		entry.Command = BufferEntry::RemoveStroke;
		entry.stroke = stroke;
		buffer.push(entry);
		return 0;
	}

	VERBOSE(*verb_out << "\tMatrixAnalyzer: removing stroke " << stroke << std::endl;);

	return groupingStrat->removeStroke(stroke);
}

int MatrixAnalyzer::updateRecognition()
{
	// Go through the buffer and process the commands
	while (!buffer.empty()) {
		BufferEntry entry = buffer.front();
		
		switch (entry.Command) {
			case BufferEntry::AddStroke:
				addStroke(entry.stroke);
				break;
			case BufferEntry::AddStrokeGroup:
				addStrokes(*entry.strokeGroup);
				break;
			case BufferEntry::RemoveStroke:
				removeStroke(entry.stroke);
				break;
		};

		buffer.pop();
	}
	return 0;
}

int MatrixAnalyzer::getConfidence(const Matrix* m, double& confidence) const
{
	int status = groupingStrat->getConfidence(m, confidence);
	VERBOSE(*verb_out << "\tMatrixAnalyzer: confidence: " << confidence << std::endl;);
	return status;
}

bool MatrixAnalyzer::hasValidLongForm()
{
	if (matrixMode == UNDER_SPECIFIED) {
		for (std::vector<Matrix*>::iterator it = underspecifiedMatrices.begin(); it != underspecifiedMatrices.end(); it++) {
			UnderspecifiedMatrix* m = dynamic_cast<UnderspecifiedMatrix*>(*it);
			m->rebuild();
			if (m->isValid()) {
				return true;
			}
		}
	}

	return false;
}

MatrixIterator* MatrixAnalyzer::createIterator()
{
	if (hasValidLongForm()) {
		return new MatrixIterator(underspecifiedMatrices);
	}
	
	return new MatrixIterator(fullySpecifiedMatrices);
}

LongFormMatrixIterator* MatrixAnalyzer::createLongFormIterator()
{
	LongFormMatrixIterator* lfmi = new LongFormMatrixIterator(underspecifiedMatrices);
	return lfmi;
}

MatrixIterator::MatrixIterator(std::vector<Matrix*> matrices)
{
	index = 0;
	for (std::vector<Matrix*>::const_iterator it = matrices.begin(); it != matrices.end(); it++) {
		MatrixIterator::matrices.push_back(*it);
	}
}

MatrixIterator::~MatrixIterator()
{
	/*while (matrices.size() > 0) {
		Matrix* m = matrices.back();
		matrices.pop_back();
		delete m;
	}*/
}

Matrix* MatrixIterator::next()
{
	if (index == matrices.size()) {
		return NULL;
	}
	Matrix* m = matrices.at(index++);
	m->setRetrievalType(true);
	return m;
}

bool MatrixIterator::hasNext()
{
	return index < matrices.size();
}

LongFormMatrixIterator::LongFormMatrixIterator(std::vector<Matrix*> matrices)
{
	index = 0;
	LongFormMatrixIterator::matrices = matrices;
	/*UnderspecifiedMatrix* um = NULL;
	for (std::vector<Matrix*>::const_iterator it = matrices.begin(); it != matrices.end(); it++) {
		Matrix* m = *it;
		if ((um = dynamic_cast<UnderspecifiedMatrix*>(m)) != 0) {
			UnderspecifiedMatrix* newUM = new UnderspecifiedMatrix(*um);
			LongFormMatrixIterator::matrices.push_back(newUM);
		}
	}*/
}

LongFormMatrixIterator::~LongFormMatrixIterator()
{
}

Matrix* LongFormMatrixIterator::next()
{
	if (index == matrices.size()) {
		return NULL;
	}
	Matrix* m = matrices.at(index++);
	m->setRetrievalType(false);
	return m;
}

bool LongFormMatrixIterator::hasNext()
{
	return index < matrices.size();
}


MatrixAnalyzer* CreateMatrixAnalyzer(MathRecognizer* parser, int caps)
{
	return new MatrixAnalyzer(parser, caps);
}

void DestroyMatrixAnalyzer(MatrixAnalyzer* ma)
{
	if (ma != NULL) {
		delete ma;
	}
}

};
