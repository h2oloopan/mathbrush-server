#include "segment.h"
#include <fstream>
#include <vector>
#include "grouping-strat.h"
#include "rect.h"
#include "matrix.h"
#include "matrix-underspecified.h"
#include "vector.h"
#include "MatrixAnalyzer.h"
#include "MathRecognizer.h"

#include <math.h>

namespace scg
{

const double LocalOptimizationGrouping::WIDTH_TOL = 2.0;
const double LocalOptimizationGrouping::HEIGHT_TOL = 0.9;
const double LocalOptimizationGrouping::MIN_OVERLAPPING_PROPORTION = 0.25;
const double LocalOptimizationGrouping::MIN_CHAR_HEIGHT_INCHES = 0.06;
long LocalOptimizationGrouping::MIN_CHAR_HEIGHT = 0;
const double LocalOptimizationGrouping::SPAN_OVERLAP_ON_ELEM_UPDATE_TOL = 0.5;

const double LocalOptimizationGrouping::Span::BUFFER_PCT = 0.3;

LocalOptimizationGrouping::Span::Span() {}
LocalOptimizationGrouping::Span::Span(MatrixElement* elem, SpanType spanType)
{
	addElement(elem, spanType);
}

void LocalOptimizationGrouping::Span::addElement(MatrixElement* elem, SpanType spanType)
{
	elems.push_back(elem);
	spanTypes.push_back(spanType);
}

void LocalOptimizationGrouping::Span::removeElement(MatrixElement* elem)
{
	size_t i;
	for (i = 0; i < elems.size() && elems[i] != elem; i++) {}
	if (i != elems.size()) {
		elems.erase(elems.begin() + i);
		spanTypes.erase(spanTypes.begin() + i);
	}
}

bool LocalOptimizationGrouping::Span::isEmpty() const
{
	return (elems.size() == 0);
}

Interval LocalOptimizationGrouping::Span::getSpanInterval() const
{
	// We must compute the interval every time this is called, in case an element has changed (added/removed strokes)
	Interval interval; //	TODO if no elements?
	bool firstStroke = true;
	double avgHeight = 0.0;
	size_t elemsInAvgHeight = 0;

	// Go over each element in the span and merge all of its strokes' intervals
	for (size_t i = 0; i < elems.size(); i++) {
		std::vector<const RawStroke*> strokes = elems[i]->getStrokes();
		long elemHeight = elems[i]->getBoundingRect().height();
		if (elemHeight > MIN_CHAR_HEIGHT) {
			avgHeight += elemHeight;
			elemsInAvgHeight++;
		}

		// TODO make better? (observer pattern on MatrixElement to observe changes?)
		
		for (std::vector<const RawStroke*>::iterator strokeIter = strokes.begin(); strokeIter != strokes.end(); strokeIter++) {
			Rect<long> rect = (*strokeIter)->bounds();
			
			// Construct the stroke interval
			Interval strokeInterval;
			if (spanTypes[i] == HEIGHT) {
				strokeInterval.start = rect.top;
				strokeInterval.end = rect.bottom;
			} else {
				strokeInterval.start = rect.left;
				strokeInterval.end = rect.right;
			}

			interval = (firstStroke) ? strokeInterval : mergeIntervals(interval, strokeInterval);
			firstStroke = false;
		}
	}

	if (elemsInAvgHeight > 0 && interval.length() < avgHeight) {
		avgHeight /= elemsInAvgHeight;
		interval.start -= static_cast<int>(avgHeight * BUFFER_PCT);
		interval.end += static_cast<int>(avgHeight * BUFFER_PCT);
	}

	return interval;
}



LocalOptimizationGrouping::LocalOptimizationGrouping(MathRecognizer* parser, std::vector<Matrix*>& matrices, int caps)
	: ElemGroupingStrat(parser, matrices, caps)
{
	MIN_CHAR_HEIGHT = (long) (MIN_CHAR_HEIGHT_INCHES * TABLETPC_DPI);
	averageStrokeHeight = 1;
	nStrokesInAvg = 0;
}

LocalOptimizationGrouping::~LocalOptimizationGrouping()
{
}

// Probability that 2 elements with 'x' average character widths between each other should be grouped
double P(double x)
{
	/* The probability that 2 elements should be grouped is modeled as a logistic
	 * function that has approximately 100% probability if they are touching, and 0%
	 * probability if they are very far. It was chosen because it is a good model of how
	 * the grouping should work: elements very to relatively close should more than likely
	 * be grouped; after a certain factor (here it is 1 character width away), the rate at which
	 * they should be grouped decreases exponentially.
	 * 
	 * The logistic function is:
	 *  P(t) = A + (K - A) / [1 + T*e^(-B*(t + M))]^(1/T)
	 * where A-lower asymptote, K-upper asymptote, T-near which asymptote maximum growth occurs,
	 * B-growth rate, M-time of max rate of growth
	 * 
	 * Here we use A=0, K=1, T=1
	 */
	const int B = 5; //growth rate
	const int M = -1; //time of maximum rate of growth (i.e. P(M) is the inflection point)
	return 1 / (1 + exp(-B * (x + M)));
}

int LocalOptimizationGrouping::getConfidence(const Matrix* m, double& confidence) const
{
	size_t nrows, ncols;

	const UnderspecifiedMatrix* underspec = dynamic_cast<const UnderspecifiedMatrix*>(m);
	if (underspec == NULL) {
		m->getDimensions(nrows, ncols);
	} else {
		underspec->getShortFormDimensions(nrows, ncols);
	}

	if (nrows == 1 && ncols == 1) {
		confidence = 1.0;
		return 0;
	}

	// Calculate the confidence of the spacing (conf_sp)
	double conf_sp = 0.0;
	unsigned comparisons = 0;
	/* Compare how well each element ligns up with the others in its row/column. If they
	 * fully encompass each other, it is a great fit. If one is overlapping, it is not as good
	 */
	for (size_t row1 = 0; row1 < nrows; row1++) {
		for (size_t col1 = 0; col1 < ncols; col1++) {
			MatrixElement* e1 = (*m)(row1, col1);
			
			if (e1 == NULL) continue; //placeholder
			Rect<long> e1Bbox = e1->getBoundingRect();
			Interval e1Height(e1Bbox.top, e1Bbox.bottom);
			Interval e1Width(e1Bbox.left, e1Bbox.right);

			// Compare height interval to rest of elements in row
			for (size_t col2 = col1 + 1; col2 < ncols; col2++) {
				MatrixElement* e2 = (*m)(row1, col2);
				if (e2 == NULL) continue; //placeholder
				Rect<long> e2Bbox = e2->getBoundingRect();
				Interval e2Height(e2Bbox.top, e2Bbox.bottom);
				conf_sp += intervalOverlapProportion(e1Height, e2Height);
				comparisons++;
			}

			// Compare width interval to rest of elements in col
			for (size_t row2 = row1 + 1; row2 < nrows; row2++) {
				MatrixElement* e2 = (*m)(row2, col1);
				if (e2 == NULL) continue; //placeholder
				Rect<long> e2Bbox = e2->getBoundingRect();
				Interval e2Width(e2Bbox.left, e2Bbox.right);
				conf_sp += intervalOverlapProportion(e1Width, e2Width);
				comparisons++;
			}
		}
	}
	conf_sp /= comparisons;

	/* Calculate the confidence of the element grouping (conf_gr). We do this by basing it off of
	 * the approximation of the probability that 2 adjacent elements should not be grouped
	 */
	double conf_gr = 0.0;
	comparisons = 0;
	for (size_t row = 0; row < nrows; row++) {
		MatrixElement* e1 = NULL;
		size_t col1 = 0;
		// Find the first cell that is not null
		for (col1 = 0; col1 < ncols; col1++) {
			e1 = (*m)(row, col1);
			if (e1 != NULL) break;
		}

		if (e1 == NULL) continue;

		// Handle first non-null column separately
		if (row < nrows - 1) {
			MatrixElement* e3 = (*m)(row + 1, col1);
			if (e3 != NULL) {
				long dist = e3->getBoundingRect().top - e1->getBoundingRect().bottom;
				double p_grouped = P((double) dist / (averageStrokeHeight / 2.0 * HEIGHT_TOL));
				conf_gr += p_grouped;
				comparisons++;
			}
		}

		for (size_t col = col1 + 1; col < ncols; col++) {
			MatrixElement* e2 = (*m)(row, col);
			if (e2 == NULL) { //placeholder
				continue;
			}

			long dist = e2->getBoundingRect().left - e1->getBoundingRect().right;
			// Calc the probability that the strokes in the element should be grouped
			double p_grouped = P((double) dist / (averageStrokeHeight / 2.0 * WIDTH_TOL));
			conf_gr += p_grouped;
			comparisons++;

			if (row < nrows - 1) {
				MatrixElement* e3 = (*m)(row + 1, col);
				if (e3 != NULL) {
					long dist = e3->getBoundingRect().top - e1->getBoundingRect().bottom;
					double p_grouped = P((double) dist / (averageStrokeHeight / 2.0 * HEIGHT_TOL));
					conf_gr += p_grouped;
					comparisons++;
				}
			}

			e1 = e2;
		}
	}
	conf_gr /= comparisons;

	confidence = (conf_sp + conf_gr) / 2.0;
	return 0;
}

void LocalOptimizationGrouping::calcScopeRect(Rect<long>& boundingRect, long offset)
{
	boundingRect.left -= static_cast<long>(offset * WIDTH_TOL);
	boundingRect.right += static_cast<long>(offset * WIDTH_TOL);
	boundingRect.top -= static_cast<long>(offset * HEIGHT_TOL);
	boundingRect.bottom += static_cast<long>(offset * HEIGHT_TOL);
}

void LocalOptimizationGrouping::removeElementFromSpan(std::vector<Span>& spanVect, MatrixElement* elem, size_t index)
{
	if (spanVect.size() <= index) return;

	spanVect[index].removeElement(elem);
	if (spanVect[index].isEmpty()) {
		spanVect.erase(spanVect.begin() + index);
	}
}

bool LocalOptimizationGrouping::determineRowCol(MatrixElement& elem, size_t& row, size_t& col, bool& newIntermRow, bool& newIntermCol)
{
	bool ret = false;
	newIntermRow = false;
	newIntermCol = false;

	const Rect<long>& rect = elem.getBoundingRect();
	Interval heightInterval(rect.top, rect.bottom);
	Interval widthInterval(rect.left, rect.right);

	// Find row
	if (capabilities & MatrixAnalyzer::RECOGNIZE_ROWS) {
		for (size_t i = 0; i < rowSpan.size(); i++) {
			Interval rowInterval = rowSpan[i].getSpanInterval();
			double overlapProportion = intervalOverlapProportion(rowInterval, heightInterval);
			if (overlapProportion > 0) {
				if (overlapProportion > MIN_OVERLAPPING_PROPORTION) {
					// Belongs to this row
					row = i;
					break;
				} else if (heightInterval.start < rowInterval.start) {
					// Need to insert before this row
					row = i;
					newIntermRow = true;
					break;
				}
			} else if (heightInterval.end < rowInterval.start) {
				// Needs to insert before this row
				row = i;
				newIntermRow = true;
				break;
			}
			
			row = i + 1;
		}
	}
	else {
		row = 0;
	}

	VERBOSE(*verb_out << "\tLocalOptimizationGrouping: elem <" << &elem << "> has row, col: " << row);
	if (newIntermRow || row == rowSpan.size()) {
		VERBOSE(*verb_out << "(new)");

		Span rSpan(&elem, Span::HEIGHT);
		sorted_insert(rSpan, rowSpan);
		ret |= true;
	} else {
		VERBOSE(*verb_out << "(append)");

		rowSpan[row].addElement(&elem, Span::HEIGHT);
		ret |= false;
	}

	// Find column
	if (capabilities & MatrixAnalyzer::RECOGNIZE_COLUMNS) {
		for (size_t i = 0; i < colSpan.size(); i++) {
			Interval colInterval = colSpan[i].getSpanInterval();
			double overlapProportion = intervalOverlapProportion(colInterval, widthInterval);
			if (overlapProportion > 0) {
				if (overlapProportion > MIN_OVERLAPPING_PROPORTION) {
					// Belongs to this column
					col = i;
					break;
				} else if (widthInterval.start < colInterval.start) {
					// Need to insert before this column
					col = i;
					newIntermCol = true;
					break;
				}
			} else if (widthInterval.end < colInterval.start) {
				// Needs to insert before this column
				col = i;
				newIntermCol = true;
				break;
			}
			
			col = i + 1;
		}
	}
	else {
		col = 0;
	}

	VERBOSE(*verb_out << ", " << col);
	if (newIntermCol || col == colSpan.size()) {
		VERBOSE(*verb_out << "(new)" << std::endl);

		Span cSpan(&elem, Span::WIDTH);
		sorted_insert(cSpan, colSpan);
		ret |= true;
	} else {
		VERBOSE(*verb_out << "(append)" << std::endl);

		colSpan[col].addElement(&elem, Span::WIDTH);
		ret |= false;
	}

	return ret;
}

void LocalOptimizationGrouping::parseRowsCols()
{
	Matrix* m = matrices->front(); // TODO analyse more than 1 matrix?

	while (elementsToProcess.size() > 0) {
		ElementAction elemAction = elementsToProcess.front();
		elementsToProcess.pop();

		ElementAction::ActionType actionType = elemAction.actionType;
		MatrixElement* elem = elemAction.elem;
		
		if (actionType == ElementAction::UPDATE) {
			VERBOSE(*verb_out << "\tLocalOptimizationStrategy: update elem " << elem << std::endl);

			// If elements changed, check that their respective rows and columns are still valid
			size_t elemRow, elemCol;
			if (m->getPosition(*elem, elemRow, elemCol) == E_NOTFOUND) {
				continue;
			}

			Rect<long> rect = elem->getBoundingRect();
			Interval heightInterval(rect.top, rect.bottom), widthInterval(rect.left, rect.right);
			size_t newRow = elemRow, newCol = elemCol;
			if (elemRow + 1 < rowSpan.size() && 
				intervalOverlapProportion(heightInterval, rowSpan[elemRow + 1].getSpanInterval()) > SPAN_OVERLAP_ON_ELEM_UPDATE_TOL) {
				newRow++;
			} else if (elemRow > 0 && 
				intervalOverlapProportion(heightInterval, rowSpan[elemRow - 1].getSpanInterval()) > SPAN_OVERLAP_ON_ELEM_UPDATE_TOL) {
				newRow--;
			}

			if (elemCol + 1 < colSpan.size() && 
				intervalOverlapProportion(widthInterval, colSpan[elemCol + 1].getSpanInterval()) > SPAN_OVERLAP_ON_ELEM_UPDATE_TOL) {
				newCol++;
			} else if (elemCol > 0 && 
				intervalOverlapProportion(widthInterval, colSpan[elemCol - 1].getSpanInterval()) > SPAN_OVERLAP_ON_ELEM_UPDATE_TOL) {
				newCol--;
			}

			// If it belongs to another span now, move it in the matrix and the spans
			if (elemRow != newRow || elemCol != newCol) {
				m->move(*elem, newRow, newCol);

				if (elemRow != newRow) {
					rowSpan[newRow].addElement(elem, Span::HEIGHT);
					removeElementFromSpan(rowSpan, elem, elemRow);
				}
				if (elemCol != newCol) {
					colSpan[newCol].addElement(elem, Span::WIDTH);
					removeElementFromSpan(colSpan, elem, elemCol);
				}
			}
		} else if (actionType == ElementAction::ADD) {
			VERBOSE(*verb_out << "\tLocalOptimizationStrategy: add elem " << elem << std::endl);

			// If there are any new elements, add them to the matrix
			const Rect<long>& rect = elem->getBoundingRect();

			Interval heightInterval(rect.top, rect.bottom);
			Interval widthInterval(rect.left, rect.right);

			// If this is the first element, default to it being the top-left
			if (m->isEmpty()) {
				m->insert(*elem, 0, 0);
				
				Span rSpan(elem, Span::HEIGHT), cSpan(elem, Span::WIDTH);
				rowSpan.push_back(rSpan);
				colSpan.push_back(cSpan);
				continue;
			}

			size_t row = 0, col = 0;
			bool newIntermRow = false, newIntermCol = false;
			bool needNewCell = determineRowCol(*elem, row, col, newIntermRow, newIntermCol);
			
			// If we are not adding a new row or column for this element (or filling a placeholder), it should be merged with another
			MatrixElement* existingElem = NULL;
			if (!needNewCell && (existingElem = (*m)(row, col)) != NULL) {
				mergeMatrixElements(existingElem, elem);
				continue;
			}
			
			m->insert(*elem, row, col, newIntermRow, newIntermCol);
		} else if (actionType == ElementAction::REMOVE) {
			VERBOSE(*verb_out << "\tLocalOptimizationStrategy: remove elem " << elem << std::endl);

			// If any elements were removed, remove them from the matrix and the spans
			size_t row, col;
			if (m->getPosition(*elem, row, col) == E_NOTFOUND) {
				continue;
			}

			// Remove element from span calculations
			removeElementFromSpan(rowSpan, elem, row);
			removeElementFromSpan(colSpan, elem, col);
			
			// Remove element from the matrix
			m->remove(*elem);
			delete elem;
		}
	}
}

bool LocalOptimizationGrouping::mergeMatrixElements(MatrixElement* existingElem, MatrixElement* obsoleteElem)
{
	// Don't do anything if they belong to the same group
	if (existingElem == NULL || obsoleteElem == NULL || existingElem == obsoleteElem) return false;

	elementsToProcess.push(ElementAction(ElementAction::REMOVE, obsoleteElem)); // mark the obsolete element to be removed
	
	// Add all the obsolete element's strokes to the new element
	std::vector<const RawStroke*> obsoleteElemStrokes = obsoleteElem->getStrokes();
	for (std::vector<const RawStroke*>::iterator otherIter = obsoleteElemStrokes.begin(); otherIter != obsoleteElemStrokes.end(); otherIter++) {
		existingElem->addStroke(*otherIter);
		strokeToElem[*otherIter] = existingElem;
	}

	return true;
}

void LocalOptimizationGrouping::determineElement(const RawStroke* stroke)
{
	const unsigned OFFSET = averageStrokeHeight / 2;

	Rect<long> boundingRect[2];
	Rect<long> scopeRect[2];
	boundingRect[0] = stroke->bounds();

	// Inflate the bounding rect of character c to get the scope of c
	scopeRect[0] = boundingRect[0];
	calcScopeRect(scopeRect[0], OFFSET);

	// Iterate over every registered stroke and see which strokes group with the new one
	ElementAction possibleUpdate(ElementAction::UPDATE, NULL);
	for (std::vector<const RawStroke*>::iterator it = registeredStrokes.begin(); it != registeredStrokes.end(); it++) {
		const RawStroke* c = *it;

		// Calculate scope rect for this registered stroke
		boundingRect[1] = strokeToElem[c]->getBoundingRect();
		scopeRect[1] = boundingRect[1];
		calcScopeRect(scopeRect[1], OFFSET);

		// If SC(stroke) intersects any other element, or SC(c) intersects 'stroke' where c is another element, group them
		if (rectangles_overlap(scopeRect[1], boundingRect[0]) || rectangles_overlap(scopeRect[0], boundingRect[1])) {
			if (strokeToElem.find(stroke) != strokeToElem.end()) {
				// If we have already found a group for the new stroke, merge the 2 groups together
				MatrixElement* existingElem = strokeToElem[stroke];
				MatrixElement* obsoleteElem = strokeToElem[c];

				VERBOSE(*verb_out << "\tLocalOptimizationStrategy: merge " << existingElem << " with obsolete element " << obsoleteElem << " for stroke " << stroke << std::endl);

				if (mergeMatrixElements(existingElem, obsoleteElem)) {
					possibleUpdate.elem = existingElem;
				}
			} else {
				// Otherwise, create a new grouping for the element
				MatrixElement* existingElem = strokeToElem[c];
				existingElem->addStroke(stroke);
				strokeToElem[stroke] = existingElem;

				possibleUpdate.elem = existingElem;

				VERBOSE(*verb_out << "\tLocalOptimizationStrategy: append stroke " << stroke << " to group " << existingElem << std::endl);
			}
		}
	}

	if (possibleUpdate.elem != NULL) {
		elementsToProcess.push(possibleUpdate);
	}

	// If the new stroke is not part of a group yet, create a new element
	if (strokeToElem.find(stroke) == strokeToElem.end()) {
		MatrixElement* elem = new MatrixElement;
		elem->addStroke(stroke);
		elementsToProcess.push(ElementAction(ElementAction::ADD, elem));

		strokeToElem[stroke] = elem;

		VERBOSE(*verb_out << "\tLocalOptimizationStrategy: create new element " << elem << " for stroke " << stroke << std::endl);
	}
}

void LocalOptimizationGrouping::segment(std::vector<const RawStroke*>& newStrokes)
{
	for (std::vector<const RawStroke*>::iterator it = newStrokes.begin(); it != newStrokes.end(); it++) {
		segment(*it);
	}
}

void LocalOptimizationGrouping::segment(const RawStroke* stroke, bool isBeingAdded)
{
	if (isBeingAdded) {
		determineElement(stroke); // determine and store the element that the stroke belongs to
		parseRowsCols(); // determine rows and columns
		registeredStrokes.push_back(stroke);
	} else {
		/* If a stroke was removed, 
		 *  1. remove the element's strokes from the registered strokes
		 *  2. remove the element from the matrix
		 *  3. add all of the element's strokes less the removed stroke back in
		 */
		MatrixElement* elem = strokeToElem[stroke];

		if (elem == NULL) return;

		std::vector<const RawStroke*> elemStrokes = elem->getStrokes();
		std::vector<const RawStroke*> strokesLessRemovedStroke(elemStrokes);

		// Remove all of the element's strokes from the registered strokes
		for (size_t i = 0; i < registeredStrokes.size() && elemStrokes.size() > 0; i++) {
			const RawStroke* regStroke = registeredStrokes[i];
			// If this is one of the element's strokes, remove it
			for (size_t j = 0; j < elemStrokes.size(); j++) {
				if (regStroke == elemStrokes[j]) {
					registeredStrokes.erase(registeredStrokes.begin() + i);
					elemStrokes.erase(elemStrokes.begin() + j);
					i--;
					break;
				}
			}
		}

		// Re-evaluate the strokes of the element (they could now be separated)
		elementsToProcess.push(ElementAction(ElementAction::REMOVE, elem));
		
		// Remove strokeToElem association for all strokes in this element (the valid ones will be added back later)
		elemStrokes = elem->getStrokes();
		for (size_t i = 0; i < elemStrokes.size(); i++) {
			std::map<const RawStroke*, MatrixElement*>::iterator strokeElem = strokeToElem.find(elemStrokes[i]);
			if (strokeElem != strokeToElem.end()) {
				VERBOSE(*verb_out << "\tLocalOptimizationStrategy: remove <stroke, element> association <" << elemStrokes[i] << ", " << strokeToElem[elemStrokes[i]] << ">" << std::endl);
				strokeToElem.erase(strokeElem);
			}
		}

		parseRowsCols();

		// If there are other strokes in the element, add them back in
		if (strokesLessRemovedStroke.size() > 1) {
			// Remove the deleted stroke from the list, and re-evaluate all of the other strokes in that element
			for (std::vector<const RawStroke*>::iterator it = strokesLessRemovedStroke.begin(); it != strokesLessRemovedStroke.end(); it++) {
				if (*it == stroke) {
					strokesLessRemovedStroke.erase(it);
					break;
				}
			}

			for (std::vector<const RawStroke*>::iterator it = strokesLessRemovedStroke.begin(); it != strokesLessRemovedStroke.end(); it++) {
				VERBOSE(*verb_out << "\tLocalOptimizationStrategy: add stroke " << *it << " back in " << std::endl);
				segment(*it); //add it back in
			}
		}
	}
}

int LocalOptimizationGrouping::addStroke(const RawStroke* stroke)
{
	// Update average stroke height
	Rect<long> rect = stroke->bounds();

	if (rect.height() > MIN_CHAR_HEIGHT) {
		averageStrokeHeight = ((averageStrokeHeight * nStrokesInAvg) + rect.height()) / (nStrokesInAvg + 1);
		nStrokesInAvg++;
	}

	segment(stroke);

	return 0;
}

int LocalOptimizationGrouping::addStrokes(const RawStrokeGroup& strokeGroup)
{
	std::vector<const RawStroke*> newStrokes;
	for (size_t i = 0; i < strokeGroup.nstrokes; ++i) {
		const RawStroke* stroke = &strokeGroup.strokes[i];

		// Update average stroke height
		Rect<long> rect = stroke->bounds();
		if (rect.height() > MIN_CHAR_HEIGHT) {
			averageStrokeHeight = ((averageStrokeHeight * nStrokesInAvg) + rect.height()) / (nStrokesInAvg + 1);
			nStrokesInAvg++;
		}

		newStrokes.push_back(stroke);
	}

	segment(newStrokes);

	return 0;
}

int LocalOptimizationGrouping::removeStroke(const RawStroke* stroke)
{
	// Update average stroke height
	Rect<long> rect = stroke->bounds();
	if (rect.height() > MIN_CHAR_HEIGHT) {
		averageStrokeHeight = ((averageStrokeHeight * registeredStrokes.size()) - rect.height()) / (registeredStrokes.size() - 1);
		nStrokesInAvg--;
	}

	segment(stroke, false);
	return 0;
}

void LocalOptimizationGrouping::forceGrouping(const RawStroke** strokes, size_t nstrokes)
{
	// Remove the strokes from their current groupings, then create and process a new element containing these strokes
	MatrixElement* elem = new MatrixElement;
	for (size_t i = 0; i < nstrokes; i++) {
		const RawStroke* stroke = strokes[i];
		removeStroke(stroke);
		
		elem->addStroke(stroke);
		strokeToElem[stroke] = elem;
		registeredStrokes.push_back(stroke);
	}

	elementsToProcess.push(ElementAction(ElementAction::ADD, elem));
	parseRowsCols();
}

// =================================================================================================== //	


/*
int ALPHA = 4;

GravitationalGrouping::GravitationalGrouping(MathRecognizer* parser, std::vector<Matrix*>& matrices, int caps) : ElemGroupingStrat(parser, matrices, caps)
{
}

GravitationalGrouping::~GravitationalGrouping()
{
}

void GravitationalGrouping::GravitationalModel::calcRowForces()
{
	// Calculate forces between elements. Group how? normalization?
	unsigned normalizer = 0;
	size_t i=0,j=0;
	for (std::vector<GravElement>::iterator it1 = elems.begin(); it1 != elems.end(); it1++) {
		j=i+1;
		for (std::vector<GravElement>::iterator it2 = it1 + 1; it2 != elems.end(); it2++) {
			Point com1 = it1->centreOfMass, com2 = it2->centreOfMass;
			Rect<long> box1 = it1->box, box2 = it2->box;
			double dist = 
				sqrt(pow((double) std::max(box1.left - box2.right, box2.left - box1.right), 2) + pow((double) com1.y - com2.y, 2));
				//(double) std::max(box1.left - box2.right, box2.left - box1.right);

			// TODO assumes no overlapping
			double expDist = pow(dist, ALPHA);
			double force = it1->mass * it2->mass / expDist;
			/*if (normalizer == 0) {
				while ((int) force % 10 == 0) {
					force *= 10;
					normalizer += 1;
				} 
			} else {
				force *= pow((double) 10, (double) normalizer);
			}* /
			VERBOSE(*verb_out << i << " to " << j << ":\t" << force << std::endl);
			j++;
		}
		i++;
	}
}

void GravitationalGrouping::GravitationalModel::addGravitationalElement(GravElement ge)
{
	elems.push_back(ge);
	calcRowForces();
}

void GravitationalGrouping::addGravitationalMass(size_t numPoints, Rect<long> box)
{
	GravElement ge;
	ge.box = box;
	ge.mass = 
		//(double) numPoints;// / (max(box.right - box.left, 1) * max(box.bottom - box.top, 1));
		(double) box.area();
	
	Point p;
	p.x = (ge.box.right + ge.box.left) / 2;
	p.y = (ge.box.bottom + ge.box.top) / 2;
	ge.centreOfMass = p;

	gravModel.addGravitationalElement(ge);
}

int GravitationalGrouping::addStroke(const RawStroke* stroke)
{
	Matrix* m = matrices->front();

	Rect<long> rect = stroke->bounds();
	addGravitationalMass(stroke->npoints, rect);

	return 0; 
}

int GravitationalGrouping::addStrokes(const RawStrokeGroup& strokeGroup)
{
	for (size_t i = 0; i < strokeGroup.nstrokes; ++i) {
		const RawStroke* stroke = &strokeGroup.strokes[i];
		addStroke(stroke);
	}

	return 0;
}

int GravitationalGrouping::removeStroke(const RawStroke* stroke)
{
	return 0;
}

int GravitationalGrouping::getConfidence(const Matrix* m, double& confidence) const
{
	confidence = 1.0;
	return 0;
}

*/

}