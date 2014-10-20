#ifndef SCG_GROUPING_STRAT_H_
#define SCG_GROUPING_STRAT_H_

#include <algorithm>
#include <vector>
#include <queue>
#include <map>

#include "stroke.h"
#include "group.h"
#include "matrix.h"
#include "interval.h"
#include "vector.h"

namespace scg
{

/* Interface for the element grouping strategy. Inheriting classes will implement a certain algorithm for grouping strokes into
 * matrix elements and parsing them into rows and columns.
 */
class ElemGroupingStrat
{
protected:
	const std::vector<Matrix*>* matrices; // From MatrixAnalyzer
	const MathRecognizer* parser;
	int capabilities;

public:
	ElemGroupingStrat(MathRecognizer* parser, std::vector<Matrix*>& matrices, int caps) : parser(parser), matrices(&matrices), capabilities(caps) {}

	void setMatrices(std::vector<Matrix*>& matrices) { ElemGroupingStrat::matrices = &matrices; }

	virtual int addStroke(const RawStroke* stroke) = 0;
	virtual int addStrokes(const RawStrokeGroup& strokeGroup) = 0;
	// This does *not* delete the stroke (that should be done in MathRecognizer)
	virtual int removeStroke(const RawStroke* stroke) = 0;

	virtual void forceGrouping(const RawStroke** strokes, size_t nstrokes) = 0;

	virtual int getConfidence(const Matrix* m, double& confidence) const = 0;
};

/* Rendition of algorithm according to Tausky et. al., "Managing Ambiguity in Mathematical Matrices"
 */
class LocalOptimizationGrouping : public ElemGroupingStrat
{
	static const double WIDTH_TOL; // tolerance for how much to increase a stroke's scope's width
	static const double HEIGHT_TOL; // tolerance for how much to increase a stroke's scope's height
	static const double MIN_OVERLAPPING_PROPORTION; // minimum overlap proportion (as a %) between an element's width/height and the row/column span for it to be part of that row/column 
	static const double MIN_CHAR_HEIGHT_INCHES; // min heigh in inches (to calculate MIN_CHAR_HEIGHT)
	static long MIN_CHAR_HEIGHT; // minimum height (in pixels) for a stroke to constitute a character (omit elements such as '-')
	static const double SPAN_OVERLAP_ON_ELEM_UPDATE_TOL; // same as MIN_OVERLAPPING_PROPORTION, except applied when an element is being updated (not created). this value should be > MIN_OVERLAPPING_PROPORTION because it should have a larger boundary to overcome if it flows into another row/column

	// Defines the height or width span of a group of MatrixElement objects (i.e. the combined union of their width or height)
	class Span {
	public:
		enum SpanType { HEIGHT, WIDTH };
	private:
		/* Add a certain percent of the average height of the span's elements in ink pixels to the start and end of the span as a buffer
		 * if the span is very narrow (e.g. for the number '1')
		 */
		static const double BUFFER_PCT;

		std::vector<MatrixElement*> elems; // elements comprised in this span
		std::vector<SpanType> spanTypes; // these two vectors are concurrently accessed i.e. elem[i] has span type spanTypes[i]
	public:
		Span();
		Span(MatrixElement* elem, SpanType spanType);

		void addElement(MatrixElement* elem, SpanType spanType);
		void removeElement(MatrixElement* elem);
		bool isEmpty() const;

		// TODO assumes that the left-hand side of the interval doesn't change after the first stroke (because of sorted_insert) is this ok?
		Interval getSpanInterval() const;

		inline bool operator==(const Span& other) const { return this->getSpanInterval() == other.getSpanInterval(); }
		inline bool operator<=(const Span& other) const { return this->getSpanInterval() <= other.getSpanInterval(); }
		inline bool operator<(const Span& other) const { return this->getSpanInterval() < other.getSpanInterval(); }
	}; // Span

	// Remove an element from the list of spans. If the resulting span doesn't have any more elements, remove the entire span from the list
	void removeElementFromSpan(std::vector<Span>& spanVect, MatrixElement* elem, size_t index);

	// Used in a FIFO queue explaining the actions to take on an element
	struct ElementAction {
		enum ActionType { ADD, UPDATE, REMOVE } actionType;
		MatrixElement* elem;

		ElementAction(ActionType actionType, MatrixElement* elem) : actionType(actionType), elem(elem) {}
	};

	unsigned averageStrokeHeight, nStrokesInAvg;
	std::vector<const RawStroke*> registeredStrokes; // cache of all strokes registered in the grouping
	std::vector<Span> rowSpan, colSpan; // spans of the rows and columns
	std::queue<ElementAction> elementsToProcess; // elements that have been added/removed/changed
	std::map<const RawStroke*, MatrixElement*> strokeToElem; // maps a stroke to its element group

public:
	LocalOptimizationGrouping(MathRecognizer* parser, std::vector<Matrix*>& matrices, int caps);
	virtual ~LocalOptimizationGrouping();

	int addStroke(const RawStroke* stroke);
	int addStrokes(const RawStrokeGroup& strokeGroup);
	int removeStroke(const RawStroke* stroke);

	void forceGrouping(const RawStroke** strokes, size_t nstrokes);

	int getConfidence(const Matrix* m, double& confidence) const;

private:
	/* Compute which row and column the element should be in, and whether it should be a new intermediate row or column
	 * Returns whether we need a new cell or not
	 */
	bool determineRowCol(MatrixElement& elem, size_t& row, size_t& col, bool& newIntermRow, bool& newIntermCol);

	/* Determine the element that a stroke belongs to, possibly creating a new element. This element is stored
	 * via the strokeToElement map. This function updates the elementsToProcess vector with the relevant
	 * elements to process, as well as the strokeToElem map for new or changed elements
	 */
	void determineElement(const RawStroke* stroke);

	/* Analyses the elementsToProcess vector to determine which elements need processing. This allows incremental parsing,
	 * as opposed to an algorithm that re-parses the whole matrix every time. elementsToProcess will be cleared at the
	 * end of this function's execution.
	 */
	void parseRowsCols();

	// Segment stroke(s) into matrix elements to be parsed into rows and columns
	void segment(std::vector<const RawStroke*>& newStrokes);
	void segment(const RawStroke* stroke, bool isBeingAdded=true);

	// Calculate the scope of a stroke's bounding rectangle
	void calcScopeRect(Rect<long>& boundingRect, long offset);

	/* Merge the obsoleteElem into the existingElem. If successful, elementsToProcess will contain a REMOVE action on
	 * obsoleteElem and the strokeToElem map will be updated for the merged strokes.
	 */
	bool mergeMatrixElements(MatrixElement* existingElem, MatrixElement* obsoleteElem);
};

/* Treat each stroke as a mass that is capable of attracting other masses (strokes) through
 * a gravitational force. The force diminishes exponentially as the radius increases. Highly
 * attracted strokes are considered a grouping.
 */
/*class GravitationalGrouping : public ElemGroupingStrat
{
	struct Point { unsigned x, y; };

	struct GravElement {
		double mass;
		Rect<long> box;
		unsigned elemNum;
		Point centreOfMass;
	};

	class GravitationalModel
	{
		std::vector<GravElement> elems;

	public:
		void addGravitationalElement(GravElement ge);

	private:
		void calcRowForces();
	};

	GravitationalModel gravModel;

public:
	GravitationalGrouping(MathRecognizer* parser, std::vector<Matrix*>& matrices);
	virtual ~GravitationalGrouping();

	int addStroke(const RawStroke* stroke);
	int addStrokes(const RawStrokeGroup& strokeGroup);
	int removeStroke(const RawStroke* stroke);

	int getConfidence(const Matrix* m, double& confidence) const;

private:
	void addGravitationalMass(size_t numPoints, Rect<long> area);
};
*/

/* As per "Analyzing Sketch Content Using In-Air Packet Information" */
class InAirGrouping : public ElemGroupingStrat
{
	// % of matrix width when an in-air stroke should be defined as a row-break
	static const double ROW_BREAK_TOL;
	// % of matrix height when an in-air stroke should be defined as a col-break
	static const double COL_BREAK_TOL;
	// alpha_row, alpha_col as defined in the paper
	static const double ALPHA_ROW;
	static const double ALPHA_COL;

	// hybrid constants
	static const double ALPHA_DIST;
	static const double SPAN_OVERLAP_ON_ELEM_UPDATE_TOL; // same as MIN_OVERLAPPING_PROPORTION, except applied when an element is being updated (not created). this value should be > MIN_OVERLAPPING_PROPORTION because it should have a larger boundary to overcome if it flows into another row/column
	static const double MIN_OVERLAPPING_PROPORTION; // minimum overlap proportion (as a %) between an element's width/height and the row/column span for it to be part of that row/column 
	static const double MIN_CHAR_HEIGHT_INCHES; // min heigh in inches (to calculate MIN_CHAR_HEIGHT)
	static long MIN_CHAR_HEIGHT; // minimum height (in pixels) for a stroke to constitute a character (omit elements such as '-')

	// temporally-ordered strokes
	std::vector<const RawStroke*> strokes;

	// Row/Col breakpoints
	struct BreakPoint {
		BreakPoint(size_t inAirIndex) : inAirIndex(inAirIndex) {}

		// the nth in air stroke that defines the break
		size_t inAirIndex;

		bool operator<(const BreakPoint& other) const { return inAirIndex < other.inAirIndex; }
	};

	// (For hybrid segmentation)
	// Defines the height or width span of a group of MatrixElement objects (i.e. the combined union of their width or height)
	class Span {
	public:
		enum SpanType { HEIGHT, WIDTH };
	private:
		/* Add a certain percent of the average height of the span's elements in ink pixels to the start and end of the span as a buffer
		 * if the span is very narrow (e.g. for the number '1')
		 */
		static const double BUFFER_PCT;

		std::vector<MatrixElement*> elems; // elements comprised in this span
		std::vector<SpanType> spanTypes; // these two vectors are concurrently accessed i.e. elem[i] has span type spanTypes[i]
	public:
		Span();
		Span(MatrixElement* elem, SpanType spanType);

		void addElement(MatrixElement* elem, SpanType spanType);
		void removeElement(MatrixElement* elem);
		bool isEmpty() const;

		// TODO assumes that the left-hand side of the interval doesn't change after the first stroke (because of sorted_insert) is this ok?
		Interval getSpanInterval() const;

		inline bool operator==(const Span& other) const { return this->getSpanInterval() == other.getSpanInterval(); }
		inline bool operator<=(const Span& other) const { return this->getSpanInterval() <= other.getSpanInterval(); }
		inline bool operator<(const Span& other) const { return this->getSpanInterval() < other.getSpanInterval(); }
	};

	std::vector<Span> rowSpan, colSpan; // spans of the rows and columns

public:
	enum StatTypes { STANDARD_DEVIATION, MEAN };

	InAirGrouping(MathRecognizer* parser, std::vector<Matrix*>& matrices, int caps);

	int addStroke(const RawStroke* stroke);
	int addStrokes(const RawStrokeGroup& strokeGroup);
	int removeStroke(const RawStroke* stroke);

	void forceGrouping(const RawStroke** strokes, size_t nstrokes);

	int getConfidence(const Matrix* m, double& confidence) const;

private:
	// for findMinorBreaks()
	enum MajorBreakType { ROW, COL, HYBRID };

	void analyse(const std::vector<const RawStroke*> &strokes);

	// Row-major row segmentation
	std::vector<BreakPoint> findMajorRowBreaks(const std::vector<long> &inAirStrokeWidths, long matrixWidth) const;
	// Column-major column segmentation
	std::vector<BreakPoint> findMajorColBreaks(const std::vector<long> &inAirStrokeHeights, long matrixHeight) const;

	// Row-major column segmentation
	std::vector<BreakPoint> findMinorColBreaks(const std::vector<long> &inAirStrokeWidths, const std::vector<long> &penDownStrokeWidths, const std::vector<BreakPoint> &rowBreakPoints) const;
	// Col-major row segmentation
	std::vector<BreakPoint> findMinorRowBreaks(const std::vector<long> &inAirStrokeHeights, const std::vector<long> &penDownStrokeHeights, const std::vector<BreakPoint> &colBreakPoints) const;

	// Generic algorithm used by findMinorRowBreaks and findMinorColBreaks
	std::vector<BreakPoint> findMinorBreaks(const std::vector<long> &inAirStrokeLengths, const std::vector<long> &penDownStrokeLengths, const std::vector<InAirGrouping::BreakPoint> &majorBreakPoints, MajorBreakType type) const;
	
	// Hybrid analysis (not row-major or column-major)
	std::vector<BreakPoint> findMinorHybridBreaks(const std::vector<long> &inAirStrokeLengths, const std::vector<long> &penDownStrokeLengths) const;

	// The "strokes" parameter in each case is the pen-down strokes (in-air strokes are implicitly calculated from these)
	std::vector<long> getStrokeWidths(const std::vector<const RawStroke*> &strokes) const;
	std::vector<long> getStrokeHeights(const std::vector<const RawStroke*> &strokes) const;
	std::vector<long> getInAirStrokeWidths(const std::vector<const RawStroke*> &strokes) const;
	std::vector<long> getInAirStrokeHeights(const std::vector<const RawStroke*> &strokes) const;

	// Returns a map with computed stats, retrievable by StatType values
	std::map<StatTypes, double> calcStats(const std::vector<long> &data) const;
	// Returns the indices in the data vector that contain the n largest data values
	std::vector<size_t> getNLargestOutlierIndices(size_t n, const std::vector<long> &data) const;
	// Compute (spatially) which row and column the element should be in, and whether it should be a new intermediate row or column
	bool determineRowCol(MatrixElement& elem, size_t& row, size_t& col, bool& newIntermRow, bool& newIntermCol);
};

}

#endif
