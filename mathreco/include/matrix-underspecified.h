#ifndef SCG_MATRIX_UNDERSPECIFIED_H
#define SCG_MATRIX_UNDERSPECIFIED_H

#include <vector>
#include <map>

#include "matrix.h"
#include "matrix-fullyspecified.h"
#include "matrix-observer.h"
#include "RangeExpander.h"
#include "matrix-graph.h"
#include "surrogate-tree.h"

namespace scg
{

/* A matrix that contains ellipses. The UnderspecifiedMatrix class will fill in the cells as needed, according to the progression
 * of the surrounding elements. If a Matrix instance is being converted to an UnderspecifiedMatrix, call the constructor that takes
 * a Matrix instance, and the UnderspecifiedMatrix will incorporate its information and delete it when its life is done.
 */
class UnderspecifiedMatrix : public Matrix, MatrixElementObserver
{
private:
	mutable std::string _str;
	FullyspecifiedMatrix* inheritedMatrix;
	bool isDirty;

	struct Region {
		enum Type {TRIANGULAR, UPPER_TRIANGULAR, LOWER_TRIANGULAR, RECTANGULAR, UNDETERMINED};

		Region() : type(UNDETERMINED) {}
		void addElem(MatrixElement* elem);
		bool operator==(const Region& other);

		Type type;

		std::vector<MatrixElement*> elems;

		inline bool operator<(const UnderspecifiedMatrix::Region& other) const { return elems < other.elems; }
	};
	std::vector<Region> regions;

	//horizontal, vertical, diagonal, anti-diagonal ellipses strokes (using maps so we can have O(1) retrieval)
	std::map<MatrixElement*, bool> hEllipses, vEllipses, dEllipses, adEllipses;
	//blank (or fill) expressions (replace '?' character for a region)
	std::map<std::pair<size_t,size_t>, SurrogateTree*> blankExprs;
	//the defined fill inside of a rectangular region, or beside a triangular region (if not defined, assume fill of 0)
	std::map<Region, SurrogateTree*> regionFill;
	//MatrixElements that actually represent a fill expression for a region
	std::map<MatrixElement*, bool> fillExprs;

	//maps the MatrixElement instance in the abstract matrix to its corresponding index in concreteMatrix
	std::map<MatrixElement*, std::pair<size_t, size_t> > abstractToConcreteIdx;
	
	//the accepted graph representing the matrix (ellipses are represented as edges joining distinct elements)
	Graph g;
	
	//the concrete matrix defined by the ellipses
	std::vector<std::vector<const ExpressionTree*> > concreteMatrix;
	//number of rows, cols in the concrete matrix
	size_t mRows, mCols;

	bool isValidMatrix;

	//used for implicitly filling a region with 0s
	SurrogateTree* zeroTree;

	bool isShortFormRetrieval;

public:
	enum EllipsisType { HORIZONTAL, VERTICAL, DIAGONAL, ANTIDIAGONAL };

	// 'm' will be deleted once the underspecified matrix's life is done
	UnderspecifiedMatrix(MathRecognizer* parser, Matrix* m = NULL);
	virtual ~UnderspecifiedMatrix();

	// Get the dimensions of the expanded matrix
	int getDimensions(size_t& rows, size_t& cols) const;
	// Get the dimensions of the short-form version
	int getShortFormDimensions(size_t& rows, size_t& cols) const;

	int getCell(size_t row, size_t col, ExpressionIterator*& results) const;
	void setRetrievalType(bool isShortForm); //are we retrieving the short form or long form of the matrix?

	// The stroke group must contain the 3 dots
	void addEllipsis(const RawStroke** ellipsis, EllipsisType type);
	// Return if the matrix thus far is a valid underspecified matrix
	bool isValid();

	const char* str() const;

	// The following methods' calls are passed on to the inherited matrix (and we observe any changes)
	void insert(MatrixElement& elem, size_t row, size_t col, bool newIntermediateRow=false, bool newIntermediateCol=false);
	void remove(MatrixElement& elem);
	void move(MatrixElement& elem, size_t newRow, size_t newCol);
	void clear();
	bool hasPlaceholders() const;
	int getPosition(const MatrixElement& elem, size_t& row, size_t& col) const;
	bool isEmpty() const;
	MatrixElement* operator() (const size_t row, const size_t col) const;

	/* Rebuild graph and concrete matrix (if the graph is valid) */
	void rebuild();

	void notify(MatrixElementEvent evt); //MatrixElementObserver

private:
	UnderspecifiedMatrix(const UnderspecifiedMatrix& other) {}
	UnderspecifiedMatrix& operator=(const UnderspecifiedMatrix& other);

	/* Create a graph that maintains concrete expressions as nodes and ellipsis as edges */
	Graph createGraph();
	/* Determine if the given graph is valid according to criteria that identifies a valid matrix construction */
	bool isValidGraph(Graph &g);
	/* Track regions (cycles) from the graph */
	std::vector<Region> findRegions(Graph& g);
	std::vector<Region> findRegions(Graph& g, std::vector<Graph::Node>& path, size_t numEdges, std::map<Graph::Edge, bool>& isInRegion);
	/* Build the matrix from graph "g" and store it for querying */
	void buildMatrix();
	/* Identify any expressions that are used to fill a range of the matrix */
	void identifyFillerExpressions(Graph& g);
	/* Store regions in the instance var, store blank expressions, and convert TRIANGULAR regions to UPPER or LOWER triangular */
	void analyseAndStoreRegions(std::vector<Region> regions, std::map<MatrixElement*, std::pair<size_t, size_t> > abstractToConcreteIdx);
	/* Get the fill for inside a rectangular region or beside a triangular region (see regionFill) */
	ExpressionTree* getRegionFill(size_t row, size_t col) const;

	/* Add an edge to 'g' according to the the previous and next elements of an ellipsis */
	bool addEdgeFromEllipsis(Graph& g, std::vector<MatrixElement*> prevElems, MatrixElement* nextElem, Graph::Edge::Type ellipsisType);
	// Getting the previous and next elements for the given ellipsis (horiz, vert, diag, anti-diag)
	std::vector<MatrixElement*> getHPrevElems(MatrixElement* ellipsis);
	std::vector<MatrixElement*> getVPrevElems(MatrixElement* ellipsis);
	std::vector<MatrixElement*> getDPrevElems(MatrixElement* ellipsis);
	std::vector<MatrixElement*> getADPrevElems(MatrixElement* ellipsis);
	MatrixElement* getHNextElem(MatrixElement* ellipsis);
	MatrixElement* getVNextElem(MatrixElement* ellipsis);
	MatrixElement* getDNextElem(MatrixElement* ellipsis);
	MatrixElement* getADNextElem(MatrixElement* ellipsis);

	/* Get the width of the matrix defined by 'g' */
	size_t getWidth(Graph &g);
	/* Get the width of the matrix defined by 'g' */
	size_t getHeight(Graph &g);
	/* Helper method for getHeight/getWidth - do not use directly */
	size_t getDimension(Graph& g, Graph::Node& node, std::map<Graph::Edge, bool>& visited, std::map<Graph::Node, size_t>& cache, bool isWidth);

	bool isEllipsis(MatrixElement* elem);
	bool isFillerExpr(MatrixElement* elem);
	const ExpressionTree* getExprTree(MatrixElement* elem);

	/* Delete the concrete matrix */
	void deleteMatrix();
};


}

#endif
