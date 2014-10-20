#include "matrix-underspecified.h"

#include <algorithm>
#include "grammar-values.h"
#include "MathRecognizer.h"
#include "interval.h"


namespace scg
{


void UnderspecifiedMatrix::Region::addElem(MatrixElement* elem)
{
	elems.push_back(elem);

	const int TRIANGLE_EDGES = 3;
	const int RECT_EDGES = 4;

	if (elems.size() == TRIANGLE_EDGES) {
		type = TRIANGULAR;
	} else if (elems.size() == RECT_EDGES) {
		type = RECTANGULAR;
	}
}

bool UnderspecifiedMatrix::Region::operator ==(const UnderspecifiedMatrix::Region &other)
{
	if (elems.size() != other.elems.size()) return false;

	for (std::vector<MatrixElement*>::iterator it = elems.begin(); it != elems.end(); it++) {
		MatrixElement* elem = *it;
		bool found = false;
		for (std::vector<MatrixElement*>::const_iterator it1 = other.elems.begin(); it1 != other.elems.end(); it1++) {
			found |= (*it1 == elem);
		}
		if (!found) return false;
	}
	return true;
}

UnderspecifiedMatrix::UnderspecifiedMatrix(MathRecognizer* parser, Matrix* m)
: Matrix(parser), isShortFormRetrieval(false), mRows(0), mCols(0), isValidMatrix(false), isDirty(false)
{
	if ((inheritedMatrix = dynamic_cast<FullyspecifiedMatrix*>(m)) == 0) {
		inheritedMatrix = new FullyspecifiedMatrix(parser);
	}

	zeroTree = new SurrogateNode(NUM_EXPR);
	SurrogateTree* terminal = new SurrogateNode(TERMINAL_EXPR);
	terminal->setTerminalValue("0");
	zeroTree->addChild(terminal);
}

UnderspecifiedMatrix::~UnderspecifiedMatrix()
{
	deleteMatrix();
	delete inheritedMatrix;
	delete zeroTree;

	for (std::map<Region, SurrogateTree*>::iterator it = regionFill.begin(); it != regionFill.end(); it++) {
		delete it->second;
	}
}

void UnderspecifiedMatrix::deleteMatrix()
{
	if (concreteMatrix.size() > 0) {
		for (size_t i = 0; i < concreteMatrix.size(); i++) {
			for (size_t j = 0; j < concreteMatrix[i].size(); j++) {
				if (concreteMatrix[i][j] != NULL) {
					delete concreteMatrix[i][j];
				}
			}
			concreteMatrix[i].clear();
		}
		concreteMatrix.clear();
	}

	regions.clear();

	for (std::map<std::pair<size_t,size_t>, SurrogateTree*>::iterator it = blankExprs.begin(); it != blankExprs.end(); it++) {
		delete it->second;
	}
	blankExprs.clear();
}

void UnderspecifiedMatrix::notify(MatrixElementEvent evt)
{
	isDirty = true;
}

void UnderspecifiedMatrix::insert(MatrixElement& elem, size_t row, size_t col, bool newIntermediateRow, bool newIntermediateCol)
{
	inheritedMatrix->insert(elem, row, col, newIntermediateRow, newIntermediateCol);
	elem.registerObserver(*this);
	isDirty = true;
}

void UnderspecifiedMatrix::remove(MatrixElement& elem)
{
	inheritedMatrix->remove(elem);
	elem.removeObserver(*this);
	isDirty = true;
}

void UnderspecifiedMatrix::move(MatrixElement& elem, size_t newRow, size_t newCol)
{
	inheritedMatrix->move(elem, newRow, newCol);
	isDirty = true;
}

void UnderspecifiedMatrix::clear()
{
	size_t rows, cols;
	inheritedMatrix->getDimensions(rows, cols);
	for (size_t r = 0; r < rows; r++) {
		for (size_t c = 0; c < cols; c++) {
			MatrixElement* elem = (*inheritedMatrix)(r,c);
			if (elem != NULL) {
				elem->removeObserver(*this);
			}
		}
	}

	deleteMatrix();

	hEllipses.clear();
	vEllipses.clear();
	dEllipses.clear();
	adEllipses.clear();
	for (std::map<Region, SurrogateTree*>::iterator it = regionFill.begin(); it != regionFill.end(); it++) {
		delete it->second;
	}
	regionFill.clear();
	fillExprs.clear();
	abstractToConcreteIdx.clear();
	g = Graph();
	mRows = 0;
	mCols = 0;
	isValidMatrix = false;

	inheritedMatrix->clear();
}

bool UnderspecifiedMatrix::hasPlaceholders() const
{
	return inheritedMatrix->hasPlaceholders();
}

int UnderspecifiedMatrix::getPosition(const MatrixElement& elem, size_t& row, size_t& col) const
{
	return inheritedMatrix->getPosition(elem, row, col);
}

MatrixElement* UnderspecifiedMatrix::operator() (const size_t row, const size_t col) const
{
	return inheritedMatrix->operator ()(row, col);
}

bool UnderspecifiedMatrix::isEmpty() const
{
	return inheritedMatrix->isEmpty();
}

int UnderspecifiedMatrix::getCell(size_t row, size_t col, ExpressionIterator*& results) const
{
	results = NULL;

	if (isShortFormRetrieval) {
		std::pair<size_t, size_t> coords = std::make_pair(row, col);
		if (blankExprs.find(coords) != blankExprs.end()) {
			SurrogateTree* fillTree = blankExprs.find(coords)->second;
			if (fillTree == NULL) {
				results = parser->CreateIteratorForExpression(CreateBlankExpression(), true, false);
			}
			else {
				results = parser->CreateIteratorForExpression(ConvertSurrogateTreeToExpressionTree(fillTree), true, false);
			}
		} else {
			inheritedMatrix->getCell(row, col, results);
		}
	} else {
		if (row < mRows && col < mCols) {
			if (concreteMatrix[row][col] != NULL) {
				//const raw_expression_node* node = new const raw_expression_node(*dynamic_cast<const raw_expression_node*>(concreteMatrix[row][col]));
				const ExpressionTree *node = concreteMatrix[row][col];
				VERBOSE(*verb_out << "\tUnderspecifiedMatrix: (" << row << "," << col << ") = " <<  node->long_str() << std::endl);
				results = parser->CreateIteratorForExpression(node, false, false);
			} else {
				ExpressionTree* t = getRegionFill(row, col);
				if (t != NULL) {
					results = parser->CreateIteratorForExpression(t, false, false);
				}
			}
		}
	}
	
	if (results == NULL) {
		VERBOSE(*verb_out << "\tUnderspecifiedMatrix: could not find range for (" << row << "," << col << "). defaulting to Matrix::getCell" << std::endl);
		inheritedMatrix->getCell(row, col, results);
	}

	return (results == NULL) ? E_NOTFOUND : 0;
}

void UnderspecifiedMatrix::setRetrievalType(bool isShortForm)
{
	isShortFormRetrieval = isShortForm;
}

int UnderspecifiedMatrix::getDimensions(size_t& rows, size_t& cols) const
{
	if (!isShortFormRetrieval) {
		rows = mRows;
		cols = mCols;
	} else {
		inheritedMatrix->getDimensions(rows, cols);
	}

	return 0;
}

int UnderspecifiedMatrix::getShortFormDimensions(size_t& rows, size_t& cols) const
{
	return inheritedMatrix->getDimensions(rows, cols);
}

const char* UnderspecifiedMatrix::str() const
{
	std::stringstream ss;

	size_t rows, cols;
	getDimensions(rows, cols);
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			ss << i << "," << j;

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

	_str = ss.str();
	const char* output = _str.c_str();
	return output;
}

void UnderspecifiedMatrix::addEllipsis(const RawStroke** ellipsis, EllipsisType type)
{
	VERBOSE(*verb_out << "\tUnderspecifiedMatrix: add ellipsis");

	size_t rows, cols;
	inheritedMatrix->getDimensions(rows, cols);
	// Store the MatrixElement that contains the ellipsis in the correct set
	for (size_t r = 0; r < rows; r++) {
		for (size_t c = 0; c < cols; c++) {
			MatrixElement* elem = (*inheritedMatrix)(r, c);
			if (elem == NULL) continue;

			const RawStroke* stk = elem->getStrokes()[0];
			if (stk == ellipsis[0] || stk == ellipsis[1] || stk == ellipsis[2]) {
				if (type == HORIZONTAL) hEllipses[elem] = true;
				else if (type == VERTICAL) vEllipses[elem] = true;
				else if (type == DIAGONAL) dEllipses[elem] = true;
				else if (type == ANTIDIAGONAL) adEllipses[elem] = true;
				
				isDirty = true;
				return;
			}
		}
	}
}

bool UnderspecifiedMatrix::isValid()
{
	if (mRows == 0 || mCols == 0) return false;
	return isValidMatrix;

	/*for (int i = 0; i < mRows; i++) {
		for (int j = 0; j < mCols; j++) {
			if (concreteMatrix[i][j] == NULL && getRegionFill(i, j) == NULL) return false;
		}
	}

	return true;*/
}

void UnderspecifiedMatrix::rebuild()
{
	if (!isDirty) {
		return;
	}

	for (std::map<Region, SurrogateTree*>::iterator it = regionFill.begin(); it != regionFill.end(); it++) {
		delete it->second;
	}
	regionFill.clear();

	Graph graph = createGraph();
	identifyFillerExpressions(graph);
	if (isValidGraph(graph)) {
		this->g = graph;
		buildMatrix();
	} else {
		VERBOSE(*verb_out << "\tNot valid graph" << std::endl);
	}


	isDirty = false;
}

void UnderspecifiedMatrix::identifyFillerExpressions(Graph& g)
{
	size_t rows, cols;
	inheritedMatrix->getDimensions(rows, cols);

	// Check to see if this is a filler element
	std::vector<Region> regions = findRegions(g);
	for (std::vector<Region>::iterator it = regions.begin(); it != regions.end(); it++) {
		Region region = *it;

		// Get upper-left and lower-right coordinates
		size_t c1 = std::numeric_limits<size_t>::max();
		size_t c2 = 0;
		size_t r1 = std::numeric_limits<size_t>::max();
		size_t r2 = 0;
		for (size_t i = 0; i < region.elems.size(); i++) {
			size_t r,c;
			inheritedMatrix->getPosition(*region.elems[i], r, c);
			c1 = std::min(c1, c);
			c2 = std::max(c2, c);
			r1 = std::min(r1, r);
			r2 = std::max(r2, r);
		}

		bool isUpperTriangular = false;
		switch (region.type) {
			case Region::TRIANGULAR:
				// Determine if it's upper or lower triangular
				for (size_t i = 0; i < region.elems.size(); i++) {
					size_t r, c;
					inheritedMatrix->getPosition(*region.elems[i], r, c);
					if ((r == r1 && c == c2) || (r == r2 && c == c1)) {
						isUpperTriangular = (r == r1);
						break;
					}
				}

				if (isUpperTriangular) {
					// it would be to the left of the triangle
					for (size_t row = 1; row < rows; row++) {
						for (size_t col = 0; col < c1 + (row - r1); col++) {
							if ((*inheritedMatrix)(row, col) != NULL) {
								MatrixElement* elem = (*inheritedMatrix)(row, col);
								const ExpressionTree* t = getExprTree(elem);
								regionFill[region] = expressionToSurrogateTree(t);
								fillExprs[elem] = true;
							}
						}
					}
				} else {
					// it would be to the right of the triangle
					for (size_t row = 0; row < rows; row++) {
						for (size_t col = c1 + (row - r1) + 1; col < cols; col++) {
							if ((*inheritedMatrix)(row, col) != NULL) {
								MatrixElement* elem = (*inheritedMatrix)(row, col);
								const ExpressionTree* t = getExprTree(elem);
								regionFill[region] = expressionToSurrogateTree(t);
								fillExprs[elem] = true;
							}
						}
					}
				}
				break;
			case Region::RECTANGULAR:
				// it would be in the middle of the rectangular region
				if ((*inheritedMatrix)(r1 + 1, c1 + 1) != NULL) {
					MatrixElement* elem = (*inheritedMatrix)(r1 + 1, c1 + 1);
					const ExpressionTree* t = getExprTree(elem);
					regionFill[region] = expressionToSurrogateTree(t);
					fillExprs[elem] = true;
				}
				break;
		}
	}
}

Graph UnderspecifiedMatrix::createGraph()
{
	Graph g;
	size_t rows, cols;
	inheritedMatrix->getDimensions(rows, cols);
	
	// Add all elements except ellipses as nodes. Ellipses will be represented as edges between nodes
	for (size_t r = 0; r < rows; r++) {
		for (size_t c = 0; c < cols; c++) {
			MatrixElement* elem = (*inheritedMatrix)(r, c);
			if (elem == NULL || isEllipsis(elem)) continue;
			
			Graph::Node node(elem);
			g.addNode(node);
		}
	}

	// We will need to analyse these after
	std::vector<std::pair<Graph::Node*, Graph::Node*> > zeroWeights;
	std::map<size_t, size_t> lineWeight;

	// For each ellipsis type, find its bounding elements
	for (std::map<MatrixElement*, bool>::iterator it = hEllipses.begin(); it != hEllipses.end(); it++) {
		MatrixElement* ellipsis = it->first;
		std::vector<MatrixElement*> prevElems = getHPrevElems(ellipsis);
		MatrixElement* next = getHNextElem(ellipsis);
		if (!addEdgeFromEllipsis(g, prevElems, next, Graph::Edge::HORIZONTAL)) {
			return Graph();
		}

		// If it's an ambiguous edge (e.g 1...1), store it to analyse below
		Graph::Node* start = g.elemToNode(prevElems.back());
		Graph::Node* end = g.elemToNode(next);
		size_t weight = g.getEdge(*start, *end).getWeight();
		if (weight == 0) {
			zeroWeights.push_back(std::make_pair(start, end));
		} else {
			size_t erow, ecol;
			inheritedMatrix->getPosition(*ellipsis, erow, ecol);
			lineWeight[ecol] = weight;
		}
	}
	for (std::vector<std::pair<Graph::Node*, Graph::Node*> >::iterator it = zeroWeights.begin(); it != zeroWeights.end(); it++) {
		Graph::Node* start = it->first;
		Graph::Node* end = it->second;
		
		size_t erow, ecol;
		inheritedMatrix->getPosition(*g.nodeToElem(*start), erow, ecol);
		ecol += 1; //the ellipsis is 1 column to the right of the starting expression
		if (lineWeight.find(ecol) != lineWeight.end()) {
			g.setEdgeWeight(*start, *end, lineWeight[ecol]);
		}
	}

	zeroWeights.clear();

	for (std::map<MatrixElement*, bool>::iterator it = vEllipses.begin(); it != vEllipses.end(); it++) {
		MatrixElement* ellipsis = it->first;
		std::vector<MatrixElement*> prevElems = getVPrevElems(ellipsis);
		MatrixElement* next = getVNextElem(ellipsis);
		if (!addEdgeFromEllipsis(g, prevElems, next, Graph::Edge::VERTICAL)) {
			return Graph();
		}

		// If it's an ambiguous edge (e.g 1...1), store it to analyse below
		Graph::Node* start = g.elemToNode(prevElems.back());
		Graph::Node* end = g.elemToNode(next);
		size_t weight = g.getEdge(*start, *end).getWeight();
		if (weight == 0) {
			zeroWeights.push_back(std::make_pair(start, end));
		} else {
			size_t erow, ecol;
			inheritedMatrix->getPosition(*ellipsis, erow, ecol);
			lineWeight[ecol] = weight;
		}
	}
	for (std::vector<std::pair<Graph::Node*, Graph::Node*> >::iterator it = zeroWeights.begin(); it != zeroWeights.end(); it++) {
		Graph::Node* start = it->first;
		Graph::Node* end = it->second;
		
		size_t erow, ecol;
		inheritedMatrix->getPosition(*g.nodeToElem(*start), erow, ecol);
		erow += 1; //the ellipsis is 1 row down from the starting expression
		if (lineWeight.find(erow) != lineWeight.end()) {
			g.setEdgeWeight(*start, *end, lineWeight[erow]);
		}
	}

	for (std::map<MatrixElement*, bool>::iterator it = dEllipses.begin(); it != dEllipses.end(); it++) {
		MatrixElement* ellipsis = it->first;
		std::vector<MatrixElement*> prevElems = getDPrevElems(ellipsis);
		MatrixElement* next = getDNextElem(ellipsis);
		if (!addEdgeFromEllipsis(g, prevElems, next, Graph::Edge::DIAGONAL)) {
			return Graph();
		}
	}

	for (std::map<MatrixElement*, bool>::iterator it = adEllipses.begin(); it != adEllipses.end(); it++) {
		MatrixElement* ellipsis = it->first;
		std::vector<MatrixElement*> prevElems = getADPrevElems(ellipsis);
		MatrixElement* next = getADNextElem(ellipsis);
		if (!addEdgeFromEllipsis(g, prevElems, next, Graph::Edge::ANTIDIAGONAL)) {
			return Graph();
		}
	}

	return g;
}

bool UnderspecifiedMatrix::addEdgeFromEllipsis(Graph& g, std::vector<MatrixElement*> prevElems, MatrixElement* nextElem, Graph::Edge::Type ellipsisType)
{
	RangeExpander rangeExpander(parser);

	if (prevElems.size() == 0 || nextElem == NULL) return false;

	for (std::vector<MatrixElement*>::iterator it = prevElems.begin(); it != prevElems.end(); it++) {
		MatrixElement* prevElem = *it;
		const ExpressionTree* startExpr = getExprTree(prevElem);
		if (startExpr != NULL) {
			rangeExpander.appendStartExpression(startExpr);
		}
	}

	const ExpressionTree* endExpr = getExprTree(nextElem);
	if (endExpr != NULL) {
		rangeExpander.setEndExpression(endExpr);
	}

	// Remove starting elements until we have a valid range (e.g. 1 3 ... 6 would remove the first element)
	while (!rangeExpander.isValid() && rangeExpander.getStartExpressions().size() > 1) {
		rangeExpander.removeStartExpression(0);
	}

	if (!rangeExpander.isValid()) return false;

	Graph::Node& first = *(g.elemToNode(prevElems.back()));
	Graph::Node& last = *(g.elemToNode(nextElem));
	g.addEdge(first, last, ellipsisType, rangeExpander);
	return true;
}

bool UnderspecifiedMatrix::isValidGraph(Graph &g)
{
	/* TODO consistency checks
	 * Rule 1: any path from the starting node returns the same height/width
	 * Rule 2: each diagonal has to have the same number of rows and columns
	 * Rule 3: no node has 2 diagonal edges
	 */
	findRegions(g); //TODO temporary - just used to analyse 0-weighted edges (e.g. 1...1)
	return !g.isEmpty() && (getWidth(g) > 0 && getHeight(g) > 0);
}

std::vector<UnderspecifiedMatrix::Region> UnderspecifiedMatrix::findRegions(Graph& g)
{
	std::vector<Region> ret;

	// We know that the cycles in the graph will only be length 3 or 4 (triangular or rectangular regions)
	// Furthermore, since this graph has been validated, we know that any node with >1 edges is part of a cycle
	// Theorem: any edge can only be part of at most 1 region

	const int TRIANGLE_EDGES = 3;
	const int RECT_EDGES = 4;

	std::vector<Graph::Node> path;
	std::map<Graph::Edge, bool> isInRegion;

	std::vector<Graph::Node> nodes = g.getNodes();
	for (std::vector<Graph::Node>::iterator it = nodes.begin(); it != nodes.end(); it++) {
		path.push_back(*it);
		std::vector<Region> regions = findRegions(g, path, RECT_EDGES, isInRegion);
		path.pop_back();
		ret.insert(ret.end(), regions.begin(), regions.end());
	}

	for (std::vector<Graph::Node>::iterator it = nodes.begin(); it != nodes.end(); it++) {
		path.push_back(*it);
		std::vector<Region> regions = findRegions(g, path, TRIANGLE_EDGES, isInRegion);
		path.pop_back();
		ret.insert(ret.end(), regions.begin(), regions.end());

		// All triangular regions must have the same sized edges
		for (std::vector<Region>::iterator it = regions.begin(); it != regions.end(); it++) {
			Region region = *it;
			Graph::Node* n1 = g.elemToNode(region.elems[0]);
			Graph::Node* n2 = g.elemToNode(region.elems[1]);
			Graph::Node* n3 = g.elemToNode(region.elems[2]);
			
			Graph::Edge triEdges[3] = { g.getEdge(*n1, *n2), g.getEdge(*n1, *n3), g.getEdge(*n2, *n3) };
			size_t weight = std::max(triEdges[0].getWeight(), std::max(triEdges[1].getWeight(), triEdges[2].getWeight()));

			if (triEdges[0].getWeight() == 0) g.setEdgeWeight(*n1, *n2, weight);
			if (triEdges[1].getWeight() == 0) g.setEdgeWeight(*n1, *n3, weight);
			if (triEdges[2].getWeight() == 0) g.setEdgeWeight(*n2, *n3, weight);
		}
	}

	return ret;
}

std::vector<UnderspecifiedMatrix::Region> UnderspecifiedMatrix::findRegions(Graph& g, std::vector<Graph::Node>& path, size_t regionEdges, std::map<Graph::Edge, bool>& isInRegion)
{
	std::vector<Region> ret;

	// If we would be at the end of a cycle, check if there is an edge connecting the first and last elements to finalize the region
	if (regionEdges == path.size()) {
		Graph::Node& startNode = path.front();
		Graph::Node& endNode = path.back();

		std::vector<Graph::Edge> edges = g.getEdges(startNode);
		for (std::vector<Graph::Edge>::iterator eit = edges.begin(); eit != edges.end(); eit++) {
			Graph::Edge edge = *eit;
			Graph::Node& otherNode = edge.otherNode(startNode);
			if (otherNode == endNode) {
				Region r;
				for (std::vector<Graph::Node>::iterator it = path.begin(); it != path.end(); it++) {
					r.addElem(g.nodeToElem(*it));
				}
				ret.push_back(r);
				isInRegion[edge] = true;
				return ret;
			}
		}

		return ret;
	}

	// Check every edge that isn't already in a region to see if it eventually forms one
	Graph::Node curNode = path.back();
	std::vector<Graph::Edge> edges = g.getEdges(curNode);
	for (std::vector<Graph::Edge>::iterator eit = edges.begin(); eit != edges.end(); eit++) {
		Graph::Edge edge = *eit;

		if (isInRegion.find(edge) != isInRegion.end()) continue;

		Graph::Node& otherNode = edge.otherNode(curNode);
		if (path.size() > 1 && otherNode == *(path.end() - 2)) continue; //don't reverse

		path.push_back(otherNode);
		std::vector<Region> regions = findRegions(g, path, regionEdges, isInRegion);
		path.pop_back();

		if (regions.size() > 0) {
			isInRegion[edge] = true;
			ret.insert(ret.end(), regions.begin(), regions.end());
		}
	}

	return ret;
}

void UnderspecifiedMatrix::buildMatrix()
{
	isValidMatrix = true;
	size_t numElemsCreated = 0;

	std::vector<Region> regions = findRegions(g);

	// Clear previous matrix
	deleteMatrix();

	// Construct new one
	mRows = getHeight(g);
	mCols = getWidth(g);
	if (mRows <= 0 || mCols <= 0) {
		isValidMatrix = false;
		return;
	}

	//concreteMatrix.assign(mRows, std::vector<const ExpressionTree*>(mCols, NULL));
	for (size_t r = 0; r < mRows; r++) {
		std::vector<const ExpressionTree*> nullCols;    
		for (size_t c = 0; c < mCols; c++) {
			nullCols.push_back(NULL);
		}
		concreteMatrix.push_back(nullCols);
	}

	size_t curRow = 0, curCol = 0;
	
	size_t rows, cols;
	inheritedMatrix->getDimensions(rows, cols);
	if (mRows < rows || mCols < cols) return; //TODO temp (while isValidGraph() isn't done)
	std::map<Graph::Edge, bool> visited;
	std::map<MatrixElement*, std::pair<size_t, size_t> > abstractToConcreteIdx; //maps the elements from the abstract matrix to the expressions in the concrete matrix
	// TODO somewhat unsafe (does not check for out-of-bounds insertion)
	for (size_t r = 0; r < rows; r++) {
		for (size_t c = 0; c < cols; c++) {
			MatrixElement* elem = (*inheritedMatrix)(r, c);
			if (elem == NULL || isEllipsis(elem) || isFillerExpr(elem)) continue;

			Graph::Node& startNode = *(g.elemToNode(elem));
			std::vector<Graph::Edge> edges = g.getEdges(startNode);

			if (abstractToConcreteIdx.find(elem) != abstractToConcreteIdx.end()) {
				std::pair<size_t, size_t> elemIndex = abstractToConcreteIdx[elem];
				curRow = elemIndex.first;
				curCol = elemIndex.second;
			} else {
				// keep track of the row and col of an expression in the concrete matrix
				concreteMatrix[curRow][curCol] = getExprTree(elem);
				numElemsCreated++;
				abstractToConcreteIdx[elem] = std::make_pair(curRow, curCol);
			}

			// Expand each edge attached to this node
			for (std::vector<Graph::Edge>::iterator it = edges.begin(); it != edges.end(); it++) {
				Graph::Edge& edge = *it;
				if (visited.find(edge) != visited.end()) continue;

				Graph::Node& endNode = edge.otherNode(startNode);
				MatrixElement* endElem = g.nodeToElem(endNode);
				const ExpressionTree* endTree = getExprTree(endElem);

				size_t end;
				size_t weight = edge.getWeight();
				switch (edge.type)
				{
				case Graph::Edge::HORIZONTAL:
					end = curCol + edge.getWeight() + 1;
					concreteMatrix[curRow][end] = endTree;
					numElemsCreated++;
					abstractToConcreteIdx[endElem] = std::make_pair(curRow, end);

					for (size_t i = 0; i < weight; i++) {
						concreteMatrix[curRow][curCol + i + 1] = edge.expand(i);
					}
					numElemsCreated += weight;
					break;

				case Graph::Edge::VERTICAL:
					end = curRow + edge.getWeight() + 1;
					concreteMatrix[end][curCol] = endTree;
					numElemsCreated++;
					abstractToConcreteIdx[endElem] = std::make_pair(end, curCol);

					for (size_t i = 0; i < weight; i++) {
						concreteMatrix[curRow + i + 1][curCol] = edge.expand(i);
					}
					numElemsCreated += weight;
					break;

				case Graph::Edge::DIAGONAL:
					end = edge.getWeight() + 1;
					concreteMatrix[curRow + end][curCol + end] = endTree;
					numElemsCreated++;
					abstractToConcreteIdx[endElem] = std::make_pair(curRow + end, curCol + end);

					for (size_t i = 0; i < weight; i++) {
						concreteMatrix[curRow + i + 1][curCol + i + 1] = edge.expand(i);
					}
					numElemsCreated += weight;
					break;

				case Graph::Edge::ANTIDIAGONAL:
					// TODO
					end = edge.getWeight() + 1;
					concreteMatrix[curRow + end][curCol - end] = endTree;
					numElemsCreated++;
					abstractToConcreteIdx[endElem] = std::make_pair(curRow + end, curCol - end);

					for (size_t i = 0; i < weight; i++) {
						concreteMatrix[curRow + i + 1][curCol - i - 1] = edge.expand(edge.getWeight() - i - 1);
					}
					numElemsCreated += weight;
					break;
				};

				visited[edge] = true;
			}

			curCol++;
		}

		curCol = 0;
		curRow++;
	}

	analyseAndStoreRegions(regions, abstractToConcreteIdx);

	if (regions.size() > 0) {
		isValidMatrix = true; //TODO assumption (since we got this far, this is a fairly strong assumption)
	} else {
		// this is a valid matrix if the number of expanded elements indeed matches the number of cells in the matrix
		isValidMatrix = (numElemsCreated == (mRows * mCols));
	}
}

void UnderspecifiedMatrix::analyseAndStoreRegions(std::vector<Region> regions, std::map<MatrixElement*, std::pair<size_t, size_t> > abstractToConcreteIdx)
{
	for (std::vector<Region>::iterator it = regions.begin(); it != regions.end(); it++) {
		Region region = *it;

		size_t c1 = std::numeric_limits<size_t>::max();
		size_t c2 = 0;
		size_t r1 = std::numeric_limits<size_t>::max();
		size_t r2 = 0;
		
		if (region.type == Region::TRIANGULAR || region.type == Region::UPPER_TRIANGULAR || region.type == Region::LOWER_TRIANGULAR) {
			for (size_t i = 0; i < region.elems.size(); i++) {
				std::pair<size_t,size_t> coords = abstractToConcreteIdx[region.elems[i]];
				c1 = std::min(c1, coords.second);
				c2 = std::max(c2, coords.second);
				r1 = std::min(r1, coords.first);
				r2 = std::max(r2, coords.first);
			}

			// If it hasn't been specified, specify it as upper or lower triangular
			if (region.type == Region::TRIANGULAR) {
				// Determine if it's upper or lower triangular
				bool isUpperTriangular = false;
				for (size_t i = 0; i < region.elems.size(); i++) {
					std::pair<size_t,size_t> coords = abstractToConcreteIdx[region.elems[i]];
					if ((coords.first == r1 && coords.second == c2) || (coords.first == r2 && coords.second == c1)) {
						isUpperTriangular = (coords.first == r1);
						break;
					}
				}

				region.type = (isUpperTriangular) ? Region::UPPER_TRIANGULAR : Region::LOWER_TRIANGULAR;
				this->regions.push_back(region);

				const size_t NUM_ROWS_IN_FILL = r2 - r1;
				c1 = r1 = std::numeric_limits<size_t>::max();
				c2 = r2 = 0;
				// determine the upper-left and lower-right coordinates
				for (size_t i = 0; i < region.elems.size(); i++) {
					MatrixElement* elem = region.elems[i];
					size_t row, col;
					inheritedMatrix->getPosition(*elem, row, col);
					c1 = std::min(c1, col);
					c2 = std::max(c2, col);
					r1 = std::min(r1, row);
					r2 = std::max(r2, row);
				}

				// track blank expressions for interface output
				if (region.type == Region::UPPER_TRIANGULAR) {
					for (size_t i = 1; i < NUM_ROWS_IN_FILL + 1; i++) {
						size_t r = r1 + i;
						for (size_t c = c1; c < r - r1; c++) {
							blankExprs[std::make_pair(r, c)] = (regionFill.find(region) != regionFill.end()) ? regionFill[region]->copy(): NULL;
						}
					}
				} else {
					// lower triangular
					for (size_t i = 0; i < NUM_ROWS_IN_FILL; i++) {
						size_t r = r1 + i;
						for (size_t c = c1 + i + 1; c <= c2; c++) {
							blankExprs[std::make_pair(r, c)] = (regionFill.find(region) != regionFill.end()) ? regionFill[region]->copy(): NULL;
						}
					}
				}
			}
		} else if (region.type == Region::RECTANGULAR) {
			this->regions.push_back(region);

			c1 = r1 = std::numeric_limits<size_t>::max();
			c2 = r2 = 0;
			for (size_t i = 0; i < region.elems.size(); i++) {
				MatrixElement* elem = region.elems[i];
				size_t row, col;
				inheritedMatrix->getPosition(*elem, row, col);
				c1 = std::min(c1, col);
				c2 = std::max(c2, col);
				r1 = std::min(r1, row);
				r2 = std::max(r2, row);
			}
			// track blank expressions for interface output
			for (size_t r = r1 + 1; r < r2; r++) {
				for (size_t c = c1 + 1; c < c2; c++) {
					blankExprs[std::make_pair(r, c)] = (regionFill.find(region) != regionFill.end()) ? regionFill[region]->copy(): NULL;
				}
			}
		}
	}

	this->abstractToConcreteIdx = abstractToConcreteIdx;
}

ExpressionTree* UnderspecifiedMatrix::getRegionFill(size_t row, size_t col) const
{
	for (std::vector<Region>::const_iterator it = regions.begin(); it != regions.end(); it++) {
		Region region = *it;

		// top left and bottom right indices of the region
		size_t c1 = std::numeric_limits<size_t>::max();
		size_t c2 = 0;
		size_t r1 = std::numeric_limits<size_t>::max();
		size_t r2 = 0;
		
		if (region.type == Region::UPPER_TRIANGULAR) {
			for (size_t i = 0; i < region.elems.size(); i++) {
				std::pair<size_t,size_t> coords = abstractToConcreteIdx.find(region.elems[i])->second;
				c1 = std::min(c1, coords.second);
				c2 = std::max(c2, coords.second);
				r1 = std::min(r1, coords.first);
				r2 = std::max(r2, coords.first);
			}

			if (row >= r1 && row <= r2 && col >= c1 && col <= c2) {
				size_t startCol = c1 + (row - r1);
				size_t endCol = c2;

				// left side of the triangle is filled with zeros
				if (col < startCol) {
					SurrogateTree* t = (regionFill.find(region) != regionFill.end()) ? regionFill.find(region)->second : zeroTree;
					return ConvertSurrogateTreeToExpressionTree(t);
				}

				RangeExpander across(parser);
				across.appendStartExpression(concreteMatrix[row][startCol]);
				across.setEndExpression(concreteMatrix[row][endCol]);

				return across.at(col - startCol);
			}
		} else if (region.type == Region::LOWER_TRIANGULAR) {
			for (size_t i = 0; i < region.elems.size(); i++) {
				std::pair<size_t,size_t> coords = abstractToConcreteIdx.find(region.elems[i])->second;
				c1 = std::min(c1, coords.second);
				c2 = std::max(c2, coords.second);
				r1 = std::min(r1, coords.first);
				r2 = std::max(r2, coords.first);
			}

			if (row >= r1 && row <= r2 && col >= c1 && col <= c2) {
				size_t startCol = c1;
				size_t endCol = c1 + (row - r1);

				// right side of the triangle is filled with zeros
				if (col > endCol) {
					SurrogateTree* t = (regionFill.find(region) != regionFill.end()) ? regionFill.find(region)->second : zeroTree;
					return ConvertSurrogateTreeToExpressionTree(t);
				}

				RangeExpander across(parser);
				across.appendStartExpression(concreteMatrix[row][startCol]);
				across.setEndExpression(concreteMatrix[row][endCol]);

				return across.at(col - startCol);
			}
		} else if (region.type == Region::RECTANGULAR) {
			for (size_t i = 0; i < region.elems.size(); i++) {
				std::pair<size_t,size_t> coords = abstractToConcreteIdx.find(region.elems[i])->second;
				c1 = std::min(c1, coords.second);
				c2 = std::max(c2, coords.second);
				r1 = std::min(r1, coords.first);
				r2 = std::max(r2, coords.first);
			}

			if (row >= r1 && row <= r2 && col >= c1 && col <= c2) {
				// If there is a concrete fill specified, use it
				if (regionFill.find(region) != regionFill.end()) {
					SurrogateTree* t = regionFill.find(region)->second;
					return ConvertSurrogateTreeToExpressionTree(t);
				}

				RangeExpander across(parser);

				for (size_t i = c1; i > 0; i--) {
					const ExpressionTree* t = concreteMatrix[row][i-1];
					if (t == NULL) break;
					across.insertStartExpression(t);
				}
				across.setEndExpression(concreteMatrix[row][c2]);

				while (!across.isValid() && across.getStartExpressions().size() > 0) {
					across.removeStartExpression(0);
				}
				
				return across.at(col - c1);
			}
		}
	}

	return NULL;
}

size_t UnderspecifiedMatrix::getWidth(Graph &g)
{
	MatrixElement* top_left = (*inheritedMatrix)(0,0);
	if (top_left == NULL) return 0;

	std::map<Graph::Edge, bool> empty;
	std::map<Graph::Node, size_t> empty_cache;
	size_t width = getDimension(g, *(g.elemToNode(top_left)), empty, empty_cache, true);
	return width + 1; //+1 to account for top node
}

size_t UnderspecifiedMatrix::getHeight(Graph &g)
{
	MatrixElement* top_left = (*inheritedMatrix)(0,0);
	if (top_left == NULL) return 0;

	std::map<Graph::Edge, bool> empty;
	std::map<Graph::Node, size_t> empty_cache;
	size_t height = getDimension(g, *(g.elemToNode(top_left)), empty, empty_cache, false);
	return height + 1; //+1 to account for top node
}

size_t UnderspecifiedMatrix::getDimension(Graph& g, Graph::Node& node, std::map<Graph::Edge, bool>& visited, std::map<Graph::Node, size_t>& cache, bool isWidth)
{
	/* Start at the top left and recursively analyse each connected node, recursively progressing towards
	 * the bottom right. If each path returns the same dimension, we can conclude that it is a valid
	 * one
	 */

	size_t ret = 0;
	Graph::Node* right = NULL;
	Graph::Node* down = NULL;
	Graph::Node* down_right = NULL;
	Graph::Node* up_right = NULL;
	Graph::Edge::Type noProgress = isWidth ? Graph::Edge::VERTICAL : Graph::Edge::HORIZONTAL;

	size_t returns[4];
	size_t numReturns = 0;

	size_t nRow, nCol;
	inheritedMatrix->getPosition(*node.elem, nRow, nCol);

	std::vector<Graph::Edge> edges = g.getEdges(node);
	for (std::vector<Graph::Edge>::iterator it = edges.begin(); it != edges.end(); it++) {
		Graph::Edge edge = *it;
		Graph::Node& otherNode = edge.otherNode(node);

		//if (visited.find(edge) != visited.end()) continue;

		Graph::Node* other = NULL;
		Rect<long> startBbox = node.elem->getBoundingRect();
		Rect<long> endBbox = otherNode.elem->getBoundingRect();
		switch (edge.type)
		{
		case Graph::Edge::HORIZONTAL:
			if (endBbox.right < startBbox.left) continue; //don't progress backwards
			other = right = &otherNode;
			break;
		case Graph::Edge::VERTICAL:
			if (endBbox.bottom < startBbox.top) continue; //don't progress backwards
			other = down = &otherNode;
			break;
		case Graph::Edge::DIAGONAL:
			if (endBbox.right < startBbox.left) continue; //don't progress backwards
			other = down_right = &otherNode;
			break;
		case Graph::Edge::ANTIDIAGONAL:
			if (endBbox.right < startBbox.left) continue; //don't progress backwards
			other = up_right = &otherNode;
			break;
		};

		visited[edge] = true;

		if (other != NULL) {
			size_t d;
			if (cache.find(*other) != cache.end()) {
				d = cache[*other];
			} else {
				// Recursively add the dimension of the connected element
				d = getDimension(g, *other, visited, cache, isWidth);
			}

			if (d < 0) { cache[node] = d; return d; }

			// If the edge vertical, don't count it towards the width (not going horizontally) (vice versa for height)
			ret = returns[numReturns++] = (edge.type == noProgress) ? d : edge.getWeight() + d + 1;
		}
 	}

	size_t rows, cols;
	inheritedMatrix->getDimensions(rows, cols);
	// If we didn't encounter a vertical ellipsis, use the element below (if exists)
	if (down == NULL && rows > nRow + 1) {
		MatrixElement* below = (*inheritedMatrix)(nRow + 1, nCol);
		if (below != NULL && !isFillerExpr(below)) {
			down = g.elemToNode(below);
			if (down != NULL) {
				size_t d = getDimension(g, *down, visited, cache, isWidth);
				if (d < 0) { cache[node] = d; return d; }
				ret = returns[numReturns++] = isWidth ? d : d + 1;
			}
		}
	}

	// Same for if we didn't find a horizontal ellipsis
	if (right == NULL && cols > nCol + 1) {
		MatrixElement* beside = (*inheritedMatrix)(nRow, nCol + 1);
		if (beside != NULL && !isFillerExpr(beside)) {
			right = g.elemToNode(beside);
			if (right != NULL) {
				size_t d = getDimension(g, *right, visited, cache, isWidth);
				if (d < 0) { cache[node] = d; return d; }
				ret = returns[numReturns++] = !isWidth ? d : d + 1;
			}
		}
	}

	for (size_t i = 1; i < numReturns; i++) {
		if (returns[i] != returns[i-1]) {
			cache[node] = static_cast<size_t>(-1);
			return static_cast<size_t>(-1); //means error
		}
	}

	cache[node] = ret;
	return ret;
}

bool UnderspecifiedMatrix::isEllipsis(MatrixElement* elem)
{
	if (elem == NULL) return false;

	std::map<MatrixElement*, bool>::iterator hIter = hEllipses.find(elem);
	std::map<MatrixElement*, bool>::iterator vIter = vEllipses.find(elem);
	std::map<MatrixElement*, bool>::iterator dIter = dEllipses.find(elem);
	std::map<MatrixElement*, bool>::iterator adIter = adEllipses.find(elem);
	return (hIter != hEllipses.end() || vIter != vEllipses.end() || 
		dIter != dEllipses.end() || adIter != adEllipses.end());
}

bool UnderspecifiedMatrix::isFillerExpr(MatrixElement* elem)
{
	if (elem == NULL) return false;

	return fillExprs.find(elem) != fillExprs.end();
}

const ExpressionTree* UnderspecifiedMatrix::getExprTree(MatrixElement* elem)
{
	if (elem == NULL) return NULL;

	ExpressionIterator* eit = parser->CreateIteratorForStrokes(&(elem->getStrokes()[0]), elem->getStrokes().size(), false);
	const ExpressionTree* tree = eit->next();
	// XXX eit->release();
	return tree;
}


std::vector<MatrixElement*> UnderspecifiedMatrix::getHPrevElems(MatrixElement* ellipsis)
{
	size_t ellipsisRow, ellipsisCol;
	inheritedMatrix->getPosition(*ellipsis, ellipsisRow, ellipsisCol);

	std::vector<MatrixElement*> prevElems;
	try {
		int i = 1;

		MatrixElement* prev = NULL;
		do {
			prev = (*inheritedMatrix)(ellipsisRow, ellipsisCol - i);
			i++;
		} while (prev == NULL || isEllipsis(prev));

		for (; prev != NULL && !isEllipsis(prev); i++) {
			prevElems.insert(prevElems.begin(), prev);
			prev = (*inheritedMatrix)(ellipsisRow, ellipsisCol - i);
		}
	} catch (InvalidBoundsException) { }
	return prevElems;
}

MatrixElement* UnderspecifiedMatrix::getHNextElem(MatrixElement* ellipsis)
{
	size_t ellipsisRow, ellipsisCol;
	inheritedMatrix->getPosition(*ellipsis, ellipsisRow, ellipsisCol);

	MatrixElement* nextElem = NULL;
	try {
		for (int i = 1; nextElem == NULL || isEllipsis(nextElem); i++) {
			nextElem = (*inheritedMatrix)(ellipsisRow, ellipsisCol + i);
		}
	} catch (InvalidBoundsException) { return NULL; }
	return nextElem;
}

std::vector<MatrixElement*> UnderspecifiedMatrix::getVPrevElems(MatrixElement* ellipsis)
{
	size_t ellipsisRow, ellipsisCol;
	inheritedMatrix->getPosition(*ellipsis, ellipsisRow, ellipsisCol);

	std::vector<MatrixElement*> prevElems;
	try {
		int i = 1;

		MatrixElement* prev = NULL;
		do {
			prev = (*inheritedMatrix)(ellipsisRow - i, ellipsisCol);
			i++;
		} while (prev == NULL || isEllipsis(prev));

		for (; prev != NULL && !isEllipsis(prev); i++) {
			prevElems.insert(prevElems.begin(), prev);
			prev = (*inheritedMatrix)(ellipsisRow - i, ellipsisCol);
		}
	} catch (InvalidBoundsException) { }
	return prevElems;
}

MatrixElement* UnderspecifiedMatrix::getVNextElem(MatrixElement* ellipsis)
{
	size_t ellipsisRow, ellipsisCol;
	inheritedMatrix->getPosition(*ellipsis, ellipsisRow, ellipsisCol);

	MatrixElement* nextElem = NULL;
	try {
		for (int i = 1; nextElem == NULL || isEllipsis(nextElem); i++) {
			nextElem = (*inheritedMatrix)(ellipsisRow + i, ellipsisCol);
		}
	} catch (InvalidBoundsException) { return NULL; }
	return nextElem;
}

// TODO diag stuff needs to be more robust e.g. line intersection (or does it? visual representation may lead the user)
std::vector<MatrixElement*> UnderspecifiedMatrix::getDPrevElems(MatrixElement* ellipsis)
{
	size_t ellipsisRow, ellipsisCol;
	inheritedMatrix->getPosition(*ellipsis, ellipsisRow, ellipsisCol);

	std::vector<MatrixElement*> prevElems;
	try {
		int i = 1;

		MatrixElement* prev = NULL;
		do {
			prev = (*inheritedMatrix)(ellipsisRow - i, ellipsisCol - i);
			i++;
		} while (prev == NULL || isEllipsis(prev));

		for (; prev != NULL && !isEllipsis(prev); i++) {
			prevElems.insert(prevElems.begin(), prev);
			prev = (*inheritedMatrix)(ellipsisRow - i, ellipsisCol - i);
		}
	} catch (InvalidBoundsException) { }
	return prevElems;
}

MatrixElement* UnderspecifiedMatrix::getDNextElem(MatrixElement* ellipsis)
{
	size_t ellipsisRow, ellipsisCol;
	inheritedMatrix->getPosition(*ellipsis, ellipsisRow, ellipsisCol);

	MatrixElement* nextElem = NULL;
	try {
		for (int i = 1; nextElem == NULL || isEllipsis(nextElem); i++) {
			nextElem = (*inheritedMatrix)(ellipsisRow + i, ellipsisCol + i);
		}
	} catch (InvalidBoundsException) { return NULL; }
	return nextElem;
}

std::vector<MatrixElement*> UnderspecifiedMatrix::getADPrevElems(MatrixElement* ellipsis)
{
	size_t ellipsisRow, ellipsisCol;
	inheritedMatrix->getPosition(*ellipsis, ellipsisRow, ellipsisCol);

	std::vector<MatrixElement*> prevElems;
	try {
		MatrixElement* prev = (*inheritedMatrix)(ellipsisRow + 1, ellipsisCol - 1);
		for (int i = 2; prev != NULL && !isEllipsis(prev); i++) {
			prevElems.insert(prevElems.begin(), prev);
			prev = (*inheritedMatrix)(ellipsisRow, ellipsisCol - i);
		}
	} catch (InvalidBoundsException) { }
	return prevElems;
}

MatrixElement* UnderspecifiedMatrix::getADNextElem(MatrixElement* ellipsis)
{
	size_t ellipsisRow, ellipsisCol;
	inheritedMatrix->getPosition(*ellipsis, ellipsisRow, ellipsisCol);

	MatrixElement* nextElem = NULL;
	try {
		for (int i = 1; nextElem == NULL || isEllipsis(nextElem); i++) {
			nextElem = (*inheritedMatrix)(ellipsisRow - i, ellipsisCol + i);
		}
	} catch (InvalidBoundsException) { return NULL; }
	return nextElem;
}


}
