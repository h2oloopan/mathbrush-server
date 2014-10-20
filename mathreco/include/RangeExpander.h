#ifndef SCG_RANGE_EXPANDER_H
#define SCG_RANGE_EXPANDER_H

#include <vector>
#include <istream>
#include <sstream>
#include <string>

#include "MathRecoTypes.h"
#include "functions.h"
#include "surrogate-tree.h"
#include "error.h"

namespace scg
{


// This class expands a range defined by an ellipsis between one or more expressions
class RangeExpander
{
	std::vector<const ExpressionTree *> startExprs;
	const ExpressionTree *endExpr;
	SurrogateTree* templateTree; //template tree of the implicit trees determined by startExpr and endExpr
	MathRecognizer* parser;

public:
	RangeExpander(MathRecognizer* parser = NULL, const ExpressionTree* startExpr = NULL, const ExpressionTree* endExpr = NULL);
	RangeExpander(const RangeExpander& other);
	RangeExpander& operator=(const RangeExpander& other);
	~RangeExpander();

	// Append an expression tree that defines the start of the range
	void appendStartExpression(const ExpressionTree* expr);
	// Insert an expression tree into the set of starting expressions
	void insertStartExpression(const ExpressionTree* expr, size_t index = 0);
	// Set the expression tree that defines the end of the range
	void setEndExpression(const ExpressionTree* expr);

	void removeStartExpression(size_t index);

	std::vector<const ExpressionTree *> getStartExpressions() const { return startExprs; }
	const ExpressionTree * getEndExpression() const { return endExpr; }

	MathRecognizer* getParser();
	void setParser(MathRecognizer* reco);

	// Get the number of elements defined by this range
	size_t size() const;

	// Get the value at a specific index of the range
	ExpressionTree* at(size_t index);

	// Clear the start and end elements
	void clear();

	// Returns if this range is well defined (i.e. has a starting and an ending expression) or not
	bool hasBoundingExpressions() { return startExprs.size() > 0 && endExpr != NULL; }

	bool isValid() { if (templateTree == NULL) constructTemplateTree(); return (templateTree != NULL); }

	/* Update the template tree if the bounds have changed externally. Otherwise, the update will be performed
	 * upon the next RangeExpander::at() function call.
	 */
	void updateTemplateTree();

private:
	/* Constructs a tree that acts as a template for an expression tree in the range. That is, the terminal nodes
	 * that define the range (e.g. for 1...3 it's '1' and '3') will have an interpolating function instead of a
	 * concrete value.
	 */
	void constructTemplateTree();

	/* Construct a surrogate tree from "trees", where the trees comprised in the list are all on the same level
	 * of an expression tree (e.g. all root notes, all terminals, etc.)
	 */
	SurrogateTree* constructSurrogateTree(std::vector<const ExpressionTree *> trees);

	// Determine the function that maps between the terminal nodes
	BoundedUnivariateFunction* determineMappingFunction(std::vector<const ExpressionTree *> terminalNodes);

	//expression_node convert(const ExpressionTree* tree);
};


}


#endif