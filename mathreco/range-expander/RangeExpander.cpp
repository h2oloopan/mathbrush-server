#include "equiv-transforms.h"
#include "MathRecognizer.h"
#include "RangeExpander.h"

#include <stdio.h>
#include <stack>
#include <queue>
#include <limits>

namespace scg
{

RangeExpander::RangeExpander(MathRecognizer* parser, const ExpressionTree* startExpr, const ExpressionTree* endExpr_)
	: parser(parser), templateTree(NULL) {	
	if (startExpr) {
		appendStartExpression(startExpr);
	}
	if (endExpr_) {
		setEndExpression(endExpr_);
	}
	else {
		endExpr = 0;
	}
}

RangeExpander::RangeExpander(const RangeExpander& other)
{
	for (std::vector<const ExpressionTree *>::const_iterator it = other.startExprs.begin(); it != other.startExprs.end(); it++) startExprs.push_back(*it);
	if (other.endExpr != NULL) {
		endExpr = other.endExpr;
	}

	parser = other.parser;
	templateTree = NULL;
	
	updateTemplateTree();
}

RangeExpander& RangeExpander::operator=(const RangeExpander& other)
{
	if (this != &other) {
		clear();

		for (std::vector<const ExpressionTree *>::const_iterator it = other.startExprs.begin(); it != other.startExprs.end(); it++) startExprs.push_back(*it);
		if (other.endExpr != NULL) {
			endExpr = other.endExpr;
		}

		parser = other.parser;
		templateTree = NULL;
		
		updateTemplateTree();
	}

	return *this;
}

RangeExpander::~RangeExpander()
{
	clear();
	/*startExprs.clear();
	if (templateTree != NULL) {
		delete templateTree;
		templateTree = NULL;
	}*/
}

void RangeExpander::clear()
{
	startExprs.clear();
	if (endExpr != NULL) {
		//endExpr->release();
		endExpr = 0;//expression_node();
	}

	if (templateTree != NULL) {
		delete templateTree;
		templateTree = NULL;
	}
}

void RangeExpander::appendStartExpression(const ExpressionTree* expr)
{
	if (expr == NULL) return;

	startExprs.push_back(expr);//convert(expr));
	if (templateTree != NULL) {
		delete templateTree;
		templateTree = NULL;
	}

	if (endExpr != NULL) constructTemplateTree();
}

void RangeExpander::insertStartExpression(const ExpressionTree* expr, size_t index)
{
	if (expr == NULL) return;

	if (index >= startExprs.size()) appendStartExpression(expr);
	else startExprs.insert(startExprs.begin() + index, expr);//convert(expr));

	if (templateTree != NULL) {
		delete templateTree;
		templateTree = NULL;
	}

	if (endExpr != NULL) constructTemplateTree();
}

void RangeExpander::removeStartExpression(size_t index)
{
	if (index >= startExprs.size()) return;

	startExprs.erase(startExprs.begin() + index);
	if (templateTree != NULL) {
		delete templateTree;
		templateTree = NULL;
	}

	if (endExpr != NULL) constructTemplateTree();
}

void RangeExpander::setEndExpression(const ExpressionTree* expr)
{
	if (expr == NULL) return;
	
	//if (endExpr != NULL) endExpr->release();
	endExpr = expr;//convert(expr);
	
	if (templateTree != NULL) {
		delete templateTree;
		templateTree = NULL;
	}

	if (startExprs.size() > 0) constructTemplateTree();
}

void RangeExpander::setParser(MathRecognizer* reco)
{
	this->parser = parser;
}

MathRecognizer* RangeExpander::getParser() 
{
	return parser;
}

size_t RangeExpander::size() const
{
	if (templateTree == NULL) {
		return 0;
	}

	// Traverse tree until we encounter a leaf node that is an unspecified terminal expression node,
	// then calculate the size of its domain
	std::stack<const SurrogateTree*> nodeStack;
	const UnspecifiedTerminalNode* generalTermNode = NULL;
	nodeStack.push(templateTree);

	while (!nodeStack.empty()) {
		const SurrogateTree* treeNode = nodeStack.top();
		nodeStack.pop();

		if (treeNode->getType() == TERMINAL_EXPR) {
			// If this terminal node is an unspecified terminal node, use its domain to calculate the size and return
			if ((generalTermNode = dynamic_cast<const UnspecifiedTerminalNode*>(treeNode)) != NULL) {
				const BoundedUnivariateFunction* f = generalTermNode->getInterpolatingFunction();
				std::pair<int, int> domain = f->getDomain();
				return domain.second - domain.first + 1;
			}
		} else {
			// Push the node's children onto the stack for analysis right-to-left
			for (size_t i = treeNode->nchildren(); i > 0; i--) {
				nodeStack.push(treeNode->child(i-1));
			}
		}
	}

	return 0;
}

ExpressionTree* RangeExpander::at(size_t index)
{
	if (parser == NULL) return NULL;

	if (templateTree == NULL) {
		constructTemplateTree();

		// If it wasn't properly constructed, fail
		if (templateTree == NULL) return NULL;
	}

	SurrogateTree* holder = templateTree->copy();

	const UnspecifiedTerminalNode* termNode = NULL;
	std::stack<const SurrogateTree*> nodeStack;
	nodeStack.push(holder);

	while (!nodeStack.empty()) {
		const SurrogateTree* treeNode = nodeStack.top();
		nodeStack.pop();

		if (treeNode->getType() == TERMINAL_EXPR) {
			// If this terminal node is an unspecified terminal node, use its domain to calculate the size and return
			if ((termNode = dynamic_cast<const UnspecifiedTerminalNode*>(treeNode)) != NULL) {
				(const_cast<UnspecifiedTerminalNode*>(termNode))->evaluate(static_cast<int>(index));
			}
		} else {
			// Push the node's children onto the stack for analysis right-to-left
			for (size_t i = treeNode->nchildren(); i > 0; i--) {
				nodeStack.push(treeNode->child(i-1));
			}
		}
	}

	// Convert the tree structure to an ExpressionTree
	ExpressionTree* tree = ConvertSurrogateTreeToExpressionTree(holder);

	return tree;
}

BoundedUnivariateFunction* RangeExpander::determineMappingFunction(std::vector<const ExpressionTree *> terminalNodes)
{
	BoundedUnivariateFunction* ret = NULL;

	const wchar_t* wstr_value = NULL;
	std::vector<int> values;
	for (std::vector<const ExpressionTree *>::iterator it = terminalNodes.begin(); it != terminalNodes.end(); it++) {
		int value;
		wstr_value = (*it)->wstr();
		std::wstringstream wss;
		wss << wstr_value;
		wss >> value;

		values.push_back(value);
	}

	if (values.size() == 2) {
		// If there are only 2 values, fit to a linear equation
		int m = (values[0] > values[1]) ? -1 : 1; //slope
		int b = static_cast<int>(values[0]); //y-intercept
		ret = new BoundedLinearFunction(m, b, 0, std::numeric_limits<int>::max()); //temporarily int max - so we can evaluate below

		// Sanity check
		double result;
		ret->evaluate(result, 0);
		if (result != values[0]) {
			delete ret;
			VERBOSE(*verb_out << "\tRangeExpander: Generated function does not fit first starting value (0," << values[0] << ")" << std::endl);
			return NULL;
		}
	} else {
		// If there are more than 2 values, fit a polynomial to the points. Also, don't use the last element, as it will
		// be out of sequence (e.g. 1 2 3 ... 10 --> last element is 10); instead, use it for the bounds
		std::vector<std::pair<int,int> > points;
		// Try fitting the values to x values starting at 1 and not including the last value (see comment above)
		for (size_t i = 0; i < values.size() - 1; i++) {
			points.push_back(std::make_pair(static_cast<int>(i), values[i]));
		}

		// TODO apply strategies if not a 'nice' function
		BoundedPolynomialFunction lagrange = createLagrangePolynomial(points);
		if (lagrange.getNumTerms() == 0 || (values[values.size() - 1] > values[0] && lagrange.getDegree() == 0)) return NULL;

		ret = new BoundedPolynomialFunction(lagrange);
		ret->setDomain(0, std::numeric_limits<int>::max()); //temporarily int max - so we can evaluate below

		// Sanity check
		std::pair<int, int> lastXY = points.back();
		double result;
		ret->evaluate(result, lastXY.first);
		if (result != lastXY.second) {
			delete ret;
			VERBOSE(*verb_out << "\tRangeExpander: Generated function does not fit last starting value (" << lastXY.first << "," << lastXY.second << ")" << std::endl);
			return NULL;
		}
	}

	// Solve for f(x) = upperValue to find the upper bound of the domain
	int upperValue = values[values.size() - 1];
	double result = (values[0] < upperValue) ? 0.0 : upperValue + 1;
	for (int x = 0; (values[0] < upperValue && result <= upperValue) || (values[0] > upperValue && result >= upperValue); x++) {
		ret->evaluate(result, x);
		if (result == upperValue) {
			ret->setDomain(0, x);
			break;
		}
	}

	if (result != upperValue) {
		VERBOSE(*verb_out << "\tRangeExpander: Could not determine proper function -- returning null mapping function" << std::endl);
		if (ret != NULL) delete ret;
		return NULL;
	}

	return ret;
}

SurrogateTree* RangeExpander::constructSurrogateTree(std::vector<const ExpressionTree *> trees)
{
	SurrogateNode* surrogateNode = NULL;

	if (parser == NULL) return NULL;

	// Determine if all of the trees are terminal trees (and/or if they are all the same type)
	bool allAreTerminal = true, allAreSameType = true;
	SemanticId prevType = trees[0]->type();
	for (std::vector<const ExpressionTree *>::iterator it = trees.begin(); it != trees.end(); it++) {
		SemanticId treeType = (*it)->type();//it->ptr()->type();
		allAreTerminal = allAreTerminal && (treeType == TERMINAL_EXPR);
		allAreSameType = allAreSameType && (treeType == prevType);
		prevType = treeType;
	}

	if (allAreTerminal) {
		// If they are all terminal trees, figure out if we're dealing with a range of values and which are
		// dealing with a constant value across the range
		bool allAreEqual = true;
		std::string val = trees[0]->long_str();
		for (std::vector<const ExpressionTree *>::iterator it = trees.begin(); it != trees.end(); it++) {
			const char* treeStr = (*it)->long_str();//it->ptr()->long_str();
			allAreEqual = allAreEqual && (treeStr == val);
			val = treeStr;
		}

		if (allAreEqual) {
			// If all the terminals have the same value, no interpolation is necessary
			surrogateNode = new SurrogateTree(TERMINAL_EXPR);
			surrogateNode->setTerminalValue(val); 
		} else {
			// If they do not have the same representation, there must exist a function mapping the values
			VERBOSE(*verb_out << "\tRangeExpander: Found variable terminal expression" << std::endl);
			surrogateNode = new UnspecifiedTerminalNode();

			// Determine the function that maps between the nodes
			BoundedUnivariateFunction* fxn = determineMappingFunction(trees);

			if (fxn == NULL) {
				VERBOSE(*verb_out << "\tRangeExpander: Could not map function" << std::endl);
				return NULL;
			}

			// Set the function to this unspecified terminal node so it can calculate what its value is
			// supposed to be, according to a specified index (see UnspecifiedTerminalNode::evaluate(int index) )
			((UnspecifiedTerminalNode*) surrogateNode)->setInterpolatingFunction(fxn);
		}
	} else if (allAreSameType) {
		// If they are not terminals, but they are the same type, link the similarity by analyzing their children
		surrogateNode = new SurrogateTree(trees[0]->type());

		// TODO assumes same children - valid?
		for (unsigned i = 0; i < trees[0]->nchildren(); i++) {
			std::vector<const ExpressionTree *> ithChildren;
			for (std::vector<const ExpressionTree *>::iterator it = trees.begin(); it != trees.end(); it++) {
				const ExpressionTree* subtree = (*it)->child(i);//it->ptr()->child(i);
				//const raw_expression_node* rawNode = dynamic_cast<const raw_expression_node*>(subtree);
				ithChildren.push_back(subtree);//expression_node(const_cast<raw_expression_node*>(rawNode)));
			}

			// Build a surrogate tree based on the next set of children
			SurrogateNode* surrogateChild = constructSurrogateTree(ithChildren);
			if (surrogateNode != NULL && surrogateChild != NULL) {
				surrogateNode->addChild(surrogateChild);
			} else {
				VERBOSE(*verb_out << "\tRangeExpander: Null surrogate node!" << std::endl);
				return NULL;
			}
		}
	} else {
		// If they are not of the same type, see if there is a rule that transforms one or more trees such that
		// all the trees have the same structure
		EquivalentTransforms& transforms = EquivalentTransforms::getInstance(*GetMathGrammar());
		SurrogateTree* tn = expressionToSurrogateTree(trees.back()/*.ptr()*/);

		for (size_t i = 0; i < trees.size() - 1; i++) {
			SurrogateTree* t = expressionToSurrogateTree(trees[i]/*.ptr()*/);
			if (t->getType() != tn->getType()) {
				SurrogateTree* equiv = transforms.getTransform(t, tn);
				if (equiv != NULL) {
					trees[i] = ConvertSurrogateTreeToExpressionTree(equiv);;
				} else {
					equiv = transforms.getTransform(tn, t);
					if (equiv == NULL) {
						delete t;
						delete tn;
						return NULL;
					}
					trees.back() = ConvertSurrogateTreeToExpressionTree(equiv);
				}

				delete equiv;
			}

			delete t;
		}

		delete tn;
		surrogateNode = constructSurrogateTree(trees);
	}

	return surrogateNode;
}

void RangeExpander::constructTemplateTree()
{
	VERBOSE(*verb_out << "\tRangeExpander: Constructing template tree" << std::endl);

	if (startExprs.size() == 0 || endExpr == NULL) {
		VERBOSE(*verb_out << "\tRangeExpander: Invalid start/end expressions" << std::endl);
		return;
	}

	std::vector<const ExpressionTree *> expressions;
	// If we have been given more than 1 start expression, find a function mapping in the start elements
	// Otherwise, find a linear mapping between the start and end elements
	expressions = startExprs;
	expressions.push_back(endExpr);

	templateTree = constructSurrogateTree(expressions);
	if (templateTree == NULL) {
		// We didn't successfully construct a surrogate tree
		VERBOSE(*verb_out << "\tRangeExpander: NULL surrogate tree (couldn't find valid range)" << std::endl);
	}
}

void RangeExpander::updateTemplateTree()
{
	if (startExprs.size() > 0 && endExpr != NULL) {
		if (templateTree != NULL) {
			delete templateTree;
			templateTree = NULL;
		}

		constructTemplateTree();
	}
}

/*expression_node RangeExpander::convert(const ExpressionTree* tree)
{
	const raw_expression_node* node = dynamic_cast<const raw_expression_node*>(tree);
	return expression_node(const_cast<raw_expression_node*>(node));
}*/


}
