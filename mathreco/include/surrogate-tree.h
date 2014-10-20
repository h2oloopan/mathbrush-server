#ifndef SCG_SURROGATE_TREE_H
#define SCG_SURROGATE_TREE_H

#include <vector>
#include <string>

#include "MathRecoTypes.h"
#include "functions.h"

namespace scg
{


/* A surrogate tree is a minimum representation of the expression_node class. This class contains the essentials for
 * constructing a minimal expression tree. It can be passed to the parser in order to obtain an actual expression tree
 * representation.
 */

class SurrogateTree {
protected:
	SemanticId type;
	std::vector<SurrogateTree*> children;
	std::string terminalValue; //used only for terminal nodes

public:
	SurrogateTree();
	SurrogateTree(SemanticId type);
	SurrogateTree(const SurrogateTree& other);
	SurrogateTree& operator=(const SurrogateTree& other);
	virtual ~SurrogateTree();

	SemanticId getType() const;
	void setType(SemanticId type);
	size_t nchildren() const;
	SurrogateTree* child(size_t i);
	const SurrogateTree* child(size_t i) const;
	
	void addChild(SurrogateTree* child);
	SurrogateTree* removeChild(size_t i);

	virtual SurrogateTree* copy() const;

	std::string getTerminalValue() const;
	void setTerminalValue(std::string value);
};


SurrogateTree* expressionToSurrogateTree(const ExpressionTree* t);

typedef SurrogateTree SurrogateNode;

/* Defines a terminal expression node that is defined not by a concrete value, but by
 * an interpolating function.
 */
class UnspecifiedTerminalNode : public SurrogateTree
{
	BoundedUnivariateFunction* f;

public:
	UnspecifiedTerminalNode();
	UnspecifiedTerminalNode(std::string value);
	UnspecifiedTerminalNode(const UnspecifiedTerminalNode& other);
	UnspecifiedTerminalNode& operator=(const UnspecifiedTerminalNode& other);
	~UnspecifiedTerminalNode();

	void setInterpolatingFunction(BoundedUnivariateFunction* f);
	const BoundedUnivariateFunction* getInterpolatingFunction() const;

	/* Evaluate the unspecified node's value according to the interpolating function. After
	 * a call to this method, getTerminalValue() will return the appropriate value
	 */
	void evaluate(int x);

	SurrogateTree* copy() const;
};

bool operator<(const SurrogateTree& left, const SurrogateTree& right);


}

#endif