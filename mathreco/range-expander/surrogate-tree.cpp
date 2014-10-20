#include "surrogate-tree.h"
#include <queue>

namespace scg
{


SurrogateTree::SurrogateTree() {}

SurrogateTree::SurrogateTree(SemanticId type) : type(type) {}

SurrogateTree::SurrogateTree(const SurrogateTree& other)
{
	this->type = other.type;
	for (std::vector<SurrogateTree*>::const_iterator it = other.children.begin(); it != other.children.end(); it++) {
		SurrogateTree* t = *it;
		this->children.push_back(t->copy());
	}
	this->terminalValue = other.terminalValue;
}

SurrogateTree& SurrogateTree::operator=(const SurrogateTree& other)
{
	if (this != &other) {
		this->type = other.type;
		for (std::vector<SurrogateTree*>::const_iterator it = other.children.begin(); it != other.children.end(); it++) {
			this->children.push_back(new SurrogateTree(**it));
		}
		this->terminalValue = other.terminalValue;
	}

	return *this;
}

SurrogateTree* SurrogateTree::copy() const
{
	SurrogateTree* copy = new SurrogateTree(*this);
	return copy;
}

SurrogateTree::~SurrogateTree()
{
	while (children.size() > 0) {
		SurrogateTree* t = children.back();
		children.pop_back();
		delete t;
	} 
}

SemanticId SurrogateTree::getType() const { return type; }

std::string SurrogateTree::getTerminalValue() const { return terminalValue; }

void SurrogateTree::setTerminalValue(std::string value) { terminalValue = value; }

void SurrogateTree::setType(SemanticId type) { this->type = type; }

void SurrogateTree::addChild(SurrogateTree* child) { children.push_back(child); }

SurrogateTree* SurrogateTree::removeChild(size_t i)
{
	if (i >= nchildren()) return NULL;
	
	SurrogateTree* ret = children[i];
	children.erase(children.begin() + i);
	return ret;
}

size_t SurrogateTree::nchildren() const { return children.size(); }

SurrogateTree* SurrogateTree::child(size_t i) { return children[i]; }
const SurrogateTree* SurrogateTree::child(size_t i) const { return children[i]; }


UnspecifiedTerminalNode::UnspecifiedTerminalNode() : SurrogateTree(TERMINAL_EXPR), f(NULL) {}

UnspecifiedTerminalNode::UnspecifiedTerminalNode(std::string value) : SurrogateTree(TERMINAL_EXPR), f(NULL) { terminalValue = value; }

UnspecifiedTerminalNode::UnspecifiedTerminalNode(const UnspecifiedTerminalNode& other)
{
	this->type = TERMINAL_EXPR;
	this->f = (other.f == NULL) ? NULL : other.f->copy();
}

UnspecifiedTerminalNode& UnspecifiedTerminalNode::operator=(const UnspecifiedTerminalNode& other)
{
	if (this != &other) {
		this->type = TERMINAL_EXPR;
		if (f != NULL) delete f;
		this->f = (other.f == NULL) ? NULL : other.f->copy();
	}
	return *this;
}

SurrogateTree* UnspecifiedTerminalNode::copy() const
{
	SurrogateTree* copy = new UnspecifiedTerminalNode(*this);
	return copy;
}

UnspecifiedTerminalNode::~UnspecifiedTerminalNode() { if (f != NULL) delete f; }

void UnspecifiedTerminalNode::setInterpolatingFunction(BoundedUnivariateFunction* f)
{
	if (this->f != NULL) delete this->f;
	this->f = f;
}

const BoundedUnivariateFunction* UnspecifiedTerminalNode::getInterpolatingFunction() const { return f; }

void UnspecifiedTerminalNode::evaluate(int x)
{
	if (f == NULL) return;

	std::stringstream ss;
	double result;
	if (f->evaluate(result, x) > 0) {
		VERBOSE(*verb_out << "\tUnspecifiedTerminalNode: interpolation out of range for index " << x << 
			" [" <<  f->getDomain().first << "," << f->getDomain().second << "]" << std::endl);
		return;
	}
	ss << result;
	terminalValue = ss.str();

	VERBOSE(*verb_out << "\tUnspecifiedTerminalNode: interpolated for index " << x << " and got " << result << std::endl);
}

SurrogateTree* expressionToSurrogateTree(const ExpressionTree* t)
{
	SurrogateTree* ret = new SurrogateTree(t->type());
	
	if (t->type() == TERMINAL_EXPR) {
		std::string value = t->long_str();
		ret->setTerminalValue(value);
	}


	for (size_t i = 0; i < t->nchildren(); i++) {
		const ExpressionTree* child = t->child(i);
		ret->addChild(expressionToSurrogateTree(child));
	}

	return ret;
}


bool operator<(const SurrogateTree& left, const SurrogateTree& right)
{
	int leftScore = 0;
	int rightScore = 0;
	std::queue<const SurrogateTree*> nodes;

	nodes.push(&left);
	while (!nodes.empty()) {
		const SurrogateTree* t = nodes.front();
		nodes.pop();
		leftScore += (int) t->getType();
		for (size_t i = 0; i < t->nchildren(); i++) {
			nodes.push(t->child(i));
		}
	}

	nodes.push(&right);
	while (!nodes.empty()) {
		const SurrogateTree* t = nodes.front();
		nodes.pop();
		rightScore += (int) t->getType();
		for (size_t i = 0; i < t->nchildren(); i++) {
			nodes.push(t->child(i));
		}
	}

	return leftScore < rightScore;
}


}