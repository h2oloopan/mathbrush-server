#include "equiv-transforms.h"

#include <string>
#include <fstream>
#include <cctype>
#include <queue>
#include "error.h"
#include "stream-defs.h"
#include "verb.h"
#include "grammar-values.h"
#include "utils.h"

namespace scg
{

const char EquivalentTransforms::COPY_TREE_CHAR = '%';

bool EquivalentTransforms::NoLeafCompare::operator()(const SurrogateTree& left, const SurrogateTree& right) const
{
	int leftScore = 0;
	int rightScore = 0;
	std::queue<const SurrogateTree*> nodes;

	// Add up the SemanticId values of the tree to determine its score
	nodes.push(&left);
	while (!nodes.empty()) {
		const SurrogateTree* t = nodes.front();
		nodes.pop();
		leftScore += (int) t->getType();
		for (size_t i = 0; i < t->nchildren(); i++) {
			const SurrogateTree* child = t->child(i);
			if (child->getType() != TERMINAL_EXPR) {
				nodes.push(child);
			}
		}
	}

	nodes.push(&right);
	while (!nodes.empty()) {
		const SurrogateTree* t = nodes.front();
		nodes.pop();
		rightScore += (int) t->getType();
		for (size_t i = 0; i < t->nchildren(); i++) {
			const SurrogateTree* child = t->child(i);
			if (child->getType() != TERMINAL_EXPR) {
				nodes.push(child);
			}
		}
	}

	return leftScore < rightScore;
}

FillerTerminalNode::FillerTerminalNode() : SurrogateTree(TERMINAL_EXPR) {}
FillerTerminalNode::~FillerTerminalNode() {}
SurrogateTree* FillerTerminalNode::copy() const
{
	return new FillerTerminalNode(*this);
}

CopyNode::CopyNode() : SurrogateTree() {}
CopyNode::~CopyNode() {}
SurrogateTree* CopyNode::copy() const
{
	return new CopyNode();
}

EquivalentTransforms& EquivalentTransforms::getInstance(const grammar& G)
{
	static EquivalentTransforms inst(G);
	return inst;
}

EquivalentTransforms::EquivalentTransforms(const grammar& G)
{
	std::string path;
	GetProfilePath(path);
	path = path + "tree-transforms";

	if (rules.empty()) {
		std::ifstream grammarFile(path.c_str());
		EquivalentTransforms::readGrammar(G, grammarFile);
		grammarFile.close();
	}
}

EquivalentTransforms::~EquivalentTransforms()
{
	while (rules.size() > 0) {
		std::map<SurrogateTree, std::vector<SurrogateTree*>, NoLeafCompare >::iterator entry = rules.begin();
		if (entry == rules.end()) break;

		std::vector<SurrogateTree*> alternates = entry->second;
		while (alternates.size() > 0) {
			SurrogateTree* t = alternates.back();
			alternates.pop_back();
			delete t;
		}

		rules.erase(entry);
	}
}

std::vector<SurrogateTree*> EquivalentTransforms::getTransforms(const SurrogateTree* src)
{
	return (rules.find(*src) != rules.end()) ? rules[*src] : std::vector<SurrogateTree*>();
}

SurrogateTree* EquivalentTransforms::getTransform(const SurrogateTree* src, const SurrogateTree* target)
{
	return getTransform(src, target, src);
}

SurrogateTree* EquivalentTransforms::getTransform(const SurrogateTree* src, const SurrogateTree* target, const SurrogateTree* original)
{
	const SurrogateTree* ref = (dynamic_cast<const CopyNode*>(src) != NULL) ? original : src;

	// If this node doesn't match the target, try to transform into an equivalent tree that has the same root type
	if (ref->getType() != target->getType()) {
		// Find equivalent trees whose root nodes match the target's root node
		std::vector<SurrogateTree*> txfms = matchRootNode(ref, target);
		for (std::vector<SurrogateTree*>::iterator it = txfms.begin(); it != txfms.end(); it++) {
			SurrogateTree* candidate = *it;
			SurrogateTree* ret = new SurrogateTree(candidate->getType());

			// Ensure subtrees match
			bool hasValidSubtrees = true;
			for (size_t i = 0; i < candidate->nchildren(); i++) {
				SurrogateTree* subtree = getTransform(candidate->child(i), target->child(i), original);
				if (subtree == NULL) {
					hasValidSubtrees = false;
					delete ret;
					break;
				}

				ret->addChild(subtree);
			}

			if (hasValidSubtrees) {
				return ret;
			}
		}
	} else {
		if (target->getType() == TERMINAL_EXPR) {
			// If they are terminal expressions, check to see if we need to fill a filler terminal
			return (dynamic_cast<const FillerTerminalNode*>(ref) != NULL) ? target->copy() : ref->copy();
		} else {
			// Ensure subtrees match
			bool hasValidSubtrees = true;
			SurrogateTree* ret = new SurrogateTree(ref->getType());
			for (size_t i = 0; i < ref->nchildren(); i++) {
				SurrogateTree* subtree = getTransform(ref->child(i), target->child(i), original);
				if (subtree == NULL) {
					hasValidSubtrees = false;
					delete ret;
					break;
				}

				ret->addChild(subtree);
			}
		
			if (hasValidSubtrees) {
				return ret;
			}
		}
	}

	return NULL;
}

std::vector<SurrogateTree*> EquivalentTransforms::matchRootNode(const SurrogateTree* src, const SurrogateTree* target)
{
	std::vector<SurrogateTree*> ret;
	if (src->getType() == target->getType()) {
		ret.push_back(src->copy());
		return ret;
	}

	std::vector<SurrogateTree*> transforms = getTransforms(src);
	for (std::vector<SurrogateTree*>::iterator it = transforms.begin(); it != transforms.end(); it++) {
		std::vector<SurrogateTree*> matches = matchRootNode(*it, target);
		ret.insert(ret.end(), matches.begin(), matches.end());
	}

	return ret;
}

void EquivalentTransforms::readGrammar(const grammar& G, std::istream& is)
{
	// Construct the trees in the production
	while (is && !is.eof() && is.peek() > -1) {
		std::string base = readTree(is);

		SurrogateTree* baseTree = parseTreeStructure(G, base);
		
		// Trees are split by '::'
		int numColons = 0;
		while (!is.eof() && numColons < 2) {
			int c = is.get();
			if (c == ':') numColons++;
		}

		CHECK_ISTREAM_BASIC(is);

		std::string transformation = readTree(is);
		SurrogateTree* transformationTree = parseTreeStructure(G, transformation);

		rules[*baseTree].push_back(transformationTree);

		VERBOSE(*verb_out << "\tEquivalentTransforms: Creating base tree " << baseTree << ". Alternate: " << transformationTree << std::endl);
	}
}

std::string EquivalentTransforms::readTree(std::istream& is)
{
	std::string treeStr;
	unsigned unclosedParens = 0;
	int c;

	// Ignore any prepending spaces
	while (std::isspace(is.peek())) {
		is.get();
		CHECK_ISTREAM_BASIC(is);
	}

	// Must start with '('
	c = is.get();
	CHECK_ISTREAM_BASIC(is);
	if (c != '(') {
		THROW_ERROR(E_INVALID, "surrogate tree cannot start with \'" << c << '\'');
	}

	unclosedParens = 1;

	// Get the whole tree string
	while (!is.eof() && unclosedParens > 0) {
		c = is.get();

		switch(c) {
			case ')':
				if (--unclosedParens > 0) treeStr.append(1, c);
				break;
			case '(':
				unclosedParens++;
				treeStr.append(1, c);
				break;
			default:
				treeStr.append(1, c);
				break;
		}
	}

	return treeStr;
}

SurrogateTree* EquivalentTransforms::parseTreeStructure(const grammar& G, const std::string& treeStructure)
{
	SurrogateTree* node = NULL;

	// If this subtree represents a copy of the original source tree, return a placeholder "CopyNode" to say so
	if (treeStructure[0] == COPY_TREE_CHAR) {
		node = new CopyNode();
		return node;
	}

	// Get the indices of the type of tree node and the beginning of its children
	unsigned i = 0, typeIndex = 0;
	for (i = 0; i < treeStructure.size() && treeStructure[i] != '('; i++) {
		if (std::isspace(treeStructure[i])) {
			typeIndex = i;
		}
	}

	if (typeIndex == 0) typeIndex = i;

	std::string type = treeStructure.substr(0, typeIndex);

	// If it has children (i.e. is non-terminal), parse the children as well
	if (type == "TERM") {
		std::string value;
		char c;
		for (size_t i = typeIndex; i < treeStructure.size(); i++) {
			c = treeStructure[i];

			switch (c) {
			case ' ':  case '\t':
			case '\r': case '\n':
				break;
			case ')':
				goto finished;
			default:
				value.push_back(c);
				break;
			}
		}
finished:

		if (value.size() > 0) {
			node = new SurrogateTree(TERMINAL_EXPR);
			node->setTerminalValue(treeStructure.substr(typeIndex + 1));
		} else {
			node = new FillerTerminalNode();
		}
	} else {
		node = new SurrogateTree(string_to_sid(type));

		unsigned unclosedParens = 0;
		std::vector<std::string> childrenStr;
		int startOfChild = -1;

		// Process everything inside the parentheses (i.e. the children)
		for (; i < treeStructure.size(); i++) {
			char c = treeStructure[i];
			switch(c) {
				case '(':
					unclosedParens++;
					break;
				case ')':
					unclosedParens--;
					if (unclosedParens == 0) {
						childrenStr.push_back(treeStructure.substr(startOfChild, i - startOfChild));
						startOfChild = -1;
					}
					break;
				case ' ': case '\t': case '\r': case '\n':
					break;
				default:
					if (startOfChild == -1) startOfChild = i;
					break;
			}
		}

		// Parse the child trees and add them to this node
		for (std::vector<std::string>::iterator it = childrenStr.begin(); it != childrenStr.end(); it++) {
			node->addChild(parseTreeStructure(G, *it));
		}
	}

	return node;
}



}
