#ifndef SCG_EQUIV_TRANSFORMS_H
#define SCG_EQUIV_TRANSFORMS_H

#include <istream>
#include <vector>
#include <map>
#include <string>

#include "grammar.h"
#include "surrogate-tree.h"

namespace scg
{


/* Defines a terminal expression node that can be anything, but should usually be treated as a filler to
 * match the target expression (e.g. 1 = (var)^0, so 1 ... x^10, the filler node is the terminal node for
 * the variable, and should be filled with "x" in this case)
 */
class FillerTerminalNode : public SurrogateTree
{
public:
	FillerTerminalNode();
	~FillerTerminalNode();

	SurrogateTree* copy() const;
};

/* Defines a tree that represents a copy of the source tree */
class CopyNode : public SurrogateTree
{
public:
	CopyNode();
	~CopyNode();

	SurrogateTree* copy() const;
};


/* This singleton class reads in a grammar defined just like portable_tree grammar files, and generates
 * production rules that map a certain tree structure to other equivalent tree structures
 *
 * For example, "x" (with the tree structure VAR(TERM x) ) is equivalent to "x^1" (with the
 * tree structure SUP(VAR(TERM x), NUM(TERM 1))
 */
class EquivalentTransforms {
	static const char COPY_TREE_CHAR;

	// Compares SurrogateTrees without considering the leaf (terminal) nodes
	class NoLeafCompare
	{
	public:
		bool operator()(const SurrogateTree& left, const SurrogateTree& right) const;
	};

	// Maps a SurrogateTree form to its equivalent transformations
	std::map<SurrogateTree, std::vector<SurrogateTree*>, NoLeafCompare > rules;

public:
	// singleton
	static EquivalentTransforms& getInstance(const grammar& G);

	// Get all the equivalent transformations for "src"
	std::vector<SurrogateTree*> getTransforms(const SurrogateTree* src);
	// Get the transformation from src to target
	SurrogateTree* getTransform(const SurrogateTree* src, const SurrogateTree* target);

private:
	EquivalentTransforms(const grammar& G);
	EquivalentTransforms(const EquivalentTransforms& other);
	EquivalentTransforms& operator=(const EquivalentTransforms& other);
	~EquivalentTransforms();

	// Read the transformation-trees transformations file to fill "rules"
	void readGrammar(const grammar& G, std::istream& is);
	std::string readTree(std::istream& is);
	SurrogateTree* parseTreeStructure(const grammar& G, const std::string& treeStructure);

	/* Helper function for getTransform(src, target). It needs "original" so that it can replace any CopyNode instances with a copy
	 * of the original
	 */
	SurrogateTree* getTransform(const SurrogateTree* src, const SurrogateTree* target, const SurrogateTree* original);

	// Return a list of trees that correspond to the transformations of "src" node to match "target" node
	std::vector<SurrogateTree*> matchRootNode(const SurrogateTree* src, const SurrogateTree* target);
};


}

#endif