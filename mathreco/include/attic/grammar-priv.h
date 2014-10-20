#ifndef GRAMMAR_PRIV_H_
#define GRAMMAR_PRIV_H_


#include "grammar.h"
#include "verb.h"

#include <list>
#include <string>
#include <vector>


namespace scg
{
namespace parser
{

typedef std::map<BoxPair, std::vector<double> > BoxScores;


struct ParseContextPrivate
{
    BoxScores box_scores;

	unsigned **S;       // map stroke index in input -> index in dim d

	const BboxCandidate ***B; // map index of box in dim d -> box in input 
	unsigned **R;             // map index of box in input -> index of box in dim d
	double *Tbase;            // terminal scores at each input box
	ParseCell **D;            // parse table

	const BboxCandidate ****Baux; // map index of box in a span of dim d sorted by d'-> box in input
	unsigned ***Raux;             // reverse map of Baux
	unsigned ***bn;               // min index in d' of boxes in each span_d/d'
	unsigned ***bx;               // max ...
	unsigned ***pbn;              // min index in d of boxes in each span_d/d'
	unsigned ***pbx;              // max ...

	ParseContextPrivate() : B(0), R(0), Tbase(0), D(0), Baux(0), Raux(0), bn(0), bx(0), pbn(0), pbx(0) { }
};



class BasicNameBuilder : public NameBuilder
{
public:
	explicit BasicNameBuilder(const SemanticId &term_) : term(term_) { }
	SemanticId build(ExpressionNode *childA, ExpressionNode *childB) const
		{ VERBOSE(*verb_out << "name string " << term << std::endl); return term; }
	std::ostream &print(std::ostream &os) const
		{ return os << "(STR " << term << ")"; }

private:
	SemanticId term;
};

class RefNameBuilder : public NameBuilder
{
public:
	explicit RefNameBuilder(unsigned ref_) : ref(ref_) { }
	SemanticId build(ExpressionNode *childA, ExpressionNode *childB) const
	{
	    VERBOSE(*verb_out << "name ref " << ref << std::endl);
	    if (!childA) {
	        VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
	        throw E_INVALID;
	    }
	    ExpressionNode *child = childA->child_node(ref-1);
	    if (!child) {
	        VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
	        throw E_INVALID;
	    }
	    return child->type(); }
	std::ostream &print(std::ostream &os) const
		{ return os << "(REF " << ref << ")"; }

private:
	unsigned ref;
};

class SelectorNameBuilder : public NameBuilder
{
public:
	SelectorNameBuilder(unsigned ref_) : ref(ref_) { }
	SemanticId build(ExpressionNode *childA, ExpressionNode *childB) const
	{
		VERBOSE(*verb_out << "name selector " << ref << std::endl);
		ExpressionNode *child = (ref == 1) ? childA : childB;
		if (!child) {
	        VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
	        throw E_INVALID;		
		}
		return child->type();
	}
	std::ostream &print(std::ostream &os) const
		{ return os << "(SEL " << ref << ")"; }

private:
	unsigned ref;
};


class SubSelectorNameBuilder : public NameBuilder
{
public:
	SubSelectorNameBuilder(unsigned ref_, NameBuilder *builder_) : ref(ref_), builder(builder_) { }
	SemanticId build(ExpressionNode *childA, ExpressionNode *childB) const
	{
		VERBOSE(*verb_out << "name subselect " << ref << std::endl);
		ExpressionNode *child = (ref == 1) ? childA : childB;
		if (!child) {
	        VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
	        throw E_INVALID;
		}
		return builder->build(child, 0);
	}
	std::ostream &print(std::ostream &os) const
		{ return os << "(SUBSEL " << ref << *builder << ")"; }

private:
	unsigned ref;
	NameBuilder *builder;
};

class SubRefNameBuilder : public NameBuilder
{
public:
	explicit SubRefNameBuilder(unsigned ref_, NameBuilder *builder_) : ref(ref_), builder(builder_) { }
	SemanticId build(ExpressionNode *childA, ExpressionNode *childB) const
	{
		VERBOSE(*verb_out << "name subref " << ref << std::endl);
		if (!childA) {
	        VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
	        throw E_INVALID;
		}
        ExpressionNode *child = childA->child_node(ref-1);
        if (!child) {
	        VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
	        throw E_INVALID;
        }
		return builder->build(child, 0);
	}
	std::ostream &print(std::ostream &os) const
		{ return os << "(SUB " << ref << " " << *builder << ")"; }

private:
	unsigned ref;
	NameBuilder *builder;
};

/*
class ConcatNameBuilder : public NameBuilder
{
public:
	ConcatNameBuilder(NameBuilder *first_, NameBuilder *second_) : first(first_), second(second_) { }
	std::string build(ExpressionNode *childA, ExpressionNode *childB) const
		{ *verb_out << "name cat " << std::endl; return first->build(childA, childB).append(second->build(childA, childB)); }
	std::ostream &print(std::ostream &os) const
		{ return os << "(CAT " << *first << " " << *second << ")"; }

private:
	NameBuilder *first;
	NameBuilder *second;
};
*/



class TerminalBuilder : public TreeBuilder
{
public:
	explicit TerminalBuilder(const std::string &term_) : term(term_) { }
	MutableExpressionNode *build(MutableExpressionNode *childA, MutableExpressionNode *childB) const
		{ VERBOSE(*verb_out << "building term " << term << std::endl); return make_terminal(term); }
	std::ostream &print(std::ostream &os) const
		{ return os << "(TERM " << term << ")"; }

private:
	std::string term;
};

class NullTreeBuilder : public TreeBuilder
{
public:
	MutableExpressionNode *build(MutableExpressionNode *, MutableExpressionNode *) const { VERBOSE(*verb_out << "null..." << std::endl); return new Expression; }
	std::ostream &print(std::ostream &os) const { return os << "NULL"; }
};


class RefTreeBuilder : public TreeBuilder
{
public:
	explicit RefTreeBuilder(unsigned ref_, bool subchild_ = false) : ref(ref_), subchild(subchild_) { }
	MutableExpressionNode *build(MutableExpressionNode *childA, MutableExpressionNode *childB) const
	{
	    VERBOSE(*verb_out << "building ref @ " << ref << std::endl);
	    MutableExpressionNode *child = childA->mchild(ref-1);
	    if (!child) {
		    VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
	        throw E_INVALID;
	    }
		if (subchild) {
		    MutableExpressionNode *node = new Expression;
		    node->add_child(child);
		    return node;
		}
		else {
		    return child;
		}
	}
	
	std::ostream &print(std::ostream &os) const
		{ return os << "(REF " << ref << " <" << subchild << ">)"; }

private:
	unsigned ref;
	bool subchild;
};


class SelectorTreeBuilder : public TreeBuilder
{
public:
	SelectorTreeBuilder(unsigned ref_, bool subchild_ = false) : ref(ref_), subchild(subchild_) { }
	MutableExpressionNode *build(MutableExpressionNode *childA, MutableExpressionNode *childB) const
	{
		VERBOSE(*verb_out << "selecting @ " << ref << std::endl);
		MutableExpressionNode *child = (ref == 1) ? childA : childB;
		if (!child) {
		    VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
		    throw E_INVALID;
		}
		if (subchild) {
		    MutableExpressionNode *node = new Expression;
		    node->add_child(child);
		    return node;
		}
		else {
		    return child;
		}
	}
	std::ostream &print(std::ostream &os) const
		{ return os << "(SEL " << ref << " <" << subchild << ">)"; }

private:
	unsigned ref;
	bool subchild;
};

class SubSelectorTreeBuilder : public TreeBuilder
{
public:
	SubSelectorTreeBuilder(unsigned ref_, TreeBuilder *builder_) : ref(ref_), builder(builder_) { }
	MutableExpressionNode *build(MutableExpressionNode *childA, MutableExpressionNode *childB) const
	{
		VERBOSE(*verb_out << "subselecting @ " << ref << std::endl);
		MutableExpressionNode *child = (ref == 1) ? childA : childB;
		if (!child) {
		    VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
		    throw E_INVALID;
		}
		return builder->build(child, 0);
	}
	std::ostream &print(std::ostream &os) const
		{ return os << "(SUBSEL " << ref << *builder << ")"; }

private:
	unsigned ref;
	TreeBuilder *builder;
};


class SubRefTreeBuilder : public TreeBuilder
{
public:
	SubRefTreeBuilder(unsigned ref_, TreeBuilder *builder_) : ref(ref_), builder(builder_) { }
	MutableExpressionNode *build(MutableExpressionNode *childA, MutableExpressionNode *childB) const
	{
		VERBOSE(*verb_out << "building sub @ " << ref << std::endl);
		MutableExpressionNode *child = childA->mchild(ref-1);
		if (!child) {
		    VERBOSE(*verb_out << "missing children! check the grammar!" << std::endl);
		    throw E_INVALID;
		}
		return builder->build(child, 0);
	}
	std::ostream &print(std::ostream &os) const
		{ return os << "(SUB " << ref << *builder << ")"; }

private:
	unsigned ref;
	TreeBuilder *builder;
};

class ConcatTreeBuilder : public TreeBuilder
{
public:
	ConcatTreeBuilder(TreeBuilder *first_, TreeBuilder *second_, bool subchild_ = false) : first(first_), second(second_), subchild(subchild_) { }
	MutableExpressionNode *build(MutableExpressionNode *childA, MutableExpressionNode *childB) const
	{
	    VERBOSE(*verb_out << "concat..." << std::endl);
	    MutableExpressionNode *left = first->build(childA, childB);
	    MutableExpressionNode *right = second->build(childA, childB);
	    std::string left_str(left->undecorated_str());
	    std::wstring left_wstr(left->undecorated_wstr());
	    MutableExpressionNode *ret = make_terminal(left_str.append(right->undecorated_str()), left_wstr.append(right->undecorated_wstr()));
	    
	    if (subchild) {
			MutableExpressionNode *node = new Expression;
			node->add_child(ret);
			ret = node;
	    }
	    
	    return ret;
	}
	std::ostream &print(std::ostream &os) const
		{ return os << "(CAT " << *first << " " << *second << ")"; }

private:
	TreeBuilder *first;
	TreeBuilder *second;
	bool subchild;
};


class SubTreeBuilder : public TreeBuilder
{
public:
	explicit SubTreeBuilder(const std::list<TreeBuilder *> &child_builders_) : child_builders(child_builders_) { }

	MutableExpressionNode *build(MutableExpressionNode *childA, MutableExpressionNode *childB) const
	{
		MutableExpressionNode *node = new Expression;
		VERBOSE(*verb_out << "building tree..." << std::endl);
		unsigned n = 0;
		for (std::list<TreeBuilder *>::const_iterator i = child_builders.begin(); i != child_builders.end(); ++i) {
			VERBOSE(*verb_out << " child " << n++ << std::endl);
			node->add_child((*i)->build(childA, childB));
		}
		VERBOSE(*verb_out << "done" << std::endl);
		return node;
	}
	std::ostream &print(std::ostream &os) const
	{
		os << "(TREE ";
		for (std::list<TreeBuilder *>::const_iterator i = child_builders.begin(); i != child_builders.end(); ++i) {
			(*i)->print(os);
			os << ' ';
		}
		return os << ")";
	}

private:
	std::list<TreeBuilder *> child_builders;
};


}
}


#endif
