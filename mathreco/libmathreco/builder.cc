#include "builder.h"
#include "relation.h"
#include "expr-node.h"
#include "grammar.h"
#include "symbols.h"

#include <cassert>
#include <set>


namespace scg
{


void
write_delimited_string(const std::string &s, std::ostream &os) {
	for (std::string::const_iterator i = s.begin(); i != s.end(); ++i) {
		if (*i == '$') {
			os << "\\$";
		}
		else {
			os << *i;
		}
	}
	os << '$';
}


std::string
ref_string_builder::build(const ExpressionTree *node) const {
	const ExpressionTree *tree = node->child(ref);
	if (!tree) {
		THROW_ERROR(E_INVALID, "no tree for ref " << ref);
	}
	return tree->getstr(type);
}

std::string
subref_string_builder::build(const ExpressionTree *node) const {
	const ExpressionTree *tree = node->child(ref);
	if (!tree) {
		THROW_ERROR(E_INVALID, "no tree for ref " << ref);
	}
	assert(subref < tree->nchildren());
	return tree->child(subref)->getstr(type);
}

std::string
join_string_builder::build(const ExpressionTree *node) const {
	if (node->nchildren() == 0) return std::string();
	std::stringstream ss;
	ss << node->child(0)->getstr(type);
	for (size_t i = 1; i < node->nchildren(); ++i) {
		ss << join;
		ss << node->child(i)->getstr(type);
	}
	return ss.str();
}

concat_string_builder::~concat_string_builder() {
	delete first;
	delete second;
}

std::string
concat_string_builder::build(const ExpressionTree *node) const {
	return first->build(node) + second->build(node);
}

basic_tree *
blank_tree_builder::build(const interpretation *src, subtree_accessor &acc) const {
	return mkblank();
}

basic_tree *
placeholder_tree_builder::build(const interpretation *src, subtree_accessor &acc) const {
	return mkplaceholder();
}

basic_tree *
terminal_builder::build(const interpretation *src, subtree_accessor &acc) const {
	return new terminal_tree(src, str, wstr, mathml, latex);
}

symbol_builder::symbol_builder(const std::string &name) {
	S = symdb_findsymbol_name(name);
	if (!S) {
		THROW_ERROR(E_NOTFOUND, "symbol_builder could not find symbol named " << name);
	}
}

basic_tree *
symbol_builder::build(const interpretation *src, subtree_accessor &acc) const {
	return new terminal_tree(src, S->name, std::wstring(1, S->unicode), S->mathml, S->latex);
}

basic_tree *
ref_tree_builder::build(const interpretation *src, subtree_accessor &acc) const {
	basic_tree *tree = acc[ref];
	acc.markused(ref);
	return tree;
}

basic_tree *
subref_tree_builder::build(const interpretation *src, subtree_accessor &acc) const {
	basic_tree *tree = acc[ref];
	if (!tree) {
		return 0;
	}
	assert(subref < tree->nchildren());
	tree->chown(subref, false);
	return tree->child(subref);
}

concat_tree_builder::~concat_tree_builder() {
	delete first;
	delete second;
}

basic_tree *
concat_tree_builder::build(const interpretation *src, subtree_accessor &acc) const {
	const basic_tree *left = first->build(src, acc);
	if (!left) {
		return 0;
	}
	const basic_tree *right = second->build(src, acc);
	if (!right) {
		delete left;
		return 0;
	}
	std::string mathml;
	mathml.append(left->long_str());
	mathml.append(right->long_str());
	std::string latex;
	latex.append(left->latex_str());
	latex.append(right->latex_str());
	basic_tree *ret = new terminal_tree(src, left->ccstr() + right->ccstr(), left->ccwstr() + right->ccwstr(), mathml, latex);
	delete left;
	delete right;
	return ret;
}

basic_tree *
join_tree_builder::build(const interpretation *src, subtree_accessor &acc) const {
	basic_tree *tree = new parsed_tree(src);
	for (size_t i = 0; i < src->nchildren(); ++i) {
		basic_tree *child = acc[i];
		acc.markused(i);
		if (!child) {
			delete tree;
			return 0;
		}
		tree->addchild(child);
	}
	return tree;
}

basic_tree *
pullup_tree_builder::build(const interpretation *src, subtree_accessor &acc) const {
	VERBOSE(*verb_out << "pullup_tree_builder::build()\n");
	assert(src->nchildren() == 2);

	basic_tree *child = acc[0];
	VERBOSE(*verb_out << " got child 0 => " << child->ccstr() << std::endl);
	acc.markused(0);
	if (!child) return 0;
	basic_tree *tree = new parsed_tree(src);
	//if (child->nchildren() > 1) {
		for (size_t i = 0; i < child->nchildren(); ++i) {
			basic_tree *pullup_child = child->child(i);
			if (!pullup_child) {
				delete tree;
				return 0;
			}
			VERBOSE(*verb_out << "  adding pullup-child " << pullup_child->ccstr() << std::endl);
			tree->addchild(pullup_child, false);
		}
	//}
	/*else {
		VERBOSE(*verb_out << "  adding pullup-child " << child->ccstr() << std::endl);
		tree->addchild(child, false);
	}*/

	child = acc[1];
	VERBOSE(*verb_out << " got child 1 => " << child->ccstr() << std::endl);
	acc.markused(1);
	if (!child) {
		delete tree;
		return 0;
	}
	tree->addchild(child);

	VERBOSE(*verb_out << "final tree is " << tree->ccstr() << std::endl);

	return tree;
}

subtree_builder::~subtree_builder() {
	for (std::list<tree_builder *>::iterator i = child_builders.begin(); i != child_builders.end(); ++i) {
		delete *i;
	}
}

basic_tree *
subtree_builder::build(const interpretation *src, subtree_accessor &acc) const {
	basic_tree *tree = new parsed_tree(src);
	size_t j = 0;
	for (std::list<tree_builder *>::const_iterator i = child_builders.begin(); i != child_builders.end(); ++i) {
		basic_tree *child = (*i)->build(src, acc);
		if (!child) {
			delete tree;
			return 0;
		}
		tree->addchild(child);
		++j;
	}
	return tree;
}


}
