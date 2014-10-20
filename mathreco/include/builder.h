#ifndef BUILDER_H_
#define BUILDER_H_


#include "expr-node.h"
#include "grammar.h"
#include "verb.h"

#include <list>
#include <string>
#include <vector>


namespace scg
{

void write_delimited_string(const std::string &s, std::ostream &os);

struct basic_string_builder  : public string_builder {
	explicit basic_string_builder(const std::string &s_) : s(s_) { }
	std::string build(const ExpressionTree *node) const { return s; }
private:
	std::string s;
};

struct ref_string_builder : public string_builder {
	ref_string_builder(size_t ref_, int type_) : ref(ref_), type(type_) { }
	std::string build(const ExpressionTree *node) const;
private:
	size_t ref;
	int type;
};


struct  subref_string_builder : public string_builder {
	subref_string_builder(size_t ref_, size_t subref_, int type_) : ref(ref_), subref(subref_), type(type_) { }
	std::string build(const ExpressionTree *node) const;
private:
	size_t ref;
	size_t subref;
	int type;
};

struct join_string_builder : public string_builder {
	join_string_builder(const std::string &join_, int type_) : join(join_), type(type_) { }
	std::string build(const ExpressionTree *node) const;
private:
	std::string join;
	int type;
};


struct concat_string_builder : public string_builder {
	concat_string_builder(string_builder *first_, string_builder *second_) : first(first_), second(second_) { }
	~concat_string_builder();
	std::string build(const ExpressionTree *node) const;
private:
	string_builder *first;
	string_builder *second;
};



struct terminal_builder : public tree_builder {
	explicit terminal_builder(const std::string &str_, const std::wstring &wstr_, const std::string &mathml_, const std::string &latex_)
		: str(str_), wstr(wstr_), latex(latex_), mathml(mathml_) { }
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
private:
	std::string str;
	std::wstring wstr;
	std::string latex;
	std::string mathml;
};

struct symbol_builder : public tree_builder {
	explicit symbol_builder(const std::string &name);
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
private:
	const symbol *S;
};

struct blank_tree_builder : public tree_builder {
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
};

struct placeholder_tree_builder : public tree_builder {
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
};

struct join_tree_builder : public tree_builder {
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
};

struct pullup_tree_builder : public tree_builder {
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
};

struct ref_tree_builder : public tree_builder {
	explicit ref_tree_builder(size_t ref_) : ref(ref_) { }
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
private:
	size_t ref;
};

struct subref_tree_builder : public tree_builder {
	subref_tree_builder(size_t ref_, size_t subref_) : ref(ref_), subref(subref_) { }
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
private:
	size_t ref;
	size_t subref;
};

struct concat_tree_builder : public tree_builder {
	concat_tree_builder(tree_builder *first_, tree_builder *second_, const production *P_) : first(first_), second(second_), P(P_) { }
	~concat_tree_builder();
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
private:
	tree_builder *first;
	tree_builder *second;
	const production *P;
};


struct subtree_builder : public tree_builder {
	explicit subtree_builder(const std::list<tree_builder *> &child_builders_)
		: child_builders(child_builders_) { }
	~subtree_builder();
	basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
private:
	std::list<tree_builder *> child_builders;
};


}


#endif
