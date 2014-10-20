#ifndef EXPR_NODE_H_
#define EXPR_NODE_H_


#include "MathRecoTypes.h"
#include "grammar-fwd.h"
#include "relation.h"
#include "verb.h"
#include "error.h"
#include "iohelp.h"
#include "reco-types.h"

#include <vector>
#include <ostream>
#include <string>
#include <cassert>
#include <iostream>

namespace scg
{

class interpreter;
class ordered_segments;
struct group;


struct expression_box : public ExpressionBox {
	inline int left() const { return box.left; }
	inline int top() const { return box.top; }
	inline int right() const { return box.right; }
	inline int bottom() const { return box.bottom; }

	Rect<long> box;


	expression_box() { }

	explicit expression_box(const Rect<long> &box_) : box(box_) { }

	expression_box(long left_, long top_, long right_, long bottom_)
		: box(left_, top_, right_, bottom_)
	{ }

	expression_box &operator=(const Rect<long> &rhs)
	{
		box = rhs;
		return *this;
	}

	bool operator<(const expression_box &rhs) const { return box < rhs.box; }
	bool operator==(const expression_box &rhs) const { return box == rhs.box; }
	bool operator==(const Rect<long> &rhs) const { return box == rhs; }
};

struct strid {
	enum {
		mathml,
		latex
	};
};

class basic_tree;

class interpretation {
public:
	static const char ID;
public:
	interpretation(interpreter *src, const production *P, const ordered_segments *segs);
	interpretation(reader &re, math_recognizer_base *rec);
	~interpretation();

	void clrstrs();

	recoscore score() const;// { return score_.score(); }
	/*double unnorm_score() const { return score_.unnorm_score(); } 
	const score_combiner &score_combo() const { return score_; }
	void set_score(const score_combiner &sc) { score_ = sc; }*/

	const Rect<long> bounds() const;

	size_t nchildren() const { return children.size(); }
	interpretation *child(size_t i) { return children[i]; }
	const interpretation *child(size_t i) const { return children[i]; }
	void addmetachild(interpretation *child, bool rm = false);
	void addchild(interpretation *child, double relsc, bool rm = false);

	math_recognizer_base *ctx() const;

	const char *str() const;
	const wchar_t *wstr() const;

	virtual bool haslongform() const;
	virtual interpreter *longformiter() const;

	const nonterminal *nt() const;

	virtual basic_tree *mktree() const;

public:
	const production *P;
	const ordered_segments *span;
	const interpretation *head;
	const interpretation *tail;
	interpreter *src;
	int rclass;

	group *grp;
	unsigned nterms;
	double relscore;
	double scorebias;
	double hint;

private:
	mutable int haslong;
	//score_combiner score_;
	std::vector<interpretation *> children;
	std::vector<bool> rmflags;

private:
	recoscore unnorm_score() const;

public:
	virtual const std::string &ccstr() const;
	virtual const std::wstring &ccwstr() const;
	virtual void write(writer &wr) const;
	virtual char id() const { return ID; }

private:
	mutable std::string raw;
	mutable std::wstring wraw;

private:
	friend interpretation *mkterminal(interpreter *, const symbol *, group *grp, const ordered_segments *span);
	friend interpretation *mkterminal(const std::string &s);
};

class fixed_interpretation : public interpretation {
public:
	static const char ID;
public:
	fixed_interpretation(const basic_tree *tree_, bool own_);
	fixed_interpretation(reader &re, math_recognizer_base *rec);
	~fixed_interpretation();
	void write(writer &wr) const;
	char id() const { return ID; }
	const std::string &ccstr() const;
	const std::wstring &ccwstr() const;
	basic_tree *mktree() const;

private:
	const basic_tree *tree;
	bool own;
};

class basic_tree : public ExpressionTree {
public:
	basic_tree() { }
	basic_tree(reader &re, math_recognizer_base *rec);
	~basic_tree();

	basic_tree *mkcopy() const;
	void replace_child(size_t i, basic_tree *child);

protected:
	virtual basic_tree *mkcopy_() const = 0;
public:
	const char *str() const;
	const wchar_t *wstr() const;
	const char *long_str() const;
	const char *latex_str() const;

	virtual const std::string &ccstr() const;
	virtual const std::wstring &ccwstr() const;

	interpreter *mkiter(math_recognizer_base *rec) const;

	virtual const production *prod() const = 0;
	virtual const nonterminal *nt() const = 0;
	virtual const ordered_segments *span() const = 0;

	size_t nchildren() const { return children.size(); }
	void addchild(basic_tree *child, bool rm = true) { children.push_back(child); rmflags.push_back(rm); }
	const basic_tree *child(size_t i) const { return children[i]; }
	basic_tree *child(size_t i) { return children[i]; }
	void chown(size_t i, bool rm) { rmflags[i] = rm; }

	void wrap() const;

	virtual void write(writer &wr) const;
	virtual char id() const = 0;
protected:
	virtual const string_builder *sbuilder() const = 0;
	virtual const string_builder *lbuilder() const = 0;

protected:
	mutable std::string raw;
	mutable std::wstring wraw;
	mutable std::string mathml;
	mutable std::string latex;

private:
	std::vector<basic_tree *> children;
	std::vector<bool> rmflags;
};

bool operator==(const basic_tree &lhs, const basic_tree &rhs);
bool operator!=(const basic_tree &lhs, const basic_tree &rhs);

class invented_tree : public basic_tree {
public:
	invented_tree(const production *P_) : basic_tree(), P(P_) { }
	invented_tree(reader &re, math_recognizer_base *rec);

	void write(writer &wr) const;

private:
	basic_tree *mkcopy_() const;

public:
	const production *prod() const { return P; }
	const nonterminal *nt() const;
	SemanticId type() const;

	const ordered_segments *span() const { return 0; }
	const ExpressionBox *box() const { return 0; }
	double score() const { return 0; }
	int lock() const { return E_INVALID; }
	int unlock() const { return E_INVALID; }
	bool is_locked() const { return false; }
	bool HasLongForm() const { return false; }
	ExpressionIterator *CreateLongFormIterator() const { return 0; }

	char id() const { return ID; }

protected:
	const string_builder *sbuilder() const;
	const string_builder *lbuilder() const;

private:
	const production *P;
public:
	static const char ID;
};


class parsed_tree : public basic_tree {
public:
	parsed_tree(const interpretation *intrp__) : basic_tree(), intrp_(intrp__) {	}
	parsed_tree(reader &re, math_recognizer_base *rec);

	void write(writer &wr) const;

private:
	basic_tree *mkcopy_() const;

public:
	const interpretation *intrp() const { return intrp_; }

	const production *prod() const { return intrp_ ? intrp_->P : 0; }
	const nonterminal *nt() const { return intrp_ ? intrp_->nt() : 0; }
	const ordered_segments *span() const { return intrp_ ? intrp_->span : 0; }

	SemanticId type() const;

	const ExpressionBox *box() const;

	double score() const;

	int lock() const;
	int unlock() const;
	bool is_locked() const;

	bool HasLongForm() const;
	ExpressionIterator *CreateLongFormIterator() const;

public:
	const std::string &ccstr() const;

	const string_builder *sbuilder() const;
	const string_builder *lbuilder() const;

	char id() const { return ID; }

protected:
	const interpretation *intrp_;
public:
	static const char ID;
};

class terminal_tree : public parsed_tree {
public:
	terminal_tree(const interpretation *intrp,
	              const std::string &str, const std::wstring &wstr,
	              const std::string &mathml_, const std::string &latex_);
	terminal_tree(reader &re, math_recognizer_base *rec);

	void write(writer &wr) const;

	SemanticId type() const { return TERMINAL_EXPR; }	
	char id() const { return ID; }

private:
	basic_tree *mkcopy_() const;

//protected:
	//const string_builder *sbuilder() const { return 0; }
	//const string_builder *lbuilder() const { return 0; }
public:
	static const char ID;
};

bool isterminal(const interpretation *intrp);
bool derivesterminal(const interpretation *intrp);

const interpretation *getsemanticintrp(const interpretation *intrp);
const interpretation *getlockintrp(const interpretation *intrp);


interpretation *mkterminal(interpreter *src, const symbol *S, group *grp, const ordered_segments *span);
interpretation *mkterminal(const std::string &s);

basic_tree *mkblank();
basic_tree *mkplaceholder();

int initcanonicalsids(grammar *G);

}


#endif
