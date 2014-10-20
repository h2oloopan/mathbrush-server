#ifndef GRAMMAR_H_
#define GRAMMAR_H_


#include "ordered-segments.h"
#include "grammar-fwd.h"

#ifndef NO_RECO_TYPES
#include "MathRecoTypes.h" // hook into external interface
#include "grammar-values.h"
#endif

#include <ostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>



namespace scg
{

struct relation;
class math_recognizer_base;

class subtree_accessor {
public:
	virtual ~subtree_accessor() { }
	virtual basic_tree *operator[](size_t i) = 0;
	virtual void markused(size_t i) = 0;
};

struct string_builder {
	virtual ~string_builder() { }
	virtual std::string build(const ExpressionTree *node) const = 0;
};

struct tree_builder {
	virtual ~tree_builder() { }
	virtual basic_tree *build(const interpretation *src, subtree_accessor &acc) const = 0;
};

struct attach_mode {
	enum {
		DEFAULT,
		GROUP = 0,
		SYMBOL,
		SYMBOL_HEAD,
		SYMBOL_TAIL,
		UNSPECIFIED
	};
	
	int from;
	int to;
	
	attach_mode() : from(UNSPECIFIED), to(UNSPECIFIED) {}
	attach_mode(int from_, int to_) : from(from_), to(to_) {}
};

struct production {
	std::vector<nonterminal *> rhs;
	int attach_in_parent;
	std::vector<attach_mode> attach_modes;
	
	std::vector<size_t> terminal_indices;
	
	std::vector<size_t> min_prefix_length;
	std::vector<size_t> max_prefix_length;
	
	size_t index;
	bool isop;

	size_t dim;
	relation *rel;
	size_t head, tail;
	double bias;

	string_builder *sbuild;
	string_builder *lbuild;
	tree_builder *tbuild;

	nonterminal *nt;
	SemanticId sid;

	int rclass;
	
	production(nonterminal *nt_, size_t index_) : dim(0), rel(0), head(0), tail(static_cast<size_t>(-1)), sbuild(0), lbuild(0), tbuild(0), nt(nt_), index(index_), sid(InvalidSemanticId), attach_in_parent(attach_mode::UNSPECIFIED), rclass(0), bias(1), isop(false) { }
	~production();
};

struct uniqproduction : public production {
	explicit uniqproduction(int id);
};

void clearuniqproductions();
const production *getuniqproduction(int id);

struct nonterminal {
	size_t index;
	std::vector<production *> productions;

	std::set<const nonterminal *> equivalent_nts;
	std::set<const nonterminal *> derives_;
	std::set<const nonterminal *> derived_by_;

	size_t min_length;
	size_t max_length;

	unsigned rootdist;

	int rclass;

	std::string name;
	
	nonterminal(const std::string &name_, size_t index_);
	~nonterminal();

	inline bool has_direct_descendent(const nonterminal *nt) const
		{ return equivalent_nts.find(nt) != equivalent_nts.end(); }

	inline bool isterminal() const {
		return name.size() > 2 && name[0] == '_' && name[1] == 'S';
	}

	inline bool ispartial() const {
		return !name.empty() && name[name.length()-1] != '_';
	}

	void addproduction(production *P);

	inline bool derives(const nonterminal *q) const { return derives_.find(q) != derives_.end(); }
	inline bool derived_by(const nonterminal *q) const { return derived_by_.find(q) != derived_by_.end(); }
};

struct grammar {
	const nonterminal *root;
	const nonterminal *partexpr_root;
	const nonterminal *compexpr_root;
	const nonterminal *single_expr_root;
	const nonterminal *multi_expr_root;
	const nonterminal *matrix_cell_root;
	std::vector<nonterminal *> nts;

	std::map<SemanticId, const production *> sidmap;

	grammar() : root(0), partexpr_root(0), compexpr_root(0), single_expr_root(0), multi_expr_root(0), matrix_cell_root(0) { }
	~grammar();
	
	
	void add_nonterminal(nonterminal *nt);
	const production *getcanonicalsidproduction(SemanticId sid) const;
	const production *production_for_tree(SemanticId sid, const std::vector<basic_tree *> &children) const;
	void setcanonicalsidproduction(SemanticId sid, const production *P);
};

void rebuild_production_lengths(grammar &G);
void read_grammar(std::istream &is, grammar &G);

std::istream &operator>>(std::istream &is, grammar &G);


}


#endif
