#ifndef PARSER_H_
#define PARSER_H_


#include "parser-fwd.h"
#include "expr-node.h"
#include "grammar.h"
#include "bitvec.h"
#include "ordered-segments.h"
#include "symbols.h"
#include "relation.h"

#include <list>
#include <map>
#include <queue>
#include <set>
#include <ostream>


namespace scg
{


struct match_score {
	const prototype *proto;
	double feature_score;
	double matcher_score;
	double score;

	match_score() :
		proto(0),
		feature_score(0),
		matcher_score(0),
		score(0)
		{ }

	static bool score_order_op(const match_score &left, const match_score &right)
		{ return left.score < right.score; }
};


struct segment;
struct group;

struct stroke {
	RawStroke input;
	segment *parent;

	stroke() : parent(0) { }
};


struct segment {
	static size_t INVALID_SEGMENT_POS;

	size_t pos;
	RawStroke stk;
	Rect<long> bounds;

	std::list<stroke *> children;
	std::list<group *> parents;

	segment() : pos(INVALID_SEGMENT_POS) { }
};


struct group_tags {
	enum {
		KNOWN_SYMBOL = 1,
		DEAD_GROUP = 2
	};
};


struct group {
	double proximity_score;
	double stack_score;
	double feature_score;
	double boosted_score;

	double weight;

	NormalizedStrokeGroup strokes;
	Rect<long> bounds;

	bitvec bits;

	long tags;

	std::list<segment *> children;

	std::vector<match_score> matches;

	std::vector<match_score> final_matches;

	group() :
		proximity_score(0),
		stack_score(0),
		feature_score(0),
		boosted_score(0),
		weight(0),
		tags(0)
		{ }
};

#if 0


template <typename T>
struct parse_result {
	expression_node tree;
	T src;

	parse_result(const expression_node &tree_, const T &src_) : tree(tree_), src(src_) { }

	parse_result(const parse_result &rhs) : tree(rhs.tree), src(rhs.src) { }

	inline bool operator<(const parse_result &rhs) const
		{
			double d = tree->score() - rhs.tree->score();
			if (std::abs(d) < std::numeric_limits<double>::epsilon()) {
				return tree->ccstr() < rhs.tree->ccstr();
			}
			return d > 0.0;
		}
};


typedef parse_result<std::vector<size_t> > links_parse_result;
typedef parse_result<size_t> cell_parse_result;

typedef std::set<links_parse_result> links_expression_heap;
typedef std::set<cell_parse_result> cell_expression_heap;


struct parse_links;


/*struct typed_relscore {
	virtual ~typed_relscore() { }
	virtual double lookup(const std::vector<raw_expression_node *> &children, const parse_links &links, std::vector<unsigned>::const_iterator curr, std::vector<unsigned>::const_iterator last) const = 0;
};


struct nt_relscore : public typed_relscore {
	nt_relscore() : typed_relscore() { }

	std::map<rclass_t, typed_relscore *> scores;

	double lookup(const std::vector<raw_expression_node *> &children, const parse_links &links, std::vector<unsigned>::const_iterator curr, std::vector<unsigned>::const_iterator last) const;
};

struct t_relscore : public typed_relscore {
	explicit t_relscore(const score_combiner &score_) : typed_relscore(), score(score_) { }

	score_combiner score;

	double lookup(const std::vector<raw_expression_node *> &children, const parse_links &links, std::vector<unsigned>::const_iterator curr, std::vector<unsigned>::const_iterator last) const;
};
*/

struct parse_cell;

struct parse_links {
	std::vector<parseref> children;

	std::vector<parse_cell *> cells;
	std::list<std::vector<size_t> > next_index_sets;
	std::set<std::vector<size_t> > expired_index_sets;
	std::vector<expression_node> ranked_trees;
	links_expression_heap cached_trees;

	std::vector<size_t> multiclass_children;
	std::vector<size_t> uniclass_children;

	//typed_relscore *score;
	double terminal_score;
	
	const production *P;

	explicit parse_links(const production *P_) : P(P_), terminal_score(0.0)
	{
		if (P_) {
			children.resize(P_->rhs.size());
		}
	}

	~parse_links()
		{ }

	//inline bool operator>(const parse_links &rhs) const
	//	{ return score > rhs.score; }
	
	inline bool operator==(const parse_links &rhs) const
		{ return children == rhs.children; }
	inline bool operator!=(const parse_links &rhs) const
		{ return children != rhs.children; }
};


class MatrixAnalyzer;
class MatrixIterator;

struct parse_cell {
	static const size_t INVALID_SRC;

	std::vector<parse_links> links;
	
	MatrixAnalyzer *matrix;
	MatrixIterator *mxiter;
	

	std::list<parse_cell *> parents;

	ordered_segments *segs;
	const nonterminal *nt;

	size_t prev_src;
	std::vector<size_t> next_best_ranks;
	std::vector<expression_node> ranked_trees;
	cell_expression_heap cached_trees;

	std::set<rclass_t> classes;

	int prepared_flags;
	
	explicit parse_cell(ordered_segments *segs_, const nonterminal *nt_);
};

typedef std::map<const nonterminal *, parse_cell *> nonterminal_parses;
typedef std::map<ordered_segments *, nonterminal_parses> parse_table;

/*
struct edit_op {
	enum {
		INSERTION,
		REMOVAL
	};

	int type;
	unsigned pos;

	edit_op(int type_, unsigned pos_) : type(type_), pos(pos_)
	{
		assert(type == INSERTION || type == REMOVAL);
	}
};
*/

struct MathRecognizer;

struct context {
	std::vector<stroke *> strokes;
	std::vector<segment *> segments;
	std::map<bitvec, group *> *groups;

	std::list<group *> newgroups;

	std::set<ordered_segments *> *locks;
	std::map<bitvec, ordered_segments *> *spans;
	parse_table table;

	std::vector<size_t> known;

	grammar G;

	MathRecognizer *wrapper;

	//std::list<edit_op> edits;

	explicit context(MathRecognizer *wrapper_) : groups(0), spans(0), wrapper(wrapper_) { }
};


struct parse_iterator {
	enum {
		DEFAULT = 0,
		ONLY_PRIMARY_SPAN_ALTERNATES = 1,
		ONLY_SEMANTIC_ALTERNATES = 2,
		LONG_FORM_ENABLED = 4
	};

	ordered_segments *start_span;
	int flags;
};

void restore_context(context *ctx, grammar &G, std::istream &is);
void save_context(context *ctx, std::ostream &os);

int initialize_context(context *ctx, const grammar &G);
void clear_context(context *ctx);


int add_strokes(context *ctx, const RawStroke *strokes, size_t nstrokes);
int remove_strokes(context *ctx, const size_t *indices, size_t nstrokes);


bool parse(context *ctx, ordered_segments *segs, const nonterminal *nt);

ExpressionIterator *make_iterator_internal(context *ctx, ordered_segments *span, const nonterminal *nt, int flags, raw_expression_node *parent);

bool prepare_for_iteration(const parse_iterator &it, context *ctx, ordered_segments *span, parse_cell *cell, bool found_locked_nt = false);
parse_cell *prepare_for_iteration(const parse_iterator &it, context *ctx, ordered_segments *segs, const nonterminal *nt);

expression_node extract_tree(const parse_iterator &it, ordered_segments *span, parse_cell *cell, size_t rank, bool wrap_string = true);
expression_node extract_tree_internal(const parse_iterator &it, ordered_segments *span, parse_cell *cell, size_t rank);

bool lock_semantics(context *ctx, ordered_segments *span, const raw_expression_node *node);
bool lock_expression(context *ctx, ordered_segments *span, const raw_expression_node *node);
//bool lock_children(context *ctx, ordered_segments *span, const raw_expression_node *parent);

enum {
	SPAN_UNLOCKED,
	SPAN_LOCKED_TYPE,
	SPAN_LOCKED_NODE,
	SPAN_LOCKED_CHILDREN
};


int is_locked(ordered_segments *span);
bool unlock(context *ctx, ordered_segments *span);
void unlock_all(context *ctx);


void dump_parse_table(std::ostream &os, const context *ctx);
void dump_span_table(std::ostream &os, const context *ctx);


//private functions

ordered_segments * lookup_span(context *ctx, const bitvec &bits, ordered_segments *src = 0, size_t dim = 0, size_t start = 0, size_t end = 0);
void append_edit_op(context *ctx, int type, size_t pos);


void dump_ink_tree(std::ostream &os, context *ctx);

double compute_membership(const relation *rel, const ordered_segments *segs1, const ordered_segments *segs2, const std::set<rclass_t> &cls1, const std::set<rclass_t> &cls2);

void remove_segment(context *ctx, segment *seg, bool kill_strokes = true);
int update_spans_on_removal(context *ctx, size_t pos);

void push_group_matches_to_parse_table(context *ctx, const group *grp);

void set_segment_positions(context *ctx);

void add_segment(context *ctx);
segment *create_default_segment(stroke *stk);
group *create_default_group(context *ctx, segment *seg);
std::map<bitvec, group *>::iterator remove_group(context *ctx, group *grp);

double match_group(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const symbol_info &info);
NormalizedStrokeGroup prepare_to_match(NormalizedStrokeGroup &input, const NormalizedStrokeGroup &model);

#endif

int recognizer_initialize();
void recognizer_shutdown();


}



#endif
