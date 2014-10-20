#ifndef MATHRECOGNIZER_PRIVATE_H_
#define MATHRECOGNIZER_PRIVATE_H_


#include "symbols.h"
#include "MathRecoTypes.h"
#include "MathRecognizer.h"
#include "intrpr.h"
#include "expr-node.h"
#include "grammar.h"

#include <istream>
#include <list>
#include <set>
#include <map>
#include <vector>
#include <queue>


namespace scg
{


class BoxLinkEstimator;

//int rebuild_global_symbols_db();
//const symbols_db *global_symbols_db();
//symbols_db *global_symbols_db_wr();
//void kill_global_symbols_db();
//symbols_db *global_autotrain_db();

typedef std::map<const production *, interpreter *> production_parser_map;
struct production_parse_table {
	interpreter *nt_parser;
	production_parser_map tab;
	production_parse_table() : nt_parser(0) { }
};
typedef std::map<const nonterminal *, production_parse_table> nt_parse_table;
typedef std::map<const ordered_segments *, nt_parse_table> parse_table;

class prodparser;
struct ptab;

int RebuildGrammarData();


ptab *get_ptab();


class math_recognizer_base : public MathRecognizer {
public:
	math_recognizer_base();
	explicit math_recognizer_base(std::istream &is);
	math_recognizer_base(const MathSymbol *symbols, size_t n);

	virtual ~math_recognizer_base();

	void release() { delete this; }

	int Clear();
	
	void Translate(long x, long y);

	int Save(std::ostream &os);

	ExpressionTree *ConvertSurrogateTreeToExpressionTree(const SurrogateTree *);

	ExpressionIterator *CreateIteratorForExpression(const ExpressionTree *expr, bool ownexpr, bool wrap_mathml);	

	const ExpressionTree *GetTopExpression();
	ExpressionIterator *CreateDefaultIterator(bool wrap_mathml = true);
	ExpressionIterator *CreateSubtreeIterator(const ExpressionTree *base, bool wrap_mathml = true);
	ExpressionIterator *CreateIterator(ExpressionTree *base, bool wrap_mathml = true);
	ExpressionIterator *CreateIteratorForStrokes(const RawStroke **strokes, size_t nstrokes, bool wrap_mathml = true);
	ExpressionIterator *CreateIteratorForStrokes(const nonterminal *nt, const RawStroke **strokes, size_t nstrokes, bool wrap_mathml = true);
	ExpressionIterator *CreateIteratorForStrokesByIndex(const size_t *strokes, size_t nstrokes, bool wrap_mathml = true);
	ExpressionIterator *CreateIteratorForStrokesByIndex(const nonterminal *nt, const size_t *strokes, size_t nstrokes, bool wrap_mathml = true);

	ExpressionIterator *CreateDefaultSemanticIterator(SemanticId sid, bool wrap_mathml = true);
	ExpressionIterator *CreateSemanticIteratorForStrokes(SemanticId sid, const RawStroke **strokes, size_t nstrokes, bool wrap_mathml = true);
	ExpressionIterator *CreateSemanticIteratorForStrokesByIndex(SemanticId sid, const size_t *strokes, size_t nstrokes, bool wrap_mathml = true);
	
	int AddKnownSymbols(const MathSymbol *symbols, size_t n);
	int RemoveKnownSymbols(const MathSymbol *symbols, size_t n);


	int RemoveStrokesByIndex(const size_t *strokes, size_t nstrokes);

	int AddStrokes(const RawStrokeGroup &group);
	int AddStrokes(const RawStroke *strokes, size_t nstrokes);

	int RemoveStroke(const RawStroke *stk);
	size_t GetStrokes(const RawStroke **strokes, size_t n);

protected:
	virtual int RemoveStrokesByIndex(const std::vector<size_t> &indices);

public:
	double avg_stroke_size;
	std::vector<stroke *> strokes;
	std::vector<segment *> segments;
	std::map<bitvec, group *> *groups;
	std::list<group *> newgroups;
	std::map<bitvec, ordered_segments *> *spans;
	std::vector<std::vector<size_t> > known;
	parse_table tab;

	/*struct resetspec {
		const nonterminal *nt;
		bitvec bits;
		size_t childi;
		resetspec(const nonterminal *nt_, const bitvec &bits_, size_t childi_) : nt(nt_), bits(bits_), childi(childi_) { }
	};*/
	//std::queue<resetspec> resets;
	std::queue<interpreter *> resets;

private:
	void init();
	void remove_segment(size_t segi, bool kill_strokes = true);
	std::map<bitvec, group *>::iterator remove_group(group *grp);
	void update_groups_on_removal(size_t i);
	int update_spans_on_removal(const segment *seg);
	void set_segment_positions();
	int remove_strokes(const size_t *indices, size_t nstrokes);
	void updategroups(size_t i);
	void handlenewgroups();
	void clear_dead_groups();
	void cleartablesfor(size_t pos);
	void rebalance_groups();
	void finalize_match_scores(group *grp);
	int add_strokes(const RawStroke *strokes, size_t nstrokes);
	void add_segment();
	segment *merge_segment(segment *seg);
	int update_proximity_groups(std::map<bitvec, group *> &testgrps, std::map<bitvec, group *> &pgroups, group *single_seg_grp, size_t segi);
	int update_stack_groups(std::list<group *> &cmpgroups, std::list<group *> &sgroups);
	group *build_stack_group(group *g1, group *g2);
	group *create_default_group(segment *seg, size_t segi);
	ordered_segments *find_span(size_t nsegments, const size_t *strokes, size_t nstrokes, bool create = false);
	void clear();

	void doresets();

	int save(std::ostream &os);
	int restore(std::istream &is);

public:
	interpreter *getparser(const nonterminal *nt, const ordered_segments *segs);
	//interpreter *getfullmux(const ordered_segments *segs);
	interpreter *mkparser(const nonterminal *nt, const ordered_segments *segs, bool ign_overlap = false);
	interpreter *mkprodintrpr(const production *P, const ordered_segments *segs, bool ign_overlap = false);
	interpreter *mktermparser(const nonterminal *nt, const ordered_segments *segs);
	const ordered_segments *getsegs(const bitvec &bits) const;
	ordered_segments *getsegs(const bitvec &bits, const ordered_segments *super, size_t d, size_t start, size_t end, bool ign_overlap);

private:
	friend class prodintrprbuilder;
	friend class external_iterator;
	friend class reader;
};

}


#endif
