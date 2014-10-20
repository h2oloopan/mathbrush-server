#ifndef RECOG_H_
#define RECOG_H_


#include <algorithm>
#include <limits>
#include <list>
#include <vector>

#include "feat.h"
#include "group.h"
#include "symbols.h"
//#include "recoparms.h"
#include "bitvec.h"
#include "ordered-segments.h"


namespace scg
{



struct match_score {
	const prototype *proto;
	double feature_score;
	double matcher_score;
	double score;

	/*inline bool operator<(const SegmentMatch &rhs) const
		{ return score < rhs.score; }*/
	static bool score_order_op(const match_score &left, const match_score &right)
		{ return left.score < right.score; }
};


struct tags {
	enum {
		DIRTY_SEGMENTATION = 1,
		DIRTY_GROUPING = (1 << 1),
		DIRTY_MATCHING = (1 << 2),
		DIRTY_SCORES = (1 << 3),
		MARK_FOR_REMOVAL = (1 << 4),
		MARK_FOR_ADDITION = (1 << 5)
	};
};


struct segment;
struct group;

struct stroke {
	RawStroke input;
	segment *parent;

	stroke() : parent(0), bits(0) { }
};


struct segment : public segment_base {
	RawStroke stk;
	Rect<long> bounds;

	std::list<stroke *> children;
	std::list<group *> parents;
};


struct group {
	double proximity_score;
	double stack_score;
	double feature_score;

	double weight;

	NormalizedStrokeGroup strokes;
	Rect<long> bounds;

	bitvec bits;

	std::list<segment *> children;

	std::vector<match_score> matches;

	std::vector<match_score> final_matches;

	group() :
		proximity_score(0),
		stack_score(0),
		feature_score(0),
		weight(0)
		{ }
};


struct context {
	std::map<bitvec, group *> groups;
	std::vector<segment *> segments;
	std::vector<stroke *> strokes;
};

void clear_context(reco_context *ctx);

double match_group(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const symbol_info &info);

segment * create_default_segment(stroke *stk);

int recognize(reco_context *ctx);


void remove_segment(reco_context *ctx, segment *seg, bool kill_strokes = true);
int add_strokes(reco_context *ctx, const RawStroke *strokes, unsigned nstrokes);
int remove_strokes(reco_context *ctx, const unsigned *indices, unsigned nstrokes);


/*
struct add_buf;
struct remove_buf;
struct preprocessed_add_buf;
struct preprocessed_remove_buf;

struct edit_visitor
{
	virtual void process(add_buf &buf) { }
	virtual void process(remove_buf &buf) { }
	virtual void process(preprocessed_add_buf &buf) { }
	virtual void process(preprocessed_remove_buf &buf) { }
};


struct edit_buf
{
    edit_buf() : next(0) { }
    virtual ~edit_buf() { delete next; }

	 virtual void process(edit_visitor &p) = 0;
    
    edit_buf *next;
};


void insert_edit_buf(edit_buf **head, edit_buf *entry);
void insert_edit_buf_at(edit_buf *inspt, edit_buf *entry);
void replace_edit_buf_at(edit_buf *inspt, edit_buf *entry);


struct remove_buf : public edit_buf
{
    std::vector<unsigned> indices;

    remove_buf() : edit_buf() { } 

	 void reset() { indices.clear(); }

	 void process(edit_visitor &p) { p.process(*this); }
};

struct preprocessed_remove_buf : public edit_buf
{
	 std::vector<unsigned> indices;
    std::vector<unsigned> segment_ids;
    std::vector<BboxCandidate *> boxes;

    preprocessed_remove_buf() : edit_buf() { }

	 void reset()
	 {
		indices.clear();
		segment_ids.clear();
		boxes.clear();
	 }

	 void process(edit_visitor &p) { p.process(*this); }
};

struct add_buf : public edit_buf
{
    RawStroke *strokes;
    unsigned nstrokes;
    long *ids;
    
    add_buf() : edit_buf(), strokes(0), nstrokes(0), ids(0) { }

	 void reset()
	 {
	 	delete[] strokes;
		strokes = 0;
		nstrokes = 0;
		delete[] ids;
		ids = 0;
	 }

	 void process(edit_visitor &p) { p.process(*this); }
};


struct preprocessed_add_buf : public edit_buf
{
	std::vector<unsigned> indices;
	std::vector<RecognitionSegment *> segments;

	preprocessed_add_buf() : edit_buf() { }

	void reset()
	{
		indices.clear();
		segments.clear();
	}

	void process(edit_visitor &p) { p.process(*this); }
};


int recognize_incremental(context *ctx,
                          remove_buf &removals, add_buf &insertions,
								  preprocessed_remove_buf &parse_removals, preprocessed_add_buf &parse_insertions);
*/

const stroke_feature_set &get_feature_set(const symbol_info &info);

NormalizedStrokeGroup prepare_to_match(NormalizedStroke *first_input, NormalizedStroke *last_input, const NormalizedStrokeGroup &model);
NormalizedStrokeGroup prepare_to_match(NormalizedStrokeGroup &input, const NormalizedStrokeGroup &model);

int recognizer_initialize();
void recognizer_shutdown();

std::map<bitvec, group *>::iterator remove_group(context *ctx, group *grp)


}


#endif

