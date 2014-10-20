#ifndef RECOG_H_
#define RECOG_H_


#include <algorithm>
#include <limits>
#include <list>
#include <vector>

#include "feat.h"
#include "group.h"
#include "symbols.h"
#include "recoparms.h"
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

const stroke_feature_set &get_feature_set(const symbol_info &info);

//NormalizedStrokeGroup prepare_to_match(NormalizedStroke *first_input, NormalizedStroke *last_input, const NormalizedStrokeGroup &model);
NormalizedStrokeGroup prepare_to_match(const NormalizedStrokeGroup &input, const NormalizedStrokeGroup &model);

int recognizer_initialize();
void recognizer_shutdown();

std::map<bitvec, group *>::iterator remove_group(context *ctx, group *grp)


}


#endif

