#ifndef RECO_TYPES_H_
#define RECO_TYPES_H_

#include "bitvec.h"
#include "stroke.h"
#include "group.h"
#include "feat.h"
#include "rect.h"
#include <list>
#include <vector>
#include <map>
#include <limits>

namespace scg {

//struct prototype;
struct symbol;
struct nonterminal;

struct recoscore {
	double hint;
	double score;
	
	recoscore() : hint(0), score(0) { }
	recoscore(double sc) : hint(0), score(sc) { }
	recoscore(double hint_, double sc) : hint(hint_), score(sc) { }
	inline void sethint(double h) { hint = h; };
	inline void addhint(double h) { hint += h; };
	inline void rmhint(double h) { hint -= h; };

	inline bool operator<(const recoscore &rhs) const { return hint < rhs.hint || (hint == rhs.hint && score < rhs.score); }
	inline bool operator>(const recoscore &rhs) const { return hint > rhs.hint || (hint == rhs.hint && score > rhs.score); }
	inline bool operator<=(const recoscore &rhs) const { return hint < rhs.hint || (hint == rhs.hint && score <= rhs.score); }
	inline bool operator>=(const recoscore &rhs) const { return hint > rhs.hint || (hint == rhs.hint && score >= rhs.score); }
	inline bool operator==(const recoscore &rhs) const { return hint == rhs.hint && score == rhs.score; }
	inline bool operator!=(const recoscore &rhs) const { return hint != rhs.hint || score != rhs.score; }
	inline bool operator<(double rhs) const { return hint == 0 && score < rhs; }
	inline bool operator>(double rhs) const { return hint > 0 || score > rhs; }
	inline bool operator<=(double rhs) const { return hint == 0 && score <= rhs; }
	inline bool operator>=(double rhs) const { return hint > 0 || score >= rhs; }
	inline bool operator==(double rhs) const { return hint == 0 && score == rhs; }
	inline bool operator!=(double rhs) const { return hint != 0 || score != rhs; }

	inline recoscore &operator+=(const recoscore &rhs) { hint += rhs.hint; score += rhs.score; return *this; }
	inline recoscore &operator-=(const recoscore &rhs) { hint -= rhs.hint; score -= rhs.score; return *this; }
	inline recoscore operator+(const recoscore &rhs) const { return recoscore(hint+rhs.hint, score+rhs.score); }
	inline recoscore operator-(const recoscore &rhs) const { return recoscore(hint-rhs.hint, score-rhs.score); }
	inline recoscore &operator*=(double rhs) { score *= rhs; return *this; }
	inline recoscore &operator/=(double rhs) { score /= rhs; return *this; }
	inline recoscore operator*(double rhs) const { return recoscore(hint, score*rhs); }
	inline recoscore operator/(double rhs) const { return recoscore(hint, score/rhs); }
};


std::ostream &operator<<(std::ostream &os, const recoscore &rhs);

struct match_score {
	symbol *S;
	double feature_score;
	double matcher_score;
	double score;
	double bag_score;
	recoscore final_score;

	match_score() :
		S(0),
		feature_score(0.0),
		matcher_score(std::numeric_limits<double>::infinity()),
		score(0.0),
		bag_score(1.0)
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

//class math_recognizer_base;

struct segment {
	static size_t INVALID_SEGMENT_POS;

	size_t pos;
	RawStroke stk;
	Rect<long> bounds;

	std::list<stroke *> children;
	std::list<group *> parents;

	//math_recognizer_base *ctx;
	segment() : pos(INVALID_SEGMENT_POS) { }
};


struct group_tags {
	enum {
		KNOWN_SYMBOL = 1,
		DEAD_GROUP = 2,
		LOADED_GROUP = 4
	};
};


struct group {
	double proximity_score;
	double stack_score;
	double feature_score;
	double container_score;
	double boosted_score;

	double weight;
	double pnil;

	RawStrokeGroup strokes;
	NormalizedStrokeGroup nstrokes;
	Rect<long> bounds;

	bitvec bits;

	long tags;

	std::list<segment *> children;

	std::map<const nonterminal *, match_score> matches;
	//std::vector<match_score> matches;
	//std::vector<match_score> final_matches;
	//group *parents[2];
	void mulscore(double x);

	group() :
		proximity_score(0),
		stack_score(0),
		feature_score(0),
		container_score(0),
		boosted_score(0),
		weight(0),
		tags(0),
		pnil(0)
		{ /*parents[0] = parents[1] = 0;*/ }
};

struct prototype;
bool prepare_to_match(group *input, const prototype *model, std::vector<size_t> &inputorder);
extern stroke_feature_set default_feature_set;
extern stroke_feature_set small_feature_set;

}

#endif
