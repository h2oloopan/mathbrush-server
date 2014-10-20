#include "relation.h"
#include "mathrecognizer-private.h"
#include "utils.h"
#include "grammar.h"
#include "expr-node.h"
#include "segment.h"
#include "rect.h"
#include "dist.h"
#include "stream-defs.h"
#include "stats.h"
#include "iohelp.h"

#include <fstream>
#include <map>
#include <vector>
#include <complex>
#include <cstring>
#include <algorithm>
#include <cmath>


namespace scg
{


static relation *relations[NUM_RELATIONS];
static relation *grprelations[NUM_RELATIONS];

static unsigned rel_orderings[] = { 0, 0, 0, 1, 0 };

static std::string box_name("_BOX");
static std::string symbol_name("_SYMBOL");
static std::string aggregate_name("_aggregate");


const rclass_t INVALID_CLASS = static_cast<rclass_t>(-1);
const rclass_t AGGREGATE_CLASS = 0;
const rclass_t SYMBOL_CLASS = 1;
const rclass_t BOX_CLASS = 2;
static const unsigned NUM_SPECIAL_CLASSES = 3;

// these following classes must match the order given in symbols.cc
const rclass_t BASELINE_CLASS = 3;
const rclass_t ASCENDER_CLASS = 4;
const rclass_t DESCENDER_CLASS = 5;
const rclass_t CENTERED_CLASS = 6;
const rclass_t HALF_ASCENDER_CLASS = 7;
const rclass_t HALF_ASCENDER_DESCENDER_CLASS = 8;
const rclass_t ROOT_CLASS = 9;
const rclass_t HORIZONTAL_CLASS = 10;
const rclass_t PUNCTUATION_CLASS = 11;
const rclass_t BIG_CENTERED_CLASS = 12;
const rclass_t EXTENDER_CLASS = 13;
const rclass_t LARGE_EXTENDER_CLASS = 14;
//const rclass_t RANGE_CLASS = 15;
const rclass_t NCLASSES = 15;
const rclass_t RELCLASS_MERGE = 31;

//const rclass_t EXTENDER_CLASS = cl_++;
//const rclass_t HALF_DESCENDER_CLASS = 8;
//const rclass_t FENCE_CLASS = 11;

//const rclass_t RANGE_OP_CLASS = 13;

rclass_t
tag_to_rclass(const std::string &tag) {
	if (tag == "_aggregate") {
		return AGGREGATE_CLASS;
	}
	else if (tag == "_box") {
		return BOX_CLASS;
	}
	else if (tag == "_symbol") {
		return SYMBOL_CLASS;
	}
	else if (tag == "baseline") {
		return BASELINE_CLASS;
	}
	else if (tag == "ascender") {
		return ASCENDER_CLASS;
	}
	else if (tag == "descender") {
		return DESCENDER_CLASS;
	}
	else if (tag == "extender") {
		return EXTENDER_CLASS;
	}
	else if (tag == "half-ascender") {
		return HALF_ASCENDER_CLASS;
	}
	//else if (tag == "half-descender") {
	//	return HALF_DESCENDER_CLASS;
	//}
	else if (tag == "half-ascender-descender") {
		return HALF_ASCENDER_DESCENDER_CLASS;
	}
	else if (tag == "centered") {
		return CENTERED_CLASS;
	}
	else if (tag == "big-centered") {
		return BIG_CENTERED_CLASS;
	}
	//else if (tag == "fence") {
	//	return FENCE_CLASS;
	//}
	else if (tag == "large-extender") {
		return LARGE_EXTENDER_CLASS;
	}
	//else if (tag == "range") {
	//	return RANGE_CLASS;
	//}
	else if (tag == "root") {
		return ROOT_CLASS;
	}
	else if (tag == "horizontal") {
		return HORIZONTAL_CLASS;
	}
	else if (tag == "punctuation") {
		return PUNCTUATION_CLASS;
	}
	return (rclass_t)-1;
}

std::string
rclass_to_str(rclass_t cls) {
	switch (cls) {
	case AGGREGATE_CLASS:
		return "aggregate";
	case BOX_CLASS:
		return "box";
	case SYMBOL_CLASS:
		return "symbol";
	default:
		return rclass_to_tag(cls);
	}
}
std::string
rclass_to_tag(rclass_t cls) {
	switch (cls) {
	case AGGREGATE_CLASS:
	case BOX_CLASS:
	case SYMBOL_CLASS:
		return "";
	case BASELINE_CLASS:
		return "baseline";
	case ASCENDER_CLASS:
		return "ascender";
	case DESCENDER_CLASS:
		return "descender";
	case EXTENDER_CLASS:
		return "extender";
	case HALF_ASCENDER_CLASS:
		return "half-ascender";
	//case HALF_DESCENDER_CLASS:
	//	return "half-descender";
	case HALF_ASCENDER_DESCENDER_CLASS:
		return "half-ascender-descender";
	case CENTERED_CLASS:
		return "centered";
	case BIG_CENTERED_CLASS:
		return "big-centered";
	//case FENCE_CLASS:
	//	return "fence";
	case LARGE_EXTENDER_CLASS:
		return "large-extender";
	//case RANGE_CLASS:
	//	return "range";
	case ROOT_CLASS:
		return "root";
	case HORIZONTAL_CLASS:
		return "horizontal";
	case PUNCTUATION_CLASS:
		return "punctuation";
	default:
		return "?rclass?";
	}
}


static double
getml(const Rect<long> &bbox, rclass_t rc) {
	if (rc & ((1 << ASCENDER_CLASS) | (1 << ROOT_CLASS))) {
		if (rc & ((1 << DESCENDER_CLASS) | (1 << HALF_ASCENDER_DESCENDER_CLASS))) {
			return bbox.top + bbox.height()/3.0;
		}
		return bbox.top + bbox.height()/2.0;
	}
	else if (rc & (1 << HALF_ASCENDER_CLASS)) {
		if (rc & ((1 << DESCENDER_CLASS) | (1 << HALF_ASCENDER_DESCENDER_CLASS))) {
			return bbox.top + bbox.height()/5.0;
		}
		return bbox.top+bbox.height()/3.0;
	}
	else if (rc & (1 << CENTERED_CLASS)) {
		if (rc & ((1 << DESCENDER_CLASS) | (1 << HALF_ASCENDER_DESCENDER_CLASS))) {
			return bbox.top + bbox.height()/5.0;
		}
		return bbox.top+bbox.height()/2.0;
	}
	else if (rc & (1 << HALF_ASCENDER_DESCENDER_CLASS)) {
		return bbox.top + bbox.height()/5.0;
	}
	else if (rc & ((1 << BASELINE_CLASS))) {
		return bbox.top;
	}
	else if (rc & (1 << HORIZONTAL_CLASS)) {
			return bbox.top + bbox.height()/5.0;
	}
	else if (rc & (1 << PUNCTUATION_CLASS)) {
		return bbox.top - (long)std::max(bbox.height()*2.0, 0.125*TABLETPC_DPI);
	}
	else {
		return bbox.top + bbox.height()/2.0;
	}
}

double
getbl(const Rect<long> &bbox, rclass_t rc) {
/*	
	for (int i = 3; i < NCLASSES; ++i) {
		if (rc & (1 << i)) {
			return bl(bbox, i);
		}
	}
	return bbox.bottom;*/
	if ((rc & (1 << DESCENDER_CLASS))) {
		if (rc & ((1 << ASCENDER_CLASS) | (1 << ROOT_CLASS))) {
			return bbox.bottom - bbox.height()/3.0;
		}
		else if (rc & (1 << BIG_CENTERED_CLASS)) {
			return bbox.bottom - 3.0*bbox.height()/8;
		}
		else if (rc & ((1 << CENTERED_CLASS) | (1 << HALF_ASCENDER_CLASS) | (1 << HALF_ASCENDER_DESCENDER_CLASS))) {
			return bbox.bottom - 2.0*bbox.height()/5;
		}
		return bbox.bottom - bbox.height()/2;
	}
	else if (rc & (1 << HALF_ASCENDER_DESCENDER_CLASS)) {
		if (rc & ((1 << ASCENDER_CLASS) | (1 << ROOT_CLASS))) {
			return bbox.bottom - bbox.height()/3.0;
		}
		else if (rc & (1 << BIG_CENTERED_CLASS)) {
			return bbox.bottom - 3.0*bbox.height()/8;
		}
		return bbox.bottom - 2.0*bbox.height()/5;
	}
	else if (rc & ((1 << BASELINE_CLASS) | (1 << ASCENDER_CLASS) | (1 << HALF_ASCENDER_CLASS) | (1 << ROOT_CLASS))) {
		return bbox.bottom;
	}
	else if (rc & (1 << BIG_CENTERED_CLASS)) {
		return bbox.bottom + bbox.height()/4.0;
	}
	else if (rc & (1 << CENTERED_CLASS)) {
		return bbox.bottom + bbox.height()/2.0;
	}
	else if (rc & (1 << HORIZONTAL_CLASS)) {
		return bbox.bottom + (long)std::max(bbox.height()*2.0, 0.125*TABLETPC_DPI);
	}
	else if (rc & (1 << PUNCTUATION_CLASS)) {
		return (bbox.top + bbox.bottom) / 2.0;
	}
	else {
		return bbox.bottom;
	}
}

static double FreeRelationDist = 0.0;
static void calc_parms() {
	FreeRelationDist = TABLETPC_DPI * GetParameterDouble("FreeRelationDist");
}
static int _e = RegisterParameterCallback(&calc_parms);

static long
inset(long d) {
	return static_cast<long>(std::min(0.25 * d, FreeRelationDist));
}



double
get_angle_y(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2)
{
	//return static_cast<double>(box2.top - box1.bottom) / std::max<double>(height(box1), FreeRelationDist);
	double x1,y1;
	double x2,y2;
	x1 = tail->span->x_center(); y1 = tail->span->y_center();
	x2 = head->span->x_center(); y2 = head->span->y_center();
	//get_from_centroid_y(box1, class1, x1, y1);
	//get_to_centroid_y(box2, class2, x2, y2);
	Rect<long> bounds1 = tail->bounds();
	/*Rect<long> bounds2 = head->bounds();
	x1 = (bounds1.left+bounds1.right)/2;
	y1 = (bounds1.top+bounds1.bottom)/2;
	x2 = (bounds2.left+bounds2.right)/2;
	y2 = (bounds2.top+bounds2.bottom)/2;*/
	VERBOSE(*verb_out << "angle from (" << x1 << "," << y1 << ")->(" << x2 << "," << y2 << ") = " << std::atan2(y2-y1, x2-x1) << std::endl);
	return std::atan2(y2-y1, x2-x1);
}

static double
getypt(const interpretation *intrp, rclass_t rc) {
	return getml(intrp->bounds(), rc);
}

double
get_angle_x(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2) {
	double x1,y1;
	double x2,y2;
	//x1 = tail->span->x_center();
	//y1 = tail->span->y_center();
	//x2 = head->span->x_center();
	//y2 = head->span->y_center();
	Rect<long> bounds1 = tail->bounds();
	Rect<long> bounds2 = head->bounds();
	x1 = (bounds1.left+bounds1.right)/2;
	y1 = (bounds1.top+bounds1.bottom)/2;
	x2 = (bounds2.left+bounds2.right)/2;
	y2 = (bounds2.top+bounds2.bottom)/2;
	//y1 = getypt(tail, class1);
	//y2 = getypt(head, class2);
	//x1 = tail->span->x_center();
	//y1 = getypt(tail, class1);
	//y1 = tail->span->y_center();
	//x2 = head->span->x_center();
	//y2 = getypt(head, class2);
	//y2 = head->span->y_center();
	//x1 = segs1->x_center(); x2 = segs2->x_center();
	//x1 = tail->span->x_center();
	//x2 = head->span->x_center();
	//y1 = getypt(tail, class1);
	//y2 = getypt(head, class2);
	//get_from_centroid_x(box1, class1, x1, y1);
	//get_to_centroid_x(box2, class2, x2, y2);
	//return static_cast<double>(y2-y1) / std::max<double>(height(box1), FreeRelationDist);//height(box2));
	VERBOSE(*verb_out << "angle from (" << x1 << "," << y1 << ")->(" << x2 << "," << y2 << ") = " << std::atan2((double)y2-y1, x2-x1) << std::endl);
	return std::atan2((double)y2-y1, x2-x1);
}

static long
distance_threshold(const Rect<long> &box1, const Rect<long> &box2)
{
	//return std::max(std::max(width(box1), height(box1)), std::max(width(box2), height(box2)));
	long b1thresh = (box1.width() + box1.height())/2;
	long b2thresh = (box2.width() + box2.height())/2;
	//return (long)std::max<double>(FreeRelationDist, b1thresh + b2thresh);
	//return std::min<double>(FreeRelationDist, std::max<long>(FreeRelationDist/2, (b1thresh+b2thresh)/4));
	long mx_thresh = std::max(b1thresh, b2thresh);
	long mn_thresh = 2 * std::min(b1thresh, b2thresh);
	return (long)std::min<double>(FreeRelationDist, std::max<double>(FreeRelationDist / 2, std::min(mx_thresh, mn_thresh)));
}

static double
distance_score(const Rect<long> &box1, const Rect<long> &box2)
{
	long thresh = distance_threshold(box1, box2);
	double dist = std::sqrt((double)dist_sq(box1, box2));
	//return std::max(0.0, 1.0 - dist/thresh);
	if (dist < thresh) {
		return 1.0;
	}
	else {
		return std::max(0.0, 1.0 - (dist - thresh) / (4 * thresh));
	}
}

static double
nominal_top(rclass_t c) {
	switch (c) {
	case AGGREGATE_CLASS:
	case SYMBOL_CLASS:
	case BOX_CLASS: return 1;
	case ASCENDER_CLASS: return 1;
	case DESCENDER_CLASS: return 0.5;
	case CENTERED_CLASS: return 0.75;
	case BIG_CENTERED_CLASS: return 1.0-1.0/6;
	case HORIZONTAL_CLASS: return 0.6;
	case BASELINE_CLASS: return 0.5;
	case HALF_ASCENDER_CLASS: return 0.75;
	case HALF_ASCENDER_DESCENDER_CLASS: return 0.75;
	case ROOT_CLASS: return 1;
	case PUNCTUATION_CLASS: return 0.1;
	case EXTENDER_CLASS: return 1.15;
	default:
		THROW_ERROR(E_INVALID, "nominal_bot(): unknown class " << c);
	}
}

static double
nominal_bot(rclass_t c) {
	switch (c) {
	case AGGREGATE_CLASS:
	case SYMBOL_CLASS:
	case BOX_CLASS: return 0;
	case ASCENDER_CLASS: return 0;
	case DESCENDER_CLASS: return -0.5;
	case CENTERED_CLASS: return 0.25;
	case BIG_CENTERED_CLASS: return 1.0/6;
	case HORIZONTAL_CLASS: return 0.4;
	case BASELINE_CLASS: return 0;
	case HALF_ASCENDER_CLASS: return 0;
	case HALF_ASCENDER_DESCENDER_CLASS: return -0.5;
	case ROOT_CLASS: return 0;
	case PUNCTUATION_CLASS: return -0.1;
	case EXTENDER_CLASS: return -0.15;
	default:
		THROW_ERROR(E_INVALID, "nominal_bot(): unknown class " << c);
	}
}
double
nominal_height(rclass_t c) {
	double top = 0, bot = 1;
	for (int i = 3; i < NCLASSES; ++i) {
		if (c & (1 << i)) {
			top = std::max(top, nominal_top(i));
			bot = std::min(bot, nominal_bot(i));
		}
	}
	if (top < bot) return 1;
	return top - bot;
}

double
gettl(const Rect<long> &bbox, rclass_t c) {
	switch (c) {
	case AGGREGATE_CLASS:
	case SYMBOL_CLASS:
	case BOX_CLASS:
		return bbox.top;
	case ASCENDER_CLASS:
		return bbox.top;
	case DESCENDER_CLASS:
		return bbox.top - 0.5*bbox.height();
	case CENTERED_CLASS:
		return bbox.top - 0.5*bbox.height();
	case HORIZONTAL_CLASS:
		return bbox.top - std::max(bbox.height()*2.0, 0.125*TABLETPC_DPI);
	case BASELINE_CLASS:
		return bbox.top - bbox.height();
	case HALF_ASCENDER_CLASS:
		return bbox.top - bbox.height()/3.0;
	case HALF_ASCENDER_DESCENDER_CLASS:
		return bbox.top - bbox.height()/5.0;
	case ROOT_CLASS:
		return bbox.top;
	case PUNCTUATION_CLASS:
		return bbox.top;
	default:
		THROW_ERROR(E_INVALID, "gettl(): unknown class " << c);
	}
}

static double
right_height_score(long h1, long h2, rclass_t c1, rclass_t c2) {
	double nh1 = nominal_height(c1);
	double nh2 = nominal_height(c2);
	double nr = nh1/nh2;
	double r = (double)h1/h2;
	double a = std::max(r/nr, nr/r);
	VERBOSE(*verb_out << "right_height_score: heights " << h1 << "/" << h2 << " and nominal heights " << nh1 << "/" << nh2 << " give ratio " << a << " for score " << 1/a << std::endl);
	return 1/std::pow(a, 1.0/4);
}

static double
angle_height_score(long h1, long h2, rclass_t c1, rclass_t c2) {
	double nh1 = nominal_height(c1);
	double nh2 = nominal_height(c2)/2;
	double nr = nh1/nh2;
	double r = (double)h1/h2;
	double a = std::max(r/nr, nr/r);
	VERBOSE(*verb_out << "angle_height_score: heights " << h1 << "/" << h2 << " and nominal heights " << nh1 << "/" << nh2 << " give ratio " << a << " for score " << 1/a << std::endl);
	return 1/std::pow(a, 1.0/4);
}


static bool
useheight(int class1, int class2) {
	return !((class1|class2) & ((1 << AGGREGATE_CLASS) | (1 << BOX_CLASS)));
	//return class1 != AGGREGATE_CLASS && class2 != AGGREGATE_CLASS;
}

struct fuzzy_relation : public relation {
	fuzzy_relation(const std::string &name_, double t0_, double t1_, double t2_, get_angle_fn t_);

	std::string name;
	double t0, t1, t2;
	
	get_angle_fn t;

	virtual double membership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const = 0;

protected:
	double basemembership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const;
};


static fuzzy_relation *fuzrelations[NUM_RELATIONS];

fuzzy_relation::fuzzy_relation(const std::string &name_, double t0_, double t1_, double t2_, get_angle_fn t_)
	: relation(name_), t0(t0_), t1(t1_), t2(t2_), t(t_)
	{ }

double
fuzzy_relation::basemembership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tailatt, int headatt) const {
	Rect<long> box1 = tail->bounds();
	Rect<long> box2 = head->bounds();
	double ct = t(tail, head, class1, class2);
	ct *= 180.0/M_PI;
	
	double angle_score;
	if (ct <= t0) {
		return 0.0;
	}
	else if (ct < t1) {
		angle_score = (ct - t0) / (t1 - t0);
	}
	else if (ct < t2) {
		angle_score = 1.0 - (ct - t1) / (t2 - t1);
	}
	else {
		return 0.0;
	}

	//std::cout << name << " " << ct << std::endl;
	double overlap = overlap_proportion(box1, box2);
	double score = angle_score * (1-0.75*overlap);
	double dist_score = 1;
	if (tailatt != attach_mode::GROUP || headatt != attach_mode::GROUP) {
		dist_score = distance_score(box1, box2);
		score *= dist_score;
	}
	VERBOSE(*verb_out << "measuring " << name << " for boxes " << box1 << " and " << box2 << "; classes " << class1 << " and " << class2 << std::endl);
	VERBOSE(*verb_out << "interpretations " << tail->str() << " and " << head->str() << std::endl);
	VERBOSE(*verb_out << " scores " << angle_score << " *? " << dist_score << " * (1-" << 0.5*overlap << ") = " << score << std::endl);
	/*
	if ((class1 == PUNCTUATION_CLASS && class2 != ROOT_CLASS) || (class2 == PUNCTUATION_CLASS && class1 != ROOT_CLASS)) {
		overlap = 0.0;
	}*/
	return score;
}

struct B_rel : public fuzzy_relation {
	B_rel() : fuzzy_relation("B", 0, 90, 180, &get_angle_y) { }
	double membership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const {
		Rect<long> box1 = tail->bounds();
		Rect<long> box2 = head->bounds();
		double score = basemembership(tail, head, class1, class2, tatt, hatt);
		long left = std::max(box1.left, box2.left);
		long right = std::min(box1.right, box2.right);
		long over = std::max<long>(right - left, 0);
		double prop;
		if (!tail->src || !head->src) { // dummy interpretations for pre-parse relation filtering
			prop = 1;
		}
		else {
			bool tailline = derivesterminal(tail) && !std::strcmp(tail->str(), "horzline");
			bool headline = derivesterminal(head) && !std::strcmp(head->str(), "horzline");
			bool taildot = derivesterminal(tail) && !std::strcmp(tail->str(), "dot");
			bool headdot = derivesterminal(head) && !std::strcmp(head->str(), "dot");
			if (tailline && !headline) {
				prop = (double)over / box2.width();
			}
			else if (headline && !tailline) {
				prop = (double)over / box1.width();
			}
			else if (headdot || taildot) {
				prop = 1.0;
			}
			else {
				prop = (double)over / std::min(box1.width(), box2.width());
			}
		}
		VERBOSE(*verb_out << " overlap proportion is " << prop << " giving final score " << score*prop*prop << std::endl);
		prop *= prop;
		score *= prop;
		return score;
	}
};
	
struct BR_rel : public fuzzy_relation {
	BR_rel() : fuzzy_relation("BR", -15, 35, 90, &get_angle_x) { }
	double membership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const {
		Rect<long> box1 = tail->bounds();
		Rect<long> box2 = head->bounds();
		double score = basemembership(tail, head, class1, class2, tatt, hatt);
		if (useheight(class1, class2)) {
			double bl1 = getbl(box1, class1);
			double bl2 = getbl(box2, class2);
			if (bl1 > 0 && bl2 > 0) {
				double free = std::max(0.25*box1.height(), 0.125*TABLETPC_DPI);
				double limit = std::max(1.5*box1.height(), 0.125*TABLETPC_DPI);
				double diff = std::min(limit+free, std::max(0.0, bl2 - bl1 + free));
				double blfac = std::pow(diff/(limit+free), 1.0/4);
				VERBOSE(*verb_out << "baselines " << bl1 << ',' << bl2 << " score " << diff/limit << " with diff " << diff << " and limit " << limit << " : " << blfac << std::endl);
				score *= blfac;
			}
			score *= angle_height_score(std::max(1l, box1.height()), std::max(1l, box2.height()), class1, class2);
		}
		VERBOSE(*verb_out << " final score: " << score << std::endl);
		return score;
	}
};

struct AR_rel : public fuzzy_relation {
	AR_rel() : fuzzy_relation("AR", -90, -35, 15, &get_angle_x) { }
	double membership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const {
		Rect<long> box1 = tail->bounds();
		Rect<long> box2 = head->bounds();
		double score = basemembership(tail, head, class1, class2, tatt, hatt);
		if (useheight(class1, class2)) {
			double bl1 = getbl(box1, class1);
			double bl2 = getbl(box2, class2);
			if (bl1 > 0 && bl2 > 0) {
				double free = std::max(0.25*box1.height(), 0.125*TABLETPC_DPI);
				double limit = std::max(1.5*box1.height(), 0.125*TABLETPC_DPI);
				double diff = std::min(limit+free, std::max(0.0, bl1 - bl2 + free));
				double blfac = std::pow(diff/(limit+free), 1.0/4);
				VERBOSE(*verb_out << "baseline " << bl1 << ',' << bl2 << " score " << diff/limit << " with diff " << diff << " and limit " << limit << " : " << blfac << std::endl);
				score *= blfac;
			}
			score *= angle_height_score(std::max(1l, box1.height()), std::max(1l, box2.height()), class1, class2);
		}
		VERBOSE(*verb_out << " final score: " << score << std::endl);
		return score;
	}
};

struct R_rel : public fuzzy_relation {
	R_rel() : fuzzy_relation("R", -60, 0, 60, &get_angle_x) { }
	double membership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const {
		Rect<long> box1 = tail->bounds();
		Rect<long> box2 = head->bounds();
		double score = basemembership(tail, head, class1, class2, tatt, hatt);
		if (useheight(class1, class2)) {
			double bl1 = getbl(box1, class1);
			double bl2 = getbl(box2, class2);
			if (bl1 > 0 && bl2 > 0) {
				double limit = std::max((double)std::max(box1.height(), box2.height()), 0.125*TABLETPC_DPI);
				double diff = std::min(limit, std::abs(bl1 - bl2));
				double blfac = std::pow(1-diff/limit, 1.0/4);
				VERBOSE(*verb_out << "baselines " << bl1 << ',' << bl2 << " score " << 1-diff/limit << " with diff " << diff << " and limit " << limit << " : " << blfac << std::endl);
				score *= blfac;
			}
			score *= right_height_score(std::max(1l, box1.height()), std::max(1l, box2.height()), class1, class2);
		}
		VERBOSE(*verb_out << " final score: " << score << std::endl);
		return score;
	}
};

struct C_rel : public fuzzy_relation {
	C_rel() : fuzzy_relation("C", 0, 0, 0, 0) { }

	double membership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const
	{
		Rect<long> box1 = tail->bounds();
		Rect<long> box2 = head->bounds();
		double overlap = intersection(box1, box2).area();
		if (box2.area() != 0) {
			overlap /= box2.area();
		}
		VERBOSE(*verb_out << "measuring C for boxes " << box1 << " and " << box2 << "; classes " << class1 << " and " << class2 << "; overlap " << overlap << std::endl);
		return overlap;
	}
};





const unsigned ndims = 7;

//static margintree<ndims,3> *mt = 0;
//static margintree<ndims,3> *grpmt = 0;
typedef std::pair<rclass_t, rclass_t> rclass_pair_t;
typedef std::map<rclass_pair_t, normaldist> distmap_t;
//typedef std::map<rclass_pair_t, double> priormap_t;
typedef std::map<unsigned, double> relpriormap_t;
typedef std::map<rclass_pair_t, relpriormap_t> priormap_t;
static distmap_t dists[scg::NUM_RELATIONS+1][ndims];
static double prior[scg::NUM_RELATIONS+1];
static priormap_t priors;//[scg::NUM_RELATIONS+1];

static const int AGGREGATE = NUM_RELATIONS;

static const double DIST_THRES = 0.07;

static void
extract_features(Rect<long> tailbox, Rect<long> headbox, double *f) {
	double xs = std::max((double)tailbox.width(), scg::TABLETPC_DPI / 16.0);
	double ys = std::max((double)tailbox.height(), scg::TABLETPC_DPI / 16.0);
	double xshead = std::max((double)headbox.width(), scg::TABLETPC_DPI / 16.0);
	double yshead = std::max((double)headbox.height(), scg::TABLETPC_DPI / 16.0);
	//xs = ys = TABLETPC_DPI;
	double scale = std::max(std::max(xs, ys), std::max(xshead, yshead));
	//double scale = (std::max(xs, ys) + std::max(xshead, yshead)) / 2.0;
	//double scale = std::max(xs, ys);
	static const double RELATIVE_MIX = 0.5;
	//xs = tailbox.width() + TABLETPC_DPI;
	//ys = tailbox.height() + TABLETPC_DPI;
	xs = ys = RELATIVE_MIX * scale + (1.0 - RELATIVE_MIX) * TABLETPC_DPI;
	f[0] = (headbox.left - tailbox.left) / xs;
	f[1] = (headbox.right - tailbox.right) / xs;
	f[2] = (headbox.left - tailbox.right) / xs;
	f[3] = (headbox.top - tailbox.top) / ys;
	f[4] = (headbox.bottom - tailbox.bottom) / ys;
	f[5] = (headbox.top - tailbox.bottom) / ys;
	f[6] = overlap_proportion(tailbox, headbox);
	VERBOSE(
		*verb_out << "F ( " << xs << ',' << ys << ") = [ ";
		for (size_t i = 0; i < ndims; ++i) *verb_out << f[i] << ' ';
		*verb_out << "]\n";
	);
}

struct dist_relation : public relation {
	size_t reli;
	std::vector<size_t> dims;

	dist_relation(const std::string &name_, size_t reli_)
		: relation(name_), reli(reli_) { }

	bool data_ok(rclass_pair_t clspair) const {
		for (std::vector<size_t>::const_iterator i = dims.begin(); i != dims.end(); ++i) {
			size_t dimi = *i;
			std::map<rclass_pair_t, normaldist>::const_iterator j = dists[reli][dimi].find(clspair);
			if (j == dists[reli][dimi].end() || !j->second.useful(DIST_THRES)) return false;
			j = dists[NUM_RELATIONS][dimi].find(clspair);
			if (j == dists[NUM_RELATIONS][dimi].end() || !j->second.useful(DIST_THRES)) return false;
		}
		return true;
	}

	double dimscore(double *f, rclass_pair_t clspair) const {
		static const rclass_pair_t GENERIC_PAIR(AGGREGATE_CLASS, AGGREGATE_CLASS);
		double p = 1.0;
		for (std::vector<size_t>::const_iterator i = dims.begin(); i != dims.end(); ++i) {
			size_t dimi = *i;
			const normaldist &clsdist = dists[reli][dimi][clspair];
			if (clspair != GENERIC_PAIR && !clsdist.useful(DIST_THRES)) {
				VERBOSE(*verb_out << " dist for rel " << reli << " in dim " << dimi << " for classes " << clspair.first << ',' << clspair.second << " not useful; canceling\n");
				return -1.0;
			}
			double clsp = std::max(1e-6, clsdist.eval(f[dimi]));
			p *= clsp;
			VERBOSE(
				*verb_out << "  dim " << dimi << ": " << f[dimi]
						  << "  and rel-dist " << clsdist << " -> " << clsp
				          << "  -> current p " << p << std::endl;
			);
		}
		return p;
	}

	double membership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const {
		VERBOSE(*verb_out << "testing relation " << name << " between " << tail->bounds() << " and " << head->bounds() << " with classes " << class1 << " and " << class2 << " : ");
		double f[ndims];
		extract_features(tail->bounds(), head->bounds(), f);
		double sc = dimscore(f, std::make_pair(class1, class2));
		VERBOSE(*verb_out << " -> score for " << name << " = " << sc << std::endl);
		return sc;
	}
};

double
get_nil_probability(const interpretation *tail, const interpretation *head) {
	Rect<long> tailbox = tail->bounds();
	Rect<long> headbox = head->bounds();

	VERBOSE(*verb_out << "nil probability between " << tailbox << " and " << headbox << std::endl);

	static const rclass_pair_t clspair(AGGREGATE_CLASS, AGGREGATE_CLASS);

	double f[ndims];
	extract_features(tailbox, headbox, f);

	double sum = 0.0;
	for (size_t reli = 0; reli < NUM_RELATIONS; ++reli) {
		dist_relation *R = (dist_relation *)get_relation(reli);
		sum += std::min(100.0, R->dimscore(f, clspair));
		VERBOSE(*verb_out << " sum -> " << sum << std::endl);
	}
	double pnil = 1.0 - sum / (sum + 1.0);
	VERBOSE(*verb_out << " -> pnil " << pnil << std::endl);
	return pnil;
}

static unsigned
link_name_to_index(const std::string &s) {
	if (s == "AR") {
		return 0;
	}
	if (s == "R") {
		return 1;
	}
	if (s == "BR") {
		return 2;
	}
	if (s == "B" || s == "BW" || s == "BN") {
		return 3;
	}
	if (s == "C" || s == "Cs") {
		return 4;
	}
	if (s == "nil") {
		return 5;
	}

	THROW_ERROR(E_INVALID, "unknown link name " << s);
}


int
initialize_relations() {
	std::string path;
	int e = GetProfilePath(path);
	if (FAILURE(e)) {
		return e;
	}
	path += "/relations.txt";
	std::ifstream in(path.c_str());
	if (!in.is_open()) {
		return E_NOTFOUND;
	}

	for (size_t i = 0; i < NUM_RELATIONS; ++i) {
		in >> prior[i];
	}
	if (!in) return E_IO;

	for (;;) {
		std::string type;
		in >> type;
		if (!in) break;
		if (type == "prior") {
			unsigned reli;
			rclass_pair_t clspair;
			in >> reli >> clspair.first >> clspair.second;
			in >> priors[clspair][reli];
		}
		else if (type == "dist") {
			unsigned reli;
			unsigned dimi;
			rclass_pair_t clspair;
			in >> reli >> dimi >> clspair.first >> clspair.second;
			in >> dists[reli][dimi][clspair];
			dists[reli][dimi][clspair].var *= 2.0;
		}
		else {
			return E_INVALID;
		}
		if (!in) return E_IO;
	}

	/*
	const char *names[NUM_RELATIONS] = {"AR", "R", "BR", "B", "C"};
	for (size_t i = 0; i < NUM_RELATIONS; ++i) {
		if (!relations[i]) {
			relations[i] = new dist_relation(names[i], i);
		}
	}*/

	if (!relations[REL_ABOVERIGHT]) {
		dist_relation *AR = new dist_relation("AR", REL_ABOVERIGHT);
		relations[REL_ABOVERIGHT] = AR;
		AR->dims.push_back(2); AR->dims.push_back(3); AR->dims.push_back(4); //AR->dims.push_back(6);
	}
	if (!relations[REL_RIGHT]) {
		dist_relation *R = new dist_relation("R", REL_RIGHT);
		relations[REL_RIGHT] = R;
		R->dims.push_back(2); R->dims.push_back(3); R->dims.push_back(4); //R->dims.push_back(6);
	}
	if (!relations[REL_BELOWRIGHT]) {
		dist_relation *BR = new dist_relation("BR", REL_BELOWRIGHT);
		relations[REL_BELOWRIGHT] = BR;
		BR->dims.push_back(2); BR->dims.push_back(3); BR->dims.push_back(4); //BR->dims.push_back(6);
	}
	if (!relations[REL_BELOW]) {
		dist_relation *B = new dist_relation("B", REL_BELOW);
		relations[REL_BELOW] = B;
		B->dims.push_back(0); B->dims.push_back(1); B->dims.push_back(5); //B->dims.push_back(6);
	}
	if (!relations[REL_CONTAINS]) {
		dist_relation *C = new dist_relation("C", REL_CONTAINS);
		relations[REL_CONTAINS] = C;
		C->dims.push_back(0); C->dims.push_back(1); C->dims.push_back(3); C->dims.push_back(4); C->dims.push_back(6);
	}

	if (!fuzrelations[REL_ABOVERIGHT]) fuzrelations[REL_ABOVERIGHT] = new AR_rel();
	if (!fuzrelations[REL_RIGHT]) fuzrelations[REL_RIGHT] = new R_rel();
	if (!fuzrelations[REL_BELOWRIGHT]) fuzrelations[REL_BELOWRIGHT] = new BR_rel();
	if (!fuzrelations[REL_BELOW]) fuzrelations[REL_BELOW] = new B_rel();
	if (!fuzrelations[REL_CONTAINS]) fuzrelations[REL_CONTAINS] = new C_rel();

	return 0;
}


void
destroy_relations() {
	for (unsigned i = 0; i < NUM_RELATIONS; i++) {
		delete relations[i];
		relations[i] = 0;
		delete fuzrelations[i];
		fuzrelations[i] = 0;
	}
}


unsigned
ordering_for_relation(const relation *rel) {
	return rel_orderings[relindex(rel)];
}

size_t
relindex(const relation *rel) {
	relation **rp = std::find(relations, relations + NUM_RELATIONS, rel);
	assert(rp != relations + NUM_RELATIONS);
	return rp - relations;
}

relation *
get_relation(int rel) {
	return relations[rel];
}

relation *
get_group_relation(int rel) {
	return fuzrelations[rel];
}


void
set_relation(int rel, relation *r) {
	relations[rel] = r;
}


}
