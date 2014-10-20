#include "relation.h"
#include "links.h"
#include "meas.h"
#include "utils.h"
#include "grammar.h"
#include "expr-node.h"
#include "segment.h"
#include "rect.h"
#include "stream-defs.h"

#include <fstream>
#include <map>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>


namespace scg
{



static relation *relations[NUM_RELATIONS];

static unsigned rel_orderings[] = { 0, 0, 1, 1, 0 };

static std::string box_name("_BOX");
static std::string symbol_name("_SYMBOL");
static std::string aggregate_name("_aggregate");


static const unsigned NUM_SPECIAL_CLASSES = 3;
const rclass_t INVALID_CLASS = static_cast<rclass_t>(-1);
const rclass_t AGGREGATE_CLASS = 0;
const rclass_t SYMBOL_CLASS = 1;
const rclass_t BOX_CLASS = 2;

// these following classes must match the order given in symbols.cc
const rclass_t BASELINE_CLASS = 3;
const rclass_t ASCENDER_CLASS = 4;
const rclass_t DESCENDER_CLASS = 5;
const rclass_t EXTENDER_CLASS = 6;
const rclass_t HALF_ASCENDER_CLASS = 7;
const rclass_t HALF_DESCENDER_CLASS = 8;
const rclass_t HALF_ASCENDER_DESCENDER_CLASS = 9;
const rclass_t CENTERED_CLASS = 10;
const rclass_t FENCE_CLASS = 11;
const rclass_t LARGE_EXTENDER_CLASS = 12;
const rclass_t RANGE_OP_CLASS = 13;
const rclass_t ROOT_CLASS = 14;
const rclass_t HORIZONTAL_CLASS = 15;
const rclass_t PUNCTUATION_CLASS = 16;


static long
get_centroid_y(const Rect<long> &box, rclass_t rclass)
{
	switch (rclass) {
	case AGGREGATE_CLASS:
	case SYMBOL_CLASS:
	case BOX_CLASS:
	case CENTERED_CLASS:
	case LARGE_EXTENDER_CLASS:
	case RANGE_OP_CLASS:
	case FENCE_CLASS:
	case ROOT_CLASS:
	case HORIZONTAL_CLASS:
	case PUNCTUATION_CLASS:
	case ASCENDER_CLASS:
	case EXTENDER_CLASS:
		return (box.top + box.bottom) / 2;
	case BASELINE_CLASS:
		return box.top + height(box) * 0.1;
	case HALF_DESCENDER_CLASS:
		return box.top + height(box) * 0.075;
	case HALF_ASCENDER_DESCENDER_CLASS:
		return box.top + height(box) * 0.25;
	case DESCENDER_CLASS:
		return box.top + height(box) * 0.05;
	case HALF_ASCENDER_CLASS:
		return box.top + height(box) * 0.33;
	default:
		THROW_ERROR(E_INVALID, "cannot obtain centroid for unknown rel-class " << rclass);
	}
}

static double FreeRelationDist = 0.0;
static void calc_parms() {
	FreeRelationDist = TABLETPC_DPI * GetParameterDouble("FreeRelationDist");
}
static int _e = RegisterParameterCallback(&calc_parms);

static long
inset(long d)
{
	return static_cast<long>(0.25 * std::min<double>(d, FreeRelationDist));
}


static void
get_from_centroid_y(const Rect<long> &box, rclass_t rclass, long &x, long &y)
{
	/*
	if (rclass == BOX_CLASS || rclass == ROOT_CLASS || rclass == HORIZONTAL_CLASS) {
		y = box.bottom - inset(height(box));
	}
	else {*/
		y = (box.top + box.bottom) / 2;
	//}
	x = (box.left + box.right) / 2;
}

static void
get_to_centroid_y(const Rect<long> &box, rclass_t rclass, long &x, long &y)
{
	/*
	if (rclass == BOX_CLASS || rclass == ROOT_CLASS || rclass == HORIZONTAL_CLASS) {
		y = box.top + inset(height(box));
	}
	else {*/
		y = (box.top + box.bottom) / 2;
	//}
	x = (box.left + box.right) / 2;
}


static void
get_from_centroid_x(const Rect<long> &box, rclass_t rclass, long &x, long &y)
{
	/*
	if (rclass == BOX_CLASS || rclass == ROOT_CLASS || rclass == HORIZONTAL_CLASS) {
		x = box.right - inset(width(box));
	}
	else {*/
		x = (box.left + box.right) / 2;
	//}
	y = get_centroid_y(box, rclass);
}

static void
get_to_centroid_x(const Rect<long> &box, rclass_t rclass, long &x, long &y)
{
	/*
	if (rclass == BOX_CLASS || rclass == ROOT_CLASS || rclass == HORIZONTAL_CLASS) {
		x = box.left + inset(width(box));
	}
	else {*/
		x = (box.left + box.right) / 2;
	//}
	y = get_centroid_y(box, rclass);
}


double
get_angle_y(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2)
{
	//return static_cast<double>(box2.top - box1.bottom) / std::max<double>(height(box1), FreeRelationDist);
	long x1,y1;
	long x2,y2;
	get_from_centroid_y(box1, class1, x1, y1);
	get_to_centroid_y(box2, class2, x2, y2);
	VERBOSE(*verb_out << "angle from (" << x1 << "," << y1 << ")->(" << x2 << "," << y2 << ") = " << std::atan2((double)y2-y1, x2-x1) << std::endl);
	return std::atan2((double)y2-y1, x2-x1);
}

double
get_angle_x(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2)
{
	long x1,y1;
	long x2,y2;
	get_from_centroid_x(box1, class1, x1, y1);
	get_to_centroid_x(box2, class2, x2, y2);
	//return static_cast<double>(y2-y1) / std::max<double>(height(box1), FreeRelationDist);//height(box2));
	VERBOSE(*verb_out << "angle from (" << x1 << "," << y1 << ")->(" << x2 << "," << y2 << ") = " << std::atan2((double)y2-y1, x2-x1) << std::endl);
	return std::atan2((double)y2-y1, x2-x1);
}

double
get_AR_angle(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2)
{
	//return static_cast<double>(box1.bottom - box2.bottom) / std::max<double>(height(box1), FreeRelationDist);//height(box2));
	//long x1 = box1.right - inset(width(box1));
	//long y1 = box1.top + inset(height(box1));
	long x1, y1;
	get_from_centroid_x(box1, class1, x1, y1);
	long x2 = box2.left + inset(width(box2));
	long y2 = box2.bottom - inset(height(box2));
	//long x2, y2;
	//get_to_centroid_x(box2, class2, x2, y2);
	return std::atan2((double)y2-y1, x2-x1);
}

double
get_BR_angle(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2)
{
	//return static_cast<double>(box2.top - box1.top) / std::max<double>(height(box1), FreeRelationDist);//std::max(height(box1), height(box2));
	//long x1 = box1.right - inset(width(box1));
	//long y1 = box1.bottom - inset(height(box1));
	long x1, y1;
	get_from_centroid_x(box1, class1, x1, y1);
	long x2 = box2.left + inset(width(box2));
	long y2 = box2.top + inset(height(box2));
	//long x2, y2;
	//get_to_centroid_x(box2, class2, x2, y2);
	return std::atan2((double)y2-y1, x2-x1);
}


static long
distance_threshold(const Rect<long> &box1, const Rect<long> &box2)
{
	long b1thresh = (width(box1) + height(box1))/2;
	long b2thresh = (width(box2) + height(box2))/2;
	long mx_thresh = std::max(b1thresh, b2thresh);
	long mn_thresh = 2 * std::min(b1thresh, b2thresh);
	return std::min<double>(FreeRelationDist, std::max<long>(FreeRelationDist / 2, std::min(mx_thresh, mn_thresh)));
}

static double
distance_score(const Rect<long> &box1, const Rect<long> &box2)
{
	long thresh = distance_threshold(box1, box2);
	double dist = std::sqrt((double)dist_sq(box1, box2));
	if (dist < thresh) {
		return 1.0;
	}
	else {
		return std::max(0.0, 1.0 - (dist - thresh) / (2 * thresh));
	}
}

relation::relation(const std::string &name_, double t0_, double t1_, double t2_, get_angle_fn t_)
	: name(name_), t0(t0_), t1(t1_), t2(t2_), t(t_)
	{ }

double
relation::membership(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2) const
{
	double ct = t(box1, box2, class1, class2);
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

	double dist_score = distance_score(box1, box2);
	double overlap = overlap_proportion(box1, box2);
	VERBOSE(*verb_out << "measuring " << name << " for boxes " << box1 << " and " << box2 << "; classes " << class1 << " and " << class2 << "; ");
	VERBOSE(*verb_out << " scores " << angle_score << " * " << dist_score << " = " << angle_score * dist_score << std::endl);

	/*
	if ((class1 == PUNCTUATION_CLASS && class2 != ROOT_CLASS) || (class2 == PUNCTUATION_CLASS && class1 != ROOT_CLASS)) {
		overlap = 0.0;
	}*/

	return angle_score * dist_score * (1 - 0.5*overlap);
}


/*
struct AR_rel : public relation {
	AR_rel() : relation("AR") { }

	double membership(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2) const
	{
		double ct = get_angle(box1, box2, class1, class2);
		ct *= 180.0/M_PI;

		double angle_score;

		VERBOSE(*verb_out << "measuring AR for boxes " << box1 << " and " << box2 << "; classes " << class1 << " and " << class2 << "; angle: " << ct << std::endl);

		if (ct < -80.0) {
			return 0.0;
		}
		else if (ct < -40.0) {
			angle_score = 1.0 + (40.0 + ct) / 40.0;
		}
		else if (ct < 0.0) {
			angle_score = -ct / 40.0;
		}
		else {
			return 0.0;
		}

		double dist_score = distance_score(box1, box2);
		double overlap = overlap_proportion(box1, box2);
		VERBOSE(*verb_out << " scores " << angle_score << " * " << dist_score << " = " << angle_score * dist_score << std::endl);
		return angle_score * dist_score * (1 - overlap);
	}
};

struct R_rel : public relation {
	R_rel() : relation("R") { }

	double membership(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2) const
	{
		double ct = get_angle(box1, box2, class1, class2);
		ct *= 180.0/M_PI;

		double angle_score;

		VERBOSE(*verb_out << "measuring R for boxes " << box1 << " and " << box2 << "; classes " << class1 << " and " << class2 << "; angle: " << ct << std::endl);
		if (ct < -40.0) {
			return 0.0;
		}
		else if (ct < 0.0) {
			angle_score = 1.0 + ct / 40.0;
		}
		else if (ct < 40.0) {
			angle_score = 1.0 - ct / 40.0;
		}
		else {
			return 0.0;
		}
 
		double dist_score = distance_score(box1, box2);
		double overlap = overlap_proportion(box1, box2);
		VERBOSE(*verb_out << " scores " << angle_score << " * " << dist_score << " = " << angle_score * dist_score << std::endl);
		return angle_score * dist_score * (1 - overlap);
	}
};

struct BR_rel : public relation {
	BR_rel() : relation("BR") { }

	double membership(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2) const
	{
		double ct = get_angle(box1, box2, class1, class2);
		ct *= 180.0/M_PI;

		double angle_score;

		VERBOSE(*verb_out << "measuring BR for boxes " << box1 << " and " << box2 << "; classes " << class1 << " and " << class2 << "; angle: " << ct << std::endl);

		if (ct < 0.0) {
			return 0.0;
		}
		else if (ct < 35.0) {
			angle_score = ct / 35.0;
		}
		else if (ct < 80.0) {
			angle_score = 1.0 - (ct - 35.0) / 45.0;
		}
		else {
			return 0.0;
		}

		double dist_score = distance_score(box1, box2);
		double overlap = overlap_proportion(box1, box2);
		VERBOSE(*verb_out << " scores " << angle_score << " * " << dist_score << " = " << angle_score * dist_score << std::endl);
		return angle_score * dist_score * (1 - overlap);
	}
};

struct B_rel : public relation {
	B_rel() : relation("B") { }

	double membership(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2) const
	{
		double ct = get_angle(box1, box2, class1, class2);
		ct *= 180.0/M_PI;

		double angle_score;

		VERBOSE(*verb_out << "measuring B for boxes " << box1 << " and " << box2 << "; classes " << class1 << " and " << class2 << "; angle: " << ct << std::endl);

		if (ct < 0.0) {	
			return 0.0;
		}
		else if (ct < 90.0) {
			angle_score = ct / 90.0;
		}
		else {
			angle_score = (180.0 - ct) / 90.0;
		}

		double dist_score = distance_score(box1, box2);
		double overlap = overlap_proportion(box1, box2);
		VERBOSE(*verb_out << " scores " << angle_score << " * " << dist_score << " = " << angle_score * dist_score << std::endl);
		return angle_score * dist_score * (1 - overlap);
	}
};
*/

struct C_rel : public relation {
	C_rel() : relation("C", 0, 0, 0, 0) { }

	double membership(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2) const
	{
		double overlap = overlap_proportion(box1, box2);
		VERBOSE(*verb_out << "measuring C for boxes " << box1 << " and " << box2 << "; classes " << class1 << " and " << class2 << "; overlap " << overlap << std::endl);
		return overlap;
	}
};

static unsigned
link_name_to_index(const std::string &s)
{
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


void
initialize_relations()
{
	std::string path;
	GetTrainingPath(path);
	std::ifstream relin((path + "/reldata").c_str());
	if (!relin.is_open()) {
		THROW_ERROR(E_NOTFOUND, "could not open reldata file");
	}

	static std::map<std::string, get_angle_fn> angle_fns;
	if (angle_fns.empty()) {
		angle_fns["AR"] = &get_angle_x;//&get_AR_angle;
		angle_fns["R"] = &get_angle_x;
		angle_fns["BR"] = &get_angle_x;//&get_BR_angle;
		angle_fns["B"] = &get_angle_y;
	}
	for (unsigned i = 0; i < NUM_RELATIONS; ++i) {
		std::string name;
		double t0, t1, t2;
		relin >> name;
		CHECK_ISTREAM(relin, "error reading relation data");
		if (name == "C") {
			relations[i] = new C_rel();
		}
		else {
			relin >> t0 >> t1 >> t2;
			CHECK_ISTREAM(relin, "error reading relation data for " << name);
			relations[i] = new relation(name, t0, t1, t2, angle_fns[name]);
		}
	}
}


void
destroy_relations()
{
	for (unsigned i = 0; i < NUM_RELATIONS; i++) {
		delete relations[i];
	}
}


unsigned
ordering_for_relation(const relation *rel) {
	relation **rp = std::find(relations, relations + NUM_RELATIONS, rel);
	assert(rp != relations + NUM_RELATIONS);
	return rel_orderings[rp - relations];
}


relation *
get_relation(int rel) {
	return relations[rel];
}


}
