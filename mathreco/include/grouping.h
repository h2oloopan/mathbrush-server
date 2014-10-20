#ifndef GROUPING_H_
#define GROUPING_H_


#include "group.h"
#include "stroke.h"
#include "normalize.h"
#include "bitvec.h"
#include "reco-types.h"
//#include "interp.h"
//#include "stroke-alg.h"

#include <list>
#include <map>


namespace scg
{

struct context;
struct group;

void calc_grouping_metrics();

//int update_grouping(reco_context *ctx);

//int update_proximity_groups(context *ctx, std::list<group *> &newgroups);
//int update_stack_groups(context *ctx, std::list<group *> &cmpgroups, std::list<group *> &newgroups);
void dilute_group_by_dist(group *grp, const std::map<bitvec, group *> &groups);


RawStroke strip_duplicate_points(const RawStroke &s);

template <typename T>
static RawStrokeGroup
build_stroke_group(T first, T last, const Rect<long> &bounds, size_t nstrokes)
{
	RawStroke *strokes = new RawStroke[nstrokes];
	RawStroke *ps = strokes;
	while (first != last) {
		RawStroke &s = *(ps++);
		RawStroke uniq = strip_duplicate_points((*(first++))->stk);
		//NormalizedStroke ns = normalize(uniq, bounds);
		//s = ns;
		s = uniq;
		//s = smooth(trim_ends(s));
		//s = subdivide(s);
	}

	return RawStrokeGroup(strokes, nstrokes);
}


}


#endif
