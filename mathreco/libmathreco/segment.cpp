#include "segment.h"
#include "mathrecognizer-private.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <vector>

#include "dist.h"
#include "error.h"
#include "group.h"
#include "grouping.h"
#include "parms.h"
#include "rect.h"
#include "stroke.h"
//#include "stroke-alg.h"
#include "verb.h"
#include "set-utils.h"


namespace scg
{



long TABLETPC_DPI = 2540; // This is hard-coded to match Microsoft's Tablet API specs


static long MinSplitDistance = 0;
static long MaxMergeDistance = 0;
static long SplitDistance = 0;
static long SplitDistanceSq = 0;
static long MaxDotSize = 0;

const double dot_distance_meas::DY_WEIGHT = 0.5;

void calc_segment_metrics() {
	MinSplitDistance = static_cast<long>(TABLETPC_DPI * GetParameterDouble("MinSplitDistance"));
	MaxMergeDistance = static_cast<long>(TABLETPC_DPI * GetParameterDouble("MaxMergeDistance"));
	SplitDistance = static_cast<long>(TABLETPC_DPI * GetParameterDouble("SplitDistance"));
	SplitDistanceSq = SplitDistance * SplitDistance;
	MaxDotSize = static_cast<long>(TABLETPC_DPI * GetParameterDouble("MaxDotSize"));
}
static int _e = RegisterParameterCallback(&calc_segment_metrics);

DLLDECL void
SetTabletResolution(size_t dpi)
{
	TABLETPC_DPI = dpi;
	calc_segment_metrics();
	calc_grouping_metrics();
}


static long
distance_threshold(const Rect<long> &bbox1, const Rect<long> &bbox2)
{
    const static double DistanceThresholdRatio = GetParameterDouble("DistanceThresholdRatio");
    //return static_cast<long>(std::sqrt(std::max<double>(width(bbox1) * height(bbox1), width(bbox2) * height(bbox2))) * 0.25);
	return std::max(SplitDistance, (long)(DistanceThresholdRatio * std::max(std::max(bbox1.width(), bbox1.height()), std::max(bbox2.width(), bbox2.height()))));
}

double
dot_distance_meas::operator()(const segment *seg) const
{
	static const double DY_WEIGHT = 0.25;
	double dx = std::abs(bounds.x_center() - seg->bounds.x_center());
	double dy = std::abs(bounds.y_center() - seg->bounds.y_center());
	return dx + DY_WEIGHT * dy;
}

double
stroke_distance_meas::operator()(const segment *seg) const {
	return dist_sq(stk, seg->stk);
}


template <typename U>
void
shift(U first, U last, int amount)
{
	while (last != first) {
		--last;
		*(last + amount) = *(last);
	}
}


static segment *
merge_segments(segment *seg1, segment *seg2)
{
	static const double ANGLE_DIFFERENCE_THRES = M_PI/8.0;
	//static const long MergeDist = MaxDotSize / 2;

	if (seg1->stk.npoints == 0 || seg2->stk.npoints == 0) {
		return 0;
	}

	segment *merged = 0;

	const RawStroke *s1 = &seg1->stk;
	const RawStroke *s2 = &seg2->stk;
	size_t last_s1_pt = s1->npoints - 1;
	size_t last_s2_pt = s2->npoints - 1;

	// kludgy checks for stitching together strokes caused by tiny jitters
	//long dss = dist_sq(s1->x[0], s1->y[0], s2->x[0], s2->y[0]);
	//long dse = dist_sq(s1->x[0], s1->y[0], s2->x[last_s2_pt], s2->y[last_s2_pt]);
	long des = dist_sq(s1->x[last_s1_pt], s1->y[last_s1_pt], s2->x[0], s2->y[0]);
	//long dee = dist_sq(s1->x[last_s1_pt], s1->y[last_s1_pt], s2->x[last_s2_pt], s2->y[last_s2_pt]);
	/*if (dss < MaxMergeDist) {
		long d = dss;
		
		for (size_t i = 1, j = 1; i < num_points(*s1) && j < num_points(*s2) && d <= dss; ++i, ++j) {
			dss = d;
			d = dist_sq(s1->x[i], s1->y[i], s2->x[j], s2->y[j]);
		}
		if (i == num_points(*s1) && j == num_points(*s2)) {
		}
		else if (i == num_points(*s1)) {
			merged = new segment;
			merged->stk = copy(*s2);
			merged->bounds = seg2->bounds;
			merged->children = seg1->children;
			merged->children.insert(merged->children.end(), seg1->children.begin(), seg1->children.end());
			VERBOSE(*verb_out << "  merging head/head based on prox; cutting s1\n");
		}
		else if (j == num_points(*s2)) {
			merged = new segment;
			merged->stk = copy(*s1);
			merged->bounds = seg1->bounds;
			merged->children = seg1->children;
			merged->children.insert(merged->children.end(), seg2->children.begin(), seg2->children.end());
			VERBOSE(*verb_out << "  merging head/head based on prox; cutting s2\n");
		}
		else {
			--i, --j;
			merged = new segment;
			merged->stk = copy(*
		}
	}
	else if (dse < MaxMergeDist) {
	}
	else*/ if (des < MaxMergeDistance*MaxMergeDistance) {
		long d = des;
		size_t i,j;
		for (i = last_s1_pt - 1, j = 1; i > 0 && j < last_s2_pt && d <= des; --i, ++j) {
			des = d;
			d = dist_sq(s1->x[i], s1->y[i], s2->x[j], s2->y[j]);
		}
		if (i == 0 && j == s2->npoints) {
			merged = new segment;
			//merged->ctx = seg1->ctx;
			RawStroke rev1 = s1->copy();
			RawStroke rev2 = s2->copy();
			RawStroke cat = concat(rev1.reverse(), rev2.reverse());
			merged->stk = cat;
			merged->bounds = merge(seg1->bounds, seg2->bounds);
			merged->children = seg1->children;
			merged->children.insert(merged->children.end(), seg2->children.begin(), seg2->children.end());
			VERBOSE(*verb_out << "  merging head/head based on prox");
		}
		else if (i == 0) {
			merged = new segment;
			//merged->ctx = seg1->ctx;
			RawStroke rev = s1->copy();
			RawStroke cat = concat(rev.reverse(), *s2);
			merged->stk = cat;
			merged->bounds = merge(seg1->bounds, seg2->bounds);
			merged->children = seg1->children;
			merged->children.insert(merged->children.end(), seg2->children.begin(), seg2->children.end());
			VERBOSE(*verb_out << "  merging head/head based on prox");
		}
		else if (j == last_s2_pt) {
			merged = new segment;
			//merged->ctx = seg1->ctx;
			RawStroke rev = s2->copy();
			RawStroke cat = concat(*s1, rev.reverse());
			merged->stk = cat;
			merged->bounds = merge(seg1->bounds, seg2->bounds);
			merged->children = seg1->children;
			merged->children.insert(merged->children.end(), seg2->children.begin(), seg2->children.end());
			VERBOSE(*verb_out << "  merging tail/tail based on prox\n");
		}
		else {
			merged = new segment;
			//merged->ctx = seg1->ctx;
			RawStroke sub1 = s1->substroke(0, i+2);
			RawStroke sub2 = s2->substroke(j-2);
			RawStroke cat = concat(sub1, sub2);
			merged->stk = cat;
			merged->bounds = merged->stk.bounds();
			merged->children = seg1->children;
			merged->children.insert(merged->children.end(), seg2->children.begin(), seg2->children.end());
			VERBOSE(*verb_out << "  merging tail/head based on prox\n");

		}
	}
	/*else if (dee < MaxMergeDist) {
	}*/
	else if (std::abs(s1->x[0] - s2->x[last_s2_pt]) < MaxMergeDistance
	 && std::abs(s1->y[0] - s2->y[last_s2_pt]) < MaxMergeDistance) {
		double dt = std::atan2(static_cast<double>(s2->y[1] - s2->y[0]), s2->x[1] - s2->x[0])
		          - std::atan2(static_cast<double>(s1->y[last_s1_pt] - s1->y[last_s1_pt-1]), s1->x[last_s1_pt] - s1->x[last_s1_pt-1]);
		if (std::abs(dt) < ANGLE_DIFFERENCE_THRES) {
			merged = new segment;
			//merged->ctx = seg1->ctx;
			RawStroke cat = concat(*s2, *s1);
			merged->stk = cat;
			merged->bounds = merge(seg1->bounds, seg2->bounds);
			merged->children = seg2->children;
			merged->children.insert(merged->children.end(), seg1->children.begin(), seg1->children.end());
			VERBOSE(*verb_out << "  merging based on endpoint match\n");
		}
	}
	else if (std::abs(s2->x[0] - s1->x[last_s1_pt]) < MaxMergeDistance
	 && std::abs(s2->y[0] - s1->y[last_s1_pt]) < MaxMergeDistance) {
		double dt = std::atan2(static_cast<double>(s2->y[last_s2_pt] - s2->y[last_s2_pt-1]), s2->x[last_s2_pt] - s2->x[last_s2_pt-1])
		          - std::atan2(static_cast<double>(s1->y[1] - s1->y[0]), s1->x[1] - s1->x[0]);
		if (std::abs(dt) < ANGLE_DIFFERENCE_THRES) {
			merged = new segment;
			//merged->ctx = seg1->ctx;
			RawStroke cat = concat(*s1, *s2);
			merged->stk = cat;
			merged->bounds = merge(seg1->bounds, seg2->bounds);
			merged->children = seg1->children;
			merged->children.insert(merged->children.end(), seg2->children.begin(), seg2->children.end());
			VERBOSE(*verb_out << "  merging based on endpoint match\n");
		}
	}


	// check for horizontal line extensions.
	// TODO: this algorithm should be generalized to work in arbitrary orientations and then
	// it can replace the shoddy implmentation above
	else {
		if (s2->x[0] < s1->x[0]) {
			std::swap(s1, s2);
			std::swap(last_s1_pt, last_s2_pt);
		}
		
		/*double t1 = std::atan2((double)s1->y[last_s1_pt] - s1->y[0], s1->x[last_s1_pt] - s1->x[0]);
		if (std::abs(t1) > M_PI/8) {
			return 0;
		}
		double t2 = std::atan2((double)s2->y[last_s2_pt] - s2->y[0], s2->x[last_s2_pt] - s2->x[0]);
		if (std::abs(t2) > M_PI/8) {
			return 0;
		}*/

		if (dist_sq(s2->x[0], s2->y[0], *s1) <= SplitDistanceSq
	    && dist_sq(s1->x[last_s1_pt], s1->y[last_s1_pt], *s2) <= SplitDistanceSq) {
			long *x1 = s1->x;
			long *x2 = s2->x;
			long *y1 = s1->y;
			long *y2 = s2->y;
			long *endx1 = x1 + s1->npoints;
			long *endx2 = x2 + s2->npoints;
			unsigned split = 0;

			while (x1 < endx1 && *x2 > *x1) {
				++x1;
				++y1;
			}

			if (x1 == endx1) {
				return 0;
			}

			while (x1 < endx1) {
				while (x2 < endx2 && *x1 > *x2) {
					long dy = std::abs(*y2 - *y1);
					if (dy > MaxMergeDistance) {
						return 0;
					}

					++x2;
					++y2;
				}

				++x1;
				++y1;
			}

			merged = new segment;
			//merged->ctx = seg1->ctx;
			RawStroke sub2 = s2->substroke(x2 - s2->x);
			RawStroke cat = concat(*s1, sub2);
			merged->stk = cat;
			merged->bounds = merged->stk.bounds();
			merged->children = seg1->children;
			merged->children.insert(merged->children.end(), seg2->children.begin(), seg2->children.end());
			VERBOSE(*verb_out << "  merging strokes with bounds " << seg1->bounds << " and " << seg2->bounds << " based on horzline extension\n");
		}
	}

	if (merged) {
		for (std::list<stroke *>::iterator pstk = merged->children.begin(); pstk != merged->children.end(); ++pstk) {
			(*pstk)->parent = merged;
		}
	}

	return merged;
}


#pragma optimize("g", off)

segment *
math_recognizer_base::merge_segment(segment *seg)
{
	unsigned n = 0;
	bool change;
	do {
		change = false;
		//for (std::vector<segment *>::iterator pseg = segments.begin() + n; pseg != segments.end(); ++pseg) {
		for (unsigned i = n; i != segments.size(); ++i) {
			segment *mseg = segments[i];
			//VERBOSE(*verb_out << "segment is " << seg << " and pseg is " << std::endl);
			segment *merged = merge_segments(mseg, seg);
			if (merged) {
				remove_segment(i, false);
				//pseg = ctx->segments.erase(pseg);
				delete seg;
				seg = merged;
				change = true;
				break;
			}
			++n;
		}
	} while (change);

	return seg;
}

#pragma optimize("g", off)


}
