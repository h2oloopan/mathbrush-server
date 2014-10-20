#include "grouping.h"
#include "stkutils.h"
#include "parms.h"
#include "segment.h"
#include "set-utils.h"
#include "verb.h"
#include "mathrecognizer-private.h"
#include "normalize.h"
#include "interp.h"

#include <algorithm>
#include <vector>


namespace scg
{


const static double GROUP_SCORE_THRESHOLD = 0.06;//3;

static size_t MaxGroupSize = (size_t)RegisterParameterUnsigned("MaxGroupSize", reinterpret_cast<unsigned *>(&MaxGroupSize));
static double SplitDistance = 0;
static double SplitDistanceSq = 0;
static double MinSplitDistance = 0;
static double DistanceThresholdRatio = RegisterParameterDouble("DistanceThresholdRatio", &DistanceThresholdRatio);
static double StackOverlapRequiredRatio = RegisterParameterDouble("StackOverlapRequiredRatio", &StackOverlapRequiredRatio);
static double StackedDotMaxHorzDist = 0;
static double MaxSplitDistance = 0;

void
calc_grouping_metrics()
{
	SplitDistance = TABLETPC_DPI * GetParameterDouble("SplitDistance");
	SplitDistanceSq = SplitDistance * SplitDistance;
	MinSplitDistance = TABLETPC_DPI * GetParameterDouble("MinSplitDistance");
	StackedDotMaxHorzDist = TABLETPC_DPI * GetParameterDouble("StackedDotMaxHorzDist");
	MaxSplitDistance = TABLETPC_DPI * GetParameterDouble("MaxSplitDistance");
}
static int _e = RegisterParameterCallback(&calc_grouping_metrics);

static inline double
inches_to_dpi(double inches)
{ return TABLETPC_DPI * inches; }


static double
distance_threshold(const Rect<long> &bbox1, const Rect<long> &bbox2) {
	long maxsz = std::max(std::max(bbox1.width(), bbox1.height()), std::max(bbox2.width(), bbox2.height()));
	double inches = (double)maxsz/TABLETPC_DPI;
	double f;
	if (inches > 0.25) {
		f = DistanceThresholdRatio;///(6*inches);
	}
	else {
		f = 1 - inches/0.25 * (1-DistanceThresholdRatio);
	}
	VERBOSE(*verb_out << "distance_threshold: maxsz " << maxsz << " is " << inches << " inches, giving threshold factor " << f << std::endl);
	return std::min((double)MaxSplitDistance, std::max((double)MinSplitDistance, f * maxsz));
}


static double
mingrpdist(const group *g1, const group *g2) {
	double mind = std::numeric_limits<double>::infinity();
	for (size_t i = 0; i < g1->strokes.nstrokes; ++i) {
		for (size_t j = 0; j < g2->strokes.nstrokes; ++j) {
			double d = dist_sq(g1->strokes.strokes[i], g2->strokes.strokes[j]);
			mind = std::min(mind, d);
		}
	}
	return std::sqrt(mind);
}

void
group::mulscore(double x) {
	proximity_score *= x;
	stack_score *= x;
	/*
	if (parents[0]) {
		parents[0]->mulscore(x);
	}
	if (parents[1]) {
		parents[1]->mulscore(x);
	}*/
}

/*
void
dilute_group_by_dist(group *grp, const std::map<bitvec, group *> &groups) {
	VERBOSE(*verb_out << "dilute_group_by_dist(" << grp->bits << ")\n");
	for (std::map<bitvec, group *>::const_iterator i = groups.begin(); i != groups.end(); ++i) {
		const group *cmpgrp = i->second;
		if (null_intersection(grp->bits, cmpgrp->bits)) {
			double thres = distance_threshold(grp->bounds, cmpgrp->bounds);
			double dist = mingrpdist(grp, cmpgrp);
			if (dist < thres/2) {
				double dfac = dist/(thres/2);
				double ofac = 1-overlap_proportion(grp->bounds, cmpgrp->bounds);
				double factor = 0.5*(dfac+ofac);
				grp->dil_factor *= factor;
				//grp->proximity_score *= factor;
				//grp->stack_score *= factor;
				VERBOSE(*verb_out << " on cmp with group " << cmpgrp->bits << ", diluting by distance factor " << dfac << " and overlap factor" << ofac << " (combined factor " << factor << ") to " << grp->dil_factor*grp->proximity_score << "," << grp->dil_factor*grp->stack_score << std::endl);
			}
		}
	}
}*/

static double
map_score_linear(double x, double parm)
{ return std::max(1.0 - std::abs(x) / parm, 0.0); }

static double
map_score_exp(double x, double mean)
{ return std::exp(-1.0/mean * std::abs(x)); }

static inline double
map_score(double x, double parm)
{ return map_score_linear(x, parm); }


struct bbox_alignment_meas : public segment_distance_meas
{
	Rect<long> bounds;

	bbox_alignment_meas(const segment *seg) : segment_distance_meas(seg->stk), bounds(seg->bounds) { }

	double operator()(const segment *seg) const
	{
		double side_scores[4];

		Rect<long> segbounds = seg->bounds;

		double thres = std::min((double)SplitDistance, std::max(distance_threshold(bounds, segbounds), MinSplitDistance));

		side_scores[0] = map_score(bounds.right - segbounds.right, thres);
		side_scores[1] = map_score(bounds.top - segbounds.top, thres);
		side_scores[2] = map_score(bounds.left - segbounds.left, thres);
		side_scores[3] = map_score(bounds.bottom - segbounds.bottom, thres);

		VERBOSE(*verb_out << " alignment measurements: " << side_scores[0] << ", " << side_scores[1] << ", " << side_scores[2] << ", " << side_scores[3] << std::endl);

		double best_x = std::max(side_scores[0], side_scores[2]);
		double best_y = std::max(side_scores[1], side_scores[3]);
		if (best_x > 0 && best_y > 0) {
			VERBOSE(*verb_out << " score: " << std::sqrt(best_x * best_y) << std::endl);
			return std::sqrt(best_x * best_y);
		}
		else if (best_x > 0) {
			double center_score = map_score(bounds.y_center() - segbounds.y_center(), thres);
			VERBOSE(*verb_out << " score: " << std::sqrt(best_x * center_score) << std::endl);
			return std::sqrt(center_score * best_x);//(center_score == 0) ? 0 : 0.5 * (center_score + best_x);
		}
		else if (best_y > 0) {
			double center_score = map_score(bounds.x_center() - segbounds.x_center(), thres);
			VERBOSE(*verb_out << " score: " << std::sqrt(center_score * best_y) << std::endl);
			return std::sqrt(center_score * best_y);//(center_score == 0) ? 0 : 0.5 * (center_score + best_y);
		}
		else {
			return 0;
		}
	}
};


static double
add_proximity_score(double sum, const group &group)
{
	return sum + group.proximity_score;
}


RawStroke
strip_duplicate_points(const RawStroke &s) {
	long *x, *y;
	size_t n = s.npoints;
	long *endx = s.x + n;
	long x0 = *s.x, y0 = *s.y;
	if (n > 0) {
		for (x = s.x + 1, y = s.y + 1; x != endx; ++x, ++y) {
			if (*x == x0 && *y == y0) {
				--n;
			}
			else {
				x0 = *x;
				y0 = *y;
			}
		}
	}

	long *xdst = new long[n];
	long *ydst = new long[n];
	unsigned long *tdst = s.time ? new unsigned long[n] : 0;
	if (n > 0) {
		long *xd = xdst, *yd = ydst;
		unsigned long *td = tdst;

		x0 = *s.x, y0 = *s.y;
		*(xd++) = x0;
		*(yd++) = y0;
		if (tdst) {
			*(td++) = *s.time;
		}

		for (x = s.x + 1, y = s.y + 1; x != endx; ++x, ++y) {
			if (*x != x0 || *y != y0) {
				x0 = *x;
				*(xd++) = x0;
				y0 = *y;
				*(yd++) = y0;
				if (tdst) {
					*(td++) = s.time[x - s.x];
				}
			}
		}
	}

	return RawStroke(xdst, ydst, tdst, n);
}


group *
build_proximity_group(group *base, group *single_seg_grp, const std::map<bitvec, group *> allgroups) {
	if ((base->tags & group_tags::KNOWN_SYMBOL)) {// || seg->stk.npoints == 0) {
		return 0;
	}

	if (base->bits.count_set_bits() >= MaxGroupSize) {
		return 0;
	}

	VERBOSE(*verb_out << "building proximity group between base group " << base->bits << " at " << base->bounds << " and segment at " << single_seg_grp->bounds << std::endl);

	/*
	double containment_score = overlap_proportion(base->bounds, seg->bounds);
	VERBOSE(*verb_out << " containment score = " << containment_score << std::endl);
	*/

	/*
	std::pair<std::list<segment *>::const_iterator, double> minpair = maximizer(base->children.begin(), base->children.end(), 0.0, bbox_alignment_meas(seg));

	double alignment_score = minpair.second;
	VERBOSE(*verb_out << " alignment score = " << alignment_score << std::endl);*/

	const segment *seg = single_seg_grp->children.front();
	/*
	std::pair<std::list<segment *>::const_iterator, double> minpair;
	minpair = minimizer(base->children.begin(), base->children.end(), std::numeric_limits<double>::infinity(), stroke_distance_meas(seg->stk));
	const segment *closest_seg = *minpair.first;
	double mindist = std::sqrt(minpair.second);*/
	double mindist = std::numeric_limits<double>::infinity();
	for (std::list<segment *>::const_iterator i = base->children.begin(); i != base->children.end(); ++i) {
		double d = dist_sq(seg->stk, (*i)->stk);
		mindist = std::min(mindist, d);
	}
	mindist = std::sqrt(mindist);

	double boxdist = std::sqrt((double)dist_sq(base->bounds, single_seg_grp->bounds));
	VERBOSE(*verb_out << "boxdist is " << boxdist << " compared to mindist " << mindist << " and upper bound " << 1.25*boxdist << std::endl);
	if (boxdist < mindist && boxdist * 1.5 > mindist) {
		mindist = boxdist;
	}

	double thres = distance_threshold(single_seg_grp->bounds, base->bounds);//std::min(SplitDistance, std::max(distance_threshold(seg->bounds, closest_seg->bounds), MinSplitDistance));
	//double thres = distance_threshold(seg->bounds, closest_seg->bounds);//std::min(SplitDistance, std::max(distance_threshold(seg->bounds, closest_seg->bounds), MinSplitDistance));
	//std::cout << thres << '\t' << MinSplitDistance << '\t' << MaxSplitDistance << std::endl;
	//thres *= thres;
/*	double autogroup_thres = thres * ForceDistanceRatio * ForceDistanceRatio;

	double distance_score;
	if (mindist < autogroup_thres) {
		distance_score = 1.0;
	}
	else {
		distance_score = std::max(1.0 - mindist / thres, 0.0);
	}
	*/
	double overlap = overlap_proportion(base->bounds, single_seg_grp->bounds);
	//thres *= 1+overlap;

	double base_score = std::max(base->proximity_score, base->stack_score);
	double distance_score = map_score(mindist, thres);

	VERBOSE(*verb_out << " dist " << mindist << " cmp. thresh. " << thres << " => distance score = " << distance_score << std::endl);

	double proximity_score = distance_score;//std::max(alignment_score, distance_score);//(alignment_score + distance_score) / 2;//std::max(alignment_score, distance_score);
	//double proximity_score = std::max(std::max(containment_score, alignment_score), distance_score);

	double container_score = std::max(base->container_score, single_seg_grp->container_score);
	const static double CONTAINER_MINSC = -1;
	const static double CONTAINER_MAXSC = 1;
	double contfac;
	if (overlap < 2*container_score) {
		contfac = CONTAINER_MINSC * (2*container_score-overlap)/(2-overlap);
	}
	else {
		contfac = CONTAINER_MAXSC * (1 - overlap/2);
	}
	VERBOSE(*verb_out << "adjusting by container/overlap " << container_score << "/" << overlap << " for factor " << contfac << std::endl);
	contfac = (contfac - CONTAINER_MINSC) / (CONTAINER_MAXSC - CONTAINER_MINSC);
	//double container_adj = 0.25 * (overlap - 0.25 - container_score);
	/*
	if (overlap < 1 - container_score) {
		proximity_score += 0.08 * (overlap - container_score);
	}
	else {
		proximity_score += 0.08 * (2*overlap-1);
	}*/
	if (proximity_score > 0) {
		contfac = 0.5*contfac + 0.75;
		VERBOSE(*verb_out << " (Factor maps to " << contfac << ")\n");
		proximity_score *= contfac;
	}
	//else {
	//	proximity_score = (1-container_score) * overlap;
	//}
	//proximity_score += container_adj;
	//proximity_score = std::min(1.0, std::max(0.0, proximity_score));
	//proximity_score *= 1-0.5*container_score*overlap;
	if (proximity_score <= GROUP_SCORE_THRESHOLD) {
		VERBOSE(*verb_out << " no group exists here\n");
		return 0;
	}
	
	// TODO: test new seg. with the next line changed back pure geom. product, not prod. over # inter-stk links
	//double final_score = std::pow(std::pow(base_score, static_cast<int>(base->children.size()-1)) * proximity_score, 1.0 / (base->children.size()));//(base->children.size() * base->proximity_score + proximity_score) / (base->children.size() + 1);
	double final_score = proximity_score * base_score;
	VERBOSE(*verb_out << " found a group with score " << proximity_score << " giving final score " << final_score << std::endl);
	if (proximity_score <= GROUP_SCORE_THRESHOLD || final_score <= GROUP_SCORE_THRESHOLD) {
		VERBOSE(*verb_out << " final score " << final_score << " is too low\n");
		return 0;
	}

	group *newgroup = new group;
	//newgroup->parents[0] = base;
	//newgroup->parents[1] = single_seg_grp;
	newgroup->bounds = merge(base->bounds, single_seg_grp->bounds);
	newgroup->bits = bitvec_union(base->bits, single_seg_grp->bits);
	//newgroup->bits.set(segi, true);
	VERBOSE(*verb_out << " constructing new proximity group " << newgroup << " at " << newgroup->bits << " with score " << final_score << std::endl);
	newgroup->children = base->children;
	newgroup->children.insert(newgroup->children.end(), single_seg_grp->children.begin(), single_seg_grp->children.end());
	RawStrokeGroup G = build_stroke_group(newgroup->children.begin(), newgroup->children.end(), newgroup->bounds, newgroup->children.size());
	newgroup->strokes = G;
	newgroup->nstrokes = normalize(newgroup->strokes);
	//NormalizedStrokeGroup ns = normalize(newgroup->strokes);
	//NormalizedStrokeGroup ss = subdivide(ns);
	//newgroup->nstrokes = subdivide(ns);

	//double totoverlap = 0;
	//final_score *= 1.0 - 0.75*std::min(1.0, totoverlap);
	
	newgroup->proximity_score = final_score;
	newgroup->container_score = container_score;
	//newgroup->stack_score = std::max(base->stack_score, final_score);
	//dilute_group_by_dist(newgroup, allgroups);

	for (std::list<segment *>::iterator i = newgroup->children.begin(); i != newgroup->children.end(); ++i) {
		(*i)->parents.push_front(newgroup);
	}

	/*
	if (base->children.size() == 1) {
		base->proximity_score *= 1.0 - 0.75*newgroup->proximity_score;
	}
	single_seg_grp->proximity_score *= 1.0 - 0.75*newgroup->proximity_score;
	*/
	//base->proximity_score = std::max(base->proximity_score - newgroup->proximity_score, 0.0);
	//base->stack_score = std::max(base->stack_score - newgroup->proximity_score, 0.0);
	//base->mulscore(1.0 - 0.75*newgroup->proximity_score);
	//base->proximity_score *= 1.0 - 0.75*newgroup->proximity_score;
	//base->stack_score *= 1.0 - 0.75*newgroup->proximity_score;
	
	//VERBOSE(*verb_out << "base group at " << base->bits << " had scores updated to " << base->proximity_score << "; " << base->stack_score << std::endl);

	return newgroup;
}


group *
math_recognizer_base::build_stack_group(group *g1, group *g2)
{
	if (g1->tags & group_tags::KNOWN_SYMBOL || g2->tags & group_tags::KNOWN_SYMBOL) {
		return 0;
	}
	
	if (g1->bits.count_set_bits() + g2->bits.count_set_bits() > MaxGroupSize) {
		return 0;
	}

	if (!null_intersection(g1->bits, g2->bits)) {
		return 0;
	}


	VERBOSE(*verb_out << "building stack group between group " << g1->bits << " at " << g1->bounds << " and group " << g2->bits << " at " << g2->bounds << std::endl);

	Rect<long> bounds = merge(g1->bounds, g2->bounds);

	for (size_t i = 0; i < segments.size(); ++i) {
		if (!(g1->bits.at(i) || g2->bits.at(i))) {
			if (overlap_proportion(segments[i]->bounds, bounds) > 0) {
				VERBOSE(*verb_out << "failed due to overlap with segment " << i << " at " << segments[i]->bounds << std::endl);
				return 0;
			}
		}
	}

	/*if (!(g1->bounds.bottom >= g2->bounds.top || g2->bounds.bottom >= g1->bounds.top)) {
		VERBOSE(*verb_out << "failed topness test\n");
		return 0;
	}*/
	
	double horz_score;
	bool dotmode = false;
	if (bbox_is_dot(g1->bounds) || bbox_is_dot(g2->bounds)) {
		if (bbox_is_dot(g1->bounds) && !bbox_is_dot(g2->bounds)) {
			if (g2->bounds.width() > std::max(10*g1->bounds.width(), TABLETPC_DPI)) {
				VERBOSE(*verb_out << "failed box shape test\n");
				return 0;
			}
		}
		else if (bbox_is_dot(g2->bounds) && !bbox_is_dot(g1->bounds)) {
			if (g1->bounds.width() > std::max(10*g2->bounds.width(), TABLETPC_DPI)) {
				VERBOSE(*verb_out << "failed box shape test\n");
				return 0;
			}
		}
		VERBOSE(*verb_out << "using dot mode\n");
		double cdist = std::abs(0.5*((g1->bounds.left + g1->bounds.right) - (g2->bounds.left + g2->bounds.right)));
		double thres = std::max(StackedDotMaxHorzDist, (double)std::max(g1->bounds.width(), g2->bounds.width()));
		VERBOSE(*verb_out << "x-center dist is " << cdist << " vs threshold " << thres << std::endl);
		//thres = StackedDotMaxHorzDist;
		horz_score = map_score(cdist, thres);
		dotmode = true;
	}
	else {
		VERBOSE(*verb_out << "using normal mode\n");
		if (g1->bounds.width() * 2 < g1->bounds.height()
		 || g2->bounds.width() * 2 < g2->bounds.height()) {
			VERBOSE(*verb_out << "failed box shape test\n");
			return 0;
		}
		if (g1->bounds.width() * 2 < g2->bounds.width()
		 || g2->bounds.width() * 2 < g1->bounds.width()) {
			VERBOSE(*verb_out << "failed box shape test 2\n");
			return 0;
		}

		long intersect_left = std::max(g1->bounds.left, g2->bounds.left);
		long intersect_right = std::min(g1->bounds.right, g2->bounds.right);
		long horz_intersect = intersect_right - intersect_left;
		VERBOSE(*verb_out << "horz_intersect=" << horz_intersect << std::endl);
		
		long minwidth = std::min(g1->bounds.width(), g2->bounds.width());
		long maxwidth = std::max(g1->bounds.width(), g2->bounds.width());
		double req_intersect_len = StackOverlapRequiredRatio * maxwidth;
		double limit = minwidth - req_intersect_len;
		if (limit > 0) {
			VERBOSE(*verb_out << "req_intersect_len=" << req_intersect_len << std::endl);
			VERBOSE(*verb_out << "limit=" << limit << std::endl);
			horz_score = map_score(limit - (horz_intersect - req_intersect_len), limit);
		}
		else {
			horz_score = 0.0;
		}
	}

	if (horz_score <= GROUP_SCORE_THRESHOLD) {		
		VERBOSE(*verb_out << "failed horizontal test\n");
		return 0;
	}

	VERBOSE(*verb_out << " horz_score = " << horz_score << std::endl);

	long vert_dist = std::min(std::abs(g1->bounds.bottom - g2->bounds.top), std::abs(g2->bounds.bottom - g1->bounds.top));
	double thres = distance_threshold(g1->bounds, g2->bounds);
	if (dotmode) {
		thres = 2*std::max(g1->bounds.height(), g2->bounds.height());
	}
	else {
		thres = std::min(0.25*TABLETPC_DPI, (double)std::max(g1->bounds.width(), g2->bounds.width()));
	}
	//double test_thres = 2 * thres * map_score_exp(thres, inches_to_dpi(0.25));
	//long vert_dist_thres = std::min(static_cast<long>(test_thres), 2 * SplitDistance);
	double vert_score = map_score(vert_dist, 2 * thres);//vert_dist_thres);
	
	VERBOSE(*verb_out << "vert_dist = " << vert_dist << std::endl);
	//VERBOSE(*verb_out << "vert_thesh = " << thres << " -> " << test_thres << " -> " << vert_dist_thres << std::endl);

	if (vert_score <= GROUP_SCORE_THRESHOLD) {
		VERBOSE(*verb_out << " failed vertical test\n");
		return 0;
	}

	VERBOSE(*verb_out << " vert_score = " << vert_score << std::endl);

	double g1score = std::max(g1->proximity_score, g1->stack_score);
	double g2score = std::max(g2->proximity_score, g2->stack_score);
	double stack_score = std::sqrt(g1score * g2score) * (0.35 * vert_score + 0.65 * horz_score);

	double overlap = overlap_proportion(g1->bounds, g2->bounds);
	VERBOSE(*verb_out << " tempering stack_score " << stack_score << " by overlap proportion " << (0.5*(1-overlap)) << std::endl);
	stack_score *= (1-0.5*overlap);

	if (stack_score <= GROUP_SCORE_THRESHOLD) {
		VERBOSE(*verb_out << " final score " << stack_score << " is too low\n");
		return 0;
	}

	VERBOSE(*verb_out << " found a group with score " << stack_score << std::endl);
	//stack_score *= std::min((1-g1->container_score), (1-g2->container_score));
	//VERBOSE(*verb_out << " reduced to " << stack_score << " by containers" << std::endl);
	double container_score = std::max(g1->container_score, g2->container_score);
	VERBOSE(*verb_out << "reducing by container score " << container_score << "*" << overlap);
	stack_score *= 1-0.5*container_score*overlap;
	VERBOSE(*verb_out << " for final score " << stack_score << std::endl);

	bitvec bits = bitvec_union(g1->bits, g2->bits);
	group *&newgroup = (*groups)[bits];
	if (!newgroup) {
		newgroup = new group;
		//newgroup->parents[0] = g1;
		//newgroup->parents[1] = g2;
		newgroup->bits = bitvec_union(g1->bits, g2->bits);
		VERBOSE(*verb_out << " constructing new stack group " << newgroup << " at " << newgroup->bits << std::endl);
		newgroup->bounds = merge(g1->bounds, g2->bounds);
		newgroup->children = g1->children;
		newgroup->children.insert(newgroup->children.end(), g2->children.begin(), g2->children.end());
		RawStrokeGroup G = build_stroke_group(newgroup->children.begin(), newgroup->children.end(), newgroup->bounds, newgroup->children.size());
		newgroup->strokes = G;
		//NormalizedStrokeGroup ns = normalize(newgroup->strokes);
		//NormalizedStrokeGroup ss = subdivide(ns);
		//newgroup->nstrokes = subdivide(ns);
		newgroup->nstrokes = normalize(newgroup->strokes);
		
		for (std::list<segment *>::iterator i = newgroup->children.begin(); i != newgroup->children.end(); ++i) {
			(*i)->parents.push_front(newgroup);
		}

		double g1score = std::max(g1->proximity_score, g1->stack_score);
		double g2score = std::max(g2->proximity_score, g2->stack_score);
		stack_score *= g1score * g2score;
		//newgroup->proximity_score = std::sqrt(g1->proximity_score*g2->proximity_score);
		newgroup->stack_score = std::max(newgroup->stack_score, stack_score);
		newgroup->container_score = std::max(g1->container_score, g2->container_score);
		//dilute_group_by_dist(newgroup, allgroups);
		
		//g1->mulscore(1.0 - newgroup->stack_score);
		//g2->mulscore(1.0 - newgroup->stack_score);
		//g1->proximity_score *= 1.0 - newgroup->stack_score;
		//g1->stack_score *= 1.0 - newgroup->stack_score;
		//g2->proximity_score *= 1.0 - newgroup->stack_score;
		//g2->stack_score *= 1.0 - newgroup->stack_score;

		//g1->proximity_score = std::max(g1->proximity_score - newgroup->stack_score, 0.0);
		//g1->stack_score = std::max(g1->stack_score - newgroup->stack_score, 0.0);
		//g2->proximity_score = std::max(g2->proximity_score - newgroup->stack_score, 0.0);
		//g2->stack_score = std::max(g2->stack_score - newgroup->stack_score, 0.0);

		//VERBOSE(*verb_out << "group 1 at " << g1->bits << " had scores updated to " << g1->proximity_score << "; " << g1->stack_score << std::endl);
		//VERBOSE(*verb_out << "group 2 at " << g2->bits << " had scores updated to " << g2->proximity_score << "; " << g2->stack_score << std::endl);

		return newgroup;
	}

	//dilute_group_by_dist(newgroup, allgroups);
	if (stack_score > newgroup->stack_score) {
		double diff = stack_score - newgroup->stack_score;
		newgroup->stack_score = stack_score;

		//g1->mulscore(1.0 - newgroup->stack_score);
		//g2->mulscore(1.0 - newgroup->stack_score);
		//g1->proximity_score *= 1.0 - newgroup->stack_score;
		//g1->stack_score *= 1.0 - newgroup->stack_score;
		//g2->proximity_score *= 1.0 - newgroup->stack_score;
		//g2->stack_score *= 1.0 - newgroup->stack_score;

		//g1->proximity_score = std::max(g1->proximity_score - diff, 0.0);
		//g1->stack_score = std::max(g1->stack_score - diff, 0.0);
		//g2->proximity_score = std::max(g2->proximity_score - diff, 0.0);
		//g2->stack_score = std::max(g2->stack_score - diff, 0.0);		

		VERBOSE(*verb_out << "group 1 at " << g1->bits << " had scores updated to " << g1->proximity_score << "; " << g1->stack_score << std::endl);
		VERBOSE(*verb_out << "group 2 at " << g2->bits << " had scores updated to " << g2->proximity_score << "; " << g2->stack_score << std::endl);
	}

	return 0;
}


static bool
group_is_dead(const group *G)
{
	return std::max(G->proximity_score, G->stack_score) <= GROUP_SCORE_THRESHOLD;
}

int
math_recognizer_base::update_proximity_groups(std::map<bitvec, group *> &testgrps, std::map<bitvec, group *> &pgroups, group *single_seg_grp, size_t segi)
{
	if (segments.empty()) {
		return 0;
	}

	VERBOSE(*verb_out << "updating groups with single seg group " << single_seg_grp->bits << std::endl);

	//segment *seg = segments[segi];

	//std::list<group *> rmgroups;
	for (std::map<bitvec, group *>::iterator i = testgrps.begin(); i != testgrps.end(); ++i) {
		group *basegroup = i->second;
		
		if (!basegroup->tags & group_tags::KNOWN_SYMBOL && !basegroup->bits.at(segi) && (basegroup->proximity_score > 0.0 || basegroup->stack_score > 0.0)) {
			bitvec newbits = basegroup->bits;
			newbits.set(segi);
			if (groups->find(newbits) == groups->end() && pgroups.find(newbits) == pgroups.end()) {
				group *newgroup = build_proximity_group(basegroup, single_seg_grp, *groups);
				if (newgroup) {
					pgroups[newgroup->bits] = newgroup;
	
					/*
					if (group_is_dead(basegroup)) {
						VERBOSE(*verb_out << "basegroup at " << basegroup->bits << " found dead after addition of proximity group " << newgroup->bits << "; dumping...\n");
						rmgroups.push_back(basegroup);
					}*/
				}
			}
		}
	}

	/*
	for (std::list<group *>::iterator i = rmgroups.begin(); i != rmgroups.end(); ++i) {
		remove_group(*i);
	}*/

	return 0;
}


int
math_recognizer_base::update_stack_groups(std::list<group *> &cmpgroups, std::list<group *> &sgroups)
{
	if (segments.empty()) {
		return 0;
	}

	std::list<group *> erase_cmpgroups;
	
	for (std::map<bitvec, group *>::iterator i = groups->begin(); i != groups->end(); ++i) {
		group *basegroup = i->second;
		if (basegroup->tags & group_tags::KNOWN_SYMBOL) continue;
		
		for (std::list<group *>::iterator j = cmpgroups.begin(); j != cmpgroups.end(); ++j) {
			group *cmpgroup = *j;
			
			if (cmpgroup != basegroup) {
				group *newgroup = build_stack_group(basegroup, cmpgroup);
				if (newgroup) {
					sgroups.push_front(newgroup);

					if (!(cmpgroup->tags & group_tags::DEAD_GROUP) && group_is_dead(cmpgroup)) {
						VERBOSE(*verb_out << "cmpgroup at " << cmpgroup->bits << " found dead after addition of stack group " << newgroup->bits << "; dumping...\n");
						erase_cmpgroups.push_back(cmpgroup);
						cmpgroup->tags = cmpgroup->tags | group_tags::DEAD_GROUP;
					}
				}
			}
		}

		/*
		if (group_is_dead(basegroup)) {
			VERBOSE(*verb_out << "basegroup at " << basegroup->bits << " found dead after stack group insertions; dumping\n");
			std::list<group *>::iterator j = std::find(cmpgroups.begin(), cmpgroups.end(), basegroup);
			if (j != cmpgroups.end()) {
				cmpgroups.erase(j);
				if (basegroup->tags & group_tags::DEAD_GROUP) {
					erase_cmpgroups.remove(basegroup);
				}
			}
			i = remove_group(basegroup);
		}
		else {
			++i;
		}*/
	}
	
	while (!erase_cmpgroups.empty()) {
		remove_group(erase_cmpgroups.front());
		erase_cmpgroups.pop_front();
	}

	return 0;
}


}
