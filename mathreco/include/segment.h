/*
 * segment.h
 *
 * This file declares the interface to the stroke segmentation and grouping module.
 */
 
#ifndef SEGMENT_H_
#define SEGMENT_H_


#include <list>
#include <vector>

#include "group.h"
#include "stroke.h"
#include "parms.h"
#include "reco-types.h"

namespace scg
{

extern long TABLETPC_DPI;



//int update_segmentation(reco_context *ctx);
//segment *merge_segment(context *ctx, segment *seg);



struct segment_distance_meas
{
	explicit segment_distance_meas(const RawStroke &stk_) : stk(stk_) { }

	virtual double operator()(const segment *seg) const = 0;

	const RawStroke &stk;
};

struct dot_distance_meas : public segment_distance_meas
{
	static const double DY_WEIGHT;

	explicit dot_distance_meas(const segment *seg) : segment_distance_meas(seg->stk), bounds(seg->bounds) { }

	double operator()(const segment *seg) const;

	Rect<long> bounds;
};

struct stroke_distance_meas : public segment_distance_meas
{
	explicit stroke_distance_meas(const RawStroke &stroke_) : segment_distance_meas(stroke_) { }

	double operator()(const segment *seg) const;
};


}

#endif

