/*
 * merge.h
 *
 * This file implements algorithms for merging together various
 * geometrical objects.  For example, rectangles and stroke groups.
 */
 
#ifndef MERGE_H_
#define MERGE_H_


#include "rect.h"
#include "group.h"
#include "stroke.h"


namespace scg
{


NormalizedStrokeGroup merge(const NormalizedStrokeGroup &g1, const NormalizedStrokeGroup &g2, double first_weight);

NormalizedStrokeGroup merge(const NormalizedStrokeGroup &g1, const NormalizedStrokeGroup &g2);


}


#endif

