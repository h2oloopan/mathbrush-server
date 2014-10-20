/*
 * interp.h
 *
 * This file defines routines for interpolation of strokes and stroke groups.
 *
 */
 
#ifndef INTERP_H_
#define INTERP_H_


#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>

#include "dist.h"
#include "error.h"
#include "parms.h"
#include "stroke.h"
#include "stkutils.h"
#include "ink-io.h"


namespace scg
{

extern double CuspLowerBound;
extern double CuspUpperBound;
extern double InterpArclengthRatio;

template <typename T>
T
round(double x)
{
	T ret = static_cast<T>(x);
	if (x - static_cast<double>(ret) >= 0.5) {
		ret++;
	}
	return ret;
}


template <typename T>
void
find_cusps(const Stroke<T> &stroke, std::vector<size_t> &cusps) {
    cusps.push_back(0);
	for (size_t i = 1; i < stroke.npoints-1; ++i) {
		T df = dist_sq(stroke.x[i], stroke.y[i], stroke.x[i+1], stroke.y[i+1]);
		T db = dist_sq(stroke.x[i], stroke.y[i], stroke.x[i-1], stroke.y[i-1]);
		T ds = dist_sq(stroke.x[i-1], stroke.y[i-1], stroke.x[i+1], stroke.y[i+1]);
		double t = std::acos((df + db - ds) / (2*std::sqrt(df*db)));
		if (t <= CuspLowerBound || t >= CuspUpperBound) {
			cusps.push_back(i);
		}
	}
    cusps.push_back(stroke.npoints - 1);    
}


/*
 * optimal_stroke_size computes a number of points that is small, but
 * sufficient to preserve the shape of a stroke.
 */
template <typename T>
size_t
optimal_stroke_size(const Stroke<T> &stroke, const Rect<T> &group_bounds)
{
	size_t sz = stroke.npoints;
	if (sz < 5) {
		return sz;
	}
    
    //return std::max<unsigned>(5, sz / 3 + 1);
    
	std::vector<size_t> cusps;
	find_cusps(stroke, cusps);

    Rect<T> bounds = stroke.bounds();
    //double arclen_threshold = std::max(InterpArclength, (InterpArclength / NormalizationScale) * std::max(width(bounds), height(bounds)));
	 double arclen_threshold = InterpArclengthRatio;
	
	double arclen = stroke.length();
	
	if (stroke_is_dot(stroke, group_bounds)) {
    	return std::min<size_t>(6, static_cast<unsigned>(arclen / arclen_threshold) + cusps.size() + 2);
	}
	else {
    	return std::max<size_t>(10, static_cast<unsigned>(arclen / arclen_threshold) + cusps.size() + 1);
	}	
}


/* subdivide a stroke so that it contains a given number of points.  This
 * function preserves end- and cusp-points, and performs linear interpolation
 * between adjacent points, spacing out 'npoints' points at equal intervals along
 * the stroke curve.
 */
template <typename T>
Stroke<T>
linear_subdivide(const Stroke<T> &stroke, size_t npoints)
{
    size_t req_npoints = npoints;
    
	T *x = new T[npoints];
	T *y = new T[npoints];
	
	std::vector<size_t> cusps;
	find_cusps(stroke, cusps);
	
	size_t num_cusps = cusps.size();
	
	size_t curr_pt = 0;

	if (cusps.size() > npoints) {
		double cusps_per_pt = (double)cusps.size() / npoints;
		for (size_t i = 0; i < npoints; i++) {
			size_t cusp = (size_t)(cusps_per_pt * i);
			x[curr_pt] = stroke.x[cusps[cusp]];
			y[curr_pt] = stroke.y[cusps[cusp]];
			curr_pt++;
		}
		return Stroke<T>(x, y, 0, npoints);
	}

	npoints -= cusps.size();
	
	size_t sz = stroke.npoints;
	
	double total_len = 0.0;
	std::vector<double> arclen;
	arclen.push_back(0.0);
	for (size_t i = 0; i < sz - 1; i++) {
		total_len += std::sqrt(static_cast<double>(dist_sq(stroke.x[i], stroke.y[i], stroke.x[i + 1], stroke.y[i + 1])));
		arclen.push_back(total_len);
	}
	
	for (std::vector<size_t>::iterator curr_cusp = cusps.begin(); curr_cusp != cusps.end() - 1; curr_cusp++) {
		// add the current cusp to the new curve
		x[curr_pt] = stroke.x[*curr_cusp];
		y[curr_pt] = stroke.y[*curr_cusp];
		curr_pt++;

		// we want to fill the new curve up to the next cusp in this loop iteration
		double curve_len = arclen[*(curr_cusp + 1)];
		double portion_of_curve = curve_len / total_len;
		
		size_t num_pts = round<size_t>(npoints * portion_of_curve);
		// by the end of this loop iteration, there should be "num_pts" points in the new curve
		
		// find the arclength distance needed to evenly distribute "num_pts" points between the current
		// cusp and the next one
		double len_between = curve_len - arclen[*curr_cusp];
		double pt_dist = len_between / (num_pts + 1);
		double desired_len = arclen[*curr_cusp] + pt_dist;
		
		npoints -= num_pts;
		size_t startpt = *curr_cusp;
		size_t currpt = *curr_cusp;
		for (size_t j = 0; j < num_pts; ++j) {
			while (arclen[currpt+1] < desired_len) {
				currpt++;
			}
			double t = (desired_len - arclen[startpt]) / (arclen[currpt+1] - arclen[startpt]);
			double px = t*stroke.x[startpt] + (1-t)*stroke.x[currpt+1];
			double py = t*stroke.y[startpt] + (1-t)*stroke.y[currpt+1];
			x[curr_pt] = static_cast<T>(px);
			y[curr_pt] = static_cast<T>(py);
			startpt = currpt;
			curr_pt++;
			desired_len += pt_dist;
		}

		/*
		while (num_pts--) {
			desired_len += pt_dist;
			
			size_t last_point = *curr_cusp + 2;
			for (size_t i = last_point; i <= *(curr_cusp + 1); i++) {
				double point_len = arclen[i];
				if (point_len >= desired_len) {
					unsigned npointshere = 1 + (unsigned)((point_len-desired_len)/pt_dist);
					// we need to add a point here...interpolate between the two closest sample points
					double real_dist = arclen[i] - arclen[i - 1];
					double first_proportion = ((point_len - desired_len) / real_dist);
					double second_proportion = 1.0 - ((point_len - desired_len) / real_dist);
					
					double px = stroke.x[i - 1] * first_proportion + stroke.x[i] * second_proportion;
					double py = stroke.y[i - 1] * first_proportion + stroke.y[i] * second_proportion;
					assert(px >= 0 && py >= 0);
					
					x[curr_pt] = static_cast<T>(px);
					y[curr_pt] = static_cast<T>(py);
					curr_pt++;
					
					last_point = i;
					break;
				}
			}
		}*/
	}
	
	x[curr_pt] = stroke.x[cusps.back()];
	y[curr_pt] = stroke.y[cusps.back()];
	++curr_pt;
	assert(curr_pt == req_npoints);
	/*
	while (curr_pt != req_npoints) {
	    x[curr_pt] = stroke.x[cusps.back()];
	    y[curr_pt] = stroke.y[cusps.back()];
	    curr_pt++;
    }*/
    
	return Stroke<T>(x, y, 0, req_npoints);
}


struct CurvatureAtPoint
{
    CurvatureAtPoint(size_t index_, double curvature_) : index(index_), curvature(curvature_) {}
    
    size_t index;
    double curvature;
    
    bool operator<(const CurvatureAtPoint &rhs) const
      { return curvature < rhs.curvature; }
};


template <typename T>
Stroke<T>
subdivide(const Stroke<T> &stroke, size_t npoints)
{
    if (stroke.npoints < npoints || stroke.npoints < 3) {
        return linear_subdivide(stroke, npoints);
    }
    
    if (stroke.npoints <= npoints) {
        return stroke.copy();
    }
    
    std::vector<size_t> selected_points;
    selected_points.reserve(npoints);
    selected_points.push_back(0);
    
    std::vector<CurvatureAtPoint> curvature;
    //std::vector<double> curvature;    
    curvature.reserve(stroke.npoints - 1);
    
    unsigned ix = 1;
    T *x, *y;
    T *endx = stroke.x + stroke.npoints - 1;
    for (x = stroke.x + 1, y = stroke.y + 1; x < endx; x++, y++) {
        double xd = static_cast<double>((*(x + 1) - *x) - (*x - *(x - 1)));
        double yd = static_cast<double>((*(y + 1) - *y) - (*y - *(y - 1)));
        double curve = (xd / std::max(std::abs(*(x + 1) - *(x - 1)), T(1)))
                     + (yd / std::max(std::abs(*(y + 1) - *(y - 1)), T(1)));
        curvature.push_back(CurvatureAtPoint(ix++, curve));
        //curvature.push_back(curve);
    }
    
    std::sort(curvature.begin(), curvature.end());
    
    //double ceiling = std::numeric_limits<double>::infinity();
    
    const Rect<T> bounds = stroke.bounds();
    double min_dist = std::max(bounds.width(), bounds.height()) / npoints;

    std::vector<CurvatureAtPoint>::iterator max_curvature = curvature.end() - 1;
    
    while (selected_points.size() != npoints - 1 && max_curvature != curvature.end()) {
        //double max = 0.0;
        //unsigned max_index = std::numeric_limits<unsigned>::max();
        const CurvatureAtPoint &curve = *max_curvature;
        size_t index = curve.index;
        //unsigned index = 1;
        std::vector<size_t>::const_iterator i;
        for (i = selected_points.begin(); i != selected_points.end(); i++) {
            if (scg::dist_sq(stroke.x[*i], stroke.y[*i], stroke.x[index], stroke.y[index]) <= min_dist * min_dist) {
                break;
            }
        }
        if (i == selected_points.end()) {
            selected_points.push_back(index);
            if (max_curvature == curvature.begin()) {
                max_curvature = curvature.erase(max_curvature);
                if (!curvature.empty()) {
                    max_curvature = curvature.end() - 1;
                }
                min_dist *= 0.5;
					 if (min_dist < std::numeric_limits<double>::epsilon()) {
					 	break;
					 }
            }
            else {
                max_curvature = curvature.erase(max_curvature);
                if (!curvature.empty()) {
                    --max_curvature;
                }
            }
        }
        else {
            if (max_curvature == curvature.begin()) {
                max_curvature = curvature.end() - 1;
                min_dist *= 0.5;
					 if (min_dist < std::numeric_limits<double>::epsilon()) {
					 	break;
					 }
            }
            else {
                --max_curvature;
            }
        }
    }
    
    std::sort(selected_points.begin(), selected_points.end());
    selected_points.push_back(stroke.npoints-1);
    
    T *newx = new T[selected_points.size()];
    T *newy = new T[selected_points.size()];
        
    unsigned j = 0;
    for (std::vector<size_t>::const_iterator i = selected_points.begin(); i != selected_points.end(); ++i) {
        newx[j] = stroke.x[*i];
        newy[j] = stroke.y[*i];
				assert(newx[j] >= 0 && newy[j] >= 0);
        j++;
    }

    return Stroke<T>(newx, newy, 0, selected_points.size());
}


template <typename T>
StrokeGroup<T>
subdivide(const StrokeGroup<T> &G) {
	Stroke<T> *strokes = new Stroke<T>[G.nstrokes];
	Rect<T> bounds = G.bounds();
	for (size_t i = 0; i < G.nstrokes; ++i) {
		strokes[i] = subdivide(G.strokes[i], optimal_stroke_size(G.strokes[i], bounds));
	}
	return StrokeGroup<T>(strokes, G.nstrokes);
}


template <typename T>
StrokeGroup<T>
subdivide(const StrokeGroup<T> &G, unsigned npoints) {
	Stroke<T> *strokes = new Stroke<T>[G.nstrokes];
	unsigned total_points = G.npoints();
	for (size_t i = 0; i < G.nstrokes; ++i) {
		strokes[i] = subdivide(G.strokes[i], round<unsigned>(npoints * ((double)G.strokes[i].npoints / total_points)));
	}
	return StrokeGroup<T>(strokes, G.nstrokes);
}

template <typename T>
StrokeGroup<T>
subdivide(const StrokeGroup<T> &input, const StrokeGroup<T> &model, const std::vector<size_t> &inputorder) {
	Stroke<T> *strokes = new Stroke<T>[model.nstrokes];
	// preserve original input order here so that recognition functions can unravel it using inputorder as well
	for (size_t i = 0; i < model.nstrokes; ++i) {
		strokes[inputorder[i]] = subdivide(input.strokes[inputorder[i]], model.strokes[i].npoints);
		//strokes[i] = subdivide(model.strokes[i], input.strokes[inputorder[i]].npoints);
	}
	return StrokeGroup<T>(strokes, model.nstrokes);
}

template <typename T>
void
subdivide(const StrokeGroup<T> &input, const StrokeGroup<T> &model, const std::vector<size_t> &inputorder, StrokeGroup<T> &inputout, StrokeGroup<T> &modelout) {
	inputout.nstrokes = input.nstrokes;
	modelout.nstrokes = model.nstrokes;
	inputout.strokes = new Stroke<T>[input.nstrokes];
	modelout.strokes = new Stroke<T>[model.nstrokes];
	for (size_t i = 0; i < model.nstrokes; ++i) {
		if (model.strokes[i].npoints > input.strokes[i].npoints) {
			modelout.strokes[i] = model.strokes[i].copy();
			inputout.strokes[i] = subdivide(input.strokes[inputorder[i]], model.strokes[i].npoints);
		}
		else {
			inputout.strokes[i] = input.strokes[i].copy();
			modelout.strokes[i] = subdivide(model.strokes[i], input.strokes[inputorder[i]].npoints);
		}
	}
}

}


#endif

