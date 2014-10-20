#ifndef SUBDIVIDE_H_
#define SUBDIVIDE_H_

#include "stroke.h"
#include "group.h"

namespace scg {

template <typename T>
T
round(double x) {
	T ret = T(x);
	return (x - ret >= 0.5) ? ret + 1 : ret;
}

template <typename T>
NormalizedStroke
subdivide(const Stroke<T> &stk, unsigned npoints) {
    if (stroke.npoints <= npoints || stroke.npoints < 3) {
        return linear_subdivide(stk, npoints);
    }
    
    std::vector<size_t> selected_points;
    selected_points.reserve(npoints);
    selected_points.push_back(0);
    
    std::vector<CurvatureAtPoint> curvature;
    //std::vector<double> curvature;    
    curvature.reserve(num_points(stroke) - 1);
    
    unsigned ix = 1;
    T *x, *y;
    T *endx = stroke.x + num_points(stroke) - 1;
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
    
    const Rect<T> bounds = bbox(stroke);
    double min_dist = std::max(width(bounds), height(bounds)) / npoints;

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
 
/*
        for (std::vector<double>::iterator ci = curvature.begin(); ci != curvature.end(); ++ci) {
            std::vector<unsigned>::const_iterator i;
            for (i = selected_points.begin(); i != selected_points.end(); i++) {
                if (scg::dist_sq(stroke.x[*i], stroke.y[*i], stroke.x[index], stroke.y[index]) <= min_dist * min_dist) {
                    break;
                }
            }
            if (i == selected_points.end()) {
                //double total = *xi + *yi;
                double total = *ci;
                if (std::abs(total) > std::abs(max) && std::abs(total) < ceiling) {
                    max = total;
                    max_index = index;
                }
            }
            index++;
        }
        
        if (max_index == std::numeric_limits<unsigned>::max()) {
            min_dist *= 0.5;
            ceiling = std::numeric_limits<double>::infinity();
        }
        else {
            selected_points.push_back(max_index);
            ceiling = std::abs(max);
        }*/
    }
    
    std::sort(selected_points.begin(), selected_points.end());
    selected_points.push_back(scg::num_points(stroke) - 1);
    
    T *newx = DEBUG_NEW T[selected_points.size()];
    if (!newx) {
        errval = E_OUTOFMEM;
        return TempStroke<T>();
    }
    
    T *newy = DEBUG_NEW T[selected_points.size()];
    if (!newy) {
        delete[] newx;
        errval = E_OUTOFMEM;
        return TempStroke<T>();
    }
        
    unsigned j = 0;
    for (std::vector<size_t>::const_iterator i = selected_points.begin(); i != selected_points.end(); ++i) {
        newx[j] = stroke.x[*i];
        newy[j] = stroke.y[*i];
        j++;
    }

    return TempStroke<T>(newx, newy, selected_points.size());
}

template <typename T>
NormalizedStrokeGroup
subdivide(const StrokeGroup<T> &G) {
	NormalizedStroke *strokes = new NormalizedStroke[G.nstrokes];
	unsigned npoints = G.npoints();
	for (size_t i = 0; i < G.nstrokes; ++i) {
		strokes[i] = subdivide(G.strokes[i], round<unsigned>(npoints * ((double)G.strokes[i].npoints / total_points)));
	}

	return StrokeGroup<T>(newstrokes, nstrokes);
}

}

#endif