// Implementation of Suzuki's direction element matcher.

#include "direlem.h"

#include <algorithm>
#include <cmath>

#include "debug.h"
#include "error.h"
#include "match.h"
#include "memory.h"
#include "norm.h"
#include "parms.h"
#include "stroke-alg.h"

namespace scg
{


static const unsigned MaxCells = GetParameterUnsigned("MaxCells");
static const unsigned NumCellsX = GetParameterUnsigned("NumCellsX");
static const unsigned NumCellsY = GetParameterUnsigned("NumCellsY");
static const unsigned NumDirections = GetParameterUnsigned("NumDirections");
static const unsigned StructureSize = NumCellsX * NumCellsY * NumDirections;

static int
extract_structure(const NormalizedStroke &stroke, double *structure, unsigned x_cells, unsigned y_cells, unsigned directions)
{
    if (num_points(stroke) <= 1) {
        return E_INVALID;
    }
    
    double *angles = angles_along_stroke(stroke);
    if (!angles) {
        return E_OUTOFMEM;
    }
    
    double direction_cell_span = 2.0 * M_PI / directions;
    
    for (unsigned i = 0; i < num_points(stroke) - 1; i++) {
        unsigned x_cell = static_cast<unsigned>((stroke.x[i]) / (NORM_MAX / x_cells));
        unsigned y_cell = static_cast<unsigned>((stroke.y[i]) / (NORM_MAX / y_cells));

        x_cell = std::max<unsigned>(0, std::min<unsigned>(x_cells - 1, x_cell));
        y_cell = std::max<unsigned>(0, std::min<unsigned>(y_cells - 1, y_cell));

        if (angles[i] < 0.0) {
            angles[i] += M_PI + M_PI;
        }
        unsigned d_cell = static_cast<unsigned>(angles[i] / direction_cell_span);
        
        d_cell = std::max<unsigned>(0, std::min<unsigned>(directions - 1, d_cell));
        
        double dist = std::sqrt(dist_sq(stroke.x[i], stroke.y[i], stroke.x[i + 1], stroke.y[i + 1]));
        double angle1 = angles[i] - (d_cell * direction_cell_span);
        double angle2 = ((d_cell + 1) * direction_cell_span) - angles[i];
        
        double contribution1 = (dist * angle2) / (angle1 + angle2);
        double contribution2 = (dist * angle1) / (angle1 + angle2);
        
        unsigned idx;
        idx = (y_cell * (x_cells * directions)) + (x_cell * directions) + d_cell;
        structure[idx] += contribution1;

        if (d_cell == directions - 1) {
            idx = (y_cell * (x_cells * directions)) + (x_cell * directions) + 0;
            structure[idx] += contribution2;
        }
        else {
            idx = (y_cell * (x_cells * directions)) + (x_cell * directions) + (d_cell + 1);
            structure[idx] += contribution1;
        }
    }
    
    delete[] angles;
    
    return 0;
}


static double *
extract_structure(const NormalizedStrokeGroup &group, unsigned x_cells, unsigned y_cells, unsigned directions)
{
    unsigned structure_size = x_cells * y_cells * directions;
    if (structure_size == 0) {
        errval = E_INVALID;
        return 0;
    }
    
    double *structure = DEBUG_NEW double[structure_size];
    if (!structure) {
        errval = E_OUTOFMEM;
        return 0;
    }
    
    std::fill(structure, structure + structure_size, 0.0);
    
    double *stroke_structure = DEBUG_NEW double[structure_size];
    if (!stroke_structure) {
        delete[] structure;
        errval = E_OUTOFMEM;
        return 0;
    }
    
    for (NormalizedStrokeGroup::const_iterator stroke = group.begin(); stroke != group.end(); ++stroke) {
        std::fill(stroke_structure, stroke_structure + structure_size, 0.0);
        errval = extract_structure(*stroke, stroke_structure, x_cells, y_cells, directions);
        if (FAILURE(errval)) {
            delete[] stroke_structure;
            delete[] structure;
            return 0;
        }
        
        for (unsigned i = 0; i < structure_size; i++) {
            structure[i] += stroke_structure[i];// / num_points(*stroke);
        }
    }
    
    delete[] stroke_structure;
    
    return structure;
}


int 
direction_element_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, Match &match)
{
    if (num_strokes(input) != num_strokes(model)) {
        match.num_strokes = 0;
        return 0;
    }
    
    const Rect<double> bounds = bbox(input);
    unsigned x_cells = NumCellsX;
    unsigned y_cells = NumCellsY;
    unsigned extra_directions;
    if (width(bounds) > height(bounds)) {
        y_cells = NumCellsY;
        x_cells = std::min<unsigned>(MaxCells, static_cast<unsigned>(y_cells * (width(bounds) / height(bounds))));
        //x_cells = std::max<unsigned>(1, static_cast<unsigned>(y_cells * (width(bounds) / height(bounds))));
        extra_directions = static_cast<unsigned>(width(bounds) / height(bounds)) - 1;
        //extra_directions = static_cast<unsigned>(height(bounds) / width(bounds)) - 1;
    }
    else {
        x_cells = NumCellsX;
        y_cells = std::min<unsigned>(MaxCells, static_cast<unsigned>(x_cells * (height(bounds) / width(bounds))));
        //y_cells = std::max<unsigned>(1, static_cast<unsigned>(x_cells * (height(bounds) / width(bounds))));
        extra_directions = static_cast<unsigned>(height(bounds) / width(bounds)) - 1;
        //extra_directions = static_cast<unsigned>(width(bounds) / height(bounds)) - 1;
    }
    
    unsigned ndir = NumDirections;
    
    if (extra_directions * 2 >= NumDirections - 4) {
        ndir = 4;
    }
    else {
        ndir = NumDirections - 2 * extra_directions;
    }
    
    
    double *model_structure = extract_structure(model, x_cells, y_cells, ndir);
    if (!model_structure) {
        int e = errval;
        errval = 0;
        return e;
    }
    
    double *input_structure = extract_structure(input, x_cells, y_cells, ndir);
    if (!input_structure) {
        int e = errval;
        errval = 0;
        return e;
    }
    
    double total_cost = 0.0;
    for (unsigned i = 0; i < x_cells * y_cells * ndir; i++) {
        total_cost += std::abs(model_structure[i] - input_structure[i]);
    }
    
    delete[] input_structure;
    delete[] model_structure;

    match.raw_score = total_cost / num_points(input);
    match.num_strokes = num_strokes(model);
        
    return 0;
}


}

