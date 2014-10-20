#include "structural.h"

#include <cmath>

#include "debug.h"
#include "elastic.h"
#include "error.h"
#include "match.h"
#include "memory.h"
#include "stroke-alg.h"


namespace scg
{


static int MaxPerturbation = 3;
static unsigned MaxShifting = 2;
static unsigned NumDirections = 16;

enum ChainCodePrimitive
{
    Line,
    ClockwiseCurve,
    AntiClockwiseCurve,
    Loop,
    Dot
};

struct ChainCodeEntry
{
    ChainCodeEntry(ChainCodePrimitive primitive_, int direction_) : primitive(primitive_), direction(direction_) {}
    
    ChainCodePrimitive primitive;
    int direction;
};


static int
direction_code(double x1, double y1, double x2, double y2)
{
    const static double divisor = (2.0 * M_PI) / NumDirections;
    
    double angle = std::atan2(y2 - y1, x2 - x1);
    if (angle < 0) {
        angle += M_PI + M_PI;
    }
    
    return static_cast<int>((angle /*+ (0.5 * divisor)*/) / divisor) % NumDirections;
}


static int *
extract_directions(const NormalizedStroke &stroke)
{
    int *structure = DEBUG_NEW int[num_points(stroke) - 1];
    if (!structure) {
        errval = E_OUTOFMEM;
        return 0;
    }
    
    int *struc = structure;
    
    double *x = stroke.x;
    double *y = stroke.y;
    double *endx = stroke.x + num_points(stroke) - 1;
    for (; x != endx; x++, y++) {
        *(struc++) = direction_code(*x, *y, *(x + 1), *(y + 1));
        //DEBUG_ONLY(debug_out << "(" << *x << "," << *y << ")->(" << *(x + 1) << "," << *(y + 1) << ") : " << *(struc - 1) << std::endl);
    }
        
    return structure;
}

/*
static std::vector<ChainCodeEntry> *
extract_structure(const NormalizedStroke &stroke)
{
    std::vector<ChainCodeEntry> *vec = DEBUG_NEW std::vector<ChainCodeEntry>;
    std::vector<ChainCodeEntry> &codes = *vec;
    
    int *directions = extract_directions(stroke);

    DEBUG_ONLY(
    debug_out << "original directions: ";
    for (unsigned i = 0; i < num_points(stroke) - 1; i++) {
        debug_out << directions[i] << " ";
    }
    debug_out << std::endl;
    )
    
    int old_direction = *directions;
    for (unsigned i = 1; i < num_points(stroke) - 1; ) {
        int direction = directions[i];
        
        if (num_points(stroke) - i <= 2) {
            codes.push_back(ChainCodeEntry(Dot, 0));
            break;
        }
        
        unsigned loop_end = stroke_self_intersection(stroke, i - 1, 10.0);
        if (loop_end < num_points(stroke)) {
            if (i <= 2) {
                codes.clear();
                codes.push_back(ChainCodeEntry(Dot, 0));
            }
            
            codes.push_back(ChainCodeEntry(Loop, 0));
            i = loop_end + 1;
            old_direction = directions[loop_end];
        }
        else if (direction == old_direction) {
            do {
                i++;
                if (i == num_points(stroke) - 1) {
                    break;
                }
                direction = directions[i];
            } while (direction == old_direction);
            
            codes.push_back(ChainCodeEntry(Line, old_direction));
        }
        else {        
            unsigned change = static_cast<unsigned>((direction - old_direction) % NumDirections);
            if (change <= NumDirections / 2) {
                // stroke is going counter-clockwise
                unsigned first_point = i - 1;
                unsigned nline = 0;
                int old_change;
                do {
                    i++;
                    if (i == num_points(stroke) - 1) {
                        break;
                    }
                    old_direction = direction;
                    old_change = change;
                    direction = directions[i];
                    change = (direction - old_direction) % NumDirections;
                } while ((change % NumDirections <= NumDirections / 2) || (change == 0 && nline++ <= 3));
                
                codes.push_back(ChainCodeEntry(AntiClockwiseCurve, direction_code(stroke.x[first_point], stroke.y[first_point], stroke.x[i - 1], stroke.y[i - 1])));
            }
            else {
                // stroke is going clockwise
                unsigned first_point = i - 1;
                unsigned nline = 0;
                int old_change;
                do {
                    i++;
                    if (i == num_points(stroke) - 1) {
                        break;
                    }
                    old_direction = direction;
                    old_change = change;
                    direction = directions[i];
                    change = (direction - old_direction) % NumDirections;
                } while ((change % NumDirections > NumDirections / 2) || (change == 0 && nline++ <= 3));
                
                codes.push_back(ChainCodeEntry(ClockwiseCurve, direction_code(stroke.x[first_point], stroke.y[first_point], stroke.x[i - 1], stroke.y[i - 1])));
            }
        }
    }
    
    delete[] directions;
    
    DEBUG_ONLY(
    debug_out << "extracted structure: {";
    for (std::vector<ChainCodeEntry>::const_iterator i = codes.begin(); i != codes.end(); ++i) {
        switch (i->primitive) {
        case Line:
            debug_out << " { Line, " << i->direction << " }";
            break;
        case ClockwiseCurve:
            debug_out << " { CW, " << i->direction << " }";
            break;
        case AntiClockwiseCurve:
            debug_out << " { ACW, " << i->direction << " }";
            break;
        case Loop:
            debug_out << " Loop";
            break;
        case Dot:
            debug_out << " Dot";;
            break;
        }
    }
    debug_out << " }" << std::endl;
    )
    
    return vec;
}
*/


static int *
extract_structure(const NormalizedStroke &stroke)
{
    int *structure = DEBUG_NEW int[num_points(stroke) - 1];
    if (!structure) {
        errval = E_OUTOFMEM;
        return 0;
    }
    
    int *struc = structure;
    
    double *angles = angles_along_stroke(stroke);
    if (!angles) {
        delete[] structure;
        return 0;
    }

    const static double divisor = (2.0 * M_PI) / NumDirections;
    
    for (const double *a = angles; a != angles + num_points(stroke) - 1; a++) {            
        int direction = static_cast<int>(*a / divisor);
        if (direction < 0) {
            direction += NumDirections;
        }
        
        *(struc++) = direction;
    }
    
    delete[] angles;
    
    return structure;
}


template <typename T, template <typename T> class G>
static int *
extract_structure(const G<T> &group)
{
    int *structure = DEBUG_NEW int[num_points(group) - 1];
    if (!structure) {
        errval = E_OUTOFMEM;
        return 0;
    }
    
    int *struc = structure;
    
    for (typename G<T>::const_iterator stroke = group.begin(); stroke != group.end(); ++stroke) {
        double *angles = angles_along_stroke(*stroke);
        if (!angles) {
            delete[] structure;
            return 0;
        }

        const static double divisor = (2.0 * M_PI) / NumDirections;
        
        for (const double *a = angles; a != angles + num_points(*stroke) - 1; a++) {            
            int direction = static_cast<int>(*a / divisor);
            if (direction < 0) {
                direction += NumDirections;
            }
            
            *(struc++) = direction;
        }
        
        delete[] angles;
        
        if (stroke != group.end() - 1) {
            T dx = (stroke + 1)->x[0] - stroke->x[num_points(*stroke) - 1];
            T dy = (stroke + 1)->y[0] - stroke->y[num_points(*stroke) - 1];
            
            double angle = std::atan2(dy, dx);
            
            int direction = static_cast<int>(angle / divisor);
            if (direction < 0) {
                direction += NumDirections;
            }
            
            *(struc++) = direction;
        }
    }
    
    return structure;
}



static double
compare_structures(const int *model, const int *input, unsigned n, int perturb)
{
	unsigned cost = 0;
	
	for (const int *m = model, *i = input; m != model + n; m++, i++) {
		int dist = std::abs(*m - (*i + perturb));
		cost += std::min(dist, static_cast<int>(NumDirections - dist));
	}
	
	return static_cast<double>(cost) / n;
}


static double
compare_structures_shifted(const int *model, const int *input, unsigned n, int perturb)
{
	if (n == 1) {
		return compare_structures(model, input, n, perturb);
	}
	
	double cost = compare_structures(model, input, n, perturb);
	for (unsigned i = 1; i <= std::min(MaxShifting, n); i++) {
	    cost = std::min(cost, compare_structures(model, input + i, n - i, perturb));
	    cost = std::min(cost, compare_structures(model + i, input, n - i, perturb));
	}
	
	return cost;
}


static int
compare_structure_elements(int a, int b)
{
    int diff = std::abs(a - b);
    return std::min(diff, 8 - diff);
}


static double
compare_chain_codes(const ChainCodeEntry &first, const ChainCodeEntry &second)
{
    static const double DirectionCost = 10.0;
    static const double CurveCost = 30.0;
    
    switch (first.primitive) {
    case Line:
        switch (second.primitive) {
        case Line:
            return DirectionCost * ((first.direction - second.direction) % NumDirections);
        case ClockwiseCurve:
        case AntiClockwiseCurve:
            return CurveCost + DirectionCost * ((first.direction - second.direction) % NumDirections);
        case Loop:
            return 5.0 * CurveCost;
        case Dot:
            return 2.0 * CurveCost;
        default:
            throw E_INVALID;
        }
    
    case ClockwiseCurve:
        switch (second.primitive) {
        case Line:
            return CurveCost + DirectionCost * ((first.direction - second.direction) % NumDirections);
        case ClockwiseCurve:
            return DirectionCost * ((first.direction - second.direction) % NumDirections);
        case AntiClockwiseCurve:
            return 4.0 * CurveCost + DirectionCost * ((first.direction - second.direction) % NumDirections);
        case Loop:
            return 3.0 * CurveCost;
        case Dot:
            return 2.0 * CurveCost;
        default:
            throw E_INVALID;
        }

    case AntiClockwiseCurve:
        switch (second.primitive) {
        case Line:
            return CurveCost + DirectionCost * ((first.direction - second.direction) % NumDirections);
        case ClockwiseCurve:
            return 4.0 * CurveCost + DirectionCost * ((first.direction - second.direction) % NumDirections);
        case AntiClockwiseCurve:
            return DirectionCost * ((first.direction - second.direction) % NumDirections);
        case Loop:
            return 3.0 * CurveCost;
        case Dot:
            return 2.0 * CurveCost;
        default:
            throw E_INVALID;
        }
    
    case Loop:
        switch (second.primitive) {
        case Line:
            return 5.0 * CurveCost;
        case ClockwiseCurve:
            return 4.0 * CurveCost;
        case AntiClockwiseCurve:
            return 4.0 * CurveCost;
        case Loop:
            return 0.0;
        case Dot:
            return 3.0 * CurveCost;
        default:
            throw E_INVALID;
        }
    
    case Dot:
        switch (second.primitive) {
        case Line:
            return 2.0 * CurveCost;
        case ClockwiseCurve:
            return 2.0 * CurveCost;
        case AntiClockwiseCurve:
            return 2.0 * CurveCost;
        case Loop:
            return 3.0 * CurveCost;
        case Dot:
            return 0.0;
        default:
            throw E_INVALID;
        }
    
    default:
        throw E_INVALID;    
    }
}


int 
chaincode_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, Match &match)
{
    if (num_strokes(input) != num_strokes(model)) {
        match.num_strokes = 0;
        return 0;
    }

    NormalizedStroke model_stroke = concat(model);
    NormalizedStroke input_stroke = concat(input);
    double worst_cost = 0.0;
    
    //NormalizedStrokeGroup::const_iterator i, j;
    //for (i = model.begin(), j = input.begin(); i != model.end(); ++i, ++j) {
        //const NormalizedStroke &model_stroke = *i;
        //const NormalizedStroke &input_stroke = *j;

        unsigned npoints;
        
        NormalizedStroke new_model;
        NormalizedStroke new_input;
        
        if (num_points(model_stroke) > num_points(input_stroke)) {
            new_model = copy(model_stroke);
            new_input = subdivide(input_stroke, num_points(model_stroke));
            npoints = num_points(model_stroke);
        }
        else if (num_points(input_stroke) > num_points(model_stroke)) {
            new_model = subdivide(model_stroke, num_points(input_stroke));
            new_input = copy(input_stroke);
            npoints = num_points(input_stroke);
        }
        else {
            new_model = copy(model_stroke);
            new_input = copy(input_stroke);
            npoints = num_points(input_stroke);
        }

        //std::vector<ChainCodeEntry> *model_structure = extract_structure(new_model);
        int *model_structure = extract_structure(new_model);
        if (!model_structure) {
            int e = errval;
            errval = 0;
            return e;
        }
        
        //std::vector<ChainCodeEntry> *input_structure = extract_structure(new_input);
        int *input_structure = extract_structure(new_input);
        if (!input_structure) {
            int e = errval;
            errval = 0;
            return e;
        }
        
	    double total_cost = std::numeric_limits<double>::infinity();
	    //elastic_match(model_structure->begin(), model_structure->end(), input_structure->begin(), input_structure->end(),
	    //              total_cost, &compare_chain_codes);
	    //elastic_match(model_structure, model_structure + npoints/*num_points(model_stroke)*/, input_structure, input_structure + npoints/*num_points(input_stroke)*/,
	    //              total_cost, &compare_structure_elements);
	    for (int i = -MaxPerturbation; i <= MaxPerturbation; i++) {
	        total_cost = std::min(total_cost, compare_structures_shifted(model_structure, input_structure, npoints - 1, i));
	    }
	    
        //total_cost /= input_structure->size();
        if (total_cost > worst_cost) {
            worst_cost = total_cost;
        }
        //worst_cost += total_cost;
        
        //delete model_structure;
        //delete input_structure;
        delete[] input_structure;
        delete[] model_structure;
    //}

/*
    double total_cost = 0.0;
    
    elastic_match(model_structure, model_structure + num_points(model) - 1,
                  input_structure, input_structure + num_points(input) - 1,
                  total_cost, std::ptr_fun(compare_structure_elements));
*/
    //double threshold = 4.0 * num_points(input);

    match.raw_score = std::max(worst_cost, 0.1);// / std::max(num_points(input), num_points(model));
    match.num_strokes = num_strokes(model);
        
    return 0;
}

}

