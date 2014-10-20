#include "elastic.h"
#include "error.h"
#include "feat-old.h"
#include "group.h"
#include "interp.h"
#include "match.h"
#include "memory.h"
#include "merge.h"
#include "norm.h"
#include "parms.h"
#include "recog.h"
#include "stroke.h"
#include "stroke-alg.h"
#include "structural.h"


namespace scg
{


#if 0
NormalizedStrokeGroupTemp
merge(const NormalizedStrokeGroup &g1, const NormalizedStrokeGroup &g2, double first_weight)
{
    bool success = false;
    
    NormalizedStrokeGroup ret;
    
    if (num_strokes(g1) != num_strokes(g2)) {
        errval = E_INVALID;
        return NormalizedStrokeGroupTemp();
    }
    
    //NormalizedStroke s1 = concat(g1.strokes(), num_strokes(g1));
    //NormalizedStroke s2 = concat(g2.strokes(), num_strokes(g2));

    double ratio1 = height_width_ratio(bbox(g1));
    double ratio2 = height_width_ratio(bbox(g2));
    
    if (std::max(ratio1 / ratio2, ratio2 / ratio1) >= 1.25) {
        errval = E_INVALID;
        return NormalizedStrokeGroupTemp();
    }
    
    NormalizedStroke *strokes = 0;
    //unsigned **match_points = 0;
    
    const NormalizedStrokeGroup *group1;
    NormalizedStrokeGroup group2;

    if (num_points(g1) > num_points(g2)) {
        group1 = &g2;
        group2 = reorder_strokes(g1, g2);
        if (!num_strokes(group2)) {
            errval = E_INVALID;
            return NormalizedStrokeGroupTemp();
        }
    }
    else {
        group1 = &g1;
        group2 = reorder_strokes(g2, g1);
        if (!num_strokes(group2)) {
            errval = E_INVALID;
            return NormalizedStrokeGroupTemp();
        }
    }
    
    Match match;
    elastic_match(g1, g2, match/*, 0, &match_points*/);
    if (!match.num_strokes || match.raw_score > 50.0) {
        errval = E_INVALID;
        goto clean_up;
    }
    
    double second_weight = 1.0 - first_weight;
    
    strokes = DEBUG_NEW NormalizedStroke[num_strokes(g1)];
    if (!strokes) {
        errval = E_OUTOFMEM;
        goto clean_up;
    }
    
    unsigned **matches;
    NormalizedStroke *s;
    NormalizedStrokeGroup::const_iterator stroke1, stroke2;
    for (matches = match_points, s = strokes, stroke1 = g1.begin(), stroke2 = g2.begin();
         stroke1 != g1.end();
         matches++, s++, stroke1++, stroke2++)
    {
        unsigned npoints = std::max(num_points(*stroke1), num_points(*stroke2));
        
        NormalizedStroke sub1 = subdivide(*stroke1, npoints);
        NormalizedStroke sub2 = subdivide(*stroke2, npoints);
        
        double cent_dx = x_center(*stroke2) - x_center(*stroke1);
        double cent_dy = y_center(*stroke2) - y_center(*stroke1);
        
        double *x = DEBUG_NEW double[num_points(*stroke2)];
        if (!x) {
            errval = E_OUTOFMEM;
            goto clean_up;
        }
        
        double *y = DEBUG_NEW double[num_points(*stroke2)];
        if (!y) {
            delete[] x;
            errval = E_OUTOFMEM;
            goto clean_up;
        }

        if (*matches) {        
            for (unsigned i = 0; i < num_points(*stroke2); i++) {
                x[i] = (stroke1->x[(*matches)[i]] + cent_dx) * first_weight + stroke2->x[i] * second_weight;
                y[i] = (stroke1->y[(*matches)[i]] + cent_dy) * first_weight + stroke2->y[i] * second_weight;
            }
        }
        else {
            double first_points_per_second = static_cast<double>(num_points(*stroke1)) / num_points(*stroke2);
            double curr_pt = 0.0;
            for (unsigned i = 0; i < num_points(*stroke2); i++) {
                x[i] = (stroke1->x[round<unsigned>(curr_pt)] + cent_dx) * first_weight + stroke2->x[i] * second_weight;
                y[i] = (stroke1->y[round<unsigned>(curr_pt)] + cent_dy) * first_weight + stroke2->y[i] * second_weight;
                
                curr_pt += first_points_per_second;
            }
        }
        
        s->set_points(x, y, num_points(*stroke2));
    }

    ret.set_strokes(strokes, num_strokes(g2));
    
    elastic_match(g1, ret, match);
    if (!match.num_strokes || match.raw_score > 50.0) {
        errval = E_INVALID;
    }
    else {
        elastic_match(g2, ret, match);
        if (!match.num_strokes || match.raw_score > 50.0) {
            errval = E_INVALID;
        }
        else {
            success = true;
        }
    }
    
clean_up:
    /*for (unsigned **points = match_points; points != match_points + num_strokes(group2); points++) {
        delete[] *points;
    }
    delete[] match_points;
*/

    if (success) {
        return normalize(ret);
    }
    else {
        if (!ret.strokes()) {
            delete[] strokes;
        }
        
        return NormalizedStrokeGroupTemp();
    }
}
#endif


NormalizedStrokeGroup
merge(const NormalizedStrokeGroup &g1, const NormalizedStrokeGroup &g2, double first_weight)
{
    static const double MergeBoxRatioThreshold = GetParameterDouble("MergeBoxRatioThreshold");
    static const double MergeMatchScoreThreshold = GetParameterDouble("MergeMatchScoreThreshold");
    
    bool success = false;
    
    NormalizedStrokeGroup ret;
    
    if (num_strokes(g1) != num_strokes(g2)) {
        errval = E_INVALID;
        return NormalizedStrokeGroup();
    }
    
    //NormalizedStroke s1 = concat(g1.strokes(), num_strokes(g1));
    //NormalizedStroke s2 = concat(g2.strokes(), num_strokes(g2));

    double ratio1 = height_width_ratio(bbox(g1));
    double ratio2 = height_width_ratio(bbox(g2));
    
    if (std::max(ratio1 / ratio2, ratio2 / ratio1) >= MergeBoxRatioThreshold) {
        errval = E_INVALID;
        return NormalizedStrokeGroup();
    }
    
    Match match;
    NormalizedStroke *strokes = 0;

    double second_weight = 1.0 - first_weight;
    
    double score = FeatureBasedMatch(g1, g2);
    if (score > 65.0) {
        errval = E_INVALID;
        goto clean_up;
    }
    
    ElasticMatch(g1, g2, match);
    if (!match.num_strokes || match.raw_score > MergeMatchScoreThreshold) {
        errval = E_INVALID;
        goto clean_up;
    }

/*
    Match match;
    structural_match(g1, g2, match);
    if (!match.num_strokes || match.raw_score > 250.0) {
        errval = E_INVALID;
        goto clean_up;
    }
*/
/*
    for (NormalizedStrokeGroup::const_iterator stroke1 = g1.begin(), stroke2 = g2.begin(); stroke1 != g1.end(); ++stroke1, ++stroke2) {
        StrokeFeatures features1;
        StrokeFeatures features2;
        ExtractFeatures(*stroke1, features1);
        ExtractFeatures(*stroke2, features2);
        
        if (abs(features1 - features2) > StrictFeatureThreshold) {
            errval = E_INVALID;
            goto clean_up;
        }
    }*/
        
    strokes = DEBUG_NEW NormalizedStroke[num_strokes(g1)];
    if (!strokes) {
        errval = E_OUTOFMEM;
        goto clean_up;
    }
    
    { // scope this loop to eliminate warnings about the goto's above
      // skipping the initialization of these variables
    NormalizedStroke *s;
    NormalizedStrokeGroup::const_iterator stroke1, stroke2;
    for (s = strokes, stroke1 = g1.begin(), stroke2 = g2.begin(); stroke1 != g1.end(); s++, stroke1++, stroke2++) {
        unsigned npoints = std::max(num_points(*stroke1), num_points(*stroke2));
        
        NormalizedStroke sub1 = subdivide(*stroke1, npoints);
        NormalizedStroke sub2 = subdivide(*stroke2, npoints);
        
        double cent_dx = x_center(*stroke2) - x_center(*stroke1);
        double cent_dy = y_center(*stroke2) - y_center(*stroke1);
        
        double *x = DEBUG_NEW double[npoints];
        if (!x) {
            errval = E_OUTOFMEM;
            goto clean_up;
        }
        
        double *y = DEBUG_NEW double[npoints];
        if (!y) {
            delete[] x;
            errval = E_OUTOFMEM;
            goto clean_up;
        }

        for (unsigned i = 0; i < npoints; i++) {
            x[i] = (sub1.x[i] + cent_dx) * first_weight + sub2.x[i] * second_weight;
            y[i] = (sub1.y[i] + cent_dy) * first_weight + sub2.y[i] * second_weight;
        }
        
        s->set_points(x, y, npoints);
    }

    } // end of loop scope
    
    ret.set_strokes(strokes, num_strokes(g2));
    /*
    elastic_match(g1, ret, match);
    if (!match.num_strokes || match.raw_score > 50.0) {
        errval = E_INVALID;
    }
    else {
        elastic_match(g2, ret, match);
        if (!match.num_strokes || match.raw_score > 50.0) {
            errval = E_INVALID;
        }
        else {
            success = true;
        }
    }
    */
    success = true;
    
clean_up:
    if (success) {
        return ret;
    }
    else {
        if (!ret.strokes) {
            delete[] strokes;
        }
        
        return NormalizedStrokeGroup();
    }
}


NormalizedStrokeGroup
merge(const NormalizedStrokeGroup &g1, const NormalizedStrokeGroup &g2)
{
    return merge(g1, g2, 0.5);
}


}

