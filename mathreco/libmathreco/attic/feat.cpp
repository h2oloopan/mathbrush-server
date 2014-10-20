#include "feat.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "error.h"
#include "group.h"
#include "interp.h"
#include "norm.h"
#include "profile.h"
#include "recog.h"
#include "stroke.h"


namespace scg
{


const double StrokeFeatures::weights[StrokeFeatures::NumFeatures] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, /*0.75/*, 20.0, 20.0*/ };

static const double StrictThresholdSpec[StrokeFeatures::NumFeatures] = { 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 15.0, /*15.0/* 1.0, 1.0/*, 2.0*/ };
static const double TightThresholdSpec[StrokeFeatures::NumFeatures] = { 25.0, 25.0, 15.0, 15.0, 20.0, 20.0, 20.0, 20.0, 40.0, /*40.0/* 1.0, 1.0/*, 2.0*/ };
//static const double DefaultThresholdSpec[StrokeFeatures::NumFeatures] = { 40.0, 40.0, 25.0, 25.0, 30.0, 30.0, 30.0, 30.0, 50.0, /*50.0/*, 0.0, 0.0/*, 2.0*/ };
static const double DefaultThresholdSpec[StrokeFeatures::NumFeatures] = { 35.0, 35.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 30.0 };
//static const double DefaultThresholdSpec[StrokeFeatures::NumFeatures] = { 50.0, 50.0, 35.0, 35.0, 40.0, 40.0, 40.0, 40.0, 200.0/*, 0.0, 0.0/*, 2.0*/ };

const StrokeFeatures StrictFeatureThreshold(StrictThresholdSpec);
const StrokeFeatures TightFeatureThreshold(TightThresholdSpec);
const StrokeFeatures DefaultFeatureThreshold(DefaultThresholdSpec);


StrokeFeatures::StrokeFeatures() { }

StrokeFeatures::StrokeFeatures(const double f[NumFeatures])
  { std::copy(f, f + NumFeatures, features); }


double &
StrokeFeatures::operator[](size_t n)
{
    if (n < NumFeatures) {
        return features[n];
    }
    else {
        throw E_INVALID;
    }
}

const double &
StrokeFeatures::operator[](size_t n) const
{
    if (n < NumFeatures) {
        return features[n];
    }
    else {
        throw E_INVALID;
    }
}


bool
StrokeFeatures::operator<(const StrokeFeatures &rhs) const
{
    for (size_t i = 0; i < NumFeatures; ++i) {
        if (features[i] >= rhs.features[i]) {
            return false;
        }
    }
    return true;
}


bool
StrokeFeatures::operator<=(const StrokeFeatures &rhs) const
{
    for (size_t i = 0; i < NumFeatures; ++i) {
        if (features[i] > rhs.features[i]) {
            return false;
        }
    }
    return true;
}

bool
StrokeFeatures::operator>(const StrokeFeatures &rhs) const
{
    for (size_t i = 0; i < NumFeatures; ++i) {
        if (features[i] <= rhs.features[i]) {
            return false;
        }
    }
    return true;

}


bool
StrokeFeatures::operator>=(const StrokeFeatures &rhs) const
{
    for (size_t i = 0; i < NumFeatures; ++i) {
        if (features[i] < rhs.features[i]) {
            return false;
        }
    }
    return true;
}


StrokeFeatures
StrokeFeatures::operator-(const StrokeFeatures &rhs) const
{
    StrokeFeatures ret;
    for (size_t i = 0; i < NumFeatures; ++i) {
        ret.features[i] = features[i] - rhs.features[i];
    }
    return ret;
}

StrokeFeatures
abs(const StrokeFeatures &features)
{
    StrokeFeatures ret;
    for (size_t i = 0; i < StrokeFeatures::NumFeatures; ++i) {
        ret.features[i] = std::abs(features[i]);
    }
    
    return ret;
}


double
FeatureNorm(const StrokeFeatures &features)
{
    double norm = 0.0;
    /*if (features[StrokeFeatures::BoundingBoxRatio] < std::numeric_limits<double>::epsilon()) {
        for (size_t i = 0; i < StrokeFeatures::NumFeatures; ++i) {
            if (i != StrokeFeatures::BoundingBoxRatio) {
                norm += std::abs(features[i]);
            }
        }
    }
    else {*/
        for (size_t i = 0; i < StrokeFeatures::NumFeatures; ++i) {
            norm += std::abs(features[i] * StrokeFeatures::weights[i]);
        }
    //}
    return norm;    
}


int
ExtractFeatures(const NormalizedStroke &stroke, const NormalizedStrokeGroup &group, StrokeFeatures &features)
{
    int err = ExtractFeatures(stroke, features);
    if (FAILURE(err)) {
        return err;
    }
    
    Rect<double> bounds = bbox(group);
    if (width(bounds) > std::numeric_limits<double>::epsilon()) {
        features[StrokeFeatures::Width] *= NORM_MAX / width(bounds);
    }

    if (height(bounds) > std::numeric_limits<double>::epsilon()) {
        features[StrokeFeatures::Height] *= NORM_MAX / height(bounds);
    }

    return 0;
}

/*
static int
FindExtremePoints(const NormalizedStroke &stroke, std::vector<unsigned> &points)
{
    double *px = stroke.x;
    double *py = stroke.y;
    double *endx = stroke.x + stroke.npoints - 1;
    double *endy = stroke.y + stroke.npoints - 1;
    
}
*/

int
ExtractFeatures(const NormalizedStroke &stroke, StrokeFeatures &features)
{
    if (num_points(stroke) == 0) {
        return E_NOTFOUND;
    }
    
    Rect<double> bounds = bbox(stroke);
    
    features[StrokeFeatures::Top] = bounds.top;
    features[StrokeFeatures::Left] = bounds.left;
    features[StrokeFeatures::Width] = width(bounds);
    features[StrokeFeatures::Height] = height(bounds);
    
    //features[StrokeFeatures::Angle] = std::atan2(stroke.y[num_points(stroke) - 1] - stroke.y[0], stroke.x[num_points(stroke) - 1] - stroke.x[0]);
    //features[StrokeFeatures::DisplacementX] = stroke.x[num_points(stroke) - 1] - stroke.x[0];
    //features[StrokeFeatures::DisplacementY] = stroke.y[num_points(stroke) - 1] - stroke.y[0];
    features[StrokeFeatures::FirstX] = stroke.x[0];
    features[StrokeFeatures::FirstY] = stroke.y[0];
    features[StrokeFeatures::LastX] = stroke.x[num_points(stroke) - 1];
    features[StrokeFeatures::LastY] = stroke.y[num_points(stroke) - 1];
    /*
    if (features[StrokeFeatures::Width] < NORM_MAX / 10.0 || features[StrokeFeatures::Height] < NORM_MAX / 10.0) {
        features[StrokeFeatures::BoundingBoxRatio] = 0.0;
    }
    else {
        features[StrokeFeatures::BoundingBoxRatio] = width_height_ratio(bounds);
    }*/
    
    double length = 0.0;
    double *px, *py;
    double *endx = stroke.x + num_points(stroke) - 1;
    for (px = stroke.x, py = stroke.y; px < endx; px++, py++) {
        length += dist_sq(*px, *py, *(px + 1), *(py + 1));
    }
    features[StrokeFeatures::Length] = std::sqrt(length);
    /*
    double accel = 0.0;
    double dx1, dy1;
    double dx2, dy2;
    px = stroke.x;
    py = stroke.y;
    dx1 = px[1] - px[0];
    dy1 = py[1] - py[0];
    for (px = stroke.x + 2, py = stroke.y + 2; px < endx; px++, py++) {
        dx2 = *px - *(px - 1);
        dy2 = *py - *(py - 1);
        accel += std::abs(dx2 - dx1) + std::abs(dy2 - dy1);
        dx1 = dx2;
        dy1 = dy2;
    }
    features[StrokeFeatures::TotalAccel] = accel;
    */
    return 0;
}


double
FeatureBasedMatch(const StrokeFeatures &from, const StrokeFeatures &to)
{
    return FeatureNorm(from - to);
}


double
FeatureBasedMatch(const NormalizedStroke &from, const NormalizedStroke &to)
{
    int e;
    StrokeFeatures from_features;
    StrokeFeatures to_features;
    e = ExtractFeatures(from, from_features);
    if (FAILURE(e)) {
        return std::numeric_limits<double>::infinity();
    }
    e = ExtractFeatures(to, to_features);
    if (FAILURE(e)) {
        return std::numeric_limits<double>::infinity();
    }

    return FeatureBasedMatch(from_features, to_features);
}


int
ComputeProfileCorrelations(const Profile &from, const Profile &to, std::vector<InterProfileCorrelation> &correlations)
{
    static const double CorrelationFeatureThreshold = GetParameterDouble("CorrelationFeatureThreshold");//FeatureNorm(TightFeatureThreshold);
    static const double CorrelationConfidenceThreshold = GetParameterDouble("CorrelationConfidenceThreshold");
    
    for (Profile::const_iterator from_symbol = from.begin(); from_symbol != from.end(); ++from_symbol) {
        for (ProfileSymbol::const_iterator from_proto = from_symbol->begin(); from_proto != from_symbol->end(); ++from_proto) {
            std::vector<std::pair<PrototypeId, double> > proto_correlations;
            
            for (Profile::const_iterator to_symbol = to.begin(); to_symbol != to.end(); ++to_symbol) {
                for (ProfileSymbol::const_iterator to_proto = to_symbol->begin(); to_proto != to_symbol->end(); ++to_proto) {
                    if (from_proto != to_proto) {
                        if (num_strokes(from_proto->strokes) > num_strokes(to_proto->strokes)) {
									 NormalizedStrokeGroup group = copy(from_proto->strokes);
									 prepare_to_match(&group[0], &group[num_strokes(to_proto->strokes)], to_proto->strokes);
									 double feature_score = FeatureBasedMatch(group, to_proto->strokes);
                            feature_score = 1.0 - feature_score / CorrelationFeatureThreshold;
                            if (feature_score >= CorrelationConfidenceThreshold) {
                                correlations.push_back(InterProfileCorrelation(from_proto->id, to_proto->id, feature_score));
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}




double FeatureBasedMatch2(const NormalizedStrokeGroup &from_group, const NormalizedStroke &from,
                         const NormalizedStrokeGroup &to_group, const NormalizedStroke &to)
{
    int e;
    StrokeFeatures from_features;
    StrokeFeatures to_features;
    e = ExtractFeatures(from, from_group, from_features);
    if (FAILURE(e)) {
        return std::numeric_limits<double>::infinity();
    }
    e = ExtractFeatures(to, to_group, to_features);
    if (FAILURE(e)) {
        return std::numeric_limits<double>::infinity();
    }
    
    return FeatureBasedMatch(from_features, to_features);
}

double FeatureBasedMatch(const NormalizedStrokeGroup &from_group, const NormalizedStroke &from,
                         const NormalizedStrokeGroup &to_group, const NormalizedStroke &to)
{
    int e;
    StrokeFeatures from_features;
    StrokeFeatures to_features;
    e = ExtractFeatures(from, from_group, from_features);
    if (FAILURE(e)) {
        return std::numeric_limits<double>::infinity();
    }
    e = ExtractFeatures(to, to_group, to_features);
    if (FAILURE(e)) {
        return std::numeric_limits<double>::infinity();
    }
    
    return FeatureBasedMatch(from_features, to_features);
}


double
FeatureBasedMatch(const NormalizedStrokeGroup &from, const NormalizedStrokeGroup &to)
{
    if (num_strokes(from) != num_strokes(to)) {
        return std::numeric_limits<double>::infinity();
    }
    
    double score = 0.0;
    NormalizedStrokeGroup::const_iterator from_stroke;
    NormalizedStrokeGroup::const_iterator to_stroke;
    for (from_stroke = from.begin(), to_stroke = to.begin();
         from_stroke != from.end(); ++from_stroke, ++to_stroke) {
        score += FeatureBasedMatch(from, *from_stroke, to, *to_stroke);
    }
    
    return score / num_strokes(from);
}


int
FeatureBasedMatch(const NormalizedStrokeGroup &from, const NormalizedStrokeGroup &to, Match &match)
{
    match.num_strokes = 0;
    
    double score = FeatureBasedMatch(from, to);
    
    if (score != std::numeric_limits<double>::infinity()) {
        match.num_strokes = num_strokes(from);
        match.raw_score = score;
    }
    
    return 0;
}



std::ostream &
operator<<(std::ostream &os, const scg::StrokeFeatures &sf)
{
    os << "{ ";
    for (size_t i = 0; i < scg::StrokeFeatures::NumFeatures - 1; i++) {
        os << sf.features[i] << ", ";
    }
    os << sf.features[scg::StrokeFeatures::NumFeatures - 1] << " }";
    return os;
}

}
