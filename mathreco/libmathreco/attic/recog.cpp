/*
 * recog.cpp
 *
 * This file implements the algorithms which control the process of recognizing
 * ink.  Much of the file consists of experimental functions which are commented
 * out pending further examination.
 */
 
#include <algorithm>
#include <cmath>
#include <fstream>
#include <istream>
#include <list>
#include <map>
#include <sstream>
#include <vector>
#include <iostream>

#include "debug.h"
#include "deform.h"
#include "ink-io.h"
#include "direlem.h"
#include "distrib.h"
#include "elastic.h"
#include "error.h"
#include "feat.h"
#include "group.h"
#include "interp.h"
#include "log.h"
#include "match.h"
#include "memory.h"
#include "norm.h"
#include "parms.h"
#include "profile.h"
#include "profman.h"
#include "recog.h"
#include "segment.h"
#include "stroke.h"
#include "stroke-alg.h"
#include "structural.h"
#include "vector.h"
#include "verb.h"


namespace scg
{

const unsigned RecognitionSegment::INVALID_ID = static_cast<unsigned>(-1);
unsigned RecognitionSegment::next_id = 0;


static const unsigned MaxGroupSize = GetParameterUnsigned("MaxGroupSize");


static const double PenliftStabilityRatio = GetParameterDouble("PenliftStabilityRatio");

typedef int (*MatcherFunction)(const NormalizedStrokeGroup &, const NormalizedStrokeGroup &, Match &);

static MatcherFunction matcher_functions[NumMatchers] = {0};

enum Segmenters
{
    HeuristicGuesser,  // heuristic method
    HeuristicScorer,   // This method doesn't work.
    FeatureCorrelator, // feature-extraction method
    TimingAnalyzer,    // This gathers statistics from the user on the fly re: elapsed time between strokes
    NumSegmenters
};

enum SegmenterBits
{
    HeuristicGuesserBit  = (1 << HeuristicGuesser),
    HeuristicScorerBit   = (1 << HeuristicScorer),
    FeatureCorrelatorBit = (1 << FeatureCorrelator),
    TimingAnalyzerBit    = (1 << TimingAnalyzer)
};


// These are the matchers whose results will be combined:
int DefaultMatchers = ElasticMatcherBit;

// These are the segmentation routines whose results will be combined:
int DefaultSegmenters = HeuristicGuesserBit | FeatureCorrelatorBit | TimingAnalyzerBit;


static const double FeaturePruningThreshold = GetParameterDouble("FeaturePruningThreshold");
static const double FeatureConsiderationThreshold = GetParameterDouble("FeatureConsiderationThreshold");
static const double FeatureDistributionCutoff = GetParameterDouble("FeatureDistributionCutoff");

struct SegmentTimingData
{
    unsigned curr_index;
    std::vector<unsigned> elapsed;
};


struct SegmentRecognitionData
{
    std::vector<double> segment_probabilities;
    std::vector<double> stacked_segment_probabilities;
    std::vector<double> noheuristic_segment_probabilities;
    std::map<unsigned, NormalizedStrokeGroup> cached_input;
    std::vector<const Prototype *> prototypes;
};


struct RecognitionCursor
{
	const RecognitionContext *ctx;
	SegmentRecognitionData segment_data;
	unsigned start_stroke;
	unsigned nstrokes;

	std::map<PrototypeId, RecognitionResult> matches;
};



struct RecognitionContextImpl
{
	RawStrokeGroup input; // stroke group from RecognitionContext after joining and sorting
	SegmentTimingData timing;
};

RecognitionContext::~RecognitionContext()
{
	delete impl;
}


static ExponentialDistribution time_between_strokes(PenliftStabilityRatio);

double matcher_weights[NumMatchers] = {2.0, 0.0, 0.0, 0.0, 1.0};

NormalizedStrokeGroup
GetInputStrokes(const RawStroke *strokes, unsigned nstrokes, unsigned max_strokes)
{
    // this function performs preprocessing on input strokes and extracts a group
    // consisting of 'nstrokes' strokes.  It resamples and normalizes the stroke group,
    // trimming sharp curves off the ends and smoothing the strokes.
    NormalizedStrokeGroup input;
    ImmutableRawStrokeSubGroup selected_raw_input(strokes, nstrokes);
    
    input = normalize(selected_raw_input);
    input = smooth(trim_ends(input));
	 input = subdivide(input);

    return input;
}


// Rescale "group" to fit into "bounds".
static NormalizedStrokeGroup
rescale(const NormalizedStrokeGroup &group, const Rect<double> &bounds)
{
    NormalizedStroke *output = DEBUG_NEW NormalizedStroke[num_strokes(group)];
    if (!output) {
        return NormalizedStrokeGroup();
    }
    
    const Rect<double> &src_bounds = bbox(group);
    double src_height = height(src_bounds);
    double src_width = width(src_bounds);
    double dst_height = height(bounds);
    double dst_width = width(bounds);

    unsigned curr_stroke = 0;
    for (NormalizedStrokeGroup::const_iterator i = group.begin(); i != group.end(); ++i) {
        const NormalizedStroke &stroke = *i;
        unsigned npoints = num_points(stroke);
        
        double *x = DEBUG_NEW double[npoints];
        if (!x) {
            delete[] output;
            return NormalizedStrokeGroup();
        }
        
        double *y = DEBUG_NEW double[npoints];
        if (!y) {
            delete[] x;
            delete[] output;
            return NormalizedStrokeGroup();
        }
        
        const double *sx = stroke.x;
        const double *sy = stroke.y;
        double *px = x;
        double *py = y;
        const double *endx = stroke.x + npoints;
        for (; sx != endx; sx++, sy++, px++, py++) {
            *px = bounds.left + ((*sx / src_width) * dst_width);
            *py = bounds.top + ((*sy / src_height) * dst_height);
        }
        
        output[curr_stroke++].set_points(x, y, npoints);
    }
    
    return NormalizedStrokeGroup(output, num_strokes(group));
}


// This function reorders and reverses the strokes in "input_", as appropriate, so as to make
// the best match with "model".  If "prune_enabled" is true, the function will return an empty
// stroke group if the difference in any feature between the groups is too large.
NormalizedStrokeGroup
PrepareInputForMatch(const NormalizedStrokeGroup &input_, const NormalizedStrokeGroup &model, bool prune_enabled)
{
    if (num_strokes(input_) != num_strokes(model)) {
        return NormalizedStrokeGroup();
    }
    
    // We must copy the input in case we reverse strokes.
    NormalizedStrokeGroup input;
    if (num_strokes(model) > 1) {
        input = rescale(input_, bbox(model));
    }
    else {
        input = copy(input_);
    }
    
    bool *used_input = DEBUG_NEW bool[num_strokes(input)];
    std::fill(used_input, used_input + num_strokes(input), false);
    
    NormalizedStrokeGroup reordered(DEBUG_NEW NormalizedStroke[num_strokes(input)], num_strokes(input));
    
    // Iterate through the model strokes.  The group returned will have strokes in the same "order" as the model.
    for (NormalizedStrokeGroup::const_iterator model_stroke = model.begin(); model_stroke != model.end(); ++model_stroke) {
        double best_score = std::numeric_limits<double>::infinity();
        unsigned best_index;
        bool best_reversed;
        
        StrokeFeatures model_features;
        int error = ExtractFeatures(*model_stroke, model_features);
        if (FAILURE(error)) {
            return NormalizedStrokeGroup();
        }
               
        double model_dx = model_features[StrokeFeatures::LastX] - model_features[StrokeFeatures::FirstX];
        double model_dy = model_features[StrokeFeatures::LastY] - model_features[StrokeFeatures::FirstY];
        
        // This inner loop finds the best-matching stroke in the input.
        unsigned stroke_index = 0;
        for (NormalizedStrokeGroup::const_iterator input_stroke = input.begin(); input_stroke != input.end(); ++input_stroke) {
            if (!used_input[stroke_index]) {                
                StrokeFeatures input_features;
                error = ExtractFeatures(*input_stroke, input_features);
                if (FAILURE(error)) {
                    return NormalizedStrokeGroup();
                }
                
                
                VERBOSE2(
                   *verb_out << "   matching model " << model_features << std::endl;
                   *verb_out << "   against  input " << input_features << std::endl;
                );
                

                bool reversed = false;
                
                double input_dx = input_features[StrokeFeatures::LastX] - input_features[StrokeFeatures::FirstX];
                double input_dy = input_features[StrokeFeatures::LastY] - input_features[StrokeFeatures::FirstY];
                
                double dot = input_dx * model_dx
                           + input_dy * model_dy;
                // If the stroke is a loop, don't compare direction from start- to end-point since they will
                // be close to eachother and dot < 0 doesn't indicate that one is backwards.
                if (!stroke_is_closed(*input_stroke) && (dot < 0)) {
                    //input_features.x_displacement = -input_features.x_displacement;
                    //input_features.y_displacement = -input_features.y_displacement;
                    std::swap(input_features[StrokeFeatures::FirstX], input_features[StrokeFeatures::LastX]);
                    std::swap(input_features[StrokeFeatures::FirstY], input_features[StrokeFeatures::LastY]);
                    reversed = true;
                }
                
                StrokeFeatures feature_diff = model_features - input_features;

                if (!prune_enabled || (FeatureNorm/*abs*/(feature_diff) <= FeatureConsiderationThreshold)) {//FeatureNorm(DefaultFeatureThreshold))) {
                    double feature_score = FeatureNorm(feature_diff);
                    if (feature_score < best_score) {
                        best_score = feature_score;
                        best_index = stroke_index;
                        best_reversed = reversed;
                    }
                }
            }
            
            ++stroke_index;
        }
        
        // If we found a match, copy it to the output stroke group.
        if (best_score < std::numeric_limits<double>::infinity()) {
            NormalizedStroke best_input;
            
            if (best_reversed) {
                best_input = reverse_copy(input[best_index]);
            }
            else {
                best_input = copy(input[best_index]);
            }

            reordered[model_stroke - model.begin()] = best_input;
            used_input[best_index] = true;
        }
        else {
            delete[] used_input;
            return NormalizedStrokeGroup();
        }
    }
    
    delete[] used_input;
    
    return reordered;
}


int
MatchSymbol(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, RecognitionResult &match)
{
    return MatchSymbol(model, input, match, DefaultMatchers);
}

int
MatchSymbol(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, RecognitionResult &match, int matchers)
{
    int e;
    
    NormalizedStrokeGroup reordered_input = PrepareInputForMatch(input, model, false);
    if (!reordered_input.empty()) {
        const Rect<double> &bounds = bbox(input);
        
        for (unsigned i = 0; i < NumMatchers; i++) {
            if (matchers & (1 << i)) {
                Match result;
                e = (*(matcher_functions[i]))(model, reordered_input, result);
                
                VERBOSE2(*verb_out << e << ", " << result.num_strokes << ", " << result.raw_score << std::endl);
                
                if (SUCCESS(e) && result.num_strokes > 0 && result.raw_score < std::numeric_limits<double>::infinity()) {
                    if (!match.valid[i] || result.raw_score < match.matcher_scores[i]) {
                        match.matcher_scores[i] = result.raw_score;
                        match.valid[i] = true;
                    }
                }
            }
        }
    }

    return 0;
}


static int
ComputeRelativeSegmentataionScores(double *scores, unsigned count)
{
    double worst = scg::max<double, double *>(scores, scores + count);
    while (count--) {
        if (*scores) {
            *scores = worst / *scores;
        }
        scores++;
    }
    
    return 0;
}


// This function implements the adjustment of feature-correlation scores for segmentation
// based on symbol correlation, as described in the Maple Conference poster.
static int
BoostSegmentationScores(const ProfileManager &profman, const std::map<PrototypeId, double> &scores, std::map<unsigned, double> &scores_out)
{
    for (std::map<PrototypeId, double>::const_iterator i = scores.begin(); i != scores.end(); ++i) {
        const Prototype *prototype = profman.GetPrototypeById(i->first);

        if (num_strokes(prototype->strokes) > 1 && i->second > 0.0) {
            double correlated_scores = 0.0;
            double avg_correlation = 0.0;
            unsigned num_correlations = 0;
            if (i->second > 0.0) {
                PrototypeCorrelations correlations = profman.GetPrototypeCorrelationsFor(prototype->id);
                for (PrototypeCorrelations::const_iterator j = correlations.begin(); j != correlations.end(); ++j) {
                    std::map<PrototypeId, double>::const_iterator correlation = scores.find(j->id);
                    if (correlation != scores.end()) {
                        const Prototype *correlated_prototype = profman.GetPrototypeById(j->id);
                        correlated_scores += correlation->second * j->correlation;
                        num_correlations++;
                    }
                }
            }
                
            double weighted_correlated_score = i->second;
            if (num_correlations) {
                weighted_correlated_score += (/*i->second * */correlated_scores) / num_correlations;
            }
            
            if (weighted_correlated_score > scores_out[num_strokes(prototype->strokes)]) {
                scores_out[num_strokes(prototype->strokes)] = weighted_correlated_score;
            }
        }
    }

    return 0;
}


// This function implements the feature-correlation segmentation method
// described in the Maple Conference poster.  It also contains (commented-out)
// experimental code that models each feature comparison of each symbol
// as an exponential distribution which is updated on the fly as the user writes.
static int
ComputeCorrelationSegmentScores(RecognitionCursor &cur,
                                double *scores_out, unsigned max_segments)
{

    std::map<PrototypeId, double> prototype_score;
    std::map<PrototypeId, double> reco_score;
    
    std::map<unsigned, double> prob;
    std::map<unsigned, double> reco_prob;
    
    std::map<double, const Prototype *> prototypes;
    
    unsigned discarded = 0;
    unsigned n = 0;
    // First, simple feature comparison is used to prune the symbol database and give
    // base probabilities.
    for (ProfileManager::const_iterator i = cur.ctx->profman->begin(); i != cur.ctx->profman->end(); ++i) {
        const ProfileSymbol &symbol = *i;
        n++;

        for (ProfileSymbol::const_iterator j = symbol.begin(); j != symbol.end(); ++j) {
            const Prototype &prototype = *j;

            unsigned model_strokes = num_strokes(prototype.strokes);
            
            if (model_strokes <= cur.nstrokes) {                
                NormalizedStrokeGroup &input = cur.segment_data.cached_input[model_strokes];
                if (input.empty()) {
                    input = GetInputStrokes(cur.ctx->impl->input.begin() + cur.start_stroke, model_strokes, cur.nstrokes);
                }
                
                const NormalizedStrokeGroup &model = prototype.strokes;

                NormalizedStrokeGroup reordered_input = PrepareInputForMatch(input, model, false);
                if (!reordered_input.empty()) {
                    double score;
                    score = FeatureBasedMatch(reordered_input, model);
                    if (score < FeatureConsiderationThreshold) {
                        score = 1.0 - (score / FeatureConsiderationThreshold);
                        VERBOSE2(*verb_out << "feature score for " << cur.ctx->profman->GetInfoBySymbol(symbol).name() << "/" << prototype.id << " = " << score << std::endl);
                        cur.segment_data.prototypes.push_back(&prototype);
                        prototype_score[j->id] = score;
                        
                        if (score > prob[model_strokes]) {
                            prob[model_strokes] = score;
                        }
                    }
                    else {
                        discarded++;
                     }
                }
            }
        }
    }
    VERBOSE2(*verb_out << "ready to boost; discarded " << discarded << " prototypes" << std::endl);
    
    // Now the probabilities of multi-stroke symbols are boosted by the probabilities of fewer-stroke symbols,
    // weighted by how much the smaller symbols resemble the first few strokes of the larger ones.
    BoostSegmentationScores(*cur.ctx->profman, prototype_score, prob);

    double baseline = 0.0;
    for (std::map<unsigned, double>::const_iterator i = prob.begin(); i != prob.end(); i++) {
        baseline = std::min(baseline, i->second);
    }
    VERBOSE2(
	 	*verb_out << "feature segmentation scores:" << std::endl;
    	for (std::map<unsigned, double>::iterator i = prob.begin(); i != prob.end(); i++) {
        *verb_out << "  " << i->first << ": " << i->second << std::endl;
    	}
	);

    for (unsigned i = 0; i < max_segments; i++) {
        scores_out[i] = prob[i + 1];//(prob[i + 1] + reco_prob[i + 1]) / 2.0;
    }
    
    return 0;
}


static double
HeuristicGuessWeight(SegmentationConfidence confidence)
{
    const static double HighHeuristicWeight = GetParameterDouble("HighHeuristicWeight");
    const static double MediumHeuristicWeight = GetParameterDouble("MediumHeuristicWeight");
    const static double LowHeuristicWeight = GetParameterDouble("LowHeuristicWeight");
    
    switch (confidence) {
    case High:
        return HighHeuristicWeight;
    case Medium:
        return MediumHeuristicWeight;
    case Low:
        return LowHeuristicWeight;
    }
    
    return 0.0;
    
}

template <typename T>
static bool
GroupContainsDot(const T &group)
{
    for (typename T::const_iterator i = group.begin(); i != group.end(); ++i) {
        if (stroke_is_dot(*i)) {
            return true;
        }
    }
    return false;
}


static int
ComputeSegmentationProbabilities(RecognitionCursor &cur, int segmenters)
{
    // This function computes the probability that the next input symbol in 'strokes' has n
    // strokes, for all n.  These probabilities are put into 'unboosted_prob'.  'prob' is the
    // same as 'unboosted_prob', except the results of the segmentation heuristics are used
    // to reweight probabilities.
    
    cur.segment_data.segment_probabilities.reserve(MaxGroupSize);
    cur.segment_data.segment_probabilities.resize(MaxGroupSize, 0.0);

    cur.segment_data.stacked_segment_probabilities.reserve(MaxGroupSize);
    cur.segment_data.stacked_segment_probabilities.resize(MaxGroupSize, 0.0);

    cur.segment_data.noheuristic_segment_probabilities.reserve(MaxGroupSize);
    cur.segment_data.noheuristic_segment_probabilities.resize(MaxGroupSize, 0.0);
    
    
    /*if (nstrokes > 1 && stroke_is_dot(strokes[1])) {
        segment_data.segment_probabilities[1] = 1.0;
        segment_data.stacked_segment_probabilities[1] = 1.0;
        segment_data.noheuristic_segment_probabilities[1] = 1.0;

        for (ProfileManager::const_iterator i = profman.begin(); i != profman.end(); ++i) {
            const ProfileSymbol &symbol = *i;
            
            for (ProfileSymbol::const_iterator j = symbol.begin(); j != symbol.end(); ++j) {
                if (GroupContainsDot(j->strokes)) {
                    segment_data.prototypes.push_back(&(*j));
                }
            }
        }

        return 0;
    }*/
    
    SegmentationGuess heuristic_guess;
    if (segmenters & HeuristicGuesserBit) {
        heuristic_guess = segment_group_size(cur.ctx->impl->input.begin() + cur.start_stroke, cur.nstrokes);

        heuristic_guess.guess.num_strokes = std::max<unsigned>(1, std::min(heuristic_guess.guess.num_strokes, MaxGroupSize));
        heuristic_guess.stacked_guess.num_strokes = std::max<unsigned>(1, std::min(heuristic_guess.stacked_guess.num_strokes, MaxGroupSize));
    
        //segmenter_scores[HeuristicGuesser][heuristic_guess.guess - 1] = 1.0;
        
        cur.segment_data.segment_probabilities[heuristic_guess.guess.num_strokes - 1] += HeuristicGuessWeight(heuristic_guess.guess.confidence);
        cur.segment_data.stacked_segment_probabilities[heuristic_guess.stacked_guess.num_strokes - 1] += HeuristicGuessWeight(heuristic_guess.stacked_guess.confidence);
    }
        

    if (segmenters & FeatureCorrelatorBit) {
        double *correlation_scores = DEBUG_NEW double[MaxGroupSize];
        ComputeCorrelationSegmentScores(cur, correlation_scores, MaxGroupSize);
        
        for (unsigned i = 0; i < MaxGroupSize; i++) {
            cur.segment_data.segment_probabilities[i] += correlation_scores[i];
            //segment_data.stacked_segment_probabilities[i] += correlation_scores[i];
            cur.segment_data.noheuristic_segment_probabilities[i] += correlation_scores[i];
        }
        
        delete[] correlation_scores;
    }
    else {
        cur.segment_data.noheuristic_segment_probabilities = cur.segment_data.segment_probabilities;
        
        for (ProfileManager::const_iterator i = cur.ctx->profman->begin(); i != cur.ctx->profman->end(); ++i) {
            const ProfileSymbol &symbol = *i;
            
            for (ProfileSymbol::const_iterator j = symbol.begin(); j != symbol.end(); ++j) {
                if ((cur.segment_data.segment_probabilities[num_strokes(j->strokes) - 1] > 0.0)
					  && (num_strokes(j->strokes) <= cur.nstrokes)) {
                    cur.segment_data.prototypes.push_back(&(*j));
                }
            }
        }
    }
    
    // This bit implements the exponential decay model for elapsed time between strokes.
    // Note that this method cannot _reduce_ the hints from other methods because the user
    // may pause for any length of time between strokes.
    VERBOSE2(*verb_out << "settled? " << time_between_strokes.IsSettled() << std::endl);
    if ((segmenters & TimingAnalyzerBit) && !cur.ctx->impl->timing.elapsed.empty() && time_between_strokes.IsSettled()) {
        static const double TimingWeight = GetParameterDouble("TimingWeight");
        static std::vector<unsigned> input_strokes;
        
        input_strokes.clear();
        input_strokes.push_back(cur.ctx->segments[cur.ctx->impl->timing.curr_index]->input_order.front());
        
        VERBOSE2(*verb_out << "timing analyzer confidence here:" << std::endl);

        //double lambda = 1.0 / average_penlift_time;
        for (size_t n = 1; n < std::min(cur.nstrokes, MaxGroupSize); n++) {
            unsigned next_stroke = cur.ctx->segments[cur.ctx->impl->timing.curr_index + n]->input_order.front();
            input_strokes.push_back(next_stroke);
            
            std::vector<unsigned>::const_iterator j;
            for (j = input_strokes.begin(); j != input_strokes.end(); ++j) {
                if (std::find(input_strokes.begin(), input_strokes.end(), *j + 1) == input_strokes.end()
                 && std::find(input_strokes.begin(), input_strokes.end(), *j - 1) == input_strokes.end()) {
                    break;
                }
            }
            
            if (j == input_strokes.end()) {
                std::sort(input_strokes.begin(), input_strokes.end());
                double prob = 1.0;
                for (j = input_strokes.begin(); j < input_strokes.end() - 1; ++j) {
                    prob *= 1.0 - time_between_strokes.Cdf(cur.ctx->impl->timing.elapsed[*j]);//std::exp(-lambda * timing_data.elapsed[*j]);
                }
                VERBOSE2(*verb_out << "  " << n + 1 << " strokes: " << prob << std::endl);
                if (cur.segment_data.noheuristic_segment_probabilities[n] > 0.0) {
                    cur.segment_data.noheuristic_segment_probabilities[n] += prob * TimingWeight;
                }
                if (cur.segment_data.noheuristic_segment_probabilities[n] > 0.0) {
                    cur.segment_data.segment_probabilities[n] += prob * TimingWeight;
                }
            }
        }    
    }
    
    double total = std::accumulate(cur.segment_data.segment_probabilities.begin(), cur.segment_data.segment_probabilities.end(), 0.0);
    for (unsigned i = 0; i < MaxGroupSize; i++) {
        cur.segment_data.segment_probabilities[i] /= total;
    }

    total = std::accumulate(cur.segment_data.stacked_segment_probabilities.begin(), cur.segment_data.stacked_segment_probabilities.end(), 0.0);
    for (unsigned i = 0; i < MaxGroupSize; i++) {
        cur.segment_data.stacked_segment_probabilities[i] /= total;
    }
    
    total = std::accumulate(cur.segment_data.noheuristic_segment_probabilities.begin(), cur.segment_data.noheuristic_segment_probabilities.end(), 0.0);
    for (unsigned i = 0; i < MaxGroupSize; i++) {
        cur.segment_data.noheuristic_segment_probabilities[i] /= total;
    }
    
    return 0;
}


static int
MatchSymbols(RecognitionCursor &cur, int matchers = DefaultMatchers)
{
    int e;
    
    for (std::vector<const Prototype *>::const_iterator i = cur.segment_data.prototypes.begin(); i != cur.segment_data.prototypes.end(); ++i) {
        const Prototype &prototype = **i;

		  const RawStroke *first_stroke = cur.ctx->impl->input.begin() + cur.start_stroke;

        unsigned model_strokes = num_strokes(prototype.strokes);
        const SymbolInfo *info = &cur.ctx->profman->GetInfoByPrototype(prototype);
        if (model_strokes <= cur.nstrokes) {
            NormalizedStrokeGroup &input = cur.segment_data.cached_input[model_strokes];
            if (input.empty()) {
                input = GetInputStrokes(first_stroke, model_strokes, cur.nstrokes);
            }
            
            const NormalizedStrokeGroup &model = prototype.strokes;
            RecognitionResult &match = cur.matches[prototype.id];
            e = MatchSymbol(model, input, match, matchers);
            if (FAILURE(e)) {
                return e;
            }
            
            match.bbox = bbox(first_stroke, first_stroke + model_strokes);
        }
    }
    
    return 0;
}


// Find ink offsets at which to start recognition based on the current start position and recognition
// results.  New start points are the current start point plus the number of strokes in each recognition
// result (discarding duplicates).
static int
find_new_start_points(RawStrokeGroup &strokes, 
                      unsigned curr_start, const std::vector<Match> &matches, std::vector<unsigned> &start_points)
{
    if (matches.empty()) {
        start_points.push_back(curr_start + 1);
        return 0;
    }
    
    unsigned nstrokes = num_strokes(strokes);
    
    for (std::vector<Match>::const_iterator i = matches.begin(); i != matches.end(); ++i) {
        const Match &match = *i;            

        unsigned new_start = curr_start + match.num_strokes;
        if (new_start < nstrokes) {
            std::vector<unsigned>::iterator j;
            for (j = start_points.begin(); j != start_points.end(); ++j) {
                if (new_start > *j) {
                    j = start_points.insert(j, new_start);
                    break;
                }
                
                if (new_start == *j) {
                    break;
                }
            }
            if (j == start_points.end()) {
                start_points.insert(j, new_start);
            }
        }
    }
    
    return 0;
}


static int
InsertMatch(const Match &match, std::vector<Match> &results, bool truncate = true)
{
    std::vector<Match>::iterator j;
    for (j = results.begin(); j != results.end(); ++j) {
        if (match.raw_score < j->raw_score) {
            j = results.insert(j, match);
            for (std::vector<Match>::iterator k = j + 1; k != results.end(); ++k) {
                if (k->symbol_info->unicode_char() == match.symbol_info->unicode_char()) {
                    results.erase(k);
                    break;
                }
            }
            break;
        }
        else if (match.symbol_info->unicode_char() == j->symbol_info->unicode_char()) {
            break;
        }                    
    }
    if (j == results.end()) { 
        results.push_back(match);
    }
    
    if (truncate) {
        if (results.size() > DefaultNumMatches) {
            results.erase(results.begin() + DefaultNumMatches, results.end());
        }
    }
    
    return 0;
}


// This function converts raw match scores to ratios below the baseline score.
// The baseline score is the worst (ie. highest) score in the vector.
// This function assumes that matches are sorted from best to worst.
static int
ComputeRelativeMatchScores(std::vector<Match> &matches)
{
    for (std::vector<Match>::iterator i = matches.begin(); i != matches.end();) {
        if (i->raw_score == std::numeric_limits<double>::infinity()) {
            i = matches.erase(i);
        }
        else {
            ++i;
        }
    }
    if (matches.size() == 1) {
        matches.front().confidence = 1.0;
    }
    else if (!matches.empty()) {
        double baseline_score = matches.back().raw_score;
        double total = 0.0;
        for (std::vector<Match>::iterator i = matches.begin(); i != matches.end(); ++i) {
            i->confidence = baseline_score / i->raw_score;
        }
    }

    return 0;
}


static int
recognize_multiple_segments(RecognitionContext &ctx)
{
    int e;
    

    std::vector<unsigned> start_points;
    start_points.push_back(0);

    unsigned old_start = 0;
    unsigned curr_start;    
 
 	const RawStrokeGroup &strokes = ctx.impl->input;
    unsigned segments = num_strokes(ctx.impl->input);

    while (!start_points.empty()) {
        // Determine at which stoke to begin recognition
        curr_start = start_points.back();
        start_points.pop_back();
        if (curr_start >= segments) {
            break;
        }
        
        ctx.impl->timing.curr_index = curr_start;
        
        VERBOSE2(*verb_out << std::endl << "beginning match at stroke " << curr_start << std::endl);

        // Recognition of a symbol consists of 5 phases:
        // 1. Calculuation of segmentation probabilities based on feature extraction.
        // (UNUSED & incorporated into step 6) 2. Pre-recognition phase for stacked and container symbols, as segmentation heuristics
        //    often miss these.
        // 3. Recognition phase for all other symbols, using normal segmentation routines.
        // (UNUSED) 4. Matching with any confusion symbols arising from phases 2 & 3.
        // (UNUSED) 5. Using prototype correlation data to reweight match results.
        // 6. Combination of individual match results into final results.  (Voting phase).
        
        //std::vector<Match> &matches = results[curr_start];
        std::map<PrototypeId, RecognitionResult> matches;
        
		  RecognitionCursor cur;
		  cur.ctx = &ctx;
		  cur.start_stroke = curr_start;
		  cur.nstrokes = segments - curr_start;

        // 1. Calculuation of segmentation probabilities based on feature extraction.  Upon return,
        // prototypes will be filled with symbols from the pruned database.
        e = ComputeSegmentationProbabilities(cur, DefaultSegmenters);
        if (FAILURE(e)) {
            return e;
        }
        
        VERBOSE2(
      	  *verb_out << "Segmentation results:" << std::endl;
       	 for (unsigned i = 0; i < cur.segment_data.segment_probabilities.size(); i++) {
         	   *verb_out << "  " << i + 1 << " strokes: " << cur.segment_data.segment_probabilities[i] << " ; " << cur.segment_data.stacked_segment_probabilities[i] << std::endl;
			}
        );
        
        
        // 3. Recognition phase for all other symbols, using normal segmentation routines.
        e = MatchSymbols(cur);
        if (FAILURE(e)) {
            return e;
        }
        
        VERBOSE2(
    	    *verb_out << "Recognition results:" << std::endl;
     	   for (std::map<PrototypeId, RecognitionResult>::const_iterator i = cur.matches.begin(); i != cur.matches.end(); ++i) {
     	       *verb_out << "  " << ctx.profman->GetInfoByPrototypeId(i->first).name() << " -> " << "(";
     	       for (unsigned j = 0; j < NumMatchers; j++) {
     	           *verb_out << i->second.matcher_scores[j];
     	           if (j != NumMatchers - 1) {
     	               *verb_out << ", ";
     	           }
     	       }
     	       *verb_out << ")" << std::endl;
     	   }
        );
        
        std::vector<Match> matcher_results[NumMatchers];
        for (std::map<PrototypeId, RecognitionResult>::const_iterator i = cur.matches.begin(); i != cur.matches.end(); ++i) {
            const RecognitionResult &match = i->second;
            Match out_match;
            out_match.bbox = match.bbox;
            out_match.symbol_id = i->first;
            out_match.num_strokes = num_strokes(ctx.profman->GetPrototypeById(i->first)->strokes);
            out_match.symbol_info = &ctx.profman->GetInfoByPrototypeId(i->first);

            double segmentation_bias;
            if ((out_match.symbol_info->get_symbol_type() & SymbolInfo::SYMBOLTYPE_CONTAINER)
             || (out_match.num_strokes > cur.segment_data.segment_probabilities.size())) {
                segmentation_bias = cur.segment_data.noheuristic_segment_probabilities[out_match.num_strokes - 1];
            }
            else if (out_match.symbol_info->get_symbol_type() & SymbolInfo::SYMBOLTYPE_STACKED) {
                segmentation_bias = cur.segment_data.stacked_segment_probabilities[out_match.num_strokes - 1];
            }
            else {
 					segmentation_bias = cur.segment_data.segment_probabilities[out_match.num_strokes - 1];
               //segmentation_bias = std::max(0.001, std::sqrt(segmentation_bias * segmentation_bias * segmentation_bias));
            }
             
            for (unsigned i = 0; i < NumMatchers; i++) {
                if (match.valid[i]) {
                    VERBOSE2(*verb_out << "matching vs. " << out_match.symbol_info->name() << "; seg. bias " << segmentation_bias << "->");
                    double seg_bias = segmentation_bias / (1.0 + (match.matcher_scores[i] / 50.0));
                    VERBOSE2(*verb_out << segmentation_bias << std::endl);
                    out_match.raw_score = match.matcher_scores[i] / seg_bias;
                    InsertMatch(out_match, matcher_results[i]);
                }
            }
        }
        
        for (unsigned i = 0; i < NumMatchers; i++) {
            ComputeRelativeMatchScores(matcher_results[i]);

            VERBOSE2(
            *verb_out << "Relative scores (" << i << "):" << std::endl;
            for (std::vector<Match>::const_iterator j = matcher_results[i].begin(); j != matcher_results[i].end(); ++j) {
                *verb_out << "  " << j->symbol_info->name() << " -> " << j->confidence << std::endl;
            }
            );
        }
        
        std::vector<Match> &match_results = ctx.segments[curr_start]->results;

        for (unsigned i = 0; i < NumMatchers; i++) {
            for (std::vector<Match>::iterator j = matcher_results[i].begin(); j != matcher_results[i].end(); ++j) {
                std::vector<Match>::iterator k;
                for (k = match_results.begin(); k != match_results.end(); ++k) {
                    if (k->symbol_info->unicode_char() == j->symbol_info->unicode_char()) {
                        Match adjusted_match = *k;
                        match_results.erase(k);
                        adjusted_match.raw_score += j->confidence * matcher_weights[i];
                        InsertMatch(adjusted_match, match_results, false);
                        k = match_results.begin(); // just prevent k == match_results.end() below; we know there's at least one entry
                        break;
                    }
                }
                
                if (k == match_results.end()) {
                    Match new_match = *j;
                    new_match.raw_score = j->confidence;
                    InsertMatch(new_match, match_results, false);
                }
            }
        }
        
        if (match_results.size() > DefaultNumMatches) {
            match_results.erase(match_results.begin(), match_results.end() - DefaultNumMatches - 1);
        }

        // Take the reciprocal of the scores so that lower is better.  
        for (std::vector<Match>::iterator i = match_results.begin(); i != match_results.end(); ++i) {
            i->raw_score = 1.0 / i->raw_score;
        }
        
        // The calls to InsertMatch() above will have placed lowest-scored (ie. worst) results
        // first.  Now that we've taken the reciprocal of all the scores and lower scores are better,
        // we need to reverse the list of results so it is sorted from best (lower) to worst (higher).
        std::reverse(match_results.begin(), match_results.end());
        ComputeRelativeMatchScores(match_results);

        double total = 0.0;
        for (std::vector<Match>::iterator i = match_results.begin(); i != match_results.end(); ++i) {
            total += i->confidence;
        }
        for (std::vector<Match>::iterator i = match_results.begin(); i != match_results.end(); ++i) {
            i->confidence = i->confidence / total;
        }

        VERBOSE2(
        *verb_out << "Final results:" << std::endl;
        if (!match_results.empty()) {
            const Rect<long> &box = match_results.front().bbox;
            *verb_out << " (" << box.left << "," << box.top << ")-(" << box.right << "," << box.bottom << ")" << std::endl;            
        }
        for (std::vector<Match>::const_iterator i = match_results.begin(); i != match_results.end(); ++i) {
            *verb_out << "  " << i->symbol_info->name() << " -> " << i->confidence << std::endl;
        }
        );
        
        find_new_start_points(ctx.impl->input, curr_start, match_results, start_points);

        old_start = curr_start;
    }

    return 0;
}


// This is the entry-point to the recognizer
int
recognize(RecognitionContext &ctx)
{
	 if (!ctx.impl) {
	 	ctx.impl = new RecognitionContextImpl;
	}

	ctx.clear_segments();

    segment_set_profile(ctx.profman);

    RawStrokeGroup joined_input;
    std::vector<unsigned> isolated_dots;

    RawStrokeGroup temp_input = join_strokes(ctx, isolated_dots);
    unsigned ndots = static_cast<unsigned>(isolated_dots.size());
    Stroke<long> *strokes_nodots = DEBUG_NEW Stroke<long>[num_strokes(temp_input) - ndots];
    unsigned j = 0;
    std::vector<unsigned>::iterator dot = isolated_dots.begin();
    for (unsigned i = 0; i < num_strokes(temp_input); i++) {
        if (dot != isolated_dots.end() && *dot == i) {
            VERBOSE2(*verb_out << "DOT AT " << *dot << std::endl);
            ++dot;
        }
        else {
            strokes_nodots[j++] = copy(temp_input[i]);
        }
    }
    ctx.impl->input.set_strokes(strokes_nodots, num_strokes(temp_input) - ndots);

    VERBOSE2(
	 	*verb_out << "NO DOT STROKES:" << std::endl;
    	unsigned i;
  	  	for (i = 0; i < num_strokes(joined_input); i++) {
  	     Rect<long> bb = bbox(joined_input[i]);
        *verb_out << "(" << bb.left << "," << bb.top << ")->(" << bb.right << "," << bb.bottom << ")" << std::endl;
    	}    
	);

    // extract timing information from the strokes (namely elapsed time between them), if available
	 SegmentTimingData &timing_data = ctx.impl->timing;
    timing_data.elapsed.clear();
    
	 unsigned num_input_strokes = num_strokes(ctx.strokes);
	 unsigned num_joined_strokes = num_strokes(ctx.impl->input);

	 RawStrokeGroup::const_iterator i;
    for (i = ctx.strokes.begin(); i != ctx.strokes.end(); ++i) {
        if (!i->time) {
            break;
        }
    }
    if (i == ctx.strokes.end()) {
		timing_data.elapsed.insert(timing_data.elapsed.end(), num_input_strokes - 1, 0);
		unsigned j = 0;
      for (RawStrokeGroup::const_iterator this_stroke = ctx.strokes.begin(); this_stroke != ctx.strokes.end() - 1; ++this_stroke) {
            RawStrokeGroup::const_iterator next_stroke = this_stroke + 1;
            timing_data.elapsed[j] = next_stroke->time[0] - this_stroke->time[num_points(*this_stroke) - 1];
            VERBOSE2(*verb_out << "elapsed time " << j << "->" << j + 1 << " is " << timing_data.elapsed[j] << "ms" << std::endl);
				++j;
        }
    }
        
    int e = 0;
    e = recognize_multiple_segments(ctx);
    if (FAILURE(e)) {
        return e;
    }
    
    
    if (!timing_data.elapsed.empty()) {
        // update the input order, as it may have been shifted during stroke joining
		  VERBOSE2(
        		for (size_t i = 0; i < num_joined_strokes; i++) {
            	*verb_out << "input_order[" << i << "] = " << ctx.segments[i]->input_order.front() << std::endl;
        		}
	 		);
            
        bool got_lift = false;
        for (std::vector<RecognitionSegment *>::iterator i = ctx.segments.begin(); i != ctx.segments.end(); ++i) {
		  		std::vector<Match> &results = (*i)->results;
            if (!results.empty()) {
                const Match &best_match = results.front();
                if (best_match.num_strokes > 1) {
                    VERBOSE2(*verb_out << "multi-stroke at segment " << i - ctx.segments.begin() << " consists of ");
                    std::vector<unsigned> input_strokes;
                    input_strokes.reserve(best_match.num_strokes);
						  for (std::vector<RecognitionSegment *>::const_iterator j = i; j != i + best_match.num_strokes; ++j) {
							const std::vector<unsigned> &input_order = (*j)->input_order;
						  	for (std::vector<unsigned>::const_iterator k = input_order.begin(); k != input_order.end(); ++k) {
                        input_strokes.push_back(*k);
                        VERBOSE2(*verb_out << *k << ' ');
							}
                    }
                    VERBOSE2(*verb_out << std::endl);
                    
                    std::vector<unsigned>::const_iterator s;
                    for (s = input_strokes.begin(); s != input_strokes.end(); ++s) {
                        if (std::find(input_strokes.begin(), input_strokes.end(), *s + 1) == input_strokes.end()
								 && std::find(input_strokes.begin(), input_strokes.end(), *s - 1) == input_strokes.end()) {
                            VERBOSE2(*verb_out << "stroke " << *s << " missing neighbouring stroke!" << std::endl);
                            break;
                        }
                    }
                    
                    if (s == input_strokes.end()) {
                        std::sort(input_strokes.begin(), input_strokes.end());
                        
                        unsigned symbol_elapsed = 0;
                        for (std::vector<unsigned>::const_iterator k = input_strokes.begin(); k < input_strokes.end() - 1; ++k) {
                            time_between_strokes.AddSample(timing_data.elapsed[*k]);
                            VERBOSE2(*verb_out << "adding sample " << *k << "; average is now " << time_between_strokes.sample_mean() << std::endl);
                        }
                        
                        got_lift = true;
                    }
                }
            }
        }
	 
	 	timing_data.elapsed.clear();
    }
    

    if (!isolated_dots.empty()) {
        static SymbolInfo isolated_dot_info;
        isolated_dot_info.set_name("dot");
        isolated_dot_info.set_unicode_char('.');
        
        Match isolated_dot_match;
        isolated_dot_match.num_strokes = 1;
        isolated_dot_match.symbol_info = &isolated_dot_info;
        isolated_dot_match.confidence = 1.0;
        
        unsigned ndots = 0;
        for (std::vector<unsigned>::const_iterator i = isolated_dots.begin(); i != isolated_dots.end(); ++i) {
			std::vector<RecognitionSegment *>::iterator segment = ctx.segments.insert(ctx.segments.begin() + *i + ndots, new RecognitionSegment);
			isolated_dot_match.bbox = bbox(ctx.strokes[*i]);
			(*segment)->results.push_back(isolated_dot_match);
        }
    }
    
    return e;
}



static void
insert_dirty_index_for_removal(std::vector<unsigned> &dirty, unsigned index)
{
	std::vector<unsigned>::iterator i;
	for (i = dirty.begin(); i != dirty.end(); ++i) {
		if (*i > index) {
			break;
		}
		else if (*i == index) {
			return;
		}
	}

	i = dirty.insert(i, index);
	++i;

	if (i != dirty.end()) {
		for (; i != dirty.end(); ++i) {
			--(*i);
		}
	}
}


static int
process_stroke_removal(RecognitionContext &ctx, std::vector<RawStroke *> &strokes, unsigned index)
{
	std::vector<RecognitionSegment *>::iterator rm_segment = ctx.segments.end();
	std::vector<unsigned>::iterator input_stroke;

	for (std::vector<RecognitionSegment *>::iterator j = ctx.segments.begin(); j != ctx.segments.end(); ++j) {
		RecognitionSegment &segment = **j;
		for (std::vector<unsigned>::iterator k = segment.input_order.begin(); k != segment.input_order.end(); ++k) {
			if (*k == index) {
				std::cout << "found stroke " << index << " at segment " << j - ctx.segments.begin() << std::endl;
				rm_segment = j;
				segment.dirty = RecognitionSegment::DIRTY_GROUPING;
				input_stroke = k;
			}
			else if (*k > index) {
				--(*k);
			}
		}
	}

	if (rm_segment != ctx.segments.end()) {
		RecognitionSegment &segment = **rm_segment;
		// remove the removed stroke from the input strokes used by this segment
		segment.input_order.erase(input_stroke);

		if (segment.input_order.empty()) {
			ctx.segments.erase(rm_segment);
			std::cout << "removing segment " << rm_segment - ctx.segments.begin() << std::endl;
		}
		else {
			unsigned prev_segments = std::min(MaxGroupSize-1, static_cast<unsigned>(rm_segment - ctx.segments.begin()));
			for (std::vector<RecognitionSegment *>::iterator k = rm_segment - prev_segments; k != rm_segment + 1; ++k) {
				RecognitionSegment &dirty_segment = **k;
				std::vector<Match> &dirty_results = dirty_segment.results;
				for (unsigned m = 0; m < dirty_results.size(); ) {
					const Match &match = dirty_results[m];
					if (match.num_strokes >= rm_segment - k) {
						dirty_results.erase(dirty_results.begin() + m);
					}
					else {
						++m;
					}
				}
			}
		}
	}
}



int
recognize_incremental(RecognitionContext &ctx, add_strokes_buf *additions, remove_strokes_buf *removals)
{
	int e;

	unsigned nstrokes = num_strokes(ctx.strokes);

	std::vector<RawStroke *> strokes;
	for (RawStrokeGroup::iterator i = ctx.strokes.begin(); i != ctx.strokes.end(); ++i) {
		strokes.push_back(&*i);
	}

	for (remove_strokes_buf *buf = removals; buf; buf = buf->next) {
		std::sort(buf->indices.begin(),buf->indices.end(), std::greater<unsigned>());
		for (std::vector<unsigned>::const_iterator i = buf->indices.begin(); i != buf->indices.end(); ++i) {
			strokes[*i] = 0;
			process_stroke_removal(ctx, strokes, *i);
		}
	}

	unsigned orig_nsegments = ctx.segments.size();
	for (add_strokes_buf *buf = additions; buf; buf = buf->next) {
		for (RawStroke *s = buf->strokes; s != buf->strokes + buf->nstrokes; ++s) {
			RecognitionSegment *seg = new RecognitionSegment;
			seg->stroke = copy(*s);
			seg->input_order.push_back(nstrokes++);
			seg->dirty = RecognitionSegment::DIRTY_SEGMENTATION;
			ctx.segments.push_back(seg);
		}
	}

	/* TODO:
	   - split segmentation (strokes-->segments) and group estimation into separate modules
		- change segmentation (join_strokes et al) to be two pass:
		  + pass 1 : sort strokes by nearness to each other
		  + pass 2 : join strokes into segments by the following algorithm:
		  
		      curr_stroke = strokes[0]
				build segment, mark dirty grouping
				do {
					change = false
					while (next segment joins with curr segment) {
						merge segments
					}
					while (prev segmemnt joins with curr segment) {
						merge segments
						change = true
					}
				} while (changes are still occurring)

		- to process stroke additions, create segments for each stroke and mark dirty segmentation, grouping
	*/

	e = update_segmentation(ctx.segments);

	return e;
}



int
recognizer_initialize()
{
    // ensure the recognizer is only initialized once
    if (matcher_functions[ElasticMatcher]) {
        return 0;
    }
    
    matcher_functions[ElasticMatcher]          = &ElasticMatch;
    matcher_functions[DeformationMatcher]      = &deformation_match;
    matcher_functions[DirectionElementMatcher] = &direction_element_match;
    matcher_functions[StrokeFeatureMatcher]    = &FeatureBasedMatch;
    matcher_functions[ChainCodeMatcher]        = &chaincode_match;
        
    return 0;
}


void
recognizer_shutdown()
{
}


}

