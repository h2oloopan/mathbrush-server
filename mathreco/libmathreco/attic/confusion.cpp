#include "confusion.h"

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "error.h"
#include "interp.h"
#include "match.h"
#include "parms.h"
#include "profile.h"
#include "profman.h"
#include "recog.h"


namespace scg
{

const static double ElasticConfusionThreshold = GetParameterDouble("ElasticConfusionThreshold");
const static double DirectionConfusionThreshold = GetParameterDouble("DirectionConfusionThreshold");
const static double ChainCodeConfusionThreshold = GetParameterDouble("ChainCodeConfusionThreshold");
const static double ConfusionRatioConfidenceThreshold = GetParameterDouble("ConfusionRatioConfidenceThreshold");


int
ComputeProfileConfusion(const ProfileManager &profman, const Profile &from, const Profile &to, std::vector<InterProfileCorrelation> *confusion)
{
    //std::map<PrototypeId, double> self_match_score;
    
    for (Profile::const_iterator from_symbol = from.begin(); from_symbol != from.end(); ++from_symbol) {
        DEBUG_ONLY(debug_out << "matching from " << profman.GetInfoBySymbol(*from_symbol).name() << std::endl);
        for (ProfileSymbol::const_iterator from_proto = from_symbol->begin(); from_proto != from_symbol->end(); ++from_proto) {
            RecognitionResult match;
            //MatchSymbol(from_proto->strokes, from_proto->strokes, match, ElasticMatcherBit);
            
            //self_match_score[from_proto->id] = match.matcher_scores[ElasticMatcher];

            for (Profile::const_iterator to_symbol = to.begin(); to_symbol != to.end(); ++to_symbol) {
                DEBUG_ONLY(debug_out << "  to " << profman.GetInfoBySymbol(*to_symbol).name() << std::endl);
                for (ProfileSymbol::const_iterator to_proto = to_symbol->begin(); to_proto != to_symbol->end(); ++to_proto) {
                    if (to_proto->id != from_proto->id) { // don't get symbols confused with themselves
                        if (num_strokes(to_proto->strokes) == num_strokes(from_proto->strokes)) {
                            match.invalidate();
                            MatchSymbol(from_proto->strokes, to_proto->strokes, match, ElasticMatcherBit | DirectionElementMatcherBit | ChainCodeMatcherBit);
                            
                            double confidence = 1.0 - (match.matcher_scores[ElasticMatcher] / ElasticConfusionThreshold);
                            if (match.valid[ElasticMatcher] && confidence > ConfusionRatioConfidenceThreshold) {
                                confusion[ElasticMatcher].push_back(InterProfileCorrelation(from_proto->id, to_proto->id, confidence));
                            }
                            confidence = 1.0 - (match.matcher_scores[DirectionElementMatcher] / DirectionConfusionThreshold);
                            if (match.valid[DirectionElementMatcher] && confidence > ConfusionRatioConfidenceThreshold) {
                                confusion[DirectionElementMatcher].push_back(InterProfileCorrelation(from_proto->id, to_proto->id, confidence));
                            }
                            confidence = 1.0 - (match.matcher_scores[ChainCodeMatcher] / ChainCodeConfusionThreshold);
                            if (match.valid[ChainCodeMatcher] && confidence > ConfusionRatioConfidenceThreshold) {
                                confusion[ChainCodeMatcher].push_back(InterProfileCorrelation(from_proto->id, to_proto->id, confidence));
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}


}

