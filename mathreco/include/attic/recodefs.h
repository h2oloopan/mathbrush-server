#ifndef RECODEFS_H_
#define RECODEFS_H_


namespace scg
{


enum MatcherIndices
{
    ElasticMatcher,
    DeformationMatcher,
    DirectionElementMatcher,
    StrokeFeatureMatcher,
    ChainCodeMatcher,
    NumMatchers
};

enum MatcherBits
{
    ElasticMatcherBit          = (1 << ElasticMatcher),
    DeformationMatcherBit      = (1 << DeformationMatcher),
    DirectionElementMatcherBit = (1 << DirectionElementMatcher),
    StrokeFeatureMatcherBit    = (1 << StrokeFeatureMatcher),
    ChainCodeMatcherBit        = (1 << ChainCodeMatcher)
};


}


#endif

