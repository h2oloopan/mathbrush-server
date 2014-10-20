#ifndef STRUCTURAL_H_
#define STRUCTURAL_H_


#include "group.h"


namespace scg
{


struct Match;


int chaincode_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, Match &match);


}


#endif

