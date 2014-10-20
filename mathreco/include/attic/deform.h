#ifndef DEFORM_H_
#define DEFORM_H_


#include "group.h"


namespace scg
{


struct Match;


int deformation_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, Match &match);


}


#endif

