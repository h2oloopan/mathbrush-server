#ifndef DIRELEM_H_
#define DIRELEM_H_


#include "group.h"


namespace scg
{


struct Match;


int direction_element_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, Match &match);


}


#endif

