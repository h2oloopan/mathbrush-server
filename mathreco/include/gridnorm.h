#ifndef _GRIDNORM_H_
#define _GRIDNORM_H_


#include "group.h"


namespace scg
{


int gridnorm_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const std::vector<size_t> &inputorder, double &score);


}


#endif
