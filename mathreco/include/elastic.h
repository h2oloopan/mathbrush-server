#ifndef ELASTIC_H_
#define ELASTIC_H_

#include "group.h"
#include <vector>

namespace scg {

void elastic_shutdown();
double elasticdist(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const std::vector<size_t> &inputorder);

}

#endif
