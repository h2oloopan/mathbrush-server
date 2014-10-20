#include "interp.h"
#include "parms.h"

namespace scg {

double CuspLowerBound = RegisterParameterDouble("CuspLowerBound", &CuspLowerBound);
double CuspUpperBound = RegisterParameterDouble("CuspUpperBound", &CuspUpperBound);
double InterpArclengthRatio = RegisterParameterDouble("InterpArclengthRatio", &InterpArclengthRatio);

}