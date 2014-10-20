#ifndef LCONV_H_
#define LCONV_H_

#include "MathRecoTypes.h"

const scg::ExpressionTree *latex2tree(const char *latex, int &stat);
void reset_relations();
void clear_relations();

#endif
