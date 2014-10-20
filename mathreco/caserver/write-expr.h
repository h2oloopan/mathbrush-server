#ifndef WRITE_EXPR_H_
#define WRITE_EXPR_H_

#include "MathRecoTypes.h"

int writeexpr(const scg::ExpressionTree *tree, char **buf, size_t len, size_t *wrlen);

#endif
