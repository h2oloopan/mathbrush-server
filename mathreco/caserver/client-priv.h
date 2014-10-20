#ifndef CLIENT_PRIV_H_
#define CLIENT_PRIV_H_

#include "MathRecognizer.h"
#include "symbols.h"

scg::ExpressionTree *exprcmd(int cmd); 
int casclient_pushexpr(const scg::ExpressionTree *tree);

#endif
