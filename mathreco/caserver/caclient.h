#ifndef CACLIENT_H_
#define CACLIENT_H_

#include "MathRecognizer.h"
#include "CASOperation.h"
#include <set>

#define CACLIENT_NETERROR -1
#define CACLIENT_LIBERROR -2

int casclient_geterror(); // gets the error code set when an operation is invoked
const char *casclient_replystring(); // gets a string describing the outcome of an attempted operation

int casclient_connect();
int casclient_disconnect();
bool casclient_isconnected();

std::set<CASOperation> casclient_getoperations(const scg::ExpressionTree *expr);
scg::ExpressionTree *casclient_dooperation(const CASOperation &op, scg::ExpressionTree *expr);

scg::ExpressionTree *casclient_evaluate(scg::ExpressionTree *expr);
scg::ExpressionTree *casclient_evalnum(scg::ExpressionTree *expr);
scg::ExpressionTree *casclient_simplify(scg::ExpressionTree *expr);

scg::ExpressionTree *casclient_expand(scg::ExpressionTree *expr);
scg::ExpressionTree *casclient_factor(scg::ExpressionTree *expr);

scg::ExpressionTree *casclient_solve(scg::ExpressionTree *expr, const scg::ExpressionTree *var);
scg::ExpressionTree *casclient_solvesystem(scg::ExpressionTree *expr);

scg::ExpressionTree *casclient_determinant(scg::ExpressionTree *expr);
scg::ExpressionTree *casclient_inverse(scg::ExpressionTree *expr);
scg::ExpressionTree *casclient_rank(scg::ExpressionTree *expr);
scg::ExpressionTree *casclient_nullspace(scg::ExpressionTree *expr);

#endif
