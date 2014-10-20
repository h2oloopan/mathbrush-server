#ifndef CASOPERATION_H_
#define CASOPERATION_H_

#include "cmdcode.h"
#include "MathRecoTypes.h"
#include <string>
#include <set>

// Encapsulates a single operation to be applied to an expression.
struct CASOperation
{
    int commandCode;
    std::string variable;
    std::string commandName() const;
    std::string text() const;
    CASOperation(); 
    CASOperation(int code);
    CASOperation(int code, std::string var);
    
    // determines the relevant operations for a given expression.
    bool operator<(const CASOperation &other) const; // used by std::set
    
    // serialization
    int serialize(char *buf, size_t *n) const;
    CASOperation(const char *buf, size_t n);
};

// Determines the relevant operations for a given expression.
void getOperations(const scg::ExpressionTree *expr, std::set<CASOperation> &ops);

// Determines the variables in an expression.
void getVariables(const scg::ExpressionTree *expr, std::set<std::string> &vars);

const scg::ExpressionTree *getVariableExpr(const scg::ExpressionTree *expr, std::string variable);

#endif
