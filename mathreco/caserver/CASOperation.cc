#include "CASOperation.h"
#include "MathRecoTypes.h"
#include "grammar-values.h"
#include "ode.h"
#include "expr.h"
#include "cmdcode.h"
#include "error.h" // for error codes used in serialization
#include <set>
#include <cstring>
#include <vector>

CASOperation::CASOperation() :  commandCode(CASCMD_NONE), variable("") {}
CASOperation::CASOperation(int code) : commandCode(code), variable("") {}
CASOperation::CASOperation(int code, std::string var) : commandCode(code), variable(var) {}

int CASOperation::serialize(char *buf, size_t *n) const
{
    if (buf && !n) return E_INVALID;
    
    // buffer to store the serialized contents
    std::vector<char> vbuf;
    
    // command code
	vbuf.push_back((char)commandCode);
    // length of variable string
    //int l = variable ? std::strlen(variable) : 0;
    int l = variable.length()+1; // include the terminal character
	vbuf.push_back((char)l);
	
    // actual variable string
    for (int i = 0; i < l; ++i)
        vbuf.push_back(variable.c_str()[i]);
    
    // write to buffer if it is ready to receive the data
    int e = E_OK;
    if (buf)
    {
        if (vbuf.size() <= *n)
            memcpy(buf, &vbuf[0], vbuf.size());
        else
            e = E_NOTREADY;
    }
    if (n)
        *n = vbuf.size();
    
    return e;
}

// initialize from buffer
CASOperation::CASOperation(const char *buf, size_t n)
{
    if (!buf || n < 1)
    {
        commandCode = CASCMD_NONE; // can't deserialize
		return;
	}
    
    // command code
	commandCode = (int)*buf;
	++buf; --n;
    
	// length of variable string
    const char l = *buf;
	++buf; --n;
    
	if (n < (unsigned)l)
    {
        commandCode = CASCMD_NONE; // can't deserialize
		return;
	}
    
	// variable string
    variable = "";
    if (l > 0)
    {
        char var[l];
        for (int i=0; i < l; ++i)
        {
            var[i] = *buf;
            buf++; n--;
        }
        variable = std::string(var);
    }
}

bool CASOperation::operator<(const CASOperation& other) const
{
    if (commandCode == other.commandCode)
        return variable < other.variable;
    
    return commandCode < other.commandCode;
}

std::string CASOperation::commandName() const
{
    switch (commandCode)
    {
        case CASCMD_EXPR:       return "Expression";
        case CASCMD_EVAL:       return "Evaluate";
        case CASCMD_EVALNUM:    return "Evaluate numerically";
        case CASCMD_SIMPLIFY:   return "Simplify";
        case CASCMD_EXPAND:     return "Expand";
        case CASCMD_SOLVE:      return "Solve";
        case CASCMD_SOLVEODE:   return "Solve ODE";
        case CASCMD_SOLVEIVP:   return "Solve IVP";
        case CASCMD_SOLVESYSTEM:return "Solve system";
        case CASCMD_FACTOR:     return "Factor";
        case CASCMD_DETERMINANT:return "Determinant";
        case CASCMD_INVERSE:    return "Inverse";
        case CASCMD_RANK:       return "Rank";
        case CASCMD_NULLSPACE:  return "Nullspace";
        default:                return "(unknown command)";
    }
}

std::string CASOperation::text() const
{
    return commandName() + (variable.empty() ? "" : (" for " + variable));
}

// Look for a node in the given expression tree corresponding to the variable given by name.
// Return a pointer to this node (subtree) if found, NULL otherwise.
const scg::ExpressionTree *getVariableExpr(const scg::ExpressionTree *expr, std::string variable)
{
    if (variable.empty()) return NULL;
    
    // if this is a variable expression, check if it's the one we're looking for
    if (expr->type() == scg::VAR_EXPR)
        if (strcmp(expr->child(scg::VAR_NAME)->str(), variable.c_str()) == 0)
            return expr;
    
    const scg::ExpressionTree *varExpr = NULL;

    // traverse the tree looking for expressions of type VAR_EXPR
    for (size_t i = 0; i < expr->nchildren(); ++i)
    {
        varExpr = getVariableExpr(expr->child(i), variable);
        if (varExpr) break;
    }
    
    return varExpr;
}

// Gets a list of the variables involved in the given expression.
void getVariables(const scg::ExpressionTree *expr, std::set<std::string> &vars)
{
    // if this is a variable expression, we're done
    if (expr->type() == scg::VAR_EXPR)
    {
        vars.insert(std::string(expr->child(scg::VAR_NAME)->str()));
        return;
    }

    // traverse the tree looking for expressions of type VAR_EXPR
    for (size_t i = 0; i < expr->nchildren(); ++i)
    {
        getVariables(expr->child(i), vars);
    }
}

// Gets a list of operations that are relevant for the given expression.
void getOperations(const scg::ExpressionTree *expr, std::set<CASOperation> &ops)
{
    scg::SemanticId id = expr->type();
    
    exprtree *t = mkexprtree(expr);
    exprgrp *grp = dynamic_cast<exprgrp*>(t);
    bool flag = false;
    if (grp)
    {
        // 1. Check ODE
        DE *de = new DE;
        if (OdeParser::isDifferentialEquation(grp, de))
        {
            ops.insert(CASOperation(CASCMD_SOLVEODE, de->dependentVariable
                                    + "(" + de->independentVariable + ")"));
            flag = true;
        }
        delete de;
        
        // 2. Check IVP
        IVP *ivp = new IVP;
        if (OdeParser::isInitialValueProblem(grp, ivp))
        {
            ops.insert(CASOperation(CASCMD_SOLVEIVP, ivp->de->dependentVariable
                                    + "(" + ivp->de->independentVariable + ")"));
            flag = true;
        }
    }
    delete t;
    
    if (flag) return; // ODEs exclude other operations
    
    // TODO: find out if expr has variables. if so, certain operations may apply
    // (e.g. differentiate (needs to be added to server), solve for variable)
    std::set<std::string> vars;
    getVariables(expr, vars);
    
    switch (id)
    {
        case scg::NUM_EXPR:
            // special case for numbers pi and e
            if (!strcmp(expr->child(scg::NUM_VALUE)->str(), "pi") ||
                !strcmp(expr->child(scg::NUM_VALUE)->str(),  "e"))
                ops.insert(CASOperation(CASCMD_EVALNUM));
            else
                ops.insert(CASOperation(CASCMD_FACTOR));
            break;
            
        case scg::NEG_EXPR:
            getOperations(expr->child(scg::UNARY_TERM), ops);
            break;
            
        case scg::REL_EXPR:
            ops.insert(CASOperation(CASCMD_EXPAND));
            for (std::set<std::string>::iterator i = vars.begin(); i != vars.end(); ++i)
                ops.insert(CASOperation(CASCMD_SOLVE, *i));
            break;
            
        case scg::ADD_EXPR:
        case scg::MULT_EXPR:
        case scg::SUP_EXPR:
        case scg::ROOT_EXPR:
        case scg::FN_EXPR:
        case scg::INTEGRAL_EXPR:
        case scg::SUM_EXPR:
        case scg::FACTORIAL_EXPR:
            ops.insert(CASOperation(CASCMD_EXPAND));
            ops.insert(CASOperation(CASCMD_EVAL));
            ops.insert(CASOperation(CASCMD_EVALNUM));
            switch (id)
            {
                case scg::ADD_EXPR:
                case scg::MULT_EXPR:
                case scg::SUP_EXPR:
                case scg::ROOT_EXPR:
                    ops.insert(CASOperation(CASCMD_SIMPLIFY));
                    ops.insert(CASOperation(CASCMD_FACTOR));
                    for (std::set<std::string>::iterator i = vars.begin(); i != vars.end(); ++i)
                        ops.insert(CASOperation(CASCMD_SOLVE, *i));
            }
            break;
            
        case scg::FRAC_EXPR:
            getOperations(expr->child(scg::FRAC_NUMER), ops);
            getOperations(expr->child(scg::FRAC_DENOM), ops);
            break;
            
        case scg::PAREN_EXPR:
            // handle regular parentheses at least
            if (!strcmp(expr->child(scg::PAREN_OPEN )->str(), "(")
             && !strcmp(expr->child(scg::PAREN_CLOSE)->str(), ")"))
                getOperations(expr->child(scg::PAREN_CONTENTS), ops);
            break;
            
        case scg::MATRIX_EXPR:
            ops.insert(CASOperation(CASCMD_DETERMINANT));
            ops.insert(CASOperation(CASCMD_INVERSE));
            ops.insert(CASOperation(CASCMD_RANK));
            ops.insert(CASOperation(CASCMD_NULLSPACE));
            break;
            
        // TODO: can both of these types be solved end-to-end?
        // enable Multi in the grammar
        case scg::LIST_EXPR:
        case scg::MULTI_EXPR:
            ops.insert(CASOperation(CASCMD_SOLVESYSTEM));
            break;
    }
    
}
