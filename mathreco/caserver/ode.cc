#include "expr.h"
#include "ode.h"
#include "symbols.h"
#include "log.h"
#include <string>
#include <vector>
//#include <stdlib.h>
#include <stdexcept>
#include <cstdio>

// this is used as a placeholder for the independent variable name in sage code
// in some cases it may not be determined until the end
#define INDEP_VAR_PLACEHOLDER "independentVariable"

std::vector<std::string> indepVarCandidates; 
//bool locked = false; // use this to protect info from being changed

// short OdeParser::order = -1; // set to 0 or -1 if no derivative found
// std::string OdeParser::dependentVariable = "";
// std::string OdeParser::independentVariable = "";
// std::string OdeParser::sage = "";

// void OdeParser::reset()
// {
// 	//locked = false;
// // 	order = -1; // set to 0 or -1 if no derivative found
// // 	dependentVariable = "";
// // 	independentVariable = "";
// // 	sage = "";
// 	indepVarCandidates.clear();
// }

// determines if the expression represents the dependent variable.
// If we have determined the dependent variable already, it must match.
// Otherwise, if the expression could be the dependent variable, we take it as such.
bool OdeParser::parseDependentVariable(exprtree* tree, DE* de) {

	exprleaf* xterm;
	
	if (tree->sid == scg::VAR_EXPR)
		xterm = (exprleaf*)((exprgrp*)tree)->children[0];
	else if (tree->sid == scg::TERMINAL_EXPR) // terminal will do too
		xterm = (exprleaf*)tree;
	else
		return false;	
	
	if (xterm->terms.size()!=1) // assume one-character variables
		return false;
	std::string xstr = scg::symdb_findsymbol_unicode(xterm->terms[0])->name;

	if (de->dependentVariable == "") { // && !locked) {
		de->dependentVariable = xstr;
	}
	else { // we had already encountered the dependent variable
		if (de->dependentVariable != xstr)
			return false; // inconsistent (multivariate not supported)
	}
	
	return true;
	
}

// Get the variable name as a string. 
// Return true if you can; false otherwise.
bool getVariableName(exprtree *tree, std::string &name) {
	
	exprleaf* term;
	
	if (tree->sid == scg::VAR_EXPR) 
		term = (exprleaf*)((exprgrp*)tree)->children[0];
	else if (tree->sid == scg::TERMINAL_EXPR)
		term = (exprleaf*)tree;
	else
		return false;
	
	if(term->terms.size() != 1) // assume one-character variable	
		return false;
	name = scg::symdb_findsymbol_unicode(term->terms[0])->name;
	return true;
}

bool OdeParser::parseIndependentVariable(exprtree* tree, DE* de) {
	
	// assume tree is the independent variable
	std::string tstr;
	if (!getVariableName(tree,tstr))
		return false;

	if (de->independentVariable == "") { // && !locked) {
		de->independentVariable = tstr;
	}
	else { // we already encountered the independent variable
		if (de->independentVariable != tstr)
			return false; // inconsistent 
	}
	
	return true;
}
	
bool OdeParser::parseLeibnizNumerator(exprgrp* dnx, int &nval, DE* de)
{
	// looking for expression of the form d^nx
	if (!(dnx->sid == scg::MULT_EXPR))
		return false;
	
	exprgrp* dn = (exprgrp*)dnx->children[0];
	
	//int nval; // order of the derivative
	
	if (dn->sid == scg::VAR_EXPR) { // first derivative (n==1)
		
		// this is the 'd' (not really a variable)
		exprleaf* dterm = (exprleaf*)((exprgrp*)dn)->children[0];
		if (!(dterm->terms.size()==1 && dterm->terms[0] == 'd'))
			return false;

		nval = 1; // implied
	}
	else if (dn->sid == scg::SUP_EXPR) { // higher order derivative (n>=2)
	
		exprgrp* d = (exprgrp*)dn->children[0];
		if (d->sid != scg::VAR_EXPR) // this is the 'd' (not really a variable)
			return false;
		exprleaf* dterm = (exprleaf*)d->children[0];
		if (!(dterm->terms.size()==1 && dterm->terms[0] == 'd'))
			return false;
		
		exprgrp* n = (exprgrp*)dn->children[1]; // this is the order of the derivative
		if (n->sid != scg::NUM_EXPR) 
			return false;
		exprleaf* nterm = (exprleaf*)n->children[0];
		if (nterm->terms.size() != 1) // todo: 2<=n<=9 is an (unnecessary?) limitation
			return false;
		std::string nstr = scg::symdb_findsymbol_unicode(nterm->terms[0])->name;
		
		nval = atoi(nstr.c_str()); // get integer from unicode value
	}
	else
		return false;
	
	if (!(nval >= 1 && nval <= 9)) // n must be in an appropriate range
		return false;
	
	//exprgrp* x = (exprgrp*)dnx->children[1];
	return parseDependentVariable(dnx->children[1], de);
		
// 	if (!(x->sid == scg::VAR_EXPR)) // this is the dependent variable 
// 		return false;
// 	exprleaf* xterm = (exprleaf*)x->children[0];
// 	if (!(xterm->terms.size()==1)) // assume one-character variable
// 		return false;
// 	std::string xstr = scg::symdb_findsymbol_unicode(xterm->terms[0])->name;
// 
// 	if (dependentVariable == "") {
// 		dependentVariable = xstr;
// 	}
// 	else { // we already encountered the dependent variable
// 		if (!(dependentVariable == xstr))
// 			return false; // inconsistent (or multivariate)
// 	}
	
	return true;
}

// nval should have been determined in processing numerator
bool OdeParser::parseLeibnizDenominator(exprgrp* dtn, int nval, DE* de) {
	
	if (dtn->sid != scg::MULT_EXPR)
		return false;

	exprgrp* d = (exprgrp*)dtn->children[0]; // this is the 'd' (not really a variable)
	
	if (d->sid != scg::VAR_EXPR) 
		return false;
	exprleaf* dterm = (exprleaf*)d->children[0];
	if (!(dterm->terms.size()==1 && dterm->terms[0] == 'd'))
		return false;

	exprgrp* tn = (exprgrp*)dtn->children[1];
	//exprgrp* t;
	exprtree* t;

	if (nval == 1)
		t = tn;
	else if (nval >= 2) {
		if (tn->sid != scg::SUP_EXPR)
			return false;
		
		//t = (exprgrp*)tn->children[0];
		t = tn->children[0];
		
		exprgrp* n = (exprgrp*)tn->children[1]; 
		if (!(n->sid == scg::NUM_EXPR)) // this is the order of the derivative (again)
			return false;
		exprleaf* nterm = (exprleaf*)n->children[0];
		if (!(nterm->terms.size()==1)) // note: 2<=n<=9 is a limitation
			return false;
		std::string nstr = scg::symdb_findsymbol_unicode(nterm->terms[0])->name;
		if (!(nval == atoi(nstr.c_str()))) // both n's must match	
			return false;
	}

	return parseIndependentVariable(t, de);
}

// x = dependent variable
// t = independent variable
// n = order of the derivative
std::string writesagederivative(std::string x, std::string t, int n) {
	
	std::string s = "diff(" + x + "," + t;
	if (n > 1) {
 		char nbuf[2]; // n is one digit (plus termination character)
 		sprintf(nbuf, "%d", n);
 		std::string nstr = nbuf;
		s += "," + nstr;
	}
	s += ")";
	return s;
}

// Determine whether or not the given fraction expression could represent a derivative (e.g. dy/dt).
// If so, update member variables accordingly.
// todo: what to do if we detect a derivative with mistmatched (in)dependent variable. define a flag.
bool OdeParser::parseLeibnizDerivative(exprgrp* tree, int &n, DE* de) {
	
	// looking for expression of the form d^nx/dt^n

	if (!(tree->sid == scg::FRAC_EXPR))
		return false;

	exprgrp* dnx = (exprgrp*)tree->children[0]; // numerator

	if (!parseLeibnizNumerator(dnx, n, de))
		return false;
	
	exprgrp* dtn = (exprgrp*)tree->children[1]; // denominator
	if (!parseLeibnizDenominator(dtn, n, de))
		return false;

	return true;
	
	/////// MOVED THIS:
//  	if (n > de->order) // && !locked)
//  		de->order = n; // this is the highest degree encountered thus far
// 		
// 	// we now have the derivative	
// 	std::string s = writesagederivative(de->dependentVariable, de->independentVariable, n);
// 	tree->setsage(s); // tell tree we found the right sage expression for it (more advanced than what tree->writesage() produces)
// 	
// 	return true;
}

// note: 
// if this returns true, it has set n to the number of primes.
// if this returns false, it has not touched n
bool OdeParser::parsePrimedDependentVariable(exprgrp* tree, int &n, DE* de) {
	
	if (tree->sid != scg::PRIME_EXPR)
		return false;
	
	exprtree* var = tree->children[0];
	if (!parseDependentVariable(var, de))
		return false;
	
	n = ((exprgrp*)tree->children[1])->children.size(); // prime count
	
	return true;
}

// expecting a function expression e.g. y'(x).
// first child: y'
// second child: (x)
bool OdeParser::parsePrimeDerivative(exprgrp* tree, DE* de) {
	
	if (tree->sid != scg::FN_EXPR) 
		return false;
	
	// primed variable, e.g. y'
	exprgrp* primedVar = (exprgrp*)tree->children[0];
	int n;
	if (!parsePrimedDependentVariable(primedVar, n, de))
		return false;
	
	// independent variable, e.g. the x in (x)
	exprtree* indepVar = tree->children[1];
	if (indepVar->sid == scg::PAREN_EXPR)
		indepVar = ((exprgrp*)indepVar)->children[1];
	if (!parseIndependentVariable(indepVar, de))
		return false;
		
	if (n > de->order) // && !locked)
		de->order = n;
	
	// we now have the derivative
	std::string s = writesagederivative(de->dependentVariable, de->independentVariable, n);
	tree->setsage(s); // tell tree we found the right sage expression for it (more advanced than what tree->writesage() produces)

// TODO what if we have the primed variable by itself? (I think this is handled fine)

}

// what about nonlinear terms?
//
// Objectives:
// (1) determine if the expression could be a valid LHS or RHS for a DE (return value)
// (2) Determine if the expression contains a derivative (containsDerivative). If so:
//    (a) Record the maximum degree of such terms
//    (b) Record the dependent variable name.
//    (c) Record the independent variable name, if explicit.
//    (d) Verify that all such terms have the same dependent, and independent (if explicit), variables.
bool OdeParser::parseDEterms(exprgrp* tree, DE* de) {
	
	//containsDerivative = false;
	//std::string discard = "";
	std::string v, i, s;
	
	// the following used in FN_EXPR case
	exprgrp* primedVar;
	exprtree* indepVar;
	int n;
	
	switch (tree->sid) {
		case scg::VAR_EXPR:
			if (getVariableName(tree, v) && v != de->dependentVariable) 				
				indepVarCandidates.push_back(v);
		case scg::NUM_EXPR:
			//tree->writesage(sage);
			return true;

		case scg::REL_EXPR:
			return false;
		
		case scg::ADD_EXPR:
		case scg::MULT_EXPR:
			return parseDEterms((exprgrp*)tree->children[0], de) && parseDEterms((exprgrp*)tree->children[1], de);

		case scg::FN_EXPR: // function name (child 0) can be prime expression
// 			if (parsePrimeDerivative(tree))
// 				return true;
			
			
			// primed variable, e.g. y'
			primedVar = (exprgrp*)tree->children[0];
			if (parsePrimedDependentVariable(primedVar,n, de)) {
				
				// independent variable, e.g. the x in (x)
				exprtree* indepVar = tree->children[1];
				if (indepVar->sid == scg::PAREN_EXPR)
					indepVar = ((exprgrp*)indepVar)->children[1];
				if (!parseIndependentVariable(indepVar, de))
					return false; // primed variable with no valid independent variable
					
				if (n > de->order) // && !locked)
					de->order = n;
				
				// we now have the derivative
				std::string s = writesagederivative(de->dependentVariable, de->independentVariable, n);
				tree->setsage(s); // tell tree we found the right sage expression for it (more advanced than what tree->writesage() produces)
				
				return true;
			}
			
			// some other function
			return parseDEterms((exprgrp*)tree->children[0], de) && parseDEterms((exprgrp*)tree->children[1], de);
		
		case scg::PAREN_EXPR:
			return parseDEterms((exprgrp*)tree->children[1], de);
		
		case scg::ROOT_EXPR: // todo: don't allow derivatives in root?
		case scg::NEG_EXPR:
			return parseDEterms((exprgrp*)tree->children[0], de);
			
		case scg::PRIME_EXPR:
			// treat prime expressions as derivatives
			//int n;
			if (!parsePrimedDependentVariable(tree, n, de)) //parsePrimeDerivative(tree);
				return false;
			// this is a primed variable by itself (independent variable not specified)
			// if we have already discovered the independentVariable, use it
			if (n > de->order) // && !locked)
				de->order = n;
			// may need to decide this at the end (by later discovery or default)
			i = de->independentVariable == "" ? INDEP_VAR_PLACEHOLDER : de->independentVariable;
			s = writesagederivative(de->dependentVariable, i, n);
			tree->setsage(s); // tell tree we found the right sage expression for it (more advanced than what tree->writesage() produces)
			return true;

		case scg::FRAC_EXPR:
			// this could be a derivative
			if (parseLeibnizDerivative(tree, n, de)) {
				if (n > de->order) // && !locked)
					de->order = n; // this is the highest degree encountered thus far
					
				// we now have the derivative	
				// tell tree we found the right sage expression for it (more advanced than what tree->writesage() produces)
				tree->setsage( writesagederivative(de->dependentVariable, de->independentVariable, n)); 
				
				return true;
			}
			// todo: don't allow derivatives in denominator?
			// note: break/return is intentionally missing here
		case scg::SUP_EXPR:
		case scg::SUBSCR_EXPR:
			// todo: don't allow derivatives in exponent?
			// todo: not sure what to do with subscripts. subscript should probably be a NUM_EXPR and nothing else.
			return parseDEterms((exprgrp*)tree->children[0], de) && parseDEterms((exprgrp*)tree->children[1], de);
			
		default:
			return false;
	}
	
	return true;


// 	if (tree->sid == scg::ADD_EXPR)
// 	{
// 		parsedeterm((exprgrp*)tree->children[0], s); // writes term
// 		// write +
// 		s += " + ";
// 		parsedeterms((exprgrp*)tree->children[1], s); //writes terms
// 	}
// 	else
// 	{
// 		parsedeterm(tree, s);
// 	}
	
}

void replace(std::string &str, const std::string &oldStr, const std::string &newStr) {
	size_t pos = 0;
	while((pos = str.find(oldStr, pos)) != std::string::npos) {
		str.replace(pos, oldStr.length(), newStr);
		pos += newStr.length();
	}
}

void OdeParser::chooseIndependentVariable(DE* de) {

	if (de->independentVariable != "") // || locked) 
		return; // already chosen, or not allowed right now
		
	// the most common choices
// 	std::vector<std::string> defaultIndepVar;
// 	defaultIndepVar.push_back("t");
// 	defaultIndepVar.push_back("x");
// 	defaultIndepVar.push_back("r");
// 	bool contains = std::find(defaultIndepVar.begin(), defaultIndepVar.end(), *i) != defaultIndepVar.end();
	
	// if we have no candidates to work with, or the only candidate is not actually a candidate, choose x or t
// 	if (indepVarCandidates.size() == 0 || indepVarCandidates[0] == dependentVariable) { 
// 		// take default
// 		independentVariable = (dependentVariable == "x") ? "t" : "x";
// 		return;
// 	}
	// iterate over candidates, hoping to find an x or a t
	for (std::vector<std::string>::const_iterator i = indepVarCandidates.begin(); i != indepVarCandidates.end(); ++i) {
		if ((*i == "x" || *i == "t") && *i != de->dependentVariable) { 
			de->independentVariable = *i;
			return;
		}
	}
	// if x or t aren't candidates, just take the first one that works
	for (std::vector<std::string>::const_iterator i = indepVarCandidates.begin(); i != indepVarCandidates.end(); ++i) {
		if (*i != de->dependentVariable) { 
			de->independentVariable = *i;
			return;
		}
	}
	// choose x or t by default
	de->independentVariable = (de->dependentVariable == "x") ? "t" : "x";
	
	indepVarCandidates.clear(); // done with this
}
// Try to parse expression tree as a DE. 
bool OdeParser::parseDE(exprgrp *tree, DE* de) {

	//std::string sage = ""; 
	//reset();
	
	// 1. expression must be an equation
	if (tree->sid != scg::REL_EXPR) 
		return false;
		
	static const scg::symbol *op = scg::symdb_findsymbol_unicode(((exprleaf*)tree->children[2])->terms[0]);
	if (!(op->name == "eq")) 
		return false;

	exprgrp* left = (exprgrp*)tree->children[0];  // terms making up LHS of DE
	//bool lhsContainsDerivative = false;
	if (!parseDEterms(left, de))//, sage, lhsContainsDerivative)) 
		return false;
	
	//sage += " == "; generate whole sage string at the end!
	
	exprgrp* right = (exprgrp*)tree->children[1];  // terms making up RHS of DE
	//bool rhsContainsDerivative = false;
	if (!parseDEterms(right, de))//, sage, rhsContainsDerivative))
		return false;
	
	// 2. there must be a derivative term somewhere (either LHS or RHS)
	//if (!(lhsContainsDerivative || rhsContainsDerivative)) 
	if (de->order <= 0)
		return false;
	
	// etc.
		
	// if we've made it this far, tree must represent a valid DE 
	//s += sage;
//if (!locked) {
	tree->writesage(de->sage);
	// TODO: if independentVariable == "" (if the indep. var. isn't explicitly stated)
	// choose x by default (unless it is the dep. variable, in which case choose t)
	// you should also keep a list of candidates: everytime a variable is encountered, add it to a vector
	// if one of them is x or t, use it. otherwise use the first one.
	chooseIndependentVariable(de);
	
//}
	// replace "independentVariable" with actual variable name in sage (requires code snippet)
	replace(de->sage,INDEP_VAR_PLACEHOLDER,de->independentVariable);
	
	return true;
}

// If successful create the corresponding DE object (or sage string?) and return true.
// Otherwise return false.
bool OdeParser::isDifferentialEquation(exprgrp *tree, DE *&de) { //, std::string &s) { //DifferentialEquation *de) {

	de = new DE();	
	bool retVal = parseDE(tree, de);
	if (!retVal) {
		delete de; de = NULL;
	}
	return retVal;
}

// Determine if the expression tree represents a basic arithmetic expression
// (where a "basic expression" is any combination of expression of 
// "basic" types: var, num, add, paren, root, neg, mult, frac, sup).
// constant indicates if we should determine if it is constant with respect to the dependent and independent variables
// (i.e. does not involve the dependent or independent variables).
bool OdeParser::parseBasicExpr(exprtree* tree, bool constant, DE* de) {

	std::string v;
	
	switch(tree->sid) {
	
		case scg::TERMINAL_EXPR:
			if (constant) {
				// if symbol matches one of the variables, return false. else return true.
				//((exprleaf*)tree)->terms
				if (getVariableName(tree,v) && (v == de->dependentVariable || v == de->independentVariable))
					return false;
			} 
			return true;
			
		case scg::VAR_EXPR:
			return parseBasicExpr(((exprgrp*)tree)->children[0], constant, de);
			
		case scg::NUM_EXPR:
			return true;
		
		case scg::ADD_EXPR:
		case scg::MULT_EXPR:
		case scg::FRAC_EXPR:
		case scg::SUP_EXPR:
			return parseBasicExpr(((exprgrp*)tree)->children[0], constant, de) 
			    && parseBasicExpr(((exprgrp*)tree)->children[1], constant, de);

		case scg::PAREN_EXPR:
			return parseBasicExpr(((exprgrp*)tree)->children[1], constant, de);

		case scg::ROOT_EXPR:
		case scg::NEG_EXPR:
			return parseBasicExpr(((exprgrp*)tree)->children[0], constant, de);
			
		default:
			return false;
	}	
}

// overloaded
// never used?
// bool OdeParser::parseBasicExpr(exprtree* tree) {
// 	return parseBasicExpr(tree, false, NULL);
// }

// parse basic expression in parentheses
bool OdeParser::parseBasicParenExpr(exprtree* tree, bool constant, DE* de) {
	return (tree->sid == scg::PAREN_EXPR) && parseBasicExpr(((exprgrp*)tree)->children[1], constant, de);
}

bool OdeParser::parseInitialConditionLHS(exprtree* tree, /*int &order,*/ DE* de, IC* ic) {

	//int order = 0;
	
	// lhs needs to be either    FN     e.g. y'(0) or      MULT  e.g. dy/dt(0) 
	//                          / \                       /  \
	//                 PRIME/VAR  PAREN     FRAC/PRIME/VAR   PAREN 
	exprtree *a, *b;
	switch (tree->sid) {
		
		case scg::FN_EXPR:
			a = ((exprgrp*)tree)->children[0];
			b = ((exprgrp*)tree)->children[1];
			if (parsePrimedDependentVariable((exprgrp*)a, ic->order, de))
				;
			else if (parseDependentVariable(a, de))
				ic->order = 0;
				//return parseBasicParenExpr(((exprgrp*)tree)->children[1], true, de);
			else
				return false;
			break;
			//return parseBasicParenExpr(b, true, de);
// 			return (parsePrimedDependentVariable((exprgrp*)((exprgrp*)tree)->children[0], order, de) 
// 			     || parseDependentVariable(((exprgrp*)tree)->children[0], de))
// 			     && parseBasicParenExpr(((exprgrp*)tree)->children[1], true, de);
			
		case scg::MULT_EXPR:
			a = ((exprgrp*)tree)->children[0];
			b = ((exprgrp*)tree)->children[1];
			if (parseLeibnizDerivative((exprgrp*)a, ic->order, de) 
			 || parsePrimedDependentVariable((exprgrp*)a, ic->order, de))
				;
			else if (parseDependentVariable(a, de))
				ic->order = 0;
			else
				return false;
			break;
			//return parseBasicParenExpr(b, true, de);
// 			return (parseLeibnizDerivative((exprgrp*)((exprgrp*)tree)->children[0], order, de) 
// 			     || parsePrimedDependentVariable((exprgrp*)((exprgrp*)tree)->children[0], order, de) 
// 			     || parseDependentVariable(((exprgrp*)tree)->children[0], de))
// 			     && parseBasicParenExpr(((exprgrp*)tree)->children[1], true, de);
			
		default:
			return false;
	}
	if (parseBasicParenExpr(b, true, de)) {
		ic->initialTimeExpr = (exprgrp*)((exprgrp*)b)->children[1];
		return true;
	}
	return false;
}


 
bool OdeParser::isInitialCondition(exprgrp* tree, /*int &order,*/ DE* de, IC* &ic) {

	// we need an equation and we first need to know the DE it relates to
	if (tree->sid != scg::REL_EXPR || de == NULL)
		return false;
	
	ic = new IC();	

	// for now: 
	// - lhs must be something like y'(0) or dy/dt(0), where 0 is a basic expression not involving variables of DE
	// - rhs must be a basic expression not involving variables of DE
	exprtree* lhs = tree->children[0];
	if (!parseInitialConditionLHS(lhs, /*order,*/ de, ic)) {
		delete ic; ic = NULL;
		return false;
	}
	
	exprtree* rhs = tree->children[1];
	if (parseBasicExpr(rhs, true, de)) {
		ic->initialValueExpr = (exprgrp*)rhs;
		return true;
	}
	delete ic; ic = NULL;
	return false;
}


// Try to parse expression tree as a DE with initial conditions. 
// If successful create the corresponding IVP object and return true.
// Otherwise return false.
bool OdeParser::isInitialValueProblem(exprgrp *tree, IVP* &ivp) { //, std::string &s) { //InitialValueProblem *ivp) {

	// First, we need a MULTI_EXPR
	if (tree->sid != scg::MULTI_EXPR || tree->children.size() != 1)
		return false;

	exprgrp* rows = (exprgrp*)tree->children[0];
	if (rows->sid != scg::MATRIXROWS_EXPR)
		return false;
	
	ivp = new IVP();
	
	ivp->deExpr = NULL;
	ivp->de = NULL;
// 	exprgrp* theDEtree = NULL;
// 	DE* theDE = NULL;
	int maxOrder = 0; // keep track of maximum order (in case more than one DE is encountered, in which case one of them should really be an IC).
	
	
	
	// iterate over all expressions
	for (size_t i = 0; i < rows->children.size(); ++i) {
		if (rows->children[i]->sid != scg::MATRIXROW_EXPR)
			return false;
		exprgrp* row = (exprgrp*)rows->children[i]; 
		for (size_t j = 0; j < row->children.size(); ++j) {
			if (row->children[j]->sid == scg::TERMINAL_EXPR) {
				delete ivp; ivp = NULL;
				return false;
			}
			exprgrp* cell = (exprgrp*)row->children[j];
			
			// "visit" expression
			DE* de;
			if (isDifferentialEquation(cell, de)) {
				if (de->order > maxOrder) {
					maxOrder = de->order;
					// this is our DE
					ivp->de = de;
					ivp->deExpr = cell; 
// 					theDE = de;
// 					theDEtree = cell; 
					
				}
				//break; // found a DE
			}
		}
		
	}		

	if (ivp->de == NULL) {
// 	if (theDE == NULL)
		delete ivp; ivp = NULL;
		return false; // didn't find a DE
	}
		
// 	int numICs = 0; // number of initial conditions	
	
	// now that we have the DE (so we have identified the variables), lock them in place
	// and look for initial conditions
	ivp->initialConds.resize(ivp->de->order, NULL); // initialize them all to NULL. index corresponds to order.
	//locked = true;
	for (size_t i = 0; i < rows->children.size(); ++i) {
		if (rows->children[i]->sid != scg::MATRIXROW_EXPR) {
			delete ivp; ivp = NULL;
			return false;
		}
		exprgrp* row = (exprgrp*)rows->children[i]; 
		for (size_t j = 0; j < row->children.size(); ++j) {
			if (row->children[j]->sid == scg::TERMINAL_EXPR) {
				delete ivp; ivp = NULL;
				return false;
			}
			exprgrp* cell = (exprgrp*)row->children[j];
			if (cell == ivp->deExpr)
				continue; // skip the DE
			
			// "visit" expression
			// use a copy of DE information in case it gets changed during parsing
			DE* decopy = ivp->de->copy();
// 			DE* decopy = theDE->copy();
			IC* ic;
			if (isInitialCondition(cell, decopy, ic)) { 
				// note: redundancy (in terms of order) is allowed, but there can't be any conditions missing in the end
				if (ic->order >= 0 && ic->order < ivp->initialConds.size())
					ivp->initialConds[ic->order] = ic;
//				numICs++;
// 				std::string s;
// 				cell->writesage(s);
// 				logmsg("Found initial condition. order %d:, sage: %s\n", ic->order, s.c_str());
				
			}
			delete decopy;
		}
		
	}
	//locked = false;
	
	// check that we have the right degrees for the initial conditions
	// and that the initial times are all the same
	for (int i = 0; i < ivp->initialConds.size(); i++) {
		if (ivp->initialConds[i] == NULL) {
			delete ivp; ivp = NULL;
			return false;
		}
		exprgrp *a = ivp->initialConds[0]->initialTimeExpr;
		exprgrp *b = ivp->initialConds[i]->initialTimeExpr;
		if (!a->equals(b)) {
			delete ivp; ivp = NULL;
			return false;
		}
	}
	
	return true;
	
	// need a list structure indexed by order, and fill it up with pointers to the initial expressions
// 	bool retVal = numICs == ivp.de->order;
// 	//delete theDE;
// 	return retVal;
}
