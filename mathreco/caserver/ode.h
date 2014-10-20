#include "expr.h"
#include "symbols.h"
#include <string>
#include <vector>

// these structures should represent just enough information to be able to generate the sage code

// struct derivative {
// 	unsigned short order;
// 	std::string dependentVariable;
// 	std::string independentVariable;
// 	// include enum for notation (leibniz/prime)? not semantically relevant
// 	int writesage(std::string &s); 
// };
// 
struct DE {
 	short order;
 	std::string dependentVariable;
 	std::string independentVariable;
 	std::string sage;
	// exprgrp* tree; ?
	DE()
	{ 
		order = -1; // set to 0 or -1 if no derivative found
		dependentVariable = "";
		independentVariable = "";
		sage = "";
	}
	DE* copy()
	{
		DE* copy = new DE;
		copy->order = order;
		copy->dependentVariable = dependentVariable;
		copy->independentVariable = independentVariable;
		copy->sage = sage;
		return copy;
	}
};

// TODO implement destructors!
struct IC {
	int order;
	exprgrp* initialTimeExpr;
	exprgrp* initialValueExpr;
};

struct IVP {
	DE* de;
	exprgrp* deExpr;
	std::vector<IC*> initialConds;
	//exprgrp* initialTimeExpr;
	//std::vector<exprgrp*> initialValueExprs;
	~IVP() {
		for (std::vector<IC*>::iterator i = initialConds.begin(); i != initialConds.end(); ++i)
			delete *i;
	}
};



class OdeParser {
	
public:
	static bool isDifferentialEquation(exprgrp *tree, DE *&de);
	static bool isInitialCondition(exprgrp* tree, DE* de, IC* &ic);
	static bool isInitialValueProblem(exprgrp *tree, IVP* &ivp);
// 	static short order; // set to 0 or -1 if no derivative found
// 	static std::string dependentVariable;
// 	static std::string independentVariable;
// 	static std::string sage;


private: 
	OdeParser(); // private constructor (no need to create an instance)
	
	//static void reset();
	
	static bool parseDEterms(exprgrp* tree, DE* de);
	static bool parseDE(exprgrp *tree, DE* de);

	static void chooseIndependentVariable(DE* de);
	static bool parseIndependentVariable(exprtree* tree, DE* de);
	static bool parseDependentVariable(exprtree* tree, DE* de);
	
	static bool parsePrimedDependentVariable(exprgrp* tree, int &nval, DE* de);
	static bool parsePrimeDerivative(exprgrp* tree, DE* de);
	
	static bool parseLeibnizDerivative(exprgrp* tree, int &order, DE* de);
	static bool parseLeibnizNumerator(exprgrp* tree, int &nval, DE* de);
	static bool parseLeibnizDenominator(exprgrp* tree, int nval, DE* de);
	
	static bool parseBasicParenExpr(exprtree* tree, bool constant, DE* de);
	static bool parseBasicExpr(exprtree* tree, bool constant, DE* de);
	//static bool parseBasicExpr(exprtree* tree);
	
	static bool parseInitialConditionLHS(exprtree* tree, DE* de, IC* ic);
};
