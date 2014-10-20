// This file was automatically generated.
// It should not be changed directly; instead, change
// the grammar file and run the "header" program.


#ifndef SCG_GRAMMAR_VALUES_H_
#define SCG_GRAMMAR_VALUES_H_


#include "dlldecl.h"

#include <string>
#include <map>


namespace scg {

// Semantic type identifiers

typedef int SemanticId;

static const SemanticId InvalidSemanticId = -1;

static const SemanticId NULL_EXPR = -3;
static const SemanticId TERMINAL_EXPR = -2;
static const SemanticId BLANK_EXPR = 0;

/*
static const SemanticId BINOP_EXPR = 1000;
static const SemanticId OP_EXPR = 1001;
static const SemanticId FUNC_EXPR = 1002;
static const SemanticId SUB_EXPR = 1003;
static const SemanticId FACTORIAL_EXPR = 1004;
static const SemanticId ARROW_EXPR = 1005;
static const SemanticId SERIES_EXPR = 1006;
static const SemanticId LIMHEAD_EXPR = 1007;
static const SemanticId FUNCNAME_EXPR = 1008;
static const SemanticId LIST_EXPR = 1009;
static const SemanticId LISTHEAD_EXPR = 1010;
static const SemanticId QUANTFORM_EXPR = 1011;
static const SemanticId QUANTLIST_EXPR = 1012;
static const SemanticId QUANT_EXPR = 1013;
static const SemanticId INTEGRATION_EXPR = 1014;
static const SemanticId FLOAT_EXPR = 1015;
static const SemanticId BAREOP_EXPR = 1016;
static const SemanticId SERHEAD_DIAG_EXPR = 1017;
static const SemanticId SERHEAD_VERT_EXPR = 1018;
static const SemanticId INTHEAD_DIAG_EXPR = 1019;
static const SemanticId INTHEAD_VERT_EXPR = 1020;
static const SemanticId NUMBER_EXPR = 1021;
*/

/*static const SemanticId BINOP_EXPR = 1;
static const SemanticId SUM_EXPR = 2;
static const SemanticId MULT_EXPR = 3;
static const SemanticId FACT_EXPR = 4;
static const SemanticId FRAC_EXPR = 5;
static const SemanticId PAREN_EXPR = 6;
static const SemanticId FUNC_EXPR = 7;
static const SemanticId SUP_EXPR = 8;
static const SemanticId SUB_EXPR = 9;
static const SemanticId ROOT_EXPR = 10;
static const SemanticId LIM_EXPR = 11;
static const SemanticId SERIES_EXPR = 12;
static const SemanticId INT_EXPR = 13;
static const SemanticId VAR_EXPR = 14;
static const SemanticId NUM_EXPR = 15;*/


const SemanticId NUM_EXPR = 1;
const SemanticId VAR_EXPR = 2;
const SemanticId NEG_EXPR = 3;
const SemanticId REL_EXPR = 4;
const SemanticId ADD_EXPR = 5;
const SemanticId MULT_EXPR = 6;
const SemanticId FRAC_EXPR = 7;
const SemanticId SUP_EXPR = 8;
const SemanticId SUBSCR_EXPR = 9;
const SemanticId PAREN_EXPR = 10;
const SemanticId ROOT_EXPR = 11;
const SemanticId FN_EXPR = 12;
const SemanticId INTEGRAL_EXPR = 13;
const SemanticId LIMIT_EXPR = 14;
const SemanticId PRIME_EXPR = 15;
const SemanticId SUM_EXPR = 16;
const SemanticId MATRIX_EXPR = 17;
const SemanticId MATRIXROWS_EXPR = 18;
const SemanticId MATRIXROW_EXPR = 19;
const SemanticId DOTS_EXPR = 20;
const SemanticId NAME_EXPR = 21;
const SemanticId LIST_EXPR = 22;
const SemanticId PLUSMINUS_EXPR = 23;
const SemanticId FACTORIAL_EXPR = 24;
const SemanticId PRIMESYMBOL_EXPR = 25;
const SemanticId PRIMESYMBOLS_EXPR = 26;
const SemanticId MULTI_EXPR = 27;

//static const SemanticId SIZE_EXPR = 18;
//static const SemanticId LIST_EXPR = 19;

const SemanticId PLACEHOLDER_EXPR = 1001;


// Children labels


const size_t BIN_OP = 2;
const size_t EXPR_LHS = 0;
const size_t EXPR_RHS = 1;
const size_t FN_ARG = 1;
const size_t FN_NAME = 0;
const size_t FRAC_DENOM = 1;
const size_t FRAC_NUMER = 0;
const size_t INTEGRAL_VAR = 3;
const size_t NAME_VALUE = 0;
const size_t UNARY_TERM = 0;
const size_t NUM_VALUE = 0;
const size_t PAREN_CLOSE = 2;
const size_t PAREN_CONTENTS = 1;
const size_t PAREN_OPEN = 0;
const size_t PRIME_COUNT = 1;
const size_t PRIME_NAME = 0;
const size_t RANGE_LOLIM = 0;
const size_t RANGE_OPERAND = 2;
const size_t RANGE_UPLIM = 1;
const size_t ROOT_CONTENTS = 0;
const size_t SUBSCR_BASE = 0;
const size_t SUBSCR_SUB = 1;
const size_t SUP_BASE = 0;
const size_t SUP_POW = 1;
const size_t VAR_NAME = 0;
const size_t MATRIX_NUMROWS = 0;
const size_t MATRIX_NUMCOLS = 1;
const size_t MATRIX_ROWS = 2;
const size_t LIMIT_VAR = 0;
const size_t LIMIT_APPROACH = 1;
const size_t LIMIT_VALUE = 2;
//const size_t MULTI_COUNT = 0;
const size_t MULTI_ROWS = 0;

SemanticId mktempsid();
SemanticId mktermsid(unsigned short unicode);
unsigned short termsid2unicode(SemanticId sid);

void initialize_grammar_typemap();
void destroy_grammar_typemap();
std::map<std::string, SemanticId> &grammar_typemap();

SemanticId string_to_sid(const std::string &s);
bool sidisterminal(SemanticId sid);
std::string sid_to_string(SemanticId sid);

bool sidistemporary(SemanticId sid);
bool sidismeaningful(SemanticId sid);

}

#endif
