// This file was automatically generated.
// It should not be changed directly; instead, change
// the grammar file and run the "header" program.


// map type strings to semantic type values

#include <string>
#include <map>

#include "grammar-values.h"
#include "mathrecognizer-private.h"

namespace scg {

std::map<std::string, SemanticId> *grammar_typemap_ = 0;

static const SemanticId TERM_BASE = 10000;
static const SemanticId TEMP_BASE = 1000;

static unsigned ntempsids;

SemanticId mktempsid() { return TEMP_BASE + ntempsids++; }
SemanticId mktermsid(unsigned short unicode) { return TERM_BASE + unicode; }
bool sidisterminal(SemanticId sid) { return sid >= TERM_BASE; }
bool sidistemporary(SemanticId sid) { return TEMP_BASE <= sid && sid < TERM_BASE; }
bool sidismeaningful(SemanticId sid) { return sid != InvalidSemanticId && !sidistemporary(sid); }

unsigned short
termsid2unicode(SemanticId sid) {
	if (sidisterminal(sid)) return sid - TERM_BASE;
	return -1;
}

std::map<std::string, SemanticId> &
grammar_typemap() {
	return *grammar_typemap_;
}

void
initialize_grammar_typemap() {
	ntempsids = 0;
	if (grammar_typemap_) return;
	grammar_typemap_ = new std::map<std::string, SemanticId>;
	(*grammar_typemap_)[std::string("ADD")] = ADD_EXPR;
	(*grammar_typemap_)[std::string("FN")] = FN_EXPR;
	(*grammar_typemap_)[std::string("FRAC")] = FRAC_EXPR;
	(*grammar_typemap_)[std::string("INTEGRAL")] = INTEGRAL_EXPR;
	(*grammar_typemap_)[std::string("MULT")] = MULT_EXPR;
	(*grammar_typemap_)[std::string("NEG")] = NEG_EXPR;
	(*grammar_typemap_)[std::string("NUM")] = NUM_EXPR;
	(*grammar_typemap_)[std::string("PAREN")] = PAREN_EXPR;
	(*grammar_typemap_)[std::string("PRIME")] = PRIME_EXPR;
	(*grammar_typemap_)[std::string("REL")] = REL_EXPR;
	(*grammar_typemap_)[std::string("ROOT")] = ROOT_EXPR;
	(*grammar_typemap_)[std::string("SUBSCR")] = SUBSCR_EXPR;
	(*grammar_typemap_)[std::string("SUM")] = SUM_EXPR;
	(*grammar_typemap_)[std::string("SUP")] = SUP_EXPR;
	(*grammar_typemap_)[std::string("VAR")] = VAR_EXPR;
	(*grammar_typemap_)[std::string("LIMIT")] = LIMIT_EXPR;
	(*grammar_typemap_)[std::string("MATRIX")] = MATRIX_EXPR;
	(*grammar_typemap_)[std::string("MATRIXROWS")] = MATRIXROWS_EXPR;
	(*grammar_typemap_)[std::string("MATRIXROW")] = MATRIXROW_EXPR;
	(*grammar_typemap_)[std::string("DOTS")] = DOTS_EXPR;
	(*grammar_typemap_)[std::string("NAME")] = NAME_EXPR;
	(*grammar_typemap_)[std::string("LIST")] = LIST_EXPR;
	(*grammar_typemap_)[std::string("PLUSMINUS")] = PLUSMINUS_EXPR;
	(*grammar_typemap_)[std::string("FACTORIAL")] = FACTORIAL_EXPR;
	(*grammar_typemap_)[std::string("PRIMESYMBOL")] = PRIMESYMBOL_EXPR;
	(*grammar_typemap_)[std::string("PRIMESYMBOLS")] = PRIMESYMBOLS_EXPR;
	
	(*grammar_typemap_)[std::string("MULTI")] = MULTI_EXPR;

	(*grammar_typemap_)[std::string("_PLACEHOLDER")] = PLACEHOLDER_EXPR;
	(*grammar_typemap_)[std::string("_BLANK")] = BLANK_EXPR;
	(*grammar_typemap_)[std::string("_TERMINAL")] = TERMINAL_EXPR;
	(*grammar_typemap_)[std::string("_NIL")] = NULL_EXPR;

	
	/*(*grammar_typemap_)[std::string("DECOR")] = DECOR_EXPR;
	(*grammar_typemap_)[std::string("SIZE")] = SIZE_EXPR;
	(*grammar_typemap_)[std::string("LIST")] = LIST_EXPR;*/

	/*(*grammar_typemap_)[std::string("OP")] = OP_EXPR;
	(*grammar_typemap_)[std::string("BINOP")] = BINOP_EXPR;
	(*grammar_typemap_)[std::string("FUNC")] = FUNC_EXPR;
	(*grammar_typemap_)[std::string("SUB")] = SUB_EXPR;
	(*grammar_typemap_)[std::string("FACTORIAL")] = FACTORIAL_EXPR;
	(*grammar_typemap_)[std::string("ARROW")] = ARROW_EXPR;
	(*grammar_typemap_)[std::string("SERIES")] = SERIES_EXPR;
	(*grammar_typemap_)[std::string("LIMHEAD")] = LIMHEAD_EXPR;
	(*grammar_typemap_)[std::string("FUNCNAME")] = FUNCNAME_EXPR;
	(*grammar_typemap_)[std::string("LIST")] = LIST_EXPR;
	(*grammar_typemap_)[std::string("LISTHEAD")] = LISTHEAD_EXPR;
	(*grammar_typemap_)[std::string("QUANTFORM")] = QUANTFORM_EXPR;
	(*grammar_typemap_)[std::string("QUANTLIST")] = QUANTLIST_EXPR;
	(*grammar_typemap_)[std::string("QUANT")] = QUANT_EXPR;
	(*grammar_typemap_)[std::string("INTEGRATION")] = INTEGRATION_EXPR;
	(*grammar_typemap_)[std::string("FLOAT")] = FLOAT_EXPR;
	(*grammar_typemap_)[std::string("BAREOP")] = BAREOP_EXPR;
	(*grammar_typemap_)[std::string("NUMBER")] = NUMBER_EXPR;
	(*grammar_typemap_)[std::string("SERHEAD_DIAG")] = SERHEAD_DIAG_EXPR;
	(*grammar_typemap_)[std::string("SERHEAD_VERT")] = SERHEAD_VERT_EXPR;
	(*grammar_typemap_)[std::string("INTHEAD_DIAG")] = INTHEAD_DIAG_EXPR;
	(*grammar_typemap_)[std::string("INTHEAD_VERT")] = INTHEAD_VERT_EXPR;
	(*grammar_typemap_)[std::string("SUBSUP")] = SUBSUP_EXPR;*/
	/*(*grammar_typemap_)[std::string("OP")] = OP_EXPR;
	(*grammar_typemap_)[std::string("BINOP")] = BINOP_EXPR;
	(*grammar_typemap_)[std::string("MULT")] = MULT_EXPR;
	(*grammar_typemap_)[std::string("PAREN")] = PAREN_EXPR;
	(*grammar_typemap_)[std::string("FUNC")] = FUNC_EXPR;
	(*grammar_typemap_)[std::string("SUP")] = SUP_EXPR;
	(*grammar_typemap_)[std::string("SUB")] = SUB_EXPR;
	(*grammar_typemap_)[std::string("ROOT")] = ROOT_EXPR;
	(*grammar_typemap_)[std::string("FRAC")] = FRAC_EXPR;
	(*grammar_typemap_)[std::string("VAR")] = VAR_EXPR;
	(*grammar_typemap_)[std::string("NUM")] = NUM_EXPR;*/

	/*(*grammar_typemap_)[std::string("BINOP")] = BINOP_EXPR;
	(*grammar_typemap_)[std::string("SUM")] = SUM_EXPR;
	(*grammar_typemap_)[std::string("MULT")] = MULT_EXPR;
	(*grammar_typemap_)[std::string("FACT")] = FACT_EXPR;
	(*grammar_typemap_)[std::string("FRAC")] = FRAC_EXPR;
	(*grammar_typemap_)[std::string("PAREN")] = PAREN_EXPR;
	(*grammar_typemap_)[std::string("FUNC")] = FUNC_EXPR;
	(*grammar_typemap_)[std::string("SUP")] = SUP_EXPR;
	(*grammar_typemap_)[std::string("SUB")] = SUB_EXPR;
	(*grammar_typemap_)[std::string("ROOT")] = ROOT_EXPR;
	(*grammar_typemap_)[std::string("LIM")] = LIM_EXPR;
	(*grammar_typemap_)[std::string("SERIES")] = SERIES_EXPR;
	(*grammar_typemap_)[std::string("INT")] = INT_EXPR;
	(*grammar_typemap_)[std::string("VAR")] = VAR_EXPR;
	(*grammar_typemap_)[std::string("NUM")] = NUM_EXPR;*/
}

void
destroy_grammar_typemap()
{
	delete grammar_typemap_;
	grammar_typemap_ = 0;
}

static const std::string NILSTR("nil");
static const std::string TERMSTR("TERM");
static const std::string TEMPSTR("_TEMP");

SemanticId
string_to_sid(const std::string &s) {
	if (s == NILSTR) return InvalidSemanticId;
	else if (s == TERMSTR) return TERMINAL_EXPR;

	std::map<std::string, SemanticId>::const_iterator i = grammar_typemap_->find(s);
	if (i != grammar_typemap_->end()) return i->second;
	const symbol *S = symdb_findsymbol_name(s);
	if (S) {
		return S->sid;
	}
	return InvalidSemanticId;
}

std::string
sid_to_string(SemanticId sid) {
	if (sid == TERMINAL_EXPR) {
		return TERMSTR;
	}
	else if (sidistemporary(sid)) {
		return TEMPSTR;
	}
	else if (sid == InvalidSemanticId) {
		return NILSTR;
	}
	for (std::map<std::string, SemanticId>::const_iterator i = grammar_typemap_->begin(); i != grammar_typemap_->end(); ++i) {
		if (i->second == sid) {
			return i->first;
		}
	}
	if (sidisterminal(sid)) {
		const symbol *S = symdb_findsymbol_sid(sid);
		if (S) {
			return S->name;
		}
	}
	return NILSTR;
}

}
