#include "write-expr.h"
#include "grammar-values.h"
#include "etype.h"
#include "expr.h"
#include "cmdcode.h"
#include "log.h"
#include "int32.h"
#include "verb.h"
#include <string.h>
using scg::verb_out;

static char
tr(scg::SemanticId id) {
	switch (id) {
	case scg::BLANK_EXPR: return OP_EMPTY;
	case scg::FN_EXPR: return OP_FUNC;
	case scg::FRAC_EXPR: return OP_FRAC;
	case scg::MULT_EXPR: return OP_MUL;
	case scg::PAREN_EXPR: return OP_PARENS;
	case scg::SUBSCR_EXPR: return OP_SUBSCR;
	case scg::SUP_EXPR: return OP_SUPSCR;
	case scg::NUM_EXPR: return OP_NUM;
	case scg::VAR_EXPR: return OP_NAME;
	case scg::NEG_EXPR: return OP_NEG;
	case scg::INTEGRAL_EXPR: return OP_INTEGRAL;
	case scg::SUM_EXPR: return OP_SUM;
	default: return 0;
	}
}

static int
writeexpr_def(const scg::ExpressionTree *tree, char **buf, size_t len, size_t *wrlen) {
	(*buf)[0] = tr(tree->type());
	++(*buf);
	++(*wrlen);
	for (size_t i = 0; i < tree->nchildren(); ++i) {
		int e = writeexpr(tree->child(i), buf, len, wrlen);
		if (e != 0) {
			return e;
		}
	}
	return 0;
}

int
writeexpr(const scg::ExpressionTree *tree, char **buf, uint32_t len, uint32_t *wrlen) {
	int e = 0;
	if (len == *wrlen) {
		return -2;
	}

	VERBOSE(*verb_out << "cas: writing expression " << tree->str() << " of type " << tree->type() << std::endl);
	switch (tree->type()) {
	case scg::BLANK_EXPR: {
		(*buf)[0] = OP_EMPTY;
		++(*buf);
		++(*wrlen);
		break;
	}
	case scg::PAREN_EXPR: {
		(*buf)[0] = OP_PARENS;
		++(*buf);
		++(*wrlen);
		e = writeexpr(tree->child(scg::PAREN_CONTENTS), buf, len, wrlen);
		break;
	}
	case scg::ADD_EXPR: {
		const scg::ExpressionTree *op = tree->child(scg::ADD_OP);
		if (!strcmp(op->str(), "plus")) {
			(*buf)[0] = OP_ADD;
		}
		else if (!strcmp(op->str(), "horzline")) {
			(*buf)[0] = OP_SUB;
		}
		else {
			return -1;
		}
		++(*buf);
		++(*wrlen);
		e = writeexpr(tree->child(scg::EXPR_LHS), buf, len, wrlen);
		if (e == 0) {
			e = writeexpr(tree->child(scg::EXPR_RHS), buf, len, wrlen);
		}
		break;
	}
	case scg::REL_EXPR: {
		const scg::ExpressionTree *op = tree->child(scg::REL_OP);
		if (!strcmp(op->str(), "eq")) (*buf)[0] = OP_EQ;
		else if (!strcmp(op->str(), "neq")) (*buf)[0] = OP_NEQ;
		else if (!strcmp(op->str(), "lt")) (*buf)[0] = OP_LT;
		else if (!strcmp(op->str(), "leq")) (*buf)[0] = OP_LEQ;
		else if (!strcmp(op->str(), "gt")) (*buf)[0] = OP_GT;
		else if (!strcmp(op->str(), "geq")) (*buf)[0] = OP_GEQ;
		else return -1;
		++(*buf);
		++(*wrlen);
		e = writeexpr(tree->child(scg::EXPR_LHS), buf, len, wrlen);
		if (e == 0) {
			e = writeexpr(tree->child(scg::EXPR_RHS), buf, len, wrlen);
		}
		break;
	}
	case scg::NUM_EXPR: {
		const scg::ExpressionTree *term = tree->child(0);
		if (!strcmp(term->str(), "infin")) {
			(*buf)[0] = OP_INFIN;
			++(*buf);
			++(*wrlen);
			break;
		}
	}
	case scg::VAR_EXPR: {
		const scg::ExpressionTree *term = tree->child(0);
		(*buf)[0] = tr(tree->type());
		++(*buf);
		++(*wrlen);
		size_t nmlen = strlen(term->str());
		if (nmlen >= EXPR_VALUELEN-1 || nmlen + 1 > (len - *wrlen)) {
			return -2;
		}
		(*buf)[0] = nmlen;
		++(*buf);
		++(*wrlen);
		strncpy(*buf, term->str(), nmlen);
		*buf += nmlen;
		*wrlen += nmlen;
		break;
	}

	case scg::MATRIX_EXPR: {
		(*buf)[0] = OP_MATRIX;
		++(*buf);
		++(*wrlen);
		e = writeexpr(tree->child(scg::MATRIX_NUMROWS), buf, len, wrlen);
		if (e == 0) {
			e = writeexpr(tree->child(scg::MATRIX_NUMCOLS), buf, len, wrlen);
			if (e == 0) {
				tree = tree->child(scg::MATRIX_ROWS);
				for (size_t i = 0; i < tree->nchildren(); ++i) {
					e = writeexpr(tree->child(i), buf, len, wrlen);
				}
			}
		}
		break;
	}
	case scg::MATRIXROW_EXPR: {
		for (size_t i = 0; i < tree->nchildren(); ++i) {
			e = writeexpr(tree->child(i), buf, len, wrlen);
		}
		break;
	}
	case scg::INTEGRAL_EXPR: {
		bool loblank = (tree->child(scg::RANGE_LOLIM)->type() == scg::BLANK_EXPR);
		bool hiblank = (tree->child(scg::RANGE_UPLIM)->type() == scg::BLANK_EXPR);
		if (loblank == hiblank) {
			e = writeexpr_def(tree, buf, len, wrlen);
		}
		else {
			e = CASREP_BADEXPR;
		}
		break;
	}
	case scg::SUM_EXPR: {
		const scg::ExpressionTree *lolim = tree->child(scg::RANGE_LOLIM);
		const scg::ExpressionTree *hilim = tree->child(scg::RANGE_UPLIM);
		bool loblank = (lolim->type() == scg::BLANK_EXPR);
		bool hiblank = (hilim->type() == scg::BLANK_EXPR);
		if (loblank || hiblank) {
			e = CASREP_BADEXPR;
			break;
		}
		if (!loblank && lolim->type() != scg::REL_EXPR) {
			e = CASREP_BADEXPR;
			break;
		}
		const scg::ExpressionTree *op = lolim->child(scg::REL_OP);
		if (op->type() != scg::TERMINAL_EXPR || strcmp(op->str(), "eq") != 0) {
			e = CASREP_BADEXPR;
			break;
		}
		const scg::ExpressionTree *var = lolim->child(scg::EXPR_LHS);
		if (var->type() != scg::VAR_EXPR) {
			e = CASREP_BADEXPR;
			break;
		}

		(*buf)[0] = OP_SUM;
		++(*buf);
		++(*wrlen);
		e = writeexpr(lolim->child(scg::EXPR_RHS), buf, len, wrlen);
		if (e != 0) break;
		e = writeexpr(hilim, buf, len, wrlen);
		if (e != 0) break;
		e = writeexpr(tree->child(scg::RANGE_OPERAND), buf, len, wrlen);
		if (e != 0) break;
		e = writeexpr(var, buf, len, wrlen);
		break;
	}
	default: {
		e = writeexpr_def(tree, buf, len, wrlen);
	}
	}
	return e;
}*/

