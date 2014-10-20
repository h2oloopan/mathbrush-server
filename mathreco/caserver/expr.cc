#include "expr.h"
#include "grammar-values.h"
#include "symbols.h"
#include "MathRecognizer.h"
#include "grammar.h"
#include "cmdcode.h"
#include "log.h"
#include "expr-node.h"
#include "int32.h"
#include <cstdio>
#include <algorithm>
#include <string>

exprgrp::~exprgrp() {
	for (std::vector<exprtree *>::iterator i = children.begin(); i != children.end(); ++i) {
		delete *i;
	}
}

bool exprgrp::equals(exprtree* other) {
	if (other == NULL || sid != other->sid)
		return false;
	if (children.size() != ((exprgrp*)other)->children.size())
		return false;
	std::vector<exprtree *>::const_iterator i,j;
	for (i = children.begin(),
	     j = ((exprgrp*)other)->children.begin(); 
	     i != children.end() && j != ((exprgrp*)other)->children.end(); ++i, ++j) {
		if (!(*i)->equals(*j))
			return false;
	}
	return true;
}

bool exprleaf::equals(exprtree* other) {
	if (other == NULL || sid != other->sid)
		return false;
	if (terms.size() != ((exprleaf*)other)->terms.size())
		return false;
	std::vector<unsigned short>::const_iterator i, j;
	for (i = terms.begin(),
	     j = ((exprleaf*)other)->terms.begin(); 
	     i != terms.end() && j != ((exprleaf*)other)->terms.end(); ++i, ++j) {
		if (*i != *j)
			return false;
	}
	return true;
}

const scg::symbol *
exprleaf::getsymbol(unsigned short unicode) {
	return scg::symdb_findsymbol_unicode(unicode);
}

int
exprleaf::writesage(std::string &s) {
	for (std::vector<unsigned short>::const_iterator i = terms.begin(); i != terms.end(); ++i) {
		const scg::symbol *S = getsymbol(*i);
		if (!S) return CASREP_BADEXPR;
		s += S->sage;
	}
	return CASREP_OK;
}

static bool
isnamereserved(const std::string &s) {
	const static std::string RESERVED[] = {
		"pi", "e", "sin", "cos", "tan", "exp", "erf", "ln", "log"
	};
	const static size_t NRESERVED = sizeof(RESERVED)/sizeof(*RESERVED);
	return std::find(RESERVED, RESERVED+NRESERVED, s) != RESERVED+NRESERVED;
}

int
exprleaf::writevardecls(std::string &s, std::set<std::string> &decls) {
	for (std::vector<unsigned short>::const_iterator i = terms.begin(); i != terms.end(); ++i) {
		const scg::symbol *S = getsymbol(*i);
		if (!S) return CASREP_BADEXPR;
		const std::string &str = S->sage;
		if (isnamereserved(str)) return CASREP_BADEXPR;
		if (decls.find(str) == decls.end()) {
			const char *cstr = str.c_str();			
			logmsg("%s = var('%s');\n", cstr, cstr);
			s += str + " = var('" + str + "');\n";
			decls.insert(str);
		}
	}
	return CASREP_OK;
}

int
exprgrp::writevardecls(std::string &s, std::set<std::string> &decls) {
	for (std::vector<exprtree *>::const_iterator i = children.begin(); i != children.end(); ++i) {
		exprtree *child = *i;
		if (child->sid == scg::TERMINAL_EXPR && sid != scg::VAR_EXPR) continue;		
		int e = (*i)->writevardecls(s, decls);
		if (e != CASREP_OK)
			return e;
	}
	return CASREP_OK;
}

static int
sageunarycmd(const char *cmd, exprtree *tree, std::string &s) {
	s += std::string(cmd) + "(";
	int e = tree->writesage(s);
	s += ")";
	return e;
}

static int
sagebinarycmd(const char *cmd, exprgrp *tree, std::string &s) {
	s += "("; 
	int e = tree->children[0]->writesage(s);
	if (e == CASREP_OK) {
		s += ")" + std::string(cmd) + "(";
		if (e == CASREP_OK)
			e = tree->children[1]->writesage(s);
		s += ")";
	}
	return e;
}

void
exprgrp::setsage(std::string s) {
	sage = s;
}


int
exprgrp::writesage(std::string &s) {
	if (sage != "") { // see setsage method
		s += sage;
		//sage = ""; // one-time use
		return CASREP_OK;
	}
	int e = CASREP_OK;
	switch (sid) {
	case scg::NEG_EXPR: 
		e = sageunarycmd("-", children[0], s);
		break;
	case scg::ROOT_EXPR:
		e = sageunarycmd("sqrt", children[0], s);
		break;
	case scg::FACTORIAL_EXPR:
		e = sageunarycmd("factorial", children[0], s);
		break;

	case scg::NUM_EXPR:
	case scg::NAME_EXPR:
	case scg::VAR_EXPR:
		e = sageunarycmd("", children[0], s);
		break;
		
	case scg::FN_EXPR:
		e = children[0]->writesage(s); // function name
		if (e == CASREP_OK) {
			s += "(";
			e = children[1]->writesage(s); // function operand
			s += ")";
		}
		break;
		
	case scg::PAREN_EXPR:
		s += "(";
		e = sageunarycmd("", children[1], s);
		s += ")";
		break;

	case scg::MULT_EXPR:
		e = sagebinarycmd("*", this, s);
		break;
	case scg::FRAC_EXPR:
		e = sagebinarycmd("/", this, s);
		break;
	case scg::SUP_EXPR:
		e = sagebinarycmd("^", this, s); 
		break;
	case scg::SUBSCR_EXPR:
		e = sagebinarycmd("_", this, s);
		break;

	case scg::ADD_EXPR:
	case scg::REL_EXPR:
		e = children[0]->writesage(s);
		if (e == CASREP_OK) {
			e = children[2]->writesage(s);
			if (e == CASREP_OK) {
				e = children[1]->writesage(s);
			}
		}
		break;

	case scg::MATRIX_EXPR: {
		if (children.size() != 3) {
			e = CASREP_BADEXPR;
			break;
		}
		exprtree *rows = children[2];
		assert(rows->sid == scg::MATRIXROWS_EXPR);
		s += "Matrix([";
		e = rows->writesage(s);
		if (e == CASREP_OK) {
			s += "])";
		}
		break;
	}
	case scg::MATRIXROWS_EXPR: {
		for (size_t i = 0; i < children.size(); ++i) {
			s += "[";
			assert(children[i]->sid == scg::MATRIXROW_EXPR);
			e = children[i]->writesage(s);
			if (e != CASREP_OK) {
				break;
			}
			s += "]";
			if (i < children.size()-1) {
				s += ",";
			}
		}
		break;
	}
	case scg::MATRIXROW_EXPR: {
		for (size_t i = 0; i < children.size(); ++i) {
			s += "(";
			e = children[i]->writesage(s);
			if (e != CASREP_OK) {
				break;
			}
			s += ")";
			if (i < children.size()-1) {
				s += ",";
			}
		}
		break;
	}
	case scg::INTEGRAL_EXPR: {
		if (children.size() != 4) {
			e = CASREP_BADEXPR;
			break;
		}
		exprtree *lolim = children[0];
		exprtree *hilim = children[1];
		exprtree *integrand = children[2];
		exprtree *var = children[3];
		int loempty, hiempty;
		assert(var->sid == scg::VAR_EXPR);
		s += "integral(";
		e = integrand->writesage(s);
		s += ", ";
		e = var->writesage(s);
		loempty = (lolim->sid == scg::BLANK_EXPR);
		hiempty = (hilim->sid == scg::BLANK_EXPR);
		if (loempty != hiempty) {
			e = CASREP_EXPR;
			break;
		}
		else if (!loempty) {
			s += ", ";
			e = lolim->writesage(s);
			s += ", ";
			e = hilim->writesage(s);
		}
		s += ")";
		break;
	}
	case scg::SUM_EXPR: { // example format: sum(n^2, n, l, 10)
		if (children.size() != 3) {
			e = CASREP_BADEXPR;
			break;
		}
		if (children[0]->sid != scg::REL_EXPR) {
			e = CASREP_INVALIDCMD; // unsupported summation expression
			break;
		}
		exprgrp *sub = (exprgrp*)children[0]; // subscript, as REL_EXPR, e.g. i=0
		exprtree *hilim = children[1];
		exprtree *summand = children[2];
		
		// identify the index variable
		exprtree *var = sub->children[0];
		exprleaf *eq = (exprleaf*)sub->children[2]; // terminal indicating relation
		if (var->sid != scg::VAR_EXPR || eq->terms[0] != '=') {
			e = CASREP_INVALIDCMD; // unsupported summation expression
			break;
		}
		exprtree *lolim = sub->children[1];
		
		// now we can proceed		
		s += "sum(";
		e = summand->writesage(s);
		s += ", ";
		e = var->writesage(s);
		int loempty, hiempty;
		loempty = (lolim->sid == scg::BLANK_EXPR);
		hiempty = (hilim->sid == scg::BLANK_EXPR);
		if (loempty || hiempty) {
			e = CASREP_BADEXPR;
			break;
		}
		s += ", ";
		e = lolim->writesage(s);
		s += ", ";
		e = hilim->writesage(s);
		s += ")";
		break;
	}
	case scg::LIST_EXPR: { // multiple expressions
		s += "[";
		for (size_t i = 0; i < children.size(); ++i) {
			s += "(";
			e = children[i]->writesage(s);
			if (e != CASREP_OK) {
				break;
			}
			s += ")";
			if (i < children.size()-1) {
				s += ",";
			}
		}
		s += "]";
		break;
	}
	case scg::MULTI_EXPR: { // multiple expressions
		if (children.size() != 1) {
			e = CASREP_BADEXPR;
			break;
		}
		exprgrp *rows = (exprgrp*)children[0];
		assert(rows->sid == scg::MATRIXROWS_EXPR);
		//e = rows->writesage(s);
		s += "[";
		for (size_t i = 0; i < rows->children.size(); ++i) {
			assert(rows->children[i]->sid == scg::MATRIXROW_EXPR);
			e = rows->children[i]->writesage(s);
			if (e != CASREP_OK) {
				break;
			}
			if (i < rows->children.size()-1) {
				s += ",";
			}
		}		
		s += "]";
		break;
	}
	default:
		e = CASREP_SYSERR;
	}
	return e;
}

exprtree *
mkexprtree(const char **stream, size_t *len) {
	if (*len < 2) {
		return 0;
	}
	char sid = **stream;
	++(*stream);
	--(*len);
	char nchildren = **stream;
	++(*stream);
	--(*len);
	if (*len < 2*(unsigned)nchildren) {
		return 0;
	}
	if ((scg::SemanticId)sid == scg::TERMINAL_EXPR) {
		exprleaf *leaf = new exprleaf;
		while (nchildren--) {
			int16_t val = *(int16_t *)(*stream);
			*stream += 2;
			*len -= 2;
			leaf->addterm((unsigned short)ntohs(val));
		}
		return leaf;
	}
	else {
		exprgrp *grp = new exprgrp(sid);
		while (nchildren--) {
			exprtree *child = mkexprtree(stream, len);
			if (!child) {
				delete grp;
				return 0;
			}
			grp->addchild(child);
		}
		return grp;
	}
}

int
exprleaf::writestream(std::vector<char> &buf) {
	/*if (terms.size() > (unsigned char)~0 || 2 + 2*terms.size() > len) {
		return CASREP_TOOLARGE;
	}*/
	if (terms.size() > (unsigned char)~0) return CASREP_TOOLARGE;
	buf.push_back((char)scg::TERMINAL_EXPR);
	buf.push_back((char)terms.size());
	for (std::vector<unsigned short>::const_iterator i = terms.begin(); i != terms.end(); ++i) {
		//((int16_t *)(*buf))[0] = htons((int16_t)*i);
		int16_t wc = htons((int16_t)*i);
		buf.push_back(((char *)&wc)[0]);
		buf.push_back(((char *)&wc)[1]);
	}

	return 0;
}

int
exprgrp::writestream(std::vector<char> &buf) {
	/*if (children.size() > (unsigned char)~0 || 2 + 2*children.size() > len) {
		return CASREP_TOOLARGE;
	}*/
	if (children.size() > (unsigned char)~0) return CASREP_TOOLARGE;
	buf.push_back((char)sid);
	buf.push_back((char)children.size());
	for (std::vector<exprtree *>::const_iterator i = children.begin(); i != children.end(); ++i) {
		int e = (*i)->writestream(buf);
		if (e != CASREP_OK) {
			return e;
		}
	}

	return 0;
}


scg::basic_tree *
exprgrp::mkscgtree() {
	const scg::grammar &G = *scg::GetMathGrammar();

	if (sid == scg::BLANK_EXPR) return scg::mkblank();
	const scg::production *P = G.getcanonicalsidproduction(sid);
	if (!P) {
		return 0;
	}
	scg::basic_tree *tree = new scg::invented_tree(P);
	for (std::vector<exprtree *>::const_iterator i = children.begin(); i != children.end(); ++i) {
		scg::basic_tree *child = (*i)->mkscgtree();
		if (!child) {
			tree->release();
			return 0;
		}
		tree->addchild(child);
	}
	return tree;
}


scg::basic_tree *
exprleaf::mkscgtree() {
	std::stringstream str;
	std::stringstream mathml;
	std::stringstream latex;
	std::wstring ws(terms.size(), 0);
	for (size_t i = 0; i < terms.size(); ++i) {
		ws[i] = (wchar_t)terms[i];
		const scg::symbol *S = getsymbol(terms[i]);
		if (!S) {
			return 0;
		}
		str << S->name;
		mathml << S->mathml;
		latex << S->latex;
	}
	return new scg::terminal_tree(0, str.str(), ws, mathml.str(), latex.str());
}

exprtree *
mkexprtree(const scg::ExpressionTree *expr) {
	if (expr->type() == scg::TERMINAL_EXPR) {
		exprleaf *leaf = new exprleaf;
		const wchar_t *s = expr->wstr();
		if (!s) {
			return 0;
		}
		while (*s) {
			leaf->addterm((unsigned short)*s);
			++s;
		}
		return leaf;
	}
	else {
		exprgrp *grp = new exprgrp(expr->type());
		for (size_t i = 0; i < expr->nchildren(); ++i) {
			exprtree *child = mkexprtree(expr->child(i));
			if (!child) {
				delete grp;
				return 0;
			}
			grp->addchild(child);
		}
		return grp;
	}
}
