#ifndef EXPR_H_
#define EXPR_H_

#include "mathrecognizer-private.h"
#include "expr-node.h"
#include "symbols.h"
#include <vector>
#include <set>
#include <cstdio>

struct exprtree {
	scg::SemanticId sid;

	exprtree() : sid(scg::InvalidSemanticId) { }
	explicit exprtree(scg::SemanticId sid_) : sid(sid_) { }
	virtual ~exprtree() { }
	virtual bool equals(exprtree* other) = 0;

 	virtual int writesage(std::string &s) = 0;
	virtual int writevardecls(std::string &s, std::set<std::string> &decls) = 0;
	virtual int writestream(std::vector<char> &buf) = 0;
	virtual scg::basic_tree *mkscgtree() = 0;
};

struct exprgrp : public exprtree {
	std::vector<exprtree *> children;

	exprgrp() : exprtree() { sage = ""; }
	explicit exprgrp(scg::SemanticId sid) : exprtree(sid) { }
	~exprgrp();
	bool equals(exprtree* other);

	void addchild(exprtree *child) { children.push_back(child); }
	
	int writesage(std::string &s);
	
	// allow outside user (in particular, OdeParser) to tell this isntance what it's sage code should be
	void setsage(std::string s); 
	
	int writevardecls(std::string &s, std::set<std::string> &decls);
	int writestream(std::vector<char> &buf);
	scg::basic_tree *mkscgtree();
	
private:
	std::string sage;
};

struct exprleaf : public exprtree {
	std::vector<unsigned short> terms;

	exprleaf() : exprtree(scg::TERMINAL_EXPR) { }

	void addterm(unsigned short t) { terms.push_back(t); }
	bool equals(exprtree* other);
	
	int writesage(std::string &s);
	int writevardecls(std::string &s, std::set<std::string> &decls);
	int writestream(std::vector<char> &buf);
	scg::basic_tree *mkscgtree();

private:
	const scg::symbol *getsymbol(unsigned short unicode);
};


exprtree *mkexprtree(const char **stream, size_t *len);
exprtree *mkexprtree(const scg::ExpressionTree *expr);

#endif

