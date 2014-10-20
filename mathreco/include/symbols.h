#ifndef SYMBOLS_H_
#define SYMBOLS_H_

#include "group.h"
#include "dlldecl.h"
#include "relation.h"
#include "grammar-values.h"

#include <istream>
#include <ostream>
#include <string>
#include <map>
#include <list>
#include <vector>


namespace scg {

typedef unsigned short unicode_char;

int symdb_init();
int symdb_shutdown();

int symdb_setuserprofile(const char *name);

int symdb_removesym(const std::string &name, bool idx = true);
int symdb_reindex();

struct symbol;
symbol *symdb_firstsymbol();
symbol *symdb_nextsymbol(symbol *cur);

symbol *symdb_findsymbol_name(const std::string &name);
symbol *symdb_findsymbol_unicode(unicode_char unicode);
symbol *symdb_findsymbol_sid(SemanticId sid);

struct prototype {
	symbol *owner;
	RawStrokeGroup strokes;
	NormalizedStrokeGroup nstrokes;
	prototype *next;

	prototype();
};

struct symbol {
	std::string name;
	std::string mathml;
	std::string latex;
	std::string sage;
	unicode_char unicode;
	SemanticId sid;
	rclass_t rclass;

	const nonterminal *nt;

	int attr;
	bool is_normal() const;
	bool is_stacked() const;
	bool is_container() const;
	bool is_small() const;
	bool is_dotted() const;

	prototype *firstproto();
	prototype *nextproto(prototype *cur);
	void clear();

	symbol();
	~symbol();
	prototype *prots;
};

}

#endif
