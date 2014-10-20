#ifndef DERIVATION_H_
#define DERIVATION_H_


#include "grammar.h"
#include "expr-node.h"
#include "parser.h"
#include "builder.h"

#include <ostream>


namespace scg
{


struct derivation_create_flags {
	enum {
		build_immediately,
		defer_build
	};
};


struct derivation {
	const production *P;
	expression_node expr;
	std::vector<derivation *> children;

public:
	explicit derivation(const production *P_) : P(P_), expr(0) { }
	derivation(const production *P_, const expression_node &expr_) : P(P_), expr(expr_) { }
	virtual ~derivation() {
		for (std::vector<derivation *>::iterator i = children.begin(); i != children.end(); ++i) {
			delete *i;
		}
	}

	void add_child(derivation *d) { children.push_back(d); }

	void write(std::ostream &os, unsigned indent = 0) const;

	derivation *copy() const;

	void build();

	bool less(const derivation *rhs) const;
	bool equal(const derivation *rhs) const;
};

std::ostream & operator<<(std::ostream &os, const derivation &d);

//derivation *read_derivation(std::istream &is, const grammar &g);


}


#endif
