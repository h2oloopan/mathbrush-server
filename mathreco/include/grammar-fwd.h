#ifndef GRAMMAR_FWD_H_
#define GRAMMAR_FWD_H_

#include <cstdlib> // size_t

namespace scg
{


typedef int SemanticId;

struct nonterminal;
struct production;
struct grammar;

struct string_builder;
struct tree_builder;


struct grammar_hook {
	virtual ~grammar_hook() { }
	virtual void start_tree(const nonterminal *, SemanticId) = 0;
	virtual void add_sid_child(SemanticId) = 0;
	virtual void add_ref_child(const nonterminal *) = 0;
	virtual void add_subref_child(const nonterminal *, size_t) = 0;
	virtual void end_tree() = 0;
};

void set_grammar_hook(grammar_hook *ghook);

}



#endif
