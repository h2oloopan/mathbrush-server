#ifndef DERIV_TOOLS_H_
#define DERIV_TOOLS_H_


#include "derivation.h"
#include "grammar.h"


int initialize_deriver();
void shutdown_deriver();

int set_derivation_symbol(const std::string &sym);
int set_term_prob_increment(double inc);
int set_derivation_seed(unsigned long int seed);

scg::derivation *random_derivation();

scg::CNFGrammar *get_derivation_grammar();


#endif
