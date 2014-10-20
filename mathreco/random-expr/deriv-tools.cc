#include "derivation.h"
#include "grammar.h"
#include "parser.h"
#include "utils.h"
#include "iter-tools.h"
#include "mathrecognizer-private.h"
#include "profile.h"
#include "grammar-values.h"
#include "MathRecognizer.h"


#include <map>
#include <string>
#include <set>
#include <vector>
#include <list>
#include <cstdlib>


struct cmp_nt_derivation {
	bool operator()(const scg::nt_derivation *lhs, const scg::nt_derivation *rhs) const
		{ return lhs->less(rhs); }
};

struct cmp_derivation {
	bool operator()(const scg::derivation *lhs, const scg::derivation *rhs) const
		{ return lhs->less(rhs); }
};

typedef std::set<scg::nt_derivation *, cmp_nt_derivation> nt_derivation_set_t;
typedef std::set<scg::derivation *, cmp_derivation> derivation_set_t;

static std::map<scg::CNFGrammar::gid_t, std::map<scg::SemanticId, nt_derivation_set_t> > sid_derivations;
static std::map<scg::CNFGrammar::gid_t, derivation_set_t> terminal_derivations;
static scg::CNFGrammar *cnf = 0;

static scg::CNFGrammar::gid_t start_gid;
static double p_start = 0.0;
static double p_inc = 0.1;


static const scg::CNFProduction *
direct_derivation(const scg::CNFNonTerminal *nt, scg::CNFGrammar::gid_t gid)
{
	for (std::list<scg::CNFProduction *>::const_iterator pp = nt->productions.begin(); pp != nt->productions.end(); ++pp) {
		const scg::CNFProduction *P = *pp;
		if (P->first == gid && P->second == scg::CNFGrammar::INVALID_GID) {
			return P;
		}
	}
	return 0;
}


static int create_flags = scg::derivation_create_flags::defer_build;


static void
identify_reachable_sids(const scg::CNFGrammar &cnf)
{
	std::set<const scg::CNFNonTerminal *> active;

	for (std::map<scg::SemanticId, std::set<scg::CNFGrammar::gid_t> >::const_iterator i = cnf.sid_to_gid_map.begin(); i != cnf.sid_to_gid_map.end(); ++i) {
		scg::SemanticId sid = i->first;
		const std::set<scg::CNFGrammar::gid_t> &gids = i->second;
		for (std::set<scg::CNFGrammar::gid_t>::const_iterator pg = gids.begin(); pg != gids.end(); ++pg) {
			const scg::CNFNonTerminal *nt = cnf.nonterminal(*pg);
			//std::cout << "sid to gid " << sid << " -> " << nt->gid << " (" << nt->name << ")\n";
			/*for (std::list<scg::CNFNonTerminal *>::const_iterator j = nt->parents.begin(); j != nt->parents.end(); ++j) {
				const scg::CNFNonTerminal *parent = *j;
				std::list<scg::CNFProduction *>::const_iterator pp;
				const scg::CNFProduction *P = direct_derivation(parent, nt->gid);
				if (P) {
					sid_derivations[parent->gid][sid].insert(new scg::nt_derivation(P, parent, 0, 0, create_flags));
					active.insert(parent);
				}
			}*/

			for (std::list<scg::CNFProduction *>::const_iterator pp = nt->productions.begin(); pp != nt->productions.end(); ++pp) {
				const scg::CNFProduction *P = *pp;
				if (P->type == sid) {
					scg::nt_derivation *d = new scg::nt_derivation(P, nt, 0, 0, create_flags);
					//std::cout << "add derivation " << *d << std::endl;
					sid_derivations[nt->gid][sid].insert(d);
					active.insert(nt);
				}
			}
		}
	}

	while (!active.empty()) {
		const scg::CNFNonTerminal *nt = *active.begin();
		active.erase(active.begin());
		
		//std::cout << "active is " << nt->name << " (" << nt->gid << ")\n";

		for (std::list<scg::CNFNonTerminal *>::const_iterator i = nt->parents.begin(); i != nt->parents.end(); ++i) {
			const scg::CNFProduction *P = direct_derivation(*i, nt->gid);
			if (P) {
				std::map<scg::SemanticId, nt_derivation_set_t> &destmap = sid_derivations[(*i)->gid];
				std::map<scg::SemanticId, nt_derivation_set_t> &src = sid_derivations[nt->gid];
				for (std::map<scg::SemanticId, nt_derivation_set_t>::const_iterator psrc = src.begin(); psrc != src.end(); ++psrc) {
					nt_derivation_set_t &dest = destmap[psrc->first];
					for (nt_derivation_set_t::const_iterator pd = psrc->second.begin(); pd != psrc->second.end(); ++pd) {
						scg::nt_derivation *d = new scg::nt_derivation(P, *i, *pd, 0, create_flags);
						//std::cout << "add derivation " << *d << std::endl;
						dest.insert(d);
						active.insert(*i);
					}
				}
			}
		}
	}
}


static void
identify_terminal_derivers(const scg::CNFGrammar &cnf)
{
	std::set<scg::CNFNonTerminal *> active;

	for (std::vector<scg::CNFTerminal>::const_iterator i = cnf.terminals.begin(); i != cnf.terminals.end(); ++i) {
		const scg::CNFTerminal &t = *i;
		for (std::list<scg::CNFNonTerminal *>::const_iterator j = t.parents.begin(); j != t.parents.end(); ++j) {
			const scg::CNFProduction *P = direct_derivation(*j, t.gid);
			if (P) {
				scg::derivation *d = new scg::nt_derivation(P, *j, new scg::term_derivation(&*i, P, create_flags), 0, create_flags);
				//std::cout << "new terminal derivation " << *d << " for " << (*j)->gid << " has tree " << *d->expr << std::endl;
				terminal_derivations[(*j)->gid].insert(d);
				active.insert(*j);
			}
		}
	}

	while (!active.empty()) {
		scg::CNFNonTerminal *nt = *active.begin();
		active.erase(active.begin());

		for (std::list<scg::CNFNonTerminal *>::const_iterator i = nt->parents.begin(); i != nt->parents.end(); ++i) {
			const scg::CNFProduction *P = direct_derivation(*i, nt->gid);
			if (P) {
				derivation_set_t &dest = terminal_derivations[(*i)->gid];
				derivation_set_t &src = terminal_derivations[nt->gid];
				for (derivation_set_t::const_iterator pd = src.begin(); pd != src.end(); ++pd) {
					dest.insert(new scg::nt_derivation(P, *i, *pd, 0, create_flags));
				}
				active.insert(*i);
			}
		}
	}
}


static void
dump_terminal_derivations(std::ostream &os, const scg::CNFGrammar &cnf)
{
	os << "terminal derivations:";

	for (std::map<scg::CNFGrammar::gid_t, derivation_set_t>::const_iterator i = terminal_derivations.begin(); i != terminal_derivations.end(); ++i) {
		os << " from " << cnf.nonterminal(i->first)->name << " (gid " << i->first << "): ";
		for (derivation_set_t::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
			//os << cnf.terminal(gid)->name << " (gid " << gid << ") ";
			os << "  " << **j << std::endl;
		}
		os << std::endl;
	}
}



static double
increment_term_prob(double p)
{
	return p + p_inc;
}


static scg::derivation * random_derivation(const scg::CNFGrammar &cnf, scg::CNFGrammar::gid_t gid, double term_prob);

static scg::derivation *
random_sid_derivation(const scg::CNFGrammar &cnf, scg::CNFGrammar::gid_t gid, double term_prob)
{
	//std::cout << "sid derivation for gid " << gid << " nt " << cnf.nonterminal(gid)->name << std::endl;

	std::map<scg::SemanticId, nt_derivation_set_t> &sids = sid_derivations[gid];
	if (sids.empty() || term_prob >= 1.0) {
		return random_derivation(cnf, gid, term_prob);
	}

	unsigned s = std::rand() % sids.size();
	//std::cout << "out of " << sids.size() << " sids available, I'm choosing " << s << std::endl;
	std::map<scg::SemanticId, nt_derivation_set_t>::iterator i = scg::inc_iter(sids.begin(), s);
	nt_derivation_set_t &derivations = i->second;
	unsigned d = std::rand() % derivations.size();
	//std::cout << "out of " << derivations.size() << " derivations available, I'm choosing " << d << std::endl;
	scg::nt_derivation *left = dynamic_cast<scg::nt_derivation *>((*scg::inc_iter(derivations.begin(), d))->copy());
	scg::nt_derivation *tail = left;
	while (tail->left) {
		tail = dynamic_cast<scg::nt_derivation *>(tail->left);
	}

	if (scg::CNFGrammar::is_terminal(tail->P->first)) {
		tail->left = new scg::term_derivation(cnf.terminal(tail->P->first), tail->P, create_flags);
	}
	else if (scg::CNFGrammar::is_nonterminal(tail->P->first)) {
		tail->left = random_derivation(cnf, tail->P->first, increment_term_prob(term_prob));
	}

	if (scg::CNFGrammar::is_terminal(tail->P->second)) {
		tail->right = new scg::term_derivation(cnf.terminal(tail->P->second), tail->P, create_flags);
	}
	else if (scg::CNFGrammar::is_nonterminal(tail->P->second)) {
		tail->right = random_derivation(cnf, tail->P->second, increment_term_prob(term_prob));
	}

	return left;
}


static scg::derivation *
random_derivation(const scg::CNFGrammar &cnf, scg::CNFGrammar::gid_t gid, double term_prob)
{
	const scg::CNFNonTerminal *nt = cnf.nonterminal(gid);
	if (!nt) {
		abort();
	}

	//std::cout << "derivation for gid " << gid << " nt " << cnf.nonterminal(gid)->name << std::endl;

	const derivation_set_t &terminals = terminal_derivations[gid];
	if (!terminals.empty()) {
		double p = (double)std::rand() / RAND_MAX;
		if (p < term_prob) {
			unsigned n = std::rand() % terminals.size();
			//std::cout << "out of " << terminals.size() << " terminals available, I'm choosing " << n << std::endl;
			derivation_set_t::const_iterator i = scg::inc_iter(terminals.begin(), n);
			scg::derivation *d = (*i)->copy();
			//std::cout << "copied derivation " << *d->expr << "(" << d->expr->long_str() << ")\n";
			return d;
		}
	}

	unsigned n = std::rand() % nt->productions.size();
	//std::cout << "out of " << nt->productions.size() << " productions available, I'm choosing " << n << std::endl;
	const scg::CNFProduction *P = *scg::inc_iter(nt->productions.begin(), n);
	
	if (!scg::CNFGrammar::is_valid_gid(P->second)) {
		if (scg::CNFGrammar::is_terminal(P->first)) {
			return new scg::term_derivation(cnf.terminal(P->first), P, create_flags);
		}
		else {
			return new scg::nt_derivation(P, nt, random_sid_derivation(cnf, P->first, term_prob), 0, create_flags);
		}
	}
	else {
		double next_p = increment_term_prob(term_prob);
		scg::derivation *d1, *d2;
		if (scg::CNFGrammar::is_terminal(P->first)) {
			d1 = new scg::term_derivation(cnf.terminal(P->first), P, create_flags);
		}
		else {
			d1 = random_sid_derivation(cnf, P->first, next_p);
		}

		//std::cout << "got left derivation " << *d1->expr << "(" << d1->expr->long_str() << ")\n";
		//std::cout << "got left derivation " << *d1 << std::endl;

		if (scg::CNFGrammar::is_terminal(P->second)) {
			d2 = new scg::term_derivation(cnf.terminal(P->second), P, create_flags);
		}
		else {
			d2 = random_sid_derivation(cnf, P->second, next_p);
		}

		//std::cout << "got right derivation " << *d2->expr << "(" << d2->expr->long_str() << ")\n";
		//std::cout << "got right derivation " << *d2 << std::endl;

		return new scg::nt_derivation(P, nt, d1, d2, create_flags);
	}

	return 0;
}


int
initialize_deriver()
{
	std::string training_path;
	scg::GetTrainingPath(training_path);
	std::ifstream ifs((training_path+"/training.grammar").c_str());
	if (!ifs.is_open() || !ifs) {
		return 0;
	}

	scg::initialize_grammar_typemap();

	scg::Grammar g;
	ifs >> g;
	cnf = g.make_cnf(scg::grammar_conversion_flags::allow_barren_symbols);

	//std::cout << *cnf << std::endl;

	start_gid = cnf->start_gid;

	identify_terminal_derivers(*cnf);
	identify_reachable_sids(*cnf);

	/*
	for (std::map<scg::CNFGrammar::gid_t, std::map<scg::SemanticId, nt_derivation_set_t> >::const_iterator i = sid_derivations.begin(); i != sid_derivations.end(); ++i) {
		std::cout << "gid " << i->first << " (" << cnf->nonterminal(i->first)->name << ")\n";
		for (std::map<scg::SemanticId, nt_derivation_set_t>::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
			std::cout << " sid " << j->first << std::endl;
			for (nt_derivation_set_t::const_iterator k = j->second.begin(); k != j->second.end(); ++k) {
				std::cout << **k << std::endl;
			}
		}
	}
	*/

	return 0;
}

void
shutdown_deriver()
{
	delete cnf;
	scg::destroy_grammar_typemap();
}

int
set_derivation_symbol(const std::string &sym)
{
	const scg::CNFNonTerminal *nt = cnf->nonterminal(sym);
	if (!nt) {
		return E_NOTFOUND;
	}
	else {
		start_gid = nt->gid;
	}
	
	return 0;
}

int
set_derivation_seed(unsigned long int seed)
{
	std::srand(seed);
	return 0;
}

int
set_term_prob_increment(double inc)
{
	p_inc = inc;
	return 0;
}


scg::derivation *
random_derivation()
{
	scg::derivation *d = random_sid_derivation(*cnf, start_gid, p_start);
	//std::cout << "generated derivation\n" << *d << std::endl;
	//scg::VerboseOutputToStream(std::cout);
	//scg::SetVerbosity(1);
	d->build();
	//scg::SetVerbosity(0);
	return d;
}

scg::CNFGrammar *
get_derivation_grammar()
{
	return cnf;
}
