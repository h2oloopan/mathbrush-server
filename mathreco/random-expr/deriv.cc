#include "derivation.h"
#include "deriv-tools.h"
#include "verb.h"
#include "MathRecognizer.h"

#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <stack>


static std::string rel_names[] = { "AR", "R", "BR", "B", "BN", "Cs", "_" };

/*
static std::map<std::string, std::vector<const scg::term_derivation *> > terminals;

static std::map<const scg::term_derivation *, std::string> aliases;

static std::map<const scg::derivation *, std::string> groups;

static std::map<const scg::derivation *, std::string> heads;
static std::map<const scg::derivation *, std::string> tails;


struct linkspec {
	std::string name1;
	std::string name2;
	std::string link;
};


std::vector<linkspec> links;
*/


struct link_constructor : public scg::derivation_visitor {
	std::string str;

	/*
	void visit(const scg::term_derivation &d)
	{
		std::vector<const scg::term_derivation *> &derivations = terminals[d.t->name];
		derivations.push_back(&d);
		std::stringstream ss;
		ss << d.t->name << '_' << derivations.size();
		std::string alias = ss.str();
		aliases[&d] = alias;
		groups[&d] = alias;
		heads[&d] = alias;
		tails[&d] = alias;
	}*/

	std::stack<bool> current_head;
	std::stack<bool> current_tail;

	link_constructor()
	{
		current_head.push(false);
		current_tail.push(false);
	}

	void visit(const scg::term_derivation &d)
	{
		if (current_head.top()) {
			str += "^";
		}
		str += d.t->name;
		if (current_tail.top()) {
			str += "$";
		}
	}

	void visit(const scg::nt_derivation &d)
	{
		if (d.P->rel != -1) {
			if (current_head.top()) {
				str += "^";
			}
			if (!d.nt->name.empty()) {
				str += "(";
			}

			current_head.push(d.P->defines_head && d.P->head == 0);
			current_tail.push(d.P->defines_tail && d.P->tail == 0);

			d.left->visit(*this);

			current_tail.pop();
			current_head.pop();

			std::stringstream ss;
			ss << " _" << rel_names[d.P->rel] << "_ ";
			str += ss.str();

			current_head.push(d.P->defines_head && d.P->head == 1);
			current_tail.push(d.P->defines_tail && d.P->tail == 1);

			d.right->visit(*this);

			current_tail.pop();
			current_head.pop();

			if (!d.nt->name.empty()) {
				str += ")";
			}
			if (current_tail.top()) {
				str += "$";
			}
		}
		else {
			if (d.left) {
				//current_head.push(d.P->defines_head && d.P->head == 0);
				//current_tail.push(d.P->defines_tail && d.P->tail == 0);
				d.left->visit(*this);
				//current_tail.pop();
				//current_head.pop();
			}
			if (d.right) {
				//current_head.push(d.P->defines_head && d.P->head == 1);
				//current_tail.push(d.P->defines_tail && d.P->tail == 1);
				d.right->visit(*this);
				//current_tail.pop();
				//current_head.pop();
			}
		}
	}


/*
	void visit(const scg::nt_derivation &d)
	{
		std::string &G = groups[&d];
		if (d.left) {
			d.left->visit(*this);

			if (G.empty()) {
				G = groups[d.left];
			}
			else {
				G += ", ";
				G += groups[d.left];
			}
		}

		if (d.right) {
			d.right->visit(*this);
			G += ", ";
			G += groups[d.right];
		}

		heads[&d] = heads[(d.P->head == 0) ? d.left : d.right];
		tails[&d] = tails[(d.P->tail == 0) ? d.left : d.right];

		if (d.P->rel != -1) {
			linkspec L;
			L.link = rel_names[d.P->rel];
			L.name1 = tails[d.left];
			L.name2 = heads[d.right];
			links.push_back(L);

			if (groups[d.right] != heads[d.right]) {
				L.name2 = groups[d.right];
				links.push_back(L);
				if (groups[d.left] != tails[d.left]) {
					L.name1 = groups[d.left];
					links.push_back(L);
				}
			}
			if (groups[d.left] != tails[d.left]) {
				L.name1 = groups[d.left];
				L.name2 = heads[d.right];
				links.push_back(L);
			}

/*
			if (d.P->tail_box == scg::box_selection::symbol) {
				L.name1 = tails[d.left];
			}
			else {
				L.name1 = groups[d.left];
			}

			if (d.P->head_box == scg::box_selection::symbol) {
				L.name2 = heads[d.right];
			}
			else {
				L.name2 = groups[d.right];
			}

			links.push_back(L);
			* /
		}
	}
	*/
};



/*
static void
write_link_table(std::ostream &os)
{
	for (std::map<const scg::term_derivation *, std::string>::const_iterator pa = aliases.begin(); pa != aliases.end(); ++pa) {
		os << pa->first->t->name << " as " << pa->second << std::endl;
	}

	os << std::endl;
	os << "__CUT__\n\n";

	for (std::vector<linkspec>::const_iterator pl = links.begin(); pl != links.end(); ++pl) {
		if (pl->name1.find_first_of(',') != std::string::npos) {
			os << '<' << pl->name1 << "> ";
		}
		else {
			os << pl->name1 << ' ';
		}
		os << pl->link << ' ';
		if (pl->name2.find_first_of(',') != std::string::npos) {
			os << '<' << pl->name2 << ">\n";
		}
		else {
			os << pl->name2 << std::endl;
		}
	}
}*/


static void
write_portable_parse_tree(std::ostream &os, const scg::expression_node *expr, scg::CNFGrammar *g)
{
	if (expr->type() == scg::TERMINAL_EXPR) {
		os << "(TERM " << expr->undecorated_str() << ')';
	}
	else {
		os << '(' << g->semantic_ids[expr->type()] << ' ';
		for (unsigned i = 0; i < expr->nchildren(); ++i) {
			write_portable_parse_tree(os, expr->mchild(i), g);
			if (i + 1 != expr->nchildren()) {
				os << ' ';
			}
		}
		os << ')';
	}
}



int
main(int argc, char *argv[])
{
	initialize_deriver();

	std::ostream *texout = &std::cout;
	std::ostream *linksout = &std::cerr;
	std::ostream *specout = &std::cerr;
	

	unsigned long int seed = std::time(0);
	if (argc > 1) {
		std::stringstream ss;
		for (unsigned i = 1; i < argc; ++i) {
			ss << argv[i] << ' ';
		}

		for (;;) {
			std::string parm;
			ss >> parm;
			if (!ss) break;

			if (parm == "seed") {
				ss >> seed;
			}
			else if (parm == "texout") {
				std::string filename;
				ss >> filename;
				texout = new std::ofstream(filename.c_str());
				if (!*texout) {
					texout = &std::cout;
				}
			}
			else if (parm == "linksout") {
				std::string filename;
				ss >> filename;
				linksout = new std::ofstream(filename.c_str());
				if (!*linksout) {
					linksout = &std::cerr;
				}
			}
			else if (parm == "specout") {
				std::string filename;
				ss >> filename;
				specout = new std::ofstream(filename.c_str());
				if (!*specout) {
					specout = &std::cerr;
				}
			}
			else if (parm == "start") {
				std::string symbol;
				ss >> symbol;
				if (set_derivation_symbol(symbol) < 0) {
					std::cerr << "deriv: " << symbol << " is not a grammar nonterminal\n";
				}
			}
			else if (parm == "pinc") {
				double inc;
				ss >> inc;
				set_term_prob_increment(inc);
			}
			else {
				std::cerr << "deriv: " << parm << " is not a valid parameter name\n";
			}

			if (ss.fail()) {
				std::cerr << "deriv: bad command line at " << parm << std::endl;
			}
		}
	}
	else {
		seed = std::time(0);
	}

	std::cerr << "deriv: using seed " << seed << std::endl;

	set_derivation_seed(seed);

	scg::derivation *d = random_derivation();
	//scg::VerboseOutputToStream(std::cerr);
	//scg::SetVerbosity(1);
	//d->build();

	link_constructor lc;
	d->visit(lc);

	std::cout << *d << std::endl;

	//write_link_table(*linksout);
	*linksout << lc.str << std::endl;

	write_portable_parse_tree(*specout, d->expr, get_derivation_grammar());
	*specout << std::endl;

	
	*texout << "\\documentclass[12pt]{article}\n";
	*texout << "\\pagestyle{empty}\n";
	*texout << "\\begin{document}\n";
	*texout << "\\begin{displaymath}\n";
	*texout << d->expr->long_str() << std::endl;
	*texout << "\\end{displaymath}\n";
	*texout << "\\end{document}\n";

	if (texout != &std::cout) {
		delete texout;
	}
	if (linksout != &std::cerr) {
		delete linksout;
	}
	if (specout != &std::cerr) {
		delete specout;
	}

	shutdown_deriver();
	
	return 0;
}
