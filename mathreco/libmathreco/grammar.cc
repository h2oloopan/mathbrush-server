#include "grammar.h"
#include "algorithm.h"
#include "builder.h"
#include "error.h"
#include "rect.h"
#include "relation.h"
#include "utils.h"
#include "symbols.h"
#include "mathrecognizer-private.h"
#include "strutils.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cctype>
#include <cassert>
#include <cstring>


namespace scg
{

grammar_hook *ghook = 0;
void set_grammar_hook(grammar_hook *ghook_) { ghook = ghook_; }


static string_builder *
make_default_string_builder(const nonterminal *nt, const production &P) {
	if (nt->isterminal()) {
		return new basic_string_builder(nt->name.substr(2));
	}
	return 0;
}


static tree_builder *
make_default_tree_builder(const nonterminal *nt, const production &P) {
	if (P.rhs.size() == 1 && P.rhs[0]->isterminal()) {
		std::string sname = P.rhs[0]->name.substr(2);
		const symbol *S = symdb_findsymbol_name(sname);
		if (!S) {
			THROW_ERROR(E_NOTFOUND, "the grammar references the unknown symbol " << sname);
		}
		std::wstring ws;
		ws.append(1, S->unicode);
		if (ghook) {
			ghook->start_tree(nt, mktermsid(S->unicode));
			ghook->end_tree();
		}
		return new terminal_builder(S->name/*value*/, ws, S->mathml, S->latex);
	}

	if (ghook) ghook->start_tree(nt, InvalidSemanticId);
	std::list<tree_builder *> builders;
	for (size_t i = 0; i < P.rhs.size(); ++i) {
		builders.push_back(new ref_tree_builder(i));
		if (ghook) ghook->add_ref_child(P.rhs[i]);
	}
	if (ghook) ghook->end_tree();
	return builders.size() == 1 ? builders.front() : new subtree_builder(builders);
}



static tree_builder *
parse_tree_builder(grammar &G, const std::string &spec, nonterminal *nt, production *P) {
	if (spec.empty()) {
		return 0;
	}

	std::string::const_iterator p = spec.begin();
	while (p != spec.end() && *p != '(') {
		do {++p;} while (std::isspace(*p));
	}

	std::string type = spec.substr(0, p-spec.begin());
	if (!type.empty()) {
		if (type == "_SYMBOL") {
			P->sid = TERMINAL_EXPR;
		}
		else if (type == "_TEMP") {
			P->sid = mktempsid();
		}
		else {
			bool canonical = false;
			if (type[type.length()-1] == '*') {
				type.erase(type.end()-1);
				canonical = true;
			}
			std::map<std::string, SemanticId>::const_iterator psid = grammar_typemap().find(type);
			if (psid == grammar_typemap().end()) {
				THROW_ERROR(E_NOTFOUND, "the semantic type " << type << " is unknown; mismatch in compiled vs. read grammar?");
			}

			P->sid = psid->second;
			if (canonical) {
				assert(!G.getcanonicalsidproduction(P->sid));
				G.setcanonicalsidproduction(P->sid, P);
			}
		}
	}

	if (ghook) ghook->start_tree(nt, P->sid);

	std::list<tree_builder *> builders;

	bool force_subchild = false;

	tree_builder *builder = 0;

	size_t curr_child = 0;

	do {++p;} while (std::isspace(*p));
	while (p != spec.end() && *p != ')') {
		if (*p == '!') {
			force_subchild = true;
			do {++p;} while (std::isspace(*p));
		}
		if (*p == '+') {
			std::string::const_iterator q = p;
			++p;
			while (*q && *q != ')') {
				++q;
			}
			std::string subspec = spec.substr(p - spec.begin(), q - p + 1);
			builder = parse_tree_builder(G, subspec, 0, P);
			if (!builder) {
				return 0;
			}
			p = q + 1;
		}
		else if (*p == '%') {
			do {++p;} while (std::isspace(*p));
			if (*p == '?') {
				if (builder) {
					builder = new concat_tree_builder(builder, new placeholder_tree_builder, P);
				}
				else {
					builder = new placeholder_tree_builder;
				}
				if (ghook) ghook->add_sid_child(PLACEHOLDER_EXPR);
				++p;
			}
			else if (*p == '#') {
				if (builder) builder = new concat_tree_builder(builder, new join_tree_builder, P);
				else builder = new join_tree_builder;
				++p;
			}
			else if (*p == '^') {
				if (builder) builder = new concat_tree_builder(builder, new pullup_tree_builder, P);
				else builder = new pullup_tree_builder;
				++p;
			}
			else {
				size_t ref = *(p++) - '0';
				if (ref == 0) {
					if (builder) {
						builder = new concat_tree_builder(builder, new blank_tree_builder(), P);
					}
					else {
						builder = new blank_tree_builder();
					}
					if (ghook) ghook->add_sid_child(BLANK_EXPR);
				}
				else if (p != spec.end() && *p == ':') {
					--ref;
					assert(ref < P->rhs.size());
					unsigned subref = *(++p) - '0' - 1;
					do {++p;} while (std::isspace(*p));
					if (builder) {
						builder = new concat_tree_builder(builder, new subref_tree_builder(ref, subref), P);
					}
					else {
						builder = new subref_tree_builder(ref, subref);
					}
					if (ghook) ghook->add_subref_child(P->rhs[ref], subref);
				}
				else {
					--ref;
					assert(ref < P->rhs.size());
					if (builder) {
						builder = new concat_tree_builder(builder, new ref_tree_builder(ref), P);
					}
					else {
						builder = new ref_tree_builder(ref);
					}
					if (ghook) ghook->add_ref_child(P->rhs[ref]);
				}
			}
		}
		else if (*p == '\'') {
			++p;
			std::string::const_iterator start = p;
			while (*p != '\'') ++p;
			do {++p;} while (std::isspace(*p));
		}
		else if (*p == ',') {
			++curr_child;
			builders.push_back(builder);
			builder = 0;
			do {++p;} while (std::isspace(*p));
		}
		else {
			std::string buf;
			std::string::const_iterator start = p;
			while (p != spec.end() && *p != ')' && *p != ',' && *p != '%' && *p != '\'' && *p != '}') {
				if (*p == '\\') {
					++p;
					if (p == spec.end()) {
						break;
					}
				}
				buf.append(1, *p);
				++p;
			}
			//std::string buf = spec.substr(start - spec.begin(), p - start);
			tree_builder *tb;
			if (type == "_SYMBOL") {
				tb = new symbol_builder(buf);
				if (ghook) {
					const symbol *S = symdb_findsymbol_name(buf);
					assert(S);
					ghook->add_sid_child(S->sid);
				}
				P->rclass |= AGGREGATE_CLASS | SYMBOL_CLASS;
			}
			else {
				tb = new terminal_builder(buf, str2wstr(buf), "", "");
				if (ghook) ghook->add_sid_child(TERMINAL_EXPR);
			}
			if (builder) {
				builder = new concat_tree_builder(builder, tb, P);
			}
			else {
				builder = tb;
			}
		}
	}

	if (builder) {
		builders.push_back(builder);
	}

	if (ghook) ghook->end_tree();

	if (!builders.empty()) {
		if (!force_subchild && builders.size() == 1) {
			return builders.front();
		}
		return new subtree_builder(builders);
	}

	return 0;
}



static string_builder *
parse_string_builder(const std::string &spec, int type) {
	string_builder *builder = 0;

	if (spec.empty()) {
		return 0;
	}

	std::string::const_iterator p = spec.begin();
	while (p != spec.end() && *p != '`') {
		if (*p == '%') {
			string_builder *ref_builder = 0;
			++p;
			size_t ref = *(p++) - '0' - 1;

			if (p != spec.end() && *p == ':') {
				size_t subref = *(++p) - '0' - 1;
				do {++p;} while (std::isspace(*p));

				ref_builder = new subref_string_builder(ref, subref, type);
			}
			else {
				ref_builder = new ref_string_builder(ref, type);
			}

			if (builder) {
				builder = new concat_string_builder(builder, ref_builder);
			}
			else {
				builder = ref_builder;
			}
		}
		else if (*p == '@') {
			std::string::const_iterator start = p + 1;
			while (p != spec.end() && *p != '(') {
				++p;
			}
			std::string cmd = spec.substr(start - spec.begin(), p - start);
			start = p + 1;
			while (p != spec.end() && *p != ')') {
				++p;
			}
			std::string arg = spec.substr(start - spec.begin(), p - start);
			++p;
			if (cmd == "join") {
				if (builder) builder = new concat_string_builder(builder, new join_string_builder(arg, type));
				else builder = new join_string_builder(arg, type);
			}
			else {
				return 0;
			}
		}
		else {
			std::string::const_iterator start = p;
			while (p != spec.end() && *p != '`' && *p != '%' && *p != '@') {
				++p;
			}
			std::string buf = spec.substr(start - spec.begin(), p - start);
			if (builder) {
				builder = new concat_string_builder(builder, new basic_string_builder(buf)); 
			}
			else {
				builder = new basic_string_builder(buf);
			}
		}
	}

	return builder;
}

static int
readrclass(const std::string &s) {
	int rclass = 0;
	rclass_t rc;
	std::string::size_type i = 3;
	std::string::size_type j = s.find_first_of(',', i);
	while (j != std::string::npos) {
		rc = tag_to_rclass(s.substr(i, j-i));
		if (rc != -1) {
			rclass = rclass | (1 << rc);
		}
		i = j+1;
		j = s.find_first_of(',', i);
	}
	rc = tag_to_rclass(s.substr(i, s.length()-1-i));
	if (rc != -1) {
		rclass = rclass | (1 << rc);
	}
	return rclass;
}

static void
parse_nonterminal(std::istream &is, grammar &G, nonterminal *nt, std::map<std::string, nonterminal *> &nt_table)
{
	VERBOSE2(*verb_out << "grammar: parsing nonterminal " << nt->name << std::endl);
	int prod_attach = attach_mode::UNSPECIFIED;
	std::string s;
	is >> s;
	if (s == ":g:") {
		prod_attach = attach_mode::GROUP;
	}
	else if (s == ":s:") {
		prod_attach = attach_mode::SYMBOL;
	}
	else if (s == ":h:") {
		prod_attach = attach_mode::SYMBOL_HEAD;
	}
	else if (s == ":t:") {
		prod_attach = attach_mode::SYMBOL_TAIL;
	}
	else if (s == "::") {
		prod_attach = attach_mode::UNSPECIFIED;
	}
	else {
		THROW_ERROR(E_INVALID, "production name must be followed by \"::\" (read " << s << ")\n");
	}

	production *P = new production(nt, nt->productions.size());
	P->attach_in_parent = prod_attach;

	for (;;) {
		is >> s;
		if (!is || is.eof()) {
			THROW_ERROR(E_IO, "input error while reading nonterminal " << nt->name);
		}

		if (s == ";") {
			if (P->tail == (size_t)-1) {
				P->tail = P->rhs.size() - 1;
			}
			if (P->sid == InvalidSemanticId && P->rhs.size() == 1 && P->rhs[0]->isterminal()) {
				P->sid = TERMINAL_EXPR;
			}
			nt->addproduction(P);
			break;
		}
		else if (s == "op") {
			P->isop = true;
		}
		else if (s.length() > 4 && s[0] == 'r' && s[1] == 'c' && s[2] == '(' && s[s.length()-1] == ')') {
			P->rclass |= readrclass(s);
		}
		else if (s.length() > 6 && !std::strncmp(s.c_str(), "bias(", 5) && s[s.length()-1] == ')') {
			std::stringstream ss;
			ss << s.substr(5, s.length()-6);
			ss >> P->bias;
		}
		else if (s == "|") {
			if (P->tail == (size_t)-1) {
				P->tail = P->rhs.size() - 1;
			}
			if (P->rhs.size() == 1 && P->rhs[0]->isterminal()) {
				P->sid = TERMINAL_EXPR;
			}
			nt->addproduction(P);
			P = new production(nt, nt->productions.size());
			P->attach_in_parent = prod_attach; 
		}
		else {
			size_t n = s.length();
			if (s[0] == '^') {
				s = s.substr(1);
				n--;
				P->head = P->rhs.size();
			}
			
			if (s[s.length()-1] == '$') {
				s = s.substr(0, s.length() - 1);
				P->tail = P->rhs.size();
			}

			if (s[0] == '<') {
				std::stringstream ss;
				ss << s;
				while (s[s.length()-1] != '>') {
					is >> s;
					if (!is || is.eof()) {
						THROW_ERROR(E_INVALID, "EOF encountered while reading relation");
					}
					ss << s;
				}
				s = ss.str();
				n = s.length();

				size_t split = s.find_first_of(':');
				if (split != std::string::npos) {
					std::string type = s.substr(split + 1, 2);
					s = s.substr(0, split);

					if (type.length() != 2) {
						THROW_ERROR(E_INVALID, "attachment mode specifier must be only two characters; found " << type);
					}
					else {
						attach_mode &mode = P->attach_modes.back();
						int *modep = &mode.from;
						for (int i = 0; i < 2; ++i) {
							if (std::tolower(type[i]) == 'g') {
								*modep = attach_mode::GROUP;
							}
							else if (std::tolower(type[i]) == 's') {
								*modep = attach_mode::SYMBOL;
							}
							else {
								THROW_ERROR(E_INVALID, "attachment mode must be g or s, not " << type[0]);
							}
							
							modep = &mode.to;
						}
					}
				}
				else {
					s[n - 1] = 0;
				}
				
				const char *link_type = s.c_str() + 1;
				if (!std::strcmp(link_type, "B")) {
					P->rel = get_relation(REL_BELOW);
				}
				else if (!std::strcmp(link_type, "R")) {
					P->rel = get_relation(REL_RIGHT);
				}
				else if (!std::strcmp(link_type, "AR")) {
					P->rel = get_relation(REL_ABOVERIGHT);
				}
				else if (!std::strcmp(link_type, "BR")) {
					P->rel = get_relation(REL_BELOWRIGHT);
				}
				else if (!std::strcmp(link_type, "C")) {
					P->rel = get_relation(REL_CONTAINS);
				}
				P->dim = ordering_for_relation(P->rel);
			}
			else if (s[0] == '[') {
				std::stringstream ss;
				ss << s;
				while (s[s.length()-1] != ']') {
					is >> s;
					if (!is || is.eof()) {
						THROW_ERROR(E_INVALID, "error: EOF encountered while reading nonterminal");
					}
					ss << s;
				}
				s = ss.str();
				n = s.length();

				const char *nt_label = s.c_str() + 1;
				s[n - 1] = 0;

				nonterminal *&rhs_nt = nt_table[nt_label];
				if (!rhs_nt) {
					rhs_nt = new nonterminal(nt_label, G.nts.size());
					G.add_nonterminal(rhs_nt);
				}
				P->rhs.push_back(rhs_nt);
				P->attach_modes.push_back(attach_mode());
			}
			else if (s[0] == '{') {
				std::stringstream ss;
				ss << s;
				while (s[s.length()-1] != '}') {
					is >> s;
					if (!is || is.eof()) {
						THROW_ERROR(E_INVALID, "error: EOF encountered while reading nonterminal");
					}
					ss << s;
				}
				s = ss.str();
				n = s.length();

				P->tbuild = parse_tree_builder(G, s.substr(1, n - 2), nt, P);
			}
			else if (s[0] == '`') {
				std::string spec;
				if (s[n-1] == '`' && n > 1) {
					spec = std::string(&s[0] + 1, n-2);
				}
				else {
					spec = std::string(&s[0] + 1, n-1);
					do {
						is >> s;
						n = s.length();
						spec += " ";
						if (s[n-1] == '`') {
							spec += std::string(&s[0], n-1);
						}
						else {
							spec += std::string(&s[0], n);
						}
					} while (s[n-1] != '`');
				}
				if (!P->sbuild) {
					P->sbuild = parse_string_builder(spec, strid::mathml);
				}
				else {
					P->lbuild = parse_string_builder(spec, strid::latex);
				}
			}
			else {
				P->terminal_indices.push_back(P->rhs.size());
				std::string term("_S");
				term.append(s);
				nonterminal *&rhs_nt = nt_table[term];
				if (!rhs_nt) {
					rhs_nt = new nonterminal(term, G.nts.size());
					symbol *S = symdb_findsymbol_name(s);
					if (!S) {
						THROW_ERROR(E_NOTFOUND, "the grammar references the unknown symbol " << s);
					}
					/*
					production *termP = new production(rhs_nt, rhs_nt->productions.size());
					termP->sid = TERMINAL_EXPR;
					S->info.P = termP;
					rhs_nt->productions.push_back(termP);
					*/
					S->nt = rhs_nt;
					//rhs_nt->rclasses.insert(S->info.classes.begin(), S->info.classes.end());
					rhs_nt->rclass = rhs_nt->rclass | S->rclass;
					G.add_nonterminal(rhs_nt);
				}
				P->rhs.push_back(rhs_nt);
				P->attach_modes.push_back(attach_mode());
			}
		}
	}
}


std::istream &
operator>>(std::istream &is, grammar &g) {
	read_grammar(is, g);
	rebuild_production_lengths(g);
	return is;
}

static void
compute_equivalent_nts(nonterminal *root) {
	if (!root->equivalent_nts.empty()) {
		return;
	}

	root->equivalent_nts.insert(root);
	for (std::vector<production *>::iterator pp = root->productions.begin(); pp != root->productions.end(); ++pp) {
		production &P = **pp;
		if (P.rhs.size() == 1) {
			nonterminal *child = P.rhs[0];
			compute_equivalent_nts(child);
			root->equivalent_nts.insert(child->equivalent_nts.begin(), child->equivalent_nts.end());
		}
	}
	
	VERBOSE2(
		*verb_out << "equivalent NTS for " << root->name << " : ";
		for (std::set<const nonterminal *>::const_iterator i = root->equivalent_nts.begin(); i != root->equivalent_nts.end(); ++i) {
			*verb_out << (*i)->name << ", ";
		}
		*verb_out << std::endl;
	);
}

static void
compute_rootdists(nonterminal *root, unsigned d) {
	if (root->rootdist != ~0) return;

	root->rootdist = d;
	for (std::vector<production *>::iterator pp = root->productions.begin(); pp != root->productions.end(); ++pp) {
		production *P = *pp;
		for (std::vector<nonterminal *>::iterator i = P->rhs.begin(); i != P->rhs.end(); ++i) {
			compute_rootdists(*i, d+1);
		}
	}
}

static std::pair<size_t, size_t>
compute_lengths(nonterminal *root) {
	if (root->min_length != std::numeric_limits<size_t>::max() && root->max_length != 0) {
		return std::make_pair(root->min_length, root->max_length);
	}
	if (root->min_length == 0 && root->max_length == 0) {
		// recursion in grammar
		return std::make_pair(1, std::numeric_limits<size_t>::max());
	}

	/*for (std::list<const ordered_segments *>::const_iterator i = root->terminal_productions.begin(); i != root->terminal_productions.end(); ++i) {
		const ordered_segments *osegs = *i;
		root->min_length = std::min(root->min_length, osegs->size());
		root->max_length = std::max(root->max_length, osegs->size());
	}*/

	if (root->name[0] == '_' && root->name[1] == 'S') {
		std::string sname = root->name.substr(2);
		symbol *S = symdb_findsymbol_name(sname);
		if (!S) {
			THROW_ERROR(E_NOTFOUND, "the grammar references the unknown symbol " << sname);
		}
		for (prototype *p = S->firstproto(); p; p = S->nextproto(p)) {
			root->min_length = std::min(root->min_length, p->strokes.nstrokes);
			root->max_length = std::max(root->max_length, p->strokes.nstrokes);
		}
	}
	else if (root->name[0] == '!') {
		root->min_length = 1;
		root->max_length = std::numeric_limits<size_t>::max();
	}
	else {
		root->min_length = 0;
		root->max_length = 0;
		size_t minlen = std::numeric_limits<size_t>::max();
		size_t maxlen = 0;
		for (std::vector<production *>::iterator pp = root->productions.begin(); pp != root->productions.end(); ++pp) {
			production &P = **pp;
			P.min_prefix_length.resize(P.rhs.size());
			P.max_prefix_length.resize(P.rhs.size());
			std::vector<size_t>::iterator pminlen = P.min_prefix_length.begin();
			std::vector<size_t>::iterator pmaxlen = P.max_prefix_length.begin();
			size_t minacc = 0;
			size_t maxacc = 0;
			for (std::vector<nonterminal *>::const_iterator i = P.rhs.begin(); i != P.rhs.end(); ++i) {
				std::pair<size_t, size_t> lens = compute_lengths(*i);
				minacc += lens.first;
				if (maxacc != std::numeric_limits<size_t>::max()) {
					if (lens.second == std::numeric_limits<size_t>::max()) {
						maxacc = lens.second;
					}
					else {
						maxacc += lens.second;
					}
				}
				*pminlen = minacc;
				*pmaxlen = maxacc;
				++pminlen;
				++pmaxlen;
			}
			minlen = std::min(minlen, P.min_prefix_length.back());
			maxlen = std::max(maxlen, P.max_prefix_length.back());
		}
		root->min_length = minlen;
		root->max_length = maxlen;
	}
	
	if (root->min_length == std::numeric_limits<size_t>::max()) {
		root->min_length = 1;
	}
	if (root->max_length == 0) {
		root->max_length = std::numeric_limits<size_t>::max();
	}
	VERBOSE2(*verb_out << "length range for " << root->name << " is (" << root->min_length << "," << root->max_length << ")\n");

	return std::make_pair(root->min_length, root->max_length);
}


void
rebuild_production_lengths(grammar &G) {
	for (std::vector<nonterminal *>::iterator i = G.nts.begin(); i != G.nts.end(); ++i) {
		for (std::vector<production *>::iterator j = (*i)->productions.end(); j != (*i)->productions.end(); ++j) {
			(*j)->min_prefix_length.clear();
			(*j)->max_prefix_length.clear();
		}
		(*i)->min_length = std::numeric_limits<size_t>::max();
		(*i)->max_length = 0;
		
		(*i)->equivalent_nts.clear();
	}

	for (std::vector<nonterminal *>::iterator i = G.nts.begin(); i != G.nts.end(); ++i) {
		compute_lengths(*i);
		compute_equivalent_nts(*i);
		/*if ((*i)->max_length == 0) {
			(*i)->max_length = std::numeric_limits<size_t>::max();
		}*/
	}
	if (G.root) compute_rootdists((nonterminal *)G.root, 0);
	if (G.single_expr_root) compute_rootdists((nonterminal *)G.single_expr_root, 0);
	if (G.multi_expr_root) compute_rootdists((nonterminal *)G.multi_expr_root, 0);
	if (G.compexpr_root) compute_rootdists((nonterminal *)G.compexpr_root, 0);
	if (G.partexpr_root) compute_rootdists((nonterminal *)G.partexpr_root, 0);
	if (G.matrix_cell_root) compute_rootdists((nonterminal *)G.matrix_cell_root, 0);
	/*for (std::vector<nonterminal *>::iterator i = G.nts.begin(); i != G.nts.end(); ++i) {
		if ((*i)->equivalent_nts.empty()) {
			compute_equivalent_nts(*i);
		}
	}*/
}

static int getrclass(nonterminal *nt);

static int
getrclass(production *P) {
	if (P->rclass == 0) {
		P->rclass = (1 << AGGREGATE_CLASS); // prevent infinite recursion
		if (P->rhs.size() == 1) return getrclass(P->rhs[0]);
		assert(P->rel);
		if (P->rel->name != "R") {
			P->rclass |= (1 << BOX_CLASS);
		}
		else {
			for (std::vector<nonterminal *>::iterator i = P->rhs.begin(); i != P->rhs.end(); ++i) {
				P->rclass |= getrclass(*i);
			}
		}
		VERBOSE2(
			*verb_out << "rclass for production " << P->nt->name << " :: ";
			for (std::vector<nonterminal *>::iterator i = P->rhs.begin(); i != P->rhs.end(); ++i) {
				*verb_out << (*i)->name << ' ';
			}
			*verb_out << "is ";
			for (int i = 0; i < NCLASSES-1; ++i) {
				if (P->rclass & (1 << i)) {
					*verb_out << rclass_to_str(i) << ' ';
				}
			}
			*verb_out << std::endl;
		);
	}
	return P->rclass;
}

static int
getrclass(nonterminal *nt) {
	if (nt->rclass == 0) {
		nt->rclass = (1 << AGGREGATE_CLASS); // prevent infinite recursion
		for (std::vector<production *>::iterator i = nt->productions.begin(); i != nt->productions.end(); ++i) {
			production *P = *i;
			nt->rclass |= getrclass(P);
		}
		VERBOSE2(
			*verb_out << "rclass for nonterminal " << nt->name << " is ";
			for (int i = 0; i < NCLASSES; ++i) {
				if (nt->rclass & (1 << i)) {
					*verb_out << rclass_to_str(i) << ' ';
				}
			}
			*verb_out << std::endl;
		);
	}
	return nt->rclass;
}

void
read_grammar(std::istream &is, grammar &G) {
	static std::string training_path;

	std::string start_state;
	std::string partstart_state;
	std::string compstart_state;
	std::string singlestart_state;
	std::string multistart_state;
	std::string matrix_cell_state;
	std::map<std::string, nonterminal *> nt_table;

	//curr_group_directions = &g.group_directions;
	for (;;) {
		std::string label;
		is >> label;

		if (is.eof()) {
			break;
		}

		while (label[0] == '#') {
			char buf[256];
			do {
				buf[255] = 0;
				is.getline(buf, sizeof(buf));
			} while (buf[255]);

			is >> label;
			if (is.eof()) {
				goto done;
			}
			if (!is) {
				THROW_ERROR(E_IO, "input error while reading comment; last line " << label);
			}
		}

		if (label == "@root") {
			is >> label;
			if (!is || is.eof()) {
				THROW_ERROR(E_IO, "nonterminal name missing in @root statement");
			}
			start_state = label;
		}
		else if (label == "@partial") {
			is >> label;
			if (!is || is.eof()) {
				THROW_ERROR(E_IO, "nonterminal name missing in @partial statement");
			}
			partstart_state = label;
		}
		else if (label == "@complete") {
			is >> label;
			if (!is || is.eof()) {
				THROW_ERROR(E_IO, "nonterminal name missing in @complete statement");
			}
			compstart_state = label;
		}
		else if (label == "@single") {
			is >> label;
			if (!is || is.eof()) {
				THROW_ERROR(E_IO, "nonterminal name missing in @single statement");
			}
			singlestart_state = label;
		}
		else if (label == "@multi") {
			is >> label;
			if (!is || is.eof()) {
				THROW_ERROR(E_IO, "nonterminal name missing in @multi statement");
			}
			multistart_state = label;
		}
		else if (label == "@matrix_cell") {
			is >> label;
			if (!is || is.eof()) {
				THROW_ERROR(E_IO, "nonterminal name missing in @matrix_cell statement");
			}
			matrix_cell_state = label;
		}
		else {
			nonterminal *&nt = nt_table[label];
			if (!nt) {
				nt = new nonterminal(label, G.nts.size());
				G.add_nonterminal(nt);
			}
			parse_nonterminal(is, G, nt, nt_table);
		}
	}

done:
	G.root = nt_table[start_state];
	if (!G.root) {
		WARNING(E_NOTFOUND, "the start state, named " << start_state << " has no definition");
	}
	G.partexpr_root = nt_table[partstart_state];
	if (!G.partexpr_root) {
		WARNING(E_NOTFOUND, "the partial start state, named " << partstart_state << " has no definition");
	}
	G.compexpr_root = nt_table[compstart_state];
	if (!G.compexpr_root) {
		WARNING(E_NOTFOUND, "the complete start state, named " << compstart_state << " has no definition");
	}
	G.single_expr_root = nt_table[singlestart_state];
	if (!G.single_expr_root) {
		WARNING(E_NOTFOUND, "the single-expr start state, named " << singlestart_state << " has no definition");
	}
	G.multi_expr_root = nt_table[multistart_state];
	if (!G.multi_expr_root) {
		WARNING(E_NOTFOUND, "the multi-expr start state, name " << multistart_state << " has no definition");
	}
	G.matrix_cell_root = nt_table[matrix_cell_state];
	if (!G.matrix_cell_root) {
		WARNING(E_NOTFOUND, "the matrix-cell root state, name " << matrix_cell_state << " has no definition");
	}

	for (std::vector<nonterminal *>::iterator i = G.nts.begin(); i != G.nts.end(); ++i) {
		nonterminal *nt = *i;
		for (std::vector<production *>::iterator j = nt->productions.begin(); j != nt->productions.end(); ++j) {
			production &P = **j;
			if (!P.sbuild) {
				P.sbuild = make_default_string_builder(nt, P);
			}
			if (!P.lbuild) {
				P.lbuild = make_default_string_builder(nt, P);
			}
			if (!P.tbuild) {
				P.tbuild = make_default_tree_builder(nt, P);
			}
		}
	}

	// XXX: hack to let the parser know that matrices derive expressions
	// (extra-gramatically)
	if (nt_table["E"]) {
		if (nt_table["!MATRIX"]) {
			production P(nt_table["!MATRIX"], 0);
			P.rclass = 1 << BOX_CLASS;
			P.rhs.push_back(nt_table["E"]);
			nt_table["!MATRIX"]->addproduction(&P);
			rebuild_production_lengths(G);
			getrclass((nonterminal *)G.root);
			nt_table["!MATRIX"]->productions.clear();
		}
		if (nt_table["!MULTI"]) {
			production P(nt_table["!MULTI"], 0);
			P.rclass = 1 << BOX_CLASS;
			P.rhs.push_back(nt_table["E"]);
			nt_table["!MULTI"]->addproduction(&P);
			rebuild_production_lengths(G);
			getrclass((nonterminal *)G.root);
			nt_table["!MULTI"]->productions.clear();
		}
	}
	else {
		rebuild_production_lengths(G);
		getrclass((nonterminal *)G.root);
	}
}

nonterminal::nonterminal(const std::string &name_, size_t index_)
	: name(name_), index(index_), min_length(std::numeric_limits<size_t>::max()), max_length(0), rclass(0), rootdist(~0) {
}

nonterminal::~nonterminal() {
	for (std::vector<production *>::iterator i = productions.begin(); i != productions.end(); ++i) {
		delete *i;
	}
}

void
nonterminal::addproduction(production *P) {
	productions.push_back(P);
	for (std::vector<nonterminal *>::iterator i = P->rhs.begin(); i != P->rhs.end(); ++i) {
		nonterminal *child = *i;
		child->derived_by_.insert(this);
		child->derived_by_.insert(derived_by_.begin(), derived_by_.end());
		for (std::set<const nonterminal *>::iterator j = child->derives_.begin(); j != child->derives_.end(); ++j) {
			nonterminal *nt = const_cast<nonterminal *>(*j);
			nt->derived_by_.insert(this);
			nt->derived_by_.insert(derived_by_.begin(), derived_by_.end());
		}
		derives_.insert(child);
		derives_.insert(child->derives_.begin(), child->derives_.end());
		for (std::set<const nonterminal *>::iterator j = derived_by_.begin(); j != derived_by_.end(); ++j) {
			nonterminal *nt = const_cast<nonterminal *>(*j);
			nt->derives_.insert(child);
			nt->derives_.insert(child->derives_.begin(), child->derives_.end());
		}
	}
}

grammar::~grammar() {
	for (std::vector<nonterminal *>::iterator i = nts.begin(); i != nts.end(); ++i) {
		delete *i;
	}
}


void
grammar::add_nonterminal(nonterminal *nt) {
	nts.push_back(nt);
}


production::~production() {
	delete sbuild;
	delete lbuild;
	delete tbuild;
}

const production *
grammar::getcanonicalsidproduction(SemanticId sid) const {
	std::map<SemanticId, const production *>::const_iterator i = sidmap.find(sid);
	return (i == sidmap.end()) ? 0 : i->second;
}

const production *
grammar::production_for_tree(SemanticId sid, const std::vector<basic_tree *> &children) const {
	if (sid == PRIMESYMBOL_EXPR) {
		const production *P = getcanonicalsidproduction(sid);
		if (!P) return 0;
		const nonterminal *nt = P->nt;
		if (!nt) return 0;
		if (children.size() != 1) return 0;
		long nprimes = strtol(children[0]->str(), 0, 10);
		if (nprimes <= 0 || nprimes >= nt->productions.size()) return 0;
		return nt->productions[nprimes-1];
	}
	else {
		return getcanonicalsidproduction(sid);
	}
}

void
grammar::setcanonicalsidproduction(SemanticId sid, const production *P) {
	sidmap[sid] = P;
}

static std::map<int, const production *> &
uniqprods() {
	static std::map<int, const production *> ups;
	return ups;
}

uniqproduction::uniqproduction(int id) : production(0, id) {
	if (uniqprods().find(id) != uniqprods().end()) {
		THROW_ERROR(E_ALREADY, "uniqproduction with id " << id << " already registered");
	}
	uniqprods()[id] = this;
}

const production *
getuniqproduction(int id) {
	return id ? uniqprods()[id] : 0;
}

void
clearuniqproductions() {
	uniqprods().clear();
}

}
