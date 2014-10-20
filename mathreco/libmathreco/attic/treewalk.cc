#include "treewalk.h"
#include "expr-node.h"
#include "grammar.h"
#include "parser-priv.h"
//#include "parser.h"
#include "fuzzy.h"
#include "mathrecognizer-private.h"
#include "profman.h"
#include "links.h"
#include "parms.h"
#include "ref.h"
#include "verb.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <errno.h>


namespace scg
{


unsigned NUM_TERMINAL_PREVIEWS = GetParameterUnsigned("NumTerminalPreviews");
unsigned NUM_NONTERMINAL_PREVIEWS = GetParameterUnsigned("NumNonTerminalPreviews");


expression_node *
make_terminal(const std::string &term, bool want_unicode)
{
	return make_terminal(term, CNFGrammar::INVALID_GID);
}


expression_node *
make_terminal(const std::string &term, int gid, bool want_unicode)
{
	VERBOSE(*verb_out << "MAKE TERM " << term << " (" << gid << ")\n");
#ifndef NO_WSTRING
	if (want_unicode && global_symbols_db()) {
		const symbol *S = global_symbols_db()->find_symbol(term);
		if (S) {
		    return make_terminal(term, std::wstring(1, S->info.unicode), gid);
		}
	}

	wchar_t *wcs;
	size_t mbsz = mbstowcs(0, term.c_str(), 0);
	wcs = DEBUG_NEW wchar_t[mbsz + 1];
	if (mbstowcs(wcs, term.c_str(), mbsz+1) != mbsz) {
		throw E_INVALID;
	}

	expression_node *expr = make_terminal(term, std::wstring(wcs), gid);
	delete[] wcs;
#else
	expression_node *expr = DEBUG_NEW expression_node;
	expr->set_type(TERMINAL_EXPR);
	expr->set_grammar_type(gid);
	expr->set_long_string(term);
	expr->set_str(term);
#endif
	return expr;
}


#ifndef NO_WSTRING
expression_node *
make_terminal(const std::string &term, const std::wstring &wterm, int gid)
{
	expression_node *node = DEBUG_NEW expression_node;
	node->set_type(TERMINAL_EXPR);
	node->set_grammar_type(gid);
	node->set_long_string(term);
	node->set_str(term);
	node->set_wstr(wterm);
	return node;
}
#endif


static expression_box
build_box(const BboxCandidate * const *boxes, const std::vector<bool> &selection)
{
	expression_box box;

	std::vector<bool>::const_iterator i;
	const BboxCandidate * const *b = boxes;
	for (i = selection.begin(); i != selection.end(); ++i, ++b) {
		if (*i) {
			box = (*b)->box;
			break;
		}
	}

	for (; i != selection.end(); ++i, ++b) {
		if (*i) {
			box = merge(box.box, (*b)->box.box);
		}
	}

	return box;
}


static expression_node *
construct_expression(const CNFProduction *P, expression_node *left, expression_node *right)
{
	SemanticId default_type;
	if (CNFGrammar::is_terminal(P->first) && !CNFGrammar::is_nonterminal(P->second)) {
		default_type = TERMINAL_EXPR;
	}
	else if (CNFGrammar::is_nonterminal(P->first) && P->second == CNFGrammar::INVALID_GID) {
		default_type = left->type();
	}
	else {
		default_type = InvalidSemanticId;
	}

	expression_node *expr;

	if (P->tree_builder) {
		expr = P->tree_builder->build(left, right);

		if (CNFGrammar::is_valid_sid(P->type)) {
			expr->set_type(P->type);
		}
		else {
			expr->set_type(default_type);
		}
	}
	else {
		expr = DEBUG_NEW expression_node;
		expr->set_type(default_type);
	}

	if (P->string_builder) {
		expr->set_long_string(P->string_builder->build(left, right));
	}

	return expr;
}


expression_iterator::expression_iterator(const ParseContext *ctx_, std::set<CNFGrammar::gid_t> gids_, const unsigned *strokes_, unsigned nstrokes_, int flags)
	: ExpressionIterator(), best_tree(0)
{
	SetFlags(flags);

	std::vector<bool> strokes(ctx_->boxes.size(), false);

	for (const unsigned *s = strokes_; s != strokes_ + nstrokes_; ++s) {
		strokes[*s] = true;
	}

	for (unsigned d = 0; d < ctx_->orderings.size(); ++d) {
		const DimData &dim = ctx_->impl->dim(d);
		unsigned ibmin = ctx_->boxes.size() - 1;
		unsigned ibmax = 0;

		for (const unsigned *s = strokes_; s != strokes_ + nstrokes_; ++s) {
			unsigned dimi = dim.parse_to_dim_index[*s];
			if (dimi < ibmin) {
				ibmin = dimi;
			}
			if (dimi > ibmax) {
				ibmax = dimi;
			}
		}

		for (std::set<CNFGrammar::gid_t>::const_iterator i = gids_.begin(); i != gids_.end(); ++i) {
			iterators.push_back(std::make_pair(DEBUG_NEW tree_iterator(shared, ctx_, strokes), 0));
			CellRef ref(*i, ibmin, ibmax, d);
			shared.cache[cell_key(ref, strokes)] = iterators.back().first;
			iterators.back().first->build(ref);
		}
	}
}


expression_iterator::~expression_iterator()
{
	for (iter_cache_t::iterator i = shared.cache.begin(); i != shared.cache.end(); ++i) {
		delete i->second;
	}
}

void
tree_iterator::extract_source_terminals(const CellSource &src, std::vector<bool> &terminals)
{
	const RangeProjectionData &rd = ctx->impl->dim(src.dim).range(src.start, src.end).projection(self_key.root.dim);
	
	std::vector<unsigned>::const_iterator first = rd.boxes.begin() + src.auxstart;
	std::vector<unsigned>::const_iterator last = rd.boxes.begin() + src.auxend + 1;
	while (first != last) {
		terminals[*first] = true;
		++first;
	}
}

void
tree_iterator::extract_source_left_terminals(const CellSource &src, std::vector<bool> &terminals)
{
	const RangeProjectionData &rd = ctx->impl->dim(src.dim).range(src.start, src.end).projection(self_key.root.dim);
	
	std::vector<unsigned>::const_iterator first = rd.boxes.begin() + src.auxstart;
	std::vector<unsigned>::const_iterator last = rd.boxes.begin() + src.split + 1;
	while (first != last) {
		terminals[*first] = true;
		++first;
	}
}

void
tree_iterator::extract_source_right_terminals(const CellSource &src, std::vector<bool> &terminals)
{
	const RangeProjectionData &rd = ctx->impl->dim(src.dim).range(src.start, src.end).projection(self_key.root.dim);
	
	std::vector<unsigned>::const_iterator first = rd.boxes.begin() + src.split + 1;
	std::vector<unsigned>::const_iterator last = rd.boxes.begin() + src.auxend + 1;
	while (first != last) {
		terminals[*first] = true;
		++first;
	}
}


tree_iterator::tree_iterator(iterator_shared_data &shared_, const ParseContext *ctx_, const std::vector<bool> &terminals_) : ctx(ctx_), shared(shared_), is_terminal(false), is_distinct_subtree(false), cell(0), last_clean(0)
{
	bounds = build_box(&ctx->boxes[0], terminals_);
	self_key.terminals = terminals_;
}

tree_iterator::~tree_iterator()
{
	for (std::vector<expression_node *>::iterator i = ranked_trees.begin(); i != ranked_trees.end(); ++i) {
		delete *i;
	}
}


tree_iterator *
tree_iterator::build_subiter(const CellRef &child, const std::vector<bool> &terminals, unsigned nbranches)
{
	if (child.gid != CNFGrammar::INVALID_GID) {
		cell_key key(child, terminals);
		iter_cache_t::iterator i = shared.cache.lower_bound(key);
		if (i != shared.cache.end() && i->first == key) {
			return i->second;
		}

		tree_iterator *iter = DEBUG_NEW tree_iterator(shared, ctx, terminals);
		shared.cache.insert(i, std::make_pair(key, iter));

		iter->build_(child, nbranches);
		return iter;
	}

	return 0;
}


const cell_key cell_key::NIL;


std::ostream &
operator<<(std::ostream &os, const cell_key &key)
{
	os << key.root << " (";
	if (key.terminals.empty()) {
		os << "none";
	}
	else {
		os << key.terminals;
	}
	return os << ")";
}


void
tree_iterator::build(const CellRef &root)
{
	build_(root, false);
}


static bool
terminals_dirty(const ParseContext *ctx, const terminal_mask_t &terminals, unsigned first_clean, unsigned last_clean)
{
	for (std::vector<terminal_mask_t>::const_iterator pdirty = ctx->dirty.begin() + first_clean; pdirty != ctx->dirty.begin() + last_clean; ++pdirty) {
		const terminal_mask_t &dirty = *pdirty;
		VERBOSE(*verb_out << " comparing dirty mask " << dirty << " vs. " << terminals << std::endl);
		assert(dirty.size() == terminals.size());
		terminal_mask_t::const_iterator d = dirty.begin();
		terminal_mask_t::const_iterator i = terminals.begin();
		for ( ; d != dirty.end(); ++d, ++i) {
			if (*d && !*i) {
				break;
			}
		}
		if (d == dirty.end()) {
			return true;
		}
	}
	return false;
}


static bool
partial_intersection(const ParseContext *ctx, const terminal_mask_t &terminals, unsigned first_clean, unsigned last_clean)
{
	for (std::vector<terminal_mask_t>::const_iterator pdirty = ctx->dirty.begin() + first_clean; pdirty != ctx->dirty.begin() + last_clean; ++pdirty) {
		const terminal_mask_t &dirty = *pdirty;
		VERBOSE(*verb_out << " comparing dirty mask " << dirty << " vs. " << terminals << std::endl);
		assert(dirty.size() == terminals.size());
		terminal_mask_t::const_iterator d = dirty.begin();
		terminal_mask_t::const_iterator i = terminals.begin();
		bool intersection = false;
		bool complete = true;
		for ( ; d != dirty.end(); ++d, ++i) {
			if (*d) {
				if (*i) {
					intersection = true;
				}
				else {
					complete = false;
				}
			}
		}
		if (intersection && !complete) {
			return true;
		}
	}
	return false;
}



static std::pair<cell_key, cell_key>
make_key_pair(const tree_iterator *left, const tree_iterator *right)
{
	std::pair<cell_key, cell_key> keys;

	keys.first = left->key();
	keys.second = right ? right->key() : cell_key::NIL;

	return keys;
}


bool
tree_iterator::next_indices_exist(const tree_iterator *left, const tree_iterator *right)
{
	std::pair<cell_key, cell_key> keys = make_key_pair(left, right);
	next_index_cache_t::iterator i = next_indices.find(keys);
	return i != next_indices.end();
}


tree_iterator::next_indices_t &
tree_iterator::get_next_indices(const tree_iterator *left, const tree_iterator *right)
{
	std::pair<cell_key, cell_key> keys = make_key_pair(left, right);
	next_index_cache_t::iterator i = next_indices.lower_bound(keys);
	if (i == next_indices.end() || i->first != keys) {
		i = next_indices.insert(i, std::make_pair(keys, next_indices_t()));
		i->second.push_back(0);
	}
	return i->second;
}


void
tree_iterator::build_(const CellRef &root, unsigned nbranches)
{
	VERBOSE(*verb_out << "tree_iterator::build() at root " << root << ", terminals " << self_key.terminals << std::endl);

	self_key.root = root;

	is_distinct_subtree = nbranches > 1;
	VERBOSE(*verb_out << " marking distinct=" << is_distinct_subtree << std::endl);

	if (CNFGrammar::is_terminal(root.gid)) {
		const CNFTerminal *term = ctx->g->terminal(root.gid);
		expression_node *tree = make_terminal(term->name, root.gid);
		double score = ctx->impl->terminal_score(root.gid, ctx->impl->dim(root.dim).dim_to_parse_index[root.start]);
		tree->set_score(score);
		tree->set_strokes(self_key.terminals);
		tree->set_box(bounds);
		tree->set_cellref(root);
		tree->set_grammar_type(root.gid);
		VERBOSE(*verb_out << " terminal tree is " << tree->str() << std::endl);
		ranked_trees.push_back(tree);

		VERBOSE(*verb_out << " marking is_terminal=1 at " << self_key << std::endl);
		is_terminal = true;

		return;
	}

	const ParseCell *cell = parse_cell(ctx, root);
	if (!cell) return;

	this->cell = cell;

	bool allow_only_transient_productions = false;
	
	std::map<terminal_mask_t, lock_descriptor>::const_iterator pldesc = ctx->locks.find(self_key.terminals);
	if (pldesc != ctx->locks.end()) {
		const lock_descriptor &ldesc = pldesc->second;
		VERBOSE(*verb_out << " build(): this node is locked; ");
		VERBOSE(
			*verb_out << " valid gids are: ";
			for (std::vector<CNFGrammar::gid_t>::const_iterator i = ldesc.valid_gids.begin(); i != ldesc.valid_gids.end(); ++i) {
				*verb_out << *i << ' ';
			}
			*verb_out << std::endl;
		);
		if (std::find(ldesc.valid_gids.begin(), ldesc.valid_gids.end(), self_key.root.gid) == ldesc.valid_gids.end()) {
			VERBOSE(*verb_out << "gid " << self_key.root.gid << " invalid here; returning null\n");
			allow_only_transient_productions = true;//return;
		}
		VERBOSE(*verb_out << "copying " << ldesc.expr->str() << std::endl);
		ranked_trees.push_back(DEBUG_NEW expression_node(*ldesc.expr));
		return;
	}

	for (std::vector<ParseRef>::const_iterator pref = cell->refs.begin(); pref != cell->refs.end(); ++pref) {
		VERBOSE(*verb_out << " parseref " << *pref << std::endl);
		if (allow_only_transient_productions
		 && (!CNFGrammar::is_nonterminal(pref->spec.prod->first)
		  || CNFGrammar::is_valid_gid(pref->spec.prod->second))) {
			continue;
		}
		for (std::set<CellSource>::const_iterator psrc = pref->sources.begin(); psrc != pref->sources.end(); ++psrc) {
			std::vector<bool> source_terminals(ctx->boxes.size(), false);
			extract_source_terminals(*psrc, source_terminals);
			VERBOSE(*verb_out << "  source " << *psrc << " gives terminals " << source_terminals << std::endl);
			if (source_terminals == self_key.terminals) {
				VERBOSE(*verb_out << "EXPLORING ref " << *pref << " and source " << *psrc << std::endl);
				double score;
				double rel_score;

				std::vector<bool> child_terminals(ctx->boxes.size(), false);
				extract_source_left_terminals(*psrc, child_terminals);

				if (!ctx->dirty.empty() && pref->parents.second.gid != CNFGrammar::INVALID_GID
				 && partial_intersection(ctx, child_terminals, 0, ctx->dirty.size())) {
					VERBOSE(*verb_out << "   this source splits a dirty region; can't use it\n");
					continue;
				}

				tree_iterator *left = build_subiter(pref->parents.first, child_terminals, nbranches + ((pref->parents.second.gid != CNFGrammar::INVALID_GID) ? 1 : 0));

				tree_iterator *right = 0;
				expression_node *left_tree = left->tree_(0, true);
				if (!left_tree) continue;

				VERBOSE(*verb_out << "BACK to ref " << *pref << " and source " << *psrc << std::endl);

				double left_score = left_tree->score();
				VERBOSE(*verb_out << "   left_score = " << left_score << std::endl);

				if (pref->parents.second.gid != CNFGrammar::INVALID_GID) {
					terminal_mask_t left_terminals = child_terminals;

					std::fill(child_terminals.begin(), child_terminals.end(), false);
					extract_source_right_terminals(*psrc, child_terminals);

					if (!ctx->dirty.empty()
					 && (partial_intersection(ctx, left_terminals, 0, ctx->dirty.size())
					  || partial_intersection(ctx, child_terminals, 0, ctx->dirty.size()))) {
						VERBOSE(*verb_out << "   this source splits a dirty region; can't use it\n");
						continue;
					}

					right = build_subiter(pref->parents.second, child_terminals, nbranches + 1);

					if (next_indices_exist(left, right)) {
						continue;
					}

					expression_node *right_tree = right->tree_(0, true);
					if (!right_tree) {
						continue;
					}

					VERBOSE(*verb_out << "BACK to ref " << *pref << " and source " << *psrc << std::endl);

					unsigned nleft = (left_tree->type() == TERMINAL_EXPR) ? NUM_TERMINAL_PREVIEWS : NUM_NONTERMINAL_PREVIEWS;
					unsigned nright = (right_tree->type() == TERMINAL_EXPR) ? NUM_TERMINAL_PREVIEWS : NUM_NONTERMINAL_PREVIEWS;

					for (unsigned ileft = 0; ileft < nleft && left_tree; ++ileft, left_tree = left->tree_(ileft, true)) {
						right_tree = right->tree_(0, true);
						for (unsigned iright = 0; iright < nleft && right_tree; ++iright, right_tree = right->tree_(iright, true)) {
							rel_score = get_relation(pref->spec.prod->rel, left_tree->rel_class(), right_tree->rel_class())->membership(left->bounds, right->bounds);
							score = cross_scores_rel(left_tree->score(), right_tree->score(), rel_score);
							score *= pref->spec.prod->scale;
							VERBOSE(*verb_out << "   got trees at indices (" << ileft << "," << iright << ") : scores " << left_tree->score() << "," << right_tree->score() << "," << rel_score << " -> " << score << std::endl);
							if (score > std::numeric_limits<double>::epsilon()) {
								get_next_indices(left, right);
								local_scores.insert(iteration_link(score, left, ileft, right, iright, 0, rel_score, pref->spec.prod));
							}
						}
					}
				}
				else {
					if (next_indices_exist(left, 0)) {
						continue;
					}
					score = left_score;
					score *= pref->spec.prod->scale;
					if (score > std::numeric_limits<double>::epsilon()) {
						get_next_indices(left, right);
						local_scores.insert(iteration_link(score, left, 0, 0, 0, 0, 0.0, pref->spec.prod));
					}
				}
			}
		}
	}

	VERBOSE(*verb_out << "after EXHAUSTING refs at " << self_key.root << ", there are " << local_scores.size() << " local score\n");

	if (!local_scores.empty()) {
		const iteration_link &best_score = *local_scores.begin();

		expression_node *tree = build_expression(best_score, true);
		VERBOSE(*verb_out << " seeding with best tree " << tree->str() << " for key " << self_key << std::endl);
		ranked_trees.push_back(tree);

		is_terminal = best_score.left->is_terminal;
		if (best_score.right) {
			is_terminal = is_terminal && best_score.right->is_terminal;
			is_terminal = is_terminal && (best_score.P->rel == StrokeOrder);
		}
		VERBOSE(*verb_out << " marking is_terminal=" << is_terminal << " at " << self_key << std::endl);
		VERBOSE(*verb_out << " marking is_distinct=" << is_distinct_subtree << " at " << self_key << std::endl);

		increment_next_indices(best_score);
	}
}



static void
insert_index_pair(std::set<std::pair<unsigned, unsigned> > &set, const std::pair<unsigned, unsigned> &idxs)
{
	std::set<std::pair<unsigned, unsigned> >::iterator pci;
	for (pci = set.begin(); pci != set.end(); ++pci) {
		if (*pci == idxs) break;
		if (pci->first == idxs.first) {
			if (idxs.second < pci->second) {
				set.erase(pci);
				set.insert(idxs);
			}
			break;
		}
		if (pci->second == idxs.second) {
			if (idxs.first < pci->first) {
				set.erase(pci);
				set.insert(idxs);
			}
			break;
		}
	}
	if (pci == set.end()) {
		set.insert(idxs);
	}
}

expression_node *
tree_iterator::build_expression(const iteration_link &score, bool building)
{
	expression_node *tree;

	expression_node *left = score.left->tree_(score.lefti, building);
	expression_node *right = score.right ? score.right->tree_(score.righti, building) : 0;

	if (!left || (score.right && !right)) {
		return 0;
	}
	tree = construct_expression(score.P, left, right);

	tree->set_context(const_cast<ParseContext *>(ctx));
	tree->set_score(score.score);
	tree->set_strokes(self_key.terminals);
	tree->set_box(bounds);
	tree->set_cellref(self_key.root);
	tree->set_grammar_type(self_key.root.gid);
	if (!score.right) {
		VERBOSE(*verb_out << " at root " << self_key.root << ", adding rec-child " << score.left->self_key.root << std::endl);
		tree->set_recursive_children(left->recursive_children());
		tree->add_recursive_child(score.left->self_key.root.gid);
		VERBOSE(
			*verb_out << "  children are now : ";
			for (std::vector<CNFGrammar::gid_t>::const_iterator i = tree->recursive_children().begin(); i != tree->recursive_children().end(); ++i) {
				*verb_out << *i << ", ";
			}
			*verb_out << std::endl;
		);
	}

	return tree;
}


void
tree_iterator::increment_next_indices(const iteration_link &score)
{
	next_index_cache_t::iterator pnxi;

	if (score.right) {
		pnxi = next_indices.find(std::make_pair(score.left->key(), score.right->key()));
	}
	else {
		pnxi = next_indices.find(std::make_pair(score.left->key(), cell_key::NIL));
	}

	assert(pnxi != next_indices.end());

	next_indices_t &nxi = pnxi->second;

	if (score.right) {
		assert(score.righti < nxi.size());
		if (score.lefti >= nxi[score.righti]) {
			VERBOSE(*verb_out << " incrementing right-index pair (" << nxi[score.righti] << "," << score.righti << ") to (");
			nxi[score.righti] = score.lefti + 1;
			VERBOSE(*verb_out << nxi[score.righti] << "," << score.righti << ")\n");
		}
		if (nxi.size() - 1 == score.righti) {
			VERBOSE(*verb_out << " inserting new right-index pair (0," << nxi.size() << ")\n");
			nxi.push_back(0);
		}
	}
	else {
		++nxi[0];
	}
}



expression_node *
tree_iterator::tree(unsigned rank)
{
	return tree_(rank, false);
}


expression_node *
tree_iterator::tree_(unsigned rank, bool building)
{
	if (self_key.root.gid == 2 && self_key.root.dim == 0 && self_key.root.start == 0 && self_key.root.end == 2) {
		int a = 0;
	}
	VERBOSE(*verb_out << "asking for tree " << rank << " at node " << self_key << "; last_clean=" << last_clean << "; terminal=" << is_terminal << "; distinct=" << is_distinct_subtree << "; flags (" << shared.prevent_terminal_alternates << "," << shared.prevent_subtree_alternates << ")" << std::endl);

	if (last_clean != ctx->dirty.size()) {
		VERBOSE(*verb_out << " this iterator may be dirty...");
		if (terminals_dirty(ctx, self_key.terminals, last_clean, ctx->dirty.size())) {
			VERBOSE(*verb_out << " it is; rebuilding...\n");
			next_indices.clear();
			local_scores.clear();
			while (!ranked_trees.empty()) {
				delete ranked_trees.back();
				ranked_trees.pop_back();
			}
	
			build(self_key.root);
		}
		else {
			VERBOSE(*verb_out << " it's not\n");
		}
		last_clean = ctx->dirty.size();
	}

	bool allow_only_transient_productions = false;
	
	std::map<terminal_mask_t, lock_descriptor>::const_iterator pldesc = ctx->locks.find(self_key.terminals);
	if (pldesc != ctx->locks.end()) {
		const lock_descriptor &ldesc = pldesc->second;
		VERBOSE(*verb_out << " this node is locked; ");
		if (std::find(ldesc.valid_gids.begin(), ldesc.valid_gids.end(), self_key.root.gid) == ldesc.valid_gids.end()) {
			VERBOSE(*verb_out << "gid " << self_key.root.gid << " invalid here; returning null\n");
			allow_only_transient_productions = true;
			//return 0;
		}
		else if (rank == 0) {
			VERBOSE(*verb_out << "returning " << ldesc.expr->str() << std::endl);
			return ldesc.expr;
		}
		else {
			VERBOSE(*verb_out << "returning null\n");
			return 0;
		}
	}

	bool valid = true;

	if (rank > 0) {
		if (is_terminal && shared.prevent_terminal_alternates) {
			if (building) {
				valid = false;
			}
			else {
				VERBOSE(*verb_out << " this node is marked terminal; cannot descend; returning null\n");
				return 0;
			}
		}
		if (is_distinct_subtree && shared.prevent_subtree_alternates) {
			if (building) {
				valid = false;
			}
			else {
				VERBOSE(*verb_out << " this node is marked distinct; cannot descend; returning null\n");
				return 0;
			}
		}
	}

	if (rank < ranked_trees.size()) {
		VERBOSE(*verb_out << " tree of rank " << rank << " has already been computed\n");
		for (std::vector<expression_node *>::iterator i = ranked_trees.begin(); i != ranked_trees.end(); ++i) {
			if ((*i)->valid()) {
				if (rank == 0) {
					VERBOSE(*verb_out << " tree exists: " << (*i)->str() << std::endl);
					return *i;
				}
				--rank;
			}
		}
		VERBOSE(*verb_out << " tree of rank " << rank << " used to exist, but flag conditions invalidated it\n");
		return 0;
	}

	VERBOSE(*verb_out << " asking for rank " << rank << "...there are " << ranked_trees.size() << " computed and " << local_scores.size() << " pending scores\n");
	assert(rank == ranked_trees.size());  // we should only be incrementing a single rank at a time

	if (CNFGrammar::is_terminal(self_key.root.gid) || local_scores.empty()) {
		return 0;
	}

	VERBOSE(*verb_out << " checking in next index pairs...");

	iteration_link prev_best = *local_scores.begin();
	local_scores.erase(local_scores.begin());

	next_index_cache_t::iterator pnxi;
	if (prev_best.right) {
		VERBOSE(*verb_out << " looking up index cache for key " << self_key << " in pair " << prev_best.left->key() << " , " << prev_best.right->key() << std::endl);
		pnxi = next_indices.find(std::make_pair(prev_best.left->key(), prev_best.right->key()));
	}
	else {
		VERBOSE(*verb_out << " looking up index cache for key " << self_key << " in single " << prev_best.left->key() << std::endl);
		pnxi = next_indices.find(std::make_pair(prev_best.left->key(), cell_key::NIL));
	}

	assert(pnxi != next_indices.end());

	const next_indices_t &nxi = pnxi->second;

	unsigned righti = 0;
	for (std::vector<unsigned>::const_iterator plefti = nxi.begin(); plefti != nxi.end(); ++plefti, ++righti) {
		VERBOSE(*verb_out << " at " << self_key << ", considering index pair " << *plefti << "," << righti << std::endl);
		expression_node *left_tree = prev_best.left->tree_(*plefti, building);
		if (!left_tree) continue;

		prev_best.lefti = *plefti;
		expression_node *right_tree = 0;
		if (left_tree && prev_best.right) {
			right_tree = prev_best.right->tree_(righti, building);
			if (!right_tree) continue;
			prev_best.righti = righti;
			prev_best.score = cross_scores_rel(left_tree->score(), right_tree->score(), prev_best.rel_score);
		}
		else {
			prev_best.score = left_tree->score();
		}

		// check whether this combination of right and left subtrees has already been found
		score_heap_t::iterator first = local_scores.lower_bound(prev_best);
		score_heap_t::iterator last = local_scores.upper_bound(prev_best);
		while (first != last) {
			--last;
			if (*last == prev_best) {
				break;
			}
		}

		// if it hasn't, add it
		if ((first == last) && (first == local_scores.end() || *first != prev_best)) {
			VERBOSE(*verb_out << " adding local score at " << self_key << " for " << prev_best.left->self_key);
			if (prev_best.right) {
				VERBOSE(*verb_out << " and " << prev_best.right->self_key);
			}
			VERBOSE(*verb_out << "; score " << prev_best.score << std::endl);
			local_scores.insert(prev_best);
		}
	}

	if (local_scores.empty()) {
		return 0;
	}

	expression_node *tree;

	do {
		const iteration_link &new_best = *local_scores.begin();

		tree = build_expression(new_best, building);
		if (tree) {
			tree->set_valid(valid);
			increment_next_indices(new_best);
			ranked_trees.push_back(tree);
			VERBOSE(*verb_out << "tree(" << rank << ") for root " << self_key.root << " gives tree " << tree->str() << std::endl);
		}
		else {
			local_scores.erase(local_scores.begin());
		}
	} while (!tree && !local_scores.empty());

	return tree;
}


const ExpressionTree *
expression_iterator::next()
{
	VERBOSE(*verb_out << "next()\n");

	std::vector<std::pair<tree_iterator *, unsigned> >::iterator best = iterators.end();
	best_tree = 0;
	for (std::vector<std::pair<tree_iterator *, unsigned> >::iterator i = iterators.begin(); i != iterators.end(); ++i) {
		expression_node *tree = i->first->tree(i->second);
		if (tree && (best == iterators.end() || tree->score() > best_tree->score())) {
			best_tree = tree;
			best = i;
		}
	}
	if (best != iterators.end()) {
		VERBOSE(*verb_out << "move_next() has top score " << best_tree->score() << " and tree " << best_tree->str() << std::endl);
		++best->second;
	}

	return best_tree ? DEBUG_NEW expression_node(*best_tree) : 0;
}



void
clear_lock(ParseContext *ctx, const terminal_mask_t &terminals)
{
	std::map<terminal_mask_t, lock_descriptor>::iterator i = ctx->locks.find(terminals);
	if (i != ctx->locks.end()) {
		clear_lock_descriptor(i->second);
		ctx->locks.erase(i);
	}

	std::vector<terminal_mask_t>::iterator j = std::find(ctx->dirty.begin(), ctx->dirty.end(), terminals);
	if (j != ctx->dirty.end()) {
		ctx->dirty.erase(j);
	}
}


void
clear_lock_descriptor(lock_descriptor &ldesc)
{
	delete ldesc.expr;
	ldesc.expr = 0;
	ldesc.valid_gids.clear();
}


void
clear_lock_table(ParseContext *ctx)
{
	while (!ctx->locks.empty()) {
		std::map<terminal_mask_t, lock_descriptor>::iterator i = ctx->locks.begin();
		clear_lock_descriptor(i->second);
		ctx->locks.erase(i);
	}

	ctx->dirty.clear();
}


}
