#include "parser.h"
#include "fuzzy.h"
#include "grammar.h"
#include "rect.h"
#include "verb.h"
#include "vector.h"
#include "relation.h"
#include "parms.h"
#include "grammar-values.h"
#include "MatrixAnalyzer.h"
#include "expr-iter.h"
#include "MathRecognizer.h"
#include "mathrecognizer-private.h"

#include <deque>
#include <cassert>
#include <cstdlib>
#include <functional>



namespace scg
{


#if 0
static double SCORE_THRESHOLD = RegisterParameterDouble("ScoreThreshold", &SCORE_THRESHOLD);


struct parsepos {
	const production *P;
	size_t pos;
	size_t segpos;
	parse_links links;

	parsepos(const production *P_, size_t pos_, size_t segpos_) : P(P_), pos(pos_), segpos(segpos_), links(P_) { }
	parsepos(const production *P_, size_t pos_, size_t segpos_, const parse_links &links_) : P(P_), pos(pos_), segpos(segpos_), links(links_) { }
};


struct termparsepos {
	const production *P;
	size_t termpos;
	size_t segpos;
	std::vector<parse_links> links;

	termparsepos(const production *P_, size_t termpos_, size_t segpos_) : P(P_), termpos(termpos_), segpos(segpos_) { }
	termparsepos(const production *P_, size_t termpos_, size_t segpos_, const std::vector<parse_links> &links_) : P(P_), termpos(termpos_), segpos(segpos_), links(links_) { }
};


static void invalidate_cell(parse_cell *cell);
static void invalidate_single_cell(parse_cell *cell);

static std::ostream &
operator<<(std::ostream &os, const parsepos &pos)
{
	os << "production: " << *pos.P << "; Ppos:" << pos.pos << "; Ipos:" << pos.segpos;
	return os;
}

static std::ostream &
operator<<(std::ostream &os, const termparsepos &pos)
{
	os << "production: " << *pos.P << "; Tpos:" << pos.termpos << "; Ipos:" << pos.segpos;
	return os;
}


static std::ostream &
operator<<(std::ostream &os, const parseref &ref)
{
	if (ref.nt) {
		os << ref.nt->name << ':' << ref.span->bits;
	}
	else {
		os << "unknown";
	}
	return os;
}

std::ostream &
operator<<(std::ostream &os, const parse_links &links)
{
	//os << "<" << ((links.P && links.P->rel) ? links.P->rel->name : "nil") << "> [ ";
	os << "<" << ((links.P && links.P->rel) ? links.P->rel->name : "nil") << "> [ ";
	for (std::vector<parseref>::const_iterator i = links.children.begin(); i != links.children.end(); ++i) {
		os << *i << ' ';
	}
	return os << ']';
	//return os << "] (" << links.score << ")";
}



static bool parse_internal(context *ctx, ordered_segments *segs, const nonterminal *nt, bool ign_overlap);

double
compute_membership(const relation *rel, const ordered_segments *segs1, const ordered_segments *segs2, const std::set<rclass_t> &cls1, const std::set<rclass_t> &cls2)
{
	//std::vector<double> M;
	//measure(b1, b2, M);

	double mem = 0.0;

	VERBOSE(
		*verb_out << "compute membership in " << rel->name << " from " << segs1->bounds() << " to " << segs2->bounds() << " with classes { ";
		for (std::set<rclass_t>::const_iterator i = cls1.begin(); i != cls1.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "} and { ";
		for (std::set<rclass_t>::const_iterator i = cls2.begin(); i != cls2.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "} ";
	);
	
	if (cls1.empty()) {
		if (cls2.empty()) {
			//return rel->membership(M, BOX_CLASS, BOX_CLASS);
			mem = rel->membership(segs1->bounds(), segs2->bounds(), AGGREGATE_CLASS, AGGREGATE_CLASS);
		}
		else {
			for (std::set<rclass_t>::const_iterator c2 = cls2.begin(); c2 != cls2.end(); ++c2) {
				//mem = std::max(mem, rel->membership(M, BOX_CLASS, *c2));
				mem = std::max(mem, rel->membership(segs1->bounds(), segs2->bounds(), AGGREGATE_CLASS, *c2));
			}
		}
	}

	else if (cls2.empty()) {
		for (std::set<rclass_t>::const_iterator c1 = cls1.begin(); c1 != cls1.end(); ++c1) {
			//mem = std::max(mem, rel->membership(M, *c1, BOX_CLASS));
			mem = std::max(mem, rel->membership(segs1->bounds(), segs2->bounds(), *c1, AGGREGATE_CLASS));
		}
	}
	
	else {
		for (std::set<rclass_t>::const_iterator c1 = cls1.begin(); c1 != cls1.end(); ++c1) {
			for (std::set<rclass_t>::const_iterator c2 = cls2.begin(); c2 != cls2.end(); ++c2) {
				//mem = std::max(mem, rel->membership(M, *c1, *c2));
				mem = std::max(mem, rel->membership(segs1->bounds(), segs2->bounds(), *c1, *c2));
			}
		}
	}

	VERBOSE(*verb_out << " : " << mem << std::endl);
	
	return mem;
}

static std::set<rclass_t> EMPTY_CLASS_SET;


ordered_segments *
lookup_span(context *ctx, const bitvec &bits, ordered_segments *src, size_t dim, size_t start, size_t end)
{
	assert(ctx->spans);

	/*
	VERBOSE(
		*verb_out << "looking for span " << bits << " in options" << std::endl;
		for (std::map<bitvec, ordered_segments *>::iterator i = ctx->spans->begin(); i != ctx->spans->end(); ++i) {
			*verb_out << "  " << i->first << std::endl;
		}
	);
	*/

	std::map<bitvec, ordered_segments *>::iterator i = ctx->spans->find(bits);
	if (i == ctx->spans->end()) {
		if (src) {
			//ordered_segments *span = new ordered_segments(ctx->segments, bits, ctx->edits.size());
			ordered_segments *span = src->slice(dim, start, end);
			VERBOSE(*verb_out << "creating new span " << span << " at " << bits << std::endl);
			(*ctx->spans)[bits] = span;
			return span;
		}
		else {
			VERBOSE(*verb_out << "span not found.\n");
			return 0;
		}
	}
	
	return i->second;
}

static bool parse_production(context *ctx, parse_cell *cell, ordered_segments *segs, parse_links &links, size_t curr, bool ign_overlap);

static bool
parse_nts(context *ctx, parse_cell *cell, parse_links &links, size_t curr, size_t term_index,
		  ordered_segments *segs, ordered_segments *postsegs,
		  ordered_segments *termsegs, const nonterminal *termnt, bool ign_overlap) {
	const production *P = links.P;
	size_t d = P->dim;

	if (curr == term_index) {
		VERBOSE(*verb_out << "finished parsing ntstring in " << *P << " from " << curr << " to " << term_index << " with links " << links << std::endl);
		if (curr == P->rhs.size()) {
			assert(!termsegs && !termnt);
			cell->links.push_back(links);
			return true;
		}
		else {
			assert(termsegs && termnt);
			links.children[curr] = parseref(termsegs, termnt);
			return parse_production(ctx, cell, postsegs, links, curr + 1, ign_overlap);
		}
	}

	VERBOSE(*verb_out << "parsing ntstring in " << *P << " from " << curr << " to " << term_index << " on " << segs->bits << " with links " << links << std::endl);

	bool suc = false;
	size_t leftover = segs->size();
	for (size_t i = curr + 1; i < term_index; ++i) {
		if (P->rhs[i]->max_length > leftover) {
			leftover = 0;
			break;
		}
		leftover -= P->rhs[i]->max_length;
	}
	size_t minlen = std::max(P->rhs[curr]->min_length, leftover);
	VERBOSE(*verb_out << P->rhs[curr]->max_length << ' ' << segs->size() << ' ' << P->min_prefix_length[term_index-1] << ' ' << P->min_prefix_length[curr] << std::endl);
	size_t maxlen = std::min<size_t>(P->rhs[curr]->max_length, segs->size() - (P->min_prefix_length[term_index-1] - P->min_prefix_length[curr]));
	VERBOSE(*verb_out << "minlen=" << minlen << "; maxlen=" << maxlen << std::endl);
	for (size_t len = minlen; len <= maxlen; ++len) {
		bitvec currbits = segs->slice_bits(d, 0, len);
		ordered_segments *currsegs = lookup_span(ctx, currbits, segs, d, 0, len);
		if (currsegs) {
			if (parse_internal(ctx, currsegs, P->rhs[curr], ign_overlap)) {
				links.children[curr] = parseref(currsegs, P->rhs[curr]);
				if (len == segs->size()) {
					assert(curr == term_index - 1);
					suc = parse_nts(ctx, cell, links, curr + 1, term_index, 0, postsegs, termsegs, termnt, ign_overlap) || suc;
				}
				else {
					bitvec rembits = segs->slice_bits(d, len, segs->size());
					ordered_segments *remsegs = lookup_span(ctx, rembits, segs, d, len, segs->size());
					if (remsegs) {
						suc = parse_nts(ctx, cell, links, curr + 1, term_index, remsegs, postsegs, termsegs, termnt, ign_overlap) || suc;
					}
				}
			}
		}
	}
	return suc;
}

static bool
parse_production(context *ctx, parse_cell *cell, ordered_segments *segs, parse_links &links, size_t curr, bool ign_overlap) {
	const production *P = links.P;
	size_t d = P->dim;

	if (curr == P->rhs.size()) {
		VERBOSE(*verb_out << "  ADDING PARSE " << links << std::endl);
		cell->links.push_back(links);
		return true;
	}

	VERBOSE(*verb_out << "parsing production " << *P << " on " << segs->bits << " at " << curr << " with links " << links << std::endl);

	std::vector<size_t>::const_iterator next_term;
	next_term = std::lower_bound(P->terminal_indices.begin(), P->terminal_indices.end(), curr);
	if (next_term == P->terminal_indices.end()) {
		return parse_nts(ctx, cell, links, curr, P->rhs.size(), segs, 0, 0, 0, ign_overlap);
	}
	else {
		bool suc = false;
		size_t term_index = *next_term;
		const nonterminal *termnt = P->rhs[term_index];

		size_t min_prelen;
		if (curr == 0) {
			if (term_index == 0) {
				min_prelen = 0;
			}
			else {
				min_prelen = P->min_prefix_length[term_index-1];
			}
		}
		else {
			min_prelen = P->min_prefix_length[term_index-1] - P->min_prefix_length[curr-1];
		}

		size_t max_postlen = 0;
		for (size_t i = term_index + 1; i < P->rhs.size(); ++i) {
			if (P->rhs[i]->max_length > segs->size()) {
				max_postlen = segs->size();
				break;
			}
			max_postlen += P->rhs[i]->max_length;
			if (max_postlen >= segs->size()) {
				max_postlen = segs->size();
				break;
			}
		}
		//max_postlen -= min_prelen + termnt->min_length;

		/*size_t leftover = segs->size() - max_postlen;
		if (leftover > termnt->max_length) {
			leftover -= termnt->max_length;
		}
		else {
			leftover = 0;
		}
		min_prelen = std::max(min_prelen, leftover);*/
		size_t min_postlen = P->min_prefix_length[P->rhs.size()-1] - P->min_prefix_length[term_index];
		size_t max_prelen = segs->size() - termnt->min_length - min_postlen;
		size_t maxstart = 0;
		for (size_t i = curr; i < term_index; ++i) {
			if (P->rhs[i]->max_length > segs->size()) {
				maxstart = segs->size();
				break;
			}
			maxstart += P->rhs[i]->max_length;
			if (maxstart >= segs->size()) {
				maxstart = segs->size();
				break;
			}
		}
		max_prelen = maxstart;//std::min(max_prelen, maxstart);

		assert(min_prelen >= 0 && min_prelen <= segs->size());
		assert(max_prelen >= 0 && max_prelen <= segs->size());
		assert(min_postlen >= 0 && min_postlen <= segs->size());
		assert(max_postlen >= 0 && max_postlen <= segs->size());

		VERBOSE(*verb_out << "term_index=" << term_index << "; min_prelen=" << min_prelen << "; min_postlen=" << min_postlen << "; max_prelen=" << max_prelen << "; max_postlen=" << max_postlen << std::endl);
		for (size_t prelen = min_prelen; prelen <= max_prelen; ++prelen) {
			//assert(min_postlen + prelen + termnt->min_length <= segs->size());
			//size_t min_termlen = std::max(termnt->min_length, segs->size() - max_postlen - prelen);
			//size_t max_termlen = std::min(termnt->max_length, segs->size() - min_postlen - prelen);
			//VERBOSE(*verb_out << "min_termlen=" << min_termlen << "; max_termlen=" << max_termlen << std::endl);
			for (size_t termlen = termnt->min_length; termlen <= termnt->max_length; ++termlen) {
				VERBOSE(*verb_out << "  trying prelength " << prelen << "; terminal length " << termlen << std::endl);
				if (prelen + termlen > segs->size()) {
					break;
				}
				if (segs->size() - (prelen + termlen) > max_postlen) {
					continue;
				}
				if (segs->size() - (prelen + termlen) < min_postlen) {
					break;
				}
				//assert(prelen + termlen + min_postlen <= segs->size());
				/*if (prelen + termlen + min_postlen > segs->size()) {
					//abort();
					VERBOSE(*verb_out << "   cannot fit in segments with total min length " << prelen + termlen + min_postlen << std::endl);
					break;
				}*/
				bitvec termbits = segs->slice_bits(d, prelen, prelen + termlen);
				ordered_segments *termsegs = lookup_span(ctx, termbits, segs, d, prelen, prelen + termlen);
				if (termsegs) {
					const parse_cell *term_cell = ctx->table[termsegs][termnt];
					/*VERBOSE(
						*verb_out << "cell is " << term_cell << std::endl;
						*verb_out << "current options in span " << termsegs << " are: ";
						const nonterminal_parses &opts = ctx->table[termsegs];
						for (nonterminal_parses::const_iterator i = opts.begin(); i != opts.end(); ++i) {
							*verb_out << "  " << i->first->name << " -> " << i->second << std::endl;
						}
					);*/

					if (term_cell) {
						VERBOSE(*verb_out << "     terminal " << termnt->name << " exists here\n");
						if (prelen > 0) {
							bitvec prebits = segs->slice_bits(d, 0, prelen);
							ordered_segments *presegs = lookup_span(ctx, prebits, segs, d, 0, prelen);
							double local_score = compute_membership(P->rel, presegs, termsegs, EMPTY_CLASS_SET, termnt->rclasses);
							if (local_score > SCORE_THRESHOLD) {
								if (min_postlen == 0) {
									suc = parse_nts(ctx, cell, links, curr, term_index, presegs, 0, termsegs, termnt, ign_overlap) || suc;
								}
								else {
									bitvec postbits = segs->slice_bits(d, prelen + termlen, segs->size());
									ordered_segments *postsegs = lookup_span(ctx, postbits, segs, d, prelen + termlen, segs->size());
									if (postsegs) {
										suc = parse_nts(ctx, cell, links, curr, term_index, presegs, postsegs, termsegs, termnt, ign_overlap) || suc;
									}
								}
							}
						}
						else {
							links.children[curr] = parseref(termsegs, termnt);
							if (min_postlen == 0) {
								suc = parse_production(ctx, cell, 0, links, curr + 1, ign_overlap) || suc;
							}
							else {
								bitvec postbits = segs->slice_bits(d, prelen + termlen, segs->size());
								ordered_segments *postsegs = lookup_span(ctx, postbits, segs, d, prelen + termlen, segs->size());
								if (postsegs) {
									suc = parse_production(ctx, cell, postsegs, links, curr + 1, ign_overlap) || suc;
								}
							}
						}
					}
				}
			}
		}
		return suc;
	}
}

static bool
parse_internal(context *ctx, ordered_segments *segs, const nonterminal *nt, bool ign_overlap) {
	parse_cell *&cell = ctx->table[segs][nt];
	if (cell) {
		VERBOSE(
			*verb_out << "already parsed at " << nt->name << ":" << segs->bits << " and there are ";
			if (cell->links.empty()) {
				*verb_out << "no ";
			}
			*verb_out << "results\n";
		);
		return !cell->links.empty();
	}

	/*if (!ign_overlap) {
		for (size_t i = 0; i < segs->bits.size(); ++i) {
			if (!segs->bits.at(i)) {
				Rect<long> exbox = ctx->segments[i]->bounds;
				Rect<long> segbox = segs->bounds();
				if (overlap_proportion(exbox, segbox) >= 0.25) {
					VERBOSE(*verb_out << "not parsing " << nt->name << " at " << segs->bits << " because of overlap with segment " << i << std::endl);
					return false;
				}
			}
		}
	}*/

	cell = new parse_cell(segs, nt);

	VERBOSE(*verb_out << "parsing " << nt->name << " at " << segs->bits << "(" << nt << " at " << segs << ")" << std::endl);

	
	if (nt->name[0] == '!')  { // indicates matrix
		VERBOSE(*verb_out << "creating matrix for strokes " << segs->bits << std::endl);
		cell->matrix = CreateMatrixAnalyzer(ctx->wrapper);
		for (std::vector<const segment *>::const_iterator i = segs->begin(TIME_ORDER); i != segs->end(TIME_ORDER); ++i) {
			cell->matrix->addStroke(&(*i)->stk);
		}
		cell->mxiter = cell->matrix->createIterator();
		bool results = (cell->mxiter->next() != 0);
		delete cell->mxiter;
		cell->mxiter = 0;
		if (!results) {
			DestroyMatrixAnalyzer(cell->matrix);
		}
		return results;
	}

	for (std::vector<production *>::const_iterator i = nt->productions.begin(); i != nt->productions.end(); ++i) {
		const production *P = *i;
		if (P->min_prefix_length[P->rhs.size()-1] <= segs->size()) {
			if (P->rhs.size() > 1 || P->rhs[0]->max_length == 0 || segs->size() <= P->rhs[0]->max_length) {
				parse_links links(P);
				links.children.resize(links.P->rhs.size());
				parse_production(ctx, cell, segs, links, 0, ign_overlap || P->rel == get_relation(REL_CONTAINS));
			}
		}
	}

	
	if (!cell->links.empty()) {
		cell->classes.insert(AGGREGATE_CLASS);
		for (std::vector<parse_links>::iterator i = cell->links.begin(); i != cell->links.end(); ++i) {
			if (i->children.size() > 1) {
				cell->classes.insert(BOX_CLASS);
			}
			else {
				const parseref &ref = i->children.front();
				parse_cell *child = ctx->table[ref.span][ref.nt];
				cell->classes.insert(child->classes.begin(), child->classes.end());
			}
		}
	}

	VERBOSE(
		*verb_out << "returning from " << nt->name << " at " << segs->bits << "; there are " << cell->links.size() << " results, as follows:\n";
		for (std::vector<parse_links>::const_iterator i = cell->links.begin(); i != cell->links.end(); ++i) {
			*verb_out << ' ' << *i << std::endl;
		}
	);

	return !cell->links.empty();
}


bool
parse(context *ctx, ordered_segments *segs, const nonterminal *nt)
{
	VERBOSE(
		*verb_out << "PARSING...existing spans are\n";
		for (parse_table::const_iterator i = ctx->table.begin(); i != ctx->table.end(); ++i) {
			*verb_out << i->first->bits << " -> " << i->first << std::endl;
		}
	);
	//rebuild_production_lengths(ctx->G, ctx->G.start);
	return parse_internal(ctx, segs, nt, false);
}


class parser_subtree_accessor : public subtree_accessor {
public:
	parser_subtree_accessor(const parse_iterator &it_, const parse_links &links_, const std::vector<size_t> &ranks_, std::vector<expression_node> &nodes_)
		: it(it_), links(links_), ranks(ranks_), nodes(nodes_) { }

	expression_node operator[](size_t i) const {
		return extract_tree_internal(it, links.children[i].span, links.cells[i], ranks[i]);
	}

	bool verify(size_t i, const expression_node &node) const {
		return links.children[i].nt == node->nt_low || links.children[i].nt->has_direct_descendent(node->nt_low);
	}

private:
	const parse_iterator &it;
	const parse_links &links;
	const std::vector<size_t> &ranks;
	std::vector<expression_node> &nodes;
};


static expression_node
build_tree(const parse_iterator &it, const parse_links &links, const std::vector<size_t> &ranks, std::vector<expression_node> &nodes, score_combiner &score)
{
	nodes.resize(links.P->rhs.size());
	parser_subtree_accessor acc(it, links, ranks, nodes);
	expression_node tree = links.P->tbuild->build(acc);
	if (tree) {
		/*
		for (size_t i = 0; i < tree->nchildren(); ++i) {
			assert(tree->mchild(i)->parent);
		}*/

		tree->ranks = ranks;

		if (tree->nchildren() == 1) {
			tree->rclasses = tree->mchild(0)->rclasses;
			tree->set_matrix(tree->mchild(0)->matrix());
		}
		if (tree->rclasses.empty()) {
			tree->rclasses.insert(BOX_CLASS);
		}
		tree->rclasses.insert(AGGREGATE_CLASS);

		/*for (size_t i = 0; i < tree->nchildren(); ++i) {
			expression_node child = tree->mchild(i);
			if ((child->type() == InvalidSemanticId || child->type() == PLACEHOLDER_EXPR)
			 && child->score() > 0) {
				score = score.add(child->score_combo());
			}
		}*/
	}

	return tree;
}


static expression_node
build_tree(const parse_iterator &it, const parse_links &links, const std::vector<size_t> &ranks, score_combiner &score)
{
	std::vector<expression_node> dummy;
	return build_tree(it, links, ranks, dummy, score);
}


expression_node
make_matrix_expression(parse_cell *mxcell, MatrixAnalyzer *mxan, Matrix *mx)
{
	expression_node tree(new raw_expression_node);
	tree->set_type(MATRIX_EXPR);

	tree->set_matrix(mxan);

	size_t nrows, ncols;
	mx->getDimensions(nrows, ncols);

	std::stringstream ss, latex;
	ss << nrows;
	tree->add_child(make_terminal(ss.str(), 1.0, false));
	ss.str("");
	ss << ncols;
	tree->add_child(make_terminal(ss.str(), 1.0, false));
	ss.str("");

	score_combiner sc;
	expression_node rows(new raw_expression_node);

	for (size_t i = 0; i < nrows; ++i) {
		std::stringstream row_ss;
		std::stringstream latex_ss;
		row_ss << "<mtr>";
		score_combiner row_sc;
		expression_node row(new raw_expression_node);
		row->set_type(MATRIXROW_EXPR);
		//bitvec bits(ctx->segments.size(), false);
		for (size_t j = 0; j < ncols; ++j) {
			ExpressionIterator *pit;
			VERBOSE(*verb_out << "make_matrix_expression: getting iterator for row " << i << ", col " << j << std::endl);
			mx->getCell(i, j, pit);
			if (!pit) {
				return NULL_EXPRESSION_NODE;
			}
			iterator_base *it = (iterator_base *)pit;
			parse_cell *parent_cell = it->get_cell();
			if (parent_cell) {
				parent_cell->parents.push_back(mxcell);
			}
			expression_node cell = it->next_node();
			if (!cell) {
				return NULL_EXPRESSION_NODE;
			}
			row->add_child(cell);
			row_sc = row_sc.add(cell->score_combo());
			//row_sc = row_sc.addterm(cell->unnorm_score());
			//sc.addterm(cell->unnorm_score());
			row_ss << "<mtd>" << cell->long_str() << "</mtd>";
			latex_ss << cell->latex_str();
			if (j != ncols - 1) {
				latex_ss << " & ";
			}
			//bits.union_insert(cell->span->bits);
		}
		row_ss << "</mtr>";
		latex_ss << " \\\\\n";

		row->set_score(row_sc);
		sc = sc.add(row_sc);
		row->set_cell(0);
		//row->set_span(lookup_span(ctx, bits));
		row->set_nt(0);
		row->nt_low = 0;
		row->prod = 0;
		row->rclasses.insert(BOX_CLASS);
		row->rclasses.insert(AGGREGATE_CLASS);
		row->set_score(row_sc);

		row->set_latex(latex_ss.str());
		row->set_long_string(row_ss.str());
		rows->add_child(row);

		ss << row_ss.str();
		latex << latex_ss.str();
	}

	rows->set_latex(latex.str());
	rows->set_long_string(ss.str());
	rows->rclasses.insert(BOX_CLASS);
	rows->rclasses.insert(AGGREGATE_CLASS);
	VERBOSE(*verb_out << "matrix score: " << sc.nrel << " rels; " << sc.nterm << " terms ; " << sc.running_rel << " ; " << sc.running_term << " = " << sc.score() << std::endl);
	rows->set_score(sc);

	tree->add_child(rows);
	
	//tree->set_long_string(ss.str());
	tree->rclasses.insert(BOX_CLASS);
	tree->rclasses.insert(AGGREGATE_CLASS);
	tree->set_score(sc);

	tree->set_span(mxcell->segs);
	tree->head_span = tree->tail_span = mxcell->segs;

	VERBOSE(*verb_out << "MATRIX expression is " << tree->long_str() << " with score " << tree->score() << std::endl);
	
	return tree;
}


static expression_node
extract_matrix_tree(context *ctx, parse_cell *cell, ordered_segments *span)
{
	Matrix *mx = cell->mxiter->next();
	if (!mx) {
		return NULL_EXPRESSION_NODE;
	}
	
	mx->rebuild();
	expression_node tree = make_matrix_expression(cell, cell->matrix, mx);
	if (!tree) {
		return NULL_EXPRESSION_NODE;
	}
	tree->set_nt(cell->nt);
	tree->nt_low = cell->nt;
	tree->prod = 0;

	expression_node rows = tree->mchild(MATRIX_ROWS);
	rows->set_nt(cell->nt);
	rows->nt_low = cell->nt;
	rows->set_context(ctx);
	rows->set_cell(cell);
	rows->set_span(span);

	for (size_t i = 0; i < rows->nchildren(); ++i) {
		expression_node row = rows->mchild(i);
		row->set_context(ctx);
	}

	double conf;
	cell->matrix->getConfidence(mx, conf);
	//score_combiner sc(conf);
	//sc.add(rows->score_combo());
	score_combiner sc = rows->score_combo();
	sc.addrel(conf);
	tree->set_score(sc);

	return tree;
}


static expression_node extract_tree(const parse_iterator &it, parse_links &links, size_t rank);


static double
measure_relation(const parse_links &links, size_t from, size_t to, const expression_node &from_tree, const expression_node &to_tree, rclass_t from_class, rclass_t to_class)
{
	assert(links.children.size() == links.P->rhs.size());
	
	const relation *rel = links.P->rel;
	const parseref &fromref = links.children[from];
	const parseref &toref = links.children[to];
	
	const attach_mode &mode = links.P->attach_modes[from];
	
	VERBOSE(*verb_out << "in measure_relation, mode is ");
	const ordered_segments *from_segs;
	if (mode.from == attach_mode::GROUP) {
		from_segs = fromref.span;
		VERBOSE(*verb_out << "GROUP " << from_segs->bounds() << " -> ");
	}
	else {
		assert(mode.from == attach_mode::SYMBOL);
		from_segs = from_tree->tail_span;
		VERBOSE(*verb_out << "SYMBOL " << from_segs->bounds() << " -> ");
	}
	
	const ordered_segments *to_segs;
	if (mode.to == attach_mode::GROUP) {
		to_segs = toref.span;
		VERBOSE(*verb_out << "GROUP " << to_segs->bounds() << " : ");
	}
	else {
		assert(mode.to == attach_mode::SYMBOL);
		to_segs = to_tree->head_span;
		VERBOSE(*verb_out << "SYMBOL " << to_segs->bounds() << " : ");
	}
	
	VERBOSE(*verb_out << std::endl);
	return rel->membership(from_segs->bounds(), to_segs->bounds(), from_class, to_class);
}


static double
getrelscore(const parse_links &links, size_t i, const expression_node &from_tree, const expression_node &to_tree) {
	std::set<rclass_t>::const_reverse_iterator from_class = from_tree->rclasses.rbegin();
	std::set<rclass_t>::const_reverse_iterator to_class = to_tree->rclasses.rbegin();
	int which = 1;
	while (from_class != from_tree->rclasses.rend() && to_class != to_tree->rclasses.rend()) {
		double s = measure_relation(links, i - 1, i, from_tree, to_tree, *from_class, *to_class);
		if (s != -1.0) {
			return s;
		}
		if (which == 1) {
			++from_class;
			which = 2;
		}
		else {
			++to_class;
			which = 1;
		}
	}
	return -1.0;
}

template <typename T>
static void
extract_multiclass_trees(const parse_iterator &it, parse_links &links, const score_combiner &score, std::vector<expression_node> children,
                         std::vector<size_t> &ranks, T first, T curr, T last,
						 ordered_segments *head_span, ordered_segments *tail_span)
{
	VERBOSE(*verb_out << "extract_multiclass_tree() at " << links << " with score " << score.score() << std::endl);

	size_t firsti;
	if (curr == first) {
		firsti = 0;
	}
	else {
		firsti = *(curr - 1) + 1;
	}

	score_combiner final_score = score;
	if (curr == last) {
		assert(children.size() == firsti);
		
		for (size_t i = firsti; i < links.children.size(); ++i) {
			VERBOSE(*verb_out << "extracting subtree ranked " << ranks[i] << " from " << links.children[i].nt->name << ':' << links.children[i].span->bits << std::endl);
			expression_node subtree = extract_tree_internal(it, links.children[i].span, links.cells[i], ranks[i]);
			if (!subtree) {
				return;
			}
			if (i == links.P->head) {
				head_span = subtree->head_span;
				assert(head_span);
			}
			if (i == links.P->tail) {
				tail_span = subtree->tail_span;
				assert(tail_span);
			}
			VERBOSE(*verb_out << "subtree is " << *subtree << std::endl);
			VERBOSE(*verb_out << " subtree score is " << subtree->score() << std::endl);
			final_score = final_score.add(subtree->score_combo());
			//final_score = final_score.addterm(subtree->unnorm_score());
			if (i > 0) {
				const expression_node &from_tree = children[i - 1];
				const expression_node &to_tree = subtree;
				double rscore = getrelscore(links, i, from_tree, to_tree);
				if (rscore <= SCORE_THRESHOLD) {
					return;
				}
				VERBOSE(*verb_out << "adding rel-score " << rscore << " at " << i << std::endl);
				//final_score = final_score.add(best_rel_score);
				final_score = final_score.addrel(rscore);
			}
			
			children.push_back(subtree);
		}

		VERBOSE(*verb_out << "parse score is " << final_score.score() << std::endl);

		VERBOSE(
			*verb_out << " at end; building tree with ranks [ ";
			for (std::vector<size_t>::const_iterator i = ranks.begin(); i != ranks.end(); ++i) {
				*verb_out << *i << ' ';
			}
			*verb_out << "]\n";
		);
		
		expression_node tree = build_tree(it, links, ranks, final_score);
		if (!tree) {
			return;
		}

		if (!links.P->sbuild && children.size() == 1 && links.P->nt->sid == InvalidSemanticId) {
			tree->set_string_builder(children[0]->get_string_builder());
			tree->set_latex_builder(children[0]->get_latex_builder());
		}
		else if (links.P->sbuild) {
			tree->set_long_string("");
			tree->set_string_builder(links.P->sbuild);
			tree->set_latex_builder(links.P->lbuild);
		}

		if (tree->type() == MATRIX_EXPR) {
			tree->set_matrix(children[0]->matrix());
		}
		else if (tree->type() == PAREN_EXPR) {
			if (children.size() == 1) {
				tree->set_matrix(children[0]->matrix());
			}
			else {
				tree->set_matrix(children[1]->matrix());
			}
		}

		//assert(tree);
		VERBOSE(*verb_out << "  tree is " << tree->long_str() << std::endl);

		VERBOSE(*verb_out << "  final_score is " << final_score.score() << std::endl);
		if (final_score.score() > SCORE_THRESHOLD) {
			if (tree->noperations == 0) {
				for (size_t i = 0; i < tree->nchildren(); ++i) {
					expression_node child = tree->mchild(i);
					VERBOSE(*verb_out << " subtree nops is " << child->noperations << std::endl);
					tree->noperations += child->noperations;
				}

				SemanticId sid = links.P->nt->sid;
				/*
				if (sid > 0 && sid != VAR_EXPR && sid != NUM_EXPR && sid != PAREN_EXPR) {
					VERBOSE(*verb_out << "the production " << links.P->nt->name << " counts as an operation!\n");
					++tree->noperations;
				}*/
				VERBOSE(*verb_out << "setting nops of " << *tree << " to " << tree->noperations << std::endl);
			}
			tree->set_score(final_score);
			if (tree->nt && tree->nt != links.P->nt) {
				tree = tree->copy_shallow();
				//tree = expression_node(new raw_expression_node(*tree));
				//tree = expression_node(tree);
			}
			if (!tree->nt_low || links.children.size() > 1) {
			//if (tree->type() > 0) {
				VERBOSE(*verb_out << "setting low NT with type " << tree->type() << " to " << links.P->nt->name << std::endl);
				tree->nt_low = links.P->nt;
				tree->parse_links = links.children;
				tree->prod = links.P;
			}
			tree->set_nt(links.P->nt);
			tree->head_span = head_span;//links.children[links.P->head].span;
			tree->tail_span = tail_span;//links.children[links.P->tail].span;
			VERBOSE(*verb_out << "setting head span to " << tree->head_span->bits << "; tail span to " << tree->tail_span->bits << std::endl);
			links.cached_trees.insert(links_parse_result(tree, ranks));
		}
	}
	else {
		assert(children.size() == firsti);
		for (size_t i = firsti; i < *curr; ++i) {
			expression_node subtree = extract_tree_internal(it, links.children[i].span, links.cells[i], ranks[i]);
			if (!subtree) {
				return;
			}
			if (i == links.P->head) {
				head_span = subtree->head_span;
				assert(head_span);
			}
			if (i == links.P->tail) {
				tail_span = subtree->tail_span;
				assert(tail_span);
			}
			final_score = final_score.add(subtree->score_combo());
			//final_score = final_score.addterm(subtree->unnorm_score());
			if (i > 0) {
				const expression_node &from_tree = children[i - 1];
				const expression_node &to_tree = subtree;
				double rscore = getrelscore(links, i, from_tree, to_tree);
				if (rscore <= SCORE_THRESHOLD) {
					return;
				}
				VERBOSE(*verb_out << "adding rel-score " << rscore << " at " << i << std::endl);
				//final_score = final_score.add(best_rel_score);
				final_score = final_score.addrel(rscore);
			}
			
			children.push_back(subtree);
		}

		size_t orig_rank = ranks[*curr];
		parse_cell *cell = links.cells[*curr];
		expression_node tree = extract_tree_internal(it, links.children[*curr].span, cell, ranks[*curr]);
		while (tree) {
			if (*curr == links.P->head) {
				head_span = tree->head_span;
				assert(head_span);
			}
			if (*curr == links.P->tail) {
				tail_span = tree->tail_span;
				assert(tail_span);
			}
			const expression_node &to_tree = tree;
			children.push_back(to_tree);
			if (*curr > 0) {
				const expression_node &from_tree = children[*curr - 1];
				double rscore = getrelscore(links, *curr, from_tree, to_tree);
				if (rscore >= SCORE_THRESHOLD) {
					score_combiner tmp_score = final_score;
					tmp_score = tmp_score.add(tree->score_combo());
					//tmp_score = tmp_score.addterm(tree->unnorm_score());
					VERBOSE(*verb_out << "adding rel-score " << rscore << " at " << *curr << std::endl);
					//tmp_score = tmp_score.add(best_rel_score);
					tmp_score = tmp_score.addrel(rscore);
					extract_multiclass_trees(it, links, tmp_score, children, ranks, first, curr + 1, last, head_span, tail_span);
				}
			}
			else {
				extract_multiclass_trees(it, links, final_score.add(to_tree->score_combo()), children, ranks, first, curr + 1, last, head_span, tail_span);
				//extract_multiclass_trees(it, links, final_score.addterm(to_tree->unnorm_score()), children, ranks, first, curr + 1, last);
			}
			children.pop_back();
			
			++ranks[*curr];
			tree = extract_tree_internal(it, links.children[*curr].span, cell, ranks[*curr]);
		}
		ranks[*curr] = orig_rank;
	}
}


static bool
fixup_tree_hack(const parse_cell *cell, expression_node node)
{
	/*
	if (node->type() == INTEGRAL_EXPR) {
		VERBOSE(*verb_out << "verify tree " << node->long_str() << " is new vs.\n");
		for (std::vector<expression_node>::const_iterator i = cell->ranked_trees.begin(); i != cell->ranked_trees.end(); ++i) {
			VERBOSE(*verb_out << "            " << (*i)->long_str() << std::endl);
			if (!strcmp((*i)->long_str(), node->long_str())) {
				return false;
			}
		}
	}
	else if (node->type() == VAR_EXPR) {
		expression_node term = node->mchild(NAME_VALUE);
		if (term->str_ == "sjn" || term->str_ == "Sin") {
			term->set_str("sin");
			term->set_long_string("sin");
			term->set_wstr(L"sin");
			node->set_long_string("<mi>sin</mi>");
			node->set_str("sin");
			node->set_wstr(L"sin");
		}
		else if (term->str_ == "Cos" || term->str_ == "cOs"
		 || term->str_ == "coS" || term->str_ == "COs"
		 || term->str_ == "CoS" || term->str_ == "cOS"
		 || term->str_ == "COS") {
			term->set_str("cos");
			term->set_long_string("cos");
			term->set_wstr(L"cos");
			node->set_long_string("<mi>cos</mi>");
			node->set_str("cos");
			node->set_wstr(L"cos");
		}
	}
	*/
	return true;
}


static bool
prepare_for_iteration(const parse_iterator &it, context *ctx, parse_cell *parent, parse_links &links, bool found_locked_nt)
{
	static std::vector<rclass_t> dummy;

	VERBOSE(
		*verb_out << "prepare_for_iteration at " << links << " with production ";
		if (links.P) *verb_out  << *links.P << std::endl;
		else *verb_out << "nil\n";
	);
	
	assert(links.ranked_trees.empty() && links.multiclass_children.empty() && links.uniclass_children.empty());

	if (!links.P) {
		assert(links.children.size() == 1);
		const nonterminal *nt = links.children.front().nt;
		assert(nt->name.length() > 2);
		assert(nt->name[0] == '_');
		assert(nt->name[1] == 'S');
		//VERBOSE(*verb_out << "creating terminal node for " << links.children.front().nt->name << " at " << links.score << std::endl);
		VERBOSE(*verb_out << "creating terminal node for " << links.children.front().nt->name << " at " << links.terminal_score << std::endl);
		links.cells.push_back(ctx->table[links.children[0].span][links.children[0].nt]);
		expression_node tree = make_terminal(links.children.front().nt->name.substr(2), links.terminal_score);
		//tree->set_score(score_combiner(links.score->lookup(tree->children, links, links.multiclass_children.begin(), links.multiclass_children.end())));
		//VERBOSE(*verb_out << " set score to " << tree->score() << std::endl);

		tree->rclasses.insert(nt->rclasses.begin(), nt->rclasses.end());

		VERBOSE(
			*verb_out << " node has classes ";
			for (std::set<rclass_t>::const_iterator i = tree->rclasses.begin(); i != tree->rclasses.end(); ++i) {
				*verb_out << *i << ' ';
			}
		);


		tree->set_span(links.children[0].span);
		tree->head_span = tree->tail_span = links.children[0].span;
		assert(tree->head_span);
		assert(tree->tail_span);
		VERBOSE(*verb_out << "setting head span to " << tree->head_span->bits << "; tail span to " << tree->tail_span->bits << std::endl);
		tree->set_context(ctx);
		tree->set_nt(links.children[0].nt);
		tree->nt_low = links.children[0].nt;
		tree->prod = 0;

		VERBOSE(*verb_out << "setting low NT to child nt " << links.children[0].nt->name << std::endl);
		links.ranked_trees.push_back(tree);
		return true;
	}


	links.cells.reserve(links.children.size());
	for (std::vector<parseref>::iterator i = links.children.begin(); i != links.children.end(); ++i) {
		links.cells.push_back(ctx->table[i->span][i->nt]);
		links.cells.back()->parents.push_back(parent);
	}

	std::vector<parse_cell *>::iterator pc = links.cells.begin();
	for (std::vector<parseref>::iterator i = links.children.begin(); i < links.children.end(); ++i, ++pc) {
		if ((*pc)->classes.size() > 1) {
			links.multiclass_children.push_back(i - links.children.begin());
		}
		else {
			links.uniclass_children.push_back(i - links.children.begin());
		}
	}

	VERBOSE(
		*verb_out << "multiclass children = [ ";
		for (std::vector<size_t>::const_iterator i = links.multiclass_children.begin(); i != links.multiclass_children.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "]\nuniclass children = [ ";
		for (std::vector<size_t>::const_iterator i = links.uniclass_children.begin(); i != links.uniclass_children.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "]\n";
	);

	score_combiner base_score;
	std::vector<parseref>::iterator first_mc;
	if (links.multiclass_children.empty()) {
		first_mc = links.children.end() - 1;
	}
	else {
		first_mc = links.children.begin() + (links.multiclass_children.front() - 1);
	}

	/*VERBOSE(*verb_out << "building scoring tree for links " << links << std::endl);
	if (links.children.size() == 1) {
		links.score = new t_relscore(score_combiner(-1.0));
	}
	else {
		build_relscore_tree(links.score, score_combiner(), links, links.multiclass_children.begin(), links.multiclass_children.begin(), links.multiclass_children.end(), INVALID_CLASS);
	}

	if (!links.score) {
		VERBOSE(*verb_out << " nothing had a high enough rel score; there will be no results here\n");
		return false;
	}*/

	pc = links.cells.begin();
	for (std::vector<parseref>::iterator i = links.children.begin(); i != links.children.end(); ++i, ++pc) {
		if (!prepare_for_iteration(it, ctx, i->span, *pc, found_locked_nt)) {
			return false;
		}
	}

	std::vector<size_t> ranks(links.children.size(), 0);
	std::vector<expression_node> children;
	extract_multiclass_trees(it, links, score_combiner(), children, ranks, links.multiclass_children.begin(), links.multiclass_children.begin(), links.multiclass_children.end(), 0, 0);

	VERBOSE(*verb_out << "after extracting trees at " << links << " there are " << links.cached_trees.size() << " options\n");

	if (links.cached_trees.empty()) {
		return false;
	}

	const links_parse_result &best = *links.cached_trees.begin();
	VERBOSE(*verb_out << "for links " << links << " the top-ranked tree is " << best.tree->long_str() << " at " << best.tree->score() << " with type " << best.tree->type() << std::endl);
	VERBOSE(*verb_out << *best.tree << std::endl);
	links.ranked_trees.push_back(best.tree);
	//cell->cached_trees.pop();
	links.cached_trees.erase(links.cached_trees.begin());


	//for (std::vector<size_t>::const_iterator j = links.uniclass_children.begin(); j != links.uniclass_children.end(); ++j) {
	for (size_t j = 0; j < links.children.size(); ++j) {
		++ranks[j];
		links.next_index_sets.push_back(ranks);
		--ranks[j];
	}

	return true;
}

const size_t parse_cell::INVALID_SRC = static_cast<size_t>(-1);


bool
prepare_for_iteration(const parse_iterator &it, context *ctx, ordered_segments *span, parse_cell *cell, bool found_locked_nt)
{
	VERBOSE(*verb_out << "prepare_for_iteration(" << cell << ")\n");

	if (cell->links.empty() && cell->nt->name[0] != '!') {
		return false;
	}

	if (span->locked_node) {
		//cell->ranked_trees.push_back(span->locked_node);
		return true;
	}
	else if (!found_locked_nt && span->locked_nt && cell->nt->sid != InvalidSemanticId && cell->nt != span->locked_nt && !cell->nt->has_direct_descendent(span->locked_nt)) {
		VERBOSE(*verb_out << " NT " << cell->nt->name << " does not match locked NT " << span->locked_nt->name << std::endl);
		return false;
	}

	if (span->locked_nt) {
		VERBOSE(*verb_out << " span " << span->bits << " locked to NT " << span->locked_nt->name << std::endl);
	}
	
	if (cell->prev_src != parse_cell::INVALID_SRC) {
		if (it.flags == cell->prepared_flags) {
			VERBOSE(*verb_out << "this cell is already initialized with " << cell->ranked_trees.size() << " ranked trees\n");
			return !cell->ranked_trees.empty();
		}
		invalidate_single_cell(cell);
		VERBOSE(*verb_out << "INVALIDATED CELL! flags " << it.flags << " vs. " << cell->prepared_flags << std::endl);
	}
	
	cell->prepared_flags = it.flags;
	VERBOSE(*verb_out << "preparing flags of " << cell->nt->name << ':' << cell->segs->bits << " to " << cell->prepared_flags << std::endl);

	cell->prev_src = 0;

	VERBOSE(*verb_out << "prepare_for_iteration at cell " << cell->nt->name << ':' << cell->segs->bits << " with " << cell->links.size() << " links\n" << std::flush);

	
	if (cell->nt->name[0] == '!') { // indicates matrix
		if (it.flags & parse_iterator::LONG_FORM_ENABLED) {
			cell->mxiter = cell->matrix->createLongFormIterator();
		}
		else {
			cell->mxiter = cell->matrix->createIterator();
		}
		expression_node tree = extract_matrix_tree(ctx, cell, span);
		if (tree) {
			tree->set_context(ctx);
			tree->set_cell(cell);
			tree->set_span(span);
			VERBOSE(*verb_out << "matrix tree is " << tree.ptr() << " or " << tree->long_str() << " at " << tree->score() << " with bounds " << tree->span->bounds() << std::endl);
			cell->ranked_trees.push_back(tree);
			return true;
		}
		return false;
	}
	else {
		for (std::vector<parse_links>::iterator i = cell->links.begin(); i != cell->links.end(); ++i) {
			VERBOSE(*verb_out << " considering links " << *i << std::endl);
			if (span->locked_nt) {
				if (i->P->nt != span->locked_nt && (i->P->rhs.size() > 1 || (i->P->rhs[0] != span->locked_nt && !i->P->rhs[0]->has_direct_descendent(span->locked_nt)))) {
					VERBOSE(*verb_out << "  (can't derive only locked NT here\n");
					continue;
				}
				if (i->P->nt == span->locked_nt && !span->locked_parse.empty() && span->locked_parse != i->children) {
					VERBOSE(*verb_out << "   (doesn't match locked parse...skipping)\n");
					continue;
				}
			}
			cell->next_best_ranks.push_back(0);
			if (found_locked_nt || !span->locked_nt || (i->P && (span->locked_nt == i->P->nt || i->P->nt->has_direct_descendent(span->locked_nt)))) {
				if (prepare_for_iteration(it, ctx, cell, *i, found_locked_nt || (i->P && span->locked_nt == i->P->nt))) {
					expression_node tree = extract_tree(it, *i, 0);
					if (tree) {
						tree->set_context(ctx);
						tree->set_cell(cell);
						tree->set_span(span);
						score_combiner sc = tree->score_combo();
						if (span->locked_nt) {
							sc.bias(100.0);
						}
						tree->set_score(sc);
						VERBOSE(*verb_out << "at " << *i << " with base index set, tree is " << tree.ptr() << " or " << tree->long_str() << " at " << tree->score() << " with bounds " << tree->span->bounds() << std::endl);
						//cell->cached_trees.push(cell_parse_result(tree, i - cell->links.begin()));
						cell->cached_trees.insert(cell_parse_result(tree, i - cell->links.begin()));
					}
				}
			}
		}
	}

	if (cell->cached_trees.empty()) {
		return false;
	}

	//const cell_parse_result &best = cell->cached_trees.top();
	cell_parse_result best = *cell->cached_trees.begin();
	
	fixup_tree_hack(cell, best.tree);
	
	VERBOSE(*verb_out << "for cell " << cell->nt->name << ':' << cell->segs->bits << " the top-ranked tree is " << best.tree->long_str() << " at " << best.tree->score() << " with type " << best.tree->type() << std::endl);
	cell->ranked_trees.push_back(best.tree);
	++cell->next_best_ranks[best.src];
	cell->prev_src = best.src;
	//cell->cached_trees.pop();
	cell->cached_trees.erase(cell->cached_trees.begin());

	return true;
}


parse_cell *
prepare_for_iteration(const parse_iterator &it, context *ctx, ordered_segments *segs, const nonterminal *nt)
{
	parse_table::iterator i = ctx->table.find(segs);
	if (i == ctx->table.end()) {
		return 0;
	}
	nonterminal_parses::iterator j = i->second.find(nt);
	if (j == i->second.end()) {
		return 0;
	}

	parse_cell *cell = j->second;
	return prepare_for_iteration(it, ctx, segs, cell, false) ? cell : 0;
}


static expression_node
extract_tree(const parse_iterator &it, parse_links &links, size_t rank)
{
	if (rank < links.ranked_trees.size()) {
		return links.ranked_trees[rank];
	}

	while (rank > links.ranked_trees.size()) {
		if (!extract_tree(it, links, links.ranked_trees.size())) {
			return NULL_EXPRESSION_NODE;
		}
	}
	assert(rank == links.ranked_trees.size());

	VERBOSE(
		*verb_out << "extract_tree(" << rank << ") at " << links << " with " << links.cached_trees.size() << " reserved trees:" << std::endl;
		for (links_expression_heap::const_iterator i = links.cached_trees.begin(); i != links.cached_trees.end(); ++i) {
			*verb_out << "  " << i->tree->score() << " : " << i->tree->long_str() << std::endl;
		}
		
		*verb_out << "next index options:\n";
		for (std::list<std::vector<size_t> >::const_iterator i = links.next_index_sets.begin(); i != links.next_index_sets.end(); ++i) {
			*verb_out << " [ ";
			for (std::vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
				*verb_out << *j << ' ';
			}
			*verb_out << "]\n";
		}
	);


	size_t last_child;
	if (links.multiclass_children.empty()) {
		last_child = links.children.size();
	}
	else {
		last_child = links.multiclass_children[0];
	}


	for (std::list<std::vector<size_t> >::iterator i = links.next_index_sets.begin(); i != links.next_index_sets.end(); ++i) {
		std::vector<size_t> &ranks = *i;
		std::vector<expression_node> children;
		extract_multiclass_trees(it, links, score_combiner(), children, ranks, links.multiclass_children.begin(), links.multiclass_children.begin(), links.multiclass_children.end(), 0, 0);
		links.expired_index_sets.insert(*i);
	}

	links.next_index_sets.clear();

	if (links.cached_trees.empty()) {
		VERBOSE(*verb_out << "there are NO TREES at " << links << "; returning 0\n");
		return NULL_EXPRESSION_NODE;
	}

	//const links_parse_result &best = links.cached_trees.top();
	const links_parse_result &best = *links.cached_trees.begin();
	links.ranked_trees.push_back(best.tree);

	std::vector<size_t> ranks = best.src;
	for (std::vector<size_t>::const_iterator i = links.uniclass_children.begin(); i != links.uniclass_children.end(); ++i) {
		++ranks[*i];
		if (links.expired_index_sets.find(ranks) == links.expired_index_sets.end()) {
			links.next_index_sets.push_back(ranks);
		}
		--ranks[*i];
	}

	//links.cached_trees.pop();
	links.cached_trees.erase(links.cached_trees.begin());

	VERBOSE(*verb_out << "extract_tree(" << rank << ") at " << links << " gives " << links.ranked_trees[rank]->long_str() << " at " << links.ranked_trees[rank]->score() << std::endl);

	return links.ranked_trees[rank];
}


parse_cell::parse_cell(ordered_segments *segs_, const nonterminal *nt_) : segs(segs_), nt(nt_), prev_src(parse_cell::INVALID_SRC), prepared_flags(-1), matrix(0), mxiter(0) { }


expression_node
extract_tree_internal(const parse_iterator &it, ordered_segments *span, parse_cell *cell, size_t rank)
{
	VERBOSE(*verb_out << "extract tree ranked " << rank << " at " << cell->nt->name << ":" << span->bits << " with " << cell->ranked_trees.size() << " already found" << std::endl);

	if (rank > 0) {
		if (it.flags & parse_iterator::ONLY_SEMANTIC_ALTERNATES && cell->nt->sid != InvalidSemanticId) {
			return NULL_EXPRESSION_NODE;
		}
		if (it.flags & parse_iterator::ONLY_PRIMARY_SPAN_ALTERNATES && span != it.start_span) {
			return NULL_EXPRESSION_NODE;
		}
	}

	if (rank == 0) {
		if (span->locked_node) {
			if (cell->nt != span->locked_node->nt_low && !cell->nt->has_direct_descendent(span->locked_node->nt_low)) {
				return NULL_EXPRESSION_NODE;
			}
			return span->locked_node;
		}
		if (span->locked_nt && cell->nt != span->locked_nt && !cell->nt->has_direct_descendent(span->locked_nt)) {
			return NULL_EXPRESSION_NODE;
		}
		//else if (span->locked_nt && !cell->nt->has_direct_descendent(span->locked_nt)) {//cell->nt != span->locked_nt && cell->nt->sid != InvalidSemanticId) {
		//	return NULL_EXPRESSION_NODE;
		//}
	}
	
	if (rank < cell->ranked_trees.size()) {
		return cell->ranked_trees[rank];
	}

	if (span->locked_node) {
		return NULL_EXPRESSION_NODE;
	}


	while (rank > cell->ranked_trees.size()) {
		if (!extract_tree(it, span, cell, cell->ranked_trees.size())) {
			return NULL_EXPRESSION_NODE;
		}
	}
	assert(rank == cell->ranked_trees.size());

	VERBOSE(
		*verb_out << "extract_tree(" << rank << ") at cell " << cell << " with " << cell->cached_trees.size() << " reserved trees" << std::endl;
		for (cell_expression_heap::const_iterator i = cell->cached_trees.begin(); i != cell->cached_trees.end(); ++i) {
			*verb_out << "  " << i->tree->score() << " : " << i->tree->long_str() << std::endl;
		}
		*verb_out << "options:\n";
		for (size_t i = 0; i < cell->links.size(); ++i) {
			*verb_out << "  rank " << cell->next_best_ranks[i] << " at " << cell->links[i] << std::endl;
		}
	);

	expression_node next_tree;

	
	if (cell->nt->name[0] == '!') { // indicates matrix
		next_tree = extract_matrix_tree(cell->ranked_trees[0]->ctx, cell, span);
	}
	else {
		next_tree = extract_tree(it, cell->links[cell->prev_src], cell->next_best_ranks[cell->prev_src]);
	}

	if (next_tree) {
		VERBOSE(*verb_out << "caching tree " << *next_tree << " in cell " << cell->nt->name << ':' << cell->segs->bits << std::endl);
		next_tree->set_context(cell->ranked_trees[0]->ctx);
		next_tree->set_cell(cell);
		next_tree->set_span(span);
		//cell->cached_trees.push(cell_parse_result(next_tree, cell->prev_src));
		cell_expression_heap::const_iterator qq = cell->cached_trees.find(cell_parse_result(next_tree, cell->prev_src));
		if (qq != cell->cached_trees.end()) {
			VERBOSE(*verb_out << " THIS WILL OVERWRITE " << qq->tree->long_str() << std::endl);
		}
		cell->cached_trees.insert(cell_parse_result(next_tree, cell->prev_src));
	}

	if (cell->cached_trees.empty()) {
		VERBOSE(*verb_out << "there are NO TREES at " << cell->nt->name << ':' << span->bits << "; returning 0\n");
		return NULL_EXPRESSION_NODE;
	}

	//const cell_parse_result &best = cell->cached_trees.top();
	cell_parse_result best = *cell->cached_trees.begin();

	bool have_usable_tree = fixup_tree_hack(cell, best.tree);

	++cell->next_best_ranks[best.src];
	cell->prev_src = best.src;

	if (have_usable_tree) {
		cell->ranked_trees.push_back(best.tree);
		//cell->cached_trees.pop();
		cell->cached_trees.erase(cell->cached_trees.begin());
		VERBOSE(*verb_out << "extract_tree(" << rank << ") at cell " << cell << " gives " << cell->ranked_trees[rank]->long_str() << " at " << cell->ranked_trees[rank]->score() << std::endl);
		return cell->ranked_trees[rank];
	}
	else {
		//cell->cached_trees.pop();
		cell->cached_trees.erase(cell->cached_trees.begin());
		return extract_tree_internal(it, span, cell, rank);
	}
}


/*
static expression_node
collapse_tree(expression_node node, const nonterminal *rootnt)
{
	if (node->nchildren() == 0) {
		raw_expression_node *copy = new raw_expression_node(*node);
		copy->nt_low = node->nt;
		return expression_node(copy);
	}

	if (node->type() == InvalidSemanticId) {
		assert(node->nchildren() == 1);
		return collapse_tree(node->mchild(0), rootnt);
	}

	expression_node newnode(new raw_expression_node);
	newnode->parent = node->parent;
	newnode->set_type(node->type());
	newnode->set_score(node->score_combo());
	newnode->set_long_string(node->long_str());
	newnode->set_span(node->span);
	newnode->set_cell(node->cell);
	newnode->set_context(node->ctx);
	newnode->set_nt(rootnt);
	newnode->nt_low = node->nt;
	if (node->nchildren() == 1) {
		newnode->add_child(collapse_tree(node->mchild(0), rootnt).ptr());
	}
	else {
		for (size_t i = 0; i < node->nchildren(); ++i) {
			newnode->add_child(collapse_tree(node->mchild(i), node->mchild(i)->nt).ptr());
		}
	}
	return newnode;
}
*/

expression_node
extract_tree(const parse_iterator &it, ordered_segments *span, parse_cell *cell, size_t rank, bool wrap_string)
{
	expression_node tree = extract_tree_internal(it, span, cell, rank);
	VERBOSE(*verb_out << "EXTRACTED tree " << (tree ? tree->long_str() : "null") << std::endl);
	if (!tree) {
		return NULL_EXPRESSION_NODE;
	}
	expression_node node(new raw_expression_node(*tree));
	return node;//tree ? collapse_tree(tree, tree->nt) : 0;
}


void
dump_span_table(std::ostream &os, const context *ctx)
{
	for (std::map<bitvec, ordered_segments *>::const_iterator i = ctx->spans->begin(); i != ctx->spans->end(); ++i) {
		os << i->first << " -> " << i->second << std::endl;
	}
}


void
dump_parse_table(std::ostream &os, const context *ctx)
{
	for (parse_table::const_iterator i = ctx->table.begin(); i != ctx->table.end(); ++i) {
		for (nonterminal_parses::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
			os << j->first->name << ':' << i->first->bits << "(" << j->first << " : " << i->first << ") - bounds " << i->first->bounds() << ", ptr " << j->second << std::endl;
			if (j->second) {
				os << " " << j->second->classes.size() << " rclasses\n";
				for (std::vector<parse_links>::const_iterator k = j->second->links.begin(); k != j->second->links.end(); ++k) {
					os << "  " << *k << std::endl;
				}
			}
			else {
				os << "  nil\n";
			}
		}
	}
}


static void
invalidate_single_cell(parse_cell *cell)
{
	VERBOSE(
		*verb_out << "invalidating cell " << cell->nt->name << ":" << cell->segs->bits << " (" << cell << ") with links:\n";
		for (std::vector<parse_links>::iterator j = cell->links.begin(); j != cell->links.end(); ++j) {
			*verb_out << " " << *j << std::endl;
		}
	);
	for (std::vector<parse_links>::iterator j = cell->links.begin(); j != cell->links.end(); ++j) {
		parse_links &links = *j;
		links.next_index_sets.clear();
		links.expired_index_sets.clear();
		links.ranked_trees.clear();
		links.cells.clear();
		links.multiclass_children.clear();
		links.uniclass_children.clear();
		while (!links.cached_trees.empty()) {
			//links.cached_trees.pop();
			links.cached_trees.erase(links.cached_trees.begin());
		}
	}
	
	cell->prev_src = parse_cell::INVALID_SRC;

	cell->ranked_trees.clear();
	while (!cell->cached_trees.empty()) {
		//cell->cached_trees.pop();
		cell->cached_trees.erase(cell->cached_trees.begin());
	}
	cell->next_best_ranks.clear();
}


static void
invalidate_cell(parse_cell *cell)
{
	if (cell->prev_src == parse_cell::INVALID_SRC) {
		return;
	}

	invalidate_single_cell(cell);
	
	for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		invalidate_cell(*i);
	}
	cell->parents.clear();
}


static void
lock_existing_semantics(context *ctx, parse_cell *cell, bool lock)
{
	if (cell->nt->sid != InvalidSemanticId) {
		lock = false;
	}
	
	if (lock) {
		VERBOSE(*verb_out << "locking existing semantics of " << cell->segs->bits << " to " << cell->nt->name << std::endl);
		cell->segs->locked_nt = cell->nt;
	}

	nonterminal_parses &parses = ctx->table[cell->segs];
	for (nonterminal_parses::iterator i = parses.begin(); i != parses.end(); ++i) {
		if (i->second) {
			invalidate_cell(i->second);
		}
	}	

	for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		lock_existing_semantics(ctx, *i, lock);
	}
}

static void
invalidate_table_to_top(context *ctx, parse_cell *cell, std::set<ordered_segments *> &cleared) {
	for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		/*if (cleared.find((*i)->segs) == cleared.end()) {
			cleared.insert((*i)->segs);
			nonterminal_parses &parses = ctx->table[(*i)->segs];
			for (nonterminal_parses::iterator j = parses.begin(); j != parses.end(); ++j) {
				if (j->second) {
					invalidate_table_to_top(ctx, j->second, cleared);
					invalidate_single_cell(j->second);
				}
			}
		}
		if (cell->parents.empty()) {
			assert(false);
		}*/
		invalidate_table_to_top(ctx, *i, cleared);
		invalidate_single_cell(*i);
	}
	cell->parents.clear();
}

static void
unlock_recursive(parse_cell *cell, std::set<parse_cell *> &seen)
{
	if (seen.find(cell) == seen.end()) {
		seen.insert(cell);
		cell->segs->locked_nt = 0;

		for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
			unlock_recursive(*i, seen);
		}
	}
}


bool
lock_semantics(context *ctx, ordered_segments *span, const raw_expression_node *node)
{
	if (span->locked_node) {
		return false;
	}

	VERBOSE(*verb_out << "locking semantics of " << span->bits << " to " << node->nt_low->name << std::endl);
	span->locked_nt = node->nt_low;
	
	span->locked_parse = node->parse_links;
	VERBOSE(
		*verb_out << "  locked parse is ";
		for (std::vector<parseref>::const_iterator i = span->locked_parse.begin(); i != span->locked_parse.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << std::endl;
	);

	/*if (IsAutoTrainingEnabled() && node->prod && node->prod->rel) {
		VERBOSE(*verb_out << "autotraining for links " << node->parse_links << std::endl);
		for (size_t i = 0; i < node->parse_links.size() - 1; ++i) {
			const parseref &crpref = node->parse_links[i];
			const parseref &nxpref = node->parse_links[i+1];
			const parse_cell *fromcell = ctx->table[crpref.span][crpref.nt];
			const parse_cell *tocell = ctx->table[nxpref.span][nxpref.nt];
			expression_node fromtree = fromcell->ranked_trees[node->ranks[i]];
			expression_node totree = tocell->ranked_trees[node->ranks[i+1]];
			for (std::set<rclass_t>::const_iterator j = fromtree->rclasses.begin(); j != fromtree->rclasses.end(); ++j) {
				for (std::set<rclass_t>::const_iterator k = totree->rclasses.begin(); k != totree->rclasses.end(); ++k) {
					VERBOSE(*verb_out << "autotrain: adding sample of link " << crpref.nt->name << ':' << crpref.span->bits << " -> " << nxpref.nt->name << ':' << nxpref.span->bits << " with classes " << *j << " and " << *k << " to relation " << node->prod->rel->name << std::endl);
					node->prod->rel->add_sample(crpref.span, nxpref.span, *j, *k); 
				}
			}
		}
	}*/

	parse_cell *cell = ctx->table[span][span->locked_nt];
	assert(cell);
	std::set<ordered_segments *> cleared;
	invalidate_table_to_top(ctx, cell, cleared);
	/*for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		//lock_existing_semantics(ctx, *i, true);
		invalidate_table_to_top(ctx, *i);
	}*/

	/*nonterminal_parses &parses = ctx->table[span];
	for (nonterminal_parses::iterator i = parses.begin(); i != parses.end(); ++i) {
		if (i->second) {
			invalidate_cell(i->second);
		}
	}*/

	ctx->locks->insert(span);
	//invalidate_cell(cell);

	return true;
}



bool
lock_expression(context *ctx, ordered_segments *span, const raw_expression_node *node)
{
	if (span->locked_nt) {
		return false;
	}

	if (span->locked_node) {
		span->locked_node = NULL_EXPRESSION_NODE;
	}

	VERBOSE(*verb_out << "locking tree at " << span->bits << " to " << node->long_str() << std::endl);

	expression_node locked_node(new raw_expression_node(*node));
	locked_node->set_score(score_combiner(100.0));
	span->locked_node = locked_node;
	//span->locked_parse = node->parse_links;

	parse_cell *cell = ctx->table[span][node->nt_low];
	assert(cell);
	std::set<ordered_segments *> cleared;
	invalidate_table_to_top(ctx, cell, cleared);
	/*for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		//lock_existing_semantics(ctx, *i, true);
		invalidate_table_to_top(ctx, *i);
	}*/

	if (IsAutoTrainingEnabled() && derives_only_terminal(node)) {
		const ExpressionTree *termnode = node;
		while (termnode->type() != TERMINAL_EXPR) {
			termnode = termnode->child(0);
		}
		std::string name = termnode->str();
		VERBOSE(*verb_out << "autotrain: adding terminal example for " << name << std::endl);
		symbols_db *autodb = global_autotrain_db();
		if (autodb) {
			symbol *S = autodb->find_symbol(name);
			if (S) {
				prototype *P = new prototype;
				P->id.pkgid = AUTOTRAIN_PKGID;
				P->id.unicode = S->info.unicode;
				P->id.seq = grab_next_prototype_seq(*autodb, P->id);
				RawStroke *stks = new RawStroke[span->size()];
				RawStroke *pstk = stks;
				for (ordered_segments::ordered_iterator i = span->begin(0); i != span->end(0); ++i, ++pstk) {
					const segment *seg = *i;
					*pstk = copy(seg->stk);
				}
				P->raw_strokes.set_strokes(stks, span->size());
				VERBOSE(*verb_out << "autotrain: example has " << span->size() << " strokes; adding to autotrain db\n");
				S->add_prototype(P);
				S->update_active_table(P);

				symbols_db *db = global_symbols_db_wr();
				S = db->find_symbol(name);
				if (S) {
					prototype *newP = new prototype;
					newP->id = P->id;
					newP->raw_strokes = copy(P->raw_strokes);
					VERBOSE(*verb_out << "autotrain: and adding to current db\n");
					S->add_prototype(newP);
					S->update_active_table(newP);
				}
				// TODO: should the correlation data be updated here as well?
			}
		}
	}
	
	/*nonterminal_parses &parses = ctx->table[span];
	for (nonterminal_parses::iterator i = parses.begin(); i != parses.end(); ++i) {
		if (i->second) {
			invalidate_cell(i->second);
		}
	}*/
	
	//invalidate_cell(cell);
	ctx->locks->insert(span);
	return true;
}


/*
bool
lock_children(context *ctx, ordered_segments *span, const raw_expression_node *parent)
{
	if (span->locked_node) {
		return false;
	}

	if (!span->locked_parse.empty()) {
		span->locked_parse.clear();
	}

	VERBOSE(*verb_out << "locking child spans at " << span->bits << " to ");
	for (size_t i = 0; i < parent->nchildren(); ++i) {
		const expression_node &child = parent->mchild(i);
		if (child->span) {
			span->locked_parse.push_back(child->span->bits);
		}
		else {
			span->locked_parse.push_back(bitvec(ctx->segments.size(), false));
		}
		VERBOSE(*verb_out << span->locked_parse.back() << ' ');
	}
	VERBOSE(*verb_out << std::endl);

	parse_cell *cell = ctx->table[span][parent->nt_low];
	assert(cell);
	/*for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		lock_existing_semantics(*i);
	}* /

	nonterminal_parses &parses = ctx->table[span];
	for (nonterminal_parses::iterator i = parses.begin(); i != parses.end(); ++i) {
		if (i->second) {
			invalidate_cell(i->second);
		}
	}

	//invalidate_cell(cell);

	return true;
}
*/


int
is_locked(ordered_segments *span)
{
	if (span->locked_node) {
		return SPAN_LOCKED_NODE;
	}

	if (span->locked_nt) {
		return SPAN_LOCKED_TYPE;
	}
	
	if (!span->locked_parse.empty()) {
		return SPAN_LOCKED_CHILDREN;
	}

	return SPAN_UNLOCKED;
}

void
unlock_all(context *ctx)
{
	while (!ctx->locks->empty()) {
		ordered_segments *span = *ctx->locks->begin();
		unlock(ctx, span);
	}
}

bool
unlock(context *ctx, ordered_segments *span)
{
	const nonterminal *nt = 0;
	if (span->locked_node) {
		assert(!span->locked_nt);
		nt = span->locked_node->nt;
		span->locked_node = NULL_EXPRESSION_NODE;
	}
	else if (span->locked_nt) {
		nt = span->locked_nt;
		span->locked_nt = 0;
	}
	else {
		VERBOSE(*verb_out << "attempted to unlock unlocked span " << span->bits << std::endl);
		return true;
	}
	
	if (!span->locked_parse.empty()) {
		span->locked_parse.clear();
	}

	nonterminal_parses &parses = ctx->table[span];
	for (nonterminal_parses::iterator i = parses.begin(); i != parses.end(); ++i) {
		if (i->second) {
			VERBOSE(*verb_out << "invalidating cell " << i->second << " for nt " << i->first->name << std::endl);
			invalidate_cell(i->second);
		}
	}

	std::set<ordered_segments *>::iterator plock = ctx->locks->find(span);
	assert(plock != ctx->locks->end());
	ctx->locks->erase(plock);

	if (nt) {
		parse_cell *cell = ctx->table[span][nt];
		assert(cell);
		std::set<parse_cell *> seen;
		unlock_recursive(cell, seen);
		invalidate_cell(cell);
	}
	
	return true;
}


static void
clear_parent_refs(parse_cell *cell)
{
	VERBOSE(*verb_out << "clearing parent refs for " << cell << " @ " << cell->nt->name << ":" << cell->segs->bits << std::endl);
	for (std::vector<parse_links>::iterator i = cell->links.begin(); i != cell->links.end(); ++i) {
		for (std::vector<parse_cell *>::iterator j = i->cells.begin(); j != i->cells.end(); ++j) {
			parse_cell *child = *j;
			if (!child) continue;
			VERBOSE(*verb_out << "child is " << child << std::endl);
			std::list<parse_cell *>::iterator parentp = std::find(child->parents.begin(), child->parents.end(), cell);
			if (parentp != child->parents.end()) {
				child->parents.erase(parentp);
			}
		}
	}
}

static void
clear_parent_refs(nonterminal_parses &tab)
{
	for (nonterminal_parses::iterator j = tab.begin(); j != tab.end(); ++j) {
		parse_cell *cell = j->second;
		if (cell) {
			clear_parent_refs(cell);
		}
	}
}

static void
clear_table(nonterminal_parses &tab)
{
	for (nonterminal_parses::iterator j = tab.begin(); j != tab.end(); ++j) {
		parse_cell *cell = j->second;
		if (cell) {
			VERBOSE(*verb_out << "deleting cell " << cell << " at " << cell->nt->name << ":" << cell->segs->bits << std::endl);
			delete cell;
		}
	}
	tab.clear();
}


void
clear_context(context *ctx)
{
	while (!ctx->segments.empty()) {
		remove_segment(ctx, ctx->segments.back());
		//ctx->segments.pop_back();
	}

	assert(ctx->strokes.empty());

	if (ctx->groups) {
		assert(ctx->groups->empty());
		delete ctx->groups;
		ctx->groups = 0;
	}

	if (ctx->spans) {
		for (std::map<bitvec, ordered_segments *>::iterator i = ctx->spans->begin(); i != ctx->spans->end(); ++i) {
			delete i->second;
		}

		delete ctx->spans;
		ctx->spans = 0;
	}
	
	if (ctx->locks) {
		delete ctx->locks;
		ctx->locks = 0;
	}

	for (parse_table::iterator i = ctx->table.begin(); i != ctx->table.end(); ++i) {
		clear_table(i->second);
	}

	ctx->table.clear();
}


int
initialize_context(context *ctx, const grammar &G)
{
	ctx->spans = new std::map<bitvec, ordered_segments *>;
	ctx->groups = new std::map<bitvec, group *>;
	ctx->locks = new std::set<ordered_segments *>;

	ctx->G = G;

	return 0;
}


int
update_spans_on_removal(context *ctx, size_t pos)
{
	std::map<bitvec, ordered_segments *> *newspans = new std::map<bitvec, ordered_segments *>;

	for (parse_table::iterator i = ctx->table.begin(); i != ctx->table.end(); ++i) {
		nonterminal_parses &nt_table = i->second;
		for (nonterminal_parses::iterator j = nt_table.begin(); j != nt_table.end(); ++j) {
			parse_cell *cell = j->second;
			if (!cell) continue;
			for (std::list<parse_cell *>::iterator k = cell->parents.begin(); k != cell->parents.end(); ) {
				const parse_cell *parent = *k;
				const bitvec &bits = parent->segs->bits;
				if (pos < bits.size() && bits.at(pos)) {
					//VERBOSE(*verb_out << "parent " << parent->nt->name << ":" << bits << " of " << cell->nt->name << ":" << cell->segs->bits << " affected by removal at " << pos << std::endl);
					k = cell->parents.erase(k);
				}
				else {
					++k;
				}
			}
		}
	}

	
	for (std::map<bitvec, ordered_segments *>::const_iterator i = ctx->spans->begin(); i != ctx->spans->end(); ++i) {
		ordered_segments *span = i->second;
		if (span) {
			bitvec &bits = span->bits;
			if (pos < bits.size()) {
				if (!bits.at(pos)) {
					bits.erase(pos);
					(*newspans)[bits] = span;
				}
				else {
					parse_table::iterator j = ctx->table.find(span);
					if (j != ctx->table.end()) {
						//clear_parent_refs(j->second);
						clear_table(j->second);
						ctx->table.erase(j);
					}
					std::set<ordered_segments *>::iterator k = ctx->locks->find(span);
					if (k != ctx->locks->end()) {
						ctx->locks->erase(k);
					}
				}
			}
			else {
				(*newspans)[bits] = span;
			}
		}
	}

	/*for (std::map<bitvec, ordered_segments *>::const_iterator i = newspans->begin(); i != newspans->end(); ++i) {
		ordered_segments *span = i->second;
		assert(span && span->bits.size() <= ctx->segments.size());
	}*/
	
	delete ctx->spans;
	ctx->spans = newspans;
		
	VERBOSE(
		*verb_out << "after span removal, span table is...\n";
		dump_span_table(*verb_out, ctx);
	);

	return 0;
}

/*
void
append_edit_op(context *ctx, int type, size_t pos)
{
	ctx->edits.push_front(edit_op(type, pos));
}
*/

#endif


static ordered_segments *
getsegs(const bitvec &bits) {
	std::map<bitvec, ordered_segments *>::iterator i = osegs.find(bits);
	return i == osegs.end() ? 0 : i->second;
}

static ordered_segments *
getsegs(const bitvec &bits, const ordered_segments *super, size_t d, size_t start, size_t end) {
	ordered_segments *&segs = osegs[bits];
	if (!segs) {
		segs = super->slice(d, start, end);
	}
	return segs;
}

static interpreter<raw_expression_node> *mkparser(const context *ctx, const nonterminal *nt, const ordered_segments *segs);

static bool
contains(const Rect<long> &super, const Rect<long> &sub) {
	return super.left <= sub.left && super.right >= sub.right && super.top <= sub.top && super.bottom >= sub.bottom;
}

template <typename T>
static T *
getnth(interpreter<T> *intrpr, size_t n) {
	size_t nk = intrpr->nknown();
	T *t;
	while (intrpr->nknown() <= n) {
		t = intrpr->next();
		if (!t) {
			return 0;
		}
	}
	return intrpr->nth(n);
}


class treeparser : public parser<raw_expression_node> {
private:
	struct basic_tree_accessor : public subtree_accessor {
		explicit basic_tree_accessor(std::vector<interpreter<raw_expression_node> *> &parsers_)
			: parsers(parsers_) {
		}

		expression_node operator[](size_t i) const {
			return expression_node(parsers[i]->next());
		}

		bool verify(size_t i, const expression_node &node) const {
			return true;
		}
	private:
		std::vector<interpreter<raw_expression_node> *> &parsers;
	};

	struct ranked_tree_accessor : public subtree_accessor {
		ranked_tree_accessor(std::vector<interpreter<raw_expression_node> *> &parsers_, const std::vector<size_t> &ranks_)
			: parsers(parsers_), ranks(ranks_) {
			std::cerr << "ranked_tree_accessor([";
			for (std::vector<size_t>::const_iterator i = ranks.begin(); i != ranks.end(); ++i) {
				std::cerr << *i;
				if (i + 1 != ranks.end()) {
					std::cerr << ',';
				}
			}
			std::cerr << "])\n";
		}

		expression_node operator[](size_t i) const {
			return expression_node(getnth(parsers[i], ranks[i]));
		}

		bool verify(size_t i, const expression_node &node) const {
			return true;
		}
	private:
		std::vector<interpreter<raw_expression_node> *> &parsers;
		const std::vector<size_t> &ranks;
	};

public:
	treeparser(const production *P_, const ordered_segments *segs_,
		       std::vector<interpreter<raw_expression_node> *> &parsers_, const score_combiner &rscores_,
			   const ordered_segments *headspan_, const ordered_segments *tailspan_)
		: P(P_), segs(segs_), rscores(rscores_), headspan(headspan_), tailspan(tailspan_) {
		assert(parsers_.size() == P->rhs.size());
		parsers.reserve(parsers_.size());
		for (std::vector<interpreter<raw_expression_node> *>::const_iterator i = parsers_.begin(); i != parsers_.end(); ++i) {
			parsers.push_back(new iterator<raw_expression_node>(*i, false));
		}
		std::vector<size_t> ranks(P->rhs.size(), 0);
		addtree(basic_tree_accessor(parsers), ranks);
	}
	~treeparser() {
		for (std::vector<interpreter<raw_expression_node> *>::iterator i = parsers.begin(); i != parsers.end(); ++i) {
			delete *i;
		}
	}

private:
	void addtree(const subtree_accessor &acc, const std::vector<size_t> &ranks) {
		std::cerr << "addtree(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			std::cerr << " " << P->rhs[i]->name;
		}
		std::cerr << ", [";
		for (std::vector<size_t>::const_iterator i = ranks.begin(); i != ranks.end(); ++i) {
			std::cerr << *i;
			if (i + 1 != ranks.end()) {
				std::cerr << ',';
			}
		}
		std::cerr << "])\n";
		expression_node node = P->tbuild->build(acc);
		if (node) {
			if (!node->nt) {
				node->nt = P->nt;
			}
			if (!node->nt_low) {
				node->nt_low = P->nt;
			}
			node->span = const_cast<ordered_segments *>(segs);
			node->set_string_builder(P->sbuild);
			node->set_latex_builder(P->lbuild);
			node->prod = P;
			node->set_score(node->score_combo().add(rscores));
			node->head_span = const_cast<ordered_segments *>(headspan);
			node->tail_span = const_cast<ordered_segments *>(tailspan);
			std::cerr << "adding known parse " << *node.ptr() << " with score " << node->score() << " making " << known.size() + 1 << std::endl;
			node->add_ref();
			known.push(treerec(node.ptr(), ranks));
			known_ranks.insert(ranks);
		}
	}

	raw_expression_node *getnext() {
		std::cerr << this << "->getnext(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			std::cerr << " " << P->rhs[i]->name;
		}
		std::cerr << ")\n";

		if (!last.ranks.empty()) {
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				++last.ranks[i];
				if (known_ranks.find(last.ranks) == known_ranks.end()) {
					addtree(ranked_tree_accessor(parsers, last.ranks), last.ranks);
				}
				--last.ranks[i];
			}
		}

		if (known.empty()) {
			std::cerr << " (none)\n";
			return 0;		
		}

		last = known.top();
		known.pop();
		std::cerr << " " << *last.node << std::endl;
		return last.node;
	}

private:
	struct treerec {
		raw_expression_node *node;
		std::vector<size_t> ranks;
		treerec() : node(0) { }
		explicit treerec(raw_expression_node *node_, const std::vector<size_t> &ranks_) : node(node_), ranks(ranks_) { }
		bool operator<(const treerec &rhs) const {
			return node->score() < rhs.node->score();
		}
	};

private:
	const production *P;
	const ordered_segments *segs;
	std::vector<interpreter<raw_expression_node> *> parsers;
	score_combiner rscores;
	const ordered_segments *headspan;
	const ordered_segments *tailspan;
	std::priority_queue<treerec> known;
	treerec last;
	std::set<std::vector<size_t> > known_ranks;
};


static size_t
minlen(const production *P, size_t start, size_t len) {
	assert(start + len <= P->rhs.size());
	if (len == 0) {
		return 0;
	}
	if (start == 0) {
		return P->min_prefix_length[len - 1];
	}
	else {
		return P->min_prefix_length[start + len - 1] - P->min_prefix_length[start - 1];
	}
}

static size_t
maxlen(const production *P, size_t start, size_t len) {
	assert(start + len <= P->rhs.size());
	if (len == 0) {
		return 0;
	}
	if (start == 0) {
		return P->max_prefix_length[len - 1];
	}
	else if (P->max_prefix_length[start + len - 1] == std::numeric_limits<size_t>::max()) {
		if (P->max_prefix_length[start - 1] == std::numeric_limits<size_t>::max()) {
			size_t tlen = 0;
			for (size_t i = start; i < start + len; ++i) {
				if (P->rhs[i]->max_length == std::numeric_limits<size_t>::max()) {
					return std::numeric_limits<size_t>::max();
				}
				tlen += P->rhs[i]->max_length;
			}
			return tlen;
		}
		else {
			return std::numeric_limits<size_t>::max();
		}
	}
	else {
		return P->max_prefix_length[start + len - 1] - P->max_prefix_length[start - 1];
	}
}

static size_t
sumormax(size_t a, size_t b) {
	if (a == std::numeric_limits<size_t>::max() || b == std::numeric_limits<size_t>::max()) {
		return std::numeric_limits<size_t>::max();
	}
	else {
		return a + b;
	}
}

struct prodlenspec {
	size_t minprelen;
	size_t maxprelen;
	size_t minxlen;
	size_t maxxlen;
	size_t minpostlen;
	size_t maxpostlen;

	prodlenspec(const production *P, size_t start, size_t xstart, size_t xlen, size_t end,
	            const ordered_segments *segs, size_t startseg) {
		size_t left = end - (xstart + xlen);
		size_t minprelen_fromprod = minlen(P, start, xstart - start);
		size_t minpostlen_fromprod = minlen(P, xstart + xlen, left);
		size_t maxprelen_fromprod = maxlen(P, start, xstart - start);
		size_t maxpostlen_fromprod = maxlen(P, xstart + xlen, left);
		size_t minxlen_fromprod = minlen(P, xstart, xlen);
		size_t maxxlen_fromprod = maxlen(P, xstart, xlen);

		size_t nsegs = segs->size() - startseg;
		size_t s = sumormax(maxxlen_fromprod, maxpostlen_fromprod);
		minprelen = std::max(minprelen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(minxlen_fromprod, minpostlen_fromprod);
		maxprelen = std::min(maxprelen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(maxprelen_fromprod, maxpostlen_fromprod);
		minxlen = std::max(minxlen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(minprelen_fromprod, minpostlen_fromprod);
		maxxlen = std::min(maxxlen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(maxprelen_fromprod, maxxlen_fromprod);
		minpostlen = std::max(minpostlen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(minprelen_fromprod, minxlen_fromprod);
		maxpostlen = std::min(maxpostlen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));

		std::cerr << "prodlenspec(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			std::cerr << " " << P->rhs[i]->name;
		}
		std::cerr << ", start " << start << ", xstart " << xstart << ", xlen " << xlen << ", end " << end << ", bits " << segs->bits << ", startseg " << startseg << ")\n";
		std::cerr << " pre=(" << minprelen << "," << maxprelen << ")  ";
		std::cerr << " x=(" << minxlen << "," << maxxlen << ")  ";
		std::cerr << " post=(" << minpostlen << "," << maxpostlen << ")\n";
	}

	inline bool valid() const {
		return minprelen <= maxprelen && minxlen <= maxxlen && minpostlen <= maxpostlen;
	}
};

static double
rscore(const production *P, size_t i,
       interpreter<raw_expression_node> *headintrp, const ordered_segments *headsegs,
	   interpreter<raw_expression_node> *tailintrp, const ordered_segments *tailsegs) {
	const ordered_segments *head;
	const ordered_segments *tail;
	switch (P->attach_modes[i].from) {
	case attach_mode::GROUP:
		head = headsegs;
		break;
	case attach_mode::SYMBOL:
		raw_expression_node *node = getfirst(headintrp);
		assert(node);
		head = node->head_span;
		break;
	}
	switch (P->attach_modes[i].to) {
	case attach_mode::GROUP:
		tail = tailsegs;
		break;
	case attach_mode::SYMBOL:
		raw_expression_node *node = getfirst(tailintrp);
		assert(node);
		tail = node->tail_span;
		break;
	}
	return P->rel->membership(head->bounds(), tail->bounds(), 0, 0);
}

class prodparser : public multiplexor<raw_expression_node> {
public:
	prodparser(const context *ctx_, const production *P_, const ordered_segments *segs_)
		: ctx(ctx_), P(P_), segs(segs_), childsegs(P->rhs.size()), parsers(P->rhs.size()), rscores(P->rhs.size() - 1) {
		if (segs->size() < P->min_prefix_length[P->rhs.size()-1]) {
			return;
		}
		std::cerr << "prodparser(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			std::cerr << " " << P->rhs[i]->name;
		}
		std::cerr << ", " << segs->bits << ")\n";
		conf(0, 0);
	}

private:
	void confnt(size_t cterm, size_t cseg, size_t crhsi) {
		std::cerr << "confnt(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			std::cerr << " " << P->rhs[i]->name;
		}
		std::cerr << ", bits " << segs->bits << ", term " << cterm << ", seg " << cseg << ", rhs " << crhsi << std::endl;

		if (crhsi == P->rhs.size()) {
			score_combiner rscom;
			for (std::vector<double>::const_iterator i = rscores.begin(); i != rscores.end(); ++i) {
				rscom = rscom.addrel(*i);
			}
			addparser(new treeparser(P, segs, parsers, rscom, headspan, tailspan));
		}
		else if (cterm < P->terminal_indices.size() && crhsi == P->terminal_indices[cterm]) {
			if (P->rel && crhsi > 0) {
				double rsc = rscore(P, crhsi - 1, parsers[crhsi - 1], childsegs[crhsi - 1], parsers[crhsi], childsegs[crhsi]);
				if (rsc <= 0.0) {
					return;
				}
				rscores[crhsi - 1] = rsc;
			}
			if (crhsi == P->head) {
				raw_expression_node *node = getfirst(parsers[crhsi]);
				headspan = node->head_span;
			}
			if (crhsi == P->tail) {
				raw_expression_node *node = getfirst(parsers[crhsi]);
				tailspan = node->tail_span;
			}
			confnt(cterm + 1, cseg + childsegs[P->terminal_indices[cterm]]->bits.size(), crhsi + 1);
		}
		else {
			prodlenspec lenspec(P, crhsi, crhsi, 1, P->rhs.size(), segs, cseg);
			if (!lenspec.valid()) {
				return;
			}
			size_t d = P->rel ? ordering_for_relation(P->rel) : 0;
			for (size_t len = lenspec.minxlen; len <= lenspec.maxxlen; ++len) {
				std::cerr << " confnt trying length " << len << std::endl;
				const ordered_segments *subsegs = getsegs(segs->slice_bits(d, cseg, cseg + len), segs, d, cseg, cseg + len);
				parsers[crhsi] = mkparser(ctx, P->rhs[crhsi], subsegs);
				if (parsers[crhsi]) {
					if (P->rel && crhsi > 0) {
						double rsc = rscore(P, crhsi - 1, parsers[crhsi - 1], childsegs[crhsi - 1], parsers[crhsi], subsegs);
						if (rsc <= 0.0) {
							return;
						}
						rscores[crhsi - 1] = rsc;
					}
					childsegs[crhsi] = subsegs;
					if (crhsi == P->head) {
						raw_expression_node *node = getfirst(parsers[crhsi]);
						headspan = node->head_span;
					}
					if (crhsi == P->tail) {
						raw_expression_node *node = getfirst(parsers[crhsi]);
						tailspan = node->tail_span;
					}
					confnt(cterm, cseg + len, crhsi + 1);
				}
			}
		}
	}

	void conf(size_t cterm, size_t cseg) {
		std::cerr << "conf(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			std::cerr << " " << P->rhs[i]->name;
		}
		std::cerr << ", bits " << segs->bits << ", term " << cterm << ", seg " << cseg << std::endl;
		if (cterm == P->terminal_indices.size()) {
			confnt(0, 0, 0);
		}
		else {
			size_t last_tpos;
			if (cterm == 0) {
				last_tpos = 0;
			}
			else {
				last_tpos = P->terminal_indices[cterm-1]+1;
			}
			size_t tpos = P->terminal_indices[cterm];
			prodlenspec lenspec(P, last_tpos, tpos, 1, P->rhs.size(), segs, cseg);
			if (!lenspec.valid()) {
				return;
			}

			size_t d = P->rel ? ordering_for_relation(P->rel) : 0;
			const nonterminal *termnt = P->rhs[tpos];

			for (size_t start = cseg + lenspec.minprelen; start <= cseg + lenspec.maxprelen; ++start) {
				for (size_t len = lenspec.minxlen; len <= lenspec.maxxlen; ++len) {
					std::cerr << " conf trying start " << start << ", len " << len << std::endl;
					const ordered_segments *subsegs = getsegs(segs->slice_bits(d, start, start + len), segs, d, start, start + len);
					std::map<bitvec, group *>::const_iterator pgrp = ctx->groups->find(subsegs->bits);
					if (pgrp != ctx->groups->end()) {
						const group *grp = pgrp->second;
						std::string tname = termnt->name.substr(2);
						for (std::vector<match_score>::const_iterator pmtch = grp->final_matches.begin(); pmtch != grp->final_matches.end(); ++pmtch) {
							if (pmtch->proto->info().name == tname) {
								std::cerr << "found terminal " << tname << " at " << subsegs->bits << std::endl;
								if (P->rel && tpos != 0) {
									const ordered_segments *presegs = getsegs(segs->slice_bits(d, cseg, start - cseg), segs, d, cseg, start - cseg);
									double rsc = P->rel->membership(presegs->bounds(), subsegs->bounds(), 0, 0);
									if (rsc <= 0.0) {
										return;
									}
								}
								staticintrpr<raw_expression_node> *tintrpr = new staticintrpr<raw_expression_node>;
								expression_node tnode = make_terminal(tname, pmtch->score);
								tnode->head_span = tnode->tail_span = const_cast<ordered_segments *>(subsegs);
								tnode->add_ref();
								tintrpr->addknown(tnode.ptr(), true);
								childsegs[tpos] = subsegs;
								parsers[tpos] = tintrpr;
								conf(cterm + 1, start + len);
								break;
							}
						}
					}
				}
			}
		}
	}

private:
	const context *ctx;
	const production *P;
	const ordered_segments *segs;

	std::vector<const ordered_segments *> childsegs;
	std::vector<interpreter<raw_expression_node> *> parsers;

	std::vector<double> rscores;

	const ordered_segments *headspan;
	const ordered_segments *tailspan;
};

class ntparser : public multiplexor<raw_expression_node> {
public:
	ntparser(const context *ctx, const nonterminal *nt, const ordered_segments *segs) {
		std::cerr << "ntparser(" << nt->name << ":" << segs->bits << ")\n";
		for (std::vector<production *>::const_iterator i = nt->productions.begin(); i != nt->productions.end(); ++i) {
			addparser(new prodparser(ctx, *i, segs));
		}
	}
};

static interpreter<raw_expression_node> *
mkparser(const context *ctx, const nonterminal *nt, const ordered_segments *segs) {
	std::cerr << "mkparser(" << nt->name << ":" << segs->bits << ")\n";
	if (segs->size() < nt->min_length || segs->size() > nt->max_length) {
		return 0;
	}
	std::map<const nonterminal *, interpreter<raw_expression_node> *> &tab = parsers[segs];
	std::map<const nonterminal *, interpreter<raw_expression_node> *>::iterator i = tab.find(nt);
	interpreter<raw_expression_node> *intrpr;
	if (i == tab.end()) {
		intrpr = new ntparser(ctx, nt, segs);
		if (intrpr->nknown() == 0) {
			raw_expression_node *test = intrpr->next();
			if (!test) {
				delete intrpr;
				intrpr = 0;
			}
		}
		tab.insert(i, std::make_pair(nt, intrpr));
	}
	else {
		intrpr = i->second;
		std::cerr << " (exists at " << intrpr << ")\n";
	}
	return intrpr;
}

}
