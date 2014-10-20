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

#include <deque>
#include <cassert>
#include <cstdlib>



namespace scg
{


static double SCORE_THRESHOLD = RegisterParameterDouble("ScoreThreshold", &SCORE_THRESHOLD);


struct parsepos {
	const production *P;
	unsigned pos;
	unsigned segpos;
	parse_links links;

	parsepos(const production *P_, unsigned pos_, unsigned segpos_) : P(P_), pos(pos_), segpos(segpos_), links(P_) { }
	parsepos(const production *P_, unsigned pos_, unsigned segpos_, const parse_links &links_) : P(P_), pos(pos_), segpos(segpos_), links(links_) { }
};


struct termparsepos {
	const production *P;
	unsigned termpos;
	unsigned segpos;
	std::vector<parse_links> links;

	termparsepos(const production *P_, unsigned termpos_, unsigned segpos_) : P(P_), termpos(termpos_), segpos(segpos_) { }
	termparsepos(const production *P_, unsigned termpos_, unsigned segpos_, const std::vector<parse_links> &links_) : P(P_), termpos(termpos_), segpos(segpos_), links(links_) { }
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



static bool parse_internal(context *ctx, ordered_segments *segs, const nonterminal *nt);

double
compute_membership(const relation *rel, const ordered_segments *from, const ordered_segments *to, const std::set<rclass_t> &rc1, const std::set<rclass_t> &rc2)
{
	//std::vector<double> M;
	//measure(b1, b2, M);

	std::set<rclass_t> cls1 = rc1;
	std::set<rclass_t> cls2 = rc2;
	if (cls1.empty()) {
		cls1.insert(BOX_CLASS);
	}
	cls1.insert(AGGREGATE_CLASS);
	if (cls2.empty()) {
		cls2.insert(BOX_CLASS);
	}
	cls2.insert(AGGREGATE_CLASS);


	VERBOSE(
		*verb_out << "compute membership in " << rel->name << " from " << from->bounds() << " to " << to->bounds() << " with classes { ";
		for (std::set<rclass_t>::const_iterator i = cls1.begin(); i != cls1.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "} and { ";
		for (std::set<rclass_t>::const_iterator i = cls2.begin(); i != cls2.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "}\n";
	);

	double mem = -1.0;

	for (std::set<rclass_t>::reverse_iterator i1 = cls1.rbegin(); i1 != cls1.rend(); ++i1) {
		for (std::set<rclass_t>::reverse_iterator i2 = cls2.rbegin(); i2 != cls2.rend(); ++i2) {
			double m = rel->membership(from, to, *i1, *i2);
			VERBOSE(*verb_out << "classes " << *i1 << " and " << *i2 << " give " << m << std::endl);
			if (m != -1.0) {
				mem = std::max(mem, m);
			}
			if (*i2 < NUM_RCLASSES && mem != -1.0) {
				break;
			}
		}
		if (*i1 < NUM_RCLASSES && mem != -1.0) {
			break;
		}
	}

	VERBOSE(*verb_out << " : " << mem << std::endl);
	return mem;

	/*
	
	if (cls1.empty()) {
		if (cls2.empty()) {
			//return rel->membership(M, BOX_CLASS, BOX_CLASS);
			mem = rel->membership(b1, b2, BOX_CLASS, BOX_CLASS);
		}
		else {
			for (std::set<rclass_t>::const_iterator c2 = cls2.begin(); c2 != cls2.end(); ++c2) {
				//mem = std::max(mem, rel->membership(M, BOX_CLASS, *c2));
				mem = std::max(mem, rel->membership(b1, b2, BOX_CLASS, *c2));
			}
		}
	}

	else if (cls2.empty()) {
		for (std::set<rclass_t>::const_iterator c1 = cls1.begin(); c1 != cls1.end(); ++c1) {
			//mem = std::max(mem, rel->membership(M, *c1, BOX_CLASS));
			mem = std::max(mem, rel->membership(b1, b2, *c1, BOX_CLASS));
		}
	}
	
	else {
		for (std::set<rclass_t>::const_iterator c1 = cls1.begin(); c1 != cls1.end(); ++c1) {
			for (std::set<rclass_t>::const_iterator c2 = cls2.begin(); c2 != cls2.end(); ++c2) {
				//mem = std::max(mem, rel->membership(M, *c1, *c2));
				mem = std::max(mem, rel->membership(b1, b2, *c1, *c2));
			}
		}
	}

	VERBOSE(*verb_out << " : " << mem << std::endl);
	
	return mem;*/
}

static std::set<rclass_t> EMPTY_CLASS_SET;


ordered_segments *
lookup_span(context *ctx, const bitvec &bits, ordered_segments *src, unsigned dim, unsigned start, unsigned end)
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
/*		
	if (i != ctx->spans->end() && i->second->edit_version == ctx->edits.size()) {
		VERBOSE(*verb_out << "found up-to-date pointer " << i->second << std::endl);
		return i->second;
	}
	else {
		bitvec scratch = bits;
		unsigned rollbacks = 0;
		for (std::list<edit_op>::const_iterator j = ctx->edits.begin(); j != ctx->edits.end(); ++j) {
			++rollbacks;

			if (j->type == edit_op::INSERTION) {
				if (scratch.at(j->pos)) {
					break;
				}
				scratch.erase(j->pos);
				VERBOSE(*verb_out << "rolling back insertion at " << j->pos << "; bits became " << scratch << std::endl);
			}
			else {
				scratch.insert(j->pos, false);
				VERBOSE(*verb_out << "rolling back removal at " << j->pos << "; bits became " << scratch << std::endl);
			}

			i = ctx->spans->find(scratch);
			if (i != ctx->spans->end() && i->second->edit_version == ctx->edits.size() - rollbacks) {
				VERBOSE(*verb_out << "found pointer at " << i->second << "; bringing it up-to-date\n");
				ordered_segments *span = i->second;
				span->bits = bits;
				span->edit_version = ctx->edits.size();
				ctx->spans->erase(i);
				(*ctx->spans)[bits] = span;
				return span;
			}
			break;
		}

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
	}*/
}


bool
parse_nt(context *ctx, const parsepos &startpos, ordered_segments *segs, unsigned stop_pos, std::vector<parse_links> &links)
{
	std::deque<parsepos> agenda;

	VERBOSE(*verb_out << "parsing nonterminals " << startpos << " to " << stop_pos << " at " << segs->bits << std::endl);

	agenda.push_back(startpos);

	unsigned N = segs->nsegments();

	while (!agenda.empty()) {
		parsepos &pos = agenda.front();
		const production *P = pos.P;
		unsigned d = P->dim;

		VERBOSE(*verb_out << " popping " << pos << std::endl);

		if (pos.pos == stop_pos) {
			links.push_back(pos.links);
		}
		else {
			unsigned min_before;
			if (pos.pos == 0) {
				min_before = 0;
			}
			else {
				if (startpos.pos == 0) {
					min_before = P->min_prefix_length[pos.pos - 1];
				}
				else {
					min_before = P->min_prefix_length[pos.pos - 1] - P->min_prefix_length[startpos.pos - 1];
				}
			}
			unsigned min_after = P->min_prefix_length[stop_pos - 1] - P->min_prefix_length[pos.pos];

			VERBOSE(*verb_out << "  min_before=" << min_before << "; min_after=" << min_after << std::endl);

			const nonterminal *currnt = P->rhs[pos.pos];
			unsigned min_length = currnt->min_length;
			if (min_after == 0) {
				min_length = std::max(min_length, segs->size() - pos.segpos);
			}

			for (unsigned len = min_length; len <= segs->size() - min_before - min_after; ++len) {
				VERBOSE(*verb_out << "   using len=" << len << std::endl);
				if (pos.segpos + len > segs->size()) {
					VERBOSE(*verb_out << "   cannot fit in segments with length " << len << std::endl);
					break;
				}

				bitvec currbits = segs->slice_bits(d, pos.segpos, pos.segpos + len);
				assert(currbits.count_set_bits() == len);

				ordered_segments *currsegs = lookup_span(ctx, currbits, segs, d, pos.segpos, pos.segpos + len);
				/*ordered_segments *&currsegs = (*ctx->spans)[currbits];
				if (!currsegs) {
					VERBOSE(*verb_out << "making new slice for " << currbits << std::endl << "(there are " << (*ctx->spans).size() << " spans now)\n");
					currsegs = segs->slice(d, pos.segpos, pos.segpos + len);
				}*/

				parseref &ref = pos.links.children[pos.pos-1];

				if (pos.pos == startpos.pos
				 || compute_membership(P->rel, ref.span, currsegs, ref.nt->rclasses, currnt->rclasses) > SCORE_THRESHOLD) {
				
					if (parse_internal(ctx, currsegs, currnt) > 0) {
						pos.links.children[pos.pos] = parseref(currsegs, currnt);
						agenda.push_back(parsepos(P, pos.pos + 1, pos.segpos + len, pos.links));
						VERBOSE(*verb_out << " pushing " << agenda.back() << std::endl);
					}
				}
				else {
					VERBOSE(*verb_out << "     local rel score is too low to continue\n");
				}
			}
		}

		agenda.pop_front();
	}

	VERBOSE(
		*verb_out << "returning from " << startpos << " at " << segs->bits << "; there are " << links.size() << " results, as follows:\n";
		for (std::vector<parse_links>::const_iterator i = links.begin(); i != links.end(); ++i) {
			*verb_out << ' ' << *i << std::endl;
		}
	);

	return !links.empty();
}


static void
cross_links(std::vector<parse_links> &dst, const std::vector<parse_links> &rhs, unsigned pos)
{
	if (rhs.empty()) {
		return;
	}

	std::vector<parse_links>::iterator i = dst.begin(); 
	unsigned n = dst.size();
	while (n--) {
		std::vector<parse_links>::const_iterator j = rhs.begin();
		std::copy(j->children.begin() + pos, j->children.end(), i->children.begin() + pos);
		if (rhs.size() > 1) {
			++j;
			for (; j != rhs.end(); ++j) {
				unsigned index = i - dst.begin();
				dst.push_back(*i);
				i = dst.begin() + index;
				parse_links &newlinks = dst.back();
				std::copy(j->children.begin() + pos, j->children.end(), newlinks.children.begin() + pos);
			}
		}
		++i;
	}
}


static bool
parse_internal(context *ctx, ordered_segments *segs, const nonterminal *nt)
{
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

	cell = new parse_cell(segs, nt);

	VERBOSE(*verb_out << "parsing " << nt->name << " at " << segs->bits << "(" << nt << " at " << segs << ")" << std::endl);
	
	/*
	if (nt->name[0] == '!')  { // indicates matrix
		VERBOSE(*verb_out << "creating matrix for strokes " << segs->bits << std::endl);
		cell->matrix = CreateMatrixAnalyzer(ctx->wrapper);
		for (std::vector<const segment *>::const_iterator i = segs->segs[0].begin(); i != segs->segs[0].end(); ++i) {
			cell->matrix->addStroke(&(*i)->stk);
		}
		cell->mxiter = cell->matrix->createIterator();
		bool results = (cell->mxiter->next() != 0);
		delete cell->mxiter;
		if (!results) {
			DestroyMatrixAnalyzer(cell->matrix);
		}
		return results;
	}*/

	std::vector<parse_links> &links = cell->links;

	std::deque<termparsepos> agenda;
	for (std::list<production>::const_iterator i = nt->productions.begin(); i != nt->productions.end(); ++i) {
		if (i->min_prefix_length[i->rhs.size()-1] <= segs->size()) {
			if (i->rhs.size() > 1 || i->rhs[0]->max_length == 0 || segs->size() <= i->rhs[0]->max_length) {
				agenda.push_back(termparsepos(&*i, 0, 0));
				VERBOSE(*verb_out << " pushing " << agenda.back() << std::endl);
			}
		}
	}

	while (!agenda.empty()) {
		termparsepos &pos = agenda.front();
		const production *P = pos.P;
		unsigned d = P->dim;

		VERBOSE(*verb_out << " popping " << pos << std::endl);

		if (pos.termpos == P->terminal_indices.size()) {
			if (pos.segpos == 0) {
				if (parse_nt(ctx, parsepos(P, 0, 0), segs, P->rhs.size(), pos.links)) {
					links.insert(links.end(), pos.links.begin(), pos.links.end());
					VERBOSE(*verb_out << "found parses of " << nt->name << " on " << segs->bits << std::endl);
				}
			}
			else {
				if (pos.segpos < segs->size()) {
					bitvec currbits = segs->slice_bits(d, pos.segpos, segs->size());
					ordered_segments *currsegs = lookup_span(ctx, currbits, segs, d, pos.segpos, segs->size());
					/*ordered_segments *&currsegs = (*ctx->spans)[currbits];
					if (!currsegs) {
						currsegs = segs->slice(d, pos.segpos, segs->size());
					}*/
					std::vector<parse_links> sublinks;
					parsepos ntpos(P, P->terminal_indices[pos.termpos-1]+1, 0);
					if (parse_nt(ctx, ntpos, currsegs, P->rhs.size(), sublinks)) {
						cross_links(pos.links, sublinks, ntpos.pos);
						links.insert(links.end(), pos.links.begin(), pos.links.end());
						VERBOSE(*verb_out << "found parses of " << nt->name << " on " << segs->bits << std::endl);
					}
				}
				else {
					links.insert(links.end(), pos.links.begin(), pos.links.end());
					VERBOSE(*verb_out << "found parses of " << nt->name << " on " << segs->bits << std::endl);
				}
			}
		}
		else {
			unsigned term_index = P->terminal_indices[pos.termpos];
			const nonterminal *termnt = P->rhs[term_index];

			unsigned min_postlength = P->min_prefix_length[P->rhs.size()-1] - P->min_prefix_length[term_index];

			unsigned min_prelength;
			if (pos.termpos == 0) {
				if (term_index == 0) {
					min_prelength = 0;
				}
				else {
					min_prelength = P->min_prefix_length[term_index-1];
				}
			}
			else {
				min_prelength = P->min_prefix_length[term_index-1] - P->min_prefix_length[P->terminal_indices[pos.termpos - 1]];
			}

			unsigned max_prelength;
			if (term_index == 0) {
				max_prelength = 0;
			}
			else {
				max_prelength = segs->size() - (P->min_prefix_length[P->rhs.size()-1] - P->min_prefix_length[term_index-1]) - pos.segpos;
			}

			VERBOSE(*verb_out << " term_index=" << term_index << "; min_prelength=" << min_prelength << "; max_prelength=" << max_prelength << std::endl);

			unsigned min_length = termnt->min_length;

			if (min_postlength == 0) {
				if (min_prelength == 0) {
					min_length = segs->size() - pos.segpos;
				}
				else {
					min_prelength = std::max(min_prelength, segs->size() - std::min(segs->size(), termnt->max_length + pos.segpos));
				}
			}
			else if (min_prelength == 0) {
				min_postlength = std::max(min_postlength, segs->size() - std::min(segs->size(), termnt->max_length + pos.segpos));
			}

			for (unsigned prelen = min_prelength; prelen <= max_prelength; ++prelen) {
				unsigned seg_term_index = pos.segpos + prelen;
				for (unsigned termlen = min_length; termlen <= termnt->max_length; ++termlen) {
					VERBOSE(*verb_out << "  trying prelength " << prelen << "; terminal length " << termlen << std::endl);
					if (pos.segpos + prelen + termlen + min_postlength > segs->size()) {
						VERBOSE(*verb_out << "   cannot fit in segments with total min length " << prelen + termlen + min_postlength << std::endl);
						break;
					}
					bitvec currbits = segs->slice_bits(d, seg_term_index, seg_term_index + termlen);
					VERBOSE(*verb_out << "    currbits=" << currbits << std::endl);
					ordered_segments *currsegs = lookup_span(ctx, currbits);
					//ordered_segments *currsegs = (*ctx->spans)[currbits];
					if (currsegs) {
						const parse_cell *term_cell = ctx->table[currsegs][termnt];
						VERBOSE(
							*verb_out << "cell is " << term_cell << std::endl;
							*verb_out << "current options in span " << currsegs << " are: ";
							const nonterminal_parses &opts = ctx->table[currsegs];
							for (nonterminal_parses::const_iterator i = opts.begin(); i != opts.end(); ++i) {
								*verb_out << "  " << i->first->name << " -> " << i->second << std::endl;
							}
						);

						if (term_cell) {
							VERBOSE(*verb_out << "     terminal " << termnt->name << " exists here\n");
							if (prelen > 0) {
								bitvec prebits = segs->slice_bits(d, pos.segpos, seg_term_index);
								ordered_segments *presegs = lookup_span(ctx, prebits, segs, d, pos.segpos, seg_term_index);
								/*ordered_segments *&presegs = (*ctx->spans)[prebits];
								if (!presegs) {
									presegs = segs->slice(d, pos.segpos, seg_term_index);
								}
								*/

								double local_score = compute_membership(P->rel, presegs, currsegs, EMPTY_CLASS_SET, termnt->rclasses);

								if (local_score > SCORE_THRESHOLD) {
									std::vector<parse_links> prefix_links;
									parsepos ntpos(P, pos.termpos == 0 ? 0 : P->terminal_indices[pos.termpos-1]+1, 0);
									if (parse_nt(ctx, ntpos, presegs, term_index, prefix_links) > 0) {
										if (pos.links.empty()) {
											assert(pos.termpos == 0);
											pos.links = prefix_links;
										}
										else {
											cross_links(pos.links, prefix_links, ntpos.pos);
										}
										for (std::vector<parse_links>::iterator i = pos.links.begin(); i != pos.links.end(); ++i) {
											i->children[term_index] = parseref(currsegs, termnt);
										}
										agenda.push_back(termparsepos(P, pos.termpos + 1, pos.segpos + prelen + termlen, pos.links));
										VERBOSE(*verb_out << " pushing " << agenda.back() << std::endl);
									}
								}
							}
							else {
								assert(P != 0);
								if (pos.links.empty()) {
									pos.links.push_back(parse_links(P));
									pos.links.back().children[0] = parseref(currsegs, termnt);
								}
								else {
									std::vector<parse_links> termlinks;
									termlinks.push_back(parse_links(P));
									termlinks.back().children[term_index] = parseref(currsegs, termnt);
									cross_links(pos.links, termlinks, term_index);
								}
								agenda.push_back(termparsepos(P, pos.termpos + 1, pos.segpos + termlen, pos.links));
								VERBOSE(*verb_out << " pushing " << agenda.back() << std::endl);
							}
						}
					}
				}
			}
		}

		agenda.pop_front();
	}

	for (std::vector<parse_links>::iterator i = links.begin(); i != links.end(); ++i) {
		if (i->children.size() > 1) {
			cell->classes.insert(BOX_CLASS);
		}
		else {
			const parseref &ref = i->children.front();
			parse_cell *child = ctx->table[ref.span][ref.nt];
			cell->classes.insert(child->classes.begin(), child->classes.end());
		}
	}

	VERBOSE(
		*verb_out << "returning from " << nt->name << " at " << segs->bits << "; there are " << links.size() << " results, as follows:\n";
		for (std::vector<parse_links>::const_iterator i = links.begin(); i != links.end(); ++i) {
			*verb_out << ' ' << *i << std::endl;
		}
	);

	return !links.empty();
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
	rebuild_production_lengths(ctx->G, ctx->G.start);
	return parse_internal(ctx, segs, nt);
}


class parser_subtree_accessor : public subtree_accessor {
public:
	parser_subtree_accessor(const parse_iterator &it_, const parse_links &links_, const std::vector<unsigned> &ranks_, std::vector<expression_node> &nodes_)
		: it(it_), links(links_), ranks(ranks_), nodes(nodes_) { }

	expression_node operator[](unsigned i) const {
		return extract_tree_internal(it, links.children[i].span, links.cells[i], ranks[i]);
	}

	bool verify(unsigned i, const expression_node &node) const {
		return links.children[i].nt == node->nt_low || links.children[i].nt->has_direct_descendent(node->nt_low);
	}

private:
	const parse_iterator &it;
	const parse_links &links;
	const std::vector<unsigned> &ranks;
	std::vector<expression_node> &nodes;
};


static expression_node
build_tree(const parse_iterator &it, const parse_links &links, const std::vector<unsigned> &ranks, std::vector<expression_node> &nodes)
{
	nodes.resize(links.P->rhs.size());
	parser_subtree_accessor acc(it, links, ranks, nodes);
	expression_node tree = links.P->tbuild->build(acc);
	if (tree) {
		/*
		for (unsigned i = 0; i < tree->nchildren(); ++i) {
			assert(tree->mchild(i)->parent);
		}*/

		if (tree->nchildren() == 1) {
			tree->rclasses = tree->mchild(0)->rclasses;
			tree->set_matrix(tree->mchild(0)->matrix());
		}
		if (tree->rclasses.empty()) {
			tree->rclasses.insert(BOX_CLASS);
		}
	}

	return tree;
}


static expression_node
build_tree(const parse_iterator &it, const parse_links &links, const std::vector<unsigned> &ranks)
{
	std::vector<expression_node> dummy;
	return build_tree(it, links, ranks, dummy);
}


/*
expression_node
make_matrix_expression(parse_cell *mxcell, MatrixAnalyzer *mxan, Matrix *mx)
{
	expression_node tree(new raw_expression_node);
	tree->set_type(MATRIX_EXPR);

	tree->set_matrix(mxan);

	unsigned nrows, ncols;
	mx->getDimensions(nrows, ncols);

	std::stringstream ss;
	ss << nrows;
	tree->add_child(make_terminal(ss.str(), 1.0, false));
	ss.str("");
	ss << ncols;
	tree->add_child(make_terminal(ss.str(), 1.0, false));
	ss.str("");

	score_combiner sc;
	expression_node rows(new raw_expression_node);

	for (unsigned i = 0; i < nrows; ++i) {
		std::stringstream row_ss;
		row_ss << "<mtr>";
		score_combiner row_sc;
		expression_node row(new raw_expression_node);
		row->set_type(MATRIXROW_EXPR);
		row->set_score(score_combiner(1.0));
		//bitvec bits(ctx->segments.size(), false);
		for (unsigned j = 0; j < ncols; ++j) {
			ExpressionIterator *pit;
			mx->getCell(i, j, pit);
			if (!pit) {
				return NULL_EXPRESSION_NODE;
			}
			iterator_base *it = (iterator_base *)pit;
			parse_cell *parent_cell = it->get_cell();
			if (parent_cell) {
				parent_cell->parents.push_back(mxcell);
			}
			expression_node cell = it->next_node(false);
			if (!cell) {
				return NULL_EXPRESSION_NODE;
			}
			row->add_child(cell);
			row_sc = row_sc.add(cell->score_combo());
			row_ss << "<mtd>" << cell->long_str() << "</mtd>";

			//bits.union_insert(cell->span->bits);
		}
		row_ss << "</mtr>";

		row->set_cell(0);
		//row->set_span(lookup_span(ctx, bits));
		row->set_nt(0);
		row->nt_low = 0;
		row->rclasses.insert(BOX_CLASS);
		row->set_score(row_sc);

		row->set_long_string(row_ss.str());
		rows->add_child(row);

		ss << row_ss.str();
	}

	sc = sc.add(rows);

	rows->set_long_string(ss.str());
	rows->rclasses.insert(BOX_CLASS);
	rows->set_score(sc);

	tree->add_child(rows);
	
	//tree->set_long_string(ss.str());
	tree->rclasses.insert(BOX_CLASS);
	tree->set_score(sc);	

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
	
	expression_node tree = make_matrix_expression(cell, cell->matrix, mx);
	if (!tree) {
		return NULL_EXPRESSION_NODE;
	}
	tree->set_nt(cell->nt);
	tree->nt_low = cell->nt;

	expression_node rows = tree->mchild(0);
	rows->set_nt(cell->nt);
	rows->nt_low = cell->nt;
	rows->set_context(ctx);
	rows->set_cell(cell);
	rows->set_span(span);

	for (unsigned i = 0; i < rows->nchildren(); ++i) {
		expression_node row = rows->mchild(i);
		row->set_context(ctx);
	}

	double conf;
	cell->matrix->getConfidence(mx, conf);
	score_combiner sc(conf);
	sc.add(rows->score_combo());
	tree->set_score(sc);

	return tree;
}*/


static expression_node extract_tree(const parse_iterator &it, parse_links &links, unsigned rank);

/*
double
nt_relscore::lookup(const std::vector<raw_expression_node *> &children, const parse_links &links, std::vector<unsigned>::const_iterator curr, std::vector<unsigned>::const_iterator last) const
{
	VERBOSE(
		*verb_out << "nt_relscore lookup in branch with rclasses ";
		for (std::map<rclass_t, typed_relscore *>::const_iterator i = scores.begin(); i != scores.end(); ++i) {
			*verb_out << i->first << ' ';
		}
		*verb_out << std::endl;
	);
	double best = 0.0;
	assert(*curr < children.size());
	const raw_expression_node *curr_child = children[*curr];
	assert(curr_child);
	VERBOSE(
		*verb_out << "child has rclasses ";
		for (std::set<rclass_t>::const_iterator i = curr_child->rclasses.begin(); i != curr_child->rclasses.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << std::endl;
	);

	for (std::set<rclass_t>::const_iterator i = curr_child->rclasses.begin(); i != curr_child->rclasses.end(); ++i) {
		std::map<rclass_t, typed_relscore *>::const_iterator j = scores.find(*i);
		if (j != scores.end()) { // scores can be pruned if they are too near zero
			assert(j->second);
			best = std::max(best, j->second->lookup(children, links, curr + 1, last));
		}
	}
	return best;
}

double
t_relscore::lookup(const std::vector<raw_expression_node *> &children, const parse_links &links, std::vector<unsigned>::const_iterator curr, std::vector<unsigned>::const_iterator last) const
{
	//assert(curr == last);
	VERBOSE(*verb_out << "t_relscore lookup at leaf with score " << score.score() << std::endl);
	return score.score();
}
*/

static double
measure_relation(const parse_links &links, unsigned from, unsigned to, const expression_node &from_tree, const expression_node &to_tree, rclass_t from_class, rclass_t to_class)
{
	assert(links.children.size() == links.P->rhs.size());
	
	const relation *rel = links.P->rel;
	const parseref &fromref = links.children[from];
	const parseref &toref = links.children[to];
	
	const ordered_segments *fromsegs = fromref.span;
	const ordered_segments *tosegs = toref.span;
	/*
	const attach_mode &mode = links.P->attach_modes[from];
	
	VERBOSE(*verb_out << "in measure_relation, mode is ");
	if (mode.from == attach_mode::GROUP) {
		fromsegs = fromref.span;
		VERBOSE(*verb_out << "GROUP " << from_box << " -> ");
	}
	else {
		assert(mode.from == attach_mode::SYMBOL);
		from_box = from_tree->tail_span;
		VERBOSE(*verb_out << "SYMBOL " << from_box << " -> ");
	}
	
	const ordered_segments *tosegs;
	if (mode.to == attach_mode::GROUP) {
		to_box = toref.span;
		VERBOSE(*verb_out << "GROUP " << to_box << " : ");
	}
	else {
		assert(mode.to == attach_mode::SYMBOL);
		to_box = to_tree->head_span;
		VERBOSE(*verb_out << "SYMBOL " << to_box << " : ");
	}*/
	
	double m = rel->membership(fromsegs, tosegs, from_class, to_class);
	VERBOSE(*verb_out << m << std::endl);
	return m;
}

/*
static score_combiner
combine_score_range(const parse_links &links, unsigned first, unsigned last, const score_combiner &curr_score)
{
	if (first == last || first + 1 == last) {
		return curr_score;
	}

	const parse_cell *curr = links.cells[first];
	assert(curr->classes.size() == 1);
	const parse_cell *next = links.cells[first + 1];
	assert(next->classes.size() == 1);

	double rel_score = measure_relation(links, first, first + 1, *curr->classes.begin(), *next->classes.begin());

	if (rel_score == SCORE_DEFER) {
		return combine_score_range(links, first + 1, last, curr_score);
	}
	
	else if (rel_score <= SCORE_THRESHOLD) {
		return score_combiner(0.0);
	}

	return combine_score_range(links, first + 1, last, curr_score.add(rel_score));
}


template <typename T>
static void
build_relscore_tree(typed_relscore *&score, const score_combiner &curr_score, parse_links &links, T first, T curr, T last, rclass_t last_rclass)
{
	VERBOSE(*verb_out << "build_relscore_tree() at " << links << " with current score " << curr_score.score() << std::endl);

	score_combiner active_score = curr_score;

	std::vector<parseref>::iterator base = links.children.begin();
	if (curr == last) {
		if (first == last) {
			active_score = combine_score_range(links, 0, links.children.size(), active_score);
		}
		else if (*(curr - 1) != links.children.size() - 1) {
			active_score = combine_score_range(links, *(curr - 1) + 1, links.children.size(), active_score);
			const parse_cell *leading_cell = links.cells[*(curr - 1) + 1];
			assert(leading_cell->classes.size() == 1);
			active_score = active_score.add(measure_relation(links, *(curr - 1), *(curr - 1) + 1, last_rclass, *leading_cell->classes.begin()));
		}

		if (active_score.score() > SCORE_THRESHOLD) {
			score = new t_relscore(active_score);
		}
		else {
			score = 0;
		}
		return;
	}


	if (curr != first && *(curr - 1) != *curr - 1) {
		active_score = combine_score_range(links, *(curr - 1) + 1, *curr - 1, active_score);
		const parse_cell *leading_cell = links.cells[*(curr - 1) + 1];
		assert(leading_cell->classes.size() == 1);
		active_score = active_score.add(measure_relation(links, *(curr - 1), *(curr - 1) + 1, last_rclass, *leading_cell->classes.begin()));

		if (active_score.score() <= SCORE_THRESHOLD) {
			score = 0;
			return;
		}
	}


	nt_relscore *nt_score = new nt_relscore;
	const parse_cell *cell = links.cells[*curr];
	if (curr != first && (*(curr - 1) == *curr - 1)) {
		for (std::set<rclass_t>::const_iterator i = cell->classes.begin(); i != cell->classes.end(); ++i) {
			score_combiner tmp_score(measure_relation(links, *curr - 1, *curr, last_rclass, *i));
			typed_relscore *&subscore = nt_score->scores[*i];
			assert(!subscore);
			build_relscore_tree(subscore, active_score.add(tmp_score), links, first, curr + 1, last, *i);
			if (!subscore) {
				nt_score->scores.erase(*i);
			}
		}
	}
	else {
		for (std::set<rclass_t>::const_iterator i = cell->classes.begin(); i != cell->classes.end(); ++i) {
			typed_relscore *&subscore = nt_score->scores[*i];
			assert(!subscore);
			if (*curr > 0) {
				const parse_cell *trailing_cell = links.cells[*curr - 1];
				assert(trailing_cell->classes.size() == 1);
				score_combiner tmp_score(measure_relation(links, *curr - 1, *curr, *trailing_cell->classes.begin(), *i));
				build_relscore_tree(subscore, active_score.add(tmp_score), links, first, curr + 1, last, *i);
			}
			else {
				build_relscore_tree(subscore, active_score, links, first, curr + 1, last, *i);
			}
			if (!subscore) {
				nt_score->scores.erase(*i);
			}
		}
	}

	if (nt_score->scores.empty()) {
		delete nt_score;
		score = 0;
	}
	else {
		VERBOSE(
			*verb_out << " nt scoring tree has entries for classes ";
			for (std::map<rclass_t, typed_relscore *>::const_iterator qq = nt_score->scores.begin(); qq != nt_score->scores.end(); ++qq) {
				*verb_out << qq->first << ' ';
			}
			*verb_out << std::endl;
		);

		score = nt_score;
	}
}
*/

template <typename T>
static void
extract_multiclass_trees(const parse_iterator &it, parse_links &links, const score_combiner &score, std::vector<expression_node> children,
                         std::vector<unsigned> &ranks, T first, T curr, T last)
{
	VERBOSE(*verb_out << "extract_multiclass_tree() at " << links << " with score " << score.score() << std::endl);

	unsigned firsti;
	if (curr == first) {
		firsti = 0;
	}
	else {
		firsti = *(curr - 1) + 1;
	}

	score_combiner final_score = score;
	if (curr == last) {
		assert(children.size() == firsti);
		
		for (unsigned i = firsti; i < links.children.size(); ++i) {
			VERBOSE(*verb_out << "extracting subtree ranked " << ranks[i] << " from " << links.children[i].nt->name << ':' << links.children[i].span->bits << std::endl);
			expression_node subtree = extract_tree_internal(it, links.children[i].span, links.cells[i], ranks[i]);
			if (!subtree) {
				return;
			}
			VERBOSE(*verb_out << "subtree is " << *subtree << std::endl);
			VERBOSE(*verb_out << " subtree score is " << subtree->score() << std::endl);
			final_score = final_score.add_struct(subtree->score());//subtree->score_combo());
			if (i > 0) {
				const expression_node &from_tree = children[i - 1];
				const expression_node &to_tree = subtree;
				double best_rel_score = 0.0;
				best_rel_score = compute_membership(links.P->rel, links.children[i-1].span, links.children[i].span, from_tree->rclasses, to_tree->rclasses);
				/*
				for (std::set<rclass_t>::const_iterator from_class = from_tree->rclasses.begin(); from_class != from_tree->rclasses.end(); ++from_class) {
					for (std::set<rclass_t>::const_iterator to_class = to_tree->rclasses.begin(); to_class != to_tree->rclasses.end(); ++to_class) {
						double local_rel_score = measure_relation(links, i - 1, i, from_tree, to_tree, *from_class, *to_class);
						if (local_rel_score >= best_rel_score) {
							best_rel_score = local_rel_score;
						}
					}
				}*/
				if (best_rel_score < SCORE_THRESHOLD) {
					return;
				}
				VERBOSE(*verb_out << "adding rel-score " << best_rel_score << " at " << i << std::endl);
				final_score = final_score.add_rel(best_rel_score);
			}
			
			children.push_back(subtree);
		}

		VERBOSE(*verb_out << "parse score is " << final_score.score() << std::endl);

		VERBOSE(
			*verb_out << " at end; building tree with ranks [ ";
			for (std::vector<unsigned>::const_iterator i = ranks.begin(); i != ranks.end(); ++i) {
				*verb_out << *i << ' ';
			}
			*verb_out << "]\n";
		);

		if (links.cells.size() ==3 && links.cells[1]->nt->name[0] == '!') {
			int a = 0;
		}
		expression_node tree = build_tree(it, links, ranks);
		if (!tree) {
			return;
		}

		if (!links.P->sbuild && children.size() == 1 && links.P->nt->sid == InvalidSemanticId) {
			tree->set_string_builder(children[0]->get_string_builder());
		}
		else if (links.P->sbuild) {
			tree->set_long_string("");
			tree->set_string_builder(links.P->sbuild);
		}

		/*
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
		}*/

		//assert(tree);
		VERBOSE(*verb_out << "  tree is " << tree->long_str() << std::endl);

		VERBOSE(*verb_out << "  final_score is " << final_score.score() << std::endl);
		if (final_score.score() > SCORE_THRESHOLD) {
			if (tree->noperations == 0) {
				for (unsigned i = 0; i < tree->nchildren(); ++i) {
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
			}
			tree->set_nt(links.P->nt);
			tree->head_span = links.children[links.P->head].span;
			tree->tail_span = links.children[links.P->tail].span;
			VERBOSE(*verb_out << "setting head span to " << tree->head_span->bits << "; tail span to " << tree->tail_span->bits << std::endl);
			links.cached_trees.insert(links_parse_result(tree, ranks));
		}
	}
	else {
		assert(children.size() == firsti);
		for (unsigned i = firsti; i < *curr; ++i) {
			expression_node subtree = extract_tree_internal(it, links.children[i].span, links.cells[i], ranks[i]);
			if (!subtree) {
				return;
			}
			final_score = final_score.add_struct(subtree->score());//->score_combo());

			if (i > 0) {
				const expression_node &from_tree = children[i - 1];
				const expression_node &to_tree = subtree;
				double best_rel_score = 0.0;
				best_rel_score = compute_membership(links.P->rel, links.children[i-1].span, links.children[i].span, from_tree->rclasses, to_tree->rclasses);

				/*for (std::set<rclass_t>::const_iterator from_class = from_tree->rclasses.begin(); from_class != from_tree->rclasses.end(); ++from_class) {
					for (std::set<rclass_t>::const_iterator to_class = to_tree->rclasses.begin(); to_class != to_tree->rclasses.end(); ++to_class) {
						double local_rel_score = measure_relation(links, i - 1, i, from_tree, to_tree, *from_class, *to_class);
						if (local_rel_score >= best_rel_score) {
							best_rel_score = local_rel_score;
						}
					}
				}*/
				if (best_rel_score < SCORE_THRESHOLD) {
					return;
				}
				VERBOSE(*verb_out << "adding rel-score " << best_rel_score << " at " << i << std::endl);
				final_score = final_score.add_rel(best_rel_score);
			}
			
			children.push_back(subtree);
		}

		unsigned orig_rank = ranks[*curr];
		parse_cell *cell = links.cells[*curr];
		expression_node tree = extract_tree_internal(it, links.children[*curr].span, cell, ranks[*curr]);
		while (tree) {
			const expression_node &to_tree = tree;
			children.push_back(to_tree);
			if (*curr > 0) {
				const expression_node &from_tree = children[*curr - 1];
				double best_rel_score = 0.0;
				best_rel_score = compute_membership(links.P->rel, links.children[*curr-1].span, links.children[*curr].span, from_tree->rclasses, to_tree->rclasses);
				/*
				for (std::set<rclass_t>::const_iterator from_class = from_tree->rclasses.begin(); from_class != from_tree->rclasses.end(); ++from_class) {
					for (std::set<rclass_t>::const_iterator to_class = to_tree->rclasses.begin(); to_class != to_tree->rclasses.end(); ++to_class) {
						double local_rel_score = measure_relation(links, *curr - 1, *curr, from_tree, to_tree, *from_class, *to_class);
						if (local_rel_score > best_rel_score) {
							best_rel_score = local_rel_score;
						}
					}
				}*/
						
				if (best_rel_score >= SCORE_THRESHOLD) {
					score_combiner tmp_score = final_score;
					tmp_score = tmp_score.add_struct(tree->score());//_combo());
					VERBOSE(*verb_out << "adding rel-score " << best_rel_score << " at " << *curr << std::endl);
					tmp_score = tmp_score.add_rel(best_rel_score);
					extract_multiclass_trees(it, links, tmp_score, children, ranks, first, curr + 1, last);
				}
			}
			else {
				extract_multiclass_trees(it, links,
				final_score.add_struct(to_tree->score())/*_combo())*/, children, ranks, first, curr + 1, last);
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
		VERBOSE(*verb_out << "setting head span to " << tree->head_span->bits << "; tail span to " << tree->tail_span->bits << std::endl);
		tree->set_context(ctx);
		tree->set_nt(links.children[0].nt);
		tree->nt_low = links.children[0].nt;
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
		for (std::vector<unsigned>::const_iterator i = links.multiclass_children.begin(); i != links.multiclass_children.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "]\nuniclass children = [ ";
		for (std::vector<unsigned>::const_iterator i = links.uniclass_children.begin(); i != links.uniclass_children.end(); ++i) {
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

	std::vector<unsigned> ranks(links.children.size(), 0);
	std::vector<expression_node> children;
	extract_multiclass_trees(it, links, score_combiner(), children, ranks, links.multiclass_children.begin(), links.multiclass_children.begin(), links.multiclass_children.end());

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


	for (std::vector<unsigned>::const_iterator j = links.uniclass_children.begin(); j != links.uniclass_children.end(); ++j) {
		++ranks[*j];
		links.next_index_sets.push_back(ranks);
		--ranks[*j];
	}

	return true;
}

const unsigned parse_cell::INVALID_SRC = static_cast<unsigned>(-1);


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
	VERBOSE(*verb_out << "preparing flags of " << cell << " to " << cell->prepared_flags << std::endl);

	cell->prev_src = 0;

	VERBOSE(*verb_out << "prepare_for_iteration at cell " << cell << " with " << cell->links.size() << " links\n" << std::flush);

	
	/*
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
	else */{
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
			/*
			if (i->P && span->locked_nt == i->P->nt && !span->locked_parse.empty() && i->P->nt->sid != InvalidSemanticId) {
				// This assertion will not be true when there are written operators
				// that are hidden as semantic types (eg plus, minus, etc.)
				// So we just ensure that all of the require children are
				// present in sequential order in the actual child
				// structure manifested in the ink.
				//assert(span->locked_parse.size() == i->children.size());
				std::list<bitvec>::const_iterator reqi = span->locked_parse.begin();
				std::vector<parseref>::const_iterator acti = i->children.begin();
				while (reqi != span->locked_parse.end()) {
					if (reqi->count_set_bits() > 0) {
						VERBOSE(*verb_out << "seeking child span " << *reqi << std::endl);
						VERBOSE(*verb_out << " comparing to " << acti->span->bits << std::endl);
						while (*reqi != acti->span->bits) {
							++acti;
							if (acti == i->children.end()) {
								goto do_not_prepare_links;
							}
							VERBOSE(*verb_out << " comparing to " << acti->span->bits << std::endl);
						}
						++acti;
					}
					++reqi;
				}
			}*/
				if (prepare_for_iteration(it, ctx, cell, *i, found_locked_nt || (i->P && span->locked_nt == i->P->nt))) {
					expression_node tree = extract_tree(it, *i, 0);
					if (tree) {
						tree->set_context(ctx);
						tree->set_cell(cell);
						tree->set_span(span);
						VERBOSE(*verb_out << "at " << *i << " with base index set, tree is " << tree.ptr() << " or " << tree->long_str() << " at " << tree->score() << " with bounds " << tree->span->bounds() << std::endl);
						//cell->cached_trees.push(cell_parse_result(tree, i - cell->links.begin()));
						cell->cached_trees.insert(cell_parse_result(tree, i - cell->links.begin()));
					}
//do_not_prepare_links:
//				;
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
	
	VERBOSE(*verb_out << "for cell " << cell << " the top-ranked tree is " << best.tree->long_str() << " at " << best.tree->score() << " with type " << best.tree->type() << std::endl);
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
extract_tree(const parse_iterator &it, parse_links &links, unsigned rank)
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
		for (std::list<std::vector<unsigned> >::const_iterator i = links.next_index_sets.begin(); i != links.next_index_sets.end(); ++i) {
			*verb_out << " [ ";
			for (std::vector<unsigned>::const_iterator j = i->begin(); j != i->end(); ++j) {
				*verb_out << *j << ' ';
			}
			*verb_out << "]\n";
		}
	);


	unsigned last_child;
	if (links.multiclass_children.empty()) {
		last_child = links.children.size();
	}
	else {
		last_child = links.multiclass_children[0];
	}


	for (std::list<std::vector<unsigned> >::iterator i = links.next_index_sets.begin(); i != links.next_index_sets.end(); ++i) {
		std::vector<unsigned> &ranks = *i;
		std::vector<expression_node> children;
		extract_multiclass_trees(it, links, score_combiner(), children, ranks, links.multiclass_children.begin(), links.multiclass_children.begin(), links.multiclass_children.end());
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

	std::vector<unsigned> ranks = best.src;
	for (std::vector<unsigned>::const_iterator i = links.uniclass_children.begin(); i != links.uniclass_children.end(); ++i) {
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
extract_tree_internal(const parse_iterator &it, ordered_segments *span, parse_cell *cell, unsigned rank)
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
		for (unsigned i = 0; i < cell->links.size(); ++i) {
			*verb_out << "  rank " << cell->next_best_ranks[i] << " at " << cell->links[i] << std::endl;
		}
	);

	expression_node next_tree;

	
	/*
	if (cell->nt->name[0] == '!') { // indicates matrix
		next_tree = extract_matrix_tree(cell->ranked_trees[0]->ctx, cell, span);
	}
	else */ {
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
		for (unsigned i = 0; i < node->nchildren(); ++i) {
			newnode->add_child(collapse_tree(node->mchild(i), node->mchild(i)->nt).ptr());
		}
	}
	return newnode;
}
*/

expression_node
extract_tree(const parse_iterator &it, ordered_segments *span, parse_cell *cell, unsigned rank, bool wrap_string)
{
	expression_node tree = extract_tree_internal(it, span, cell, rank);
	VERBOSE(*verb_out << "EXTRACTED tree " << (tree ? tree->long_str() : "null") << std::endl);
	if (!tree) {
		return NULL_EXPRESSION_NODE;
	}
	const char *s = tree->long_str();
	if (strncmp(s, MATHML_HEADER.c_str(), MATHML_HEADER.size()) == 0) {
		wrap_string = false;
	}
	expression_node node(new raw_expression_node(*tree, wrap_string));
	return node;//tree ? collapse_tree(tree, tree->nt) : 0;
}


void
dump_span_table(std::ostream &os, context *ctx)
{
	for (std::map<bitvec, ordered_segments *>::const_iterator i = ctx->spans->begin(); i != ctx->spans->end(); ++i) {
		os << i->first << " -> " << i->second << std::endl;
	}
}


void
dump_parse_table(std::ostream &os, context *ctx)
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

/*
static void
lock_existing_semantics(parse_cell *cell)
{
	VERBOSE(*verb_out << "locking existing semantics of " << cell->segs->bits << " to " << cell->nt->name << std::endl);

	cell->segs->locked_nt = cell->nt;

	for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		lock_existing_semantics(*i);
	}
}*/

static void
unlock_recursive(parse_cell *cell)
{
	cell->segs->locked_nt = 0;

	for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		unlock_recursive(*i);
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

	nonterminal_parses &parses = ctx->table[span];
	/*parse_cell *cell = parses[span->locked_nt];
	assert(cell);
	for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		lock_existing_semantics(*i);
	}*/

	for (nonterminal_parses::iterator i = parses.begin(); i != parses.end(); ++i) {
		if (i->second) {
			invalidate_cell(i->second);
		}
	}	

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
	locked_node->set_score(score_combiner(1.0));//(1.0));
	span->locked_node = locked_node;
	//span->locked_parse = node->parse_links;

	nonterminal_parses &parses = ctx->table[span];
	//parse_cell *cell = parses[node->nt_low];
	//assert(cell);
	/*for (std::list<parse_cell *>::iterator i = cell->parents.begin(); i != cell->parents.end(); ++i) {		
		lock_existing_semantics(*i);
	}*/
	
	for (nonterminal_parses::iterator i = parses.begin(); i != parses.end(); ++i) {
		if (i->second) {
			invalidate_cell(i->second);
		}
	}
	
	//invalidate_cell(cell);

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
	for (unsigned i = 0; i < parent->nchildren(); ++i) {
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

	if (nt) {
		parse_cell *cell = ctx->table[span][nt];
		assert(cell);
		unlock_recursive(cell);
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

	ctx->G = G;

	return 0;
}


int
update_spans_on_removal(context *ctx, unsigned pos)
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
append_edit_op(context *ctx, int type, unsigned pos)
{
	ctx->edits.push_front(edit_op(type, pos));
}
*/

}
