#include "binfmt.h"
#include "grammar.h"
#include "builder.h"
#include "stream-defs.h"
#include "error.h"
#include "ink-io.h"
#include "mathrecognizer-private.h"
#include "grouping.h"
#include "rect.h"
#include "MatrixAnalyzer.h"
#include "reco-types.h"

#include <string>
#include <istream>
#include <ostream>
#include <sstream>
#include <cstdlib>
#include <cmath>



namespace scg
{


const static long RECOGNIZER_VERSION = 8;
const static char VERSION_STRING[] = "LMRV";

/*
static void
save_prototype_id(writer &wr, const prototype_id &id) {
	wr.write((long)id.pkgid);
	wr.write(id.unicode);
	wr.write((long)id.seq);
}*/

static void
save_group(writer &wr, const group *grp) {
	wr.write(grp->bits);
	wr.write(grp->proximity_score);
	wr.write(grp->stack_score);
	wr.write(grp->feature_score);
	wr.write(grp->container_score);
	wr.write(grp->boosted_score);
	wr.write(grp->weight);
	wr.write(grp->pnil);
	wr.write(grp->tags);
	wr.write(grp->children.size());
	for (std::list<segment *>::const_iterator j = grp->children.begin(); j != grp->children.end(); ++j) {
		wr.write((*j)->pos);
	}
	/*if (grp->parents[0]) {
		assert(grp->parents[1]);
		bool hasparents = true;
		wr.write(hasparents);
		// NB: because of bitvec ordering, parents will already be written
		wr.write(grp->parents[0]->bits);
		wr.write(grp->parents[1]->bits);
	}
	else {
		assert(!grp->parents[1]);
		bool hasparents = false;
		wr.write(hasparents);
	}*/
	/*wr.write(grp->matches.size());
	for (std::vector<match_score>::const_iterator j = grp->matches.begin(); j != grp->matches.end(); ++j) {
		save_prototype_id(wr, j->proto->id);
		wr.write(j->feature_score);
		wr.write(j->matcher_score);
		wr.write(j->score);
	}*/
	//std::cout << "writing group for bits " << grp->bits << std::endl;
	wr.write(grp->matches.size());
	for (std::map<const nonterminal *, match_score>::const_iterator j = grp->matches.begin(); j != grp->matches.end(); ++j) {
		const match_score &match = j->second;
		//save_prototype_id(wr, match.proto->id);
		wr.write(match.S->unicode);
		wr.write(match.feature_score);
		wr.write(match.matcher_score);
		wr.write(match.score);
		wr.write(match.final_score);
		//std::cout << "  symbol " << match.S->info.name << " score " << match.final_score << std::endl;
	}
}

static void
save_stroke(writer &wr, const RawStroke &stk) {
	wr.write(stk.npoints);
	for (size_t j = 0; j < stk.npoints; ++j) {
		wr.write(stk.x[j]);
		wr.write(stk.y[j]);
		wr.write((long)(stk.time ? stk.time[j] : 0));
	}
}

static void
save_rect(writer &wr, const Rect<long> &r) {
	wr.write(r.left);
	wr.write(r.top);
	wr.write(r.right);
	wr.write(r.bottom);
}

int
math_recognizer_base::save(std::ostream &os) {
	writer wr(os);
	os.write(VERSION_STRING, sizeof(VERSION_STRING));
	wr.write(RECOGNIZER_VERSION);
	wr.write(segments.size());
	for (std::vector<segment *>::const_iterator i = segments.begin(); i != segments.end(); ++i) {
		const segment *seg = *i;
		save_stroke(wr, seg->stk);
	}

	wr.write(strokes.size());
	for (std::vector<stroke *>::const_iterator i = strokes.begin(); i != strokes.end(); ++i) {
		const stroke *stk = *i;
		save_stroke(wr, stk->input);
		wr.write(stk->parent->pos);
	}

	wr.write(known.size());
	for (std::vector<std::vector<size_t> >::const_iterator i = known.begin(); i != known.end(); ++i) {
		save_rect(wr, segments[i->front()]->bounds);
		wr.write(i->size());
		for (std::vector<size_t>::const_iterator j = i->begin(); j != i->end(); ++j) {
			wr.write(*j);
		}
	}

	if (groups) {
		wr.write(groups->size());
		for (std::map<bitvec, group *>::const_iterator i = groups->begin(); i != groups->end(); ++i) {
			const group *grp = i->second;
			save_group(wr, grp);
		}
	}
	else {
		wr.write((size_t)0);
	}

	wr.write(newgroups.size());
	for (std::list<group *>::const_iterator i = newgroups.begin(); i != newgroups.end(); ++i) {
		save_group(wr, *i);
	}

	const ordered_segments *allsegs = getsegs(bitvec(segments.size(), true));
	interpreter *rootintrpr = allsegs ? mkparser(GetMathGrammar()->root, getsegs(bitvec(segments.size(), true))) : 0;

	if (spans) {
		wr.write(spans->size());
		for (std::map<bitvec, ordered_segments *>::const_iterator i = spans->begin(); i != spans->end(); ++i) {
			wr.write(i->first);
		}
	}
	else {
		wr.write((size_t)0);
	}

	/*
	wr.write(updates.size());
	for (std::vector<lockspec>::const_iterator i = updates.begin(); i != updates.end(); ++i) {
		wr.write((long)i->type);
		wr.write((long)i->rev);
		wr.write((long)i->cancelat);
		if (i->cancelat == 0) {
			wr.write(i->intrpr);
		}
		else {
			wr.write((interpreter *)0);
		}
	}*/

	wr.write(rootintrpr);
	if (GetMathGrammar()->compexpr_root) {
		rootintrpr = allsegs ? mkparser(GetMathGrammar()->compexpr_root, getsegs(bitvec(segments.size(), true))) : 0;
		wr.write(rootintrpr);
	}
	if (GetMathGrammar()->partexpr_root) {
		rootintrpr = allsegs ? mkparser(GetMathGrammar()->partexpr_root, getsegs(bitvec(segments.size(), true))) : 0;
		wr.write(rootintrpr);
	}
	if (GetMathGrammar()->single_expr_root) {
		rootintrpr = allsegs ? mkparser(GetMathGrammar()->single_expr_root, getsegs(bitvec(segments.size(), true))) : 0;
		wr.write(rootintrpr);
	}
	if (GetMathGrammar()->multi_expr_root) {
		rootintrpr = allsegs ? mkparser(GetMathGrammar()->multi_expr_root, getsegs(bitvec(segments.size(), true))) : 0;
		wr.write(rootintrpr);
	}

	wr.write(tab.size());
	for (parse_table::const_iterator i = tab.begin(); i != tab.end(); ++i) {
		const nt_parse_table &nttab = i->second;
		wr.write(i->first->bits);
		assert(!spans || (*spans)[i->first->bits] == i->first);
		wr.write(nttab.size());
		for (nt_parse_table::const_iterator j = nttab.begin(); j != nttab.end(); ++j) {
			const nonterminal *nt = j->first;
			const production_parse_table &prodtab = j->second;
			wr.write(nt->index);
			wr.writeptr(prodtab.nt_parser);
			wr.write(prodtab.tab.size());
			for (production_parser_map::const_iterator k = prodtab.tab.begin(); k != prodtab.tab.end(); ++k) {
				size_t prod_index = std::find(nt->productions.begin(), nt->productions.end(), k->first) - nt->productions.begin();
				assert(prod_index < nt->productions.size());
				wr.write(prod_index);
				wr.writeptr(k->second);
			}
		}
	}
	return 0;
}

static void
restore_stroke(reader &re, RawStroke &stk) {
	re.read(stk.npoints);
	stk.x = new long[stk.npoints];
	stk.y = new long[stk.npoints];
	stk.time = new unsigned long[stk.npoints];
	for (size_t i = 0; i < stk.npoints; ++i) {
		re.read(stk.x[i]);
		re.read(stk.y[i]);
		re.read((long &)stk.time[i]);
	}
	//std::cout << stk << std::endl;
}

static void
restore_rect(reader &re, Rect<long> &rc) {
	re.read(rc.left);
	re.read(rc.top);
	re.read(rc.right);
	re.read(rc.bottom);
}

/*
static void
restore_prototype_id(reader &re, prototype_id &id) {
	long prox;
	re.read(prox);
	id.pkgid = prox;
	re.read(id.unicode);
	re.read(prox);
	id.seq = prox;
}*/

static void
restore_match_score(reader &re, match_score &ms) {
	//prototype_id id;
	//restore_prototype_id(re, id);
	unicode_char unicode;
	re.read(unicode);
	ms.S = symdb_findsymbol_unicode(unicode);
	re.read(ms.feature_score);
	re.read(ms.matcher_score);
	re.read(ms.score);
	re.read(ms.final_score);
	/*const symbols_db *db = global_symbols_db();
	ms.proto = db->find_prototype(id);
	if (!ms.proto) {
		THROW_ERROR(E_NOTFOUND, "restored match score refers to unknown prototype " << id);
	}
	ms.S = ms.proto->owner;*/
	//std::cout << "score for symbol " << ms.S->info.name << " = " << ms.final_score << std::endl;
}

static void
restore_group(reader &re, math_recognizer_base *rec, group *grp, std::vector<segment *> &segments) {
	re.read(grp->bits);
	re.read(grp->proximity_score);
	re.read(grp->stack_score);
	re.read(grp->feature_score);
	re.read(grp->container_score);
	re.read(grp->boosted_score);
	re.read(grp->weight);
	re.read(grp->pnil);
	re.read(grp->tags);
	grp->tags |= group_tags::LOADED_GROUP;

	size_t n;
	re.read(n);
	assert(n == grp->bits.count_set_bits());
	while (n--) {
		size_t pos;
		re.read(pos);
		if (pos >= segments.size()) {
			THROW_ERROR(E_INVALID, "while reading group, child pos " << pos << " is too large");
		}
		segment *child = segments[pos];
		grp->bounds = grp->children.empty() ? child->bounds : merge(grp->bounds, child->bounds);
		grp->children.push_back(child);
		child->parents.push_back(grp);
	}

	RawStrokeGroup G = build_stroke_group(grp->children.begin(), grp->children.end(), grp->bounds, grp->children.size());
	grp->strokes = G;

	/*re.read(n);
	grp->matches.resize(n);
	for (size_t i = 0; i < n; ++i) {
		restore_match_score(re, grp->matches[i]);
	}*/
	/*bool hasparents;
	re.read(hasparents);
	if (hasparents) {
		// NB: because of bitvec ordering, parents will already be loaded
		bitvec bits;
		re.read(bits);
		grp->parents[0] = (*rec->groups)[bits];
		assert(grp->parents[0]);
		re.read(bits);
		grp->parents[1] = (*rec->groups)[bits];
		assert(grp->parents[1]);
	}*/

	//std::cout << "group for bits " << grp->bits << std::endl;
	re.read(n);
	for (size_t i = 0; i < n; ++i) {
		match_score match;
		restore_match_score(re, match);
		grp->matches[match.S->nt] = match;
	}
}

int
math_recognizer_base::restore(std::istream &os) {
	for (size_t i = 0; i < sizeof(VERSION_STRING); ++i) {
		if (os.peek() != VERSION_STRING[i]) {
			ERR(E_OUTOFDATE, "");
			return ERROR_CODE();
		}
		os.get();
	}
	reader re(os);
	long version;
	re.read(version);
	if (version < RECOGNIZER_VERSION) {
		ERR(E_OUTOFDATE, "");
		return ERROR_CODE();
	}
	size_t n;
	re.read(n);
	segments.resize(n);
	avg_stroke_size = 0;
	for (size_t i = 0; i < n; ++i) {
		segments[i] = new segment;
		segments[i]->pos = i;
		restore_stroke(re, segments[i]->stk);
		segments[i]->bounds = segments[i]->stk.bounds();
		double sz = std::max(segments[i]->bounds.width(), segments[i]->bounds.height());
		avg_stroke_size += sz;
	}
	avg_stroke_size /= segments.size();

	re.read(n);
	//std::cout << "SCG_INK\n" << n << std::endl;
	strokes.resize(n);
	for (size_t i = 0; i < n; ++i) {
		strokes[i] = new stroke;
		restore_stroke(re, strokes[i]->input);
		size_t parpos;
		re.read(parpos);
		strokes[i]->parent = segments[parpos];
		segments[parpos]->children.push_back(strokes[i]);
		//std::cout << strokes[i]->input << std::endl;
		VERBOSE(*verb_out << "restore stroke with bounds " << strokes[i]->input.bounds() << std::endl);
	}
	//std::abort();

	re.read(n);
	known.resize(n);
	for (size_t i = 0; i < n; ++i) {
		Rect<long> bounds;
		restore_rect(re, bounds);
		size_t nsegs;
		re.read(nsegs);
		known[i].resize(nsegs);
		for (size_t j = 0; j < nsegs; ++j) {
			size_t idx;
			re.read(idx);
			segments[idx]->bounds = bounds;
			known[i][j] = idx;
		}
	}

	re.read(n);
	assert(groups);
	while (n--) {
		group *grp = new group;
		restore_group(re, this, grp, segments);
		(*groups)[grp->bits] = grp;
	}

	re.read(n);
	while (n--) {
		group *grp = new group;
		restore_group(re, this, grp, segments);
		newgroups.push_back(grp);
	}
	
	assert(spans);
	re.read(n);
	while (n--) {
		bitvec bits;
		re.read(bits);
		ordered_segments *&span = (*spans)[bits];
		assert(!span);
		span = new ordered_segments(segments, bits);
	}

	/*
	re.read(n);
	//updates.resize(n);
	for (size_t i = 0; i < n; ++i) {
		// TODO: updates may not be valid in case
		// of stroke deletions. Read in a lockspec
		// and add a new applylockspec() function
		// which will apply it or not based on cancelat,
		// but add it to the updates list regardless.
		lockspec ls;
		re.read(ls.type);
		re.read((long &)ls.rev);
		re.read((long &)ls.cancelat);
		re.read(&ls.intrpr, this);
		updates.push_back(ls);
		if (ls.type == LS_LOCK && ls.cancelat == 0) {
			assert(ls.intrpr);
			locks[ls.intrpr->span()] = ls;
		}
	}*/

	interpreter *rootintrpr;
	re.read(&rootintrpr, this); // root
	if (GetMathGrammar()->compexpr_root) {
		re.read(&rootintrpr, this); // complete expr root
	}
	if (GetMathGrammar()->partexpr_root) {
		re.read(&rootintrpr, this); // partial expr root
	}
	if (GetMathGrammar()->single_expr_root) {
		re.read(&rootintrpr, this); // single expr root
	}
	if (GetMathGrammar()->multi_expr_root) {
		re.read(&rootintrpr, this); // multi expr root
	}

	re.read(n);
	while (n--) {
		bitvec bits;
		re.read(bits);
		const ordered_segments *span = (*spans)[bits];
		assert(span);
		nt_parse_table &nttab = tab[span];
		size_t m;
		re.read(m);
		while (m--) {
			size_t ntindex;
			re.read(ntindex);
			if (ntindex >= GetMathGrammar()->nts.size()) {
				ERR(E_INVALID, "while reading table, nonterminal index " << ntindex << " is too large");
				return E_INVALID;
			}
			const nonterminal *nt = GetMathGrammar()->nts[ntindex];
			interpreter *intrpr;
			re.readptr(&intrpr);
			if (intrpr != (interpreter *)1) {
				nttab[nt].nt_parser = intrpr;
			}

			size_t nprodintrprs;
			re.read(nprodintrprs);
			while (nprodintrprs--) {
				size_t prod_index;
				re.read(prod_index);
				assert(prod_index < nt->productions.size());
				const production *P = nt->productions[prod_index];
				re.readptr(&intrpr);
				if (intrpr != (interpreter *)1) {
					nttab[nt].tab[P] = intrpr;
				}
			}
		}
	}

	return 0;
}

}
