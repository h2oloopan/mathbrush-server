#include "MathRecognizer.h"
#include "grammar.h"
#include "error.h"
#include "utils.h"
#include "parms.h"
#include "verb.h"
#include "relation.h"
#include "mathrecognizer-private.h"
#include "group.h"
#include "symbols.h"
#include "ptab.h"
#include "tpc-group.h"
#include "grammar-values.h"
#include "ordered-segments.h"
#include "expr-iter.h"
#include "surrogate-tree.h"
#include "intrpr.h"
#include "extern-iter.h"
#include "reco-types.h"
#include "stream-defs.h"
#include "strutils.h"
#include "replacements.h"

#include <algorithm>
#include <fstream>
#include <cassert>
#include <locale.h>

namespace scg
{

MathSymbol::MathSymbol() : symbol(0)
	{ }

MathSymbol::MathSymbol(unicode_char symbol_, const Rect<long> &bounds_)
	: symbol(symbol_), bounds(bounds_)
	{ }

inline bool
MathSymbol::operator==(const MathSymbol &rhs) const
	{ return symbol == rhs.symbol && bounds == rhs.bounds; }

inline bool
MathSymbol::operator!=(const MathSymbol &rhs) const
	{ return !(*this == rhs); }

static grammar *global_grammar = 0;
static ptab *global_ptab = 0;
ptab *get_ptab() { return global_ptab; }

int
RebuildGrammarData() {
	if (global_grammar) {
		rebuild_production_lengths(*global_grammar);
	}
	return 0;
}

#ifdef USING_TABLETPC

enum StkId_Typ {
	UNKSTK,
	TPCSTK,
	INVSTK,
	KNWNSTK
};

struct StrokeId {
	int typ;
	ptrdiff_t id;

	StrokeId() : typ(UNKSTK), id(0) { }
	StrokeId(int typ_, ptrdiff_t id_) : typ(typ_), id(id_) { }
};

static StrokeId
tablet_stkid(long id) {
	return StrokeId(TPCSTK, (ptrdiff_t)id);
}

static StrokeId
invalid_stkid() {
	return StrokeId(INVSTK, 0);
}

static StrokeId
known_stkid(const MathSymbol *ms) {
	return StrokeId(KNWNSTK, (ptrdiff_t)ms);
}

std::ostream &
operator<<(std::ostream &os, const StrokeId &sid) {
	os << sid.typ;
	switch (sid.typ) {
	case UNKSTK:
	case INVSTK:
	case KNWNSTK:
		break;
	case TPCSTK:
		os << ' ' << sid.id;
		break;
	}
	return os;
}

std::istream &
operator>>(std::istream &is, StrokeId &sid) {
	is >> sid.typ;
	switch (sid.typ) {
	case UNKSTK:
	case INVSTK:
	case KNWNSTK:
		break;
	case TPCSTK:
		is >> sid.id;
		break;
	}
	return is;
}

struct tablet_math_recognizer : public math_recognizer_base {
	explicit tablet_math_recognizer(std::istream &is);
	explicit tablet_math_recognizer(const scg::RawStrokeGroup &strokes_);
	tablet_math_recognizer(const MathSymbol *symbols, size_t n);
	explicit tablet_math_recognizer(scg::TPC_StrokeGroup &tpc_strokes);
	~tablet_math_recognizer() { ReleaseInk(); }

	ExpressionIterator *CreateIteratorForStrokesById(const long *strokes, size_t nstrokes, bool wrap_mathml = false);

	int Clear();

	int Save(std::ostream &os);

	int AddKnownSymbols(const MathSymbol *symbols, size_t n);
	int RemoveKnownSymbols(const MathSymbol *symbols, size_t n);


	int AddStrokes(IInkStrokeDisp **strokes, size_t nstrokes);
	int AddStrokes(IInkDisp *ink);
	int AddStrokes(const RawStrokeGroup &group);
	int AddStrokes(const RawStroke *strokes, size_t nstrokes);

	int RemoveStrokesById(const long *strokes, size_t nstrokes);
	int RemoveStrokesByIndex(const size_t *strokes, size_t nstrokes);

	int GetInk(IInkDisp **ink);
	void ReleaseInk();

public:
	std::vector<StrokeId> stroke_ids;
	IInkDisp *msink;
};


#else


struct scg_math_recognizer : public math_recognizer_base {
	explicit scg_math_recognizer(const scg::RawStrokeGroup &strokes_);
	explicit scg_math_recognizer(std::istream &is);
	scg_math_recognizer(const MathSymbol *symbols, size_t n);

	int Save(std::ostream &os);
};

#endif

static bool glob_autotrain = false;

DLLDECL void
EnableAutoTraining() {
	glob_autotrain = true;
}

DLLDECL void
DisableAutoTraining() {
	glob_autotrain = false;
}

DLLDECL bool
IsAutoTrainingEnabled() {
	return glob_autotrain;
}

extern int recognizer_initialize();
extern void recognizer_shutdown();

static std::string grammar_name;
const static std::string DEFAULT_GRAMMAR = "partial-expr";
DLLDECL void
SetGrammar(const char *gmr) {
	grammar_name = gmr;
}

void
math_recognizer_base::init() {
	avg_stroke_size = 0;
	groups = new std::map<bitvec, group *>;
	spans = new std::map<bitvec, ordered_segments *>;
}

math_recognizer_base::~math_recognizer_base() {
	clear();
}

math_recognizer_base::math_recognizer_base() {
	init();
}

math_recognizer_base::math_recognizer_base(std::istream &is) {
	init();
	int e = restore(is);
	if (e != 0) {
		throw e;
	}
}


math_recognizer_base::math_recognizer_base(const MathSymbol *symbols, size_t n) {
	init();
	AddKnownSymbols(symbols, n);
}



void
math_recognizer_base::clear() {
	while (!segments.empty()) {
		remove_segment(segments.size() - 1);
	}

	assert(strokes.empty());
	if (groups) {
		assert(groups->empty());
		delete groups;
		groups = 0;
	}

	while (!newgroups.empty()) {
		delete newgroups.front();
		newgroups.erase(newgroups.begin());
	}
	
	if (spans) {
		for (std::map<bitvec, ordered_segments *>::iterator i = spans->begin(); i != spans->end(); ++i) {
			delete i->second;
		}
		delete spans;
		spans = 0;
	}

	for (parse_table::iterator i = tab.begin(); i != tab.end(); ++i) {
		nt_parse_table &nttab = i->second;
		for (nt_parse_table::iterator j = nttab.begin(); j != nttab.end(); ++j) {
			delete j->second.nt_parser;
			production_parser_map &prodtab = j->second.tab;
			for (production_parser_map::iterator k = prodtab.begin(); k != prodtab.end(); ++k) {
				delete k->second;
			}
		}
	}
	tab.clear();
}

void
math_recognizer_base::Translate(long x, long y)
{
	for (std::vector<stroke *>::iterator i = strokes.begin(); i != strokes.end(); ++i) {
		(*i)->input.translate(x, y);
	}

	for (std::vector<segment *>::iterator i = segments.begin(); i != segments.end(); ++i) {
		segment *seg = *i;
		seg->stk.translate(x, y);
		seg->bounds.left += x;
		seg->bounds.top += y;
		seg->bounds.right += x;
		seg->bounds.bottom += y;
		VERBOSE(*verb_out << "translated segment " << seg << " to new bounds " << seg->bounds << std::endl);
	}

	if (groups) {
		for (std::map<bitvec, group *>::iterator i = groups->begin(); i != groups->end(); ++i) {
			group *grp = i->second;
			grp->bounds.left += x;
			grp->bounds.top += y;
			grp->bounds.right += x;
			grp->bounds.bottom += y;
			VERBOSE(*verb_out << "translated group " << grp << " to new bounds " << grp->bounds << std::endl);
		}
	}

	if (spans) {
		for (std::map<bitvec, ordered_segments *>::iterator i = spans->begin(); i != spans->end(); ++i) {
			i->second->translate(x, y);
		}
	}
}

int
math_recognizer_base::Save(std::ostream &os)
{
	try {
		return save(os);
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return e;
	}
}


ExpressionTree *
ConvertSurrogateTreeToExpressionTree(const SurrogateTree *st) {
	if (st->getType() == TERMINAL_EXPR) {
		std::string s = st->getTerminalValue();
		symbol *S = symdb_findsymbol_name(s);
		if (S) {
			wchar_t wc = (wchar_t)S->unicode;
			return new terminal_tree(0, S->name, std::wstring(&wc, 1), S->mathml, S->latex);
		}
		else {
			return new terminal_tree(0, s, str2wstr(s), s, s);
		}
	}

	const production *P = GetMathGrammar()->getcanonicalsidproduction(st->getType());
	if (!P) {
		return 0;
	}
	invented_tree *tree = new invented_tree(P);
	for (size_t i = 0; i < st->nchildren(); ++i) {
		ExpressionTree *ch = ConvertSurrogateTreeToExpressionTree(st->child(i));
		if (!ch) {
			delete tree;
			return 0;
		}
		tree->addchild(dynamic_cast<basic_tree *>(ch));
	}
	return tree;
}

static bool
iscomplete(const ExpressionTree *tree) {
	if (tree->type() == PLACEHOLDER_EXPR) {
		return false;
	}
	for (size_t i = 0; i < tree->nchildren(); ++i) {
		if (!iscomplete(tree->child(i))) {
			return false;
		}
	}
	return true;
}

void
math_recognizer_base::doresets() {
	while (!resets.empty()) {
		/*resetspec rs = resets.front();
		resets.pop();
		interpreter *intrpr = getparser(rs.nt, getsegs(rs.bits));
		assert(intrpr);
		if (rs.childi != ~0) {
			assert(rs.childi < intrpr->nchildren());
			intrpr = intrpr->child(rs.childi);
		}*/
		interpreter *intrpr = resets.front();
		resets.pop();
		intrpr->resetparents();
	}
}


const ExpressionTree *
math_recognizer_base::GetTopExpression() {
	bitvec bits(segments.size(), true);
	ordered_segments *&segs = (*spans)[bits];
	if (!segs) {
		segs = new ordered_segments(segments, bits);
	}

	doresets();

	interpreter *intrpr = mkparser(GetMathGrammar()->root, segs);
	if (!intrpr) {
		return 0;
	}
	external_iterator *it = new external_iterator(intrpr, false, true);
	const ExpressionTree *tree = it->next();
	delete it;
	/*if (iscomplete(tree)) {
		std::cout << "locking " << tree->str() << std::endl;
		tree->lock();
	}*/
	return tree;
}

size_t
math_recognizer_base::GetStrokes(const RawStroke **stkarr, size_t n)
{
	if (!stkarr) {
		return strokes.size();
	}

	std::vector<stroke *>::const_iterator end = strokes.begin() + std::min(n, strokes.size());
	const RawStroke **p = stkarr;
	for (std::vector<stroke *>::const_iterator i = strokes.begin(); i != end; ++i) {
		*(p++) = &(*i)->input;
	}

	return end - strokes.begin();
}

int
math_recognizer_base::RemoveStroke(const RawStroke *stk) {
	std::vector<size_t> indices;
	for (std::vector<stroke *>::const_iterator i = strokes.begin(); i != strokes.end(); ++i) {
		if (stk == &(*i)->input) {
			indices.push_back(i - strokes.begin());
			return RemoveStrokesByIndex(indices);
		}
	}
	return E_NOTFOUND;
}

ExpressionTree *
CreatePlaceholderExpression() {
	return mkplaceholder();
}

ExpressionTree *
CreateBlankExpression() {
	return mkblank();
}


ExpressionIterator *
math_recognizer_base::CreateDefaultSemanticIterator(SemanticId sid, bool wrap_mathml) {
	std::vector<size_t> s;
	for (size_t i = 0; i < segments.size(); i++) {
		s.push_back(i);
	}
	return CreateSemanticIteratorForStrokesByIndex(sid, &s[0], s.size(), wrap_mathml);
}

ExpressionIterator *
math_recognizer_base::CreateSemanticIteratorForStrokes(SemanticId sid, const RawStroke **strokes, size_t nstrokes, bool wrap_mathml) {
	if (!strokes) {
		return 0;
	}

	std::vector<size_t> indices;
	indices.reserve(nstrokes);

	for (const RawStroke **p = strokes; p != strokes + nstrokes; ++p) {
		std::vector<segment *>::const_iterator pseg;
		for (pseg = segments.begin(); pseg != segments.end(); ++pseg) {
			const segment *seg = *pseg;
			if (&seg->stk == *p) {
				indices.push_back(pseg - segments.begin());
				break;
			}
		}
		if (pseg == segments.end()) {
			ERR(E_NOTFOUND, "Some strokes could not be found in the recognition context.");
			return 0;
		}
	}

	return CreateSemanticIteratorForStrokesByIndex(sid, &indices[0], indices.size(), wrap_mathml);
}

ExpressionIterator *
math_recognizer_base::CreateSemanticIteratorForStrokesByIndex(SemanticId sid, const size_t *strokes, size_t nstrokes, bool wrap_mathml) {
	if (!strokes) {
		return 0;
	}

	ordered_segments *span = find_span(segments.size(), strokes, nstrokes, true);
	if (!span) {
		return 0;
	}

	doresets();
	multiplexor *mux = new multiplexor(this, 0, span, 2);
	const grammar *G = GetMathGrammar();
	for (std::vector<nonterminal *>::const_iterator i = G->nts.begin(); i != G->nts.end(); ++i) {
		const nonterminal *nt = *i;
		for (std::vector<production *>::const_iterator j = nt->productions.begin(); j != nt->productions.end(); ++j) {
			const production *prod = *j;
			if (prod->sid == sid) {
				interpreter *intrpr = mkprodintrpr(prod, span);
				if (intrpr) {
					mux->addparser(intrpr, false);
				}
			}
		}
	}
	if (mux->nchildren() == 0) {
		delete mux;
		return 0;
	}
	return new external_iterator(mux, true, wrap_mathml);
}



ExpressionIterator *
math_recognizer_base::CreateDefaultIterator(bool wrap_mathml) {
	std::vector<size_t> s;
	for (size_t i = 0; i < segments.size(); i++) {
		s.push_back(i);
	}
	return CreateIteratorForStrokesByIndex(&s[0], s.size(), wrap_mathml);
}


/*
void
math_recognizer_base::update(interpreter *intrpr, rev_t upto) const {
	rev_t maxrev = upto == 0 ? rev() : std::min(upto, rev());
	rev_t origrev = intrpr->rev();
	assert(maxrev >= origrev);
	bool upd = true;
	size_t actrev = intrpr->rev() + 1;
	while (actrev <= maxrev && upd) {
		lockspec ls = getupdate(actrev);
		switch (ls.type) {
		case LS_LOCK:
			if (ls.cancelat == 0) {
				upd = intrpr->lockat(this, ls.intrpr, ls.rev);
			}
			/*else {
				VERBOSE(*verb_out << "update(" << intrpr->nt()->name << ":" << intrpr->span()->bits << ") skipping cancelled lock " << ls.intrpr->nt()->name << ":" << ls.intrpr->span()->bits << std::endl);
				intrpr->advancerev(ls.rev);
			}* /
			break;
		case LS_UNLOCK:
			if (ls.cancelat == 0 || origrev >= ls.cancelat) {
				upd = intrpr->unlockat(this, ls.intrpr, ls.rev);
			}
			/*else {
				VERBOSE(*verb_out << "update(" << intrpr->nt()->name << ":" << intrpr->span()->bits << ") skipping cancelled unlock " << ls.intrpr->nt()->name << ":" << ls.intrpr->span()->bits << std::endl);
				intrpr->advancerev(ls.rev);
			}* /
			break;
		default:
			assert(false);
			break;
		}
		actrev++;
	}
}*/

external_iterator::external_iterator(interpreter *src_, bool rm_, bool wrap_mathml_)
	: src(src_), at(0), rm(rm_), wrap_mathml(wrap_mathml_) {
	//rec->update(src);
}

external_iterator::~external_iterator() {
	for (std::vector<replacement_entry>::iterator i = replacement_cache.begin(); i != replacement_cache.end(); ++i) {
		delete i->replacee;
	}
	if (rm) delete src;
}

void
external_iterator::release() {
	delete this;
}

const ExpressionTree *
external_iterator::next() {
	const basic_tree *tree = tree_from_interpretation(getnth(src, at));
	if (tree) ++at;
	return tree;
}

const ExpressionTree *
external_iterator::nth(size_t i) {
	if (i >= at) {
		return 0;
	}
	return tree_from_interpretation(src->nth(i));
}

const basic_tree *
external_iterator::tree_from_interpretation(interpretation *intrp) {
	if (!intrp) return 0;
	VERBOSE(*verb_out << "external_iterator got interpretation " << intrp->str() << std::endl);
	const basic_tree *tree = intrp->mktree();
	VERBOSE(*verb_out << "external_iterator got basic_tree " << tree->str() << std::endl);
	tree = postprocess_tree(tree);
	if (wrap_mathml) {
		tree->wrap();
	}
	return tree;
}

const basic_tree *
external_iterator::postprocess_tree(const basic_tree *tree) {
	for (std::vector<replacement_entry>::iterator i = replacement_cache.begin(); i != replacement_cache.end(); ++i) {
		if (*tree == *i->replacement) {
			tree = i->replacee;
			replacement_cache.erase(i);
			return tree;
		}
	}
	const basic_tree *repl = replace_tree(tree);
	if (repl) {
		VERBOSE(*verb_out << "replacing tree " << tree->ccstr() << " with tree " << repl->ccstr() << std::endl);
		replacement_cache.push_back(replacement_entry(tree, repl));
	}
	return repl ? repl : tree;
}

size_t
external_iterator::count() const {
	return at;
}

interpreter *
external_iterator::getsrc() {
	return src;
}

ExpressionIterator *
math_recognizer_base::CreateSubtreeIterator(const ExpressionTree *base, bool wrap_mathml) {
	const basic_tree *tree = reinterpret_cast<const basic_tree *>(base);
	/*const interpretation *intrp = tree->intrp();
	if (!intrp->P) {
		return 0;
	}
	interpreter *intrpr = mkparser(intrp->P->nt, intrp->span);*/
	interpreter *intrpr = tree->mkiter(this);
	return intrpr ? new external_iterator(intrpr, false, wrap_mathml) : 0;
}

ordered_segments *
math_recognizer_base::find_span(size_t nsegments, const size_t *strokes, size_t nstrokes, bool create)
{
	bitvec bits(nsegments, false);
	for (const size_t *s = strokes; s != strokes + nstrokes; ++s) {
		bits.set(*s, true);
	}

	std::map<bitvec, ordered_segments *>::iterator i;
	i = spans->find(bits);
	if (i != spans->end()) {
		return i->second;
	}
	
	if (create) {
		ordered_segments *span = new ordered_segments(segments, bits);
		spans->insert(i, std::make_pair(bits, span));
		return span;
	}
	return 0;
}

ExpressionIterator *
math_recognizer_base::CreateIteratorForExpression(const ExpressionTree *expr, bool ownexpr, bool wrap_mathml) {
	const basic_tree *stree = dynamic_cast<const basic_tree *>(expr);
	interpretation *intrp = new fixed_interpretation(stree, ownexpr);
	staticintrpr *intrpr = new staticintrpr(this, stree->nt(), stree->span());
	intrpr->addknown(intrp, true);
	return new external_iterator(intrpr, true, wrap_mathml);
}


ExpressionIterator *
math_recognizer_base::CreateIteratorForStrokes(const RawStroke **strokes, size_t nstrokes, bool wrap_mathml) {
	return CreateIteratorForStrokes(GetMathGrammar()->compexpr_root, strokes, nstrokes, wrap_mathml);
}

ExpressionIterator *
math_recognizer_base::CreateIteratorForStrokes(const nonterminal *nt, const RawStroke **strokes, size_t nstrokes, bool wrap_mathml) {
	if (!strokes) {
		return 0;
	}

	std::vector<size_t> indices;
	indices.reserve(nstrokes);

	VERBOSE(*verb_out << "CreateIteratorForStrokes: " << nstrokes << " strokes resolve to indices ");
	for (const RawStroke **p = strokes; p != strokes + nstrokes; ++p) {
		std::vector<segment *>::const_iterator pseg;
		for (pseg = segments.begin(); pseg != segments.end(); ++pseg) {
			const segment *seg = *pseg;
			if (&seg->stk == *p) {
				indices.push_back(pseg - segments.begin());
				VERBOSE(*verb_out << (pseg - segments.begin()) << ' ');
				break;
			}
		}
		if (pseg == segments.end()) {
			ERR(E_NOTFOUND, "Some strokes could not be found in the recognition context.");
			return 0;
		}
	}

	VERBOSE(*verb_out << std::endl);

	return CreateIteratorForStrokesByIndex(nt, &indices[0], indices.size(), wrap_mathml);
}


ExpressionIterator *
math_recognizer_base::CreateIteratorForStrokesByIndex(const nonterminal *nt, const size_t *strokes, size_t nstrokes, bool wrap_mathml) {
	if (!strokes) {
		return 0;
	}
   
	ordered_segments *span = find_span(segments.size(), strokes, nstrokes, true);
	if (!span) {
		return 0;
	}

	doresets();
   
	interpreter *intrpr = mkparser(nt, span);
	return intrpr ? new external_iterator(intrpr, false, wrap_mathml) : 0;
}


ExpressionIterator *
math_recognizer_base::CreateIteratorForStrokesByIndex(const size_t *strokes, size_t nstrokes, bool wrap_mathml) {
	return CreateIteratorForStrokesByIndex(GetMathGrammar()->compexpr_root, strokes, nstrokes, wrap_mathml);
}


extern segment *create_default_segment(stroke *stk);

int
math_recognizer_base::AddKnownSymbols(const MathSymbol *symbols, size_t n) {
	std::list<group *> newgroups;
	for (const MathSymbol *ms = symbols; ms != symbols + n; ++ms) {
		symbol *S = symdb_findsymbol_unicode(ms->symbol);
		if (!S || !S->nt) {// || !S->info.P) {
			THROW_ERROR(E_NOTFOUND, "a known symbol references the unknown or unused Unicode symbol " << std::hex << ms->symbol << std::dec);
		}

		const scg::nonterminal *nt = S->nt;
		
		size_t nstrokes = nt->min_length;
		std::vector<size_t> segs(nstrokes);
		group *grp = new group;
		grp->bounds = ms->bounds;
		bitvec bits(segments.size()+nstrokes, false);
		for (size_t i = 0; i < nstrokes; ++i) {
			add_segment();
			stroke *stk = new stroke;
			strokes.push_back(stk);
			segment *seg = create_default_segment(stk);
			seg->bounds = ms->bounds;
			segments.push_back(seg);
			segs.push_back(segments.size()-1);
			grp->children.push_back(seg);
			seg->parents.push_back(grp);
			bits.set(segments.size()-1);
		}
		
		grp->bits = bits;
		grp->proximity_score = 1.0;
		known.push_back(segs);
		grp->tags = group_tags::KNOWN_SYMBOL;

		match_score match;
		match.S = S;
		match.score = 10;
		grp->matches[match.S->nt] = match;//.push_back(match);
		(*groups)[grp->bits] = grp;
		//newgroups.push_back(grp);
	}
	
	/*while (!newgroups.empty()) {
		group *grp = newgroups.front();
		VERBOSE(*verb_out << "pushing results at new group " << grp << " at " << grp->bits << " to parse table\n");
		//push_group_matches_to_parse_table(&ctx, grp);
		newgroups.pop_front();
	}*/

	set_segment_positions();

	//VERBOSE(dump_ink_tree(*verb_out, &ctx));

	/*
	if (update_table) {
		return UpdateParseTable();
	}*/

	return 0;
}


int
math_recognizer_base::RemoveKnownSymbols(const MathSymbol *symbols, size_t n) {
	std::vector<size_t> indices;
	indices.reserve(n);

	for (const MathSymbol *ms = symbols; ms != symbols + n; ++ms) {
		for (std::vector<segment *>::const_iterator pseg = segments.begin(); pseg != segments.end(); ++pseg) {
			const segment *seg = *pseg;
			if (seg->bounds == ms->bounds) {
				indices.push_back(pseg - segments.begin());
			}
		}
	}

	std::sort(indices.begin(), indices.end(), std::greater<size_t>());
	return RemoveStrokesByIndex(indices);
}


int
math_recognizer_base::Clear()
{
	clear();
	groups = new std::map<bitvec, group *>;
	spans = new std::map<bitvec, ordered_segments *>;
	return 0;
}


int
math_recognizer_base::AddStrokes(const RawStrokeGroup &group) {
  return AddStrokes(group.strokes, group.nstrokes);
}


int
math_recognizer_base::RemoveStrokesByIndex(const size_t *strokes, size_t nstrokes) {
	std::vector<size_t> icopy(strokes, strokes + nstrokes);
	return RemoveStrokesByIndex(icopy);
}

int
math_recognizer_base::RemoveStrokesByIndex(const std::vector<size_t> &indices) {
	remove_strokes(&indices[0], indices.size());
	return 0;
}


int
math_recognizer_base::AddStrokes(const RawStroke *strokes, size_t nstrokes) {
	add_strokes(strokes, nstrokes);
  return 0;
}


#ifdef USING_TABLETPC

DLLDECL MathRecognizer *
CreateMathRecognizer() {
	TPC_StrokeGroup empty;
	try {
		return new tablet_math_recognizer(empty);
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return 0;
	}
}

int
tablet_math_recognizer::AddKnownSymbols(const MathSymbol *symbols, size_t n) {
	int e = math_recognizer_base::AddKnownSymbols(symbols, n);
	if (e == 0) {
		for (const scg::MathSymbol *ms = symbols; ms != symbols + n; ++ms) {
			stroke_ids.push_back(known_stkid(ms));
		}		
	}
	return e;
}

int
tablet_math_recognizer::RemoveKnownSymbols(const scg::MathSymbol *symbols, size_t n) {
	int e = math_recognizer_base::RemoveKnownSymbols(symbols, n);
	if (e == 0) {
		for (const scg::MathSymbol *ms = symbols; ms != symbols + n; ++ms) {
			for (std::vector<StrokeId>::iterator i = stroke_ids.begin(); i != stroke_ids.end(); ) {
				if (i->typ == KNWNSTK && i->id == (ptrdiff_t)ms) {
					i = stroke_ids.erase(i);
				}
				else {
					++i;
				}
			}
		}
		assert(strokes.size() == stroke_ids.size());
	}
	return e;
}

int
tablet_math_recognizer::RemoveStrokesByIndex(const size_t *indices, size_t n) {
	std::vector<size_t> icopy(indices, indices + n);
	std::sort(icopy.begin(), icopy.end(), std::greater<size_t>());
	for (std::vector<size_t>::const_iterator pi = icopy.begin(); pi != icopy.end(); ++pi) {
		stroke_ids.erase(stroke_ids.begin() + *pi);
	}

	return math_recognizer_base::RemoveStrokesByIndex(icopy);
}


tablet_math_recognizer::tablet_math_recognizer(const RawStrokeGroup &strokes_)
	: math_recognizer_base(), msink(0) {
	math_recognizer_base::AddStrokes(strokes_);
}

tablet_math_recognizer::tablet_math_recognizer(std::istream &is) : math_recognizer_base(is), msink(0) {
	stroke_ids.insert(stroke_ids.end(), strokes.size(), invalid_stkid());
}

tablet_math_recognizer::tablet_math_recognizer(const MathSymbol *symbols, size_t n) : math_recognizer_base(symbols, n), msink(0) { }

tablet_math_recognizer::tablet_math_recognizer(scg::TPC_StrokeGroup &tpc_strokes) : math_recognizer_base(), msink(0)
{
	if (tpc_strokes.nstrokes > 0) {
		AddStrokes(&tpc_strokes.strokes[0], tpc_strokes.nstrokes);
	}
}


int
tablet_math_recognizer::Clear()
{
	stroke_ids.clear();
	return math_recognizer_base::Clear();
}


int
tablet_math_recognizer::Save(std::ostream &os)
{
	int e = math_recognizer_base::Save(os);
	/*if (e == 0) {
		os << "STROKE_IDS " << stroke_ids.size() << ' ';
		for (std::vector<StrokeId>::const_iterator i = stroke_ids.begin(); i != stroke_ids.end(); ++i) {
			os << *i << ' ';
		}
		os << std::endl;
	}*/
	return e;
}


DLLDECL MathRecognizer *
CreateMathRecognizer(const MathSymbol *symbols, size_t n)
{
	try {
		return new tablet_math_recognizer(symbols, n);
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return 0;
	}
}


ExpressionIterator *
tablet_math_recognizer::CreateIteratorForStrokesById(const long *strokes, size_t nstrokes, bool wrap_mathml)
{
  if (stroke_ids.empty()) {
    return 0;
  }
    
  std::vector<size_t> indices;
    
  for (size_t i = 0; i < stroke_ids.size(); i++) {
	  const StrokeId &id = stroke_ids[i];
	  if (id.typ == TPCSTK
	   && std::find(strokes, strokes + nstrokes, (long)id.id) != strokes + nstrokes) {
          indices.push_back(i);
      }
  }
    
  return CreateIteratorForStrokesByIndex(&indices.front(), indices.size(), wrap_mathml);
}


int
tablet_math_recognizer::RemoveStrokesById(const long *strokes, size_t nstrokes) {
	VERBOSE(
		*verb_out << "removing ids ";
		for (const long *s = strokes; s != strokes + nstrokes; ++s) {
			*verb_out << *s << ' ';
		}
		*verb_out << " from ids ";
		for (std::vector<StrokeId>::const_iterator i = stroke_ids.begin(); i != stroke_ids.end(); ++i) {
			*verb_out << i->typ << ' ' << i->id;
		}
		*verb_out << std::endl;
	);
	
	std::vector<size_t> indices;
	indices.reserve(nstrokes);
	for (std::vector<StrokeId>::iterator i = stroke_ids.begin(); i != stroke_ids.end(); ) {
		if (i->typ == TPCSTK && std::find(strokes, strokes + nstrokes, (long)i->id) != strokes + nstrokes) {
			indices.push_back(i - stroke_ids.begin());
			i = stroke_ids.erase(i);
		}
		else {
			++i;
		}
	}
	
	return math_recognizer_base::RemoveStrokesByIndex(indices);
}


int
tablet_math_recognizer::AddStrokes(IInkStrokeDisp **strokes, size_t nstrokes) {
  TPC_StrokeGroup tpc_group(strokes, nstrokes, TPC_StrokeGroup::WEAK);
  RawStrokeGroup raw_group;
    
  int e = convert(raw_group, tpc_group);
  if (FAILURE(e)) {
    throw e;
  }

  for (size_t i = 0; i < nstrokes; i++) {
		long tmp;
    strokes[i]->get_ID(&tmp);
		stroke_ids.push_back(tablet_stkid(tmp));
  }
    
	return math_recognizer_base::AddStrokes(raw_group);
}

int
tablet_math_recognizer::AddStrokes(IInkDisp *ink) {
    HRESULT hr;
    IInkStrokes *strokes;
    if (FAILED(hr = ink->get_Strokes(&strokes))) {
        return MAKE_API_ERROR(hr);
    }
    
    long nstrokes;
    if (FAILED(hr = strokes->get_Count(&nstrokes))) {
        return MAKE_API_ERROR(hr);
    }

    IInkStrokeDisp **tpc_strokes = new IInkStrokeDisp *[nstrokes];
    
    for (long i = 0; i < nstrokes; i++) {
        if (FAILED(hr = strokes->Item(i, tpc_strokes + i))) {
            delete[] tpc_strokes;
		        return MAKE_API_ERROR(hr);
        }
    }

    return AddStrokes(tpc_strokes, nstrokes);
}

int
tablet_math_recognizer::AddStrokes(const RawStrokeGroup &group) {
	int e = math_recognizer_base::AddStrokes(group);
	if (e == 0 && strokes.size() != stroke_ids.size()) {
		stroke_ids.insert(stroke_ids.end(), group.nstrokes, invalid_stkid());
		assert(strokes.size() == stroke_ids.size());
	}
	return e;
}

int
tablet_math_recognizer::AddStrokes(const RawStroke *stkarr, size_t nstrokes) {
	int e = math_recognizer_base::AddStrokes(stkarr, nstrokes);
	if (e == 0 && strokes.size() != stroke_ids.size()) {
		stroke_ids.insert(stroke_ids.end(), nstrokes, invalid_stkid());
		assert(strokes.size() == stroke_ids.size());
	}
	return e;
}


int
tablet_math_recognizer::GetInk(IInkDisp **ink)
{
	*ink = make_empty_ink_object();
	IInkDisp *real_ink = *ink;
	if (!real_ink) {
		return E_SYSTEM;
	}
	
	stroke_ids.clear();
	for (std::vector<stroke *>::const_iterator i = strokes.begin(); i != strokes.end(); ++i) {
		IInkStrokeDisp *stroke;
		int e = convert((*i)->input, real_ink, &stroke);
		if (FAILURE(e)) {
			real_ink->Release();
			return e;
		}
		
		long id;
		HRESULT hr = stroke->get_ID(&id);
		if (FAILED(hr)) {
			return MAKE_API_ERROR(hr);
		}
	
		stroke_ids.push_back(tablet_stkid(id));
	}
	
	msink = real_ink;
	
	return 0;
}


void
tablet_math_recognizer::ReleaseInk()
{
	if (msink) {
		msink->Release();
		msink = 0;
	}
}


DLLDECL MathRecognizer *
CreateMathRecognizer(IInkDisp *ink)
{
	TPC_StrokeGroup empty;
	try {
		MathRecognizer *rec = new tablet_math_recognizer(empty);
		rec->AddStrokes(ink);
		return rec;
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return 0;
	}
}

DLLDECL MathRecognizer *
CreateMathRecognizer(std::istream &is) {
	try {
		return new tablet_math_recognizer(is);
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return 0;
	}
}



#else

DLLDECL MathRecognizer *
CreateMathRecognizer() {
	try {
		return new scg_math_recognizer(RawStrokeGroup());
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return 0;
	}
}


DLLDECL MathRecognizer *
CreateMathRecognizer(std::istream &is) {
	MathRecognizer *rec;
	try {
		rec = new scg_math_recognizer(is);
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return 0;
	}

	std::string label;
	size_t n;
	is >> label >> n;
	CHECK_ISTREAM_BASIC(is);

	if (label != "STROKE_IDS") {
		ERR(E_INVALID, "saved context is invalid");
		return 0;
	}

	long dummy;
	while (n--) {
		is >> dummy;
		CHECK_ISTREAM_BASIC(is);
	}

	return rec;
}


DLLDECL MathRecognizer *
CreateMathRecognizer(const MathSymbol *symbols, size_t n) {
	try {
		return new scg_math_recognizer(symbols, n);
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return 0;
	}
}





scg_math_recognizer::scg_math_recognizer(std::istream &is) : math_recognizer_base(is) { }

scg_math_recognizer::scg_math_recognizer(const MathSymbol *symbols, size_t n) : math_recognizer_base(symbols, n) { }

scg_math_recognizer::scg_math_recognizer(const scg::RawStrokeGroup &strokes_) : math_recognizer_base()
{
	AddStrokes(strokes_);
}



int
scg_math_recognizer::Save(std::ostream &os)
{
	int e = math_recognizer_base::Save(os);
	os << "STROKE_IDS 0\n";
	return e;
}


#endif

DLLDECL MathRecognizer *
CreateMathRecognizer(const RawStrokeGroup &strokes)
{
	try {
#ifdef USING_TABLETPC
		return new tablet_math_recognizer(strokes);
#else
		return new scg_math_recognizer(strokes);
#endif
	}
	catch (int e) {
		ENSURE_ERROR(e);
		return 0;
	}
}

DLLDECL const grammar *
GetMathGrammar()
{
	return global_grammar;
}


static size_t initialization_counter = 0;

// helper struct for auto-decrementing refcount on initialization errors
struct counter_wrapper {
	counter_wrapper(size_t *count) : count(count) { ++(*count); }
	~counter_wrapper() { if (count) --(*count); }
	void release() { count = 0; }
	size_t *count;
};

DLLDECL int
InitializeRecognizer() {
	if (initialization_counter > 0) {
		++initialization_counter;
		return 0;
	}
	counter_wrapper refcount(&initialization_counter);

	int e;
	initialize_grammar_typemap();
	if (FAILURE(e = RestoreRegisteredParms())) return e;
	if (FAILURE(e = initialize_relations())) return e;
	//if (FAILURE(e = rebuild_global_symbols_db())) throw e;
	if (FAILURE(e = symdb_init())) return e;
	if (FAILURE(e = recognizer_initialize())) return e;
	std::string path;
	if (FAILURE(e = GetProfilePath(path))) return e;
	std::string gfile;
	if (grammar_name.empty()) {
		grammar_name = DEFAULT_GRAMMAR;
	}
	gfile = path + "/" + grammar_name + ".grammar";
	std::ifstream ifs(gfile.c_str());
	if (ifs.fail() || ifs.bad() || !ifs.is_open()) {
		return E_NOTFOUND;
	}
	
	//ptab_staticinit();

	global_grammar = new grammar;
	try {
		ifs >> *global_grammar;
	}
	catch (error &e) {
		return e.code;
	}
	catch (int e) {
		return e;
	}

	/*ptab_rectify(global_grammar);

	ifs.close();
	ifs.open((path + "/ptab").c_str());
	if (ifs.is_open()) {
		global_ptab = mkptab(global_grammar);
		e = ptab_read(global_ptab, ifs);
		if (FAILURE(e)) {
			throw e;
		}
	}*/

	if (FAILURE(e = InitializeMatrixRecognizer(global_grammar))) return e;
	if (FAILURE(e = initcanonicalsids(global_grammar))) return e;

	/*
	std::vector<std::string> rm;
	for (symbol *S = symdb_firstsymbol(); S; S = symdb_nextsymbol(S)) {
		if (!S->nt) {
			rm.push_back(S->name);
			if (pt) {
				ptab_remove(pt, S->sid);
			}
		}
	}

	for (std::vector<std::string>::const_iterator i = rm.begin(); i != rm.end(); ++i) {
		symdb_removesym(*i, false);
	}*/
	refcount.release(); // don't decrement the ref count now, we're all initialized
	return 0;
}


DLLDECL void
ShutdownRecognizer() {
	if (initialization_counter == 0) return;
	if (--initialization_counter > 0) return;

	if (global_ptab) {
		rmptab(global_ptab);
		global_ptab = 0;
	}
	ShutdownMatrixRecognizer();
	clearuniqproductions();
	delete global_grammar;
	global_grammar = 0;
	recognizer_shutdown();
	symdb_shutdown();
	destroy_relations();
	ReleaseParameters();
	destroy_grammar_typemap();
}


}
