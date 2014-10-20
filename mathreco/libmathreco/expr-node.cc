#include "expr-node.h"
#include "grammar.h"
#include "verb.h"
#include "mathrecognizer-private.h"
#include "builder.h"
#include "ordered-segments.h"
#include "grammar.h"
#include "expr-iter.h"
#include "extern-iter.h"
#include "reco-types.h"
#include "strutils.h"

namespace scg
{

interpretation::interpretation(interpreter *src_, const production *P_, const ordered_segments *segs)
	: src(src_), P(P_), span(segs), head(this), tail(this), haslong(0), rclass(1 << AGGREGATE_CLASS),
	  nterms(0), relscore(0.0), scorebias(0), hint(0), grp(0) {
	//VERBOSE(*verb_out << "new " << this << std::endl);
	//std::cout << "new " << this << std::endl;
}

interpretation::~interpretation() {
	/*
	VERBOSE(*verb_out << "del " << this << " which is " << ccstr() << std::endl);
	for (size_t i = 0; i < children.size(); ++i) {
		if (rmflags[i]) {
			delete children[i];
		}
	}*/
	//std::cout << "del " << this << std::endl;
}

recoscore
interpretation::unnorm_score() const {
	recoscore sc;
	if (grp) {
		assert(nchildren() == 0);
		std::map<const nonterminal *, match_score>::const_iterator i = grp->matches.find(nt());
		if (i != grp->matches.end()) {
			sc = i->second.final_score;
			if (grp->pnil > 0) {
				sc.score = std::log(1.0 - grp->pnil) + sc.score - std::log(grp->pnil);
				sc.score *= 2.0;
				VERBOSE(*verb_out << "terminal " << nt()->name << ":" << span->bits << " has score " << i->second.final_score.score << " and pnil " << grp->pnil << " -> " << sc.score << std::endl);
			}
		}
		//sc = (i == grp->matches.end()) ? 0 : i->second.final_score;
		VERBOSE3(*verb_out << "  terminal: " << sc << std::endl);
		return sc;
	}
	for (size_t i = 0; i < nchildren(); ++i) {
		recoscore chsc = child(i)->unnorm_score();
		VERBOSE3(*verb_out << "  child score " << chsc << std::endl);
		//sc *= chsc;
		sc += chsc;
	}
	VERBOSE3(*verb_out << " * relscore " << relscore << std::endl);
	//sc *= relscore;
	sc += relscore;
	sc += scorebias;
	sc.hint += hint;
	return sc;
}

recoscore
interpretation::score() const {
	//assert(nterms > 0); // breaks dummy interpretations in rscore pre-screening
	VERBOSE3(*verb_out << "interpretation::score() @ " << nt()->name << ":" << span->bits << std::endl);
	recoscore sc = unnorm_score();
	VERBOSE3(*verb_out << " unnorm score " << sc << " with " << nterms << " terminals -> ");
	//sc = std::pow(sc, 1.0/(2*nterms-1));
	//sc = sc/(2*nterms-1);
	VERBOSE3(*verb_out << sc << std::endl);
	return sc;
}

math_recognizer_base *
interpretation::ctx() const {
	return src->ctx();
}

const nonterminal *
interpretation::nt() const {
	return P ? P->nt : (src ? src->nt() : 0);//P ? P->nt : 0;
}

const Rect<long>
interpretation::bounds() const {
	return span->bounds();
}

void
interpretation::addchild(interpretation *child, double relsc, bool rm) {
	addmetachild(child, rm);
	//score_ = score_.add(child->score_);
	//relscore *= relsc;
	relscore += relsc;
	nterms += child->nterms;
}

void
interpretation::addmetachild(interpretation *child, bool rm) {
	assert(child != this);
	children.push_back(child);
	rmflags.push_back(rm);
}

void
interpretation::clrstrs() {
	raw.clear();
	wraw.clear();
}

const char *
interpretation::str() const {
	return ccstr().c_str();
}

const std::string &
interpretation::ccstr() const {
	if (raw.empty()) {
		std::stringstream ss;
		ss.precision(2);
		bool printable = P && P->sid != InvalidSemanticId && P->sid != TERMINAL_EXPR;
		if (printable) {
			//ss << '(' << (P && P->nt ? P->nt->name : "nil") << " ";
			ss << '(' << sid_to_string(P->sid) << ' ';
		}
		for (std::vector<interpretation *>::const_iterator i = children.begin(); i != children.end(); ++i) {
			if (i != children.begin()) {
				ss << ' ';
			}
			ss << (*i)->str();
		}
		if (printable) {
			ss << ")";
		}
		/*if (isterminal(this) || nchildren() > 1) {
			ss << ":" << score();
		}*/
		raw = ss.str();
	}
	return raw;
}

const wchar_t *
interpretation::wstr() const {
	return ccwstr().c_str();
}

const std::wstring &
interpretation::ccwstr() const {
	if (wraw.empty()) {
		//std::wstringstream ss;
		//ss << L"(" << wraw;
		for (std::vector<interpretation *>::const_iterator i = children.begin(); i != children.end(); ++i) {
			/*if (i != children.begin()) {
				ss << L" ";
			}*/
			//ss << (*i)->wstr();
			wraw.append((*i)->ccwstr());
		}
		//ss << L")";
		//wraw = ss.str();
	}
	return wraw;
}

static const production *
readproduction(reader &re) {
	const production *P = 0;
	size_t ntindex, Pindex;
	re.read(ntindex);
	re.read(Pindex);
	if (ntindex) {
		const grammar &G = *GetMathGrammar();
		if (ntindex - 1 >= G.nts.size()) {
			THROW_ERROR(E_INVALID, "while loading interpretation, nonterminal index " << ntindex - 1 << " is too large");
		}
		const nonterminal *nt = G.nts[ntindex - 1];
		if (Pindex >= nt->productions.size()) {
			THROW_ERROR(E_INVALID, "while loading interpretation, production index " << Pindex << " is too large for " << nt->name);
		}
		P = nt->productions[Pindex];
	}
	else {
		P = getuniqproduction((int)Pindex);
	}
	return P;
}

interpretation::interpretation(reader &re, math_recognizer_base *rec) : haslong(0) {
	re.read(&src, rec);
	P = readproduction(re);
	if (!P) {
		re.read(raw);
		re.read(wraw);
	}
	bool grpprox;
	re.read(grpprox);
	if (grpprox) {
		bitvec grpbits;
		re.read(grpbits);
		grp = (*rec->groups)[grpbits];
		assert(grp);
	}
	long prox;
	re.read(prox);
	nterms = prox;
	re.read(relscore);
	re.read(hint);
	re.read(scorebias);
	size_t nchil;
	re.read(nchil);
	children.resize(nchil);
	rmflags.resize(nchil);
	for (size_t i = 0; i < nchil; ++i) {
		re.read(&children[i], rec);
		bool prox;
		re.read(prox);
		rmflags[i] = prox;
	}
	bitvec bits;
	re.read(bits);
	if (bits.size() > 0) {
		span = rec->getsegs(bits);
		if (!span) {
			THROW_ERROR(E_INVALID, "while loading interpretation, span not found at " << bits);
		}
	}
	re.read((interpretation **)&head, rec);
	if (!head) head = this;
	re.read((interpretation **)&tail, rec);
	if (!tail) tail = this;
	re.read(prox);
	rclass = prox;
}

void
interpretation::write(writer &wr) const {
	//VERBOSE(*verb_out << "wr " << this << std::endl);
	wr.write(src);
	wr.write((P && P->nt) ? P->nt->index + 1 : 0);
	wr.write(P ? P->index : 0);
	if (!P) {
		wr.write(raw);
		wr.write(wraw);
	}
	if (grp) {
		wr.write(true);
		wr.write(grp->bits);
	}
	else {
		wr.write(false);
	}
	wr.write((long)nterms);
	wr.write(relscore);
	wr.write(hint);
	wr.write(scorebias);
	wr.write(nchildren());
	for (size_t i = 0; i < nchildren(); ++i) {
		wr.write(children[i]);
		wr.write(rmflags[i]);
	}
	wr.write(span ? span->bits : bitvec::EMPTY);
	if (head == this) wr.write((long)0);
	else wr.write(head);
	if (tail == this) wr.write((long)0);
	else wr.write(tail);
	wr.write((long)rclass);
}

bool
isterminal(const interpretation *intrp) {
	return (!intrp->P || intrp->P->sid == TERMINAL_EXPR);
}

bool
derivesterminal(const interpretation *intrp) {
	return isterminal(intrp) || (intrp->nchildren() == 1 && derivesterminal(intrp->child(0)));
}

const interpretation *
getsemanticintrp(const interpretation *intrp) {
	if (!intrp) {
		return 0;
	}
	if (!intrp->P || intrp->P->sid != InvalidSemanticId || intrp->P->rhs.size() > 1 || isterminal(intrp)) {
		return intrp;
	}
	if (intrp->nchildren() != 1) {
		abort();
	}
	return getsemanticintrp(intrp->child(0));
}

const interpretation *
getlockintrp(const interpretation *intrp) {
	if (!intrp->P || (intrp->nchildren() == 0 && intrp->P->sid == TERMINAL_EXPR) || intrp->nchildren() > 1) {
		return intrp;
	}
	return getlockintrp(intrp->child(0));
}


struct tree_accessor : public subtree_accessor {
public:
	explicit tree_accessor(const interpretation *root_) : root(root_), used(root->nchildren(), false), subtrees(root->nchildren(), 0) {	}
	~tree_accessor() {
		for (size_t i = 0; i < root->nchildren(); ++i) {
			if (!used[i]) {
				delete subtrees[i];
			}
		}
	}

	void markused(size_t i) { used[i] = true; }

	basic_tree *operator[](size_t i) {
		if (!subtrees[i]) {
			subtrees[i] = root->child(i)->mktree();
		}
		return subtrees[i];
	}

private:
	const interpretation *root;
	std::vector<basic_tree *> subtrees;
	std::vector<bool> used;
};

basic_tree *
interpretation::mktree() const {
	const interpretation *sintrp = getsemanticintrp(this);
	basic_tree *tree;
	if (!sintrp->P) {
		tree = new terminal_tree(this, sintrp->ccstr(), sintrp->ccwstr(), sintrp->ccstr(), sintrp->ccstr());
	}
	else {
		assert(sintrp->P->tbuild);
		tree_accessor acc(sintrp);
		tree = sintrp->P->tbuild->build(this, acc);
	}
	return tree;
}

bool
interpretation::haslongform() const {
	if (haslong) {
		return haslong > 0;
	}
	for (std::vector<interpretation *>::const_iterator i = children.begin(); i != children.end(); ++i) {
		if ((*i)->haslongform()) {
			haslong = 1;
			return true;
		}
	}
	haslong = -1;
	return false;
}

interpreter *
interpretation::longformiter() const {
	if (haslongform()) {
		std::vector<interpreter *> subintrprs;
		for (std::vector<interpretation *>::const_iterator i = children.begin(); i != children.end(); ++i) {
			subintrprs.push_back((*i)->longformiter());
		}
		return new treeparser(ctx(), P, span, subintrprs);
	}
	else {
		staticintrpr *intrpr = new staticintrpr(ctx(), nt(), span);
		intrpr->addknown(const_cast<interpretation *>(this), false);
		return intrpr;
	}
}


interpretation *
mkterminal(interpreter *src, const symbol *S, group *grp, const ordered_segments *segs) {
	interpretation *intrp = new interpretation(src, /*S->P*/0, segs);
	//intrp->set_score(score_combiner(score));
	intrp->nterms = 1;
	intrp->grp = grp;
	intrp->raw = S->name;//value;
	intrp->wraw = std::wstring(1, S->unicode);
	//intrp->headspan = intrp->tailspan = segs;
	intrp->rclass = S->rclass;
	
	VERBOSE(*verb_out << "created TERMINAL NODE " << intrp << " for terminal " << S->name/*value*/ << std::endl);
	return intrp;
}

interpretation *
mkterminal(const std::string &s) {
	interpretation *intrp = new interpretation(0, 0, 0);
	//intrp->set_score(score_combiner(score));
	intrp->raw = s;
	intrp->wraw = str2wstr(s);
	return intrp;

}

const char interpretation::ID = 0;
const char fixed_interpretation::ID = 1;

fixed_interpretation::fixed_interpretation(const basic_tree *tree_, bool own_)
	: interpretation(0, 0, 0), tree(tree_), own(own_) {
	//set_score(score_combo().addterm(tree->score()));
}

fixed_interpretation::fixed_interpretation(reader &re, math_recognizer_base *rec)
	: interpretation(0, 0, 0) {
	re.read(const_cast<basic_tree **>(&tree), rec);
	re.read(own);
}

fixed_interpretation::~fixed_interpretation() {
	if (own) {
		delete tree;
	}
}

void
fixed_interpretation::write(writer &wr) const {
	wr.write(tree);
	wr.write(own);
}

const std::string &
fixed_interpretation::ccstr() const { return tree->ccstr(); }
const std::wstring &
fixed_interpretation::ccwstr() const { return tree->ccwstr(); }
basic_tree *
fixed_interpretation::mktree() const { return tree->mkcopy(); }

const char parsed_tree::ID = 1;
const char invented_tree::ID = 2;
const char terminal_tree::ID = 3;

bool
operator==(const basic_tree &lhs, const basic_tree &rhs) {
	if (lhs.type() != rhs.type()) return false;
	if (lhs.type() == TERMINAL_EXPR) return lhs.ccstr() == rhs.ccstr();
	if (lhs.nchildren() != rhs.nchildren()) return false;
	for (size_t i = 0; i < lhs.nchildren(); ++i) {
		if (*lhs.child(i) != *rhs.child(i)) return false;
	}
	return true;
}
bool operator!=(const basic_tree &lhs, const basic_tree &rhs) { return !(lhs == rhs); }

basic_tree::~basic_tree() {
	for (size_t i = 0; i < children.size(); ++i) {
		if (rmflags[i]) {
			delete children[i];
		}
	}
}

void
basic_tree::replace_child(size_t i, basic_tree *child) {
	if (i >= children.size()) {
		throw E_INVALID;
	}
	if (children[i]) children[i]->release();
	children[i] = child;
}

basic_tree::basic_tree(reader &re, math_recognizer_base *rec) {
	size_t n;
	re.read(n);
	children.resize(n);
	rmflags.resize(n);
	for (size_t i = 0; i < n; ++i) {
		re.read(const_cast<basic_tree **>(&children[i]), rec);
		bool rm;
		re.read(rm);
		rmflags[i] = rm;
	}
}

void
basic_tree::write(writer &wr) const {
	wr.write(children.size());
	for (size_t i = 0; i < children.size(); ++i) {
		wr.write(children[i]);
		wr.write(rmflags[i]);
	}
}

basic_tree *
basic_tree::mkcopy() const {
	basic_tree *cp = mkcopy_();
	for (std::vector<basic_tree *>::const_iterator i = children.begin(); i != children.end(); ++i) {
		cp->addchild((*i)->mkcopy());
	}
	return cp;
}

const char *
basic_tree::str() const {
	return ccstr().c_str();
}

const std::string &
basic_tree::ccstr() const {
	if (raw.empty()) {
		std::stringstream ss;
		ss << '(' << sid_to_string(type()) << ' ';
		for (size_t i = 0; i < nchildren(); ++i) {
			if (i != 0) {
				ss << ' ';
			}
			ss << child(i)->str();
		}
		ss << ')';
		raw = ss.str();
	}	    	    
	return raw;
}

const wchar_t *
basic_tree::wstr() const {
	return ccwstr().c_str();
}

const std::wstring &
basic_tree::ccwstr() const {
	if (wraw.empty()) {
		/*std::wstringstream ss;
		ss << L"(" << type();
		for (size_t i = 0; i < nchildren(); ++i) {
			if (i < nchildren() - 1) {
				ss << L" ";
			}
			ss << child(i)->wstr();
		}
		ss << L")";
		wraw = ss.str();*/
		for (size_t i = 0; i < nchildren(); ++i) {
			std::wstring ws = child(i)->ccwstr();
			wraw.append(ws);
		}
	}
	return wraw;
}

const char *
basic_tree::long_str() const {
	if (mathml.empty()) {
		VERBOSE(*verb_out << "basic_tree: building MathML string for tree " << ccstr() << std::endl);
		const string_builder *sb = sbuilder();
		if (sb) {
			mathml = sb->build(this);
		}
		VERBOSE(*verb_out << " got " << mathml << std::endl);
	}
	return mathml.c_str();
}

/*std::string
ExpressionTree::portable_str() const {
	if (type() == TERMINAL_EXPR) {
		return str();
	}
	std::stringstream ss;
	ss << '(' << sid_to_string(type());
	for (size_t i = 0; i < nchildren(); ++i) {
		ss << ' ' << child(i)->portable_str();
	}
	ss << ')';
	return ss.str();
}*/

const char *
basic_tree::latex_str() const {
	if (latex.empty()) {
		const string_builder *lb = lbuilder();
		if (lb) {
			latex = lb->build(this);
		}
	}
	return latex.c_str();
}

const char *
ExpressionTree::getstr(int id) const {
	if (id == strid::mathml) {
		return long_str();
	}
	else if (id == strid::latex) {
		return latex_str();
	}
	return 0;
}

void
basic_tree::wrap() const {
	std::stringstream ss;
	ss << "<math xmlns='http://www.w3.org/1998/Math/MathML'>" << long_str() << "</math>";
	mathml = ss.str();
}

interpreter *
basic_tree::mkiter(math_recognizer_base *rec) const {
	return (nt() && span()) ? rec->mkparser(nt(), span())/*rec->getfullmux(span())*/ : 0;
}

invented_tree::invented_tree(reader &re, math_recognizer_base *rec)
	: basic_tree(re, rec), P(0) {
	P = readproduction(re);
}

void
invented_tree::write(writer &wr) const {
	basic_tree::write(wr);
	if (!P) {
		wr.write((size_t)0);
		wr.write((size_t)0);
	}
	else {
		wr.write(P->nt ? P->nt->index + 1: 0);
		wr.write(P->index);
	}
}

basic_tree *
invented_tree::mkcopy_() const {
	return new invented_tree(P);
}

const nonterminal *
invented_tree::nt() const {
	return P ? P->nt : 0;
}

SemanticId
invented_tree::type() const {
	return P ? P->sid : InvalidSemanticId;
}

const string_builder *
invented_tree::sbuilder() const {
	return P ? P->sbuild : 0;
}

const string_builder *
invented_tree::lbuilder() const {
	return P ? P->lbuild : 0;
}


parsed_tree::parsed_tree(reader &re, math_recognizer_base *rec)
	: basic_tree(re, rec), intrp_(0) {
	re.read(const_cast<interpretation **>(&intrp_), rec);
}

void
parsed_tree::write(writer &wr) const {
	basic_tree::write(wr);
	wr.write(intrp_);
}

basic_tree *
parsed_tree::mkcopy_() const {
	return new parsed_tree(intrp_);
}

const string_builder *
parsed_tree::sbuilder() const {
	const interpretation *sintrp = getsemanticintrp(intrp_);
	return (sintrp && sintrp->P) ? sintrp->P->sbuild : ((intrp_->P && intrp_->P->sbuild) ? intrp_->P->sbuild : 0);
}

const string_builder *
parsed_tree::lbuilder() const {
	const interpretation *sintrp = getsemanticintrp(intrp_);
	return (sintrp && sintrp->P) ? sintrp->P->lbuild : ((intrp_->P && intrp_->P->lbuild) ? intrp_->P->lbuild : 0);
}

double
parsed_tree::score() const {
	return intrp_ ? intrp_->score().score : 0;
}

const ExpressionBox *
parsed_tree::box() const {
	return new expression_box(intrp_->bounds());
}

SemanticId
parsed_tree::type() const {
	const interpretation *pnode = getsemanticintrp(intrp_);
	return pnode->P ? pnode->P->sid : InvalidSemanticId;
}

const std::string &
parsed_tree::ccstr() const {
	return basic_tree::ccstr();
	if (raw.empty()) {
		const interpretation *pnode = getsemanticintrp(intrp_);
		if (type() == TERMINAL_EXPR && !pnode->P) {
			raw = pnode->str();
		}
		else {
			std::stringstream ss;
			ss << '(' << (pnode->P && pnode->P->nt ? pnode->P->nt->name : "nil") << " ";
			for (size_t i = 0; i < nchildren(); ++i) {
				const basic_tree *ch = child(i);
				if (i != 0) {
					ss << ' ';
				}
				ss << ch->ccstr();
			}
			ss << ')';
			raw = ss.str();
		}
	}	    	    
	return raw;
}

/*
const std::wstring &
parsed_tree::ccwstr() const {
	if (wraw.empty()) {
		const interpretation *pnode = getsemanticintrp(intrp_);
		if (type() == TERMINAL_EXPR) {
			wraw = pnode->wstr();
		}
		else {
			std::wstringstream ss;
			if (pnode->P && pnode->P->nt) {
				size_t mbsz = mbstowcs(0, pnode->P->nt->name.c_str(), 0);
				wchar_t *buf = new wchar_t[mbsz+1];
				mbstowcs(buf, pnode->P->nt->name.c_str(), mbsz+1);
				delete[] buf;
			}
			ss << L"(" << ss.str() << L" ";
			for (size_t i = 0; i < nchildren(); ++i) {
				const basic_tree *ch = child(i);
				if (i != 0) {
					ss << L" ";
				}
				ss << ch->ccwstr();
			}
			ss << L")";
			wraw = ss.str();
		}
	}
	return wraw;
}*/

int
parsed_tree::lock() const {
	if (!intrp_ || !intrp_->src) {
		ERR(E_INVALID, "This expression did not emerge from a recognition context.");
		return get_error().code;
	}
	const interpretation *intrp = getlockintrp(intrp_);
	//assert(intrp->P);
	//return rec->lock(intrp->src);
	assert(intrp->src);
	//std::cout << "LOCKING interpreter " << intrp->src << std::endl;
	if (intrp->grp) {
		assert(intrp->nchildren() == 0);
		group *grp = intrp->grp;
		std::map<const nonterminal *, match_score>::iterator i = grp->matches.find(intrp->nt());
		if (i != grp->matches.end()) {
			match_score &match = i->second;
			if (match.final_score.hint == 0) {
				match.final_score.hint = 1;
			}
		}
		/*else {
			std::cout << "COULD NOT FIND TERMINAL SCORE\n";
		}*/
	}
	return intrp->src->lock();
}

int
parsed_tree::unlock() const {
	if (!intrp_ || !intrp_->src) {
		ERR(E_INVALID, "This expression did not emerge from a recognition context.");
		return get_error().code;
	}
	//return rec->unlock(intrp_->src);
	const interpretation *intrp = getlockintrp(intrp_);
	assert(intrp->src);
	//std::cout << "UNLOCKING interpreter " << intrp->src << std::endl;
	if (intrp->grp) {
		assert(intrp->nchildren() == 0);
		group *grp = intrp->grp;
		std::map<const nonterminal *, match_score>::iterator i = grp->matches.find(intrp->nt());
		if (i != grp->matches.end()) {
			match_score &match = i->second;
			if (match.final_score.hint == 1) {
				match.final_score.hint = 0;
			}
		}
		/*else {
			std::cout << "COULD NOT FIND TERMINAL SCORE\n";
		}*/
	}
	return intrp->src->unlock();
}

bool
parsed_tree::is_locked() const {
	/*if (!rec || !intrp_ || !intrp_->src) {
		//THROW_ERROR(E_INVALID, "This expression did not emerge from a recognition context.");
		return false;
	}
	return rec->islocked(intrp_->src);*/
	const interpretation *intrp = getlockintrp(intrp_);
	assert(intrp->src);
	//std::cout << "ISLOCKED? interpreter " << intrp->src << std::endl;
	if (intrp->grp) {
		assert(intrp->nchildren() == 0);
		group *grp = intrp->grp;
		std::map<const nonterminal *, match_score>::iterator i = grp->matches.find(intrp->nt());
		if (i != grp->matches.end()) {
			match_score &match = i->second;
			return (match.final_score.hint == 1);
		}
		/*else {
			std::cout << "COULD NOT FIND TERMINAL SCORE\n";
		}*/
	}
	return intrp->src->islocked();
}

bool
parsed_tree::HasLongForm() const {
	return getsemanticintrp(intrp_)->haslongform();
}

ExpressionIterator *
parsed_tree::CreateLongFormIterator() const {
	const interpretation *sintrp = getsemanticintrp(intrp_);
	return sintrp->haslongform() ? new external_iterator(sintrp->longformiter(), true, true) : 0;
}


terminal_tree::terminal_tree(const interpretation *intrp,
							 const std::string &str, const std::wstring &wstr, const std::string &mathml_, const std::string &latex_)
	: parsed_tree(intrp) {
	raw = str;
	wraw = wstr;
	if (!intrp || !intrp->P || !intrp->P->sbuild) {
		mathml = mathml_;
	}
	if (!intrp || !intrp->P || !intrp->P->lbuild) {
		latex = latex_;
	}
}

terminal_tree::terminal_tree(reader &re, math_recognizer_base *rec)
	: parsed_tree(re, rec) {
	re.read(raw);
	re.read(wraw);
	re.read(mathml);
	re.read(latex);
}

void
terminal_tree::write(writer &wr) const {
	parsed_tree::write(wr);
	wr.write(raw);
	wr.write(wraw);
	wr.write(mathml);
	wr.write(latex);
}

basic_tree *
terminal_tree::mkcopy_() const {
	return new terminal_tree(intrp_, raw, wraw, mathml, latex);
}

struct blank {
	blank() : P(-2) {
		P.sid = BLANK_EXPR;
		P.tbuild = new blank_tree_builder;
		P.sbuild = new basic_string_builder("");
		P.lbuild = new basic_string_builder("");
	}

	uniqproduction P;
};
static blank blankspec;

basic_tree *
mkblank() {
	return new invented_tree(&blankspec.P);
}

struct placeholder {
	placeholder() : P(-3) {
		P.sid = PLACEHOLDER_EXPR;
		P.tbuild = new placeholder_tree_builder;
		P.sbuild = new basic_string_builder("<mi>?</mi>");
		P.lbuild = new basic_string_builder("?");
	}

	uniqproduction P;
};
static placeholder phspec;

basic_tree *
mkplaceholder() {
	return new invented_tree(&phspec.P);
}

int
initcanonicalsids(grammar *G) {
	G->setcanonicalsidproduction(PLACEHOLDER_EXPR, &phspec.P);
	G->setcanonicalsidproduction(BLANK_EXPR, &blankspec.P);
	return 0;
}

}
