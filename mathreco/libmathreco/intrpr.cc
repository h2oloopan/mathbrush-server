#include "intrpr.h"
#include "expr-node.h"
#include "grammar.h"
#include "mathrecognizer-private.h"
#include "verb.h"
#include "error.h"
#include <cassert>
#include <algorithm>
#include <iostream>

namespace scg {

const char staticintrpr::ID = 1;
//const char iterator::ID = 2;
const char multiplexor::ID = 3;
const char sequence::ID = 4;

interpreter::interpreter(math_recognizer_base *rec_, const nonterminal *nt, const ordered_segments *span, int priority)
	: rec(rec_), nt_(nt), span_(span), hint(0), resetpriority(priority), nlinks(0) {
	//std::cout << "new " << this << std::endl;
}

interpreter::~interpreter() {
	//std::cout << "del " << this << " for " << nt()->name << ":" << span()->bits << std::endl;
	for (size_t i = 0; i < children.size(); ++i) {
		interpreter *child = children[i];
		if (rmflags[i]) {
			delete child;
		}
		else if (child) { // child might have already been deleted
			//std::cout << "cleaning up child " << child << std::endl;
			std::set<interpreter *>::iterator j = child->parents.find(this);
			assert(j != child->parents.end());
			child->parents.erase(j);
		}
	}
	for (std::set<interpreter *>::iterator i = parents.begin(); i != parents.end(); ++i) {
		interpreter *parent = *i;
		parent->rmchild(this);
	}
}

void
interpreter::addchild(interpreter *child, bool rm) {
	//std::cout << "addchild " << child << " (" << child->nt()->name << ":" << child->span()->bits << ") to " << this << " (" << nt()->name << ":" << span()->bits << ")\n";
	assert(child);
	assert(child != this);
	//assert(child->span()->bits.subset_of(span()->bits));
	//assert(std::find(children.begin(), children.end(), child) == children.end());
	children.push_back(child);
	rmflags.push_back(rm);
	child->parents.insert(this);
}

void
interpreter::rmchild(interpreter *child) {
	assert(child);
	assert(child != this);
	for (std::vector<interpreter *>::iterator i = children.begin(); i != children.end(); ++i) {
		if (*i == child) {
			*i = 0; // don't delete; this is called from child's destructor
		}
	}
}

void
interpreter::clearchildren() {
	for (std::vector<interpreter *>::iterator i = children.begin(); i != children.end(); ++i) {
		interpreter *child = *i;
		std::set<interpreter *>::iterator j = child->parents.find(this);
		assert(j != child->parents.end());
		child->parents.erase(j);		
	}
	children.clear();
}

struct intrprcmp_for_reset {
	inline bool checknts(const interpreter *lhs, const interpreter  *rhs) const {
		const nonterminal *lnt = lhs->nt();
		const nonterminal *rnt = rhs->nt();
		return rnt->derives(lnt);
	}

	bool operator()(const interpreter *lhs, const interpreter *rhs) const {
		const bitvec &lbits = lhs->span()->bits;
		const bitvec &rbits = rhs->span()->bits;
		const nonterminal *lnt = lhs->nt();
		const nonterminal *rnt = rhs->nt();
		VERBOSE(*verb_out << "cmp " << lhs << " " << lnt->name << "/" << lhs->resetpriority << ":" << lbits << "\n    " << rhs << " " << rnt->name << "/" << rhs->resetpriority << ":" << rbits << " === ");
		bool r;
		if (lbits == rbits) {
			if (lnt == rnt) {
				if (lhs->resetpriority == rhs->resetpriority) {
					r = lhs < rhs;
				}
				else {
					r = lhs->resetpriority < rhs->resetpriority;
				}
			}
			else if (rnt->has_direct_descendent(lnt)) r = true;
			else if (lnt->has_direct_descendent(rnt)) r = false;
			else if (lnt->rootdist == rnt->rootdist) r = lnt->index < rnt->index;
			else r = lnt->rootdist > rnt->rootdist;
		}
		else if (lbits.subset_of(rbits)) {
			r = true;
		}
		else if (rbits.subset_of(lbits)) {
			r = false;
		}
		else r = lbits < rbits;
		VERBOSE(*verb_out << (r ? "yes" : "no") << std::endl);
		return r;
	}
};

void
interpreter::countlinks() {
	++nlinks;
	if (nlinks == 1) {
		for (std::set<interpreter *>::iterator i = parents.begin(); i != parents.end(); ++i) {
			(*i)->countlinks();
		}
	}
}

void
interpreter::resetparents() {
	std::set<interpreter *> Q;
	std::set<interpreter *> done;
	assert(nlinks == 0);
	countlinks();
	nlinks = 0;
	Q.insert(this);
	while (!Q.empty()) {
		interpreter *intrpr = *Q.begin();
		Q.erase(Q.begin());
		assert(done.find(intrpr) == done.end() && (done.insert(intrpr), true)); 
		intrpr->reset();
		for (std::set<interpreter *>::iterator i = intrpr->parents.begin(); i != intrpr->parents.end(); ++i) {
			interpreter *p = *i;
			if (--p->nlinks == 0) {
				Q.insert(p);
			}
		}
	}
}

bool
interpreter::islocked() const {
	return hint == 1;
}

void
interpreter::addreset() {
	rec->resets.push(this);
	/*interpreter *resetintrpr = rec->getparser(nt(), span());
	size_t childi = ~0;
	if (resetintrpr != this) {
		for (childi = 0; childi < resetintrpr->nchildren(); ++childi) {
			if (resetintrpr->child(childi) == this) {
				break;
			}
		}
		if (childi == resetintrpr->nchildren()) std::abort();
	}
	rec->resets.push(math_recognizer_base::resetspec(nt(), span()->bits, childi));*/
}

int
interpreter::lock() {
	if (hint == 0) {
		hint = 1;
		addreset();
	}
	return 0;
}

int
interpreter::unlock() {
	if (hint == 1) {
		hint = 0;
		addreset();
	}
	return 0;
}

interpreter::interpreter()
	: rec(0), nt_(0), span_(0), resetpriority(0), nlinks(0) {
}

void
interpreter::read(reader &re, math_recognizer_base *rec_) {
	rec = rec_;
	size_t ntindex;
	re.read(ntindex);
	if (ntindex > 0) {
		--ntindex;
		const grammar &G = *GetMathGrammar();
		if (ntindex >= G.nts.size()) {
			THROW_ERROR(E_INVALID, "while reading interpreter, nt index " << ntindex << " is too large");
		}
		nt_ = G.nts[ntindex];
	}
	else {
		nt_ = 0;
	}
	bitvec bits;
	re.read(bits);
	if (bits.size() == 0) {
		span_ = 0;
	}
	else {
		span_ = rec->getsegs(bits);
		if (!span_) {
			THROW_ERROR(E_INVALID, "while reading interpreter, found unknown span " <<  bits);
		}
	}
	//re.read((long &)atrev);
	re.read(hint);
	long prox;
	re.read(prox);
	resetpriority = prox;
	size_t nchil;
	re.read(nchil);
	for (size_t i = 0; i < nchil; ++i) {
		interpreter *child;
		re.read(&child, rec);
		bool rm;
		re.read(rm);
		addchild(child, rm);
	}
	//interpreter *prox;
	//re.read(&prox, rec);
	//blocking_lock = prox;
}

void
interpreter::write(writer &wr) const {
	//((interpreter *)this)->update();
	/*
	if (atrev != rec->rev()) {
		((interpreter *)this)->reset();
	}*/
	wr.write(nt() ? nt()->index+1 : 0);
	wr.write(span() ? span()->bits : bitvec::EMPTY);
	wr.write(hint);
	wr.write((long)resetpriority);
	//wr.write((long)atrev);
	wr.write(nchildren());
	for (size_t i = 0; i < nchildren(); ++i) {
		wr.write(children[i]);
		wr.write(rmflags[i]);
	}
	//wr.write(blocking_lock);
}


staticintrpr::staticintrpr(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span)
	: interpreter(rec, nt, span, 0), at(0) {
}

staticintrpr::~staticintrpr() {
	for (size_t i = 0; i < ts.size(); ++i) {
		if (rmflags[i]) {
			delete ts[i];
		}
	}
}

size_t
staticintrpr::nknown() const {
	return at;
}

interpretation *
staticintrpr::nth(size_t i) const {
	return (i < ts.size()) ? ts[i] : 0;
}

size_t
staticintrpr::nadded() const {
	return ts.size();
}

interpretation *
staticintrpr::next() {
	if (at == ts.size()) {
		return 0;
	}
	else {
		ts[at]->hint += hint;
		return ts[at++];
	}
}

void
staticintrpr::reset() {
	at = 0;
}

staticintrpr::staticintrpr() : interpreter(), at(0) { resetpriority = 0; }

void
staticintrpr::read(reader &re, math_recognizer_base *rec) {
	interpreter::read(re, rec);
	size_t nts;
	re.read(nts);
	ts.resize(nts);
	rmflags.resize(nts);
	for (size_t i = 0; i < ts.size(); ++i) {
		re.read(&ts[i], rec);
		bool prox;
		re.read(prox);
		rmflags[i] = prox;
	}
	//re.read(at);
}

void
staticintrpr::write(writer &wr) const {
	interpreter::write(wr);
	wr.write(ts.size());
	for (size_t i = 0; i < ts.size(); ++i) {
		wr.write(ts[i]);
		wr.write(rmflags[i]);
	}
	//wr.write(at);
}

bool
cmpscored(const interpretation *lhs, const interpretation *rhs) {
	return lhs->score() > rhs->score();
}

void
staticintrpr::addknown(interpretation *t, bool rm) {
	std::vector<interpretation *>::iterator i = std::upper_bound(ts.begin(), ts.end(), t, &cmpscored);
	size_t n = i - ts.begin();
	ts.insert(i, t);
	rmflags.insert(rmflags.begin() + n, rm);
}

parser::parser(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span, int priority)
	: interpreter(rec, nt, span, priority) { }

void
parser::freeknown() {
	for (std::vector<interpretation *>::iterator i = known.begin(); i != known.end(); ++i) {
		delete *i;
	}
	known.clear();
}

void
parser::clearknown() {
	known.clear();
}

size_t
parser::nknown() const {
	return known.size();
}

interpretation *
parser::nth(size_t i) const {
	if (i >= known.size()) {
		THROW_ERROR(E_INVALID, "interpretation index " << i << " is too large for parser " << this);
	}
	return known[i];
}

interpretation *
parser::next() {
	interpretation *intrp = getnext();
	if (!intrp) {
		return 0;
	}
	//VERBOSE(*verb_out << "parser::next() adding hint " << hint << " to interpretation " << intrp << std::endl);
	intrp->hint += hint;
	intrp->clrstrs();
	//VERBOSE(*verb_out << "parser::next() result: " << intrp->str() << std::endl);
	known.push_back(intrp);
	return intrp;
}

parser::parser() : interpreter() { }

void
parser::read(reader &re, math_recognizer_base *rec) {
	interpreter::read(re, rec);
	/*
	size_t nknown;
	re.read(nknown);
	known.resize(nknown);
	for (size_t i = 0; i < nknown; ++i) {
		re.read(&known[i], rec);
		VERBOSE(*verb_out << "read(): known to " << this << " (" << nt()->name << ":" << span()->bits << ") is " << known[i]->ccstr() << std::endl);
	}*/
}

void
parser::write(writer &wr) const {
	interpreter::write(wr);
	/*
	wr.write(known.size());
	for (std::vector<interpretation *>::const_iterator i = known.begin(); i != known.end(); ++i) {
		VERBOSE(*verb_out << "write(): known to " << this << " (" << nt()->name << ":" << span()->bits << ") is " << (*i)->ccstr() << std::endl);
		wr.write(*i);
	}*/
}

/*
iterator::iterator(interpreter *src_, bool ownptr)
	: interpreter(src_->ctx(), src_->nt(), src_->span(), src_->rev()), curr(0) {
	addchild(src_, ownptr);
}

size_t
iterator::nknown() const {
	return curr;
}

interpretation *
iterator::nth(size_t i) const {
	return child(0)->nth(i);
}

interpretation *
iterator::next() {
	if (curr < child(0)->nknown()) {
		return child(0)->nth(curr++);
	}
	interpretation *nx = child(0)->next();
	if (nx) {
		++curr;
	}
	return nx;
}

void
iterator::reset() {
	curr = 0;
}

iterator::iterator() : interpreter(), curr(0) { }

void
iterator::read(reader &re, math_recognizer_base *rec) {
	interpreter::read(re, rec);
	//re.read(curr);
}

void
iterator::write(writer &wr) const {
	interpreter::write(wr);
	//wr.write(curr);
}*/

multiplexor::multiplexor(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span, int priority)
	: parser(rec, nt, span, priority), prevsrc(0) {
}

void
multiplexor::addparser(interpreter *intrpr, bool rm, int prior) {
	addchild(intrpr, rm);
	intrprpos.push_back(0);
	priority.push_back(prior);
	interpretation *intrp = getfirst(intrpr);
	if (intrp) {
		addintrp(intrp, nchildren() - 1, prior);
	}
}

interpretation *
multiplexor::getnext() {
	VERBOSE(*verb_out << "mux:" << this << "->getnext() with " << known.size() << " known\n");
	if (prevsrc) {
		size_t idx = prevsrc - 1;
		interpreter *intrpr = child(idx);
		++intrprpos[idx];
		interpretation *intrp = getnth(intrpr, intrprpos[idx]);
		if (intrp) {
			addintrp(intrp, idx, priority[idx]);
		}
	}
	intrp_t best;
	if (!known.empty()) {
		best = known.back();
		known.pop_back();
		prevsrc = best.src;
	}

	return best.intrp;
}

multiplexor::intrp_t::intrp_t(interpretation *intrp_, size_t src_, int priority_) : intrp(intrp_), src(src_), priority(priority_) { }

bool
multiplexor::intrp_t::operator<(const multiplexor::intrp_t &rhs) const {
	recoscore sc = intrp->score();
	recoscore rhssc = rhs.intrp->score();
	if (sc < rhssc) {
		return true;
	}
	else if (sc == rhssc) {
		return priority < rhs.priority;
	}
	return false;
}

void
multiplexor::addintrp(interpretation *intrp, size_t srcidx, int priority) {
	intrp_t val(intrp, srcidx + 1, priority);
	std::vector<intrp_t>::iterator i = std::upper_bound(known.begin(), known.end(), val);
	known.insert(i, val);
}

void
multiplexor::reset() {
	VERBOSE(*verb_out << this << "->reset() begins with " << known.size() << "known\n");
	parser::clearknown();
	known.clear();
	prevsrc = 0;
	for (size_t i = 0; i < nchildren(); ++i) {
		interpreter *intrpr = child(i);
		intrprpos[i] = 0;
		interpretation *intrp = getfirst(intrpr);
		if (intrp) {
			addintrp(intrp, i, priority[i]);
		}
	}
	VERBOSE(*verb_out << this << "->reset() gives " << known.size() << " known\n");
}

multiplexor::multiplexor()
	: parser(), prevsrc(0), intrprpos(nchildren(), 0) { }

void
multiplexor::read(reader &re, math_recognizer_base *rec) {
	parser::read(re, rec);
	priority.resize(nchildren());
	for (size_t i = 0; i < nchildren(); ++i) {
		long prox;
		re.read(prox);
		priority[i] = (int)prox;
	}
	intrprpos.insert(intrprpos.end(), nchildren(), 0);
	reset();
	/*
	intrprpos.resize(nchildren());
	priority.resize(nchildren());
	for (size_t i = 0; i < nchildren(); ++i) {
		re.read(intrprpos[i]);
		long prox;
		re.read(prox);
		priority[i] = (int)prox;
	}
	//reset();
	re.read(prevsrc);
	size_t nknown;
	re.read(nknown);
	known.resize(nknown);
	for (size_t i = 0; i < nknown; ++i) {
		re.read(&known[i].intrp, rec);
		re.read(known[i].src);
		known[i].priority = priority[known[i].src-1];
	}*/
}

void
multiplexor::write(writer &wr) const {
	parser::write(wr);
	for (size_t i = 0; i < nchildren(); ++i) {
		wr.write((long)priority[i]);
	}
	/*
	for (size_t i = 0; i < nchildren(); ++i) {
		wr.write(intrprpos[i]);
		wr.write((long)priority[i]);
	}
	wr.write(prevsrc);
	wr.write(known.size());
	for (std::vector<intrp_t>::const_iterator i = known.begin(); i != known.end(); ++i) {
		wr.write(i->intrp);
		wr.write(i->src);
	}*/
}


sequence::sequence(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span, int priority)
	: parser(rec, nt, span, priority), curintrpr(~0), curpos(0) { }

void
sequence::reset() {
	parser::clearknown();
	curintrpr = (nchildren() == 0) ? ~0 : 0;
	curpos = 0;
}

void
sequence::addparser(interpreter *intrpr, bool rm) {
	addchild(intrpr, rm);
	if (curintrpr == ~0) {
		curintrpr = 0;
	}
}

void
sequence::read(reader &re, math_recognizer_base *rec) {
	parser::read(re, rec);
	reset();
	/*long dummy;
	re.read(dummy);
	curintrpr = (size_t)dummy;
	re.read(dummy);
	curpos = (size_t)dummy;*/
}

void
sequence::write(writer &wr) const {
	parser::write(wr);
	//wr.write((long)curintrpr);
	//wr.write((long)curpos);
}

interpretation *
sequence::getnext() {
	if (curintrpr == ~0) return 0;
	interpretation *intrpr = getnth(child(curintrpr), curpos++);
	if (!intrpr) {
		curpos = 0;
		++curintrpr;
		if (curintrpr == nchildren()) {
			curintrpr = ~0;
			return 0;
		}
		return getnext();
	}
	return intrpr;
}


interpretation *
getnth(interpreter *intrpr, size_t i) {
	if (i >= intrpr->nknown()) {
		interpretation *intrp;
		do {
			intrp = intrpr->next();
			if (!intrp) {
				return 0;
			}
		} while (i >= intrpr->nknown());
		return intrp;
	}
	return intrpr->nth(i);
}

interpretation *
getfirst(interpreter *intrpr) {
	return (intrpr->nknown() > 0) ? intrpr->nth(0) : intrpr->next();
}

}
