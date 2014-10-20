#include "ptab.h"
#include "mathrecognizer-private.h"
#include "symbols.h"
#include "grammar-values.h"
#include "expr-node.h"
#include <vector>
#include <cassert>
#include <ostream>
#include <cctype>

namespace scg {

typedef std::vector<std::set<SemanticId> > sidslots_t; // map tree child index -> applicable sids
std::map<SemanticId, sidslots_t> sidslots; // map sid to its child slots
std::map<SemanticId, sidslots_t> termsidslots; // map sid to its TERMINAL child slots
std::map<const nonterminal *, std::set<SemanticId> > ntsids; // map nt to all its sids



template <typename T>
static void
ensureindex(std::vector<T> &v, size_t i) {
	assert(i <= v.size());
	if (i == v.size()) v.push_back(T());
}

struct hook : public grammar_hook {
	SemanticId cursid;
	sidslots_t *curslots;
	size_t curchild;

	struct refdeferral {
		SemanticId sid;
		size_t childi;
		const nonterminal *srcnt;
		refdeferral(SemanticId sid_, size_t childi_, const nonterminal *srcnt_) : sid(sid_), childi(childi_), srcnt(srcnt_) { }
		void resolve() {
			sidslots_t &slots = sidslots[sid];
			assert(!slots.empty());
			assert(childi < slots.size());

			//std::cout << "resolving ref deferral for slot " << childi << " in " << sid_to_string(sid) << " from nonterminal " << srcnt->name << ": adding \n";
			if (srcnt->isterminal()) {
				slots[childi].insert(string_to_sid(srcnt->name.substr(2)));
			}
			std::set<SemanticId> &srcsids = ntsids[srcnt];
			/*for (std::set<SemanticId>::const_iterator i = srcsids.begin(); i != srcsids.end(); ++i) {
				std::cout << "  " << sid_to_string(*i) << std::endl;
			}*/
			slots[childi].insert(srcsids.begin(), srcsids.end());
		}
	};

	std::vector<refdeferral> refdefs;


	struct subrefdeferral {
		SemanticId sid;
		size_t childi;
		const nonterminal *srcnt;
		size_t srcchildi;
		subrefdeferral(SemanticId sid_, size_t childi_, const nonterminal *srcnt_, size_t srcchildi_) : sid(sid_), childi(childi_), srcnt(srcnt_), srcchildi(srcchildi_) { }
		void resolve() {
			sidslots_t &slots = sidslots[sid];
			assert(!slots.empty());
			assert(childi < slots.size());
			const std::set<SemanticId> &srcsids = ntsids[srcnt];
			for (std::set<SemanticId>::const_iterator i = srcsids.begin(); i != srcsids.end(); ++i) {
				sidslots_t &srcslots = sidslots[*i];
				assert(!srcslots.empty());
				assert(srcchildi < srcslots.size());
				const std::set<SemanticId> &childsids = srcslots[srcchildi];
				slots[childi].insert(childsids.begin(), childsids.end());
			}
		}
	};

	std::vector<subrefdeferral> subrefdefs;

	bool canresolve(const subrefdeferral &def) {
		const std::set<SemanticId> &srcsids = ntsids[def.srcnt];
		for (std::vector<subrefdeferral>::const_iterator i = subrefdefs.begin(); i != subrefdefs.end(); ++i) {
			if (srcsids.find(i->sid) != srcsids.end() && i->childi == def.srcchildi) {
				return false;
			}
		}
		return true;
	}

	hook() : cursid(InvalidSemanticId), curslots(0), curchild(0) { }

	void start_tree(const nonterminal *nt, SemanticId sid) {
		//if (sidismeaningful(sid)) {
		ntsids[nt].insert(sid);
		//}
		assert(!curslots);
		curslots = &sidslots[sid];
		cursid = sid;
	}

	void add_sid_child(SemanticId sid) {
		ensureindex(*curslots, curchild);
		(*curslots)[curchild].insert(sid);
		++curchild;
	}

	void add_ref_child(const nonterminal *nt) {
		ensureindex(*curslots, curchild);
		refdefs.push_back(refdeferral(cursid, curchild, nt));
		++curchild;
	}

	void add_subref_child(const nonterminal *nt, size_t i) {
		ensureindex(*curslots, curchild);
		subrefdefs.push_back(subrefdeferral(cursid, curchild, nt, i));
		++curchild;
	}

	void end_tree() {
		assert(curslots);
		curslots = 0;
		curchild = 0;
		cursid = InvalidSemanticId;
	}

	void finalize() {
		for (std::vector<refdeferral>::iterator i = refdefs.begin(); i != refdefs.end(); ++i) {
			i->resolve();
		}
		refdefs.clear();
		bool progress = true;
		while (!subrefdefs.empty() && progress) {
			progress = false;
			for (std::vector<subrefdeferral>::iterator i = subrefdefs.begin(); i != subrefdefs.end(); ) {
				if (canresolve(*i)) {
					i->resolve();
					i = subrefdefs.erase(i);
					progress = true;
				}
				else {
					++i;
				}
			}
		}
	}
};

static hook ghook;

void
ptab_staticinit() {
	set_grammar_hook(&ghook);
}

void
ptab_rectify(const grammar *G) {
	for (std::vector<nonterminal *>::const_iterator i = G->nts.begin(); i != G->nts.end(); ++i) {
		std::set<SemanticId> &sids = ntsids[*i];
		const std::set<const nonterminal *> &equivs = (*i)->equivalent_nts;
		for (std::set<const nonterminal *>::const_iterator j = equivs.begin(); j != equivs.end(); ++j) {
			const std::set<SemanticId> &equivsids = ntsids[*j];
			sids.insert(equivsids.begin(), equivsids.end());
		}
	}
	ghook.finalize();
	/*
	for (std::map<const nonterminal *, std::set<SemanticId> >::const_iterator i = ntsids.begin(); i != ntsids.end(); ++i) {
		std::cout << i->first->name << " derives:\n";
		const std::set<SemanticId> &sids = i->second;
		for (std::set<SemanticId>::const_iterator j = sids.begin(); j != sids.end(); ++j) {
			std::cout << "  " << sid_to_string(*j) << std::endl;
		}
	}*/

	for (std::map<SemanticId, sidslots_t>::iterator i = sidslots.begin(); i != sidslots.end(); ) {
		SemanticId sid = i->first;
		if (sid == InvalidSemanticId) {
			sidslots_t &slots = i->second;
			sidslots_t &termslots = termsidslots[sid];
			termslots.clear();
			if (slots.size() > 1) {
				slots.erase(slots.begin() + 1, slots.end());
			}
			++i;
		}
		else if (!sidismeaningful(sid)) {
			sidslots.erase(i);
			i = sidslots.upper_bound(sid);
		}
		else {
			sidslots_t &slots = i->second;
			sidslots_t &termslots = termsidslots[sid];
			termslots.resize(slots.size());
			for (size_t j = 0; j < slots.size(); ++j) {
				std::set<SemanticId> &sidset = slots[j];
				for (std::set<SemanticId>::iterator k = sidset.begin(); k != sidset.end(); ) {
					SemanticId slotsid = *k;
					if (sidisterminal(slotsid) || !sidismeaningful(slotsid)) {
						if (sidisterminal(slotsid)) {
							termslots[j].insert(slotsid);
						}
						sidset.erase(k);
						k = sidset.upper_bound(sid);
					}
					else {
						++k;
					}
				}
			}
			++i;
		}
	}

	VERBOSE(
		*verb_out << "PTAB slot mapping::\n";
		for (std::map<SemanticId, sidslots_t>::const_iterator i = sidslots.begin(); i != sidslots.end(); ++i) {
			*verb_out << sid_to_string(i->first) << " has slots:\n";
			const sidslots_t &slots = i->second;
			for (size_t j = 0; j < slots.size(); ++j) {
				*verb_out << "  " << j << ":\n";
				for (std::set<SemanticId>::const_iterator k = slots[j].begin(); k != slots[j].end(); ++k) {
					*verb_out << "    " << sid_to_string(*k) << std::endl;
				}
			}
		}
		for (std::map<SemanticId, sidslots_t>::const_iterator i = termsidslots.begin(); i != termsidslots.end(); ++i) {
			*verb_out << sid_to_string(i->first) << " has terminal slots:\n";
			const sidslots_t &slots = i->second;
			for (size_t j = 0; j < slots.size(); ++j) {
				*verb_out << "  " << j << ":\n";
				for (std::set<SemanticId>::const_iterator k = slots[j].begin(); k != slots[j].end(); ++k) {
					*verb_out << "    " << sid_to_string(*k) << std::endl;
				}
			}
		}
	);
}


static void
postprocess_tree(typetree &tt) {
	/*std::vector<typetree>::const_iterator i;
	for (i = tt.children.begin(); i != tt.children.end(); ++i) {
		if (!sidisterminal(i->sid)) {
			break;
		}
	}
	if (i == tt.children.end()) {
		typetree termchild(TERMINAL_EXPR);
		termchild.children = tt.children;
		termchild.children.push_back(typetree());
		tt.children.clear();
		tt.children.push_back(termchild);
	}*/

	if (tt.sid == REL_EXPR || tt.sid == ADD_EXPR || tt.sid == MULT_EXPR) {
		typetree &chtt = tt.children.back();
		if (chtt.sid == TERMINAL_EXPR) {
			chtt.sid = chtt.children.front().sid;
			chtt.children.clear();
		}
	}
}

int
mktypetree(typetree &tt, const ExpressionTree *tree) {
	/*
	intrp = getsemanticintrp(intrp);
	if (intrp->P && intrp->P->sid != InvalidSemanticId && intrp->P->sid != TERMINAL_EXPR) {
		tt.sid = intrp->P->sid;
	}

	if (intrp->nchildren() == 1 && tt.sid == InvalidSemanticId) {
		return mktypetree(tt, intrp->child(0));
	}
	for (size_t i = 0; i < intrp->nchildren(); ++i) {
		const interpretation *chintrp = intrp->child(i);
		typetree chtree;
		int e = mktypetree(chtree, chintrp);
		if (FAILURE(e)) {
			return e;
		}
	}

	postprocess_tree(tt);
	return 0;*/
	
	tt.sid = tree->type();
	if (tt.sid == TERMINAL_EXPR) {
		const wchar_t *ws = tree->wstr();
		while (*ws) {
			const symbol *S = symdb_findsymbol_unicode(*ws);
			if (!S) {
				return E_NOTFOUND;
			}
			tt.children.push_back(typetree());
			tt.children.back().sid = S->sid;
			ws++;
		}
	}
	else {
		for (size_t i = 0; i < tree->nchildren(); ++i) {
			tt.children.push_back(typetree());
			int e = mktypetree(tt.children.back(), tree->child(i));
			if (FAILURE(e)) {
				return e;
			}
		}
		/*
		if (tt.sid == REL_EXPR || tt.sid == ADD_EXPR || tt.sid == MULT_EXPR) {
			const ExpressionTree *ch = tree->child(2);
			assert(ch);
			if (tt.sid == MULT_EXPR && ch->type() == BLANK_EXPR) {
				tt.sid = BLANK_EXPR;
			}
			else {
				assert(ch->type() == TERMINAL_EXPR);
				const symbol *S = symdb_findsymbol_name(ch->str());
				assert(S);
				tt.sid = S->sid;
			}
			tt.children.pop_back();
		}
		else if (tt.sid == PAREN_EXPR) {
			const ExpressionTree *ch = tree->child(0);
			assert(ch && ch->type() == TERMINAL_EXPR);
			const symbol *S = symdb_findsymbol_name(ch->str());
			assert(S);
			tt.sid = S->sid;
			tt.children.pop_front();
			tt.children.pop_back();
		}*/
	}
	postprocess_tree(tt);
	return 0;
}

int
mktypetree(typetree &tt, std::istream &is) {
	while (std::isspace(is.peek())) { is.get(); }
	if (!is) return 0;

	int c = is.peek();
	if (!is) return E_EOF;

	if (c == '(') {
		is.get();
		while (std::isspace(is.peek())) { is.get(); }
		if (!is) return E_EOF;

		std::string type;
		while (is.peek() != ')' && !std::isspace(is.peek())) {
			type += is.get();
		}

		tt.sid = string_to_sid(type);
		if (tt.sid == InvalidSemanticId) {
			return E_INVALID;
		}

		while (c != ')') {
			c = is.peek();

			if (!is) return E_EOF;
			if (c == ')' || std::isspace(c)) {
				is.get();
			}
			else {
				typetree child;
				int e = mktypetree(child, is);
				if (FAILURE(e)) {
					return e;
				}
				tt.children.push_back(child);
			}
		}

		postprocess_tree(tt);
	}
	else {
		tt.sid = TERMINAL_EXPR;
		std::string symval;
		while (!std::isspace(is.peek()) && is.peek() != ')') {
			if (is.peek() == '\\') {
				is.get();
			}
			symval += is.get();
			if (!is) return E_EOF;
		}
		SemanticId sid = string_to_sid(symval);
		/*
		if (sid == InvalidSemanticId) {
			return E_INVALID;
		}
		tt.sid = sid;*/
		if (sid != InvalidSemanticId) {
			tt.children.push_back(typetree(sid));
		}
		else {
			for (size_t i = 0; i < symval.length(); ++i) {
				SemanticId sid = string_to_sid(std::string(1,symval[i]));
				if (sid != InvalidSemanticId) {
					return E_INVALID;
				}
				tt.children.push_back(typetree(sid));
			}
		}
		tt.children.push_back(typetree());
		/*
		char buf[MB_CUR_MAX+1];
		size_t n = 0;
		while (n < symval.length()) {
			wchar_t wc;
			int c = mbtowc(&wc, &symval[n], symval.length()-n);
			if (c == -1) {
				return E_INVALID;
			}
			n += c;
			c = wctomb(buf, wc);
			if (c == -1) {
				return E_INVALID;
			}
			buf[c] = '\0';
			typetree child;
			child.sid = string_to_sid(buf);
			if (child.sid == InvalidSemanticId) {
				return E_INVALID;
			}
			tt.children.push_back(child);
		}*/
	}
	return 0;
}

std::ostream &
operator<<(std::ostream &os, const typetree &tt) {
	if (tt.children.empty()) {
		return os << sid_to_string(tt.sid);
	}
	else {
		os << '(' << sid_to_string(tt.sid);
		for (std::vector<typetree>::const_iterator i = tt.children.begin(); i != tt.children.end(); ++i) {
			os << ' ' << *i;
		}
		os << ')';
	}
	return os;
}



struct typedist {
	std::map<SemanticId, unsigned> D;
	std::map<unsigned, unsigned> freqmap;

	unsigned noptions;
	unsigned nsamples;
	typedist(unsigned noptions_) : noptions(noptions_), nsamples(0) {
		freqmap[0] = noptions;
	}

	void inc(SemanticId sid) {
		unsigned &n = D[sid];
		if (--freqmap[n] == 0 && n != 0) {
			freqmap.erase(n);
		}
		++n;
		++nsamples;
		++freqmap[n];
	}

	void remove(SemanticId sid) {
		std::map<SemanticId, unsigned>::iterator j = D.find(sid);
		if (j != D.end()) {
			--noptions;
			if (--freqmap[j->second] == 0 && j->second != 0) {
				freqmap.erase(j->second);
			}
			nsamples -= j->second;
			D.erase(j);
		}
	}

	bool distisbad() const {
		std::map<unsigned, unsigned>::const_iterator zero = freqmap.begin();
		return 2*zero->second > noptions;
	}

	double P(SemanticId sid) const {
		if (nsamples == 0) return -1.0;
		std::map<unsigned, unsigned>::const_iterator zero = freqmap.begin();
		std::map<unsigned, unsigned>::const_iterator low = zero;
		assert(zero->first == 0);
		++low;
		assert(low != freqmap.end());

		if (distisbad()) return -1.0;
		unsigned nsmooth = low->first * low->second;
		std::map<SemanticId, unsigned>::const_iterator n = D.find(sid);
		unsigned nsid = (n == D.end()) ? 0 : n->second;

		if (nsid == 0) {
			// smooth # samples with minimal non-zero frequency over unseen classes
			return (double)nsmooth / (nsmooth + nsamples) / zero->second;
		}
		else {
			return (double)n->second / (nsmooth + nsamples);
		}
	}
};

static std::ostream &
operator<<(std::ostream &os, const typedist &d) {
	//os << d.ntotal << ' ' << d.nempty << ' ' << d.nsingle << ' ' << d.D.size() << std::endl;
	os << d.D.size() << std::endl;
	for (std::map<SemanticId, unsigned>::const_iterator i = d.D.begin(); i != d.D.end(); ++i) {
		os << '\t' << sid_to_string(i->first) << ' ' << i->second << std::endl;
	}
	return os;
}

static int
read_typedist(typedist *d, std::istream &is) {
	size_t n;
	//is >> d.ntotal >> d.nempty >> d.nsingle >> n;
	is >> n;
	while (n--) {
		std::string sidname;
		is >> sidname;
		unsigned &freq = d->D[string_to_sid(sidname)];
		is >> freq;
		++d->freqmap[freq];
		d->nsamples += freq;
	}
	if (!is) {
		return E_IO;
	}
	return 0;
}

struct ptab {
	struct ntcond {
		SemanticId parent;
		size_t i;
		ntcond(SemanticId parent_, size_t i_) : parent(parent_), i(i_) { }
		inline bool operator<(const ntcond &rhs) const { return parent < rhs.parent || (parent == rhs.parent && i < rhs.i); }
	};
	struct tcond {
		SemanticId parent;
		size_t childi;
		SemanticId pred;
		tcond(SemanticId parent_, size_t childi_, SemanticId pred_) : parent(parent_), childi(childi_), pred(pred_) { }
		inline bool operator<(const tcond &rhs) const { return parent < rhs.parent || (parent == rhs.parent && (childi < rhs.childi || (childi == rhs.childi && pred < rhs.pred))); }
	};
	std::map<ntcond, typedist *> ntD;
	std::map<tcond, typedist *> tD;

	// selections[sid] / chances[sid] gives a fallback probability
	mutable std::map<SemanticId, unsigned> selections;
	mutable std::map<SemanticId, unsigned> chances;

	~ptab() {
		for (std::map<ntcond, typedist *>::iterator i = ntD.begin(); i != ntD.end(); ++i) {
			delete i->second;
		}
		for (std::map<tcond, typedist *>::iterator i = tD.begin(); i != tD.end(); ++i) {
			delete i->second;
		}
	}

	double ntP(SemanticId parent, size_t i, SemanticId child) const {
		const typedist *D = ntdist(ntcond(parent, i));
		VERBOSE(*verb_out << "P(" << sid_to_string(child) << " | " << sid_to_string(parent) << ':' << i << ") = ");
		double p = -1.0;
		if (D) p = D->P(child);
		if (p == -1.0) {
			VERBOSE(*verb_out << "falling back to context-free distribution...");
			p = (double)selections[child] / chances[child];
		}
		VERBOSE(*verb_out << p << std::endl);
		return p;
	}

	double termP(SemanticId parent, size_t i, SemanticId pred, SemanticId child) const {
		const typedist *D = tdist(tcond(parent, i, pred));
		VERBOSE(*verb_out << "P(" << sid_to_string(child) << " | " << sid_to_string(parent) << ':' << i << ", " << sid_to_string(pred) << ") = ");
		double p = -1.0;
		if (D) p = D->P(child);
		if (p == -1.0) {
			VERBOSE(*verb_out << "falling back to context-free distribution...");
			p = (double)selections[child] / chances[child];
		}
		VERBOSE(*verb_out << p << std::endl);
		return p;
	}

	const typedist *ntdist(ntcond c) const {
		std::map<ntcond, typedist *>::const_iterator i = ntD.find(c);
		return i == ntD.end() ? 0 : i->second;
	}
	const typedist *tdist(tcond c) const {
		std::map<tcond, typedist *>::const_iterator i = tD.find(c);
		return i == tD.end() ? 0 : i->second;
	}

	typedist *ntdist(ntcond c) {
		typedist *&dist = ntD[c];
		if (!dist) {
			const std::set<SemanticId> &sids = sidslots[c.parent][c.i];
			dist = new typedist(sids.size());
		}
		return dist;
	}

	typedist *tdist(tcond c) {
		typedist *&dist = tD[c];
		if (!dist) {
			const std::set<SemanticId> &sids = termsidslots[c.parent][c.childi];
			//VERBOSE(*verb_out << "tdist " << sid_to_string(c.parent) << ':' << c.childi << " has " << sids.size() << " sid options\n");
			dist = new typedist(sids.size());
		}
		return dist;
	}

	void update_defaults(SemanticId sid, size_t childi, const typedist *D) {
		for (std::map<SemanticId, unsigned>::const_iterator i = D->D.begin(); i != D->D.end(); ++i) {
			selections[i->first] += i->second;
		}
		const std::set<SemanticId> &sids = sidslots[sid][childi];
		for (std::set<SemanticId>::const_iterator i = sids.begin(); i != sids.end(); ++i) {
			chances[*i] += D->nsamples;
		}
	}
};

ptab *
mkptab(const grammar *G) {
	ptab *pt = new ptab;
	return pt;
}

void
rmptab(ptab *pt) {
	delete pt;
}

static int
processterm(ptab *pt, const typetree &tt, SemanticId parentsid, size_t childi) {
	ptab::tcond cond(parentsid, childi, InvalidSemanticId);
	for (std::vector<typetree>::const_iterator i = tt.children.begin(); i != tt.children.end(); ++i) {
		typedist *dist = pt->tdist(cond);
		dist->inc(i->sid);
		cond.pred = i->sid;
	}
	return 0;
}

static int
processtree(ptab *pt, const typetree &tt, SemanticId parentsid, size_t childi) {
	if (tt.sid == TERMINAL_EXPR) {
		return processterm(pt, tt, parentsid, childi);
	}
	else {
		size_t childi = 0;
		for (std::vector<typetree>::const_iterator i = tt.children.begin(); i != tt.children.end(); ++i, ++childi) {
			ptab::ntcond c(tt.sid, childi);
			pt->ntdist(c)->inc(i->sid);
			int e = processtree(pt, *i, tt.sid, childi);
			if (FAILURE(e)) return e;
		}
	}
	return 0;
}

int
ptab_insert(ptab *pt, const typetree &tt) {
	typetree wrap;
	wrap.children.push_back(tt);
	return processtree(pt, wrap, InvalidSemanticId, 0);
}

void
ptab_remove(ptab *pt, SemanticId sid) {
	for (std::map<ptab::ntcond, typedist *>::iterator i = pt->ntD.begin(); i != pt->ntD.end(); ++i) {
		ptab::ntcond cnd = i->first;
		if (cnd.parent == sid) {
			delete i->second;
			pt->ntD.erase(i);
			i = pt->ntD.upper_bound(cnd);
		}
		else {
			i->second->remove(sid);
		}
	}
	for (std::map<ptab::tcond, typedist *>::iterator i = pt->tD.begin(); i != pt->tD.end(); ++i) {
		ptab::tcond cnd = i->first;
		if (cnd.parent == sid || cnd.pred == sid) {
			delete i->second;
			pt->tD.erase(i);
			i = pt->tD.upper_bound(cnd);
		}
		else {
			i->second->remove(sid);
		}
	}
}

int
ptab_read(ptab *pt, std::istream &is) {
	size_t n;
	is >> n;
	if (is.eof() || !is) {
		return E_IO;
	}
	while (n--) {
		std::string s;
		is >> s;
		if (is.eof()) break;
		if (is.fail() || is.bad()) return E_IO;

		SemanticId sid = string_to_sid(s);
		unsigned childi;
		is >> childi;
		typedist *D = pt->ntdist(ptab::ntcond(sid, childi));
		int e = read_typedist(D, is);
		if (FAILURE(e)) {
			return e;
		}
		pt->update_defaults(sid, childi, D);
	}

	is >> n;
	if (is.eof() || !is) {
		return E_IO;
	}
	while (n--) {
		std::string s;
		is >> s;
		if (is.eof()) break;
		if (is.fail() || is.bad()) return E_IO;
		SemanticId parentsid = string_to_sid(s);
		unsigned childi;
		is >> childi;
		std::string p;
		is >> p;
		SemanticId predsid = string_to_sid(p);
		typedist *D = pt->tdist(ptab::tcond(parentsid, childi, predsid));
		int e = read_typedist(D, is);
		if (FAILURE(e)) {
			return e;
		}
		pt->update_defaults(parentsid, childi, D);
	}
	return 0;
}

void
ptab_write(ptab *pt, std::ostream &os) {
	os << pt->ntD.size() << std::endl;
	for (std::map<ptab::ntcond, typedist *>::const_iterator i = pt->ntD.begin(); i != pt->ntD.end(); ++i) {
		os << sid_to_string(i->first.parent) << ' ' << i->first.i << std::endl;
		os << *i->second << std::endl;
	}
	os << pt->tD.size() << std::endl;
	for (std::map<ptab::tcond, typedist *>::const_iterator i = pt->tD.begin(); i != pt->tD.end(); ++i) {
		os << sid_to_string(i->first.parent) << ' ' << i->first.childi << ' ' << sid_to_string(i->first.pred) << std::endl;
		os << *i->second << std::endl;
	}
}

double
ptab_measureterms(const ptab *pt, const typetree &tt, SemanticId parentsid, size_t childi) {
	double sc = 1.0;
	SemanticId predsid = InvalidSemanticId;
	for (size_t i = 0; i < tt.children.size(); ++i) {
		const typetree &child = tt.children[i];
		sc *= pt->termP(parentsid, childi, predsid, child.sid);
		predsid = child.sid;
	}
	return sc;
}

double
ptab_measuretop(const ptab *pt, const ExpressionTree *tree) {
	typetree tt;
	mktypetree(tt, tree);
	if (tt.sid == TERMINAL_EXPR) {
		//return ptab_measureterms(pt, tt, InvalidSemanticId, 0);
		return 1.0;
	}
	double sc = 1.0;
	for (size_t i = 0; i < tt.children.size(); ++i) {
		const typetree &child = tt.children[i];
		if (child.sid == TERMINAL_EXPR) {
			sc *= ptab_measureterms(pt, child, tt.sid, i);
		}
		else {
			sc *= pt->ntP(tt.sid, i, child.sid);
		}
	}
	return sc;
}

}
