#include <vector>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include "lconv.h"
#include "cmdcode.h"
#include "MathRecognizer.h"
#include "mathrecognizer-private.h"
#include "symbols.h"
#include "relation.h"
#include "order.h"
#include "rect.h"

#include "log.h"
#include "expr.h"

const int SPLITEXPR = 42;

struct fixed_relation : public scg::relation {
	std::set<std::pair<scg::Rect<long>, scg::Rect<long> > > members;

	explicit fixed_relation(const std::string &name) : relation(name, 0, 0, 0, 0) { }

	double membership(const scg::interpretation *tail, const scg::interpretation *head, 
	                  scg::rclass_t, scg::rclass_t, int, int) const {
		/*std::cout << "TESTING " << name << " on ";
		scg::Rect<long> h = head->bounds();
		std::cout << "  (" << h.left << ',' << h.top << ")-(" << h.right << ',' << h.bottom;
		scg::Rect<long> t = tail->bounds();
		std::cout << ")  (" << t.left << ',' << t.top << ")-(" << t.right << ',' << t.bottom;*/
		double val = (members.find(std::make_pair(tail->bounds(), head->bounds())) == members.end()) ? 0 : 1;
		//std::cout << ") == " << val << std::endl;
		return val;
	}

	void write(std::ostream &os) const {
		os << name << ":\n";
		for (std::set<std::pair<scg::Rect<long>, scg::Rect<long> > >::const_iterator i = members.begin(); i != members.end(); ++i) {
			scg::Rect<long> head = i->first;
			os << "  (" << head.left << ',' << head.top << ")-(" << head.right << ',' << head.bottom;
			scg::Rect<long> tail = i->second;
			os << ")  (" << tail.left << ',' << tail.top << ")-(" << tail.right << ',' << tail.bottom << ")\n";
		}
	}
};

static fixed_relation AR("AR"), R("R"), BR("BR"), B("B"), C("C");
static fixed_relation *rels[] = { &AR, &R, &BR, &B, &C };


struct treenode {
	static size_t NN;
	static void resetN() { NN = 0; }

	int rel;
	// only one of the next two lines.
	std::string name; size_t n;

	std::vector<treenode> children; size_t head, tail;

	explicit treenode(int rel_ = scg::REL_RIGHT) : rel(rel_), n(~0),  head(~0), tail(~0) { }
	void addchild(const treenode &tn) {
		//assert(tn.isvalid());
		if (!tn.isvalid()) return;
		if (children.empty()) {
			head = 0;
		}
		tail = children.size();
		children.push_back(tn);
	}
	bool isvalid() const {
		return ((!name.empty() && n != ~0) || (!children.empty() && head != ~0 && tail != ~0));
	}
};
size_t treenode::NN = 0;

std::ostream &
operator<<(std::ostream &os, const treenode &tn) {
	if (tn.children.empty()) {
		os << tn.name;
	}
	else {
		os << scg::get_relation(tn.rel)->name << '(';
		for (std::vector<treenode>::const_iterator i = tn.children.begin(); i != tn.children.end(); ++i) {
			if (i != tn.children.begin()) os << ' ';
			os << *i;
		}
		os << ')';
	}
	return os;
}

static int readrow(const char *latex, size_t &i, treenode &tn);
static int readroworterm(const char *latex, size_t &i, treenode &tn);


static treenode
mktermnode(const std::string &name) {
	treenode tn;
	tn.name = name;
	tn.n = treenode::NN;
	scg::symbol *S = scg::symdb_findsymbol_name(name);
	treenode::NN += 1;
	return tn;
}


struct cmdspec {
	std::string cmd;
	int (*fn)(const char *, size_t &, treenode &);
};

static int ignore(const char *, size_t &, treenode &tn) { return 0; }
static int
ignorearg(const char *latex, size_t &i, treenode &tn) {
	treenode arg;
	int e = readroworterm(latex, i, arg);
	if (e != 0) return e;
	return 0;
}

static int mkalpha(const char *, size_t &, treenode &tn) { tn = mktermnode("alpha"); return 0; }
static int mkbeta(const char *, size_t &, treenode &tn) { tn = mktermnode("beta"); return 0; }
static int mkgamma(const char *, size_t &, treenode &tn) { tn = mktermnode("gamma"); return 0; }
static int mkdelta(const char *, size_t &, treenode &tn) { tn = mktermnode("delta"); return 0; }
static int mkepsilon(const char *, size_t &, treenode &tn) { tn = mktermnode("epsilon"); return 0; }
static int mkzeta(const char *, size_t &, treenode &tn) { tn = mktermnode("zeta"); return 0; }
static int mktheta(const char *, size_t &, treenode &tn) { tn = mktermnode("theta"); return 0; }
static int mklambda(const char *, size_t &, treenode &tn) { tn = mktermnode("lambda"); return 0; }
static int mkmu(const char *, size_t &, treenode &tn) { tn = mktermnode("mu"); return 0; }
static int mkxi(const char *, size_t &, treenode &tn) { tn = mktermnode("xi"); return 0; }
static int mkpi(const char *, size_t &, treenode &tn) { tn = mktermnode("pi"); return 0; }
static int mkrho(const char *, size_t &, treenode &tn) { tn = mktermnode("rho"); return 0; }
static int mksigma(const char *, size_t &, treenode &tn) { tn = mktermnode("sigma"); return 0; }
static int mktau(const char *, size_t &, treenode &tn) { tn = mktermnode("tau"); return 0; }
static int mkphi(const char *, size_t &, treenode &tn) { tn = mktermnode("phi"); return 0; }
static int mkpsi(const char *, size_t &, treenode &tn) { tn = mktermnode("psi"); return 0; }
static int mkomega(const char *, size_t &, treenode &tn) { tn = mktermnode("omega"); return 0; }
static int mkDelta(const char *, size_t &, treenode &tn) { tn = mktermnode("Delta"); return 0; }
static int mkGamma(const char *, size_t &, treenode &tn) { tn = mktermnode("Gamma"); return 0; }
static int mkSigma(const char *, size_t &, treenode &tn) { tn = mktermnode("Sigma"); return 0; }
static int mkOmega(const char *, size_t &, treenode &tn) { tn = mktermnode("Omega"); return 0; }
static int mkPi(const char *, size_t &, treenode &tn) { tn = mktermnode("Pi"); return 0; }

static int mkd(const char *, size_t &, treenode &tn) { tn = mktermnode("d"); return 0; }
static int mklbrace(const char *, size_t &, treenode &tn) { tn = mktermnode("lbrace"); return 0; }
static int mkrbrace(const char *, size_t &, treenode &tn) { tn = mktermnode("rbrace"); return 0; }

static int mkleq(const char *, size_t &, treenode &tn) { tn = mktermnode("leq"); return 0; }
static int mkgeq(const char *, size_t &, treenode &tn) { tn = mktermnode("geq"); return 0; }
static int mkneq(const char *, size_t &, treenode &tn) { tn = mktermnode("neq"); return 0; }
static int mkell(const char *, size_t &, treenode &tn) { tn = mktermnode("l"); return 0; }
static int mkprime(const char *, size_t &, treenode &tn) { tn = mktermnode("prime"); return 0; }
static int mktimes(const char *, size_t &, treenode &tn) { tn = mktermnode("times"); return 0; }
static int mkIntegral(const char *, size_t &, treenode &tn) { tn = mktermnode("Integral"); return 0; }
static int mkcdot(const char *, size_t &, treenode &tn) { tn = mktermnode("cdot"); return 0; }
static int mkinfin(const char *, size_t &, treenode &tn) { tn = mktermnode("infin"); return 0; }
static int mkplusorminus(const char *latex, size_t &i, treenode &tn) { tn = mktermnode("plusorminus"); return 0; }
static int mkarrow(const char *latex, size_t &i, treenode &tn) { tn = mktermnode("arrow"); return 0; }
static int mkequiv(const char *latex, size_t &i, treenode &tn) { tn = mktermnode("equiv"); return 0; }
static int mkapprox(const char *latex, size_t &i, treenode &tn) { tn = mktermnode("approx"); return 0; }

static int
mkfrac(const char *latex, size_t &i, treenode &tn) {
	treenode frac(scg::REL_BELOW);
	treenode numer, denom;
	int e = readroworterm(latex, i, numer);
	if (e != 0) return e;
	e = readroworterm(latex, i, denom);
	if (e != 0) return e;
	frac.addchild(numer);
	frac.addchild(mktermnode("horzline"));
	frac.addchild(denom);
	frac.head = frac.tail = 1;
	tn.addchild(frac);
	return 0;
}

static int
mkodiff(const char *latex, size_t &i, treenode &tn) {
	treenode frac(scg::REL_BELOW);
	treenode numer, denom;
	int e = readroworterm(latex, i, numer);
	if (e != 0) return e;
	e = readroworterm(latex, i, denom);
	if (e != 0) return e;
	treenode numer_wrap, denom_wrap;
	numer_wrap.addchild(mktermnode("d"));
	numer_wrap.addchild(numer);
	denom_wrap.addchild(mktermnode("d"));
	denom_wrap.addchild(denom);
	frac.addchild(numer_wrap);
	frac.addchild(mktermnode("horzline"));
	frac.addchild(denom_wrap);
	frac.head = frac.tail = 1;
	tn.addchild(frac);
	return 0;
}


static int
mksqrt(const char *latex, size_t &i, treenode &tn) {
	treenode root(scg::REL_CONTAINS);
	root.addchild(mktermnode("sqrt"));
	while (std::isspace(latex[i])) ++i;
	if (latex[i] == '[') {
		std::cerr << "mktrees: cannot parse decorated root signs\n";
		return CASREP_CONVERR;
	}
	treenode arg;
	int e = readroworterm(latex, i, arg);
	if (e != 0) return e;
	root.addchild(arg);
	tn.addchild(root);
	return 0;
}

static int
mkname(const char *n, treenode &tn) {
	while (*n) {
		treenode child = mktermnode(std::string(1, *n));
		tn.addchild(child);
		++n;
	}
	return 0;
}

static int extractarg(const char *latex, size_t &i, treenode &tn) { return readroworterm(latex, i, tn); }

static int
mkcdots(const char *latex, size_t &i, treenode &tn) {
	treenode *node;
	if (tn.rel == scg::REL_RIGHT) {
		node = &tn;
	}
	else {
		tn.addchild(treenode());
		node = &tn.children.back();
	}
	node->addchild(mktermnode("dot"));
	node->addchild(mktermnode("dot"));
	node->addchild(mktermnode("dot"));
	return 0;
}

static int mksin(const char *latex, size_t &i, treenode &tn) { return mkname("sin", tn); }
static int mkcos(const char *latex, size_t &i, treenode &tn) { return mkname("cos", tn); }
static int mktan(const char *latex, size_t &i, treenode &tn) { return mkname("tan", tn); }
static int mksec(const char *latex, size_t &i, treenode &tn) { return mkname("sec", tn); }
static int mkcsc(const char *latex, size_t &i, treenode &tn) { return mkname("csc", tn); }
static int mkcot(const char *latex, size_t &i, treenode &tn) { return mkname("cot", tn); }
static int mkarcsin(const char *latex, size_t &i, treenode &tn) { return mkname("arcsin", tn); }
static int mkarccos(const char *latex, size_t &i, treenode &tn) { return mkname("arccos", tn); }
static int mkarctan(const char *latex, size_t &i, treenode &tn) { return mkname("arctan", tn); }
static int mkerf(const char *latex, size_t &i, treenode &tn) { return mkname("erf", tn); }
static int mklog(const char *latex, size_t &i, treenode &tn) { return mkname("log", tn); }
static int mkexp(const char *latex, size_t &i, treenode &tn) { return mkname("exp", tn); }
static int mkln(const char *latex, size_t &i, treenode &tn) { return mkname("ln", tn); }
static int mklim(const char *latex, size_t &i, treenode &tn) { return mkname("lim", tn); }
static int mkarg(const char *latex, size_t &i, treenode &tn) { return mkname("arg", tn); }
static int mkRe(const char *latex, size_t &i, treenode &tn) { return mkname("Re", tn); }
static int mkIm(const char *latex, size_t &i, treenode &tn) { return mkname("Im", tn); }
static int mkmin(const char *latex, size_t &i, treenode &tn) { return mkname("min", tn); }
static int mkmax(const char *latex, size_t &i, treenode &tn) { return mkname("max", tn); }
static int mksign(const char *latex, size_t &i, treenode &tn) { return mkname("sign", tn); }
static int mkdet(const char *latex, size_t &i, treenode &tn) { return mkname("det", tn); }
static int mkst(const char *latex, size_t &i, treenode &tn) { return mkname("st", tn); }

static int
splitexpr(const char *latex, size_t &i, treenode &tn) {
	return SPLITEXPR;
}

static int
readspec(const char *latex, size_t &i, treenode &tn) {
	const static cmdspec specs[] = {
		{ "alpha", &mkalpha }, { "beta", &mkbeta }, { "gamma", &mkgamma },
		{ "delta", &mkdelta }, { "epsilon", &mkepsilon }, { "varepsilon", &mkepsilon },
		{ "zeta", &mkzeta }, { "theta", &mktheta }, { "vartheta", &mktheta },
		{ "lambda", &mklambda }, { "mu", &mkmu }, { "xi", &mkxi },
		{ "pi", &mkpi }, { "varpi", &mkpi }, { "psi", &mkpsi },
		{ "phi", &mkphi }, { "varphi", &mkphi }, { "Gamma", &mkGamma },
		{ "Delta", &mkDelta }, { "omega", &mkomega}, { "sigma", &mksigma },
		{ "Pi", &mkPi }, { "Omega", &mkOmega }, { "Sigma", &mkSigma },
		{ "leq", &mkleq }, { "geq", &mkgeq }, { "neq", &mkneq },
		{ "ne", &mkneq }, { "le", &mkleq }, { "ge", &mkgeq },
		{ "equiv", &mkequiv }, { "der", &mkd }, { "st", &mkst },
		{ "ell", &mkell }, { "prime", &mkprime }, { "zu", &extractarg },
		{ "frac", &mkfrac }, { "tfrac", &mkfrac }, { "dfrac", &mkfrac },
		{ "odiff", &mkodiff }, { "approx", &mkapprox },
		{ "sum", &mkSigma }, { "int", &mkIntegral }, { "times", &mktimes },
		{ "{", &mklbrace }, { "}", &mkrbrace },
		{ "sin", &mksin }, { "cos", &mkcos }, { "tan", &mktan },
		{ "sec", &mksec }, { "csc", &mkcsc }, { "cot", &mkcot },
		{ "arcsin", &mkarcsin }, { "arccos", &mkarccos }, { "arctan", &mkarctan },
		{ "log", &mklog }, { "exp", &mkexp }, { "ln", &mkln }, { "erf", &mkerf },
		{ "Re", &mkRe }, { "Im", &mkIm }, { "min", &mkmin }, { "max", &mkmax },
		{ "mn", &mkmin }, { "mx", &mkmax },
		{ "sqrt", &mksqrt }, { "cdot", &mkcdot }, { "cdots", &mkcdots },
		{ "dd", &mkd }, { "infty", &mkinfin }, { "pm", &mkplusorminus },
		{ "to", &mkarrow }, { "nearrow", &mkarrow }, { "searrow", &mkarrow },
		{ "lim", &mklim }, { "arg", &mkarg },
		{ "sign", &mksign }, { "det", &mkdet },
		{ "mathrm", &extractarg }, { "mathcal", &extractarg },
		{ "texttrm", &extractarg }, { "vphantom", &extractarg },
		
		{ "iff", &splitexpr }, { "implies", &splitexpr }, { "isdef", &splitexpr },

		{ "l", &ignore }, { "r", &ignore }, { "bigl", &ignore }, { "bigr", &ignore },
		{ "Bigl", &ignore }, { "Bigr", &ignore }, { "ddagger", &mkd },
		{ "quad", &ignore }, { "qquad", &ignore }, { ",", &ignore }, { ";", &ignore }, { ":", &ignore },
		{ "left", &ignore }, { "right", &ignore },
		{ "hline", &ignore }, { "left.", &ignore }, { "right.", &ignore },
		{ "big", &ignore }, { "bigl.", &ignore }, { "bigr.", &ignore },
		{ "DS", &ignore }, { "displaystyle", &ignore }, { "ds", &ignore },
		{ "Ti", &ignore }, { "notag", &ignore },
		{ "texttrm", &ignorearg }, { "eqno", &ignorearg }, { "tag", &ignorearg },
	};
	const size_t NSPECS = sizeof(specs)/sizeof(*specs);

	std::vector<bool> avail(NSPECS, true);
	size_t navail = NSPECS;
	size_t match = ~0;
	size_t at = 0;
	while (navail > 0 && latex[i+at] != '\0') {
		for (size_t j = 0; j < NSPECS; ++j) {
			if (avail[j]) {
				if (latex[i+at] != specs[j].cmd[at]) {
					avail[j] = false;
					--navail;
				}
				else if (at+1 == specs[j].cmd.length()) {
					match = j;
				}
			}
		}
		++at;
	}

	if (match == ~0) {
		std::cerr << "mktrees: unknown latex command at " << i << std::endl;
		return CASREP_CONVERR;
	}

	i += specs[match].cmd.length();
	return (*specs[match].fn)(latex, i, tn);
}

static int
readterm(const char *latex, size_t &i, treenode &tn) {
	if (latex[i] == '\\') {
		++i;
		return readspec(latex, i, tn);
	}
	else if (std::isdigit(latex[i]) || std::isalpha(latex[i])) {
		tn = mktermnode(std::string(1, latex[i]));
		++i;
		return 0;
	}
	else if (latex[i] == ';') { ++i; return 0; } // XXX: is this okay?
	else if (latex[i] == '\'') { ++i; tn = mktermnode("prime"); return 0; }
	else if (latex[i] == '(') { ++i; tn = mktermnode("lparen"); return 0; }
	else if (latex[i] == ')') { ++i; tn = mktermnode("rparen"); return 0; }
	else if (latex[i] == '[') { ++i; tn = mktermnode("lbracket"); return 0; }
	else if (latex[i] == ']') { ++i; tn = mktermnode("rbracket"); return 0; }
	else if (latex[i] == '!') { ++i; tn = mktermnode("exclaim"); return 0; }
	else if (latex[i] == '=') { ++i; tn = mktermnode("eq"); return 0; }
	else if (latex[i] == '<') { ++i; tn = mktermnode("lt"); return 0; }
	else if (latex[i] == '>') { ++i; tn = mktermnode("gt"); return 0; }
	else if (latex[i] == '+') { ++i; tn = mktermnode("plus"); return 0; }
	else if (latex[i] == '-') { ++i; tn = mktermnode("horzline"); return 0; }
	else if (latex[i] == '/') { ++i; tn = mktermnode("slash"); return 0; }
	else if (latex[i] == '.') { ++i; tn = mktermnode("dot"); return 0; }
	else {
		std::cerr << "mktrees: unsupported term char " << latex[i] << " at " << i << std::endl;
		return CASREP_CONVERR;
	}
}

static int
readroworterm(const char *latex, size_t &i, treenode &tn) {
	while (std::isspace(latex[i])) ++i;
	if (latex[i] == '{') {
		++i;
		int e = readrow(latex, i, tn);
		if (e != 0) return e;
		if (latex[i] != '}') {
			std::cerr << "mktrees: unmatched brace: expected \'}\', read " << latex[i] << std::endl;
			return CASREP_CONVERR;
		}
		++i;
		return 0;
	}
	else if (latex[i] == '\0') {
		std::cerr << "mkterms: unexpected end of string at pos " << i << std::endl;
		return CASREP_CONVERR;
	}
	else {
		return readterm(latex, i, tn);
	}
}

static int
readrow(const char *latex, size_t &i, treenode &tn) {
	while (latex[i] != '\0' && latex[i] != '}') {
		if (latex[i] == '{') {
			++i;
			treenode child;
			int e = readrow(latex, i, child);
			if (e != 0) return e;
			tn.addchild(child);
			if (latex[i] != '}') {
				std::cerr << "mktrees: unmatched brace: expected \'}\', read " << latex[i] << std::endl;
				return CASREP_CONVERR;
			}
			++i;
		}
		else if (latex[i] == '_') {
			if (tn.children.empty()) {
				std::cerr << "mktrees: no head on sub at pos " << i << std::endl;
				return CASREP_CONVERR;
			}
			++i;
			treenode base = tn.children.back();
			tn.children.back() = treenode(scg::REL_BELOWRIGHT);
			treenode &sub = tn.children.back();
			sub.addchild(base);
			treenode child;
			int e = readroworterm(latex, i, child);
			if (e != 0) return e;
			sub.addchild(child);
			sub.tail = 0;
		}
		else if (latex[i] == '^') {
			if (tn.children.empty()) {
				std::cerr << "mktrees: no head on sup at pos " << i << std::endl;
				return CASREP_CONVERR;
			}
			++i;
			treenode base = tn.children.back();
			tn.children.back() = treenode(scg::REL_ABOVERIGHT);
			treenode &sup = tn.children.back();
			sup.addchild(base);
			treenode child;
			int e = readroworterm(latex, i, child);
			if (e != 0) return e;
			sup.addchild(child);
			sup.tail = 0;
		}
		else if (!std::isspace(latex[i])) {
			treenode child;
			int e = readterm(latex, i, child);
			if (e != 0) return e;
			tn.addchild(child);
		}
		else {
			++i;
		}
	}
}

void
reset_relations() {
	for (size_t i = 0; i < scg::NUM_RELATIONS; ++i) {
		scg::set_relation(i, rels[i]);
	}
}

void
clear_relations() {
	for (size_t i = 0; i < scg::NUM_RELATIONS; ++i) {
		rels[i]->members.clear();
		scg::set_relation(i, 0);
	}
}

static int
mktreenode(const char *s, treenode &tn) {
	int e;
	size_t i = 0;
	while (s[i] != '\0' && std::isspace(s[i])) ++i;
	if (s[i] == '\0') {
		std::cerr << "empty latex string\n";
		return CASREP_CONVERR;
	}
	e = readrow(s, i, tn);
	// No splits allowed in sage expressions.
	// Each command can only return a single expression.
	/*while (e == SPLITEXPR) {
		e = writetree(tn);
		if (e != 0) {
			return 0;
		}
		tn = treenode();
		treenode::resetN();
		e = readrow(s, i, tn);
	}*/
	return e;
}

static const treenode & 
gettail(const treenode &tn) {
	if (tn.children.empty()) {
		assert(tn.n != ~0);
		return tn;
	}
	assert(tn.tail != ~0);
	return gettail(tn.children[tn.tail]);
}

static const treenode &
gethead(const treenode &tn) {
	if (tn.children.empty()) {
		assert(tn.n != ~0);
		return tn;
	}
	assert(tn.head != ~0);
	return gethead(tn.children[tn.head]);
}

static int
mksymbols(const treenode &tn, long &x, long &y, scg::MathSymbol *symb, scg::Rect<long> *boxes) {
	if (tn.children.empty()) {
		assert(!tn.name.empty());
		scg::symbol *S = scg::symdb_findsymbol_name(tn.name);
		if (!S) {
			return CASREP_CONVERR;
		}
		boxes[tn.n].left = boxes[tn.n].right = x;
		boxes[tn.n].top = boxes[tn.n].bottom = y;
		symb[tn.n].symbol = S->unicode;
		symb[tn.n].bounds = boxes[tn.n];
	}
	else {
		assert(!tn.children.empty());
		for (size_t i = 0; i < tn.children.size(); ++i) {
			long x0 = x;
			long y0 = y;
			mksymbols(tn.children[i], x, y, symb, boxes);
			if (scg::ordering_for_relation(scg::get_relation(tn.rel)) == scg::Y_ORDER) {
				x = x0;
				++y;
			}
			else {
				++x;
				y = y0;
			}
		}
		for (size_t i = 0; i < tn.children.size()-1; ++i) {
			const treenode &tail = gettail(tn.children[i]);
			const treenode &head = gethead(tn.children[i+1]);
			rels[tn.rel]->members.insert(std::make_pair(boxes[tail.n], boxes[head.n]));
		}
	}
	return 0;
}

/*
exprtree *
tree2stree(const scg::ExpressionTree *tree) {
	if (tree->type() == scg::TERMINAL_EXPR) {
		exprleaf *leaf = new exprleaf;
		const wchar_t *s = tree->wstr();
		while (*s) {
			leaf->addterm((unsigned short)*s);
			++s;
		}
	}
	else {
		exprgrp *grp = new exprgrp(tree->type());
		for (size_t i = 0; i < tree->nchildren(); ++i) {
			exprtree *child = tree2stree(tree->child(i));
			if (!child) {
				delete grp;
				return 0;
			}
			grp->addchild(child);
		}
		return grp;
	}
}*/

const scg::ExpressionTree *
latex2tree(const char *s, int &stat) {
	reset_relations();
	treenode::resetN();

	treenode tn;
	stat = mktreenode(s, tn);
	logmsg("created treenode with return %d\n", stat);
	if (stat != 0) return 0;

	scg::MathSymbol *syms = new scg::MathSymbol[treenode::NN];
	scg::Rect<long> *boxes = new scg::Rect<long>[treenode::NN];
	long x = 0, y = 0;
	mksymbols(tn, x, y, syms, boxes);
	logmsg("created %u symbols and boxes\n", treenode::NN);

	scg::MathRecognizer *rec = scg::CreateMathRecognizer(syms, treenode::NN);
	scg::ExpressionIterator *it = rec->CreateDefaultIterator();
	const scg::ExpressionTree *tree = 0;
	if (it) {
		tree = it->next();
		stat = tree ? 0 : CASREP_CONVERR;
	}
	else {
		stat = CASREP_CONVERR;
	}

	delete[] syms;
	delete[] boxes;
	if (tree) logmsg("lconv gave tree %s\n", tree->str());
	else logmsg("lconv failed\n");

	clear_relations();

	stat = CASREP_OK;
	return tree;
}
