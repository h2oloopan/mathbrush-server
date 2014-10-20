#ifndef INTRPR_H_
#define INTRPR_H_

#include <vector>
#include <map>
#include <queue>
#include <set>
#include "grammar-fwd.h"
#include "ordered-segments.h"
#include "iohelp.h"

typedef unsigned rev_t;

namespace scg {

class interpretation;
class math_recognizer_base;

class interpreter {
public:
	interpreter(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span, int priority);
	interpreter();
	virtual void read(reader &re, math_recognizer_base *rec);
	virtual ~interpreter();
	
	virtual char id() const = 0;

	virtual size_t nknown() const = 0;
	virtual interpretation *nth(size_t i) const = 0;
	virtual interpretation *next() = 0;

public:
	inline math_recognizer_base *ctx() const { return rec; }
	inline const nonterminal *nt() const { return nt_; }
	inline const ordered_segments *span() const { return span_; }

	virtual void write(writer &wr) const;

	void sethint(double h) { hint = h; }
	void addhint(double h) { hint += h; }
	void rmhint(double h) { hint -= h; }

	int lock();
	int unlock();
	bool islocked() const;

	int resetpriority;
	void resetparents();

protected:
	virtual void reset() = 0;
	void addreset();

	inline interpreter *child(size_t i) { return children[i]; }
	inline const interpreter *child(size_t i) const { return children[i]; }
	inline size_t nchildren() const { return children.size(); }
	void addchild(interpreter *child, bool rm);
	void rmchild(interpreter *child);

private:
	friend class math_recognizer_base;
	math_recognizer_base *rec;
	const nonterminal *nt_;
	const ordered_segments *span_;
	std::vector<interpreter *> children;
	std::vector<bool> rmflags;

	size_t nlinks;
	void countlinks();
	
public:
	std::set<interpreter *> parents;
protected:
	double hint;
	void clearchildren();
};

class staticintrpr : public interpreter {
public:
	staticintrpr(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span);
	staticintrpr();
	void read(reader &re, math_recognizer_base *rec);
	~staticintrpr();
	char id() const { return ID; }
	size_t nknown() const;
	interpretation *nth(size_t i) const;
	interpretation *next();
	void addknown(interpretation *t, bool rm);
	size_t nadded() const;

	void write(writer &wr) const;

protected:
	void reset();
private:
	std::vector<interpretation *> ts;
	std::vector<bool> rmflags;
	size_t at;
public:
	static const char ID;
};

class parser : public interpreter {
public:
	parser(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span, int priority);
	parser();
	void read(reader &re, math_recognizer_base *rec);
	size_t nknown() const;
	interpretation *nth(size_t i) const;
	interpretation *next();

	void write(writer &wr) const;

protected:
	virtual interpretation *getnext() = 0;
	void freeknown();
	void clearknown();
private:
	std::vector<interpretation *> known;
};

/*
class iterator : public interpreter {
public:
	iterator(interpreter *src, bool ownptr = false);
	iterator();
	void read(reader &re, math_recognizer_base *rec);
	char id() const { return ID; }
	size_t nknown() const;
	interpretation *nth(size_t i) const;
	interpretation *next();

	void write(writer &wr) const;

protected:
	void reset();
private:
	size_t curr;
public:
	static const char ID;
};*/

class multiplexor : public parser {
public:
	multiplexor(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span, int priority);
	multiplexor();
	void read(reader &re, math_recognizer_base *rec);
	char id() const { return ID; }
	void addparser(interpreter *prsr, bool rm = true, int priority = 0);

	void write(writer &wr) const;

protected:
	void reset();
	interpretation *getnext();
private:
	void addintrp(interpretation *intrp, size_t idx, int priority);
	std::vector<size_t> intrprpos;
	std::vector<int> priority;
	size_t prevsrc;

	struct intrp_t {
		interpretation *intrp;
		size_t src;
		int priority;

		intrp_t() : intrp(0), src(0), priority(0) { }
		intrp_t(interpretation *intrp, size_t src, int priority_);
		bool operator<(const intrp_t &rhs) const;
	};
	std::vector<intrp_t> known;
public:
	static const char ID;
};

class sequence : public parser {
public:
	sequence(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span, int priority);
	void addparser(interpreter *intrpr, bool rm = true);
	char id() const { return ID; }

	void read(reader &re, math_recognizer_base *rec);
	void write(writer &wr) const;

protected:
	interpretation *getnext();
	void reset();

private:
	size_t curintrpr;
	size_t curpos;
public:
	static const char ID;
};

interpretation *getnth(interpreter *intrpr, size_t n);
interpretation *getfirst(interpreter *intrpr);

}

#endif
