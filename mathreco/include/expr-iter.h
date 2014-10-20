#ifndef EXPR_ITER_H_
#define EXPR_ITER_H_

#include "grammar.h"
#include "intrpr.h"
#include "ordered-segments.h"
#include "mathrecognizer-private.h"
#include "MatrixAnalyzer.h"

#include <set>
#include <vector>

namespace scg
{

class external_iterator;

class treeparser : public parser {
public:
	treeparser(math_recognizer_base *rec, const production *P_, const ordered_segments *segs_, std::vector<interpreter *> &parsers);
	treeparser();
	void read(reader &re, math_recognizer_base *rec);
	~treeparser();

	char id() const { return ID; }

private:
	bool addtree(const std::vector<size_t> &ranks);
	interpretation *getnext();

private:
	struct treerec {
		interpretation *node;
		std::vector<size_t> ranks;
		treerec() : node(0) { }
		explicit treerec(interpretation *node_, const std::vector<size_t> &ranks_) : node(node_), ranks(ranks_) { }
		bool operator>(const treerec &rhs) const {
			return node->score() > rhs.node->score();
		}
	};

private:
	void freeknown();
	void reset();

public:
	void write(writer &wr) const;

private:
	static void read_treerec(reader &re, treerec &tr, math_recognizer_base *rec) {
		re.read(&tr.node, rec);
		read_ranks(re, tr.ranks);
	}

	static void write_treerec(writer &wr, const treerec &tr) {
		wr.write(tr.node);
		write_ranks(wr, tr.ranks);
	}

	static void read_ranks(reader &re, std::vector<size_t> &ranks) {
		size_t n;
		re.read(n);
		ranks.resize(n);
		for (size_t i = 0; i < n; ++i) {
			re.read(ranks[i]);
		}
	}

	static void write_ranks(writer &wr, const std::vector<size_t> &ranks) {
		wr.write(ranks.size());
		for (std::vector<size_t>::const_iterator i = ranks.begin(); i != ranks.end(); ++i) {
			wr.write(*i);
		}
	}

private:
	const production *P;
	typedef std::set<treerec, std::greater<treerec> > known_set;
	known_set known;
	treerec last;
	std::set<std::vector<size_t> > known_ranks;
public:
	static const char ID;
};


class matrixbuilder : public parser {
public:
	matrixbuilder(math_recognizer_base *rec_, MatrixAnalyzer *mxan_, MatrixIterator *mxiter_,
	              const nonterminal *nt, const ordered_segments *segs);
	matrixbuilder();
	~matrixbuilder();

	MatrixAnalyzer *mxan() { return mxan_; }

protected:
	interpretation *getnext();
	void reset() { abort(); }

public:
	char id() const { abort(); return 0; }

protected:
	MatrixAnalyzer *mxan_;
	MatrixIterator *mxiter;
	std::set<external_iterator *> extchildren;

private:
	struct matrix_rows_tbuilder : public tree_builder {
		basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
	};

	struct matrix_row_tbuilder : public tree_builder {
		basic_tree *build(const interpretation *src, subtree_accessor &acc) const;
	};

	struct matrix_rows_mathml : public string_builder {
		std::string build(const ExpressionTree *node) const;
	};

	struct matrix_row_mathml : public string_builder {
		std::string build(const ExpressionTree *node) const;
	};

	struct matrix_rows_latex : public string_builder {
		std::string build(const ExpressionTree *node) const;
	};

	struct matrix_row_latex : public string_builder {
		std::string build(const ExpressionTree *node) const;
	};

	static production *matrix_prod;
	static production *matrix_rows_prod;
	static production *matrix_row_prod;

private:
	friend int InitializeMatrixRecognizer(grammar *G);
	friend int ShutdownMatrixRecognizer();
};


class matrixparser : public matrixbuilder {
public:
	matrixparser(math_recognizer_base *rec_, const nonterminal *nt, const ordered_segments *segs, int caps = MatrixAnalyzer::RECOGNIZE_FULL);
	matrixparser();
	void read(reader &re, math_recognizer_base *rec_);
	~matrixparser();

	void write(writer &wr) const;

protected:
	void initmatrixreco(int caps);
	
public:
	char id() const { return ID; }
	static const char ID;

protected:
	void reset();

private:
	int capabilities;
};


class matrixintrp : public interpretation {
public:
	static const char ID;
	char id() const { return ID; }

public:
	matrixintrp(matrixbuilder *src_, const production *P, const ordered_segments *segs);
	matrixintrp(reader &re, math_recognizer_base *rec);
	~matrixintrp();
	bool haslongform() const;
	interpreter *longformiter() const;
	void write(writer &wr) const;

private:
	mutable matrixbuilder *src;
	mutable LongFormMatrixIterator *mxit;
	mutable interpreter *intrpr;
};

int InitializeMatrixRecognizer(grammar *G);
int ShutdownMatrixRecognizer();

}


#endif
