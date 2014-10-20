#include "expr-iter.h"
#include "builder.h"
#include "extern-iter.h"
#include "reco-types.h"

namespace scg {

const static double OP_FACTOR = std::log(0.95);

matrixintrp::matrixintrp(matrixbuilder *src_, const production *P, const ordered_segments *segs)
	: interpretation(src_, P, segs), src(src_), mxit(0), intrpr(0) {
}
matrixintrp::matrixintrp(reader &re, math_recognizer_base *rec)
	: interpretation(re, rec), mxit(0), intrpr(0) {
	interpreter *prox;
	re.read(&prox, rec);
	src = dynamic_cast<matrixbuilder *>(prox);
}

matrixintrp::~matrixintrp() {
	delete intrpr;
	delete mxit;
}

void
matrixintrp::write(writer &wr) const {
	interpretation::write(wr);
	wr.write(src);
}

bool
matrixintrp::haslongform() const {
	return src->mxan()->hasValidLongForm();
}

interpreter *
matrixintrp::longformiter() const {
	assert(src->mxan()->hasValidLongForm());
	// XXX: this should not be necessary!
	// possibly the problem with nested subtree
	// fixed_interpretation pointers getting screwed
	// up is related to the unspecified matrix object
	// rebuilding itself in the HasLongForm() query
	// that precedes iterator construction
	/*
	if (intrpr) {
		delete intrpr;
		intrpr = 0;
	}*/
	if (!intrpr) {
		mxit = src->mxan()->createLongFormIterator();
		assert(mxit);
		intrpr = new matrixbuilder(src->ctx(), src->mxan(), mxit, P->nt, span);
	}
	return intrpr;
}

const char matrixintrp::ID = 10;

matrixbuilder::matrixbuilder(math_recognizer_base *rec_, MatrixAnalyzer *mxan__, MatrixIterator *mxiter_,
                             const nonterminal *nt, const ordered_segments *segs)
	: parser(rec_, nt, segs, 3), mxan_(mxan__), mxiter(mxiter_) {
}

matrixbuilder::matrixbuilder() : parser() { }

matrixbuilder::~matrixbuilder() {
	for (std::set<external_iterator *>::iterator i = extchildren.begin(); i != extchildren.end(); ++i) {
		(*i)->release();
	}
}

interpretation *
matrixbuilder::getnext() {
	if (!mxiter->hasNext()) {
		return 0;
	}
	Matrix *mx = mxiter->next();
	assert(mx);
	double mxconf;
	mxan_->getConfidence(mx, mxconf);

	size_t nrows, ncols;
	mx->getDimensions(nrows, ncols);
	if (mxconf == 0 || (nrows == 1 && ncols == 1)) {
		return 0;
	}

	VERBOSE(*verb_out << "matrixparser: making interpretation with " << nrows << " rows and " << ncols << " cols\n");
	interpretation *intrp = new matrixintrp(this, matrix_prod, span());
	std::stringstream ss;
	ss << nrows;
	intrp->addmetachild(mkterminal(ss.str()), true);
	ss.str("");
	ss << ncols;
	intrp->addmetachild(mkterminal(ss.str()), true);
	interpretation *rows = new interpretation(this, matrix_rows_prod, 0);
	for (size_t i = 0; i < nrows; ++i) {
		interpretation *row = new interpretation(this, matrix_row_prod, 0);
		for (size_t j = 0; j < ncols; ++j) {
			ExpressionIterator *it;
			mx->getCell(i, j, it);
			if (!it) {
				delete intrp;
				return 0;
			}
			
			external_iterator *ext = dynamic_cast<external_iterator *>(it);
			interpreter *intrpr = ext->getsrc();
			addchild(intrpr, false);

			interpretation *intrp = getfirst(intrpr);
			VERBOSE(*verb_out << "matrixparser: cell (" << i << ',' << j << ") is " << intrp->ccstr() << std::endl);
			if (!intrp) {
				delete row;
				delete rows;
				delete intrp;
				return 0;
			}
			row->addchild(intrp, std::log(j > 0 ? mxconf : 1), false);
			row->scorebias += OP_FACTOR;
			extchildren.insert(ext);
		}
		rows->addchild(row, true);
	}
	intrp->addchild(rows, true);
	//intrp->set_score(intrp->score_combo().addop());
	intrp->rclass = 1 << AGGREGATE_CLASS;
	VERBOSE(*verb_out << "matrixparser:" << span()->bits << " found interpretation " << intrp->str() << " with score " << intrp->score() << std::endl);
	return intrp;
}

basic_tree *
matrixbuilder::matrix_rows_tbuilder::build(const interpretation *src, subtree_accessor &acc) const {
	assert(src->P == matrix_rows_prod);
	parsed_tree *tree = new parsed_tree(src);
	for (size_t i = 0; i < src->nchildren(); ++i) {
		tree->addchild(acc[i]);
		acc.markused(i);
	}
	VERBOSE(*verb_out << "matrix_rows_tbuilder built tree " << tree->str() << std::endl);
	return tree;
}

basic_tree *
matrixbuilder::matrix_row_tbuilder::build(const interpretation *src, subtree_accessor &acc) const {
	assert(src->P == matrix_row_prod);
	parsed_tree *tree = new parsed_tree(src);
	for (size_t i = 0; i < src->nchildren(); ++i) {
		tree->addchild(acc[i]);
		acc.markused(i);
	}
	VERBOSE(*verb_out << "matrix_row_tbuilder built tree " << tree->str() << std::endl);
	return tree;
}

std::string
matrixbuilder::matrix_rows_mathml::build(const ExpressionTree *node) const {
	std::stringstream ss;
	for (size_t i = 0; i < node->nchildren(); ++i) {
		ss << "<mtr>" << node->child(i)->long_str() << "</mtr>";
	}
	return ss.str();
}

std::string
matrixbuilder::matrix_row_mathml::build(const ExpressionTree *node) const {
	std::stringstream ss;
	for (size_t i = 0; i < node->nchildren(); ++i) {
		ss << "<mtd>" << node->child(i)->long_str() << "</mtd>";
	}
	return ss.str();
}

std::string
matrixbuilder::matrix_rows_latex::build(const ExpressionTree *node) const {
	std::stringstream ss;
	if (node->nchildren() > 0) {
		ss << node->child(0)->latex_str();
		for (size_t i = 1; i < node->nchildren(); ++i) {
			ss << " \\\\ " << node->child(i)->latex_str();
		}
	}
	return ss.str();
}

std::string
matrixbuilder::matrix_row_latex::build(const ExpressionTree *node) const {
	std::stringstream ss;
	for (size_t i = 0; i < node->nchildren(); ++i) {
		ss << node->child(i)->latex_str() << " ";
		if (i != node->nchildren() - 1) {
			ss << "& ";
		}
	}
	return ss.str();
}

int
InitializeMatrixRecognizer(scg::grammar *G) {
	matrixbuilder::matrix_prod = new uniqproduction(1);
	std::list<tree_builder *> builders;
	builders.push_back(new ref_tree_builder(0));
	builders.push_back(new ref_tree_builder(1));
	builders.push_back(new ref_tree_builder(2));
	matrixbuilder::matrix_prod->tbuild = new subtree_builder(builders);
	matrixbuilder::matrix_prod->sbuild = new ref_string_builder(2, strid::mathml);
	matrixbuilder::matrix_prod->lbuild = new ref_string_builder(2, strid::latex);
	matrixbuilder::matrix_prod->sid = MATRIX_EXPR;
	G->setcanonicalsidproduction(MATRIX_EXPR, matrixbuilder::matrix_prod);

	matrixbuilder::matrix_rows_prod = new uniqproduction(2);
	matrixbuilder::matrix_rows_prod->tbuild = new matrixbuilder::matrix_rows_tbuilder;
	matrixbuilder::matrix_rows_prod->sbuild = new matrixbuilder::matrix_rows_mathml;
	matrixbuilder::matrix_rows_prod->lbuild = new matrixbuilder::matrix_rows_latex;
	matrixbuilder::matrix_rows_prod->sid = MATRIXROWS_EXPR;

	G->setcanonicalsidproduction(MATRIXROWS_EXPR, matrixbuilder::matrix_rows_prod);

	matrixbuilder::matrix_row_prod = new uniqproduction(3);
	matrixbuilder::matrix_row_prod->tbuild = new matrixbuilder::matrix_row_tbuilder;
	matrixbuilder::matrix_row_prod->sbuild = new matrixbuilder::matrix_row_mathml;
	matrixbuilder::matrix_row_prod->lbuild = new matrixbuilder::matrix_row_latex;
	matrixbuilder::matrix_row_prod->sid = MATRIXROW_EXPR;
	G->setcanonicalsidproduction(MATRIXROW_EXPR, matrixbuilder::matrix_row_prod);
	return 0;
}

int
ShutdownMatrixRecognizer() {
	delete matrixbuilder::matrix_prod;
	matrixbuilder::matrix_prod = 0;
	delete matrixbuilder::matrix_rows_prod;
	matrixbuilder::matrix_rows_prod = 0;
	delete matrixbuilder::matrix_row_prod;
	matrixbuilder::matrix_row_prod = 0;
	return 0;
}


production *matrixbuilder::matrix_prod = 0;
production *matrixbuilder::matrix_rows_prod = 0;
production *matrixbuilder::matrix_row_prod = 0;


matrixparser::matrixparser(math_recognizer_base *rec_, const nonterminal *nt, const ordered_segments *segs, int caps)
	: matrixbuilder(rec_, 0, 0, nt, segs) {
	VERBOSE(*verb_out << "matrixparser(" << segs->bits << ")\n");
	initmatrixreco(caps);
}

matrixparser::matrixparser() : matrixbuilder() { }

void
matrixparser::read(reader &re, math_recognizer_base *rec_) {
	long caps;
	re.read(caps);
	matrixbuilder::read(re, rec_);
	initmatrixreco(caps);
}

matrixparser::~matrixparser() {
	delete mxiter;
	delete mxan_;
	freeknown();
}

void
matrixparser::initmatrixreco(int caps) {
	capabilities = caps;
	mxan_ = CreateMatrixAnalyzer(ctx(), caps);
	for (std::vector<const segment *>::const_iterator i = span()->begin(TIME_ORDER); i != span()->end(TIME_ORDER); ++i) {
		mxan_->addStroke(&(*i)->stk);
	}
	mxan_->updateRecognition();
	mxiter = mxan_->createIterator();
}
const char matrixparser::ID = 20;

void
matrixparser::write(writer &wr) const {
	wr.write((long)capabilities);
	parser::write(wr);
	//wr.write((size_t)0);
}

void
matrixparser::reset() {
	delete mxiter;
	delete mxan_;
	mxiter = 0;
	mxan_ = 0;
	clearchildren();
	freeknown();
	initmatrixreco(capabilities);
}

}
