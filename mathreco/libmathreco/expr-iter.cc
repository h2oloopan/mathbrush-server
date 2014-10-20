#include "expr-iter.h"
#include "mathrecognizer-private.h"
#include "intrpr.h"
#include "ordered-segments.h"
#include "bitvec.h"
#include "segment.h"
#include "ptab.h"
#include "rect.h"
#include "expr-node.h"
#include "relation.h"
#include "builder.h"
#include "reco-types.h"
#include "verb.h"

namespace scg {

const static double REL_SCORE_THRES = std::log(1e-6);
const static double INTERP_THRES = std::log(0.01);
const static double OP_FACTOR = std::log(0.95);
const static double MAX_SPAN_OVERLAP = 1;

static bool FILTER_SUBTREES = true;
static bool FILTER_RELATIONS = true;

void
DisableSubtreeFiltering() {
	FILTER_SUBTREES = false;
}

void
DisableRelationFiltering() {
	FILTER_RELATIONS = false;
}

std::ostream &
operator<<(std::ostream &os, const recoscore &rhs) {
	return os << rhs.hint << ':' << rhs.score;
}

const ordered_segments *
math_recognizer_base::getsegs(const bitvec &bits) const {
	std::map<bitvec, ordered_segments *>::iterator i = spans->find(bits);
	return i == spans->end() ? 0 : i->second;
}

ordered_segments *
math_recognizer_base::getsegs(const bitvec &bits, const ordered_segments *super, size_t d, size_t start, size_t end, bool ign_overlap) {
	std::map<bitvec, ordered_segments *>::iterator i = spans->find(bits);
	if (i != spans->end()) {
		return i->second;
	}
	ordered_segments *segs = super->slice(d, start, end);
	if (!ign_overlap) {
		bitvec grpbits(segments.size(), false);
		Rect<long> segbox = segs->bounds();
		double totol = 0;
		for (size_t j = 0; j < segs->bits.size(); ++j) {
			if (!segs->bits.at(j)) {
				double ol = overlap_proportion(segments[j]->bounds, segbox);
				if (ol > 0) {
					double rootsc = 0;
					grpbits.set(j, true);
					std::map<bitvec, group *>::const_iterator k = groups->find(grpbits);
					if (k != groups->end()) {
						const group *grp = k->second;
						ol *= (1-grp->container_score);
						/*
						for (std::vector<match_score>::const_iterator mm = grp->final_matches.begin(); mm != grp->final_matches.end(); ++mm) {
							if (mm->proto->info().name == "sqrt") {
								rootsc = mm->matcher_score;
								break;
							}
						}*/
					}
					grpbits.set(j, false);

					totol += ol;
					if (totol >= MAX_SPAN_OVERLAP) break;
				}
			}
		}
		if (totol >= MAX_SPAN_OVERLAP) {
			VERBOSE(*verb_out << "ignoring span " << segs->bits << " because of overlap\n");
			delete segs;
			segs = 0;
		}
		else {
			segs->overlap = totol;
			spans->insert(i, std::make_pair(bits, segs));
		}
	}
	else {
		spans->insert(i, std::make_pair(bits, segs));
	}
	return segs;
}

static bool
contains(const Rect<long> &super, const Rect<long> &sub) {
	return super.left <= sub.left && super.right >= sub.right && super.top <= sub.top && super.bottom >= sub.bottom;
}

static int
attach_in_parent(const interpretation *intrp) {
	if (intrp->P && intrp->P->attach_in_parent != attach_mode::UNSPECIFIED) {
		return intrp->P->attach_in_parent;
	}
	else if (intrp->nchildren() == 1) {
		return attach_in_parent(intrp->child(0));
	}
	return attach_mode::UNSPECIFIED;
}

const char *
modestr(int mode) {
	switch (mode) {
	case attach_mode::GROUP: return "GROUP";
	case attach_mode::SYMBOL: return "SYMBOL";
	case attach_mode::SYMBOL_HEAD: return "SYMBOL-HEAD";
	case attach_mode::SYMBOL_TAIL: return "SYMBOL-TAIL";
	case attach_mode::UNSPECIFIED: return "UNSPECIFIED";
	}
	return "UNKNOWN";
}


static const interpretation *
attaching_intrp(const interpretation *intrp, bool head, int &mode) {//const interpretation *symintrp, int &mode) {
	VERBOSE2(*verb_out << " starting with intrp " << intrp->str() << "; Pmode " << modestr(mode) << std::endl);
	while (mode != attach_mode::GROUP) {
		const interpretation *nextintrp = head ? intrp->head : intrp->tail;
		if (nextintrp == intrp || !nextintrp/* || !nextintrp->P*/) {
			VERBOSE2(*verb_out << " got terminal with next " << (nextintrp ? nextintrp->str() : "") << "; quitting with SYMBOL\n");
			mode = attach_mode::SYMBOL;
			break;
		}
		if (head && mode == attach_mode::SYMBOL_TAIL) {
			mode = attach_mode::GROUP;
			break;
		}
		if (!head && mode == attach_mode::SYMBOL_HEAD) {
			mode = attach_mode::GROUP;
			break;
		}
		intrp = nextintrp;
		if (intrp->P) {
			mode = intrp->P->attach_in_parent;
		}
		VERBOSE2(*verb_out << " got intrp " << intrp->str() << "; Pmode " << modestr(mode) << std::endl);
	}
	VERBOSE2(*verb_out << " RETURNING with intrp " << intrp->str() << " and mode " << modestr(mode) << std::endl);
	return intrp;
	/*switch (Pmode) {
	case attach_mode::GROUP:
		mode = Pmode;
		return intrp;
	case attach_mode::SYMBOL:
		mode = Pmode;
		return symintrp;
	/*case attach_mode::UNSPECIFIED:
		if (!intrp->P) {
			mode = attach_mode::SYMBOL;
			return intrp;
		}
		int parentmode = attach_in_parent(intrp);
		switch (parentmode) {
		case attach_mode::GROUP:
		case attach_mode::UNSPECIFIED:
			mode = attach_mode::GROUP;
			return intrp;
		case attach_mode::SYMBOL:
			mode = attach_mode::SYMBOL;
			return symintrp;
		}* /
	default:
		THROW_ERROR(E_INVALID, "invalid Pmode combination");
	}*/
}

struct relscore {
	double score;
	int tailrc;
	int headrc;
	relscore() : score(-1.0), tailrc(0), headrc(0) { }
	relscore(double score_, int tailrc_, int headrc_) : score(score_), tailrc(tailrc_), headrc(headrc_) { }
};


static double
rscore(const production *P, const ordered_segments *tail, const ordered_segments *head) {
	interpretation taildummy(0, P, tail);
	interpretation headdummy(0, P, head);
	VERBOSE(*verb_out << "prelim rscore between " << tail->bits << " and " << head->bits << std::endl);
	size_t i = relindex(P->rel);
	relation *rel = get_group_relation(i);
	return std::log(rel->membership(&taildummy, &headdummy, AGGREGATE_CLASS, AGGREGATE_CLASS, attach_mode::GROUP, attach_mode::GROUP));
	/*
	double sc = 0;
	for (int j = NCLASSES-1; j >= 0; --j) {
		for (int k = NCLASSES-1; k >= 0; --k) {
			double clsc = P->rel->membership(tail, head, j, k, attach_mode::DEFAULT, attach_mode::DEFAULT);
			if (clsc > 0) return clsc;
		}
	}
	return 0;*/
}

const symbol *
getsymbol(const interpretation *intrp) {
	if (!intrp->P) {
		return symdb_findsymbol_name(intrp->str());
	}
	if (intrp->nchildren() != 1) {
		return 0;
	}
	return getsymbol(intrp->child(0));
}

static relscore
rscore_internal(const relation *rel, const attach_mode &att, const interpretation *tail, const interpretation *head) {
	int tailmode = att.from;
	int headmode = att.to;
	int tailclass, headclass;
	const interpretation *tailintrp = tail;
	const interpretation *headintrp = head;
	VERBOSE2(*verb_out << "rscore_internal for " << rel->name << " between " << tail->span->bits << " and " << head->span->bits << std::endl);
	VERBOSE2(*verb_out << "tail is " << tail->str() << " ; real tail " << tail->tail->str() << " ; rel tail " << tailintrp->str() << std::endl);
	VERBOSE2(*verb_out << "head is " << head->str() << " ; real head " << head->head->str() << " ; rel head " << headintrp->str() << std::endl);
	VERBOSE2(
		*verb_out << "rscore between " << tailintrp->str() << " and " << headintrp->str() << " with classes ";
		for (int i = 0; i < NCLASSES; ++i) {
			if (tailintrp->rclass & (1 << i)) {
				*verb_out << rclass_to_str(i) << ' ';
			}
		}
		*verb_out << " and ";
		for (int i = 0; i < NCLASSES; ++i) {
			if (headintrp->rclass & (1 << i)) {
				*verb_out << rclass_to_str(i) << ' ';
			}
		}
	);
	if (tailmode == attach_mode::GROUP && (tailintrp->rclass & (1 << BOX_CLASS)) && (!tailintrp->P || tailintrp->P->rclass == 0)) {
		tailclass = (1 << BOX_CLASS) | (1 << AGGREGATE_CLASS);
	}
	else {
		tailclass = tailintrp->rclass;
	}
	if (headmode == attach_mode::GROUP && (headintrp->rclass & (1 << BOX_CLASS)) && (!headintrp->P || headintrp->P->rclass == 0)) {
		headclass = (1 << BOX_CLASS) | (1 << AGGREGATE_CLASS);
	}
	else {
		headclass = headintrp->rclass;
	}


	VERBOSE2(
		*verb_out << "rscore between " << tail->str() << " and " << head->str() << " with classes ";
		for (int i = 0; i < NCLASSES; ++i) {
			if (tail->rclass & (1 << i)) {
				*verb_out << rclass_to_str(i) << ' ';
			}
		}
		*verb_out << " and ";
		for (int i = 0; i < NCLASSES; ++i) {
			if (head->rclass & (1 << i)) {
				*verb_out << rclass_to_str(i) << ' ';
			}
		}
		*verb_out << std::endl;
		*verb_out << "testing interps " << tailintrp->str() << " (" << tailintrp->score() << ") and " << headintrp->str() << " (" << headintrp->score() << ") with classes ";
		for (int i = 0; i < NCLASSES; ++i) {
			if (tailintrp->rclass & (1 << i)) {
				*verb_out << rclass_to_str(i) << ' ';
			}
		}
		*verb_out << " and ";
		for (int i = 0; i < NCLASSES; ++i) {
			if (headintrp->rclass & (1 << i)) {
				*verb_out << rclass_to_str(i) << ' ';
			}
		}
		*verb_out << std::endl;
		*verb_out << "connection modes " << modestr(tailmode) << " and " << modestr(headmode) << std::endl;
	);

	static std::vector<rclass_t> tailclasses;
	static std::vector<rclass_t> headclasses;
	tailclasses.reserve(32); headclasses.reserve(32);
	tailclasses.clear(); headclasses.clear();

	const symbol *S;
	//if (S = getsymbol(tailintrp)) tailclasses.push_back(S->sid);
	//if (S = getsymbol(headintrp)) headclasses.push_back(S->sid);
	
	for (int j = NCLASSES-1; j >= 0; --j) {
		if (tailclass & (1 << j)) {
			tailclasses.push_back(j);
		}
		if (headclass & (1 << j)) {
			headclasses.push_back(j);
		}
	}

	assert(tailclasses.size() >= 1);
	assert(headclasses.size() >= 1);

	VERBOSE(
		*verb_out << "rscore_internal from " << tailintrp->span->bounds() << " to " << headintrp->span->bounds() << " with classes [ ";
		for (std::vector<rclass_t>::const_iterator i = tailclasses.begin(); i != tailclasses.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "] and [ ";
		for (std::vector<rclass_t>::const_iterator i = headclasses.begin(); i != headclasses.end(); ++i) {
			*verb_out << *i << ' ';
		}
		*verb_out << "]\n";
	);
	relscore rs;
	for (size_t i = 0; i <= tailclasses.size() + headclasses.size() - 2; ++i) {
		double score = 0.0;
		unsigned n = 0;
		rclass_t meastail = 0, meashead = 0;
		for (size_t j = 0; j <= std::min(i, tailclasses.size() - 1); ++j) {
			rclass_t tailclass = tailclasses[j];
			rclass_t headclass = headclasses[i-j];
			double clsc = rel->membership(tailintrp, headintrp, tailclass, headclass, tailmode, headmode);
			//if (rel) clsc = rel->membership(tailintrp, headintrp, tailclass, headclass, tailmode, headmode);
			//else clsc = get_nil_probability(tailintrp, headintrp, tailclass, headclass);
			//rs.score = std::max(rs.score, clsc);
			if (clsc != -1.0) {
				score += clsc;
				meastail |= 1 << tailclass;
				meashead |= 1 << headclass;
				++n;
			}
		}
		if (n > 0) {
			rs.score = score / n;
			rs.tailrc = meastail;
			rs.headrc = meashead;
			break;
		}
	}

foundscore:
	rs.headrc = headclass;
	rs.tailrc = tailclass;
	//rs.score = std::log(rs.score);
	//VERBOSE(*verb_out << "reported log score is " << rs.score << std::endl);
	return rs;
}

static relscore
rscore(const relation *rel, const attach_mode &att, interpretation *tail, interpretation *head) {
	VERBOSE(*verb_out << "rscore for " << rel->name << " from " << tail->span->bits << " to " << head->span->bits << std::endl);
	//if (rel->name == "R") {// && tail->span->bits[0] == 49 && head->span->bits[0] == 6) {
		//int a = 0;
	//}

	int tailmode = att.from;
	int headmode = att.to;
	const interpretation *tailintrp = attaching_intrp(tail, false, tailmode);//tail->tail, tailmode);
	const interpretation *headintrp = attaching_intrp(head, true, headmode);//head->head, headmode);

	attach_mode real_attach(tailmode, headmode);

	relscore rs;
	double sum = 0.0;
	for (size_t i = 0; i < NUM_RELATIONS; ++i) {
		const relation *R = get_relation(i);
		relscore subrs = rscore_internal(R, real_attach, tailintrp, headintrp);
		VERBOSE(*verb_out << " subscore for rel " << R->name << " = " << subrs.score << std::endl);
		if (R == rel) {
			rs = subrs;
		}
		sum += subrs.score;
	}
	
	double pnil = get_nil_probability(tailintrp, headintrp);
	//relscore pnil_score = rscore_internal(0, real_attach, tailintrp, headintrp);
	//double pnil = pnil_score.score;
	//double pnil = 1.0 - sum / (sum + 1.0);
	VERBOSE(*verb_out << " final score = " << rs.score << " / " << sum << " = " << rs.score / sum << "; pnil " << pnil << " = " << rs.score / pnil);
	//rs.score = (1.0 / pnil - 1.0) * rs.score / sum;
	rs.score = (1.0 - pnil) * rs.score / sum;
	rs.score = std::log(rs.score / pnil);
	//rs.score = std::log(rs.score / sum);
	VERBOSE(*verb_out << " -> log " << rs.score << std::endl);
	return rs;
}

/*
struct bldata {
	std::vector<long> baselines;
	std::vector<long> heights;
	std::vector<int> rclasses;
	std::vector<double> scores;
};

/*
static double blscore(const interpretation *intrp, bldata &dat);
static void
gatherbldata(const interpretation *intrp, bldata &dat) {
	for (size_t i = 0; i < intrp->nchildren(); ++i) {
		const interpretation *child = intrp->child(i);
		if (derivesterminal(child)) {
		}
		else if (child->P->rel == get_relation(REL_RIGHT)) {
			gatherbldata(intrp, dat);
		}
		else {
			double childsc = blscore(child, dat);
			dat.scores.push_back(childsc);
		}
	}
}* /

static double
blscore(const interpretation *intrp, bldata &dat) {
	VERBOSE(*verb_out << "blscore(" << intrp->str() << ")\n");
	if (derivesterminal(intrp)) {
		for (int i = NCLASSES-1; i >= 0; --i) {
			if (intrp->rclass & (1 << i)) {
				dat.baselines.push_back(getbl(intrp->bounds(), i));
				dat.heights.push_back(std::max(intrp->bounds().height(), (long)(0.125*TABLETPC_DPI)));
				dat.rclasses.push_back(i);
				VERBOSE(*verb_out << " gathered bl baseline " << dat.baselines.back() << " and height " << dat.heights.back() << " with rclass " << rclass_to_str(i) << std::endl);
				break;
			}
		}
		return 1;
	}
	else if (intrp->P->rel != get_relation(REL_RIGHT)) {
		double sc = 1;
		for (size_t i = 0; i < intrp->nchildren(); ++i) {
			bldata d;
			sc *= blscore(intrp->child(i), d);
		}
		dat.scores.push_back(std::pow(sc, 1.0/intrp->nchildren()));
		return 1;
	}
	else {
		for (size_t i = 0; i < intrp->nchildren(); ++i) {
			blscore(intrp->child(i), dat);
		}
		VERBOSE(*verb_out << "finalizing blscore for " << intrp->str() << std::endl);
		double blsc = 1;
		double hsc = 1;
		double sc = 1;
		unsigned n = 0;
		for (size_t i = 0; i < dat.baselines.size(); ++i) {
			for (size_t j = i+1; j < dat.baselines.size(); ++j) {
				//long hgt = std::max(dat.heights[i], dat.heights[j]);
				//double db = (double)(dat.baselines[i]-dat.baselines[j])/hgt;
				//double free = std::max(0.25*dat.heights[j], 0.125*TABLETPC_DPI);
				//double limit = std::max(1.5*dat.heights[j], 0.125*TABLETPC_DPI);
				//double diff = std::min((double)limit+free, std::max(0.0, (double)dat.baselines[j] - dat.baselines[i] + free));
				//double bsc = std::pow(diff/(limit+free), 1.0/4);
				//double bsc = std::exp(-2*std::log(2.0) * std::abs(db));
				//VERBOSE(*verb_out << "  between items " << i << " and " << j << ", baseline diff " << db << " gives score " << bsc << std::endl);

				double dhreal = (double)dat.heights[i]/dat.heights[j];
				double dhnominal = (double)nominal_height(dat.rclasses[i])/nominal_height(dat.rclasses[j]);
				double r = std::max(dhreal/dhnominal, dhnominal/dhreal);
				double hsc = 1/std::pow(r, 1.0/4);
				double ns = hsc;//std::sqrt(bsc*hsc);
				sc *= ns;
				VERBOSE(*verb_out << "  and height ratio " << r << " gives score " << hsc << " for " << ns << "final\n");
				++n;
			}
		}
		sc = std::pow(sc, 1.0/n);
		VERBOSE(*verb_out << " local score normalized to " << sc << " and global score is ");
		dat.scores.push_back(sc);
		sc = 1;
		for (std::vector<double>::const_iterator i = dat.scores.begin(); i != dat.scores.end(); ++i) {
			sc *= *i;
		}
		sc = std::pow(sc, 1.0/dat.scores.size());
		VERBOSE(*verb_out << sc << std::endl);
		return sc;
	}
}*/

treeparser::treeparser(math_recognizer_base *rec, const production *P_, const ordered_segments *segs_, std::vector<interpreter *> &parsers)
	: parser(rec, P_->nt, segs_, 1), P(P_) {
	VERBOSE(*verb_out << "treeparser(" << P->nt->name << ":" << segs_->bits << " =" << (P->rel ? P->rel->name : std::string("")) << "=";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			*verb_out << " " << P->rhs[i]->name;
		}
		*verb_out << std::endl;
	);
	assert(parsers.size() == P->rhs.size());
	for (std::vector<interpreter *>::iterator i = parsers.begin(); i != parsers.end(); ++i) {
		addchild(*i, false);
	}
	reset();
}

treeparser::~treeparser() {
	freeknown();
}

bool
treeparser::addtree(const std::vector<size_t> &ranks) {
	VERBOSE(
		*verb_out << "addtree(" << P->nt->name << ":" << span()->bits << " =" << (P->rel ? P->rel->name : std::string("")) << "=";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			*verb_out << " " << P->rhs[i]->name;
		}
		*verb_out << ", [";
		for (std::vector<size_t>::const_iterator i = ranks.begin(); i != ranks.end(); ++i) {
			*verb_out << *i;
			if (i + 1 != ranks.end()) {
				*verb_out << ',';
			}
		}
		*verb_out << "])\n";
	);
	interpretation *intrp = new interpretation(this, P, span());
	for (size_t i = 0; i < nchildren(); ++i) {
		/*
		if (!child(i)->isactive()) {
			delete intrp;
			intrp = 0;
			break;
		}*/
		interpretation *chintrp = getnth(child(i), ranks[i]);
		if (!chintrp) {
			delete intrp;
			intrp = 0;
			break;
		}
		VERBOSE(*verb_out << " treeparser child " << i << " is " << chintrp << " or " << chintrp->str() << std::endl);
		if (i == P->head) {
			intrp->head = chintrp;//->head;
		}
		if (i == P->tail) {
			intrp->tail = chintrp;//->tail;
		}
		double rsc = 0;
		if (i > 0) {
			relscore rs = rscore(P->rel, P->attach_modes[i-1], intrp->child(i-1), chintrp);
			if (rs.score <= REL_SCORE_THRES) {
				VERBOSE(*verb_out << "CUTTING INTRP (REL " << rs.score << ")!\n");
				if (intrp->child(i-1)->score().hint == 0 && chintrp->score().hint == 0) {
					//VERBOSE(*verb_out << "  failed rel-score thres with hints " << intrp->child(i-1)->str() << " ==> " << intrp->child(i-1)->hint << " and " << chintrp->str() << " ==> " << chintrp->hint << std::endl);
					delete intrp;
					intrp = 0;
					break;
				}
			}
			VERBOSE(*verb_out << " GOT RSCORE " << rs.score << std::endl);
			rsc = rs.score;
			if (P->rel == get_relation(REL_RIGHT)) {
				// XXX: note rclass may vary between neighbouring rels. that is weird.
				intrp->rclass |= rs.tailrc | rs.headrc | (1 << RELCLASS_MERGE);
			}
			intrp->rclass |= (1 << BOX_CLASS);
			//intrp->rclass &= ~(1 << SYMBOL_CLASS);
			//intrp->set_score(intrp->score_combo().addrel(rs.score));
		}
		intrp->addchild(chintrp, rsc);
	}

	if (intrp) {
		//typetree tt;
		//mktypetree(tt, intrp->mktree(ctx()));

		
		/*ExpressionTree *tree = intrp->mktree();
		double prob = ptab_measuretop(get_ptab(), tree);
		VERBOSE(*verb_out << "got probability " << prob << " from ptab for tree " << tree->str() << std::endl);
		intrp->scorebias += std::log(prob);
		tree->release();*/

		if (intrp->score() < (2*intrp->nterms-1) * INTERP_THRES) {
			VERBOSE(*verb_out << "Pruned expr with score " << intrp->score() << " vs. thres " << (2*intrp->nterms-1) * INTERP_THRES << std::endl);
			delete intrp;
			known_ranks.insert(ranks);
			return false;
		}
		if (intrp->nchildren() > 1) {
			/*bldata dat;
			double blsc = blscore(intrp, dat);
			blsc = 0.75 + 0.25 * std::exp(-4.0/3*std::log(2)*(1-blsc));
			intrp->set_score(intrp->score_combo().bias(blsc).addop());*/
			//intrp->set_score(intrp->score_combo().bias(intrp->P->bias));
			if (intrp->P->sid != InvalidSemanticId) {// && !intrp->P->nt->ispartial()) {
				//intrp->set_score(intrp->score_combo().addop());
				intrp->scorebias += OP_FACTOR;
			}
		}

		intrp->rclass |= intrp->P->rclass;

		if (intrp->score() != intrp->score()) {
			recoscore d = intrp->score();
		}
		VERBOSE(
			*verb_out << "to " <<P->nt->name << ":" << span()->bits << " =" << (P->rel ? P->rel->name : std::string("")) << "=";
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				*verb_out << " " << P->rhs[i]->name;
			}
			*verb_out << " adding known parse " << intrp->str() << " of type " << P->sid << " with score " << intrp->score() << " making " << known.size() + 1 << std::endl;
		);
		if (!P->rel) {
			intrp->rclass = intrp->child(0)->rclass;
		}
		/*else if (P->rel != get_relation(REL_RIGHT)) {
			intrp->rclass = (1 << AGGREGATE_CLASS) | (1 << BOX_CLASS);
		}*/
		VERBOSE(
			*verb_out << " expr has classes ";
			for (int i = 0; i < NCLASSES-1; ++i) {
				if (intrp->rclass & (1 << i)) {
					*verb_out << rclass_to_str(i) << ' ';
				}
			}
			*verb_out << std::endl;
		);
		//intrp->hint += hint;
		known.insert(treerec(intrp, ranks));
		return true;
	}
	known_ranks.insert(ranks);
	return false;
}

interpretation *
treeparser::getnext() {
	VERBOSE(
		*verb_out << "treeparser:" << this << "->getnext(" << P->nt->name << ":" << span()->bits << " =" << (P->rel ? P->rel->name : std::string("")) << "=";
		for (size_t i = 0; i < P->rhs.size(); ++i) {
			*verb_out << " " << P->rhs[i]->name;
		}
		*verb_out << ")\n";
	);
	if (!FILTER_SUBTREES || nchildren() == 1) {
		if (!last.ranks.empty() && ((P->sid != InvalidSemanticId && !FILTER_SUBTREES) || derivesterminal(last.node) || P->sid == InvalidSemanticId)) {
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				++last.ranks[i];
				if (known_ranks.find(last.ranks) == known_ranks.end()) {
					addtree(last.ranks);
				}
				--last.ranks[i];
			}
		}
	}

	if (known.empty()) {
		VERBOSE(*verb_out << " (none)\n");
		return 0;		
	}

	last = *known.begin();
	known.erase(known.begin());
	if (FILTER_SUBTREES && nchildren() > 1) {
		freeknown();
	}
	VERBOSE(*verb_out << " " << last.node->str() << std::endl);
	return last.node;
}

void
treeparser::freeknown() {
	parser::freeknown();
	while (!known.empty()) {
		const treerec &best = *known.begin();
		delete best.node;
		known.erase(known.begin());
	}
}

static bool
ismatrix(const production *P) {
	return P->nt->name == "MATRIX" || P->nt->name == "MATRIX_" || P->nt->name == "MX";
}

void
treeparser::reset() {
	freeknown();
	last.ranks.clear();
	last.node = 0;
	known_ranks.clear();
	std::vector<size_t> ranks(P->rhs.size(), 0);
	//addtree(ranks);
	
	if (ismatrix(P)) {
		addtree(ranks);
	}
	else {
		if (!addtree(ranks)) {
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				++ranks[i];
				if (known_ranks.find(ranks) == known_ranks.end()) {
					addtree(ranks);
				}
				--ranks[i];
			}
		}
		
		if (!known.empty()) {
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				interpretation *intrp = getfirst(child(i));
				assert(intrp);
				if (derivesterminal(intrp)) {
					++ranks[i];
					addtree(ranks);
					--ranks[i];
				}
			}
		}
	}
}

treeparser::treeparser() : parser(), P(0) { }

void
treeparser::read(reader &re, math_recognizer_base *rec) {
	parser::read(re, rec);
	size_t Pindex;
	re.read(Pindex);
	if (Pindex >= nt()->productions.size()) {
		THROW_ERROR(E_INVALID, "while reading treeparser, production index " << Pindex << " was too large for nonterminal " << nt()->name);
	}
	P = nt()->productions[Pindex];
	/*size_t n;
	re.read(n);
	while (n--) {
		treerec tr;
		read_treerec(re, tr, rec);
		known.insert(tr);
	}
	read_treerec(re, last, rec);
	re.read(n);
	while (n--) {
		std::vector<size_t> ranks;
		read_ranks(re, ranks);
		known_ranks.insert(ranks);
	}*/
	reset();
}

void
treeparser::write(writer &wr) const {
	parser::write(wr);
	wr.write(P->index);
	/*
	wr.write(known.size());
	for (known_set::const_iterator i = known.begin(); i != known.end(); ++i) {
		write_treerec(wr, *i);
	}
	write_treerec(wr, last);
	wr.write(known_ranks.size());
	for (std::set<std::vector<size_t> >::const_iterator i = known_ranks.begin(); i != known_ranks.end(); ++i) {
		const std::vector<size_t> &r = *i;
		write_ranks(wr, r);
	}*/
}

const char treeparser::ID = 10;


static size_t
minlen(const production *P, size_t start, size_t len) {
	assert(start + len <= P->rhs.size());
	if (len == 0) {
		return 0;
	}
	if (start == 0) {
		return P->min_prefix_length[len - 1];
	}
	else {
		return P->min_prefix_length[start + len - 1] - P->min_prefix_length[start - 1];
	}
}

static size_t
maxlen(const production *P, size_t start, size_t len) {
	assert(start + len <= P->rhs.size());
	if (len == 0) {
		return 0;
	}
	if (start == 0) {
		return P->max_prefix_length[len - 1];
	}
	else if (P->max_prefix_length[start + len - 1] == std::numeric_limits<size_t>::max()) {
		if (P->max_prefix_length[start - 1] == std::numeric_limits<size_t>::max()) {
			size_t tlen = 0;
			for (size_t i = start; i < start + len; ++i) {
				if (P->rhs[i]->max_length == std::numeric_limits<size_t>::max()) {
					return std::numeric_limits<size_t>::max();
				}
				tlen += P->rhs[i]->max_length;
			}
			return tlen;
		}
		else {
			return std::numeric_limits<size_t>::max();
		}
	}
	else {
		return P->max_prefix_length[start + len - 1] - P->max_prefix_length[start - 1];
	}
}

static size_t
sumormax(size_t a, size_t b) {
	if (a == std::numeric_limits<size_t>::max() || b == std::numeric_limits<size_t>::max()) {
		return std::numeric_limits<size_t>::max();
	}
	else {
		return a + b;
	}
}

struct prodlenspec {
	size_t minprelen;
	size_t maxprelen;
	size_t minxlen;
	size_t maxxlen;
	size_t minpostlen;
	size_t maxpostlen;

	prodlenspec(const production *P, size_t start, size_t xstart, size_t xlen, size_t end,
	            const ordered_segments *segs, size_t startseg) {
		size_t left = end - (xstart + xlen);
		size_t minprelen_fromprod = minlen(P, start, xstart - start);
		size_t minpostlen_fromprod = minlen(P, xstart + xlen, left);
		size_t maxprelen_fromprod = maxlen(P, start, xstart - start);
		size_t maxpostlen_fromprod = maxlen(P, xstart + xlen, left);
		size_t minxlen_fromprod = minlen(P, xstart, xlen);
		size_t maxxlen_fromprod = maxlen(P, xstart, xlen);

		assert(segs->size() > startseg);
		size_t nsegs = segs->size() - startseg;
		size_t s = sumormax(maxxlen_fromprod, maxpostlen_fromprod);
		minprelen = std::max(minprelen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(minxlen_fromprod, minpostlen_fromprod);
		maxprelen = std::min(maxprelen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(maxprelen_fromprod, maxpostlen_fromprod);
		minxlen = std::max(minxlen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(minprelen_fromprod, minpostlen_fromprod);
		maxxlen = std::min(maxxlen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(maxprelen_fromprod, maxxlen_fromprod);
		minpostlen = std::max(minpostlen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));
		s = sumormax(minprelen_fromprod, minxlen_fromprod);
		maxpostlen = std::min(maxpostlen_fromprod, ((s >= nsegs) ? 0 : nsegs - s));

		VERBOSE(
			*verb_out << "prodlenspec(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				*verb_out << " " << P->rhs[i]->name;
			}
			*verb_out << ", start " << start << ", xstart " << xstart << ", xlen " << xlen << ", end " << end << ", bits " << segs->bits << ", startseg " << startseg << ")\n";
			*verb_out << " pre=(" << minprelen << "," << maxprelen << ")  ";
			*verb_out << " x=(" << minxlen << "," << maxxlen << ")  ";
			*verb_out << " post=(" << minpostlen << "," << maxpostlen << ")\n";
		);
	}

	inline bool valid() const {
		return minprelen <= maxprelen && minxlen <= maxxlen && minpostlen <= maxpostlen;
	}
};


class prodintrprbuilder {
public:
	bool ign_overlap;
	prodintrprbuilder(math_recognizer_base *rec_, const production *P_, const ordered_segments *segs_, bool ign_overlap_)
		: rec(rec_), P(P_), segs(segs_), intrpr(0), mux(0),
		  termstarts(P_->terminal_indices.size()), childsegs(P_->rhs.size()), parsers(P_->rhs.size()),
			ign_overlap(ign_overlap_) {
		ign_overlap = ign_overlap || (P->rel && P->rel == get_relation(REL_CONTAINS));
		assert(segs->size() >= P->min_prefix_length[P->rhs.size()-1]);
		VERBOSE(
			*verb_out << "prodintrprbuilder(" << P->nt->name << ":" << segs->bits << " =" << (P->rel ? P->rel->name : std::string("")) << "=";
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				*verb_out << " " << P->rhs[i]->name;
			}
			*verb_out << ")\n";
		);
		parseterms(0, 0);
	}

private:
	void parsents(size_t cterm, size_t cseg, size_t crhsi, const ordered_segments *subspan, size_t csubseg) {
		VERBOSE(
			*verb_out << "parsents(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				*verb_out << " " << P->rhs[i]->name;
			}
			*verb_out << ", bits " << segs->bits << ", term " << cterm << ", seg " << cseg << ", rhs " << crhsi << ", subspan ";
			if (subspan) {
				*verb_out << subspan->bits;
			}
			else {
				*verb_out << "nil";
			}
			*verb_out << ", subseg " << csubseg << std::endl;
		);
		if (crhsi == P->rhs.size()) {
			if (!intrpr) {
				intrpr = new treeparser(rec, P, segs, parsers);
			}
			else {
				if (!mux) {
					mux = new multiplexor(rec, P->nt, segs, 2);
					mux->addparser(intrpr);
					intrpr = mux;
				}
				mux->addparser(new treeparser(rec, P, segs, parsers));
			}
		}
		else if (cterm < P->terminal_indices.size() && crhsi == P->terminal_indices[cterm]) {
			if (FILTER_RELATIONS && P->rel && crhsi > 0) {
				double rsc = rscore(P, childsegs[crhsi - 1], childsegs[crhsi]);
				if (rsc <= REL_SCORE_THRES) {
					return;
				}
			}
			size_t d = P->rel ? ordering_for_relation(P->rel) : 0;
			size_t nextbreak = (cterm + 1 == P->terminal_indices.size()) ? segs->size() : termstarts[cterm + 1];
			cseg += childsegs[P->terminal_indices[cterm]]->size();
			const ordered_segments *nextspan;
			if (cseg == nextbreak) {
				nextspan = 0;
			}
			else {
				nextspan = rec->getsegs(segs->slice_bits(d, cseg, nextbreak), segs, d, cseg, nextbreak, ign_overlap);
				if (!nextspan) return;
			}
			return parsents(cterm + 1, cseg, crhsi + 1, nextspan, 0);
		}
		else {
			prodlenspec lenspec(P, crhsi, crhsi, 1, cterm == P->terminal_indices.size() ? P->rhs.size() : P->terminal_indices[cterm], subspan, csubseg);
			if (!lenspec.valid()) {
				return;
			}
			size_t d = P->rel ? ordering_for_relation(P->rel) : 0;
			size_t minlen = lenspec.minxlen;
			size_t sz = subspan->size();
			if (lenspec.maxpostlen <= sz && lenspec.maxpostlen + csubseg <= sz) {
				minlen = std::max(minlen, sz - csubseg - lenspec.maxpostlen);
			}
			size_t maxlen = std::min(sz - lenspec.minpostlen - csubseg, lenspec.maxxlen);
			for (size_t len = minlen; len <= maxlen; ++len) {
				VERBOSE(*verb_out << " confnt trying length " << len << std::endl);
				const ordered_segments *subsegs = rec->getsegs(subspan->slice_bits(d, csubseg, csubseg + len), subspan, d, csubseg, csubseg + len, ign_overlap);
				if (!subsegs) continue;
				VERBOSE(*verb_out << "  got segs " << subsegs->bits << std::endl);
				if (FILTER_RELATIONS && P->rel && crhsi > 0) {
					double rsc = rscore(P, childsegs[crhsi - 1], subsegs);
					if (rsc <= REL_SCORE_THRES) {
						continue;
					}
				}
				parsers[crhsi] = rec->mkparser(P->rhs[crhsi], subsegs, ign_overlap);
				if (parsers[crhsi]) {
					// TODO: insert pre-check for specific relation score > 0 ? 
					childsegs[crhsi] = subsegs;
					parsents(cterm, cseg + len, crhsi + 1, subspan, csubseg + len);
				}
			}
		}
	}

	void parseterms(size_t cterm, size_t cseg) {
		VERBOSE(
			*verb_out << "parseterms(" << P->nt->name << " :" << (P->rel ? P->rel->name : std::string("")) << ":";
			for (size_t i = 0; i < P->rhs.size(); ++i) {
				*verb_out << " " << P->rhs[i]->name;
			}
			*verb_out << ", bits " << segs->bits << ", term " << cterm << ", seg " << cseg << std::endl;
		);

		if (cterm == P->terminal_indices.size()) {
			size_t d = P->rel ? ordering_for_relation(P->rel) : 0;
			const ordered_segments *subspan;
			if (P->terminal_indices.empty()) {
				subspan = segs;
			}
			else {
				if (termstarts[0] == 0) {
					subspan = 0;
				}
				else {
					subspan = rec->getsegs(segs->slice_bits(d, 0, termstarts[0]), segs, d, 0, termstarts[0], ign_overlap);
					if (!subspan) return;
				}
			}
			parsents(0, 0, 0, subspan, 0);
		}
		else {
			size_t last_tpos;
			if (cterm == 0) {
				last_tpos = 0;
			}
			else {
				last_tpos = P->terminal_indices[cterm-1]+1;
			}
			size_t tpos = P->terminal_indices[cterm];
			prodlenspec lenspec(P, last_tpos, tpos, 1, P->rhs.size(), segs, cseg);
			if (!lenspec.valid()) {
				return;
			}

			size_t d = P->rel ? ordering_for_relation(P->rel) : 0;
			const nonterminal *termnt = P->rhs[tpos];

			for (size_t start = cseg + lenspec.minprelen; start <= cseg + lenspec.maxprelen; ++start) {
				assert(segs->size() > lenspec.minpostlen + start);
				size_t minlen = lenspec.minxlen;
				if (lenspec.maxpostlen <= segs->size() && lenspec.maxpostlen + start <= segs->size()) {
					minlen = std::max(minlen, segs->size() - start - lenspec.maxpostlen);
				}
				size_t maxlen = std::min(segs->size() - lenspec.minpostlen - start, lenspec.maxxlen);
				for (size_t len = minlen; len <= maxlen; ++len) {
					VERBOSE(*verb_out << " conf trying start " << start << ", len " << len << std::endl);
					bool ignol = ign_overlap;
					if (P && P->rel == get_relation(REL_BELOW) && termnt->name == "_Shorzline") {
						ignol = true;
					}
					const ordered_segments *subsegs = rec->getsegs(segs->slice_bits(d, start, start + len), segs, d, start, start + len, ignol);
					if (!subsegs) continue;
					if (P->rel && tpos != 0) {
						const ordered_segments *presegs;
						if (last_tpos == tpos) {
							presegs = childsegs[last_tpos-1];
						}
						else {
							presegs = rec->getsegs(segs->slice_bits(d, cseg, start), segs, d, cseg, start, ign_overlap);
							if (!presegs) continue;
						}
						if (FILTER_RELATIONS) {
							double rsc = rscore(P, presegs, subsegs);
							if (rsc <= REL_SCORE_THRES) {
								continue;
							}
						}
					}
					interpreter *tref = rec->mkparser(termnt, subsegs, ign_overlap);
					if (tref) {
						termstarts[cterm] = start;
						childsegs[tpos] = subsegs;
						parsers[tpos] = tref;
						parseterms(cterm + 1, start + len);
					}
				}
			}
		}
	}


public:
	interpreter *intrpr;
private:
	math_recognizer_base *rec;
	const production *P;
	const ordered_segments *segs;
	multiplexor *mux;
	std::vector<const ordered_segments *> childsegs;
	std::vector<size_t> termstarts;
	std::vector<interpreter *> parsers;
};

static interpreter *
mkprodintrpr(math_recognizer_base *rec, production_parser_map &prod_tab, const production *P, const ordered_segments *segs, bool ign_overlap) {
	interpreter *&intrpr = prod_tab[P];
	if (!intrpr) {
		intrpr = prodintrprbuilder(rec, P, segs, ign_overlap).intrpr;
	}
	return intrpr;
}

interpreter *
math_recognizer_base::mkprodintrpr(const production *P, const ordered_segments *segs, bool ign_overlap) {
	if (segs->size() < P->min_prefix_length[P->rhs.size()-1] || segs->size() > P->max_prefix_length[P->rhs.size()-1]) {
		return 0;
	}
	assert((*spans)[segs->bits] == segs);
	nt_parse_table &nttab = tab[segs];
	production_parse_table &prodtab = nttab[P->nt];
	return scg::mkprodintrpr(this, prodtab.tab, P, segs, ign_overlap);
}

static interpreter *
build_ntintrpr(math_recognizer_base *rec, production_parse_table &nttab, const nonterminal *nt, const ordered_segments *segs, bool ign_overlap) {
	VERBOSE(*verb_out << "build_ntintrpr(" << nt->name << ":" << segs->bits << ")\n");
	multiplexor *mux = 0;
	interpreter *intrpr = 0;
	//if (nt->productions.size() > 1) {
		mux = new multiplexor(rec, nt, segs, 3);
		intrpr = mux;
	//}
	for (std::vector<production *>::const_iterator i = nt->productions.begin(); i != nt->productions.end(); ++i) {
		const production *P = *i;
		VERBOSE(*verb_out << " " << segs->size() << " segs with minlen " << P->min_prefix_length[P->rhs.size() - 1] << std::endl);
		if (segs->size() >= P->min_prefix_length[P->rhs.size() - 1]) {
			interpreter *locintrpr;
			int priority = 0;
			locintrpr = mkprodintrpr(rec, nttab.tab, P, segs, ign_overlap);
			if (P->rhs.size() == 1) {
				const std::string &name = P->rhs[0]->name;
				if (name[name.length()-1] == '_') {
					priority = -1;
				}
			}
			/*if (!intrpr) {
				intrpr = locintrpr;
			}
			else*/ if (locintrpr) {
				mux->addparser(locintrpr, false, priority);
			}
		}
	}
	return intrpr;
}

interpreter *
math_recognizer_base::mktermparser(const nonterminal *nt, const ordered_segments *segs) {
	std::map<bitvec, group *>::iterator pgrp = groups->find(segs->bits);
	if (pgrp != groups->end()) {
		group *grp = pgrp->second;
		std::map<const nonterminal *, match_score>::const_iterator i = grp->matches.find(nt);
		if (i != grp->matches.end()) {
			const match_score &match = i->second;
			VERBOSE(*verb_out << "found terminal " << nt->name << " at " << segs->bits << std::endl);
			staticintrpr *tintrpr = new staticintrpr(this, nt, segs);
			interpretation *tnode = mkterminal(tintrpr, match.S, grp, segs);
			tintrpr->addknown(tnode, true);
			if (tnode->score() < INTERP_THRES) {
				VERBOSE(*verb_out << "terminal score is too low - pruning!\n");
				delete tintrpr;
				return 0;
			}
			return tintrpr;
		}
	}
	return 0;
}

interpreter *
math_recognizer_base::getparser(const nonterminal *nt, const ordered_segments *segs) {
	nt_parse_table &nttab = tab[segs];
	nt_parse_table::iterator i = nttab.find(nt);
	return (i == nttab.end()) ? 0 : i->second.nt_parser;
}

/*
interpreter *
math_recognizer_base::getfullmux(const ordered_segments *segs) {
	nt_parse_table &nttab = tab[segs];
	interpreter *&intrp = nttab[0];
	if (!intrp) {
		multiplexor *mux = new multiplexor(this, 0, segs, 2);
		const grammar *G = GetMathGrammar();
		for (std::vector<nonterminal *>::const_iterator i = G->nts.begin(); i != G->nts.end(); ++i) {
			interpreter *child = getparser(*i, segs);
			if (child) {
				mux->addparser(child);
			}
		}
		intrp = mux;
	}
	return intrp;
}*/

static interpreter *
wrap(interpreter *intrpr, const nonterminal *nt) {
	multiplexor *mux = new multiplexor(intrpr->ctx(), nt, intrpr->span(), 2);
	mux->addparser(intrpr, false);
	return mux;
}

static interpreter *
build_multiline_interpreter(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *span, bool ign_overlap) {
	assert(nt->productions.size() == 2);

	VERBOSE(*verb_out << "build_multiline_interpreter(nt=" << nt->name << ", span=" << span->bits << ")\n");
	interpreter *single_expr_intrpr = rec->mkprodintrpr(nt->productions[0], span, ign_overlap);//rec->mkparser(nt, span);
	interpreter *multi_expr_intrpr = rec->mkprodintrpr(nt->productions[1], span, ign_overlap);
	
	if (!single_expr_intrpr && !multi_expr_intrpr) return 0;
	if (!single_expr_intrpr) return wrap(multi_expr_intrpr, nt);
	if (!multi_expr_intrpr) return wrap(single_expr_intrpr, nt);

	interpretation *top_single = getfirst(single_expr_intrpr);
	if (!top_single) {
		VERBOSE(*verb_out << " multiline interpreter is ONLY multi-line 2\n");
		return wrap(multi_expr_intrpr, nt);
	}
	basic_tree *single_tree = top_single->mktree();
	double single_score = single_tree->score();
	VERBOSE(*verb_out << " top single tree is " << single_tree->ccstr() << " with score " << single_score << std::endl);
	single_tree->release();

	scg::bitvec bits(span->bits.size(), true);
	int hibit = bits.highest_set_bit();
	if (hibit == -1) {
		VERBOSE(*verb_out << " multiline interpreter is ONLY single-line (wrapped) 1\n");
		return wrap(single_expr_intrpr, nt);
	}
	bits.set(hibit, false);
	const ordered_segments *prev_span = rec->getsegs(bits);
	if (!prev_span) {
		VERBOSE(*verb_out << " multiline interpreter is ONLY single-line (wrapped) 2\n");
		return wrap(single_expr_intrpr, nt);
	}

	VERBOSE(*verb_out << " prev span is " << bits << std::endl);
	interpreter *last_multi_intrpr = rec->mkprodintrpr(nt->productions[1], span, ign_overlap);
	/*if (last_multi_intrpr) {
		interpretation *last_multi_intrp = getfirst(last_multi_intrpr);
		if (last_multi_intrp) {
			basic_tree *multi_tree = last_multi_intrp->mktree();
			double multi_score = multi_tree->score();
			VERBOSE(*verb_out << " last multi interp is " << multi_tree->ccstr() << " with score " << multi_score << std::endl);
			multi_tree->release();

			if (multi_score > single_score) {
				VERBOSE(*verb_out << " multiline interpreter is sequence iterator for multi/single expressions\n");
				sequence *seq = new sequence(rec, root_nt, span, 3);
				seq->addparser(multi_intrpr);
				seq->addparser(single_expr_intrpr);
				return seq;
			}
		}
	}*/

	VERBOSE(*verb_out << " multiline interpreter is multiplexor for multi/single expressions\n");
	multiplexor *mux = new multiplexor(rec, nt, span, 3);
	mux->addparser(single_expr_intrpr, false);
	mux->addparser(multi_expr_intrpr, false);
	return mux;
}


interpreter *
mkntparser(math_recognizer_base *rec, const nonterminal *nt, const ordered_segments *segs, bool ign_overlap) {
	VERBOSE(*verb_out << "mkparser(" << nt->name << ":" << segs->bits << ")\n");
	if (segs->size() < nt->min_length || segs->size() > nt->max_length) {
		return 0;
	}

	assert((*rec->spans)[segs->bits] == segs);
	nt_parse_table &nttab = rec->tab[segs];
	nt_parse_table::iterator i = nttab.find(nt);
	if (i == nttab.end()) {
		production_parse_table &prod_tab = nttab[nt];
		interpreter *&intrpr = prod_tab.nt_parser;
		if (nt->isterminal()) {
			intrpr = rec->mktermparser(nt, segs);
		}
		else if (nt->name == "!MATRIX") {
			intrpr = new matrixparser(rec, nt, segs);
		}
		else if (nt->name == "!MULTI") {
			intrpr = new matrixparser(rec, nt, segs, MatrixAnalyzer::RECOGNIZE_ROWS);
		}
		else if (nt->name[0] == '&') {
			assert(nt->productions.size() == 2);
			intrpr = build_multiline_interpreter(rec, nt, segs, ign_overlap);
		}
		else {
			intrpr = build_ntintrpr(rec, prod_tab, nt, segs, ign_overlap);
			if (intrpr && intrpr->nknown() == 0) {
				interpretation *test = intrpr->next();
				if (!test) {
					delete intrpr;
					intrpr = 0;
				}
			}
		}
		return intrpr;
	}
	else {
		VERBOSE(*verb_out << " (exists at " << i->second.nt_parser << ")\n");
		return i->second.nt_parser;
	}
}


interpreter *
math_recognizer_base::mkparser(const nonterminal *nt, const ordered_segments *segs, bool ign_overlap) {
	return mkntparser(this, nt, segs, ign_overlap);
}

}
