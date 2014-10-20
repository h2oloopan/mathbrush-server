#include "grouping.h"
#include "segment.h"
#include "stroke.h"
#include "stkutils.h"
//#include "stroke-alg.h"
//#include "interp.h"
#include "group.h"
#include "symbols.h"
#include "symbol-bag.h"
#include "parms.h"
#include "normalize.h"
#include "interp.h"
#include "elastic.h"
#include "lsp.h"
#include "gridnorm.h"
#include "feat.h"
#include "set-utils.h"
#include "verb.h"
#include "mathrecognizer-private.h"
#include "utils.h"
#include "ink-io.h"

#include <vector>
#include <list>
#include <deque>
#include <fstream>
#include <cstring>


namespace scg {


static double SCORE_THRESHOLD = std::log(1e-1);//RegisterParameterDouble("ScoreThreshold", &SCORE_THRESHOLD);
/*static double SMALL_SYMBOL_THRESHOLD = 0;
static void calc_parms() {
	SMALL_SYMBOL_THRESHOLD = GetParameterDouble("SmallSymbolThreshold") * TABLETPC_DPI;
}
static int _e = RegisterParameterCallback(&calc_parms);*/

//const unsigned MaxGroupSize = GetParameterUnsigned("MaxGroupSize");
static double FeatureConsiderationThreshold = RegisterParameterDouble("FeatureConsiderationThreshold", &FeatureConsiderationThreshold);
static unsigned NumMatchesPerSegment = RegisterParameterUnsigned("NumMatchesPerSegment", &NumMatchesPerSegment);

stroke_feature_set default_feature_set;
stroke_feature_set small_feature_set;

static symbolbag *symbag = 0;

static std::map<std::string, unsigned> lsp_ncor;
static std::map<std::string, unsigned> dtw_ncor;
static std::map<std::string, unsigned> feat_ncor;

static unsigned writer_pkgid = -1;

void
SetWriterPkgid(unsigned id) {
	writer_pkgid = id;
}

const stroke_feature_set &
get_feature_set(const symbol *s) {
	if (s->is_small()) {
		return small_feature_set;
	}
	return default_feature_set;
}


static bool
compare_stroke_cmp_pair(const std::pair<double, bool> &left, const std::pair<double, bool> &right)
{
	return left.first < right.first;
}

/*
static int
rescale(NormalizedStroke *first, NormalizedStroke *last, const Rect<double> &bounds)
{
	Rect<double> src_bounds = bbox(first, last);

	double xscale = width(bounds) / width(src_bounds);
	double yscale = height(bounds) / height(src_bounds);

	while (first != last) {
		double *x = first->x;
		double *y = first->y;
		const double *endx = x + num_points(*first);
		while (x != endx) {
			*x = bounds.left + *x * xscale;
			*y = bounds.top + *y * yscale;
			++x;
			++y;
		}
		++first;
	}

	return 0;
}


static int
rescale(NormalizedStrokeGroup &strokes, const Rect<double> &bounds)
{
	return rescale(&*strokes.begin(), &*strokes.end(), bounds);
}*/


bool
prepare_to_match(group *input, const prototype *model, std::vector<size_t> &inputorder) {
	assert(input->strokes.nstrokes >= model->strokes.nstrokes);
	Rect<long> input_bounds = input->bounds;
	/*std::cout << "prepare_to_match:\n";
	std::cout << " model: " << model->strokes << std::endl;
	std::cout << " input: " << input->strokes << std::endl;
	*/
	//size_t nstrokes = last_input - first_input;
	std::vector<size_t> available_strokes;
	for (size_t i = 0; i < input->strokes.nstrokes; ++i) {
		available_strokes.push_back(i);
	}

	inputorder.resize(model->strokes.nstrokes);
	for (size_t mi = 0; mi < model->strokes.nstrokes; ++mi) {
		const RawStroke &mstroke = model->strokes.strokes[mi];
		double model_dx = mstroke.x[mstroke.npoints - 1] - mstroke.x[0];
		double model_dy = mstroke.y[mstroke.npoints - 1] - mstroke.y[0];

		double best_feature_score = std::numeric_limits<double>::infinity();
		std::vector<size_t>::iterator best_matching_stroke = available_strokes.end();
		for (std::vector<size_t>::iterator i = available_strokes.begin(); i != available_strokes.end(); ++i) {
			RawStroke &istroke = input->strokes.strokes[*i];

			double input_dx = istroke.x[istroke.npoints - 1] - istroke.x[0];
			double input_dy = istroke.y[istroke.npoints - 1] - istroke.y[0];

			double dot = input_dx * model_dx + input_dy * model_dy;
			if (dot < 0.0 && !stroke_is_closed(istroke, input_bounds)) {
				input_dx = -input_dx;
				input_dy = -input_dy;
				//std::cout << "reversing stroke " << *i << std::endl;
				istroke.reverse();
				input->nstrokes.strokes[*i].reverse();
			}

			double score = feature_match(default_feature_set, model->nstrokes.strokes[mi], input->nstrokes.strokes[*i]);
			if (score < best_feature_score) {
				best_feature_score = score;
				best_matching_stroke = i;
			}
		}

		if (best_matching_stroke == available_strokes.end()) {
			return false;
		}

		inputorder[mi] = *best_matching_stroke;
		available_strokes.erase(best_matching_stroke);
	}
	/*std::cout << "revised input: " << input->strokes;
	std::cout << "input order: ";
	for (size_t i = 0; i < model->strokes.nstrokes; ++i) {
		std::cout << inputorder[i] << ' ';
	}
	std::cout << std::endl;*/
	return true;
}


/*
struct correlation_matches_prototype_pred
{
	explicit correlation_matches_prototype_pred(const prototype *P_) : P(P_) { }

	inline bool operator()(const prototype_correlation &correlation) const
		{ return correlation.from == P; }

private:
	const prototype *P;
};*/

static double
pmap(double d, double d0, double p0, double d1, double p1, double dmax, double pmax) {
	if (d < d0) {
		return d/d0 * p0;
	}
	else if (d < d1) {
		return p0 + (d-d0)/(d1-d0) * (p1-p0);
	}
	else {
		return p1 + (d-d1)/(dmax-d1) * (pmax-p1);
	}
}

	const static double LSP_D0 = 0.85;
	const static double LSP_P0 = 0.1;
	const static double LSP_D1 = 1.65;
	const static double LSP_P1 = 0.97;
	const static double DTW_D0 = 0.35;
	const static double DTW_P0 = 0.035;
	const static double DTW_D1 = 0.95;
	const static double DTW_P1 = 0.98;
	const static double LSPM_D0 = 0.3;//0.1;
	const static double LSPM_P0 = 0.1;//0.03;//0.01;
	const static double LSPM_D1 = 1.2;//1.05;//0.68;
	const static double LSPM_P1 = 0.93;//0.85;//0.83;
	const static double DTWM_D0 = 0.1;//0.1;
	const static double DTWM_P0 = 0.01;
	const static double DTWM_D1 = 0.5;//0.5;
	const static double DTWM_P1 = 0.87;//0.87;
	const static double FEATM_D0 = 0.01;
	const static double FEATM_P0 = 0.01;
	const static double FEATM_D1 = 0.24;
	const static double FEATM_P1 = 0.93;
static double
blendscores(double lspd, double dtwd, double featd, unsigned npoints, unsigned nstrokes, const symbol *S) {
	VERBOSE2(*verb_out << "blending scores " << lspd << " and " << dtwd << " and " << featd << "...");
	featd = pmap(featd, FEATM_D0, FEATM_P0, FEATM_D1, FEATM_P1, 1, 1);
	lspd = pmap(lspd, LSPM_D0, LSPM_P0, LSPM_D1, LSPM_P1, 2, 1);//std::sqrt(LSPDEG/*+1*/));
	dtwd = pmap(dtwd, DTWM_D0, DTWM_P0, DTWM_D1, DTWM_P1, 2, 1);
	VERBOSE2(*verb_out << "map to " << lspd << " and " << dtwd << " and " << featd << std::endl);
	VERBOSE2(*verb_out << " (npoints=" << npoints << ", nstrokes=" << nstrokes << ")\n");
	unsigned lspnc = lsp_ncor[S->name];
	unsigned dtwnc = dtw_ncor[S->name];
	unsigned featnc = dtw_ncor[S->name];
	double lspfac, dtwfac, featfac;
	if (lspnc == 0 && dtwnc == 0 && featnc == 0) {
		lspfac = dtwfac = featfac = 1.0/3;
	}
	else {
		lspfac = (double)lspnc / (lspnc + dtwnc + featd);
		dtwfac = (double)dtwnc / (lspnc + dtwnc + featd);
		featfac = (double)featnc / (lspnc + dtwnc + featd);
	}
	//lspfac = 1; dtwfac = 0;
	//dtwfac = 1; lspfac = 0;
	//lspfac = dtwfac = 0.5;
	//featfac = 0;
	//if (npoints < 32) {
		lspfac = dtwfac = 0.4;
	/*}
	else if (npoints > 64) {
		lspfac = 0;
		dtwfac = 0.8;
	}
	else {
		lspfac = 0.4-(npoints-32.0)/32.0;
		dtwfac = 0.8-lspfac;
	}*/
	/*if (nstrokes == 1) {
		lspfac = 0.5;
		dtwfac = 0.3;
	}
	else if (nstrokes == 2) {
		lspfac = dtwfac = 0.4;
	}
	else {
		lspfac = 0.1;
		dtwfac = 0.7;
	}*/
	featfac = 0.2;
	//featfac = 0;
	//dtwfac = std::exp(-std::log(2.0)/24 * npoints);
	//lspfac = 1 - dtwfac;
	//dtwfac *= 0.8; lspfac *= 0.8;
	double d = lspfac * lspd + dtwfac * dtwd + featfac * featd;
	VERBOSE2(*verb_out << " combine with factors " << lspfac << " and " << dtwfac << " and " << featfac << " for score " << d << std::endl);
	return d;//std::sqrt(d);
}

match_score
prune_and_match_symbol(symbol *S, group *grp) {
	double grplen = 0;
	for (size_t i = 0; i < grp->nstrokes.nstrokes; ++i) {
		grplen = std::max(grplen, grp->nstrokes.strokes[i].length());
	}
	static std::vector<size_t> inputorder(1, 0);
	match_score match;
	if (S->is_small()) {
		double best_score = std::numeric_limits<double>::infinity();
		double best_fscore = 0;
		const prototype *best_proto = 0;
		for (prototype *P = S->firstproto(); P; P = S->nextproto(P)) {
			const RawStrokeGroup &model = P->strokes;
			if (model.nstrokes != grp->children.size()) {
				continue;
			}
			if (model.nstrokes > inputorder.size()) {
				size_t s = inputorder.size();
				inputorder.resize(model.nstrokes);
				for (size_t i = s; i < inputorder.size(); ++i) {
					inputorder[i] = i;
				}
			}
			double fs = feature_match(small_feature_set, grp->nstrokes, P->nstrokes, inputorder);//match_group(model, grp->strokes, S.info);
			fs = std::min(fs, 2 * FeatureConsiderationThreshold);
			fs = 1.0 - fs / (2*FeatureConsiderationThreshold);
			VERBOSE2(*verb_out << " matching vs " << S->name << "...\n");
			double s;
			gridnorm_match(P->nstrokes, grp->nstrokes, inputorder, s);
			if (s < best_score) {
				best_fscore = fs;
				best_score = s;
				best_proto = P;
			}
			VERBOSE2(*verb_out << "   match score " << s << std::endl);
		}
		if (best_proto) {
			match.feature_score = best_fscore;
			match.matcher_score = best_score;
			match.S = S;
		}
	}
	else {
		static const size_t NSCORES = 2;
		double lspscore = std::numeric_limits<double>::infinity();
		double dtwscore = std::numeric_limits<double>::infinity();
		double featscore = std::numeric_limits<double>::infinity();
		double dtwds[NSCORES], lspds[NSCORES], featds[NSCORES];
		std::fill(dtwds, dtwds+NSCORES, std::numeric_limits<double>::infinity());
		std::fill(featds, featds+NSCORES, std::numeric_limits<double>::infinity());
		std::fill(lspds, lspds+NSCORES, std::numeric_limits<double>::infinity());
		double fscore = 0;
		unsigned n = 0;
		const prototype *best_proto = 0;
		unsigned best_npts = 0;
		unsigned nwriterprotos = 0;
		/*
		if (writer_pkgid != -1) {
			for (symbol::const_iterator pproto = S.begin(); pproto != S.end(); ++pproto) {
				const prototype *P = *pproto;
				const RawStrokeGroup &model = P->strokes;
				if (model.nstrokes == grp->children.size() && P->id.pkgid == writer_pkgid) {
					if (++nwriterprotos == 5) {
						break;
					}
				}
			}
		}*/
		for (prototype *P = S->firstproto(); P; P = S->nextproto(P)) {
			//if (P->meanfeatd == -1 || P->varfeatd == 0 || P->varlspd == 0 || P->vardtwd == 0) continue;
			const RawStrokeGroup &model = P->strokes;
			if (model.nstrokes != grp->children.size()) {
				continue;
			}
			if (model.nstrokes > inputorder.size()) {
				inputorder.resize(model.nstrokes);
			}
			if (prepare_to_match(grp, P, inputorder)) {
				double modlen = 0;
				for (size_t i = 0; i < P->nstrokes.nstrokes; ++i) {
					modlen = std::max(modlen, P->nstrokes.strokes[i].length());
				}
				VERBOSE2(
					*verb_out << "  matching vs " << S->name << "...grplen " << grplen << " ; modlen " << modlen << std::endl;
					*verb_out << "inputorder = [ ";
					for (std::vector<size_t>::const_iterator kk = inputorder.begin(); kk != inputorder.end(); ++kk) {
						*verb_out << *kk << ' ';
					}
					*verb_out << "]\n";
				);
				if (std::max(modlen/grplen, grplen/modlen) > 2) {
					continue;
				}
				/*std::vector<size_t> modelorder(model.nstrokes);
				for (size_t i = 0; i < model.nstrokes; ++i) {
					modelorder[inputorder[i]] = i;
				}*/
				double featd = feature_match(default_feature_set, grp->nstrokes, P->nstrokes, inputorder);
				//VERBOSE(*verb_out << " feature score for " << info.name << " (" << P->id << ") is " << featd << std::endl);
				double fs = 1.0 - featd / FeatureConsiderationThreshold;
				//if (fs > 0) {
					fscore = std::max(fscore, fs);
					double lspd = lspdist(P->strokes, grp->strokes, inputorder);
					
					NormalizedStrokeGroup subdiv = subdivide(grp->nstrokes, P->nstrokes, inputorder);
					double dtwd = elasticdist(P->nstrokes, subdiv, inputorder);
					VERBOSE2(*verb_out << "  feat " << featd << ", LSP " << lspd << ", DTW " << dtwd << std::endl);
					if (lspd < lspscore) {
						//lspscore = lspd;
						best_proto = P;
						best_npts = grp->nstrokes.strokes[0].npoints;
						for (size_t i = 1; i < grp->nstrokes.nstrokes; ++i) {
							best_npts = std::max(best_npts, (unsigned)grp->nstrokes.strokes[i].npoints);
						}
					}
					/*if (dtwd < dtwscore) {
						dtwscore = dtwd;
						//best_npts = ((subdiv.npoints()/subdiv.nstrokes) + (P->nstrokes.npoints()/P->nstrokes.nstrokes)) / 2;
					}
					if (featd < featscore) {
						featscore = featd;
					}*/
					for (size_t i = 0; i < NSCORES; ++i) {
						if (dtwd < dtwds[i]) {
							std::memmove(dtwds + i+1, dtwds+i, sizeof(*dtwds)*(NSCORES-i-1));
							dtwds[i] = dtwd;
							break;
						}
					}
					for (size_t i = 0; i < NSCORES; ++i) {
						if (lspd < lspds[i]) {
							std::memmove(lspds + i+1, lspds+i, sizeof(*lspds)*(NSCORES-i-1));
							lspds[i] = lspd;
							break;
						}
					}
					for (size_t i = 0; i < NSCORES; ++i) {
						if (featd < featds[i]) {
							std::memmove(featds + i+1, featds+i, sizeof(*featds)*(NSCORES-i-1));
							featds[i] = featd;
							break;
						}
					}
				/*}
				else {
					VERBOSE4(*verb_out << "   pruned\n");
				}*/
			}
		}
		if (best_proto) {
			VERBOSE2(*verb_out << " matching vs " << S->name << std::endl);
			double dtwsc = 0, lspsc = 0, featsc = 0;
			size_t i;
			VERBOSE2(*verb_out << " DTW scores ");
			for (i = 0; i < NSCORES; ++i) {
				if (dtwds[i] == std::numeric_limits<double>::infinity()) {
					break;
				}
				dtwsc += dtwds[i];
				VERBOSE2(*verb_out << dtwds[i] << ' ');
			}
			if (i > 1) {
				dtwsc /= i;
				VERBOSE2(*verb_out << " -> " << dtwsc << std::endl);
				VERBOSE2(*verb_out << " Feat scores ");
				for (i = 0; i < NSCORES; ++i) {
					if (featds[i] == std::numeric_limits<double>::infinity()) {
						break;
					}
					featsc += featds[i];
					VERBOSE2(*verb_out << featds[i] << ' ');
				}
				if (i > 1) featsc /= i;
				VERBOSE2(*verb_out << " -> " << featsc << std::endl);
				VERBOSE2(*verb_out << " LSP scores ");
				for (i = 0; i < NSCORES; ++i) {
					if (lspds[i] == std::numeric_limits<double>::infinity()) {
						break;
					}
					lspsc += lspds[i];
					VERBOSE2(*verb_out << lspds[i] << ' ');
				}
				if (i > 1) lspsc /= i;
				VERBOSE2(*verb_out << " -> " << lspsc << std::endl);
				match.feature_score = fscore;
				match.matcher_score = blendscores(lspsc, dtwsc, featsc, best_npts, grp->strokes.nstrokes, S);
				match.S = S;
			}
		}
	}
	return match;
}

static int
prune_and_match_group(group *grp) {
	int e = 0;
	VERBOSE(*verb_out << "prune_and_match on group " << grp->bits << std::endl);
	for (symbol *S = symdb_firstsymbol(); S; S = symdb_nextsymbol(S)) {
		match_score match = prune_and_match_symbol(S, grp);
		if (match.S) {// && match.matcher_score < 0.8) {
			/*
			size_t n = match.proto->strokes.nstrokes;
			if (n > 1) {
				double msc = std::pow(0.95, n-1);
				match.matcher_score *= msc;
			}*/
			//grp->feature_score = std::max(grp->feature_score, match.feature_score);
			grp->feature_score = std::max(grp->feature_score, 1-match.matcher_score);
			grp->matches[S->nt] = match;//.push_back(match);
		}
	}
	VERBOSE(*verb_out << "final feature score for group " << grp->bits << " is " << grp->feature_score << std::endl);

	return 0;
}

std::map<bitvec, group *>::iterator
math_recognizer_base::remove_group(group *grp)
{
	VERBOSE(*verb_out << "remove_group " << grp->bits << std::endl);
	for (std::list<segment *>::iterator i = grp->children.begin(); i != grp->children.end(); ++i) {
		segment *child = *i;
		std::list<group *>::iterator j = std::find(child->parents.begin(), child->parents.end(), grp);
		assert(j != child->parents.end());
		child->parents.erase(j);
	}

	std::map<bitvec, group *>::iterator i = groups->find(grp->bits);
	assert(i != groups->end());
	groups->erase(i);
	i = groups->upper_bound(grp->bits);

	std::list<group *>::iterator j = std::find(newgroups.begin(), newgroups.end(), grp);
	if (j != newgroups.end()) {
		newgroups.erase(j);
	}

	delete grp;
	return i;
}

static void
clear_table(nt_parse_table &tab)
{
	for (nt_parse_table::iterator j = tab.begin(); j != tab.end(); ++j) {
		delete j->second.nt_parser;
		production_parser_map &prod_parsers = j->second.tab;
		for (production_parser_map::iterator k = prod_parsers.begin(); k != prod_parsers.end(); ++k) {
			delete k->second;
		}
	}
	tab.clear();
}

int
math_recognizer_base::update_spans_on_removal(const segment *seg) {
	std::map<bitvec, ordered_segments *> *newspans = new std::map<bitvec, ordered_segments *>;

	for (std::map<bitvec, ordered_segments *>::const_iterator i = spans->begin(); i != spans->end(); ++i) {
		ordered_segments *span = i->second;
		if (span) {
			bitvec &bits = span->bits;
			//VERBOSE(*verb_out << "span " << span << " (" << bits << ") -> ");
			if (seg->pos < bits.size()) {
				std::list<group *>::const_iterator j;
				for (j = seg->parents.begin(); j != seg->parents.end(); ++j) {
					if (!null_intersection(bits, (*j)->bits)) break;
				}
				//if (!bits.at(pos)) {
				if (j == seg->parents.end()) {
					//VERBOSE(*verb_out << bits << " : erase bit " << seg->pos << std::endl);
					bits.erase(seg->pos);
					(*newspans)[bits] = span;
				}
				else {
					//VERBOSE(*verb_out << bits << " : kill parse table\n");
					parse_table::iterator j = tab.find(span);
					if (j != tab.end()) {
						clear_table(j->second);
						tab.erase(j);
					}
					delete i->second;
				}
			}
			else {
				//VERBOSE(*verb_out << bits << " : no change\n");
				(*newspans)[bits] = span;
			}
		}
	}

	delete spans;
	spans = newspans;
		
	return 0;
}

void
math_recognizer_base::update_groups_on_removal(size_t pos)
{
	std::map<bitvec, group *> *newgroups = new std::map<bitvec, group *>;

	for (std::map<bitvec, group *>::const_iterator i = groups->begin(); i != groups->end(); ++i) {
		group *grp = i->second;
		if (grp) {
			bitvec &bits = grp->bits;
			bits.erase(pos);
			(*newgroups)[bits] = grp;
		}
	}

	delete groups;
	groups = newgroups;

	VERBOSE2(
		*verb_out << "new groups are:\n";
		for (std::map<bitvec, group *>::const_iterator i = groups->begin(); i != groups->end(); ++i) {
			*verb_out << "  " << i->first << " -> " << i->second << " at " << i->second->bits << std::endl;
		}
	);
}


size_t segment::INVALID_SEGMENT_POS;

void
math_recognizer_base::remove_segment(size_t segi, bool kill_strokes) {
	VERBOSE(*verb_out << "remove_segment " << segi << std::endl);
	segment *seg = segments[segi];
	if (kill_strokes) {
		for (std::list<stroke *>::iterator i = seg->children.begin(); i != seg->children.end(); ++i) {
			std::vector<stroke *>::iterator j = std::find(strokes.begin(), strokes.end(), *i);
			if (j != strokes.end()) {
				delete *i;
				strokes.erase(j);
			}
		}
	}
	
	update_spans_on_removal(seg);

	while (!seg->parents.empty()) {
		remove_group(seg->parents.front());
	}

	update_groups_on_removal(seg->pos);

	//append_edit_op(ctx, edit_op::REMOVAL, seg->pos);

	segments.erase(segments.begin() + segi);

	delete seg;

	set_segment_positions();
}


group *
math_recognizer_base::create_default_group(segment *seg, size_t segi)
{
	group *grp = new group;
	grp->children.push_back(seg);
	seg->parents.push_back(grp);
	grp->bounds = seg->bounds;
	bitvec bits(segments.size(), false);
	bits.set(segi, true);
	grp->bits = bits;
	grp->feature_score = 0.0;
	grp->proximity_score = 1.0;
	if (seg->stk.npoints > 0) {
		grp->strokes = build_stroke_group(&seg, &seg + 1, seg->bounds, 1);
		//NormalizedStrokeGroup ns = normalize(grp->strokes);
		//std::cout << "default group gives n stroke " << ns << std::endl;
		//NormalizedStrokeGroup ss = subdivide(ns);
		grp->nstrokes = normalize(grp->strokes);
		//grp->nstrokes = subdivide(ns);
		//std::cout << "default group gives ns stroke " << grp->nstrokes << std::endl;
	}
	return grp;
}


segment *
create_default_segment(stroke *stk)
{
	segment *seg = new segment;
	if (stk) {
		stk->parent = seg;
		seg->children.push_back(stk);
		seg->stk = stk->input.copy();
		seg->bounds = seg->stk.bounds();
	}
	return seg;
}

void
math_recognizer_base::rebalance_groups() {
	VERBOSE2(*verb_out << "rebalancing...\n");
	//double complprod = 1.0;
	std::map<group *, double> maxsupgrp;
	for (std::map<bitvec, group *>::iterator i = groups->begin(); i != groups->end(); ++i) {
		//i->second->weight = 1.0;
		/*double sc = std::max(i->second->proximity_score, i->second->stack_score);
		if (sc < 1) {
			complprod *= 1.0 - std::max(i->second->proximity_score, i->second->stack_score);
		}*/
		VERBOSE2(*verb_out << " finding max for grp " << i->first << std::endl);
		group *grp = i->second;
		for (std::map<bitvec, group *>::iterator j = groups->begin(); j != groups->end(); ++j) {
			group *tst = j->second;
			if (grp == tst) continue;
			if (grp->bits.subset_of(tst->bits)) {
				VERBOSE2(*verb_out << "  found hit " << tst->bits << " with score " << std::max(tst->proximity_score, tst->stack_score) << std::endl);
				maxsupgrp[grp] = std::max(maxsupgrp[grp], /*tst->feature_score* */std::max(tst->proximity_score, tst->stack_score));
			}
		}
	}
	for (std::map<bitvec, group *>::iterator i = groups->begin(); i != groups->end(); ++i) {
		group *grp = i->second;
		grp->weight = 1.0 - maxsupgrp[grp];
		double gscore = std::max(grp->stack_score, grp->proximity_score);
		gscore = std::min(1.0, gscore);
		double gs = gscore;
		for (size_t j = 1; j < grp->bits.count_set_bits(); ++j) {
			gscore *= gs;
		}
		grp->pnil = 1.0 - gscore * grp->weight * grp->feature_score;
		grp->pnil = std::max(0.0, std::min(1.0, grp->pnil));
		VERBOSE2(*verb_out << "weight for grp " << i->first << " is " << grp->weight << " and pnil " << grp->pnil << std::endl);
		//grp->proximity_score = 1.0 - 0.75*maxsupgrp[grp];
		//grp->stack_score = 1.0 - 0.75*maxsupgrp[grp];
	}
}


void
math_recognizer_base::finalize_match_scores(group *grp) {
	VERBOSE(*verb_out << "finalizing scores for group " << grp->bits << std::endl);
	if (grp->tags & group_tags::KNOWN_SYMBOL) {
		return;
	}

	VERBOSE(*verb_out << "average size is " << avg_stroke_size << /*"; weighted avg is " << wavg_size << */std::endl);

	double SMALL_SYMBOL_THRESHOLD = 0.125 * TABLETPC_DPI;
	double avg_size = std::max(avg_stroke_size, SMALL_SYMBOL_THRESHOLD);
	VERBOSE(*verb_out << " -> " << avg_size << std::endl);
	// weight small symbols by exponential decay:
	//  height = 0 => weight 1
	//  height = alpha * avg => weight 1/2
	double alpha = 0.125;
	double L = std::log(2.0) / (alpha * avg_size);
	double sz = std::max(grp->bounds.width(), grp->bounds.height());
	double small_weight = std::exp(-L * sz);

	grp->weight = std::min(1.0, std::max(0.0, grp->weight));
		
	VERBOSE(*verb_out << "  begin group " << grp->bits << " of size " << grp->children.size() << " and weight " << grp->weight << " at " << grp->bounds << ":" << std::endl);
	//std::cout << "  begin group " << grp->bits << " of size " << grp->children.size() << " and weight " << grp->weight << " at " << grp->bounds << ":" << std::endl;
	VERBOSE(*verb_out << "scores: prox " << grp->proximity_score << " stack " << grp->stack_score << " feat " << grp->feature_score << " boost " << grp->boosted_score << std::endl);

	//double small_weight = 1.0 - std::min(1.0, (double)sz/(2*avg_size));
	//double small_weight = 1.0 - std::min((double)sz / SMALL_SYMBOL_THRESHOLD, 1.0);//avg_size, 1.0);
	VERBOSE(*verb_out << "size here is " << sz << " vs. average " << avg_size << std::endl);
	VERBOSE(*verb_out << "   small symbol weight here is " << small_weight << std::endl);
	VERBOSE(*verb_out << "   there are " << grp->matches.size() << " matches\n");

	static const size_t NCANDIDATES = 5;
	match_score bestmatches[NCANDIDATES];
	/*for (size_t j = 0; j < NCANDIDATES; ++j) {
		bestmatches[j].S = 0;
		bestmatches[j].final_score = 0;
		bestmatches[j].matcher_score = std::numeric_limits<double>::infinity();
	}*/
	for (std::map<const nonterminal *, match_score>::iterator j = grp->matches.begin(); j != grp->matches.end(); ++j) {
		match_score &match = j->second;
		match.matcher_score = 0.1 + 0.9 * match.matcher_score;//std::max(match.matcher_score, 0.1);
		VERBOSE(*verb_out << match.S->name << ":  ");
		double w = (match.S->is_small()) ? small_weight : 1.0 - small_weight;
		match.score = w / (match.matcher_score * match.matcher_score);
		VERBOSE(*verb_out << " matcher " << match.matcher_score);
		VERBOSE(*verb_out << " -> small " << w);
		if (symbag) {
			match.bag_score = symbag->query(match.S);
			VERBOSE(*verb_out << " -> bag " << match.bag_score);
		}
		VERBOSE(*verb_out << " -> " << match.score * match.bag_score << std::endl);

		match_score *sc;
		for (sc = bestmatches; sc != bestmatches + NCANDIDATES; ++sc) {
			if (!sc->S || match.score * match.bag_score > sc->score * sc->bag_score) break;
			/*double scw = 1;
			if (sc->S) scw = sc->S->is_small() ? small_weight : 1.0 - small_weight;
			if (w / match.matcher_score > scw / sc->matcher_score) break;*/
		}
		if (sc != bestmatches + NCANDIDATES) {
			memmove(sc + 1, sc, sizeof(*bestmatches) * (bestmatches + NCANDIDATES - sc - 1));
			*sc = match;
		}
	}

	double scoresum = 0;
	double bagscores[NCANDIDATES] = {1.0};
	double bagscoresum = 0;
	for (size_t j = 0; bestmatches[j].S && j < NCANDIDATES; ++j) {
		scoresum += bestmatches[j].score;
		bagscoresum += bestmatches[j].score * bestmatches[j].bag_score;
		/*
		double smallscore = (bestmatches[j].S->is_small()) ? small_weight : 1.0 - small_weight;
		scoresum += smallscore / (bestmatches[j].matcher_score * bestmatches[j].matcher_score);
		if (symbag) {
			bagscores[j] = symbag->query(bestmatches[j].S);
			bagscoresum += bagscores[j];
		}*/
	}

	VERBOSE(*verb_out << "FINAL matches for group " << grp->bits << " with weight " << grp->weight << " and pnil " << grp->pnil << std::endl);
	grp->matches.clear();
	for (size_t j = 0; bestmatches[j].S && j < NCANDIDATES; ++j) {
		const symbol *S = bestmatches[j].S;
		if (symbag) {
			symbag->addsymbol(bestmatches[j].S, (1.0 - grp->pnil) * bestmatches[j].score / scoresum);
		}
		double prob = bestmatches[j].score * bestmatches[j].bag_score / bagscoresum;
		bestmatches[j].final_score.score = std::log(prob);
		VERBOSE(*verb_out << "  " << S->name << " -> " << prob << " -> " << bestmatches[j].final_score.score << std::endl);
		/*double smallscore = (S->is_small()) ? small_weight : 1.0 - small_weight;
		VERBOSE(*verb_out << "   for best match to " << S->name << ", matcher score " << bestmatches[j].matcher_score << " and small score " << smallscore << " -> ");
		bestmatches[j].final_score.score = smallscore / (bestmatches[j].matcher_score * bestmatches[j].matcher_score * scoresum);
		if (symbag) symbag->addsymbol(S, (1.0 - grp->pnil) * bestmatches[j].final_score.score);
		VERBOSE(*verb_out << bestmatches[j].final_score.score << " -> ");
		bestmatches[j].final_score.score = std::log(bestmatches[j].final_score.score);
		if (symbag) {
			VERBOSE(*verb_out << "symbag score " << bagscores[j] / bagscoresum << " -> ");
			bestmatches[j].final_score.score += std::log(bagscores[j] / bagscoresum);
		}
		VERBOSE(*verb_out << bestmatches[j].final_score.score << std::endl);*/

		if (bestmatches[j].final_score.score > SCORE_THRESHOLD) {
			//VERBOSE(*verb_out << "    final score = " << 1.0 / (bestmatches[j].matcher_score * scoresum) << " / " << pnil << " = LOG " << bestmatches[j].final_score.score << std::endl);
			grp->matches[bestmatches[j].S->nt] = bestmatches[j];
		}
	}

	VERBOSE2(*verb_out << "updated symbol bag to:\n" << *symbag << std::endl);
}

template <typename T>
std::map<bitvec, T *> *
append_segment(std::map<bitvec, T *> &base) {
	std::map<bitvec, T *> *newts = new std::map<bitvec, T *>;

	for (typename std::map<bitvec, T *>::const_iterator i = base.begin(); i != base.end(); ++i) {
		T *t = i->second;
		if (t) {
			bitvec &bits = t->bits;
			bits.insert(bits.size(), false);
			(*newts)[bits] = t;
		}
	}

	return newts;
}

void
math_recognizer_base::add_segment() {
	std::map<bitvec, group *> *newgroups = append_segment(*groups);
	delete groups;
	groups = newgroups;

	//std::map<bitvec, ordered_segments *> *newspans = append_segment(*ctx->spans);
	//delete ctx->spans;
	//ctx->spans = newspans;

	VERBOSE(
		*verb_out << "new groups are:\n";
		for (std::map<bitvec, group *>::const_iterator i = groups->begin(); i != groups->end(); ++i) {
			*verb_out << "  " << i->first << " -> " << i->second << " at " << i->second->bits << std::endl;
		}
	);

	set_segment_positions();
}



void
math_recognizer_base::set_segment_positions() {
	size_t n = 0;
	for (std::vector<segment *>::iterator i = segments.begin(); i != segments.end(); ++i) {
		/*
		for (std::list<group *>::iterator j = (*i)->parents.begin(); j != (*i)->parents.end(); ++j) {
			if (ctx->groups->find((*j)->bits) == ctx->groups->end()) {
				std::abort();
			}
		}*/
		(*i)->pos = n++;
	}
}


/*void
dump_ink_tree(std::ostream &os, context *ctx)
{
	os << "INK TREE:\n";
	for (std::vector<stroke *>::const_iterator i = ctx->strokes.begin(); i != ctx->strokes.end(); ++i) {
		os << "stroke " << i - ctx->strokes.begin() << " : " << *i << std::endl;
		os << "  parent : " << (*i)->parent << std::endl;
	}

	os << std::endl;

	for (std::vector<segment *>::const_iterator i = ctx->segments.begin(); i != ctx->segments.end(); ++i) {
		const segment *seg = *i;
		os << "segment " << i - ctx->segments.begin() << " : " << seg << " at " << seg->bounds << std::endl;
		os << "  children : ";
		for (std::list<stroke *>::const_iterator j = seg->children.begin(); j != seg->children.end(); ++j) {
			os << *j << ' ';
		}
		os << std::endl << "  parents : ";
		for (std::list<group *>::const_iterator j = seg->parents.begin(); j != seg->parents.end(); ++j) {
			os << *j << ' ';
		}
		os << std::endl;
	}

	os << std::endl;

	for (std::map<bitvec, group *>::const_iterator i = ctx->groups->begin(); i != ctx->groups->end(); ++i) {
		const group *grp = i->second;
		os << "group " << grp << " at " << grp->bits << " / " << grp->bounds << " with weight " << grp->weight << std::endl;
		os << "scores: prox " << grp->proximity_score << " stack " << grp->stack_score << " feat " << grp->feature_score << " boost " << grp->boosted_score << std::endl;
		os << "  children : ";
		for (std::list<segment *>::const_iterator j = grp->children.begin(); j != grp->children.end(); ++j) {
			os << *j << ' ';
		}
		/*
		os << std::endl << "  matches : ";
		for (std::vector<match_score>::const_iterator j = grp->matches.begin(); j != grp->matches.end(); ++j) {
			os << "    " << j->proto->info().name << " -> " << j->score << std::endl;
		}* /
		os << std::endl << "  final matches : ";
		for (std::vector<match_score>::const_iterator j = grp->final_matches.begin(); j != grp->final_matches.end(); ++j) {
			os << "    " << j->proto->info().name << " -> " << j->score << std::endl;
		}
	}
}*/


/*
static void
swap_cases(group *grp) {
	struct swap_t {
		std::string s1;
		std::string s2;
		swap_t(const std::string &s1_, const std::string &s2_) : s1(s1_), s2(s2_) { }
	};

	static swap_t swaps[] = { swap_t("C","c"), swap_t("K","k"), swap_t("O","o"), swap_t("P","p"),
		                      swap_t("S","s"), swap_t("V","v"), swap_t("W","w"), swap_t("X","x"),
							  swap_t("Y","y"), swap_t("Z","z") };
	static size_t nswaps = sizeof(swaps) / sizeof(*swaps);
	
	for (std::vector<match_score>::iterator i = grp->final_matches.begin(); i != grp->final_matches.end(); ++i) {
		const prototype *proto = i->proto;
		for (swap_t *s = swaps; s != swaps + nswaps; ++s) {
			if (proto->info().name == s->s2) {
				for (std::vector<match_score>::iterator j = grp->final_matches.begin(); j != i; ++j) {
					const prototype *proto2 = j->proto;
					if (proto2->info().name == s->s1) {
						VERBOSE(*verb_out << "For group " << grp->bits << " found swap between " << s->s1 << " and " << s->s2 << " with scores " << j->score << " vs. " << i->score << std::endl);
						j->proto = proto;
						i->proto = proto2;
						goto swapped;
					}
				}
			}
		}
swapped:
		;
	}
}*/

int
math_recognizer_base::add_strokes(const RawStroke *strokes, size_t nstrokes) {
	if (nstrokes == 0) {
		return 0;
	}

	VERBOSE(*verb_out << "adding " << nstrokes << " strokes\n");

	size_t newsegstart = segments.size();

	avg_stroke_size *= segments.size();

	for (const RawStroke *s = strokes; s != strokes + nstrokes; ++s) {
		if (s->npoints == 1) {
			continue;
		}

		add_segment();
		//VERBOSE(*verb_out << "adding stroke " << *s << std::endl);
	 	stroke *stk = new stroke;
#ifdef IPAD_RECOGNIZER
		stk->input = triple(*s);
#else
		stk->input = s->copy();
#endif
		this->strokes.push_back(stk);

		segment *seg = create_default_segment(stk);
		//seg->ctx = this;
		VERBOSE(*verb_out << "created segment " << seg);
		seg = merge_segment(seg);
		VERBOSE(*verb_out << " merged to " << seg << std::endl);

		segments.push_back(seg);
		double sz = std::max(seg->bounds.width(), seg->bounds.height());
		avg_stroke_size += sz;

		updategroups(segments.size()-1);
	}

	avg_stroke_size /= segments.size();

	VERBOSE(*verb_out << "after add, groups are\n";
		for (std::map<bitvec, group *>::const_iterator i = groups->begin(); i != groups->end(); ++i) {
			*verb_out << i->first << std::endl;
		}
	);
	handlenewgroups();

	/*
	for (size_t j = newsegstart; j < segments.size(); ++j) {
		const segment *cmpseg = segments[j];
		for (std::map<const ordered_segments *, lockspec>::iterator i = locks.begin(); i != locks.end(); ++i) {
			const ordered_segments *lockspan = i->first;
			const lockspec &ls = i->second;
			const interpreter *lockintrpr = ls.intrpr;
			const segment *minseg[2];
			const segment *maxseg[2];
			minseg[0] = lockspan->min(0);
			minseg[1] = lockspan->min(1);
			maxseg[0] = lockspan->max(0);
			maxseg[1] = lockspan->max(1);
			if ((*orders[0])(minseg[0], cmpseg) && (*orders[0])(cmpseg, maxseg[0])) {
				break;
			}
			if ((*orders[1])(minseg[1], cmpseg) && (*orders[1])(cmpseg, maxseg[1])) {
				break;
			}
			if (j < segments.size()) {
				// new stroke is in the middle of locked span, so transfer the lock
				bitvec bits(segments.size(), false);
				bits.union_insert(lockintrpr->span()->bits);
				bits.set(j);
				const ordered_segments *newsegs = getsegs(bits);
				interpreter *newintrpr = mkparser(lockintrpr->nt(), newsegs);
				unlock(lockintrpr);
				lock(newintrpr);
			}
		}
	}*/

	return 0;
}

void
math_recognizer_base::handlenewgroups() {
	/*
	for (std::list<group *>::iterator i = newgroups.begin(); i != newgroups.end(); ++i) {
		dilute_group_by_dist(*i, *groups);
	}
	for (std::list<group *>::iterator i = newgroups.begin(); i != newgroups.end(); ++i) {
		group *grp = *i;
		grp->proximity_score *= grp->dil_factor;
		grp->stack_score *= grp->dil_factor;
	}*/

	rebalance_groups();

	while (!newgroups.empty()) {
		group *grp = newgroups.front();
		newgroups.pop_front();
		std::map<bitvec, group *>::const_iterator i = groups->find(grp->bits);
		assert(i != groups->end() && i->second == grp);
		prune_and_match_group(grp);
		if (grp->matches.empty()) {
			remove_group(grp);
		}
		else {
			double gscore = std::max(grp->stack_score, grp->proximity_score);
			gscore = std::min(1.0, gscore);
			double gs = gscore;
			for (size_t j = 1; j < grp->bits.count_set_bits(); ++j) {
				gscore *= gs;
			}
			grp->pnil = 1.0 - gscore * grp->weight * grp->feature_score;
			grp->pnil = std::max(0.0, std::min(1.0, grp->pnil));
			finalize_match_scores(grp);
		}
	}

	//incorporate_match_correlations(*global_symbols_db());

	/*for (std::list<group *>::iterator i = newgroups.begin(); i != newgroups.end(); ++i) {
		swap_cases(*i);
	}*/

	//VERBOSE(dump_ink_tree(*verb_out, ctx));

	//push_matches_to_parse_table(ctx);

	set_segment_positions();
}

void
math_recognizer_base::clear_dead_groups() {
	std::list<group *> rmgroups;
	for (std::map<bitvec, group *>::iterator i = groups->begin(); i != groups->end(); ++i) {
		group *grp = i->second;
		if (grp->proximity_score == 0 && grp->stack_score == 0) {
			rmgroups.push_back(grp);
		}
	}

	for (std::list<group *>::iterator i = rmgroups.begin(); i != rmgroups.end(); ++i) {
		remove_group(*i);
	}
}

void
math_recognizer_base::updategroups(size_t segi) {
	std::map<bitvec, group *> prox_groups;

	segment *seg = segments[segi];

	group *single_seg_group = create_default_group(seg, segi);
	symbol *sqrt = symdb_findsymbol_name("sqrt");
	//symbol *horzline = symdb_findsymbol_name("horzline");
	match_score match = prune_and_match_symbol(sqrt, single_seg_group);
	//match_score hmatch = prune_and_match_symbol(horzline, single_seg_group);
	if (match.S) {
		/*
		if (hmatch.proto) {
			single_seg_group->container_score = std::max(0.0, hmatch.matcher_score-match.matcher_score);
			VERBOSE(*verb_out << "container score for seg.group " << single_seg_group->bits << " is " << match.matcher_score << " - " << hmatch.matcher_score << " = " << single_seg_group->container_score << std::endl);
		}
		else*/ {
			single_seg_group->container_score = 1-match.matcher_score;
			VERBOSE(*verb_out << "container score for seg.group " << single_seg_group->bits << " is " << single_seg_group->container_score << std::endl);
		}
	}
	//single_seg_group->proximity_score = 1;

	double best_score = 0.0;
	double allscores = 1.0;
	std::list<group *> all_new_groups;
	update_proximity_groups(*groups, prox_groups, single_seg_group, segi);
	std::map<bitvec, group *> new_prox_groups;
	while (!prox_groups.empty()) {
		for (std::map<bitvec, group *>::iterator i = prox_groups.begin(); i != prox_groups.end(); ++i) {
			group *insgroup = i->second;
		
			best_score = std::max(best_score, insgroup->proximity_score);
			allscores *= insgroup->proximity_score;
			VERBOSE(*verb_out << "adding group " << insgroup << " to group map at " << insgroup->bits << std::endl);
			(*groups)[insgroup->bits] = insgroup;
			newgroups.push_back(i->second);
			all_new_groups.push_back(i->second);
		}
		//all_prox_groups.insert(all_prox_groups.end(), prox_groups.begin(), prox_groups.end());
		new_prox_groups = prox_groups;
		prox_groups.clear();
		for (size_t j = 0; j < segments.size(); ++j) {
			bitvec bits(segments.size(), false);
			bits.set(j);
			std::map<bitvec, group *>::iterator k = groups->find(bits);
			if (k != groups->end()) {
				group *grp = k->second;
				update_proximity_groups(new_prox_groups, prox_groups, grp, j);
			}
		}
	}
	
	//if (best_score < 1.0) {
		single_seg_group->proximity_score = 1.0;//-0.75*best_score;//1.0 - (all_new_groups.empty() ? 0.0 : std::pow(allscores, 1.0/all_new_groups.size()));//1.0 - best_score;
		//dilute_group_by_dist(single_seg_group, *groups);
		(*groups)[single_seg_group->bits] = single_seg_group;
		all_new_groups.push_back(single_seg_group);
		newgroups.push_back(single_seg_group);
		VERBOSE(*verb_out << "created default group " << single_seg_group << " at " << single_seg_group->bits << " with score " << single_seg_group->proximity_score << std::endl);		
	//}

	for (std::list<group *>::const_iterator i = all_new_groups.begin(); i != all_new_groups.end(); ++i) {
		group *grp = *i;
		if (grp->container_score == 0) {
			match_score match = prune_and_match_symbol(sqrt, grp);
			if (match.S) {
				/*double csc;
				if (hmatch.proto) {
					csc = std::max(0.0, match.matcher_score-hmatch.matcher_score);
				}
				else {*/
					grp->container_score = 1-match.matcher_score;
				//}
				//grp->container_score = std::max(grp->container_score, csc);
			}
			VERBOSE(*verb_out << "container score for group " << grp->bits << " is " << grp->container_score << std::endl);
		}
	}

	/*std::list<group *> stack_groups;
	update_stack_groups(all_new_groups, stack_groups);
	while (!stack_groups.empty()) {
		newgroups.insert(newgroups.end(), stack_groups.begin(), stack_groups.end());
		for (std::list<group *>::iterator i = stack_groups.begin(); i != stack_groups.end(); ++i) {
			group *insgroup = *i;
			VERBOSE(*verb_out << "adding group " << insgroup << " to group map at " << insgroup->bits << std::endl);
			(*groups)[insgroup->bits] = insgroup;
		}
		all_new_groups = stack_groups;
		stack_groups.clear();
		update_stack_groups(all_new_groups, stack_groups);
	}*/

	clear_dead_groups();
}

void
math_recognizer_base::cleartablesfor(size_t pos) {
#ifdef WIN32
	for (parse_table::iterator i = tab.begin(); i != tab.end();) {
		if (i->first->bits.at(pos)) {
			clear_table(i->second);
			i = tab.erase(i);
		}
		else {
			++i;
		}
	}
#else
	std::list<const ordered_segments *> rmlist;
	for (parse_table::iterator i = tab.begin(); i != tab.end(); ++i) {
		if (i->first->bits.at(pos)) {
			clear_table(i->second);
			rmlist.push_back(i->first);
		}
	}
	while (!rmlist.empty()) {
		tab.erase(tab.find(rmlist.front()));
		rmlist.pop_front();
	}
#endif
}

int
math_recognizer_base::remove_strokes(const size_t *indices, size_t nstrokes)
{
	VERBOSE(*verb_out << "removing " << nstrokes << " strokes\n");
	//VERBOSE(dump_ink_tree(*verb_out, ctx));

	std::vector<size_t> sindices(indices, indices + nstrokes);
	std::sort(sindices.begin(), sindices.end(), std::greater<size_t>());
	for (std::vector<size_t>::const_iterator pi = sindices.begin(); pi != sindices.end(); ++pi) {
		VERBOSE(*verb_out << "removing stroke " << *pi << std::endl);
		if (*pi >= strokes.size()) {
			VERBOSE(*verb_out << " actually that index is too large. Not removing it.\n");
			continue;
		}
		//VERBOSE(*verb_out << "stroke is " << strokes[*pi]->input << std::endl);
		std::vector<segment *>::iterator pseg = std::find(segments.begin(), segments.end(), strokes[*pi]->parent);
		assert(pseg != segments.end());
		remove_segment(pseg - segments.begin());
		//ctx->segments.erase(pseg);
		for (std::vector<std::vector<size_t> >::iterator i = known.begin(); i != known.end(); ++i) {
			std::vector<size_t>::iterator j = std::find(i->begin(), i->end(), *pi);
			if (j != i->end()) {
				known.erase(i);
				break;
			}
		}
	}

	for (size_t i = 0; i < segments.size(); ++i) {
		if (segments[i]->parents.empty()) {
			cleartablesfor(i);
			updategroups(i);
		}
	}
	
	handlenewgroups();
	return 0;
}

int
recognizer_initialize()
{
	default_feature_set.add_feature(new feat_stroke_top, 1.0);
	default_feature_set.add_feature(new feat_stroke_left, 1.0);
	default_feature_set.add_feature(new feat_stroke_width, 1.0);
	default_feature_set.add_feature(new feat_stroke_height, 1.0);
	default_feature_set.add_feature(new feat_stroke_first_x, 1.0);
	default_feature_set.add_feature(new feat_stroke_first_y, 1.0);
	default_feature_set.add_feature(new feat_stroke_last_x, 1.0);
	default_feature_set.add_feature(new feat_stroke_last_y, 1.0);
	default_feature_set.add_feature(new feat_stroke_length, 0.5);
	default_feature_set.prepare_for_matching();

	small_feature_set.add_feature(new feat_aspect_ratio, 1.0);
	small_feature_set.add_feature(new feat_endtoend_dist, 1.0);
	small_feature_set.prepare_for_matching();

	/*
	std::string training_path;
	GetTrainingPath(training_path);
	std::ifstream ncin((training_path+"/ncor.dat").c_str());
	if (!ncin.is_open()) {
		return 0;
	}
	for (;;) {
		std::string nm;
		ncin >> nm;
		if (ncin.eof()) {
			break;
		}
		ncin >> lsp_ncor[nm] >> dtw_ncor[nm] >> feat_ncor[nm];
	}*/


	std::string path;
	GetProfilePath(path);
	std::ifstream bagin((path+"/symfreq.bag").c_str());
	if (bagin.is_open()) {
		symbag = new symbolbag;
		int e = symbag->read(bagin, 0.75);
		if (FAILURE(e)) {
			delete symbag;
			symbag = 0;
			return e;
		}
		VERBOSE(*verb_out << "loaded symbol bag:\n" << *symbag << std::endl);
	}

	return 0;
}


void
recognizer_shutdown()
{
	//elastic_shutdown();
	delete symbag;
	symbag = 0;
}


}
