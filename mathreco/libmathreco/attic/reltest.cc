#include "annotate.h"
#include "symbols.h"
#include "error.h"
#include "group.h"
#include "stroke.h"
#include "relation.h"
#include "MathRecognizer.h"
#include "mathrecognizer-private.h"
#include "utils.h"
#include "relation.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


static double
bessel0(double x) {
	double t = x/3.75;
	if (std::abs(x) <= 3.75) {
		double t2 = t*t;
		return 1 + t2*(3.5156229 + t2*(3.0899424 + t2*(1.2067492 + t2*(0.2659732 + t2*(0.0360758 + t2*0.0045813)))));
	}
	t = 3.75/x;
	return std::exp(x)/std::sqrt(x) * (0.39894228 + t*(0.01328592+t*(0.00225319-t*(0.00157565-t*(0.00916281-t*(0.02057706-t*(0.02635537-t*(0.01647633-t*0.00392377))))))));
}

struct classdata {
	size_t n;
	double dmu;
	double dlambda;
	double txmean;
	double tymean;
	double tvar;
	bool duseful;
	bool tuseful;
	classdata() : dmu(0.0), dlambda(0.0), txmean(0.0), tymean(0.0), tvar(0.0), duseful(false), tuseful(false) { }
};

struct blclassdata {
	size_t n;
	double mean;
	double var;
	bool useful;
	double mean2;
	double var2;
	bool useful2;
	blclassdata() : n(0), mean(0.0), var(0.0), useful(false), mean2(0), var2(0), useful2(false) { }
};

struct bclassdata {
	size_t n;
	double xmean;
	double ymean;
	double var;
	double dxmean;
	double dxvar;
	bool useful;
	bool dxuseful;
	bclassdata() : n(0), xmean(0.0), ymean(0.0), var(0.0), useful(false), dxmean(0), dxvar(0), dxuseful(false) { }
};

struct cclassdata {
	size_t n;
	double mean;
	bool useful;
	cclassdata() : n(0), mean(0.0), useful(false) { }
};

std::map<std::pair<unsigned, unsigned>, bclassdata> Rdata;
std::map<std::pair<unsigned, unsigned>, bclassdata> ARdata;
std::map<std::pair<unsigned, unsigned>, bclassdata> BRdata;
std::map<std::pair<unsigned, unsigned>, bclassdata> Bdata;
std::map<std::pair<unsigned, unsigned>, cclassdata> Cdata;

const scg::symbols_db *db;

static void
find_classes(const scg::AnnotatedStrokeGroup &G, const scg::stroke_collection &stks, std::set<unsigned> &cls) {
	for (scg::AnnotatedStrokeGroup::const_symbol_iterator j = G.symbols_begin(); j != G.symbols_end(); ++j) {
		const scg::symbol_annotation &sa = *j;
		if (sa.strokes == stks) {
			const scg::symbol *S = db->find_symbol(sa.name);
			if (!S) {
				cls.insert(scg::INVALID_CLASS);
			}
			else {
				cls.insert(S->info.classes.begin(), S->info.classes.end());
			}
			break;
		}
	}
	if (cls.empty()) {
		cls.insert(scg::BOX_CLASS);
	}
	cls.insert(scg::AGGREGATE_CLASS);
}

static double
get_centroid_y(const scg::Rect<long> &box, unsigned rclass) {
	/*
	if (rclass == scg::AGGREGATE_CLASS
	 || rclass == scg::SYMBOL_CLASS
	 || rclass == scg::BOX_CLASS
	 || rclass == scg::CENTERED_CLASS
	 || rclass == scg::LARGE_EXTENDER_CLASS
	 || rclass == scg::RANGE_OP_CLASS
	 || rclass == scg::FENCE_CLASS
	 || rclass == scg::ROOT_CLASS
	 || rclass == scg::HORIZONTAL_CLASS
	 || rclass == scg::PUNCTUATION_CLASS
	 || rclass == scg::ASCENDER_CLASS
	 || rclass == scg::EXTENDER_CLASS) {
		return (box.top + box.bottom) / 2.0;
	}
	else*/ if (rclass == scg::BASELINE_CLASS) {
		return box.top + height(box) * 0.1;
	}
	else if (rclass == scg::HALF_DESCENDER_CLASS) {
		return box.top + height(box) * 0.1;
	}
	else if (rclass == scg::HALF_ASCENDER_DESCENDER_CLASS) {
		return box.top + height(box) * 0.25;
	}
	else if (rclass == scg::DESCENDER_CLASS) {
		return box.top + height(box) * 0.05;
	}
	else if (rclass == scg::HALF_ASCENDER_CLASS) {
		return box.top + height(box) * 0.33;
	}
	else {
		return (box.top + box.bottom) / 2.0;
	}
}

static scg::Rect<long>
make_bounds(const scg::RawStrokeGroup &G, const scg::stroke_collection &stks) {
	scg::stroke_collection::const_iterator i = stks.begin();
	scg::Rect<long> bounds = scg::bbox(G.strokes[*i]);
	for (++i; i != stks.end(); ++i) {
		bounds = scg::merge(bounds, scg::bbox(G.strokes[*i]));
	}
	return bounds;
}

static std::string
standardize(const std::string &s) {
	if (s == "BN" || s == "BW") {
		return "B";
	}
	if (s == "Cs") {
		return "C";
	}
	return s;
}

static std::pair<double, double>
fromcentroid(const scg::AnnotatedStrokeGroup &G, const scg::stroke_collection &stks, const std::string &rel, unsigned rclass) {
	std::pair<double, double> pt(0,0);
	scg::Rect<long> bounds = make_bounds(G, stks);
	pt.first = 0.5 * (bounds.left + bounds.right);
	pt.second = 0.5 * (bounds.top + bounds.bottom);
	return pt;

	if (rel == "R" || rel == "AR" || rel == "BR" || rel == "C") {
		/*if (rel == "R") {
			pt.second = 0.5 * (bounds.top + bounds.bottom);
		}
		if (rel == "AR") {
			pt.second = bounds.bottom;
		}
		if (rel == "BR") {
			pt.second = bounds.top;
		}

		pt.first = bounds.right;
		if (rel == "C") {
			pt.first = bounds.right;
			pt.second = bounds.top;
		}*/
		/*else*/ {
			for (scg::stroke_collection::const_iterator i = stks.begin(); i != stks.end(); ++i) {
				const scg::RawStroke &stk = G.strokes[*i];
				bounds = scg::bbox(stk);
				pt.first = std::max(pt.first, (double)bounds.left);
			}
		}
	}
	else {
		pt.first = 0.5 * (bounds.left + bounds.right);
		//pt.second = 0.5 * (bounds.top + bounds.bottom);
		for (scg::stroke_collection::const_iterator i = stks.begin(); i != stks.end(); ++i) {
			const scg::RawStroke &stk = G.strokes[*i];
			bounds = scg::bbox(stk);
			pt.second = std::max(pt.second, (double)bounds.top);
		}
	}
	return pt;
}

static std::pair<double, double>
tocentroid(const scg::AnnotatedStrokeGroup &G, const scg::stroke_collection &stks, const std::string &rel, unsigned rclass) {
	std::pair<double, double> pt(0,0);
	scg::Rect<long> bounds = make_bounds(G, stks);
	pt.first = 0.5 * (bounds.left + bounds.right);
	pt.second = get_centroid_y(bounds, rclass);//0.5 * (bounds.top + bounds.bottom);
	return pt;
	if (rel == "R" || rel == "AR" || rel == "BR" || rel == "C") {
		if (rel == "R") {
			pt.second = 0.5 * (bounds.top + bounds.bottom);
		}
		if (rel == "AR") {
			pt.second = bounds.bottom;
		}
		if (rel == "BR") {
			pt.second = bounds.top;
		}
		if (rel == "C") {
			pt.second = bounds.top;
		}
		pt.first = bounds.left;
		pt.second = get_centroid_y(bounds, rclass);//0.5 * (bounds.top + bounds.bottom);
		//pt.second = 0.5 * (bounds.top + bounds.bottom);
	}
	else {
		pt.first = 0.5 * (bounds.left + bounds.right);
		//pt.second = 0.5 * (bounds.top + bounds.bottom);
		pt.second = bounds.top;
	}
	return pt;
}

static unsigned globaltests = 0, globalcorr = 0;

static double
blmeas(const scg::AnnotatedStrokeGroup &G, const scg::stroke_collection &from, const scg::stroke_collection &to, const std::set<unsigned> &from_classes, const std::set<unsigned> &to_classes, const std::string &rel) {
	std::map<std::pair<unsigned, unsigned>, blclassdata> *data;
	/*if (rel == "R") data = &Rdata;
	if (rel == "AR") data = &ARdata;
	if (rel == "BR") data = &BRdata;*/
	double mem = -1.0;
	for (std::set<unsigned>::reverse_iterator i1 = from_classes.rbegin(); i1 != from_classes.rend(); ++i1) {
		if (*i1 < scg::NUM_RCLASSES && mem != -1.0) {
			break;
		}
		for (std::set<unsigned>::reverse_iterator i2 = to_classes.rbegin(); i2 != to_classes.rend(); ++i2) {
			if (*i2 < scg::NUM_RCLASSES && mem != -1.0) {
				break;
			}
			const blclassdata &cd = (*data)[std::make_pair(*i1, *i2)];
			double loc = -1.0;
			if (cd.useful) {
				scg::Rect<long> frombox = make_bounds(G, from);
				scg::Rect<long> tobox = make_bounds(G, to);
				double d = std::abs(d-cd.mean)/std::sqrt(cd.var);
				double t = 1/(1 + 0.2316419*d);
				double p = 1/std::sqrt(2*M_PI)*std::exp(-(d*d)/2);
				double m = 1 - p*t*(0.319381530-t*(0.356563782-t*(1.781477937-t*(1.821255978-t*1.330274429))));
				m = 2*(1.0 - m);
				loc = m;
			}
			if (cd.useful2) {
				scg::Rect<long> frombox = make_bounds(G, from);
				scg::Rect<long> tobox = make_bounds(G, to);
				double d = (tobox.bottom - frombox.bottom) / 2540.0;
				d = std::abs(d-cd.mean2)/std::sqrt(cd.var2);
				double t = 1/(1 + 0.2316419*d);
				double p = 1/std::sqrt(2*M_PI)*std::exp(-(d*d)/2);
				double m = 1 - p*t*(0.319381530-t*(0.356563782-t*(1.781477937-t*(1.821255978-t*1.330274429))));
				m = 2*(1.0 - m);
				loc = (loc == -1.0) ? m : std::sqrt(loc*m);
			}
			if (loc > mem) {
				std::cout << "  taking score " << loc << " for " << rel << " with " << *i1 << ',' << *i2 << std::endl;
			}
			mem = std::max(mem, loc);
		}
	}
	return mem;
}

static double
cmeas(const scg::AnnotatedStrokeGroup &G, const scg::stroke_collection &from, const scg::stroke_collection &to, const std::set<unsigned> &from_classes, const std::set<unsigned> &to_classes) {
	double mem = -1.0;
	for (std::set<unsigned>::reverse_iterator i1 = from_classes.rbegin(); i1 != from_classes.rend(); ++i1) {
		for (std::set<unsigned>::reverse_iterator i2 = to_classes.rbegin(); i2 != to_classes.rend(); ++i2) {
			const cclassdata &cd = Cdata[std::make_pair(*i1, *i2)];;;;
			if (cd.useful) {
				scg::Rect<long> frombox = make_bounds(G, from);
				scg::Rect<long> tobox = make_bounds(G, to);
				double ol = 1.0 - scg::overlap_proportion(frombox, tobox);
				double m = 1.0 - std::exp(-1.0/cd.mean * ol);
				m = 1.0 - m;
				mem = std::max(mem, 1.0/cd.mean * std::exp(-ol/cd.mean));
			}
			if (*i2 < scg::NUM_RCLASSES && mem != -1.0) {
				break;
			}
		}
		if (*i1 < scg::NUM_RCLASSES && mem != -1.0) {
			break;
		}
	}
	return mem;
}

static double
bmeas(const scg::AnnotatedStrokeGroup &G, const scg::stroke_collection &from, const scg::stroke_collection &to, const std::set<unsigned> &from_classes, const std::set<unsigned> &to_classes, const std::string &rel) {
	std::map<std::pair<unsigned, unsigned>, bclassdata> *data;
	if (rel == "R") data = &Rdata;
	if (rel == "AR") data = &ARdata;
	if (rel == "BR") data = &BRdata;
	if (rel == "B") data = &Bdata;

	double mem = -1.0;
	for (std::set<unsigned>::reverse_iterator i1 = from_classes.rbegin(); i1 != from_classes.rend(); ++i1) {
		for (std::set<unsigned>::reverse_iterator i2 = to_classes.rbegin(); i2 != to_classes.rend(); ++i2) {
			const bclassdata &cd = (*data)[std::make_pair(*i1, *i2)];
			if (cd.useful) {
				std::pair<double, double> frompt = fromcentroid(G, from, rel, *i1);
				std::pair<double, double> topt = tocentroid(G, to, rel, *i2);
				double dx = topt.first - frompt.first;
				double dy = topt.second - frompt.second;
				double t = std::atan2(dy,dx);
				double dt = std::abs(t - std::atan2(cd.ymean, cd.xmean));
				dt = std::min(dt, 2*M_PI - dt);
				dt = dt/std::sqrt(cd.var);
				//double K = 1.0 / cd.tvar;
				//double m = std::exp(K * std::cos(dt)) / (2*M_PI*bessel0(K));
				t = 1/(1 + 0.2316419*dt);
				double p = 1/std::sqrt(2*M_PI)*std::exp(-(dt*dt)/2);
				double m = 1 - p*t*(0.319381530-t*(0.356563782-t*(1.781477937-t*(1.821255978-t*1.330274429))));
				m = 2*(1.0 - m);
				//std::cout << "  (for " << *i1 << "," << *i2 << ", t=" << t << ", mt=" << std::atan2(cd.ymean, cd.xmean) << ", dt=" << dt << ", p=" << p << ',' << m << ")\n";
				//std::cout << "   (score for " << *i1 << "," << *i2 << " is " << m << std::endl;

				if (cd.dxuseful) {
					scg::Rect<long> frombox = make_bounds(G, from);
					scg::Rect<long> tobox = make_bounds(G, to);
					dx = (tobox.left - frombox.right) / 2540.0;
					dx = (dx - cd.dxmean) / std::sqrt(cd.dxvar);
					t = 1/(1 + 0.2316419*dx);
					double p2 = 1/std::sqrt(2*M_PI)*std::exp(-(dx*dx)/2);
					double m2 = 1 - p*t*(0.319381530-t*(0.356563782-t*(1.781477937-t*(1.821255978-t*1.330274429))));
					m2 = 2*(1.0 - m2);
					m = std::sqrt(m*m2);
					std::cout << "p: " << p <<  " -> ";
					p = std::sqrt(p*p2);
					std::cout << p << std::endl;
				}
				mem = std::max(mem, p);
			}
			if (*i2 < scg::NUM_RCLASSES && mem != -1.0) {
				break;
			}
		}
		if (*i1 < scg::NUM_RCLASSES && mem != -1.0) {
			break;
		}
	}
	return mem;
}

static void
test_ink(const scg::AnnotatedStrokeGroup &G) {
	unsigned numtests = 0, numcorrect = 0;
	for (scg::AnnotatedStrokeGroup::const_link_iterator i = G.links_begin(); i != G.links_end(); ++i) {
		const scg::link_annotation &annot = *i;
		if (annot.from.empty() || annot.to.empty()) continue;
		std::set<unsigned> from_classes, to_classes;
		std::string relname = standardize(annot.rel);
		find_classes(G, annot.from, from_classes);
		find_classes(G, annot.to, to_classes);
		std::cout << "Link is " << relname << " between " << make_bounds(G, annot.from) << " and " << make_bounds(G, annot.to) << " with classes ";
		for (std::set<unsigned>::const_iterator q = from_classes.begin(); q != from_classes.end(); ++q) {
			std::cout << *q << ',';
		}
		std::cout << " and ";
		for (std::set<unsigned>::const_iterator q = to_classes.begin(); q != to_classes.end(); ++q) {
			std::cout << *q << ',';
		}
		std::cout << std::endl;
		std::map<double, std::string, std::greater<double> > results;
		double c = cmeas(G, annot.from, annot.to, from_classes, to_classes);
		double b = bmeas(G, annot.from, annot.to, from_classes, to_classes, "B");
		double r = bmeas(G, annot.from, annot.to, from_classes, to_classes, "R");
		double ar = bmeas(G, annot.from, annot.to, from_classes, to_classes, "AR");
		double br = bmeas(G, annot.from, annot.to, from_classes, to_classes, "BR");
		results[b] = "B";
		results[c] = "C";
		results[r] = "R";
		results[ar] = "AR";
		results[br] = "BR";

		std::string bestrel = results.begin()->second;
		for (std::map<double, std::string, std::greater<double> >::iterator j = results.begin(); j != results.end(); ++j) {
			std::cout << "  " << j->second << " = " << j->first << std::endl;
		}
		++numtests;
		if (bestrel == relname) {
			++numcorrect;
		}
		else {
			std::cout << "INCORRECT\n";
		}
	}
	std::cout << "On this file, " << numcorrect << " correct of " << numtests << std::endl;
	globaltests += numtests;
	globalcorr += numcorrect;
}

static void
loadbl(std::istream &in, const std::string &rel) {
	std::map<std::pair<unsigned, unsigned>, blclassdata> *data;
	/*if (rel == "R") data = &Rdata;
	if (rel == "AR") data = &ARdata;
	if (rel == "BR") data = &BRdata;*/
	unsigned fromcl, tocl;
	in >> fromcl >> tocl;
	blclassdata &cd = (*data)[std::make_pair(fromcl, tocl)];
	in >> cd.n >> cd.mean >> cd.var >> cd.useful;
	in >> cd.mean2 >> cd.var2 >> cd.useful2;
}

static void
loadb(std::istream &in, const std::string &rel) {
	std::map<std::pair<unsigned, unsigned>, bclassdata> *data;
	if (rel == "R") data = &Rdata;
	if (rel == "AR") data = &ARdata;
	if (rel == "BR") data = &BRdata;
	if (rel == "B") data = &Bdata;
	unsigned fromcl, tocl;
	in >> fromcl >> tocl;
	bclassdata &cd = (*data)[std::make_pair(fromcl, tocl)];
	in >> cd.n >> cd.xmean >> cd.ymean >> cd.var >> cd.useful;
	in >> cd.dxmean >> cd.dxvar >> cd.dxuseful;
}

static void
loadc(std::istream &in) {
	unsigned fromcl, tocl;
	in >> fromcl >> tocl;
	cclassdata &cd = Cdata[std::make_pair(fromcl, tocl)];
	in >> cd.n >> cd.mean >> cd.useful;
}

int
main(int argc, char *argv[])
{
	scg::RecognizerHandle rh = scg::InitializeRecognizer();
	db = scg::global_symbols_db();

	std::string path;
	scg::GetTrainingPath(path);
	std::ifstream relin((path + "/reldb").c_str());
	if (!relin.is_open()) {
		THROW_ERROR(E_NOTFOUND, "could not open reldata file");
	}

	for (size_t i = 0; i < scg::NUM_RELATIONS; ++i) {
		std::string name;
		size_t n;
		relin >> name >> n;
		while (n--) {
			if (name == "R" || name == "AR" || name == "BR") {
				loadb(relin, name);
			}
			if (name == "B") {
				loadb(relin, "B");
			}
			if (name == "C") {
				loadc(relin);
			}
		}
	}
	while (std::cin) {
		std::string line;
		std::getline(std::cin, line);
		std::ifstream inkfile(line.c_str());
		if (inkfile.is_open() && inkfile) {
			scg::AnnotatedStrokeGroup G;
			int e = scg::import_annotated_ink(inkfile, G);
			if (!FAILURE(e)) {
				test_ink(G);
			}
		}
	}
	std::cout << "On this set, " << globalcorr << " correct of " << globaltests << "(" << 100.0*globalcorr/globaltests << "%)\n";
	scg::ShutdownRecognizer(rh);
	return 0;
}
