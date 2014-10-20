#include "relation.h"
#include "links.h"
#include "meas.h"
#include "utils.h"
#include "grammar.h"
#include "expr-node.h"
#include "segment.h"
#include "rect.h"
#include "stream-defs.h"
#include "ordered-segments.h"

#include <fstream>
#include <map>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>


namespace scg
{

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

static std::map<std::pair<rclass_t, rclass_t>, bclassdata> Rdata;
static std::map<std::pair<rclass_t, rclass_t>, bclassdata> ARdata;
static std::map<std::pair<rclass_t, rclass_t>, bclassdata> BRdata;
static std::map<std::pair<rclass_t, rclass_t>, bclassdata> Bdata;
static std::map<std::pair<rclass_t, rclass_t>, cclassdata> Cdata;


static relation *relations[NUM_RELATIONS];

static unsigned rel_orderings[] = { 0, 0, 0, 1, 0 };

static std::string box_name("_BOX");
static std::string symbol_name("_SYMBOL");
static std::string aggregate_name("_aggregate");


static const unsigned NUM_SPECIAL_CLASSES = 3;
const rclass_t INVALID_CLASS = static_cast<rclass_t>(-1);
const rclass_t AGGREGATE_CLASS = 0;
const rclass_t SYMBOL_CLASS = 1;
const rclass_t BOX_CLASS = 2;
const size_t NUM_SPECIAL_RCLASSES = 3;

// these following classes must match the order given in symbols.cc
const rclass_t BASELINE_CLASS = 3;
const rclass_t ASCENDER_CLASS = 4;
const rclass_t DESCENDER_CLASS = 5;
const rclass_t EXTENDER_CLASS = 6;
const rclass_t HALF_ASCENDER_CLASS = 7;
const rclass_t HALF_DESCENDER_CLASS = 8;
const rclass_t HALF_ASCENDER_DESCENDER_CLASS = 9;
const rclass_t CENTERED_CLASS = 10;
const rclass_t FENCE_CLASS = 11;
const rclass_t LARGE_EXTENDER_CLASS = 12;
const rclass_t RANGE_OP_CLASS = 13;
const rclass_t ROOT_CLASS = 14;
const rclass_t HORIZONTAL_CLASS = 15;
const rclass_t PUNCTUATION_CLASS = 16;
const size_t NUM_RCLASSES = 17 - NUM_SPECIAL_RCLASSES;


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

static double
get_centroid_y(const scg::Rect<long> &box, unsigned rclass) {
	if (rclass == scg::BASELINE_CLASS) {
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

static std::pair<double, double>
fromcentroid(const ordered_segments *segs, rclass_t rclass) {
	std::pair<double, double> pt(0,0);
	scg::Rect<long> bounds = segs->bounds();
	pt.first = 0.5 * (bounds.left + bounds.right);
	pt.second = 0.5 * (bounds.top + bounds.bottom);
	return pt;
}

static std::pair<double, double>
tocentroid(const ordered_segments *segs, rclass_t rclass) {
	std::pair<double, double> pt(0,0);
	scg::Rect<long> bounds = segs->bounds();
	pt.first = 0.5 * (bounds.left + bounds.right);
	pt.second = get_centroid_y(bounds, rclass);//0.5 * (bounds.top + bounds.bottom);
	return pt;
}

struct brelation : public relation {
	brelation(const std::string &name_, const std::map<std::pair<rclass_t, rclass_t>, bclassdata> &data);

	double membership(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) const;

	mutable std::map<std::pair<rclass_t, rclass_t>, bclassdata> data;
};

brelation::brelation(const std::string &name_, const std::map<std::pair<rclass_t, rclass_t>, bclassdata> &data_)
	: relation(name_), data(data_)
{
	/*VERBOSE(*verb_out << "relation " << name << " has data\n";
	for (std::map<std::pair<rclass_t, rclass_t>, classdata>::const_iterator i = data.begin(); i != data.end(); ++i) {
		*verb_out << i->first.first << ' ' << i->first.second << ' ' << i->second.tuseful << ' ' << i->second.txmean << ' ' << i->second.tymean << ' ' << i->second.tvar << std::endl;
	});*/
}

double
brelation::membership(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) const
{
	const bclassdata &cd = data[std::make_pair(class1, class2)];
	if (!cd.useful) {
		return -1.0;
	}
	std::pair<double, double> frompt = fromcentroid(from, class1);
	std::pair<double, double> topt = tocentroid(to, class2);
	double dx = topt.first - frompt.first;
	double dy = topt.second - frompt.second;
	double t = std::atan2(dy,dx);
	double dt = std::abs(t - std::atan2(cd.ymean, cd.xmean));
	dt = std::min(dt, 2*M_PI - dt);
	dt = dt/std::sqrt(cd.var);
	//double K = 1.0 / cd.tvar;
	//double m = std::exp(K * std::cos(dt)) / (2*M_PI*bessel0(K));
	//t = 1/(1 + 0.2316419*dt);
	double p = 1/std::sqrt(2*M_PI)*std::exp(-(dt*dt)/2);
	VERBOSE(*verb_out << "membership is " << p << " in " << name << " between " << from->bounds() << " and " << to->bounds() << " with classes " << class1 << " and " << class2 << std::endl);
	return p;
	//double m = 1 - p*t*(0.319381530-t*(0.356563782-t*(1.781477937-t*(1.821255978-t*1.330274429))));
	//m = 2*(1.0 - m);
}


struct crelation : public relation {
	crelation(const std::map<std::pair<rclass_t, rclass_t>, cclassdata> data_) : relation("C"), data(data_) { }

	mutable std::map<std::pair<rclass_t, rclass_t>, cclassdata> data;

	double membership(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) const
	{
		const cclassdata &cd = Cdata[std::make_pair(class1, class2)];
		if (!cd.useful) {
			return -1.0;
		}
		scg::Rect<long> frombox = from->bounds();
		scg::Rect<long> tobox = to->bounds();
		double ol = 1.0 - scg::overlap_proportion(frombox, tobox);
		return 1.0/cd.mean * std::exp(-ol/cd.mean);
		//double m = 1.0 - std::exp(-1.0/cd.mean * ol);
		//m = 1.0 - m;
		//mem = std::max(mem, 1.0/cd.mean * std::exp(-ol/cd.mean));
	}
};

static unsigned
link_name_to_index(const std::string &s)
{
	if (s == "AR") {
		return 0;
	}
	if (s == "R") {
		return 1;
	}
	if (s == "BR") {
		return 2;
	}
	if (s == "B" || s == "BW" || s == "BN") {
		return 3;
	}
	if (s == "C" || s == "Cs") {
		return 4;
	}
	if (s == "nil") {
		return 5;
	}

	THROW_ERROR(E_INVALID, "unknown link name " << s);
}


void
initialize_relations()
{
	std::string path;
	GetTrainingPath(path);
	std::ifstream relin((path + "/reldb").c_str());
	if (!relin.is_open()) {
		THROW_ERROR(E_NOTFOUND, "could not open reldata file");
	}

	for (size_t i = 0; i < NUM_RELATIONS; ++i) {
		std::string name;
		size_t n;
		relin >> name >> n;
		while (n--) {
			rclass_t fromcl, tocl;
			relin >> fromcl >> tocl;
			if (name == "C") {
				cclassdata &cd = Cdata[std::make_pair(fromcl, tocl)];
				relin >> cd.n >> cd.mean >> cd.useful;
			}
			else {
				std::map<std::pair<rclass_t, rclass_t>, bclassdata> *data;
				if (name == "R") {
					data = &Rdata;
				}
				else if (name == "AR") {
					data = &ARdata;
				}
				else if (name == "BR") {
					data = &BRdata;
				}
				else if (name == "B") {
					data = &Bdata;
				}
				bclassdata &cd = (*data)[std::make_pair(fromcl, tocl)];
				relin >> cd.n >> cd.xmean >> cd.ymean >> cd.var >> cd.useful;
			}
		}
		CHECK_ISTREAM(relin, "error reading relation data for " << name);
	}
	relations[0] = new brelation("AR", ARdata);
	relations[1] = new brelation("R", Rdata);
	relations[2] = new brelation("BR", BRdata);
	relations[3] = new brelation("B", Bdata);
	relations[4] = new crelation(Cdata);
}


void
destroy_relations()
{
	for (unsigned i = 0; i < NUM_RELATIONS; i++) {
		delete relations[i];
	}
}


unsigned
ordering_for_relation(const relation *rel) {
	relation **rp = std::find(relations, relations + NUM_RELATIONS, rel);
	assert(rp != relations + NUM_RELATIONS);
	return rel_orderings[rp - relations];
}


relation *
get_relation(int rel) {
	return relations[rel];
}


}
