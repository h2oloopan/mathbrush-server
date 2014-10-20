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
#include "distrib.h"

#include <fstream>
#include <map>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>


namespace scg
{

struct classdata {
	distrib *D[4];
	classdata() {
		D[0] = 0;
		D[1] = 0;
		D[2] = 0;
		D[3] = 0;
	}
	~classdata() {
		delete D[0];
		delete D[1];
		delete D[2];
		delete D[3];
	}
	bool useful() const {
		return (D[0] && D[1] && D[2] && D[3]) && D[0]->useful() && D[1]->useful() && D[2]->useful() && D[3]->useful();
	}
};

struct cclassdata {
	exp_dist D;
	cclassdata() { }
	bool useful() const { return D.useful(); }
};

static std::map<std::pair<rclass_t, rclass_t>, classdata> Rdata;
static std::map<std::pair<rclass_t, rclass_t>, classdata> ARdata;
static std::map<std::pair<rclass_t, rclass_t>, classdata> BRdata;
static std::map<std::pair<rclass_t, rclass_t>, classdata> Bdata;
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

rclass_t
tag_to_rclass(const std::string &tag) {
	if (tag == "_aggregate") {
		return AGGREGATE_CLASS;
	}
	else if (tag == "_box") {
		return BOX_CLASS;
	}
	else if (tag == "_symbol") {
		return SYMBOL_CLASS;
	}
	else if (tag == "baseline") {
		return BASELINE_CLASS;
	}
	else if (tag == "ascender") {
		return ASCENDER_CLASS;
	}
	else if (tag == "descender") {
		return DESCENDER_CLASS;
	}
	else if (tag == "extender") {
		return EXTENDER_CLASS;
	}
	else if (tag == "half-ascender") {
		return HALF_ASCENDER_CLASS;
	}
	else if (tag == "half-descender") {
		return HALF_DESCENDER_CLASS;
	}
	else if (tag == "half-ascender-descender") {
		return HALF_ASCENDER_DESCENDER_CLASS;
	}
	else if (tag == "centered") {
		return CENTERED_CLASS;
	}
	else if (tag == "fence") {
		return FENCE_CLASS;
	}
	else if (tag == "large-extender") {
		return LARGE_EXTENDER_CLASS;
	}
	else if (tag == "range-op") {
		return RANGE_OP_CLASS;
	}
	else if (tag == "root") {
		return ROOT_CLASS;
	}
	else if (tag == "horizontal") {
		return HORIZONTAL_CLASS;
	}
	else if (tag == "punctuation") {
		return PUNCTUATION_CLASS;
	}
	return (rclass_t)-1;
}

std::string
rclass_to_tag(rclass_t cls) {
	switch (cls) {
	case AGGREGATE_CLASS:
	case BOX_CLASS:
	case SYMBOL_CLASS:
		return "";
	case BASELINE_CLASS:
		return "baseline";
	case ASCENDER_CLASS:
		return "ascender";
	case DESCENDER_CLASS:
		return "descender";
	case EXTENDER_CLASS:
		return "extender";
	case HALF_ASCENDER_CLASS:
		return "half-ascender";
	case HALF_DESCENDER_CLASS:
		return "half-descender";
	case HALF_ASCENDER_DESCENDER_CLASS:
		return "half-ascender-descender";
	case CENTERED_CLASS:
		return "centered";
	case FENCE_CLASS:
		return "fence";
	case LARGE_EXTENDER_CLASS:
		return "large-extender";
	case RANGE_OP_CLASS:
		return "range-op";
	case ROOT_CLASS:
		return "root";
	case HORIZONTAL_CLASS:
		return "horizontal";
	case PUNCTUATION_CLASS:
		return "punctuation";
	default:
		return "?rclass?";
	}
}


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

static double
inset(long x) {
	return 0.25 * std::max((double)x, 2540.0/3.0);
}

static std::pair<double, double>
fromcentroid(const ordered_segments *segs, rclass_t rclass, const relation *rel) {
	std::pair<double, double> pt(0,0);
	scg::Rect<long> bounds = segs->bounds();
	if (rel == relations[REL_BELOW]) {
		pt.first = 0.5 * (bounds.left + bounds.right);
	}
	else {
		pt.first = bounds.right - inset(width(bounds));
	}
	pt.second = 0.5 * (bounds.top + bounds.bottom);
	return pt;
/*	std::pair<double, double> pt(0,0);
	scg::Rect<long> bounds = segs->bounds();
	pt.first = 0.5 * (bounds.left + bounds.right);
	pt.second = 0.5 * (bounds.top + bounds.bottom);
	return pt;*/
}

static std::pair<double, double>
tocentroid(const ordered_segments *segs, rclass_t rclass, const relation *rel) {
	std::pair<double, double> pt(0,0);
	scg::Rect<long> bounds = segs->bounds();
	if (rel == relations[REL_BELOW]) {
		pt.first = 0.5 * (bounds.left + bounds.right);
	}
	else {
		pt.first = bounds.left + inset(width(bounds));
	}
	//pt.second = get_centroid_y(bounds, rclass);//0.5 * (bounds.top + bounds.bottom);
	pt.second = 0.5 * (bounds.top + bounds.bottom);
	return pt;
/*	std::pair<double, double> pt(0,0);
	scg::Rect<long> bounds = segs->bounds();
	pt.first = 0.5 * (bounds.left + bounds.right);
	pt.second = get_centroid_y(bounds, rclass);//0.5 * (bounds.top + bounds.bottom);
	return pt;*/
}

struct defrelation : public relation {
	defrelation(const std::string &name_, const std::map<std::pair<rclass_t, rclass_t>, classdata> &data);

	double membership(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) const;
	void add_sample(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2);
	void write(std::ostream &os) const;

	mutable std::map<std::pair<rclass_t, rclass_t>, classdata> data;
};

defrelation::defrelation(const std::string &name_, const std::map<std::pair<rclass_t, rclass_t>, classdata> &data_)
	: relation(name_), data(data_)
{
	/*VERBOSE(*verb_out << "relation " << name << " has data\n";
	for (std::map<std::pair<rclass_t, rclass_t>, classdata>::const_iterator i = data.begin(); i != data.end(); ++i) {
		*verb_out << i->first.first << ' ' << i->first.second << ' ' << i->second.tuseful << ' ' << i->second.txmean << ' ' << i->second.tymean << ' ' << i->second.tvar << std::endl;
	});*/
}

double
defrelation::membership(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) const
{
	const classdata &cd = data[std::make_pair(class1, class2)];
	// aggregate data is always useful since it is used as the last resort
	if (!cd.useful() && (class1 != AGGREGATE_CLASS || class2 != AGGREGATE_CLASS)) {
		VERBOSE(*verb_out << "membership is not useful in " << name << " between " << from->bounds() << " and " << to->bounds() << " with classes " << class1 << " and " << class2 << std::endl);
		return -1.0;
	}
	std::pair<double, double> frompt = fromcentroid(from, class1, this);
	std::pair<double, double> topt = tocentroid(to, class2, this);
	double dx = topt.first - frompt.first;
	double dy = topt.second - frompt.second;
	dx /= TABLETPC_DPI;
	dy /= TABLETPC_DPI;
	/*if (name == "AR") {
		dy *= -1;
	}*/
	double p = std::sqrt(cd.D[2]->score(std::atan2(dy,dx)) * cd.D[3]->score(std::sqrt(dx*dx+dy*dy)));//cd.D[0]->score(dx) * cd.D[1]->score(dy);
	VERBOSE(*verb_out << "membership is " << p << " in " << name << " between " << from->bounds() << " and " << to->bounds() << " with classes " << class1 << " and " << class2 << std::endl);
	return p;
}

void
defrelation::add_sample(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) {
	classdata &cd = data[std::make_pair(class1, class2)];
	std::pair<double, double> frompt = fromcentroid(from, class1, this);
	std::pair<double, double> topt = tocentroid(to, class2, this);
	double dx = topt.first - frompt.first;
	double dy = topt.second - frompt.second;
	dx /= TABLETPC_DPI;
	dy /= TABLETPC_DPI;
	/*if (name == "AR") {
		dy *= -1;
	}*/
	cd.D[0]->update(dx);
	cd.D[1]->update(dy);
	cd.D[2]->update(std::atan2(dy,dx));
	cd.D[3]->update(std::sqrt(dx*dx+dy*dy));
}

void
defrelation::write(std::ostream &os) const {
	size_t n = 0;
	for (std::map<std::pair<rclass_t, rclass_t>, classdata>::const_iterator i = data.begin(); i != data.end(); ++i) {
		const classdata &cd = i->second;
		if (cd.D[0] && cd.D[1] && cd.D[2] && cd.D[3]) {
			++n;
		}
	}
	os << name << ' ' << n << std::endl;
	for (std::map<std::pair<rclass_t, rclass_t>, classdata>::const_iterator i = data.begin(); i != data.end(); ++i) {
		const classdata &cd = i->second;
		if (cd.D[0] && cd.D[1] && cd.D[2] && cd.D[3]) {
			os << i->first.first << ' ' << i->first.second << ' ' << *cd.D[0] << ' ' << *cd.D[1] << ' ' << *cd.D[2] << ' ' << *cd.D[3] << std::endl;
		}
	}
}


struct crelation : public relation {
	crelation(const std::map<std::pair<rclass_t, rclass_t>, cclassdata> data_) : relation("C"), data(data_) { }

	mutable std::map<std::pair<rclass_t, rclass_t>, cclassdata> data;

	double membership(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) const
	{
		const cclassdata &cd = Cdata[std::make_pair(class1, class2)];
		// aggregate data is always useful since it is used as the last resort
		if (!cd.useful() && (class1 != AGGREGATE_CLASS || class2 != AGGREGATE_CLASS)) {
			return -1.0;
		}
		scg::Rect<long> frombox = from->bounds();
		scg::Rect<long> tobox = to->bounds();
		double ol = 1.0 - scg::overlap_proportion(frombox, tobox);
		return cd.D.score(ol);
	}

	void add_sample(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) {
		cclassdata &cd = Cdata[std::make_pair(class1, class2)];
		scg::Rect<long> frombox = from->bounds();
		scg::Rect<long> tobox = to->bounds();
		double ol = 1.0 - scg::overlap_proportion(frombox, tobox);
		cd.D.update(ol);
	}

	void write(std::ostream &os) const {
		os << "C " << data.size() << std::endl;
		for (std::map<std::pair<rclass_t, rclass_t>, cclassdata>::const_iterator i = data.begin(); i != data.end(); ++i) {
			os << i->first.first << ' ' << i->first.second << ' ' << i->second.D << std::endl;
		}
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
				relin >> cd.D;
			}
			else {
				std::map<std::pair<rclass_t, rclass_t>, classdata> *data;
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
				classdata &cd = (*data)[std::make_pair(fromcl, tocl)];
				if (name == "R") {
					cd.D[0] = new normal_dist;
					cd.D[1] = new normal_dist;
					cd.D[2] = new vonmises_dist;
					cd.D[3] = new gamma_dist;
				}
				else if (name == "AR") {
					cd.D[0] = new normal_dist;
					cd.D[1] = new normal_dist;
					cd.D[2] = new vonmises_dist;
					cd.D[3] = new gamma_dist;
				}
				else if (name == "BR") {
					cd.D[0] = new normal_dist;
					cd.D[1] = new normal_dist;
					cd.D[2] = new vonmises_dist;
					cd.D[3] = new gamma_dist;
				}
				else if (name == "B") {
					cd.D[0] = new normal_dist;
					cd.D[1] = new normal_dist;
					cd.D[2] = new vonmises_dist;
					cd.D[3] = new gamma_dist;
				}
				VERBOSE(*verb_out << "reading relation data for " << name << std::endl);
				relin >> *cd.D[0] >> *cd.D[1] >> *cd.D[2] >> *cd.D[3];
			}
		}
		CHECK_ISTREAM(relin, "error reading relation data for " << name);
	}
	relations[0] = new defrelation("AR", ARdata);
	relations[1] = new defrelation("R", Rdata);
	relations[2] = new defrelation("BR", BRdata);
	relations[3] = new defrelation("B", Bdata);
	relations[4] = new crelation(Cdata);
}


void
destroy_relations()
{
	std::string path;
	GetTrainingPath(path);
	std::ofstream relout((path + "/reldb").c_str());
	if (!relout.is_open()) {
		THROW_ERROR(E_NOTFOUND, "could not open reldata file");
	}

	for (unsigned i = 0; i < NUM_RELATIONS; i++) {
		relations[i]->write(relout);
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
