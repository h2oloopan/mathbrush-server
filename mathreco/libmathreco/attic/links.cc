#include "links.h"
#include "meas.h"
#include "memory.h"
#include "verb.h"
#include "error.h"

#include <cmath>
#include <iostream>
#include <complex>


#ifdef WIN32
#define ISNAN _isnan
#else
#define ISNAN std::isnan
#endif



namespace scg
{


static double
compute_membership(double x, double mean, double var)
{
	static double n = 10.0;
	static double log2 = std::log(2.0);
	static std::complex<double> negone(-1.0, 0.0);
	static std::complex<double> c(0, M_PI);

	std::complex<double> split(n * var / (1 + n * var));
	std::complex<double> z((x - mean) / (1 + std::abs(x - mean)), 0.0);
	std::complex<double> a = log2 / (c - std::log(split + 1.0) - std::log(split - 1.0));
	double l = 1.0 / std::real(std::pow(negone, a));
	std::complex<double> p(std::pow(z + 1.0, a) * std::pow(z - 1.0, a));
	return l * std::real(p);
}



int
BoxLinkEstimator::read_statistics(std::istream &is)
{
	unsigned D = link_codes.size();
	unsigned M = measurer_indices.size();

	scratch.clear();
	scratch.insert(scratch.end(), M, 0.0);

	stats.clear();
	stats.insert(stats.end(), D, std::vector<Statistics>());
	for (unsigned i = 0; i < D; i++) {
		stats[i].insert(stats[i].end(), M, Statistics());
	}

	Statistics dummy;

	std::string link;
	is >> link;
	for (unsigned i = 0; i < D; i++) {
		std::map<std::string, unsigned>::const_iterator lci = link_codes.find(link);
		unsigned li;
		if (lci == link_codes.end()) {
			std::cerr << "error cherching link " << link << std::endl;
			li = (unsigned)-1;
		}
		else {
			li = lci->second;
		}
		
		std::string name;
		for (;;) {
			is >> name;
			if (is.eof()) break;
			if (link_codes.find(name) != link_codes.end()) {
				link = name;
				break;
			}

			unsigned j;
			lci = measurer_indices.find(name);
			if (lci == measurer_indices.end()) {
				std::cerr << "error cherching name " << name << std::endl;
				j = (unsigned)-1;
			}
			else {
				j = lci->second;
			}

			Statistics *S;
			if (li == (unsigned)-1 || j == (unsigned)-1) {
				S = &dummy;
			}
			else {
				S = &stats[li][j];
			}

			is >> S->N >> S->mean >> S->sum_of_squares >> S->var;
			S->var = std::sqrt(S->var);
		}
	}
	
	return 0;
}

int
BoxLinkEstimator::write_statistics(std::ostream &os) const
{
	for (std::map<std::string, unsigned>::const_iterator link = link_codes.begin(); link != link_codes.end(); ++link) {
		os << link->first << std::endl;
		unsigned li = link->second;
		for (std::map<std::string, unsigned>::const_iterator meas = measurer_indices.begin(); meas != measurer_indices.end(); ++meas) {
			os << meas->first << " ";
			unsigned mi = meas->second;
			const Statistics &s = stats[li][mi];
			os << s.N << ' ' << s.mean << ' ' << s.sum_of_squares << ' ' << s.var * s.var << std::endl;
		}
	}

	return 0;
}


double
BoxLinkEstimator::combine_scratch_scores(unsigned rel) const
{
	const std::vector<unsigned> &sel_meas = selected_measurers[rel];
	std::vector<unsigned>::const_iterator i = sel_meas.begin();
	if (i == sel_meas.end()) {
		return 0.0;
	}

	const static double p = 5.0;
	double r = scratch[*(i++)];
	for (; i != sel_meas.end(); ++i) {
		//Yager's intersection operator
		r *= scratch[*i];
		//r = 1.0 - std::min(1.0, std::pow(std::pow(1.0 - r, p) + std::pow(1.0 - scratch[*i], p), 1.0/p));//(D[j] * membership[i][j]) / (r + (1 - r) * (D[j] + membership[i][j] - D[j] * membership[i][j]));
	}

	return std::pow(r, 1.0 / sel_meas.size());
}


double
BoxLinkEstimator::estimate_link(const std::vector<double> &m, int rel) const
{
	if (rel < 0 || rel >= selected_measurers.size()) {
		throw E_INVALID;
	}

	const std::vector<unsigned> &sel_meas = selected_measurers[rel];
	for (std::vector<unsigned>::const_iterator i = sel_meas.begin(); i != sel_meas.end(); ++i) {
		scratch[*i] = compute_membership(m[*i], stats[rel][*i].mean, stats[rel][*i].var);
		if (scratch[*i] < 0.00001) {
			return 0.0;
		}
	}

	return combine_scratch_scores(rel);
}


int
BoxLinkEstimator::train_link(const std::vector<double> &m, int rel)
{
	if (rel < 0 || rel >= selected_measurers.size()) {
		return E_INVALID;
	}

	unsigned i = 0;
	for (std::vector<scg::Measurement *>::const_iterator pm = measurers.begin(); pm != measurers.end(); ++pm) {
		const scg::Measurement *meas = *pm;
		Statistics &stat = stats[rel][i];
		double x = m[i];
		double new_mean = (stat.N * stat.mean + x) / (stat.N + 1);
		double dx = x - new_mean;
		stat.var *= stat.var;
		stat.var = (stat.N * new_mean * new_mean + stat.sum_of_squares - 2 * stat.N * new_mean * stat.mean + dx * dx) / (stat.N + 1);
		stat.sum_of_squares += x * x;
		stat.var = std::sqrt(stat.var);
		stat.mean = new_mean;
		++stat.N;
		++i;
	}

	return 0;
}


void
BoxLinkEstimator::estimate_links(const std::vector<double> &m, std::vector<double> &D) const
{
	unsigned num_directions = link_codes.size();

	D.clear();
	D.insert(D.end(), num_directions, 0.0);

	for (unsigned j = 0; j < num_directions - 1; j++) {
		D[j] = estimate_link(m, j);
	}

	VERBOSE(
		*verb_out << "\ndistribution:\n";
		for (unsigned i = 0; i < num_directions; i++) {
			*verb_out << i << " : " << D[i] << std::endl;
		}
	);
}


}



std::ostream &
operator<<(std::ostream &os, const scg::BboxLinkScore &ls)
{
	os << ls.from<< " ";
	switch (ls.type) {
	case scg::AboveRight:
		os << "UP-RIGHT-TO";
		break;
	case scg::Right:
		os << "RIGHT-TO";
		break;
	case scg::BelowRight:
		os << "DOWN-RIGHT-TO";
		break;
	case scg::BelowWide:
		os << "DOWN-TO";
		break;
	case scg::BelowNarrow:
		os << "DOWN-N-TO";
		break;
	case scg::Contains:
		os << "CONTAINS";
		break;
	}
	os << " " << ls.to << " ";
	os << "(" << ls.p << ")";

	return os;
}

