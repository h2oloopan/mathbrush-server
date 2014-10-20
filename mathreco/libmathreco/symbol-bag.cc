#include "symbol-bag.h"
#include "error.h"
#include "verb.h"
#include <cassert>


namespace scg {


//unsigned symbolbag::nsamples = 0;
//std::map<const symbol *, symbolbag::symboldata> symbolbag::freqtab;


int
symbolbag::read(std::istream &is, double uniformity_) {
	uniformity = uniformity_;
	is >> nsamples;
	if (!is) return E_IO;

	unsigned totfreq = 0;
	for (;;) {
		std::string name;
		is >> name;
		if (!is) break;
		const symbol *baseS = symdb_findsymbol_name(name);
		/*if (!baseS) {
			return E_INVALID;
		}*/
		symboldata *symd = 0;
		unsigned symbol_freq;
		unsigned num_co_symbols;
		is >> symbol_freq >> num_co_symbols;
		if (!is) return E_IO;
		if (baseS) {
			symd = &freqtab[baseS];
			symd->freq = symbol_freq;
			nextp[baseS] = symbol_freq;
			freqsum += symbol_freq;
			totfreq += symbol_freq;
		}

		while (num_co_symbols--) {
			is >> name;
			if (!is) return E_IO;
			const symbol *coS = symdb_findsymbol_name(name);
			/*if (!coS) {
				return E_INVALID;
			}*/
			unsigned cofreq;
			is >> cofreq;
			if (!is) return E_IO;
			if (baseS) {
				if (coS) {
					symd->cofreq[coS] = cofreq;
				}
			}
			else {
			}
		}
	}

	/*double uniform = (double)totfreq / freqtab.size();
	VERBOSE(*verb_out << "uniform factor = " << totfreq << " / " << freqtab.size() << " = " << uniform << std::endl);
	//double sum = 0.0;
	//freqsum = 0.0;
	for (std::map<const symbol *, double>::iterator i = nextp.begin(); i != nextp.end(); ++i) {
		double updated = uniformity * uniform + (1.0 - uniformity) * i->second;
		VERBOSE(*verb_out << "modding freq for " << i->first->name << " from " << i->second << " to " << uniformity << " * " << uniform << " + " << 1.0 - uniformity << " * " << i->second << " = " << updated << std::endl);
		//i->second = uniformity * uniform + (1.0 - uniformity) * i->second;
		//freqsum += i->second;
	}*/

	//VERBOSE(*verb_out << "initialized symbol bag to\n" << *this << std::endl);
	return 0;
}


int
symbolbag::write(std::ostream &os) const {
	os << nsamples << std::endl;
	if (!os) return E_IO;
	for (std::map<const symbol *, symboldata>::const_iterator i = freqtab.begin(); i != freqtab.end(); ++i) {
		const symboldata &symd = i->second;
		os << i->first->name << ' ' << symd.freq << ' ' << symd.cofreq.size() << std::endl;
		if (!os) return E_IO;
		for (std::map<const symbol *, unsigned>::const_iterator j = symd.cofreq.begin(); j != symd.cofreq.end(); ++j) {
			os << j->first->name << ' ' << j->second << std::endl;
			if (!os) return E_IO;
		}
	}
	return 0;
}


int
symbolbag::import(const symbolbag &other) {
	nsamples += other.nsamples;
	for (std::map<const symbol *, symboldata>::const_iterator i = other.freqtab.begin(); i != other.freqtab.end(); ++i) {
		symboldata &symd = freqtab[i->first];
		const symboldata &other_symd = i->second;
		symd.freq += other_symd.freq;
		for (std::map<const symbol *, unsigned>::const_iterator j = other_symd.cofreq.begin(); j != other_symd.cofreq.end(); ++j) {
			symd.cofreq[j->first] += j->second;
		}
	}
	return 0;
}


symbolbag::symbolbag() : freqsum(0), uniformity(0.5) {
}


int
symbolbag::addsymbol(const symbol *S, double p) {
	double &w = symweights[S];
	//p = (1.0 / nextp.size());//uniformity * (1.0 / nextp.size()) + (1.0 - uniformity) * p;
	if (p > w) {
		double adj = p - w;
		w = p;
		const symboldata &symd = freqtab[S];
		for (std::map<const symbol *, unsigned>::const_iterator i = symd.cofreq.begin(); i != symd.cofreq.end(); ++i) {
			double freq = adj * i->second;
			nextp[i->first] += freq;
			freqsum += freq;
		}
	}

	return 0;
}

static const symbol *
remap_symbol(const symbol *S) {
	if (S->name == "iota") {
		S = symdb_findsymbol_name("i");
	}
	else if (S->name == "jota") {
		S = symdb_findsymbol_name("j");
	}
	return S;
}

double
symbolbag::query(const symbol *S) {
	S = remap_symbol(S);
	//VERBOSE(*verb_out << "query " << S->name << ": " << nextp[S] << " / " << freqsum << " = " << nextp[S] / freqsum << std::endl);
	double uniform_part = 1.0 / nextp.size();
	double symbol_part = nextp[S] / freqsum;
	//VERBOSE(*verb_out << "query " << S->name << ": " << uniformity << " * " << (1.0 / nextp.size()) << " + " << 1.0 - uniformity << " * " << nextp[S] / freqsum << " = " << uniformity * uniform_part + (1.0 - uniformity) * symbol_part << std::endl);
	return uniformity * uniform_part + (1.0 - uniformity) * symbol_part;
}


std::ostream &
operator<<(std::ostream &os, const symbolbag &bag) {
	for (std::map<const symbol *, double>::const_iterator i = bag.nextp.begin(); i != bag.nextp.end(); ++i) {
		os << i->first->name << " : " << i->second / bag.freqsum << std::endl;
	}
	return os;
}


}
