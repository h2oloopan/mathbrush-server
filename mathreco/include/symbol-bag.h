#ifndef SYMBOL_BAG_H_
#define SYMBOL_BAG_H_

#include "symbols.h"
#include <istream>
#include <ostream>
#include <map>


namespace scg {


class symbolbag {
private:
	struct symboldata {
		unsigned freq;
		std::map<const symbol *, unsigned> cofreq;
	};

	unsigned nsamples;
	double freqsum;
	std::map<const symbol *, symboldata> freqtab;
	std::map<const symbol *, double> symweights;
	std::map<const symbol *, double> nextp;
	double uniformity;
public:
	symbolbag();

	int addsymbol(const symbol *S, double p = 1.0);

	double query(const symbol *S);

	int read(std::istream &is, double uniformity = 0.0);
	int write(std::ostream &os) const;

	int import(const symbolbag &other);

private:
	friend std::ostream &operator<<(std::ostream &os, const symbolbag &bag);
};

std::ostream &operator<<(std::ostream &os, const symbolbag &bag);

}

#endif
