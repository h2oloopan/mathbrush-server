#ifndef GAMMA_H_
#define GAMMA_H_

#include <cstdlib>
#include <istream>
#include <ostream>

namespace scg {

struct distribution {
	virtual ~distribution() { }
	virtual bool useful(double thres) const = 0;
	virtual double maxval() const = 0;
	virtual double eval(double x) const = 0;
	virtual double logeval(double x) const = 0;

	virtual void write(std::ostream &out) const = 0;
	virtual void read(std::istream &in) = 0;
};

struct uniformdist : public distribution {
	double mn, mx;
	uniformdist() { }
	uniformdist(double min_, double max_);
	bool useful(double) const { return true; }
	double maxval() const;
	double eval(double x) const;
	double logeval(double x) const;
	void write(std::ostream &out) const;
	void read(std::istream &in);
};

/*
struct diracmixdist : public distribution {
	double d;
	distribution *primary;
	double p;

	diracmixdist();
	diracmixdist(double d_, distribution *primary_, double p_);
	bool useful(double) const;
	double maxval() const;
	double eval(double x) const;
	double logeval(double x) const;
	void write(std::ostream &out) const;
	void read(std::istream &in);
};*/

struct diracdist : public distribution {
	double d;

	diracdist() { }
	diracdist(double d_);
	bool useful(double) const { return true; }
	double maxval() const;
	double eval(double x) const;
	double logeval(double x) const;

	void write(std::ostream &out) const;
	void read(std::istream &in);
};

struct gammadist : public distribution {
	double mu;
	double lambda;
	double c;

	double sstd;
	unsigned ssize;

	gammadist() { }
	gammadist(double mu_, double lambda_);
	bool useful(double thres) const;
	double maxval() const;
	double eval(double x) const;
	double logeval(double x) const;
	double mkc();

	void write(std::ostream &out) const;
	void read(std::istream &in);
};
gammadist estimategamma(double *sample, size_t n, double thres = 9.61);

struct betadist : public distribution {
	double alpha;
	double beta;
	double c;

	double sstd;
	unsigned ssize;

	betadist() { }
	betadist(double alpha_, double beta_);
	bool useful(double thres) const;
	double maxval() const;
	double eval(double x) const;
	double logeval(double x) const;
	double mkc();
	double gen() const;

	void write(std::ostream &out) const;
	void read(std::istream &in);
};
betadist estimatebeta(double *samples, size_t n);

struct powerdist : public distribution {
	double p;
	double sstd;
	unsigned ssize;

	powerdist() { }
	powerdist(double p_);
	bool useful(double thres) const;
	double maxval() const;
	double eval(double x) const;
	double logeval(double x) const;
	void write(std::ostream &out) const;
	void read(std::istream &in);
};
powerdist estimatepower(double *samples, size_t n);

struct expdist : public distribution {
	double lambda;
	expdist() { }
	expdist(double lambda_);
	double maxval() const;
	double eval(double x) const;
	double logeval(double x) const;

	double sstd;
	unsigned ssize;

	void write(std::ostream &out) const;
	void read(std::istream &in);
	bool useful(double thres) const;
};
expdist estimateexp(double *samples, size_t n);

struct normaldist : public distribution {
	double mean;
	double var;

	double sstd;
	unsigned ssize;

	normaldist() : ssize(0) { }
	normaldist(double mean, double var);
	bool useful(double thres) const;
	double maxval() const;
	double eval(double x) const;
	double logeval(double x) const;

	void write(std::ostream &out) const;
	void read(std::istream &in);
};
normaldist estimatenormal(double *samples, size_t n);

/*
struct misesdist : public distribution {
	double loc;
	double conc;
	double C;

	double sstd;
	unsigned ssize;
	
	misesdist() { }
	misesdist(double loc_, double conc_);
	bool useful(double thres) const;
	double eval(double x) const;
	double logeval(double x) const;
	double mkc();
};
*/

/*
std::ostream &operator<<(std::ostream &os, const expdist &D);
std::istream &operator>>(std::istream &is, expdist &D);
std::ostream &operator<<(std::ostream &os, const betadist &D);
std::istream &operator>>(std::istream &is, betadist &D);
std::ostream &operator<<(std::ostream &os, const gammadist &D);
std::istream &operator>>(std::istream &is, gammadist &D);
std::ostream &operator<<(std::ostream &os, const normaldist &D);
std::istream &operator>>(std::istream &is, normaldist &D);
std::ostream &operator<<(std::ostream &os, const misesdist &D);
std::istream &operator>>(std::istream &is, misesdist &D);
*/

std::ostream &operator<<(std::ostream &os, const distribution &D);
std::istream &operator>>(std::istream &is, distribution &D);

}

#endif
