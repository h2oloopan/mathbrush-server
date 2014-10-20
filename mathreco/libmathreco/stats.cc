#include "stats.h"
#include <limits>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>

namespace scg {

static bool
isnan(double d) {
	return d != d;
}


std::ostream &
operator<<(std::ostream &os, const distribution &D) {
	D.write(os);
	return os;
}

std::istream &
operator>>(std::istream &is, distribution &D) {
	D.read(is);
	return is;
}

static double
mean(double *sample, size_t n, size_t *nused_ = 0) {
	std::sort(sample, sample + n);
	double sum = 0.0;
	size_t nused = 0;
	double Q1 = sample[n/4];
	double Q3 = sample[3*n/4];
	double IQR = 1.5 * (Q3 - Q1);
	for (size_t i = 0; i < n; ++i) {
		if (Q1 - IQR < sample[i] && sample[i] < Q3 + IQR) {
			sum += sample[i];
			++nused;
		}
	}
	if (nused_) *nused_ = nused;
	return sum/nused;
}

static double
variance(double mean, double *sample, size_t n) {
	double Q1 = sample[n/4];
	double Q3 = sample[3*n/4];
	double IQR = 1.5 * (Q3 - Q1);
	double var = 0.0;
	size_t nused = 0;
	for (double *x = sample; x != sample + n; ++x) {
		if (Q1 - IQR < *x && *x < Q3 + IQR) {
			double dx = (*x-mean);
			var += dx*dx;
			++nused;
		}
	}
	return var / (nused - 1);
}


static double
gam(double x) {
	int i,k,m;
	double ga,gr,r,z;

	static double g[] = {
		1.0,
		0.5772156649015329,
		-0.6558780715202538,
		-0.420026350340952e-1,
		0.1665386113822915,
		-0.421977345555443e-1,
		-0.9621971527877e-2,
		0.7218943246663e-2,
		-0.11651675918591e-2,
		-0.2152416741149e-3,
		0.1280502823882e-3,
		-0.201348547807e-4,
		-0.12504934821e-5,
		0.1133027232e-5,
		-0.2056338417e-6,
		0.6116095e-8,
		0.50020075e-8,
		-0.11812746e-8,
		0.1043427e-9,
		0.77823e-11,
		-0.36968e-11,
		0.51e-12,
		-0.206e-13,
		-0.54e-14,
		0.14e-14
	};

	if (x > 171.0) {
		ga = std::numeric_limits<double>::infinity();
	}
	if (x == (int)x) {
		if (x > 0.0) {
			ga = 1.0;               // use factorial
			for (i=2;i<x;i++) {
				ga *= i;
			}
		}
		else {
			ga = std::numeric_limits<double>::infinity();
		}
	}
	else {
		if (std::abs(x) > 1.0) {
			z = std::abs(x);
			m = (int)z;
			r = 1.0;
			for (k=1;k<=m;k++) {
				r *= (z-k);
			}
			z -= m;
		}
		else {
			z = x;
		}
		gr = g[24];
		for (k=23;k>=0;k--) {
			gr = gr*z+g[k];
		}
		ga = 1.0/(gr*z);
		if (std::abs(x) > 1.0) {
			ga *= r;
			if (x < 0.0) {
				ga = -M_PI/(x*ga*std::sin(M_PI*x));
			}
		}
	}
	return ga;
}

uniformdist::uniformdist(double min_, double max_) : mn(min_), mx(max_) { }

double uniformdist::maxval() const { return eval(mn); }
double
uniformdist::eval(double x) const {
	return (x < mn || x > mx) ? 0.0 : 1.0 / (mx-mn);
}

double
uniformdist::logeval(double x) const {
	return std::log(eval(x));
}

void
uniformdist::write(std::ostream &out) const {
	out << mn << ' ' << mx;
}

void
uniformdist::read(std::istream &in) {
	in >> mn >> mx;
}

diracdist::diracdist(double d_) : d(d_) { }

double diracdist::maxval() const { return eval(d); }
double
diracdist::eval(double x) const {
	return (x == d) ? 1.0 : 0.0;
}

double
diracdist::logeval(double x) const {
	return std::log(eval(x));
}

void
diracdist::write(std::ostream &out) const {
	out << d;
}

void
diracdist::read(std::istream &in) {
	in >> d;
}

gammadist::gammadist(double mu_, double lambda_)
	: mu(mu_), lambda(lambda_), c(mkc()) {
}

bool
gammadist::useful(double thres) const {
	if (ssize == 0) return false;
	return 1.96 * sstd / std::sqrt((double)ssize) <= thres;
}

double
gammadist::mkc() {
	return std::pow(lambda, mu) / gam(mu);
}

double gammadist::maxval() const { return eval((mu-1)/lambda); }
double
gammadist::eval(double x) const {
	return c * std::pow(x, mu-1) * std::exp(-lambda*x);
}

double
gammadist::logeval(double x) const {
	return std::log(c) + (mu-1) * std::log(x) - lambda*x;
}

gammadist
estimategamma(double *sample, size_t n, double thres) {
	thres = 0.5 * std::exp(-0.5*thres); // 2-dof chi-square pdf at thres
	//std::cerr << "thres is " << thres << std::endl;

	std::sort(sample, sample + n);
	size_t nused = 0;
	double Q1 = sample[n/4];
	double Q3 = sample[3*n/4];
	double IQR = 1.5 * (Q3 - Q1);

	double sum = 0.0, logsum = 0.0;
	for (double *x = sample; x != sample + n; ++x) {
		if (Q1 - IQR < *x && *x < Q3 + IQR) {
			sum += *x;
			logsum += std::log(*x);
			++nused;
		}
	}

	double mean = sum/nused;
	double var = variance(mean, sample, n);

	double r = mean * (1-mean) / var - 1;
	//return gammadist(mean * r, (1-mean) * r);

	double logmean = logsum/nused;
	double M = std::log(mean) - logmean;

	double mu = 1.0 / (2*M);
	double lambda = mu / mean;

	double DG = std::log(mu) - (1.0 + (1.0 - (0.1 - 1.0/(21*mu*mu))/(mu*mu))/(6*mu))/(2*mu);
	double m = n * (std::log(lambda) - DG) + logsum;
	double g = (n*mu)/lambda - sum;
	double TG = (1.0 + (1.0 + (1.0 - (0.2 - 1.0/(7*mu*mu))/(mu*mu))/(3*mu))/(2*mu))/mu;
	double D = n * (mu*TG - 1.0);
	double varmu = mu/D;
	double varlambda = (lambda*lambda*TG)/D;
	double cov = lambda/D;
	double mm1 = varmu * m + cov * g;
	double mm2 = cov * m + varlambda * g;
	double S = m * mm1 + g * mm2;
	while (S > thres) {
		//std::cerr << "mu, lambda = " << mu << " , " << lambda << std::endl;
		//std::cerr << "S = " << S << std::endl;
		mu += mm1;
		lambda += mm2;
		DG = std::log(mu) - (1.0 + (1.0 - (0.1 - 1.0/(21*mu*mu))/(mu*mu))/(6*mu))/(2*mu);
		m = n * (std::log(lambda) - DG) + logsum;
		g = (n*mu)/lambda - sum;
		TG = (1.0 + (1.0 + (1.0 - (0.2 - 1.0/(7*mu*mu))/(mu*mu))/(3*mu))/(2*mu))/mu;
		D = n * (mu*TG - 1.0);
		varmu = mu/D;
		varlambda = (lambda*lambda*TG)/D;
		cov = lambda/D;
		mm1 = varmu * m + cov * g;
		mm2 = cov * m + varlambda * g;
		S = m * mm1 + g * mm2;
	}
	//std::cerr << "final mu, lambda = " << mu << " , " << lambda << std::endl;
	//std::cerr << "final S = " << S << std::endl;

	gammadist G;
	G.mu = mu;
	G.lambda = lambda;
	G.sstd = std::sqrt(var);
	G.ssize = n;
	return G;
}

betadist::betadist(double alpha_, double beta_)
	: alpha(alpha_), beta(beta_), c(mkc()), sstd(0), ssize(0) {
}

bool
betadist::useful(double thres) const {
	if (ssize == 0) return false;
	return 1.96 * sstd / std::sqrt((double)ssize) <= thres;
}

double
betadist::mkc() {
	return gam(alpha+beta) / (gam(alpha) * gam(beta));
}

static double
gengamma1(int k) {
	double g = 0;
	while (k--) {
		g += -std::log(std::rand()/(double)RAND_MAX);
	}
	return g;
}

double
betadist::gen() const {
	double a = gengamma1((int)alpha);
	double b = gengamma1((int)beta);
	return a/(a+b);
}

double betadist::maxval() const { return eval((alpha-1.0)/(alpha+beta-2.0)); }
double
betadist::eval(double x) const {
	double A = std::min(1000.0, std::pow(x, alpha-1));
	double B = std::min(1000.0, std::pow(1-x, beta-1));
	return c*A*B;
	//return c * std::pow(x, alpha-1)*std::pow(1-x, beta-1);
}

double
betadist::logeval(double x) const {
	return std::log(c) + (alpha-1)*std::log(x) + (beta-1)*std::log(1-x);
}

betadist
estimatebeta(double *sample, size_t n) {
	betadist D;
	double x = mean(sample, n);
	double var = variance(x, sample, n);
	double c = x*(1-x)/var - 1;
	D.alpha = x*c;
	D.beta = (1-x)*c;
	D.sstd = std::sqrt(var);
	D.ssize = n;
	return D;
}

expdist::expdist(double lambda_) : lambda(lambda_) { }

double expdist::maxval() const { return eval(0); }

double
expdist::eval(double x) const {
	return std::min(1000.0, lambda * std::exp(-lambda * x));
}

double
expdist::logeval(double x) const {
	return std::log(lambda) - lambda*x;
}

normaldist::normaldist(double mean_, double var_) : mean(mean_), var(var_), sstd(0), ssize(0) { }
bool
normaldist::useful(double thres) const {
	if (ssize < 5 || isnan(mean) || isnan(var) || isnan(sstd)) return false;
	return 1.96 * sstd / std::sqrt((double)ssize) <= thres;
}

double normaldist::maxval() const { return eval(mean); }

double
normaldist::eval(double x) const {
	double dx = x - mean;
	dx *= dx;
	return std::exp(-dx/(2*var))/std::sqrt(var*2*M_PI);
}

double
normaldist::logeval(double x) const {
	double dx = x - mean;
	dx *= dx;
	return -dx/(2*var) - std::log(std::sqrt(var*2*M_PI));
}

normaldist
estimatenormal(double *sample, size_t n) {
	size_t nused;
	normaldist D;
	D.mean = scg::mean(sample, n, &nused);
	D.var = variance(D.mean, sample, n);
	D.ssize = nused;
	D.sstd = std::sqrt(D.var);
	return D;
}


expdist
estimateexp(double *sample, size_t n) {
	expdist D;
	double x = mean(sample, n);
	D.lambda = 1.0 / x;
	D.ssize = n;
	D.sstd = std::sqrt(variance(x, sample, n));
	return D;
	//return expdist(std::min(1.0 / mean(sample, n), 1000.0));
}

bool
expdist::useful(double thres) const {
	if (ssize == 0) return false;
	return 1.96 * sstd / std::sqrt((double)ssize) <= thres;
}


/*
misesdist::misesdist(double loc_, double conc_) : loc(loc_), conc(conc_), C(mkc()) { } 

double
misesdist::mkc() {
	return 0;
	//return 2*M_PI * bessel0(conc);
}

bool
misesdist::useful(double thres) const {
	if (ssize == 0) return false;
	return 1.96 * sstd / std::sqrt(ssize) <= thres;
}

double
misesdist::eval(double x) const {
	return std::exp(conc * std::cos(x - loc)) / C;
}

double
misesdist::logeval(double x) const {
	return conc * std::cos(x - loc) - std::log(C);
}

misesdist
estimatemises(double *sample, size_t n) {
	double xs = 0, ys = 0;
	for (double *x = sample; x != sample + n; ++x) {
		xs += std::cos(*x);
		ys += std::sin(*x);
	}
	misesdist D;
	D.loc = std::atan2(ys,xs);
	//D.conc = impossible;
	D.C = D.mkc();
	return D;
}

std::ostream &
operator<<(std::ostream &os, const misesdist &D) {
	return os << D.loc << ' ' << D.conc << ' ' << D.sstd << ' ' << D.ssize;
}

std::istream &
operator>>(std::istream &is, misesdist &D) {
	is >> D.loc >> D.conc >> D.sstd >> D.ssize;
	D.C = D.mkc();
	return is;
}*/


void
expdist::write(std::ostream &os) const {
	os << lambda << ' ' << ssize << ' ' << sstd;
}

void
expdist::read(std::istream &is) {
	is >> lambda >> ssize >> sstd;
}

void
gammadist::write(std::ostream &os) const {
	if (isnan(mu) || mu == std::numeric_limits<double>::infinity()) os << "0 ";
	else os << mu << ' ';
	if (isnan(lambda) || lambda == std::numeric_limits<double>::infinity()) os << "0 ";
	else os << lambda << ' ';
	os << ssize << ' ';
	if (isnan(sstd)) os << "1";
	else os << sstd;
}

void
gammadist::read(std::istream &is) {
	is >> mu >> lambda >> ssize >> sstd;
	c = mkc();
}

void
betadist::write(std::ostream &os) const {
	if (isnan(alpha)) os << "0 ";
	else os << alpha << ' ';
	if (isnan(beta)) os << "0 ";
	else os << beta << ' ';
	os << ssize << ' ';
	if (isnan(sstd)) os << "1";
	else os << sstd;
}

void
betadist::read(std::istream &is) {
	is >> alpha >> beta >> ssize >> sstd;
	c = mkc();
}

void
normaldist::write(std::ostream &os) const {
	if (isnan(mean)) os << "0 ";
	else os << mean << ' ';
	if (isnan(var)) os << "0 ";
	else os << var << ' ';
	os << ssize << ' ';
	if (isnan(sstd)) os << "1";
	else os << sstd;
}

void
normaldist::read(std::istream &is) {
	is >> mean >> var >> ssize >> sstd;
	var = std::max(var, 0.0001);
}

/*
diracmixdist::diracmixdist() : primary(0) {
}

diracmixdist::diracmixdist(double d_, distribution *primary_, double p_)
	: d(d_), primary(primary_), p(p_) {
}

bool
diracmixdist::useful(double thres) const {
	return primary && primary->useful(thres);
}

double
diracmixdist::eval(double x) const {
	return (x == d) ? p : (1-p)*primary->eval(x);
}

double
diracmixdist::logeval(double x) const {
	return (x == d) ? std::log(p) : (std::log(1-p) + primary->logeval(x));
}

void
diracmixdist::write(std::ostream &out) const {
	out << d << ' ' << p << ' ' << *primary;
}

void
diracmixdist::read(std::istream &in) {
	in >> d >> p >> *primary;
}
*/

powerdist::powerdist(double p_) : p(p_), sstd(0), ssize(0) { }

bool
powerdist::useful(double thres) const {
	if (ssize == 0) return false;
	return 1.96 * sstd / std::sqrt((double)ssize) <= thres;
}

double powerdist::maxval() const { return eval(0); }

double
powerdist::eval(double x) const {
	return (p-1) * std::pow(1+x, -p);
}

double
powerdist::logeval(double x) const {
	return std::log(p-1) - p * std::log(1+x);
}

void
powerdist::write(std::ostream &out) const {
	out << p << ' ' << ssize << ' ' << sstd;
}

void
powerdist::read(std::istream &in) {
	in >> p >> ssize >> sstd;
}


powerdist
estimatepower(double *x, size_t n) {
	double N = 0;
	double mean = 0;
	for (size_t i = 0; i < n; ++i) {
		N += std::log(1+x[i]);
		mean += x[i];
	}
	mean /= n;
	powerdist D((n+N) / N);
	D.ssize = n;
	D.sstd = variance(mean, x, n);
	return D;
}

}
