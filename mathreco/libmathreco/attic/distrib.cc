#include "distrib.h"
#include "verb.h"

#include <cmath>
#include <numeric>
#include <iostream>
#include <cassert>

namespace scg {

exp_dist::exp_dist()
	: mean(0), sum(0), n(0), useful_(false) { }
exp_dist::exp_dist(double mean_)
	: mean(mean_), sum(mean_), n(1), useful_(true) { }
	
double
exp_dist::score(double x) const {
	if (x < 0) {
		return 0;
	}
	return std::exp(-x/mean)/mean;
}
	
void
exp_dist::estimate(const std::vector<double> &data) {
	sum = 0.0;
	n = 0;
	for (std::vector<double>::const_iterator i = data.begin(); i != data.end(); ++i) {
		if (*i > 0) {
			sum += *i;
			++n;
		}
	}
	mean = sum / n;
	useful_ = (mean * 1.96 / std::sqrt((double)n)) < 0.05;
}

void
exp_dist::update(double x) {
	if (x > 0) {
		sum += x;
		n++;
		mean = sum / n;
		useful_ = (mean * 1.96 / std::sqrt((double)n)) < 0.05;
	}
}

bool
exp_dist::useful() const {
	return useful_;
}

void
exp_dist::read(std::istream &in) {
	in >> n >> sum >> useful_;
	mean = sum / n;
}

void
exp_dist::write(std::ostream &out) const {
	out << n << ' ' << sum << ' ' << useful_;
}


normal_dist::normal_dist()
	: mean(0), var(1), sum(0), sumsq(0), n(0), useful_(false) { }
normal_dist::normal_dist(double mean_, double var_)
	: mean(mean_), var(var_), sum(mean), sumsq(mean*mean), n(1), useful_(true) { }

bool
normal_dist::useful() const {
	return useful_;
}

void
normal_dist::estimate(const std::vector<double> &data) {
	sum = 0;
	sumsq = 0;
	n = data.size();
	for (std::vector<double>::const_iterator i = data.begin(); i != data.end(); ++i) {
		sum += *i;
		sumsq += *i * (*i);
	}
	mean = sum / n;
	var = sumsq / n - mean*mean;
	useful_ = 1.96 * std::sqrt(var/n) <= 1.0/16;
}

double
normal_dist::score(double x) const {
	double dx = x - mean;
	return std::exp(-dx*dx / (2*var));
}

void
normal_dist::update(double x) {
	sum += x;
	sumsq += x*x;
	n++;
	mean = sum / n;
	var = sumsq / n - mean*mean;
	useful_ = 1.96 * std::sqrt(var/n) <= 1.0/16;
}

void
normal_dist::read(std::istream &in) {
	in >> n >> sum >> sumsq >> useful_;
	mean = sum / n;
	var = sumsq / n - mean*mean;
}

void
normal_dist::write(std::ostream &out) const {
	out << n << ' ' << sum << ' ' << sumsq << ' ' << useful_;
}

gamma_dist::gamma_dist()
	: mu(0), lambda(0), mu_bias(0), lambda_bias(0), sum(0), logsum(0), n(0) { }
gamma_dist::gamma_dist(double mu_, double lambda_)
	: mu(mu_), lambda(lambda_), mu_bias(0), lambda_bias(0), sum(mu_/lambda_), logsum(std::log(mu_/lambda_)), n(1) { }

void
gamma_dist::read(std::istream &in) {
	in >> n >> sum >> logsum >> useful_;
	updateparms();
}

void
gamma_dist::write(std::ostream &out) const {
	out << n << ' ' << sum << ' ' << logsum << ' ' << useful_;
}

bool
gamma_dist::useful() const {
	return useful_;
}

double
gamma_dist::score(double x) const {
	if (x < 0) {
		return 0;
	}
	return std::pow(x*lambda, mu - 1) * lambda * std::exp(-x * lambda) / gamma(mu);//std::pow(x*lambda/(mu-1), mu-1) * std::exp(mu - x*lambda - 1);
}

void
gamma_dist::update(double x) {
	if (x > 0) {
		sum += x;
		logsum += std::log(x);
		n++;
		updateparms();
	}
}

void
gamma_dist::estimate(const std::vector<double> &data) {
	n = 0;
	sum = 0.0;
	logsum = 0.0;
	for (std::vector<double>::const_iterator i = data.begin(); i != data.end(); ++i) {
		if (*i > 0) {
			sum += *i;
			logsum += std::log(*i);
			++n;
		}
	}
	updateparms();
}

void
gamma_dist::updateparms() {
	static double thres = 0.5 * std::exp(-0.5 * 9.21); // 99 percentile
	double mean = sum/n;
	double logmean = logsum/n;
	
	double M = std::log(mean) - logmean;

	mu = 1.0 / (2*M);
	lambda = mu / mean;

	double DG = digamma(mu);
	double m = n * (std::log(lambda) - DG) + logsum;
	double g = (n*mu)/lambda - sum;
	double TG = trigamma(mu);
	double D = n * (mu*TG - 1.0);
	double varmu = mu/D;
	double varlambda = (lambda*lambda*TG)/D;
	double cov = lambda/D;
	double mm1 = varmu * m + cov * g;
	double mm2 = cov * m + varlambda * g;
	double S = m * mm1 + g * mm2;
	while (S > thres) {
		mu += mm1;
		lambda += mm2;
		DG = digamma(mu);
		m = n * (std::log(lambda) - DG) + logsum;
		g = (n*mu)/lambda - sum;
		TG = trigamma(mu);
		D = n * (mu*TG - 1.0);
		varmu = mu/D;
		varlambda = (lambda*lambda*TG)/D;
		cov = lambda/D;
		mm1 = varmu * m + cov * g;
		mm2 = cov * m + varlambda * g;
		S = m * mm1 + g * mm2;
	}

	TG = trigamma(mu);
	double tG = tetragamma(mu);
	double adj = mu * TG - 1;
	double den = 2*n * adj*adj; 
	mu_bias = (mu * (TG - mu*tG) - 2) / den;
	lambda_bias = lambda * (2*mu*TG*TG - 3*TG - mu*tG) / den;

	double var = mu / lambda / lambda;
	useful_ = 1.96 * std::sqrt(var/n) <= 1.0/16;
	useful_ = n > 4;
	VERBOSE(*verb_out << "gamma dist has mu " << mu << " and lambda " << lambda << "; useful " << useful_ << std::endl);
}

double
gamma_dist::gamma(double x) {
	if (x < 0.5) {
		return M_PI / std::sin(M_PI*x) / gamma(1.0 - x);
	}
	double b = 0.99999999999980993 + 676.5203681218851/(x+1) - 1259.1392167224028/(x+2)
	         + 771.32342877765313/(x+3) - 176.61502916214059/(x+4) + 12.507343278686905/(x+5)
			 - 0.13857109526572012/(x+6) + 9.9843695780195716e-6/(x+7) + 1.5056327351493116e-7/(x+8);
	double t = x + 7.5;
	return std::sqrt(2*M_PI) * std::pow(t,x+0.5) * std::exp(-t) * b;
}

double
gamma_dist::digamma(double x) {
	double x2 = x*x;
	return std::log(x) - (1 + (1 - (0.1 - 1/(21*x2))/x2)/(6*x))/(2*x);
}

double
gamma_dist::trigamma(double x) {
	double x2 = x*x;
	return (1 + (1 + (1 - (0.2 - 1/(7*x2))/x2)/(3*x))/(2*x))/x;
}

double
gamma_dist::tetragamma(double x) {
	double x2 = x*x;
	return -(1 + (1 + (1 - (1 - 1/x2)/(3*x2))/(2*x))/x)/x2;
}


vonmises_dist::vonmises_dist()
	: mu(0), K(1), sumcos(0), sumsin(0), sumcos2(0), sumsin2(0), n(0), useful_(false) {
}

vonmises_dist::vonmises_dist(double mu_, double K_)
	: mu(mu_), K(K_), sumcos(0), sumsin(0), sumcos2(0), sumsin2(0), n(0), useful_(true) {
}

double
vonmises_dist::score(double x) const {
	return std::exp(K * std::cos(x - mu)) / (2*M_PI * bessel0(K));//std::exp(K*(std::cos(x - mu) - 1));
}

void
vonmises_dist::update(double x) {
	sumsin += std::sin(x);
	sumcos += std::cos(x);
	sumsin2 += std::sin(2*x);
	sumcos2 += std::cos(2*x);
	++n;
	update_parms();
}

double
vonmises_dist::bessel0(double x) {
	double t = x/3.75;
	if (std::abs(x) <= 3.75) {
		double t2 = t*t;
		return 1.0 + t2*(3.5156229 + t2*(3.0899424 + t2*(1.2067492 + t2*(0.2659732 + t2*(0.0360768 + t2*0.0045713)))));
	}
	else if (x > 0) {
		double tinv = 1.0/t;
		double b = 0.39894228 + tinv*(0.01328592 + tinv*(0.00225319 - tinv*(0.00157565 - tinv*(0.00916281 - tinv*(0.02057706 - tinv*(0.02635537 - tinv*(0.01647633 - tinv*0.00392377)))))));
		return b * std::exp(x) / std::sqrt(x);
	}
	else {
		assert(false);
	}
}

double
vonmises_dist::bessel1(double x) {
	double t = x/3.75;
	if (std::abs(x) <= 3.75) {
		double t2 = t*t;
		return x*(0.5 + t2*(0.87890594 + t2*(0.51498869 + t2*(0.15084934 + t2*(0.02658733 + t2*(0.00301532 + t2*0.00032411))))));
	}
	else if (x > 0) {
		double tinv = 1.0/t;
		double b = 0.39894228 - tinv*(0.03988024 + tinv*(0.00362018 - tinv*(0.00163801 - tinv*(0.01031555 - tinv*(0.02282967 - tinv*(0.02895312 - tinv*(0.01787654 - tinv*(0.00420059))))))));
		return b * std::exp(x) / std::sqrt(x);
	}
	else {
		assert(false);
	}
}

double
vonmises_dist::f(double R, double x) {
	return R * bessel0(x) - bessel1(x);
}

double
vonmises_dist::df(double R, double x) {
	return bessel1(x) * (R + 1/x) - bessel0(x);
}

void
vonmises_dist::update_parms() {
	double mcos = sumcos / n;
	double msin = sumsin / n;
	double mcos2 = sumcos2 / n;
	double msin2 = sumsin2 / n;
	double R = std::sqrt(mcos*mcos + msin*msin);
	mu = std::atan2(msin, mcos);
	
	if (n < 2) {
		useful_ = false;
		VERBOSE(*verb_out << "von-mises dist has mu " << mu << " but insufficient data to try approximating K\n");
		//std::cerr << "von-mises dist has mu " << mu << " but insufficient data to try approximating K\n";
	}
	else {
		K = 2.0;
		while (f(R, K) > 0.0) {
			K *= 2;
		}
		double k0 = 0.0;
		while (std::abs(K - k0) > 0.000001) {
			k0 = K;
			K = k0 - f(R, k0)/df(R, k0); 
		}
		useful_ = std::asin(1.96 / std::sqrt(n*R*K)) < 8*M_PI/180.0;
		useful_ = n > 4;
		//std::cout << "von-mises 95% range is " << std::asin(1.96 / std::sqrt(n*R*K))*180/M_PI << " with n=" << n << std::endl;
		VERBOSE(*verb_out << "von-mises dist has mu " << mu << " and K " << K << "; useful " << useful_ << std::endl);
		//std::cerr << "von-mises dist has mu " << mu << " and K " << K << "; useful " << useful_ << std::endl;
	}
}

void
vonmises_dist::estimate(const std::vector<double> &data) {
	sumsin = sumcos = sumsin2 = sumcos2 = 0;
	n = data.size();
	for (std::vector<double>::const_iterator i = data.begin(); i != data.end(); ++i) {
		sumsin += std::sin(*i);
		sumcos += std::cos(*i);
		sumsin2 += std::sin(2* (*i));
		sumcos2 += std::cos(2* (*i));
	}
	update_parms();
}

bool
vonmises_dist::useful() const {
	return useful_;
}

void
vonmises_dist::read(std::istream &in) {
	in >> n >> sumcos >> sumsin >> sumcos2 >> sumsin2 >> useful_;
	update_parms();
}

void
vonmises_dist::write(std::ostream &out) const {
	out << n << ' ' << sumcos << ' ' << sumsin << ' ' << sumcos2 << ' ' << sumsin2 << ' ' << useful_;
}

std::istream &
operator>>(std::istream &in, distrib &D) {
	D.read(in);
	return in;
}

std::ostream &
operator<<(std::ostream &out, const distrib &D) {
	D.write(out);
	return out;
}

}