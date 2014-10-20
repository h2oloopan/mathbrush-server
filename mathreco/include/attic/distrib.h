#ifndef DISTRIB_H_
#define DISTRIB_H_

#include <ostream>
#include <istream>
#include <vector>

namespace scg
{

class distrib {
public:
	virtual void estimate(const std::vector<double> &data) = 0;
	virtual double score(double x) const = 0;
	virtual void update(double x) = 0;
	virtual bool useful() const = 0;
protected:
	virtual void read(std::istream &in) = 0;
	virtual void write(std::ostream &out) const = 0;
	friend std::istream &operator>>(std::istream &, distrib &);
	friend std::ostream &operator<<(std::ostream &, const distrib &);
};

std::istream &operator>>(std::istream &in, distrib &D);
std::ostream &operator<<(std::ostream &out, const distrib &D);

class exp_dist : public distrib {
public:
	exp_dist();
	explicit exp_dist(double mean_);
	
	double score(double x) const;
	void estimate(const std::vector<double> &data);
	void update(double x);
	bool useful() const;

private:
	void read(std::istream &in);
	void write(std::ostream &out) const;

private:
	double mean;
	double sum;
	size_t n;
	bool useful_;
};

class normal_dist : public distrib {
public:
	normal_dist();
	normal_dist(double mean_, double var_);
	void estimate(const std::vector<double> &data);
	double score(double x) const;
	void update(double x);
	bool useful() const;

private:
	void read(std::istream &in);
	void write(std::ostream &out) const;

private:
	double mean;
	double var;
	double sum;
	double sumsq;
	size_t n;
	bool useful_;
};

class gamma_dist : public distrib {
public:
	gamma_dist();
	gamma_dist(double mu_, double lambda_);

	double score(double x) const;
	void update(double x);
	void estimate(const std::vector<double> &data);
	bool useful() const;

private:
	void updateparms();

private:
	static double gamma(double x);
	static double digamma(double x);
	static double trigamma(double x);
	static double tetragamma(double x);

private:
	void read(std::istream &in);
	void write(std::ostream &out) const;

private:
	double mu;
	double lambda;
	double mu_bias;
	double lambda_bias;
	double sum;
	double logsum;
	size_t n;
	bool useful_;
};

class vonmises_dist : public distrib {
public:
	vonmises_dist();
	vonmises_dist(double mu, double K);

	double score(double x) const;
	void update(double x);
	void estimate(const std::vector<double> &data);
	bool useful() const;

private:
	void update_parms();
	static double bessel0(double x);
	static double bessel1(double x);
	static double f(double R, double x);
	static double df(double R, double x);

private:
	void read(std::istream &in);
	void write(std::ostream &out) const;

private:
	double mu;
	double K;
	double sumcos;
	double sumsin;
	double sumcos2;
	double sumsin2;
	size_t n;
	bool useful_;
};

}


#endif

