#ifndef MEAS_H_
#define MEAS_H_


#include "dist.h"
#include "expr.h"
#include "vector-io.h"

#include <map>
#include <string>
#include <vector>


namespace scg
{

extern long MinStrokeSize;


typedef Rect<long> Bbox;

void find_tl_point(const Bbox &box, long *pt);
void find_tr_point(const Bbox &box, long *pt);
void find_bl_point(const Bbox &box, long *pt);
void find_br_point(const Bbox &box, long *pt);

template <typename T>
double minnorm(const Rect<T> &box1, const Rect<T> &box2)
{
	return (double)std::min(std::max(width(box1), height(box1)), std::max(width(box2), height(box2)));
}

template <typename T>
double norm(const Rect<T> &box1, const Rect<T> &box2)
{
	return (double)std::max(std::max(width(box1), height(box1)), std::max(width(box2), height(box2)));
	//return (double)std::max(height(box1), height(box2));
}

template <typename T>
double xnorm(const Rect<T> &box1, const Rect<T> &box2)
{
	return norm(box1, box2);//std::max(T(1), width(box2));//(double)std::max(width(box1), width(box2)); 
}

template <typename T>
double ynorm(const Rect<T> &box1, const Rect<T> &box2)
{
	return norm(box1, box2);//std::max(T(1), height(box2));//(double)std::max(height(box1), height(box2));
}

template <typename T>
double ynorm(const Rect<T> &box)
{
	return height(box) + MinStrokeSize;
	return std::max<long>(height(box), MinStrokeSize);
	//return height(box);
}
template <typename T>
double xnorm(const Rect<T> &box)
{
	return width(box) + MinStrokeSize;
	return std::max<long>(width(box), MinStrokeSize);
	//return width(box);
}

template <typename T>
double norm(const Rect<T> &box)
{
	return std::max<long>(xnorm(box), ynorm(box));
	/*
	double w = width(box);
	double h = height(box);
	return std::sqrt(w*w + h*h);
	*/
}


struct Measurement
{
	virtual double measure(const Bbox &box1, const Bbox &box2) const = 0;
	virtual std::string name() const = 0;
};


static double
diag_length(const Bbox &b)
{
	double w = width(b);
	double h = height(b);
	return std::sqrt(w*w + h*h);
}

struct horz_align : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)(std::max(box2.left - box1.right, (long)0) + std::max(box1.left - box2.right, (long)0)) / norm(box1, box2); }
	std::string name() const { return "horz_align"; }
};

struct b2y2 : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.bottom - box1.top) / ynorm(box1, box2); }
	std::string name() const { return "b2y2"; }
};

struct b2x1a : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.left - box1.right) / xnorm(box1, box2); }
	std::string name() const { return "b2x1a"; }
};

struct angle_br_tl : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.bottom, box2.left - box1.right); }
	std::string name() const { return "angle_br_tl"; }
};

struct angle_br_bl : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.bottom - box1.bottom, box2.left - box1.right); }
	std::string name() const { return "angle_br_bl"; }
};

struct angle_tl_tr : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.top, box2.right - box1.left); }
	std::string name() const { return "angle_tl_tr"; }
};

struct angle_tr_tl : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.top, box2.left - box1.right); }
	std::string name() const { return "angle_tr_tl"; }
};


struct cx : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)((box2.right + box2.left) / 2 - (box1.right + box1.left) / 2) / diag_length(box1); }
	std::string name() const { return "cx"; }
};
struct cy : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)((box2.bottom + box2.top) / 2 - (box1.bottom + box1.top) / 2) / diag_length(box1); }
	std::string name() const { return "cy"; }
};
struct ct : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2(((double)box2.bottom + box2.top) / 2 - (box1.bottom + box1.top) / 2, (box2.right + box2.left) / 2 - (box1.right + box1.left) / 2); }
	std::string name() const { return "ct"; }
};


struct dx : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)(box2.left - box1.right) / xnorm(box1); }
	std::string name() const { return "dx"; }
};

struct dy : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.top - box1.bottom) / ynorm(box1); }
	std::string name() const { return "dy"; }
};

struct dright : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)(box2.right - box1.right) / xnorm(box1); }
	std::string name() const { return "dright"; }
};

struct dtop : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)(box2.top - box1.top) / std::max(height(box1), height(box2)); }
	std::string name() const { return "dtop"; }
};

struct dmid : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)((box2.top + box2.bottom)/2.0 - (box1.top + box1.bottom)/2.0) / std::max(height(box1), height(box2)); }
	std::string name() const { return "dmid"; }
};

struct dbot : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)(box2.bottom - box1.bottom) / std::max(height(box1), height(box2)); }
	std::string name() const { return "dbot"; }
};

struct dbaseline : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (double)(box2.bottom - box1.bottom) / ynorm(box1); }
	std::string name() const { return "dbaseline"; }
};

struct midangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((box2.top + (double)height(box2)/2) - (box1.top + (double)height(box1)/2), (double)box2.left - box1.right); }
	std::string name() const { return "midangle"; }
};

struct bltrmag : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
	{
		double dx = box2.right - box1.left;
		double dy = box2.top - box1.bottom;
		return std::sqrt(dx*dx + dy*dy) / std::max(xnorm(box1), xnorm(box2));
	}
	std::string name() const { return "bltrmag"; }
};

struct brtlangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.bottom, box2.left - box1.right); }
	std::string name() const { return "brtlangle"; }
};

struct trblangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.bottom - box1.top, box2.left - box1.right); }
	std::string name() const { return "trblangle"; }
};

struct blblangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.bottom - box1.bottom, box2.left - box1.left); }
	std::string name() const { return "blblangle"; }
};

struct brbrangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.bottom - box1.bottom, box2.right - box1.right); }
	std::string name() const { return "brbrangle"; }
};

struct brcangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2(((double)box2.bottom + box2.top)/2 - box1.bottom, ((double)box2.right + box2.left)/2 - box1.right); }
	std::string name() const { return "brcangle"; }
};

struct trcangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2(((double)box2.bottom + box2.top)/2 - box1.top, ((double)box2.right + box2.left)/2 - box1.right); }
	std::string name() const { return "trcangle"; }
};

struct trtrangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.top, box2.right - box1.right); }
	std::string name() const { return "trtrangle"; }
};

struct bltrangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.bottom, box2.right - box1.left); }
	std::string name() const { return "bltrangle"; }
};

struct tltlangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.top, box2.left - box1.left); }
	std::string name() const { return "tltlangle"; }
};

struct tltrangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.top, box2.right - box1.left); }
	std::string name() const { return "tltrangle"; }
};

struct trtlangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.top - box1.top, box2.left - box1.right); }
	std::string name() const { return "trtlangle"; }
};

struct brblangle : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return std::atan2((double)box2.bottom - box1.bottom, box2.left - box1.right); }
	std::string name() const { return "brblangle"; }
};

struct overlap : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return overlap_proportion(box1, box2); }
	std::string name() const { return "overlap"; }
};

struct drl : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.left - box1.right) / xnorm(box1); }
	std::string name() const { return "drl"; }
};

struct drr : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.right - box1.right) / xnorm(box1); }
	std::string name() const { return "drr"; }
};

struct dll : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.left - box1.left) / xnorm(box1); }
	std::string name() const { return "dll"; }
};

struct dlr : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.right - box1.left) / xnorm(box1); }
	std::string name() const { return "dlr"; }
};

struct dbt : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.top - box1.bottom) / ynorm(box1); }
	std::string name() const { return "dbt"; }
};

struct dbb : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.bottom - box1.bottom) / ynorm(box1); }
	std::string name() const { return "dbb"; }
};

struct dtt : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.top - box1.top) / ynorm(box1); }
	std::string name() const { return "dtt"; }
};

struct dtb : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.bottom - box1.top) / ynorm(box1); }
	std::string name() const { return "dtb"; }
};


struct topray : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return ((double)box2.bottom - box1.top) / ynorm(box2); }
	std::string name() const { return "topray"; }
};

struct midray : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return (box2.bottom - ((double)box1.bottom + box1.top)/2) / ynorm(box2); }
	std::string name() const { return "midray"; }
};

struct botray : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return ((double)box2.bottom - box1.bottom) / ynorm(box2); }
	std::string name() const { return "botray"; }
};

struct leftray : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return ((double)box1.left - box2.left) / xnorm(box2); }
	std::string name() const { return "leftray"; }
};

struct rightray : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return ((double)box1.right - box2.left) / xnorm(box2); }
	std::string name() const { return "rightray"; }
};

struct relheight : public Measurement {
	double measure(const Bbox &box1, const Bbox &box2) const
		{ return ((double)ynorm(box2)/ynorm(box1)); }
	std::string name() const { return "relheight"; }
};


typedef std::vector<double> LinkMeasurements;

LinkMeasurements measure(const std::vector<Measurement *> &measurers, const Bbox &box1, const Bbox &box2);
LinkMeasurements measure(const std::map<std::string, Measurement *> &measurers, const Bbox &box1, const Bbox &box2);

}


#endif
