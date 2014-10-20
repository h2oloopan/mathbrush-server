#ifndef FEATURES_H_
#define FEATURES_H_


#include "stroke.h"
#include "group.h"
#include "rect.h"

#include <vector>

namespace scg
{


struct stroke_feature {
	virtual ~stroke_feature() { }
	virtual double measure(const NormalizedStroke &s) const = 0;

	virtual double diff(double lhs, double rhs) const;
};



struct stroke_feature_set {
	~stroke_feature_set()
	{
		while (!features.empty()) {
			delete features.back().f;
			features.pop_back();
		}
	}

	// protocol for use:
	//  1. call add_feature once for each feature you wish to use
	//  2. call prepare_for_matching()
	//  3. use the object with feature_match
	inline void add_feature(stroke_feature *f, double weight)
		{ features.push_back(feature_spec(f, weight)); }

	void prepare_for_matching();
	

private:
	struct feature_spec {
		stroke_feature *f;
		double w;

		feature_spec(stroke_feature *f_, double w_) : f(f_), w(w_) { }
	};
	std::vector<feature_spec> features;


public:
	typedef std::vector<feature_spec>::iterator iterator;
	typedef std::vector<feature_spec>::const_iterator const_iterator;

	inline iterator begin() { return features.begin(); }
	inline iterator end() { return features.end(); }
	inline const_iterator begin() const { return features.begin(); }
	inline const_iterator end() const { return features.end(); }
};


double feature_match(const stroke_feature_set &features, const NormalizedStroke &model, const NormalizedStroke &input);
double feature_match(const stroke_feature_set &features, const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const std::vector<size_t> &inputorder);



struct feat_stroke_top : public stroke_feature {
	double measure(const NormalizedStroke &s) const
		{ return s.bounds().top; }
};


struct feat_stroke_left : public stroke_feature {
	double measure(const NormalizedStroke &s) const
		{ return s.bounds().left; }
};

struct feat_stroke_width : public stroke_feature {
	double measure(const NormalizedStroke &s) const
		{ return s.bounds().width(); }
};

struct feat_stroke_height : public stroke_feature {
	double measure(const NormalizedStroke &s) const
		{ return s.bounds().height(); }
};

struct feat_stroke_first_x : public stroke_feature {
	double measure(const NormalizedStroke &s) const
		{ return s.x[0]; }
};

struct feat_stroke_first_y : public stroke_feature {
	double measure(const NormalizedStroke &s) const
		{ return s.y[0]; }
};

struct feat_stroke_last_x : public stroke_feature {
	double measure(const NormalizedStroke &s) const
		{ return s.x[s.npoints - 1]; }
};

struct feat_stroke_last_y : public stroke_feature {
	double measure(const NormalizedStroke &s) const
		{ return s.y[s.npoints - 1]; }
};

struct feat_stroke_length : public stroke_feature {
	double measure(const NormalizedStroke &s) const;
};


struct feat_aspect_ratio : public stroke_feature {
	double measure(const NormalizedStroke &s) const;
};

struct feat_endtoend_dist : public stroke_feature {
	double measure(const NormalizedStroke &s) const;
};



}


#endif
