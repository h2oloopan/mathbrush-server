#include "feat.h"
#include "dist.h"

#include <cassert>
#include <cmath>


namespace scg
{


void
stroke_feature_set::prepare_for_matching() {
	// normalize the weights
	double s = 0.0;
	for (const_iterator i = begin(); i != end(); ++i) {
		s += i->w;
	}

	for (iterator i = begin(); i != end(); ++i) {
		i->w /= s;
	}
}



double
feature_match(const stroke_feature_set &features, const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const std::vector<size_t> &inputorder) {
	assert(model.nstrokes == input.nstrokes);

	double s = 0;
	for (size_t i = 0; i < model.nstrokes; ++i) {
		s += feature_match(features, model.strokes[i], input.strokes[inputorder[i]]);
	}

	return s / model.nstrokes;
}


double
feature_match(const stroke_feature_set &features, const NormalizedStroke &model, const NormalizedStroke &input) {
	double s = 0;
	for (stroke_feature_set::const_iterator feat = features.begin(); feat != features.end(); ++feat) {
		s += feat->w * feat->f->diff(feat->f->measure(model), feat->f->measure(input));
	}
	return s;
}


double
stroke_feature::diff(double lhs, double rhs) const {
	return std::abs(lhs - rhs);
}


double
feat_stroke_length::measure(const NormalizedStroke &s) const {
	if (s.npoints < 2) {
		return 0;
	}

	double L = 0;
	double *x, *y;
	for (x = s.x, y = s.y; x != s.x + s.npoints - 1; ++x, ++y) {
		L += std::sqrt((double)dist_sq(*x, *y, *(x + 1), *(y + 1)));
	}

	return L;
}


double
feat_aspect_ratio::measure(const NormalizedStroke &s) const {
	Rect<double> bounds = s.bounds();
	return (double)bounds.width() / bounds.height();
}


double
feat_endtoend_dist::measure(const NormalizedStroke &s) const {
	double dx = s.x[s.npoints - 1] - s.x[0];
	double dy = s.y[s.npoints - 1] - s.y[0];
	return std::sqrt(dx * dx + dy * dy);
}

}
