#include "gridnorm.h"
#include "error.h"
#include "parms.h"

#include <numeric>


namespace scg
{

static unsigned DefaultGridRows = RegisterParameterUnsigned("DefaultGridRows", &DefaultGridRows);
static unsigned DefaultGridCols = RegisterParameterUnsigned("DefaultGridCols", &DefaultGridCols);

static int
compute_weights(const NormalizedStrokeGroup &strokes, unsigned rows, unsigned cols, unsigned *w, unsigned &t, const std::vector<size_t> *order) {
	for (size_t i = 0; i < strokes.nstrokes; ++i) {
		const NormalizedStroke &s = strokes.strokes[order ? (*order)[i] : i];
		double *x = s.x, *y = s.y;
		for (; x != s.x + s.npoints; ++x, ++y) {
			unsigned c = (unsigned)(*x * cols);
			unsigned r = (unsigned)(*y * rows);
			w[r*cols+c] += 1;
		}
	}
	t = std::accumulate(w, w + rows*cols, 0);
	return 0;
}

static int
gridnorm_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, unsigned rows, unsigned cols, const std::vector<size_t> &inputorder, double &score) {
	static unsigned mweight[64];
	static unsigned iweight[64];
	unsigned mtot = 0;
	unsigned itot = 0;

	if (rows * cols >= 64) {
		return E_INVALID;
	}
	std::fill(mweight, mweight + 64, 0);
	std::fill(iweight, iweight + 64, 0);
	compute_weights(model, rows, cols, mweight, mtot, 0);
	compute_weights(input, rows, cols, iweight, itot, &inputorder);

	score = 0;
	unsigned *m = mweight, *i = iweight;
	for (; m != mweight + rows*cols; ++m, ++i) {
		double nm = (double)*m / mtot;
		double ni = (double)*i / itot;
		nm -= ni;
		score += nm * nm;
	}

	score = std::sqrt(score);// / (DefaultGridRows * DefaultGridCols));
	//score = std::max(0.0, 2.0 * (score - 0.5));
	return 0;
}


int
gridnorm_match(const NormalizedStrokeGroup &model, const NormalizedStrokeGroup &input, const std::vector<size_t> &inputorder, double &score) {
	unsigned r = DefaultGridRows;
	unsigned c = DefaultGridCols;
	int e = gridnorm_match(model, input, r, c, inputorder, score);
	if (FAILURE(e)) {
		return e;
	}
	if (score < 0.01) {
		score = 0.01;
	}
	return 0;
}


}
