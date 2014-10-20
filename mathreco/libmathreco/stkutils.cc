#include "stkutils.h"
#include "segment.h"
#include "parms.h"

namespace scg {

static long MaxDotSize = 0;
static void calc_parms() {
	MaxDotSize = static_cast<long>(TABLETPC_DPI * GetParameterDouble("MaxDotSize"));
}
static int _e = RegisterParameterCallback(&calc_parms);

bool
bbox_is_dot(const Rect<long> &bbox) {
	if (bbox.width() == 0 && bbox.height() == 0) {
		return true;
	}
	VERBOSE(*verb_out << "MaxDotSize is " << MaxDotSize << " vs " << bbox.width() << " and " << bbox.height() << std::endl);
	if (bbox.width() > MaxDotSize || bbox.height() > MaxDotSize) {
		return false;
	}
	double hgt = std::max<double>(bbox.height(), TABLETPC_DPI/32.0);
	double wid = std::max<double>(bbox.height(), TABLETPC_DPI/32.0);
	if (std::max(wid,hgt) / std::min(wid,hgt) > 3) {
		return false;
	}
	return true;
}

bool
stroke_is_dot(const RawStroke &stroke) {
	return bbox_is_dot(stroke.bounds());
}


RawStroke
triple(const RawStroke &stk) {
	size_t n = stk.npoints;
	n = 3*(n-1)+1;
	long *x = new long[n];
	long *y = new long[n];
	unsigned long *t = 0;
	if (stk.time) {
		t = new unsigned long[n];
	}

	long *px = x, *py = y;
	unsigned long *pt = t;
	long *sx = stk.x, *sy = stk.y;
	unsigned long *st = stk.time;
	while (sx < stk.x + stk.npoints - 1) {
		long dx = sx[1] - sx[0];
		long dy = sy[1] - sy[0];
		unsigned long dt;
		if (pt) {
			dt = st[1] - st[0];
		}
		for (long i = 0; i <= 2; ++i) {
			*px = *sx + (long)(dx * (i / 3.0));
			*py = *sy + (long)(dy * (i / 3.0));
			if (pt) {
				*pt = *st + (unsigned long)(dt * (i / 3.0));
			}
			++px; ++py;
			if (pt) {
				++pt;
			}
		}
		++sx; ++sy;
		if (pt) {
			++st;
		}
	}
	*px = *sx;
	*py = *sy;
	if (pt) {
		*pt = *st;
	}
	++px, ++py;
	assert(px - x == n);
	return RawStroke(x, y, t, n);
}

RawStrokeGroup
triple(const RawStrokeGroup &stks) {
	RawStroke *outstks = new RawStroke[stks.nstrokes];
	for (size_t i = 0; i < stks.nstrokes; ++i) {
		outstks[i] = triple(stks.strokes[i]);
	}
	return RawStrokeGroup(outstks, stks.nstrokes);
}

}
