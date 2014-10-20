#ifndef ORDERED_SEGMENTS_H_
#define ORDERED_SEGMENTS_H_


#include "bitvec.h"
#include "order.h"
#include "rect.h"
#include "grammar-values.h"
#include "expr-node.h"

#include <vector>
#include <iostream>

namespace scg
{


struct segment;
struct nonterminal;


class ordered_segments {
private:
	typedef std::vector<const segment *> segvec;

public:
	typedef segvec::const_iterator ordered_iterator;

public:
	ordered_segments(const std::vector<segment *> &input_, const bitvec &bits_);
	//ordered_segments(const ordered_segments &boxes_);
	inline bool valid() const
		{ return bits.size() > 0; }
	
	inline void invalidate()
		{ bits.clear(); }
	
	inline size_t nsegments() const
		{ return input.size(); }
	inline size_t size() const
		{ return bits.count_set_bits(); }

	bitvec slice_bits(size_t dim, size_t start, size_t end) const;
	ordered_segments *slice(size_t dim, size_t start, size_t end) const;
	ordered_segments &slice_inplace(size_t dim, size_t start, size_t end);

	inline ordered_iterator begin(size_t d) const
	{
		if (segs[d].empty()) {
			build_order(d);
		}
		return segs[d].begin();
	}

	inline ordered_iterator end(size_t d) const
	{
		if (segs[d].empty()) {
			build_order(d);
		}
		return segs[d].end();
	}

	inline const segment *min(size_t d) const {
		return *begin(d);
	}
	inline const segment *max(size_t d) const {
		ordered_iterator i = end(d);
		--i;
		return *i;
	}

	const Rect<long> &bounds() const;

	inline void translate(int x, int y)
	{
		if (bbox.top != -1) {
			bbox.left += x;
			bbox.top += y;
			bbox.right += x;
			bbox.bottom += y;
		}
	}

	double x_center() const {
		long sx = 0;
		size_t npts = 0;
		for (ordered_iterator i = begin(0); i != end(0); ++i) {
			const RawStroke &stk = (*i)->stk;
			for (size_t j = 0; j < stk.npoints; ++j) {
				sx += stk.x[j];
			}
			npts += stk.npoints;
		}
		return (double)sx/npts;
	}
	double y_center() const {
		long sy = 0;
		size_t npts = 0;
		for (ordered_iterator i = begin(0); i != end(0); ++i) {
			const RawStroke &stk = (*i)->stk;
			for (size_t j = 0; j < stk.npoints; ++j) {
				sy += stk.y[j];
			}
			npts += stk.npoints;
		}
		return (double)sy/npts;
	}


private:
	ordered_segments(const std::vector<segment *> &input_);

private:
	void build_order(size_t d) const;

public:
	const std::vector<segment *> &input;
	double overlap;

	mutable Rect<long> bbox;
	mutable segvec segs[NUM_ORDERS];

	bitvec bits;
};


}


#endif
