#include "ordered-segments.h"
#include "reco-types.h"
#include "verb.h"

#include <algorithm>
#include <cassert>
#include <cstddef>


namespace scg
{


ordered_segments::ordered_segments(const std::vector<segment *> &input_, const bitvec &bits_/*, unsigned edit_version_*/)
	: input(input_), bits(bits_), overlap(0)
{
	bbox.top = -1;
}


ordered_segments::ordered_segments(const std::vector<segment *> &input_/*, unsigned edit_version_*/)
	: input(input_), bits(input_.size(), false), overlap(0)
{
	bbox.top = -1;
}

bitvec
ordered_segments::slice_bits(size_t dim, size_t start, size_t end) const
{
	if (end == start) {
		THROW_ERROR(E_INVALID, "taking slice of segments of size 0");
	}

	if (segs[dim].empty()) {
		build_order(dim);
	}

	assert((ptrdiff_t)(end - start) <= this->end(dim) - begin(dim));
	bitvec dstbits(bits.size(), false);
	for (ordered_iterator i = begin(dim) + start; i != begin(dim) + end; ++i) {
		dstbits.set((*i)->pos, true);
	}
		
	return dstbits;
}

ordered_segments *
ordered_segments::slice(size_t dim, size_t start, size_t end) const
{
	if (end == start) {
		THROW_ERROR(E_INVALID, "taking slice of segments of size 0");
	}

	if (segs[dim].empty()) {
		build_order(dim);
	}

	ordered_segments *oseg = new ordered_segments(input);
	oseg->segs[dim].resize(end - start);
	segvec::iterator dstbox = oseg->segs[dim].begin();
	bitvec &dstbits = oseg->bits;
	for (ordered_iterator i = begin(dim) + start; i != begin(dim) + end; ++i) {
		*dstbox = *i;
		++dstbox;
		dstbits.set((*i)->pos, true);
	}

	return oseg;
}

ordered_segments &
ordered_segments::slice_inplace(size_t dim, size_t start, size_t end)
{
	if (segs[dim].empty()) {
		build_order(dim);
	}

	for (ordered_iterator i = begin(dim); i != begin(dim) + start; ++i) {
		bits.erase((*i)->pos);
	}
	for (ordered_iterator i = begin(dim) + end; i != this->end(dim); ++i) {
		bits.erase((*i)->pos);
	}

	segs[dim].erase(segs[dim].begin(), segs[dim].begin() + start);
	segs[dim].erase(segs[dim].begin() + end, segs[dim].end());

	for (segvec *j = segs; j != segs + NUM_ORDERS; ++j) {
		if (j != segs + dim) {
			j->clear();
		}
	}

	return *this;
}


void
ordered_segments::build_order(size_t d) const
{
	segvec &v = segs[d];
	for (size_t i = 0; i < bits.size(); ++i) { 
		if (bits.at(i)) {
			v.push_back(input[i]);
		}
	}
	std::sort(v.begin(), v.end(), orders[d]);
}


const Rect<long> &
ordered_segments::bounds() const {
	if (bbox.top == -1) {
		for (size_t i = 0; i < bits.size(); ++i) {
			if (bits.at(i)) {
				if (bbox.top == -1) {
					bbox = input[i]->bounds;
				}
				else {
					bbox = merge(bbox, input[i]->bounds);
				}
			}
		}
	}

	return bbox;
}


}
