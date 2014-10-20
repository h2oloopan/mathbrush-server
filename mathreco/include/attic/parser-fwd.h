#ifndef PARSER_FWD_H_
#define PARSER_FWD_H_


#include "grammar-fwd.h"


namespace scg
{


class ordered_segments;


struct parseref {
	ordered_segments *span;
	const nonterminal *nt;

	parseref() : span(0), nt(0) { }
	parseref(ordered_segments *span_, const nonterminal *nt_) : span(span_), nt(nt_) { }

	inline bool operator==(const parseref &rhs) const
		{ return span == rhs.span && nt == rhs.nt; }
	inline bool operator!=(const parseref &rhs) const
		{ return span != rhs.span || nt != rhs.nt; }
};


struct parse_cell;
struct context;

struct segment;
struct group;
struct stroke;


}


#endif
