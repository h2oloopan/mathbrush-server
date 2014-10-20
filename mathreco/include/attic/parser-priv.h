#ifndef PARSER_PRIV_H_
#define PARSER_PRIV_H_


#include "grammar.h"
#include "parser-fwd.h"
#include "parser.h"
#include "verb.h"

#include <list>
#include <queue>
#include <string>
#include <vector>


namespace scg
{


struct SpanData
{
	inline unsigned min(unsigned d) const { return dim_min[d]; }
	inline unsigned max(unsigned d) const { return dim_max[d]; }

	expression_box bounds;

	SpanData(unsigned D)
	{
		dim_min = new unsigned[D];
		dim_max = new unsigned[D];
		std::fill(dim_min, dim_min + D, 0);
		std::fill(dim_max, dim_max + D, 0);
	}

	~SpanData()
	{
		delete[] dim_min;
		delete[] dim_max;
	}

	unsigned *dim_min;
	unsigned *dim_max;
};


struct RangeProjectionData
{
	std::vector<std::vector<SpanData *> *> data;

	std::vector<unsigned> boxes;

	bool clean_boxes;

	inline SpanData &span(unsigned from, unsigned to) { return *(*data[from])[to-from]; }
	inline const SpanData &span(unsigned from, unsigned to) const { return *(*data[from])[to-from]; }

	RangeProjectionData() : clean_boxes(false) { }

	~RangeProjectionData() {
		for (std::vector<std::vector<SpanData *> *>::iterator i = data.begin(); i != data.end(); ++i) {
			for (std::vector<SpanData *>::iterator j = (*i)->begin(); j != (*i)->end(); ++j) {
				delete *j;
			}
			delete *i;
		}
	}
};


struct RangeData
{
	RangeProjectionData *data;
	std::map<int, ParseCell> cells;
	const CNFGrammar *g;

	bool clean_cells;

	RangeData(const CNFGrammar *g_, unsigned D, unsigned N) : g(g_), clean_cells(false)
	{
		data = new RangeProjectionData[D];
	}

	inline ParseCell *cell(unsigned gid)
	{
		std::map<int, ParseCell>::iterator i = cells.find(gid);
		return (i == cells.end()) ? 0 : &i->second;
	}

	inline const ParseCell *cell(unsigned gid) const
	{
		std::map<int, ParseCell>::const_iterator i = cells.find(gid);
		return (i == cells.end()) ? 0 : &i->second;
	}

	inline RangeProjectionData &projection(unsigned adim) { return data[adim]; }
	inline const RangeProjectionData &projection(unsigned adim) const { return data[adim]; }

	~RangeData() {
		delete[] data;
	}
};


struct DimData
{
	std::vector<unsigned> dim_to_parse_index;
	std::vector<unsigned> parse_to_dim_index;

	std::vector<std::vector<RangeData *> *> data;

	inline RangeData &range(unsigned from, unsigned to) { return *(*data[from])[to-from]; }
	inline const RangeData &range(unsigned from, unsigned to) const { return *(*data[from])[to-from]; }

	void clear() {
		for (std::vector<std::vector<RangeData *> *>::iterator i = data.begin(); i != data.end(); ++i) {
			for (std::vector<RangeData *>::iterator j = (*i)->begin(); j != (*i)->end(); ++j) {
				delete *j;
			}
			delete *i;
		}
	}
};

struct TerminalData
{
	std::map<int, double> scores;
};


struct DirtySpan
{
	unsigned start;     // start of super-range
	unsigned length; // length of super-range
	unsigned auxstart;  // start of parse range
	unsigned auxlength;    // length of parse range

	unsigned pdim;
	unsigned adim;

	bool operator<(const DirtySpan &rhs) const
	{
		if (auxlength < rhs.auxlength) {
			return true;
		}
		else if (auxlength == rhs.auxlength) {
			if (auxstart < rhs.auxstart) {
				return true;
			}
			else if (auxstart == rhs.auxstart) {
				if (length < rhs.length) {
					return true;
				}
				else if (length == rhs.length) {
					if (start < rhs.start) {
						return true;
					}
					else if (start == rhs.start) {
						if (adim < rhs.adim) {
							return true;
						}
						else if (adim == rhs.adim) {
							return pdim < rhs.pdim;
						}
					}
				}
			}
		}

		return false;
	}
};

std::ostream &operator<<(std::ostream &os, const DirtySpan &ds);


struct ParseContextPrivate
{
	std::vector<TerminalData> terminals;
	DimData *data;

	std::vector<std::set<DirtySpan> > dirty;

	void mark_span_dirty(const DirtySpan &ds);

	inline double terminal_score(unsigned gid, unsigned box)
	{
		std::map<int, double>::const_iterator i = terminals[box].scores.find(gid);
		if (i == terminals[box].scores.end()) {
			return 0.0;
		}
		return i->second;
	}

	inline void set_terminal_score(unsigned gid, unsigned box, double score) { terminals[box].scores[gid] = score; }

	inline DimData &dim(unsigned dim) { return data[dim]; }
	inline const DimData &dim(unsigned dim) const { return data[dim]; }
};


}



#endif
