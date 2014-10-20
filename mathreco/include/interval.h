#ifndef SCG_INTERVAL_H
#define SCG_INTERVAL_H

namespace scg
{


// Defines an arbitrary interval having a start and an end value.
struct Interval {
	Interval() : start(-1), end(-1) {}
	Interval(int start, int end) : start(start), end(end) {}
	int start, end;

	inline size_t length() const { return static_cast<size_t>(std::max(end - start, 1)); }

	inline bool operator==(const Interval& other) const { return this->start == other.start; }
	inline bool operator<=(const Interval& other) const { return this->start <= other.start; }
	inline bool operator<(const Interval& other) const { return this->start < other.start; }
};

inline Interval intervalIntersection(const Interval& i1, const Interval& i2)
{
	if (i1.start < 0 && i2.start >= 0) {
		return i2;
	} else if (i1.start >= 0 && i2.start < 0) {
		return i1;
	} else if (i1.start < 0 && i2.start < 0) {
		return Interval();
	}

	unsigned intersectStart = std::max(i1.start, i2.start);
	unsigned intersectEnd = std::min(i1.end, i2.end);
	return Interval(intersectStart, intersectEnd);
}

inline Interval intervalUnion(const Interval& i1, const Interval& i2)
{
	if (i1.start < 0 && i2.start >= 0) {
		return i2;
	} else if (i1.start >= 0 && i2.start < 0) {
		return i1;
	} else if (i1.start < 0 && i2.start < 0) {
		return Interval();
	}

	unsigned unionStart = std::min(i1.start, i2.start);
	unsigned unionEnd = std::max(i1.end, i2.end);
	return Interval(unionStart, unionEnd);
}

inline bool intervalsOverlap(const Interval& int1, const Interval& int2)
{
	return (int1.end >= int2.start && int1.start <= int2.end);
}

inline double intervalOverlapProportion(const Interval& i1, const Interval& i2)
{
    if (!intervalsOverlap(i1, i2)) { return 0.0; }

    Interval intersection = intervalIntersection(i1, i2);
    return (double) intersection.length() / std::min(i1.length(), i2.length());
}

inline bool intervalsOverlapCompletely(const Interval& i1, const Interval& i2)
{
	return intervalOverlapProportion(i1, i2) == 1.0;
}

inline Interval mergeIntervals(const Interval& int1, const Interval& int2)
{
	unsigned minStart = std::min(int1.start, int2.start);
	unsigned maxEnd = std::max(int1.end, int2.end);
	return Interval(minStart, maxEnd);
}

}

#endif