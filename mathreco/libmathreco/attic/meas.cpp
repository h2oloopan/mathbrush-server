#include "meas.h"
#include "segment.h"
#include "parms.h"

namespace scg
{

long MinStrokeSize = 0;
void calc() {
	MinStrokeSize = static_cast<long>(TABLETPC_DPI * GetParameterDouble("MinStrokeSize"));
}
int _e = RegisterParameterCallback(&calc);

/*void
find_tl_point(const Bbox &box, long *pt)
{
	int tl_dist = std::numeric_limits<int>::max();
	
	for (RawStrokeGroup::const_iterator i = box.strokes.begin(); i != box.strokes.end(); ++i) {
		RawStroke::point_type *px = i->x;
		RawStroke::point_type *py = i->y;
		for (; px != i->x + i->npoints; px++,py++) {
			RawStroke::point_type x = *px - box.left;
			RawStroke::point_type y = *py - box.top;
			RawStroke::point_type d = x * x + y * y;
			if (d < tl_dist) {
				tl_dist = d;
				pt[0] = *px;
				pt[1] = *py;
			}
		}
	}
}

void
find_tr_point(const Bbox &box, long *pt)
{
	int tr_dist = std::numeric_limits<int>::max();

	for (RawStrokeGroup::const_iterator i = box.strokes.begin(); i != box.strokes.end(); ++i) {
		RawStroke::point_type *px = i->x;
		RawStroke::point_type *py = i->y;
		for (; px != i->x + i->npoints; px++,py++) {
			RawStroke::point_type x = box.right - *px;
			RawStroke::point_type y = *py - box.top;
			RawStroke::point_type d = x * x + y * y;
			if (d < tr_dist) {
				tr_dist = d;
				pt[0] = *px;
				pt[1] = *py;
			}
		}
	}
}

void
find_bl_point(const Bbox &box, long *pt)
{
	int bl_dist = std::numeric_limits<int>::max();

	for (RawStrokeGroup::const_iterator i = box.strokes.begin(); i != box.strokes.end(); ++i) {
		RawStroke::point_type *px = i->x;
		RawStroke::point_type *py = i->y;
		for (; px != i->x + i->npoints; px++,py++) {
			RawStroke::point_type x = *px - box.left;
			RawStroke::point_type y = box.bottom - *py;
			RawStroke::point_type d = x * x + y * y;
			if (d < bl_dist) {
				bl_dist = d;
				pt[0] = *px;
				pt[1] = *py;
			}
		}
	}
}

void
find_br_point(const Bbox &box, long *pt)
{
	int br_dist = std::numeric_limits<int>::max();
	
	for (RawStrokeGroup::const_iterator i = box.strokes.begin(); i != box.strokes.end(); ++i) {
		RawStroke::point_type *px = i->x;
		RawStroke::point_type *py = i->y;
		for (; px != i->x + i->npoints; px++,py++) {
			RawStroke::point_type x = box.right - *px;
			RawStroke::point_type y = box.bottom - *py;
			RawStroke::point_type d = x * x + y * y;
			if (d < br_dist) {
				br_dist = d;
				pt[0] = *px;
				pt[1] = *py;
			}
		}
	}
}
*/

LinkMeasurements
measure(const std::vector<Measurement *> &measurers, const Bbox &box1, const Bbox &box2)
{
	LinkMeasurements meas;
	for (std::vector<Measurement *>::const_iterator i = measurers.begin(); i != measurers.end(); ++i) {
		meas.push_back((*i)->measure(box1, box2));
	}

	return meas;
}


LinkMeasurements
measure(const std::map<std::string, Measurement *> &measurers, const Bbox &box1, const Bbox &box2)
{
	LinkMeasurements meas;
	for (std::map<std::string, Measurement *>::const_iterator i = measurers.begin(); i != measurers.end(); ++i) {
		meas.push_back(i->second->measure(box1, box2));
	}

	return meas;
}

}
