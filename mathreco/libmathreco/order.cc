#include "order.h"
#include "segment.h"
//#include "mathrecognizer-private.h"

namespace scg
{


const unsigned X_ORDER = 0;
const unsigned Y_ORDER = 1;
const unsigned TIME_ORDER = 2;

/*static bool
r_order(const segment *lhs, const segment *rhs) {
	const segment *base = lhs->ctx->segments[0];
	long basey = (base->bounds.top + base->bounds.bottom) / 2;
	long lhsy = (lhs->bounds.top + lhs->bounds.bottom) / 2;
	long rhsy = (rhs->bounds.top + rhs->bounds.bottom) / 2;
	long dx = lhs->bounds.left - base->bounds.left;
	long dy = lhsy - basey;
	long lhsr = dx*dx + dy*dy;
	dx = rhs->bounds.left - base->bounds.left;
	dy = rhsy - basey;
	long rhsr = dx*dx + dy*dy;
	return lhsr < rhsr;
}

static bool
th_order(const segment *lhs, const segment *rhs) {
	const segment *base = lhs->ctx->segments[0];
	long basey = (base->bounds.top + base->bounds.bottom) / 2;
	long lhsy = (lhs->bounds.top + lhs->bounds.bottom) / 2;
	long rhsy = (rhs->bounds.top + rhs->bounds.bottom) / 2;
	long dx = lhs->bounds.left - base->bounds.left;
	long dy = lhsy - basey;
	double lhsth = std::atan2((double)dy,dx);
	dx = rhs->bounds.left - base->bounds.left;
	dy = rhsy - basey;
	double rhsth = std::atan2((double)dy,dx);
	return lhsth < rhsth;
}*/

static bool
x_order(const segment *lhs, const segment *rhs)
	//{ return (3*lhs->bounds.left + lhs->stk.x_center()) < (3*rhs->bounds.left + rhs->stk.x_center()); }
	//{ return (7*lhs->bounds.left + lhs->bounds.right) < (7*rhs->bounds.left + rhs->bounds.right); }
	{ return lhs->bounds.left < rhs->bounds.left; }
	/*{
		double a = std::max(0.0, 0.1*(1 - lhs->bounds.width()/TABLETPC_DPI));
		double lp = (lhs->bounds.left + a*lhs->bounds.right)/(1+a);
		a = std::max(0.0, 0.1*(1 - rhs->bounds.width()/TABLETPC_DPI));
		double rp = (rhs->bounds.left + a*rhs->bounds.right)/(1+a);
		return lp < rp;
	}*/

static bool
y_order(const segment *lhs, const segment *rhs)
	//{ return (3*lhs->bounds.top + lhs->stk.y_center()) < (3*rhs->bounds.top + rhs->stk.y_center()); }
	//{ return (7*lhs->bounds.top + lhs->bounds.bottom) < (7*rhs->bounds.top + rhs->bounds.bottom); }
	{ return lhs->bounds.top < rhs->bounds.top; }
	/*{
		double a = std::max(0.0, 0.1*(1 - lhs->bounds.height()/TABLETPC_DPI));
		double lp = (lhs->bounds.top + a*lhs->bounds.bottom)/(1+a);
		a = std::max(0.0, 0.1*(1 - rhs->bounds.height()/TABLETPC_DPI));
		double rp = (rhs->bounds.top + a*rhs->bounds.bottom)/(1+a);
		return lp < rp;
	}*/

static bool
t_order(const segment *lhs, const segment *rhs)
	{ return lhs->pos < rhs->pos; }

order_fn orders[NUM_ORDERS] = { x_order, y_order, t_order };


}
