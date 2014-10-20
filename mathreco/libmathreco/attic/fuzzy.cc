#include "fuzzy.h"

#include <algorithm>


namespace scg
{


double
cross_scores(double lhs, double rhs)
{
	double n = std::min(lhs, rhs);
	double x = std::max(lhs, rhs);
	return n * (1.0 + x * (x - n));
}


double
cross_scores_rel(double lhs, double rhs, double rel)
{
	//*verb_out << "cross " << lhs << " x " << rhs << " wrt " << rel << std::endl;
	//return std::min(rel, std::min(lhs, rhs));

	// P --> A0 r1 A1 r2 ... rn An
	// P~ = ((((((A0~ x A1~) n r1~) x A2~) n r2~) ... x An~) n rn~)
	// But we have P --> A r B
	// so P~ = (A~ x B~) n r~

	//Yager's intersection
	//double c = 1.0 - std::min(1.0, std::pow(std::pow(1.0 - cross, p) + std::pow(1.0 - rel, p), 1.0/p));
	return cross_scores(cross_scores(lhs, rhs), rel);
}



}
