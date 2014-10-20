#ifndef ITER_TOOLS_H_
#define ITER_TOOLS_H_


namespace scg
{


template <typename T>
T
inc_iter(T it, int n)
{
	if (n < 0) {
		while (n++ != 0) {
			--it;
		}
	}
	else if (n > 0) {
		while (n-- != 0) {
			++it;
		}
	}
	return it;
}


}


#endif

