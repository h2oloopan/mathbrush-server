#include "bitvec.h"
#include "stream-defs.h"
#include "binfmt.h"

#include <algorithm>
#include <cstring>


namespace scg
{

const bitvec bitvec::EMPTY(0);

bitvec::bitvec() : bits(0), n(0), cap(0) { }

bitvec::bitvec(size_t N, bool v) : bits(0), n(0), cap(0)
{
	reallocate(N);
	std::fill(bits, bits + cap / 32, v ? ~0 : 0);
	int m = (1 << (N % 32)) - 1;
	*(bits + cap / 32 - 1) &= m;
	n = N;
}

bitvec::bitvec(const bitvec &vec) : bits(0), n(0), cap(0)
{
	reallocate(vec.size());
	std::copy(vec.bits, vec.end(), bits);
	n = vec.size();
}

bitvec &
bitvec::operator=(const bitvec &rhs)
{
	if (this != &rhs) {
		delete[] bits;
		bits = 0;
		reallocate(rhs.size());
		std::copy(rhs.bits, rhs.end(), bits);
		n = rhs.size();
	}
	return *this;
}


bitvec::~bitvec()
{ clear(); }

void
bitvec::write(writer &wr) const {
	wr.write(n);
	wr.write(cap);
	for (int *i = bits; i != bits + cap / 32; ++i) {
		wr.write((long)*i);
	}
}

void
bitvec::read(reader &re) {
	clear();
	re.read(n);
	re.read(cap);
	assert(cap % 32 == 0);
	bits = new int[cap / 32];
	for (int *i = bits; i != bits + cap / 32; ++i) {
		long prox;
		re.read(prox);
		*i = (int)prox;
	}
}

#ifdef WIN32
const unsigned char bitvec::lookup[256] =
  { 0, 1, 1, 2, 1, 2, 2, 3,
    1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    
    1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    
    2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7,
    5, 6, 6, 7, 6, 7, 7, 8 };
#endif

size_t
bitvec::count_set_bits() const
{
	if (n == 0) {
		return 0;
	}

	size_t c = 0;
	for (int *b = bits; b != end(); ++b) {
#ifdef WIN32
		register int t = *b;
		c += lookup[t & 0xff] + lookup[(t >> 8) & 0xff] + lookup[(t >> 16) & 0xff] + lookup[(t >> 24) & 0xff];
#else
		c += __builtin_popcount(*b);
#endif
	}
	/*

	int b = *(bits + (n / 32));
	int fm = 1 << (n % 32);
	for (int m = 1; m < fm; m <<= 1) {
		if (b & m) {
			++c;
		}
	}
	*/

	return c;
}


void
bitvec::reallocate(size_t N)
{
	//if (N > 0 && N % 32 == 0) --N;

	size_t newcap = (N / 32) + 1;
	int *newbits = new int[newcap];
	std::fill(newbits, newbits + newcap, 0);
	if (bits) {
		std::copy(bits, bits + std::min(cap, 32 * newcap) / 32, newbits);
		delete[] bits;
	}
	bits = newbits;
	cap = newcap * 32;
}


void
bitvec::push_bits_down(bool v, size_t i, size_t j)
{
	if ((i+1) * 32 <= n) {
		push_bits_down(bit(i, 31), i + 1, 0);
	}

	int *b = bits + i;
	int m = (1 << j) - 1;
	int m2 = ~m ^ (1 << j);
	*b = (*b & m) | ((*b & m2) << 1);
	if (v) {
		*b = *b | (1 << j);
	}
	else {
		*b = *b & ~(1 << j);
	}
}


void
bitvec::pull_bits_up(size_t i, size_t j)
{
	int *b = bits + i;
	int m = (1 << j) - 1;    // mask for bits lower than j
	int m2 = ~m ^ (1 << j);  // mask for bits higher than j
	*b = (*b & m) | ((*b & m2) >> 1);

	if ((i+1) * 32 <= n) {
		bool v = bit(i + 1, 0);
		if (v) {
			*b = *b | (1 << 31);
		}
		else {
			*b = *b & ~(1 << 31);
		}

		pull_bits_up(i + 1, 0);
	}
}


std::ostream &
operator<<(std::ostream &os, const bitvec &vec)
{
	size_t n = 0;
	for (int *b = vec.bits; b != vec.end(); ++b) {
		for (int m = 1; m != 0 && ++n <= vec.size(); m <<= 1) {
			os << ((*b & m) ? '1' : '0');
		}
	}
	return os;
}

std::istream &
operator>>(std::istream &is, bitvec &vec)
{
	vec.clear();
	
	std::string s;
	is >> s;
	CHECK_ISTREAM(is, "while reading bitvec");
	vec.reallocate(s.length());
	vec.n = s.length();

	for (size_t i = 0; i < s.length(); ++i) {
		if (s[i] == '1') {
			vec.set(i, true);
		}
		else {
			assert(s[i] == '0');
			vec.set(i, false);
		}
	}
	return is;
}


bool
null_intersection(const bitvec &lhs, const bitvec &rhs)
{
	int *b1, *b2;
	for (b1 = lhs.bits, b2 = rhs.bits; b1 != lhs.end() && b2 != rhs.end(); ++b1, ++b2) {
		if (*b1 & *b2) {
			return false;
		}
	}
	if (b1 == lhs.end() && b2 == rhs.end()) {
		return true;
	}
	else if (b1 == lhs.end()) {
		for (; b2 != rhs.end(); ++b2) {
			if (*b2) {
				return false;
			}
		}
	}
	else if (b2 == rhs.end()) {
		for (; b1 != lhs.end(); ++b1) {
			if (*b1) {
				return false;
			}
		}
	}
	return true;
}


bitvec
bitvec_union(const bitvec &lhs, const bitvec &rhs) {
	bitvec U(std::max(lhs.cap, rhs.cap), false);
	size_t lsz = lhs.cap / 32;
	size_t rsz = rhs.cap / 32;
	int *p = U.bits;
	int *lp, *rp;
	for (lp = lhs.bits, rp = rhs.bits; lp != lhs.bits + lsz && rp != rhs.bits + rsz; ++lp, ++rp, ++p) {
		*p = *lp | *rp;
	}
	if (lsz > rsz) {
		memcpy(p, lp, lsz - rsz);
	}
	else if (rsz > lsz) {
		memcpy(p, rp, rsz - lsz);
	}
	return U;
}

}
