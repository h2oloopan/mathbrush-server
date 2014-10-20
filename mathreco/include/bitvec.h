#ifndef BITVEC_H_
#define BITVEC_H_


#include "iohelp.h"
#include "error.h"
#include "verb.h"

#include <ostream>
#include <cassert>

#include <iostream>

namespace scg
{


class bitvec {
private:
	struct bit_proxy {
		bit_proxy(int *b_, int m_) : b(b_), m(m_) { }

		operator bool() const { return (*b & m) != 0; }

		inline bit_proxy &operator=(bool v)
		{
			if (v) {
				*b = *b | m;
			}
			else {
				*b = *b & ~m;
			}
			return *this;
		}

	private:
		int *b;
		int m;
	};


public:
	bitvec();
	explicit bitvec(size_t N, bool v = true);
	bitvec(const bitvec &vec);

	~bitvec();


	bitvec &operator=(const bitvec &rhs);


	inline size_t size() const
		{ return n; }

	inline int highest_set_bit() const {
		int base = (int)(cap - 32);
		for (int i = (int)(cap / 32 - 1); i >= 0; --i) {
			int word = bits[i];
			if (word != 0) {
				for (int j = 31; j >= 0; --j) {
					if (word & (1 << j)) {
						return base + j;
					}
				}
			}
			base -= 32;
		}
		return -1;
	}


	size_t count_set_bits() const;

	void write(writer &wr) const;
	void read(reader &re);

	inline bool operator[](size_t i) const
		{ return bit(i >> 5, (1 << (i & 0x1f))); }

	inline bit_proxy operator[](size_t i)
		{ return bit_proxy(bits + (i / 32), 1 << (i % 32)); }

	inline bool at(size_t i) const
	{
		if (i >= n) {
			return false;
			//THROW_ERROR(E_INVALID, "index " << i << " too large for bitvec");
		}
		return bit(i / 32, i % 32);
	}


	inline void clear()
	{
		delete[] bits;
		bits = 0;
		n = 0;
		cap = 0;
	}
	
	inline void union_insert(const bitvec &rhs)
	{
		if (n < rhs.n) {
			THROW_ERROR(E_INVALID, "size mismatch in unioned bitvectors");
		}
		for (int *b1 = bits, *b2 = rhs.bits; b2 != rhs.end(); ++b1, ++b2) {
			*b1 = (*b1 | *b2);
		}
	}

	inline bool subset_of(const bitvec &rhs) const {
		const int *b1 = bits;
		const int *b2 = rhs.bits;
		const int *end1 = end();
		const int *end2 = rhs.end();
		for (; b1 != end1 && b2 != end2; ++b1, ++b2) {
			if ((*b1 | *b2) != *b2) return false;
		}
		for (; b1 != end1; ++b1) {
			if (*b1) return false;
		}
		return true;
	}


	inline bool set(size_t i, bool v = true)
	{
		if (i >= n) {
			THROW_ERROR(E_INVALID, "index " << i << " too large for bitvec");
		}
		int *b = bits + (i / 32);
		if (v) {
			*b = *b | (1 << (i % 32));
		}
		else {
			*b = *b & ~(1 << (i % 32));
		}
		return v;
	}


	inline void insert(size_t i, bool v)
	{
		if (n + 1 == cap) {
			reallocate(cap * 2);
		}
		
		push_bits_down(v, i / 32, i % 32);

		++n;
	}


	inline void erase(size_t i)
	{
		pull_bits_up(i / 32, i % 32);

		--n;
		if (2 * n < cap) {
			reallocate(cap / 2);
		}
	}


	inline bool operator==(const bitvec &vec) const
	{
			return equals(vec);
		if (size() <= vec.size()) {
			return equals(vec);
		}
		else {
			return vec.equals(*this);
		}
	}
	
public:
	static const bitvec EMPTY;

private:
	inline bool equals(const bitvec &vec) const {
		int *b1 = bits, *b2 = vec.bits;
		int *end1 = end(), *end2 = vec.end();
		for ( ; b1 != end1 && b2 != end2; ++b1, ++b2) {
			if (*b1 != *b2) {
				return false;
			}
		}

		for ( ; b1 != end1; ++b1) {
			if (*b1) return false;
		}
		for ( ; b2 != end2; ++b2) {
			if (*b2) return false;
		}

		return true;
	}


public:
	inline bool operator!=(const bitvec &vec) const
		{ return !(*this == vec); }

	bool operator<(const bitvec &vec) const {
		return less(vec);
		/*
		bool eq = size() <= vec.size() ? equals(vec) : vec.equals(*this);
		if (eq) return false;
		if (subset_of(vec)) return true;
		if (vec.subset_of(*this)) return false;*/

		//if (size() < vec.size()) return true;
		//if (size() > vec.size()) return false;
		//return less(vec);
		if (size() <= vec.size()) {
			//VERBOSE(*verb_out << *this << " < " << vec << " = " << less(vec) << std::endl);
			return less(vec);
		}
		else {
			//VERBOSE(*verb_out << *this << " < " << vec << " = " << (!vec.less(*this) && !vec.equals(*this)) << std::endl);
			return !vec.less(*this) && !vec.equals(*this);
		}
	}

private:
	inline bool less(const bitvec &vec) const {
		int *b1 = bits, *b2 = vec.bits;
		int *end1 = end(), *end2 = vec.end();
		for ( ; b1 != end1 && b2 != end2; ++b1, ++b2) {
			if (*b1 > *b2) return false;
			if (*b1 < *b2) return true;
		}

		for ( ; b1 != end1; ++b1) {
			if (*b1) return false;
		}
		for ( ; b2 != end2; ++b2) {
			if (*b2) return true;
		}

		return false;
	}


private:
	inline bool bit(size_t i, size_t j) const
		{ return (bits[i] & (1 << j)) != 0; }


	void reallocate(size_t N);

	void push_bits_down(bool v, size_t i, size_t j);
	void pull_bits_up(size_t i, size_t j);


private:
	inline int * const end() const
		{ return bits + n/32 + 1; }
	inline int *end()
		{ return bits + n/32 + 1; }

public:
	int *bits;
	size_t n;
	size_t cap;

private:
	friend std::ostream &operator<<(std::ostream &os, const bitvec &vec);
	friend std::istream &operator>>(std::istream &is, bitvec &vec);
	friend bool null_intersection(const bitvec &lhs, const bitvec &rhs);
	friend bitvec bitvec_union(const bitvec &lhs, const bitvec &rhs);
	
#ifdef WIN32
	static const unsigned char lookup[256];
#endif
};


std::ostream &operator<<(std::ostream &os, const bitvec &vec);
std::istream &operator>>(std::istream &is, bitvec &vec);

bool null_intersection(const bitvec &lhs, const bitvec &rhs);
bitvec bitvec_union(const bitvec &lhs, const bitvec &rhs);


}


#endif
