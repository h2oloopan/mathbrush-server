#include "binfmt.h" // this must be the first include so windows' namespace is clean
#include "iohelp.h"
#include "intrpr.h"
#include "expr-iter.h"
#include "expr-node.h"
#include "verb.h"
#include <ostream>
#include <numeric>
#include <cstddef>

#ifdef WIN32
#include <float.h>
#endif

namespace scg {

bool writer::isok() const { return (bool)os; }
bool reader::isok() const { return (bool)is; }

static bin32_t
ptr32(const void *ptr) {
	if (sizeof(void *) == 8) {
		return (bin32_t)((ptrdiff_t)ptr & 0xffffffff);
	}
	else {
		return (bin32_t)(ptrdiff_t)ptr;
	}
}

void
writer::writeptr(const interpreter *intrpr) {
	if (intrprs.find(intrpr) != intrprs.end()) {
		bin32_t id = htonl(ptr32(intrpr));
		os.write((char *)&id, sizeof(id));
	}
	else {
		// XXX: this is nasty
		bin32_t unk = htonl((bin32_t)1);
		os.write((char *)&unk, sizeof(unk));
	}
}

void
writer::write(const interpreter *intrpr) {
	bin32_t id = htonl(ptr32(intrpr));
	os.write((char *)&id, sizeof(id));
	if (intrpr) {
		std::set<const interpreter *>::iterator i = intrprs.find(intrpr);
		if (i == intrprs.end()) {
			//static unsigned N = 0;
			//VERBOSE(*verb_out << "writei " << id << " for " << N++ << std::endl);
			intrprs.insert(i, intrpr);
			write(intrpr->id());
			intrpr->write(*this);
		}
	}
}

void
reader::readptr(interpreter **intrpr) {
	bin32_t id;
	is.read((char *)&id, sizeof(id));
	id = ntohl(id);
	if (!id) {
		*intrpr = 0;
	}
	else if (id == 1) {
		*intrpr = (interpreter *)1;
	}
	else {
		std::map<bin32_t, interpreter *>::iterator i = intrprs.find(id);
		assert(i != intrprs.end());
		*intrpr = i->second;
	}
}

void
reader::read(interpreter **intrpr, math_recognizer_base *rec) {
	bin32_t id;
	is.read((char *)&id, sizeof(id));
	id = ntohl(id);
	if (!id) {
		*intrpr = 0;
	}
	else {
		interpreter *&newintrpr = intrprs[id];
		char typ;
		if (!newintrpr) {
			//static unsigned N = 0;
			read(typ);
			//VERBOSE(*verb_out << "readi " << htonl(id) << " for " << N++ << std::endl);
			if (typ == staticintrpr::ID) {
				newintrpr = new staticintrpr;
			}
			/*else if (typ == iterator::ID) {
				newintrpr = new iterator;
			}*/
			else if (typ == multiplexor::ID) {
				newintrpr = new multiplexor;
			}
			else if (typ == treeparser::ID) {
				newintrpr = new treeparser;
			}
			else if (typ == matrixparser::ID) {
				newintrpr = new matrixparser;
			}
			else {
				THROW_ERROR(E_INVALID, "read invalid interpreter " << id);
			}
			newintrpr->read(*this, rec);
		}
		*intrpr = newintrpr;
	}
}

void
writer::write(const interpretation *intrp) {
	bin32_t id = htonl(ptr32(intrp));
	os.write((char *)&id, sizeof(id));
	if (intrp) {
		std::set<const interpretation *>::iterator i = intrps.find(intrp);
		if (i == intrps.end()) {
			intrps.insert(i, intrp);
			write(intrp->id());
			intrp->write(*this);
		}
	}
}

void
reader::read(interpretation **intrp, math_recognizer_base *rec) {
	bin32_t id;
	is.read((char *)&id, sizeof(id));
	id = ntohl(id);
	if (!id) {
		*intrp = 0;
	}
	else {
		interpretation *&newintrp = intrps[id];
		if (!newintrp) {
			char typ;
			read(typ);
			if (typ == interpretation::ID) {
				newintrp = new interpretation(*this, rec);
			}
			else if (typ == fixed_interpretation::ID) {
				newintrp = new fixed_interpretation(*this, rec);
			}
			else if (typ == matrixintrp::ID) {
				newintrp = new matrixintrp(*this, rec);
			}
			else {
				THROW_ERROR(E_INVALID, "read invalid interpreter " << id);
			}
		}
		*intrp = newintrp;
	}
}

void
writer::write(const basic_tree *tree) {
	bin32_t id = htonl(ptr32(tree));
	os.write((char *)&id, sizeof(id));
	if (tree) {
		std::set<const basic_tree *>::iterator i = trees.find(tree);
		if (i == trees.end()) {
			trees.insert(i, tree);
			write(tree->id());
			tree->write(*this);
		}
	}
}

void
reader::read(basic_tree **tree, math_recognizer_base *rec) {
	bin32_t id;
	is.read((char *)&id, sizeof(id));
	id = ntohl(id);
	if (!id) {
		*tree = 0;
	}
	else {
		basic_tree *&newtree = trees[id];
		if (!newtree) {
			char typ;
			read(typ);
			if (typ == parsed_tree::ID) {
				newtree = new parsed_tree(*this, rec);
			}
			else if (typ == invented_tree::ID) {
				newtree = new invented_tree(*this, rec);
			}
			else if (typ == terminal_tree::ID) {
				newtree = new terminal_tree(*this, rec);
			}
			else {
				THROW_ERROR(E_INVALID, "read invalid interpreter " << id);
			}

		}
		*tree = newtree;
	}
}

void
writer::write(char c) {
	os.write(&c, 1);
}

void
reader::read(char &c) {
	is.read(&c, 1);
}

void
writer::write(size_t i) {
	write((long)i);
}

void
reader::read(size_t &i) {
	read((long &)i);
}

void
writer::write(const bitvec &bits) {
	bits.write(*this);
}

void
reader::read(bitvec &bits) {
	bits.read(*this);
}

void
writer::write(const recoscore &rs) {
	write(rs.hint);
	write(rs.score);
}

void
reader::read(recoscore &rs) {
	read(rs.hint);
	read(rs.score);
}

void
writer::write(bool b) {
	write((char)(b ? 1 : 0));
}

void
reader::read(bool &b) {
	char c;
	read(c);
	b = (c != 0);
}

void
writer::write(long l) {
	bin32_t c = htonl((bin32_t)(l & 0xffffffff));
	os.write((char *)&c, sizeof(c));
}

void
reader::read(long &l) {
	bin32_t c;
	is.read((char *)&c, sizeof(c));
	c = ntohl(c);
	l = c;
}

void
writer::write(unsigned short s) {
	short c = htons((short)s);
	os.write((char *)&c, sizeof(c));
}

void
reader::read(unsigned short &s) {
	is.read((char *)&s, sizeof(s));
	s = ntohs(s);
}

void
writer::write(double d) {
	char s;
#ifdef WIN32
	switch (_fpclass(d)) {
	case _FPCLASS_NZ:
	case _FPCLASS_PZ:
		s = 0;
		break;
	case _FPCLASS_SNAN:
	case _FPCLASS_QNAN:
		s = 1;
		break;
	case _FPCLASS_PINF:
		s = 2;
		break;
	case _FPCLASS_NINF:
		s = 3;
		break;
	case _FPCLASS_PD:
		s = 4;
		break;
	case _FPCLASS_ND:
		s = 5;
		break;
	default:
		s = 6;
		break;
	}
#else
	switch (std::fpclassify(d)) {
	case FP_ZERO:
		s = 0;
		break;
	case FP_NAN:
		s = 1;
		break;
	case FP_INFINITE:
		if (d > 0) {
			s = 2;
		}
		else {
			s = 3;
		}
		break;
	case FP_SUBNORMAL:
		if (d > 0) {
			s = 4;
		}
		else {
			s = 5;
		}
		break;
	case FP_NORMAL:
		s = 6;
		break;
	}
#endif
	os.write(&s, 1);
	if (s == 6) {
		int exp;
		double fr = std::frexp(d, &exp);
		bin32_t b = htonl((bin32_t)exp);
		os.write((char *)&b, sizeof(b));
		b = htonl((bin32_t)((int)std::ldexp(fr, 31)));
		os.write((char *)&b, sizeof(b));
	}
}

void
reader::read(double &d) {
	bin32_t v;
	int exp;
	int fr;
	char t;
	is.read(&t, 1);
	switch (t) {
	case 0:
		d = 0.0;
		break;
	case 1:
		d = std::numeric_limits<double>::quiet_NaN();
		break;
	case 2:
		d = std::numeric_limits<double>::infinity();
		break;
	case 3:
		d = -std::numeric_limits<double>::infinity();
		break;
	case 4:
		d = std::numeric_limits<double>::denorm_min();
		break;
	case 5:
		d = -std::numeric_limits<double>::denorm_min();
		break;
	case 6:
		is.read((char *)&v, sizeof(v));
		exp = (int)ntohl(v);
		is.read((char *)&v, sizeof(v));
		fr = (int)ntohl(v);
		d = std::ldexp((double)fr, -31);
		d = std::ldexp(d, exp);
		break;
	default:
		abort();
	}
}

void
writer::write(const std::string &s) {
	write(s.length());
	os.write(&s[0], s.length());
}

void
reader::read(std::string &s) {
	size_t n;
	read(n);
	s.resize(n);
	is.read(&s[0], n);
	s[n] = '\0';
}

void
writer::write(const std::wstring &ws) {
	write(ws.length());
	for (size_t i = 0; i < ws.length(); ++i) {
		write((unsigned short)ws[i]);
	}
}

void
reader::read(std::wstring &ws) {
	size_t n;
	read(n);
	ws.resize(n);
	for (size_t i = 0; i < n; ++i) {
		read((unsigned short &)ws[i]);
	}
	ws[n] = L'\0';
}

}
