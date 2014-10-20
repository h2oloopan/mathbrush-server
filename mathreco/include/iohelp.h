#ifndef WRITER_H_
#define WRITER_H_

#include "binfmt.h"
#include <ostream>
#include <set>
#include <map>

namespace scg {

class interpreter;
class interpretation;
class basic_tree;
class bitvec;
class math_recognizer_base;
struct recoscore;

class reader {
public:
	explicit reader(std::istream &is_) : is(is_) { }
	bool isok() const;
	void read(interpreter **intrpr, math_recognizer_base *rec);
	void read(interpretation **intrp, math_recognizer_base *rec);
	void read(basic_tree **tree, math_recognizer_base *rec);
	void read(char &c);
	void read(size_t &i);
	void read(bitvec &bits);
	void read(recoscore &rs);
	void read(bool &b);
	void read(long &l);
	void read(double &d);
	void read(unsigned short &s);
	void read(std::string &s);
	void read(std::wstring &ws);
	void readptr(interpreter **intrpr);

private:
	std::istream &is;
	std::map<bin32_t, interpreter *> intrprs;
	std::map<bin32_t, interpretation *> intrps;
	std::map<bin32_t, basic_tree *> trees;
};


class writer {
public:
	explicit writer(std::ostream &os_) : os(os_) { }
	bool isok() const;
	void write(const interpreter *intrpr);
	void write(const interpretation *intrp);
	void write(const basic_tree *tree);
	void write(char c);
	void write(size_t i);
	void write(const bitvec &bits);
	void write(const recoscore &rs);
	void write(bool b);
	void write(long l);
	void write(double d);
	void write(unsigned short s);
	void write(const std::string &s);
	void write(const std::wstring &ws);
	void writeptr(const interpreter *intrpr);

private:
	std::ostream &os;
	std::set<const interpreter *> intrprs;
	std::set<const interpretation *> intrps;
	std::set<const basic_tree *> trees;
};



};


#endif
