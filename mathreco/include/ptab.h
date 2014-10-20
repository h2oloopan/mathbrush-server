#ifndef PTAB_H_
#define PTAB_H_

#include "grammar-fwd.h"
#include "grammar-values.h"
#include <istream>
#include <ostream>
#include <map>
#include <vector>
#include <string>

namespace scg {

void ptab_staticinit();
void ptab_rectify(const grammar *G);

struct ExpressionTree;
class interpretation;

struct typetree {
	SemanticId sid;
	std::vector<typetree> children;
	typetree() : sid(InvalidSemanticId) { }
	explicit typetree(SemanticId sid_) : sid(sid_) { }
};


int mktypetree(typetree &tt, std::istream &is);
int mktypetree(typetree &tt, const ExpressionTree *tree);
std::ostream &operator<<(std::ostream &os, const typetree &tt);

struct ptab;
ptab *mkptab(const grammar *G);
void rmptab(ptab *pt);
int ptab_insert(ptab *pt, const typetree &tt);
void ptab_remove(ptab *pt, SemanticId sid);
int ptab_read(ptab *pt, std::istream &is);
void ptab_write(ptab *pt, std::ostream &os);
double ptab_measurectx(ptab *pt, const typetree &tt);
double ptab_measuretop(const ptab *pt, const ExpressionTree *tree);
double ptab_measure(ptab *pt, const typetree &tt);
double ptab_measuretop_log(ptab *pt, const typetree &tt);
double ptab_measure_log(ptab *pt, const typetree &tt);

}

#endif
