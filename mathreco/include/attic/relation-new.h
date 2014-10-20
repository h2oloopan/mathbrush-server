#ifndef RELATION_H_
#define RELATION_H_


#include "grammar-fwd.h"
#include "rect.h"
#include "links.h"

#include <string>
#include <vector>


namespace scg
{

class ordered_segments;
typedef unsigned rclass_t;

extern const rclass_t INVALID_CLASS;
extern const rclass_t AGGREGATE_CLASS;
extern const rclass_t SYMBOL_CLASS;
extern const rclass_t BOX_CLASS;
extern const rclass_t BASELINE_CLASS;
extern const rclass_t ASCENDER_CLASS;
extern const rclass_t DESCENDER_CLASS;
extern const rclass_t EXTENDER_CLASS;
extern const rclass_t HALF_ASCENDER_CLASS;
extern const rclass_t HALF_DESCENDER_CLASS;
extern const rclass_t HALF_ASCENDER_DESCENDER_CLASS;
extern const rclass_t CENTERED_CLASS;
extern const rclass_t FENCE_CLASS;
extern const rclass_t LARGE_EXTENDER_CLASS;
extern const rclass_t RANGE_OP_CLASS;
extern const rclass_t ROOT_CLASS;
extern const rclass_t HORIZONTAL_CLASS;
extern const rclass_t PUNCTUATION_CLASS;
extern const size_t NUM_RCLASSES;
extern const size_t NUM_SPECIAL_RCLASSES;


struct expression_box;


enum {
	REL_ABOVERIGHT,
	REL_RIGHT,
	REL_BELOWRIGHT,
	REL_BELOW,
	REL_CONTAINS,
	//REL_NONE,
	NUM_RELATIONS
};


struct classdata {
	size_t n;
	double dmu;
	double dlambda;
	double txmean;
	double tymean;
	double tvar;
	bool duseful;
	bool tuseful;
	classdata() : n(0), dmu(0.0), dlambda(0.0), txmean(0.0), tymean(0.0), tvar(0.0), duseful(false), tuseful(false) { }
};


struct relation
{
	relation(const std::string &name_) : name(name_) { }

	std::string name;
	//mutable std::map<std::pair<rclass_t, rclass_t>, classdata> data;
	virtual double membership(const ordered_segments *from, const ordered_segments *to, rclass_t class1, rclass_t class2) const = 0;
};

void initialize_relations();
void destroy_relations();

unsigned ordering_for_relation(const relation *rel);

relation *get_relation(int rel);

}


#endif
