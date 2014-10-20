#ifndef RELATION_H_
#define RELATION_H_


#include "grammar-fwd.h"
#include "rect.h"
#include "links.h"

#include <string>
#include <vector>


namespace scg
{

typedef unsigned rclass_t;

extern const rclass_t INVALID_CLASS;
extern const rclass_t AGGREGATE_CLASS;
extern const rclass_t SYMBOL_CLASS;
extern const rclass_t BOX_CLASS;


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



typedef double (*get_angle_fn)(const Rect<long> &, const Rect<long> &, rclass_t, rclass_t);


struct relation
{
	relation(const std::string &name_, size_t n_, double t0_, double t1_, double t2_, get_angle_fn t_);

	std::string name;
	size_t n;
	double t0, t1, t2;
	
	get_angle_fn t;

	virtual double membership(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2) const;
	void add_sample(double t);
};

void initialize_relations();
void destroy_relations();

unsigned ordering_for_relation(const relation *rel);

relation *get_relation(int rel);

double get_angle_x(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2);
double get_AR_angle(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2);
double get_BR_angle(const Rect<long> &box1, const Rect<long> &box2, rclass_t class1, rclass_t class2);

}


#endif
