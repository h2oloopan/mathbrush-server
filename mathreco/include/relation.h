#ifndef RELATION_H_
#define RELATION_H_


#include "grammar-fwd.h"
#include "rect.h"

#include <string>
#include <vector>


namespace scg
{

class math_recognizer_base;
typedef size_t rclass_t;

extern const rclass_t INVALID_CLASS;
extern const rclass_t AGGREGATE_CLASS;
extern const rclass_t SYMBOL_CLASS;
extern const rclass_t BOX_CLASS;
extern const rclass_t RELCLASS_MERGE;
extern const rclass_t NCLASSES;


class interpretation;


enum {
	REL_ABOVERIGHT,
	REL_RIGHT,
	REL_BELOWRIGHT,
	REL_BELOW,
	REL_CONTAINS,
	//REL_NONE,
	NUM_RELATIONS
};



typedef double (*get_angle_fn)(const interpretation *, const interpretation *, rclass_t, rclass_t);

struct relation {
	std::string name;

	relation(const std::string &name_) : name(name_) { }
	virtual ~relation() { }
	virtual double membership(const interpretation *tail, const interpretation *head, rclass_t class1, rclass_t class2, int tatt, int hatt) const = 0;
};


int initialize_relations();
void destroy_relations();

unsigned ordering_for_relation(const relation *rel);

relation *get_relation(int rel);
relation *get_group_relation(int rel);
size_t relindex(const relation *rel);

rclass_t tag_to_rclass(const std::string &tag);
std::string rclass_to_tag(rclass_t cls);
std::string rclass_to_str(rclass_t cls);

void set_relation(int rel, relation *r);
double nominal_height(rclass_t c);
double getbl(const Rect<long> &bbox, rclass_t c);

double get_nil_probability(const interpretation *tail, const interpretation *head);

}


#endif
