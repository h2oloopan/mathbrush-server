#ifndef ANNOTATE_H_
#define ANNOTATE_H_


#include "group.h"
#include "grammar.h"

#include <algorithm>
#include <istream>
#include <ostream>
#include <list>
#include <string>


namespace scg
{


typedef std::list<size_t> stroke_collection;

std::ostream &operator<<(std::ostream &, const stroke_collection &);

struct symbol_annotation {
	 stroke_collection strokes;
    std::string name;

	symbol_annotation() { }
	 symbol_annotation(const stroke_collection &strokes_, const std::string &name_)
	 	: strokes(strokes_), name(name_) { }
};

std::ostream &operator<<(std::ostream &, const symbol_annotation &);


struct portable_tree {
	SemanticId sid;
	std::list<std::string> content;
	std::list<portable_tree> children;
	size_t nterminals;

	portable_tree() : sid(InvalidSemanticId), nterminals(0) { }


	size_t height() const;
};

int read_portable_tree(std::istream &is, portable_tree &spec);
int write_portable_tree(std::ostream &is, const portable_tree &spec);


struct link_annotation {
	stroke_collection from;
	std::string rel;
	stroke_collection to;
	int type;

	inline bool operator==(const link_annotation &rhs) const
		{ return rel == rhs.rel && type == rhs.type && from == rhs.from && to == rhs.to; }

	link_annotation() { }
	link_annotation(const stroke_collection &from_, const std::string &rel_, const stroke_collection &to_, int type_)
		: from(from_), rel(rel_), to(to_), type(type_) { }
};


struct link_type {
	enum {
		box_box,
		box_head,
		tail_box,
		tail_head,
		default_type
	};
};


struct AnnotatedStrokeGroup : public RawStrokeGroup {
	std::string source;
	std::string writer;
	std::list<symbol_annotation> symbol_annotations;
	std::list<link_annotation> link_annotations;
	portable_tree spec;

	std::vector<size_t> erased;
    
   AnnotatedStrokeGroup() { }

	void add_symbol_annotation(const stroke_collection &strokes,
	                           const std::string &terminal)
	 	{ symbol_annotations.push_back(symbol_annotation(strokes, terminal)); }

	void add_link_annotation(const stroke_collection &from_strokes,
	                         const std::string rel,
									 const stroke_collection &to_strokes,
									 int type = link_type::box_box)
	{
		link_annotation annot(from_strokes, rel, to_strokes, type);
		std::list<link_annotation>::const_iterator i = std::find(link_annotations.begin(), link_annotations.end(), annot);
		if (i == link_annotations.end()) {
	 		link_annotations.push_back(annot);
		}
	}

public:
	typedef std::list<symbol_annotation>::iterator symbol_iterator;
	typedef std::list<symbol_annotation>::const_iterator const_symbol_iterator;

	symbol_iterator symbols_begin() { return symbol_annotations.begin(); }
	symbol_iterator symbols_end() { return symbol_annotations.end(); }
	const_symbol_iterator symbols_begin() const { return symbol_annotations.begin(); }
	const_symbol_iterator symbols_end() const { return symbol_annotations.end(); }

public:
	typedef std::list<link_annotation>::iterator link_iterator;
	typedef std::list<link_annotation>::const_iterator const_link_iterator;

	link_iterator links_begin() { return link_annotations.begin(); }
	link_iterator links_end() { return link_annotations.end(); }
	const_link_iterator links_begin() const { return link_annotations.begin(); }
	const_link_iterator links_end() const { return link_annotations.end(); }
};

size_t num_strokes(const AnnotatedStrokeGroup &g);

int import_annotated_ink(std::istream &is, AnnotatedStrokeGroup &strokes, const grammar *G = 0);
int export_annotated_ink(std::ostream &os, const AnnotatedStrokeGroup &strokes);


}


#endif

