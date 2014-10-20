#include "annotate.h"
#include "ink-io.h"
#include "error.h"
#include "memory.h"
#include "stream-defs.h"

#include <functional>
#include <sstream>
#include <cctype>

namespace scg
{


std::istream &
operator>>(std::istream &is, stroke_collection &strokes)
{
	char c;
	is >> c;
	if (c != '<') {
		throw E_INVALID;
	}

	for (;;) {
		unsigned n;
		is >> n;
		is >> c;
		if (!is) {
			throw E_INVALID;
		}

		strokes.push_back(n);

		if (c == '>') {
			break;
		}
		// otherwise c should be a comma
	}

	return is;
}

std::istream &
operator>>(std::istream &is, symbol_annotation &sa) {
	is >> sa.strokes >> sa.name;
	if (sa.name == "pm") {
		sa.name = "plusorminus";
	}
	/*if (sa.name == "dot") {
		sa.name = ".";
	}*/
	return is;
}

std::istream &
operator>>(std::istream &is, link_annotation &la)
	{ return is >> la.from >> la.rel >> la.to; }


size_t
portable_tree::height() const
{
	if (children.empty()) {
		return 1;
	}

	size_t H = 0;
	for (std::list<portable_tree>::const_iterator i = children.begin(); i != children.end(); ++i) {
		H = std::max(H, i->height());
	}
	return H + 1;
}

int
write_portable_tree(std::ostream &os, const portable_tree &tree) {
	os << '(';
	os << sid_to_string(tree.sid);
	if (tree.sid == TERMINAL_EXPR) {
		for (std::list<std::string>::const_iterator i = tree.content.begin(); i != tree.content.end(); ++i) {
			os << ((i == tree.content.begin()) ? ' ' : '/');
			const std::string &s = *i;
			for (size_t j = 0; j < s.size(); ++j) {
				if (s[j] == ')') {
					os << "\\)";
				}
				else if (s[j] == '(') {
					os << "\\(";
				}
				else {
					os << s[j];
				}
			}
		}
	}
	for (std::list<portable_tree>::const_iterator i = tree.children.begin(); i != tree.children.end(); ++i) {
		os << ' ';
		write_portable_tree(os, *i);
	}
	os << ')';
	return 0;
}

static int
read_portable_parse_tree(std::istream &is, portable_tree &spec, const std::string &prev) {
	int c;

	while (std::isspace(is.peek())) {
		is.get();
		CHECK_ISTREAM_BASIC(is);
	}

	c = is.get();
	CHECK_ISTREAM_BASIC(is);
	if (c != '(') {
		THROW_ERROR(E_INVALID, "portable parse tree cannot start with \'" << c << '\'');
	}

	while (std::isspace(is.peek())) {
		is.get();
		CHECK_ISTREAM_BASIC(is);
	}

	std::string name;

	for (;;) {
		int c = is.peek();
		CHECK_ISTREAM_BASIC(is);
		if (std::isspace(c) || c == ')') {
			break;
		}
		name.append(1, c);
		is.get();
	}

	if (name == "TERM") {
		spec.sid = TERMINAL_EXPR;
		while (std::isspace(is.peek())) {
			is.get();
			CHECK_ISTREAM_BASIC(is);
		}

		spec.nterminals = 1;
		std::string curr;
		for (;;) {
			c = is.get();
			CHECK_ISTREAM_BASIC(is);

			switch (c) {
			case ' ':  case '\t':
			case '\r': case '\n':
				break;
			case ')':
				spec.content.push_back(curr);
				goto finished;
			case '/':
				++spec.nterminals;
				spec.content.push_back(curr);
				curr.clear();
				break;
			default:
				curr.append(1, c);
				break;
			}
		}
finished:
		if (!spec.content.empty() && spec.content.front() == "...") {
			ERR(E_NOTFOUND, "terminal dots not yet supported internally");
			return E_NOTFOUND;
		}
		if (spec.content.size() > 1 && prev != "NUM") {
			return E_INVALID;
		}
	}
	else {
		if (name == "NAME") {
			name = "VAR";
		}
		else if (name == "MUL") {
			name = "MULT";
		}
		else if (name == "INT") {
			name = "INTEGRAL";
		}
		else if (name == "LIM") {
			name = "LIMIT";
		}
		else if (name == "ADD") {
			spec.children.push_back(portable_tree());
			spec.children.back().sid = TERMINAL_EXPR;
			spec.children.back().content.push_back("plus");
		}
		else if (name == "SUB") {
			spec.children.push_back(portable_tree());
			spec.children.back().sid = TERMINAL_EXPR;
			spec.children.back().content.push_back("horzline");
			name = "ADD";
		}
		else if (name == "PM") {
			spec.children.push_back(portable_tree());
			spec.children.back().sid = TERMINAL_EXPR;
			spec.children.back().content.push_back("plusorminus");
			name = "ADD";
		}

		spec.sid = string_to_sid(name);
		if (spec.sid == InvalidSemanticId && !name.empty()) {
			ERR(E_NOTFOUND, "unknown semantics " << name << " found in annotation tree");
			std::cerr << "unknown semantics " << name << " found in annotation tree\n";
			return E_NOTFOUND;
		}
		for (;;) {
			while (std::isspace(is.peek())) {
				is.get();
				CHECK_ISTREAM_BASIC(is);
			}

			c = is.peek();
			if (c == ')') {
				is.get();
				break;
			}
			else if (c == '(') {
				spec.children.push_back(portable_tree());
				int e = read_portable_parse_tree(is, spec.children.back(), name);
				if (FAILURE(e)) {
					return e;
				}
				spec.nterminals += spec.children.back().nterminals;
			}
			else {
				std::string term;
				while (!std::isspace(is.peek())) {
					term += is.get();
					CHECK_ISTREAM_BASIC(is);
				}
				portable_tree pt;
				pt.sid = TERMINAL_EXPR;
				std::string::size_type i = 0;
				while (i != std::string::npos) {
					i = term.find_first_of('/', i);
				}
				//THROW_ERROR(E_INVALID, "unexpected character while reading portable parse tree: " << c);
			}
		}
	}

	if (spec.sid == ADD_EXPR) {
		spec.children.push_back(spec.children.front());
		spec.children.pop_front();
	}
	else if (spec.sid == MULT_EXPR) {
		if (spec.children.back().sid == NUM_EXPR) {
			ERR(E_NOTFOUND, "do not support RHS numbers");
			return E_NOTFOUND;
		}
		if (spec.children.front().sid == NUM_EXPR && prev == "MULT") {
			ERR(E_NOTFOUND, "do not support RHS numbers");
			return E_NOTFOUND;
		}
	}
	else if (spec.sid == VAR_EXPR) {
		portable_tree tmp = spec.children.front();
		if (tmp.content.size() > 1 && prev != "FN") {
			ERR(E_NOTFOUND, "do not support names outside of functions");
			return E_NOTFOUND;
		}
	}
	else if (spec.sid == REL_EXPR) {
		portable_tree tmp = spec.children.front();
		if (tmp.sid == TERMINAL_EXPR) {
			spec.children.pop_front();
			spec.children.push_back(tmp);
		}
	}
	else if (spec.sid == PAREN_EXPR) {
		portable_tree tmp = spec.children.front();
		if (tmp.content.size() == 1 && tmp.content.front() == "ROUND") {
			spec.children.front().content.front() = "(";
			spec.children.push_back(portable_tree());
			spec.children.back().sid = TERMINAL_EXPR;
			spec.children.back().content.push_back(")");
		}
		else {
			ERR(E_NOTFOUND, "do not support non-round parentheses");
			return E_NOTFOUND;
		}
		/*else if (tmp.content.size() == 1 && tmp.content.front() == "SQUARE") {
			spec.children.front().content.front() = "[";
			spec.children.push_back(portable_tree());
			spec.children.back().sid = TERMINAL_EXPR;
			spec.children.back().content.push_back("]");
		}
		else if (tmp.content.size() == 1 && tmp.content.front() == "CURLY") {
			spec.children.front().content.front() = "{";
			spec.children.push_back(portable_tree());
			spec.children.back().sid = TERMINAL_EXPR;
			spec.children.back().content.push_back("}");
		}*/
	}

	/*if (spec.sid == TERMINAL_EXPR) {
		if (spec.content.size() == 3) {
			std::list<std::string>::const_iterator i = std::find_if(spec.content.begin(), spec.content.end(),
			                                                        std::bind1st(std::not_equal_to<std::string>(), "."));
			if (i == spec.content.end()) {
				spec.content.clear();
				spec.content.push_back("hdots");
			}
		}
	}*/

	return 0;
}

int
read_portable_tree(std::istream &is, portable_tree &pt) {
	return read_portable_parse_tree(is, pt, "");
}


static void
validate_strokes(const AnnotatedStrokeGroup &G, stroke_collection &strokes)
{
	for (stroke_collection::iterator i = strokes.begin(); i != strokes.end(); ) {
		if (std::find(G.erased.begin(), G.erased.end(), *i) != G.erased.end()) {
			//std::cerr << "erasing stroke " << *i << " from annotation\n";
			i = strokes.erase(i);
		}
		else {
			std::vector<size_t>::const_iterator j = std::upper_bound(G.erased.begin(), G.erased.end(), *i);
			size_t offset = j - G.erased.begin();
			*i -= offset;
			++i;
		}
	}
}


int
import_annotated_ink(std::istream &is, AnnotatedStrokeGroup &strokes, const grammar *G) {
	 while (std::isspace(is.peek())) is.get();
	 while (is.peek() == '#' || is.peek() == '%') {
	 	std::string meta;
		is >> meta;
		if (meta == "#writer:") {
			is >> strokes.writer;
		}
		else if (meta == "#source:") {
			is >> strokes.source;
		}
		else if (meta[0] != '%') {
			THROW_ERROR(E_INVALID, "unknown ink meta-tag " << meta);
		}
		while (std::isspace(is.peek())) is.get();
	 }

	 try {
    	is >> static_cast<RawStrokeGroup &>(strokes);
	}
	catch (int e) {
        strokes.clear();
        return e;
    }


	RawStroke *newstrokes;
	size_t n_newstrokes = 0;
	for (size_t i = 0; i < strokes.nstrokes; ++i) {
		const RawStroke &s = strokes.strokes[i];
		if (s.npoints > 1) {
			++n_newstrokes;
		}
		else {
			strokes.erased.push_back(i);
		}
	}

	newstrokes = new RawStroke[n_newstrokes];
	n_newstrokes = 0;
	for (size_t i = 0; i < strokes.nstrokes; ++i) {
		const RawStroke &s = strokes.strokes[i];
		if (s.npoints > 1) {
			RawStroke S = s.copy();
			newstrokes[n_newstrokes] = S;
			++n_newstrokes;
		}
	}

	strokes.set_strokes(newstrokes, n_newstrokes);
    
	 std::string key;
	 is >> key;

    if (!is) {
	 	// no annotations
        return 0;
    }
    
	// the string ANNOTATIONS indicates new-style symbol and link annotations
	if (key == "ANNOTATIONS") {
		for (;;) {
			is >> key;
			if (!is) {
				break;
			}

			if (key == "SYMBOL") {
				symbol_annotation sa;
				is >> sa;
				validate_strokes(strokes, sa.strokes);
				strokes.symbol_annotations.push_back(sa);
			}
			else if (key == "SYMBOLMAP") {
				stroke_collection s;
				is >> s;
				validate_strokes(strokes, s);
				size_t d;
				is >> d;
			}
			else if (key == "TREE") {
				int e = read_portable_parse_tree(is, strokes.spec, "");
				if (FAILURE(e)) {
					return e;
				}
			}
			else if (key == "LINK") {
				link_annotation la;
				is >> la;
				validate_strokes(strokes, la.from);
				validate_strokes(strokes, la.to);
				la.type = 0;
				for (std::list<symbol_annotation>::const_iterator i = strokes.symbol_annotations.begin(); i != strokes.symbol_annotations.end(); ++i) {
					if (i->strokes == la.from) {
						la.type += 2;
					}
					if (i->strokes == la.to) {
						la.type += 1;
					}
				}
				strokes.link_annotations.push_back(la);
			}
			else if (key == "LINKBB") {
				link_annotation la;
				is >> la;
				validate_strokes(strokes, la.from);
				validate_strokes(strokes, la.to);
				la.type = link_type::box_box;
				strokes.link_annotations.push_back(la);
			}
			else if (key == "LINKTH") {
				link_annotation la;
				is >> la;
				validate_strokes(strokes, la.from);
				validate_strokes(strokes, la.to);
				la.type = link_type::tail_head;
				strokes.link_annotations.push_back(la);
			}
			else if (key == "LINKTB") {
				link_annotation la;
				is >> la;
				validate_strokes(strokes, la.from);
				validate_strokes(strokes, la.to);
				la.type = link_type::tail_box;
				strokes.link_annotations.push_back(la);
			}
			else if (key == "LINKBH") {
				link_annotation la;
				is >> la;
				validate_strokes(strokes, la.from);
				validate_strokes(strokes, la.to);
				la.type = link_type::box_head;
				strokes.link_annotations.push_back(la);
			}
			else {
				return E_INVALID;
			}
		}
	}
	else {
		// else we need to convert old-style indexed annotations to
		// new-style

    	size_t nannot;

		std::stringstream ss;
		ss << key;
		ss >> nannot;
		if (ss.fail()) {
			return E_INVALID;
		}
    
    	while (nannot--) {
     	  if (!is) {
            return E_IO;
        }
        
        size_t start;
		  size_t nstrokes;
		  std::string name;
        is >> start;
        is >> name;
        is >> nstrokes;
		  if (!is || is.fail()) {
		  	return E_INVALID;
		  }

		  std::list<size_t> S;
		  for (size_t i = start; i < start + nstrokes; ++i) {
		  	S.push_back(i);
		  }

			validate_strokes(strokes, S);
		  strokes.add_symbol_annotation(S, name);
    }

		is >> nannot;
		if (is.eof()) {
			return 0;
		}
		if (!is) {
			return E_IO;
		}

		while (nannot--) {
			if (!is) {
				return E_IO;
			}

			size_t g1, g2;
			std::string rel;
			is >> g1 >> g2 >> rel;
			if (!is) {
				return E_IO;
			}

			link_annotation la;
			if (rel == "L") {
				la.rel = "R";
			}
			else if (rel == "BR") {
				la.rel = "BR";
			}
			else if (rel == "BL") {
				la.rel = "AR";
			}
			else if (rel == "A") {
				la.rel = "B";
			}
			else if (rel == "C") {
				la.rel = "C";
			}
			else {
				continue;
			}

			std::list<symbol_annotation>::const_iterator from = strokes.symbol_annotations.begin();
			std::list<symbol_annotation>::const_iterator to = strokes.symbol_annotations.begin();
			while (g1--) { ++from; if (from == strokes.symbol_annotations.end()) return E_INVALID; }
			while (g2--) { ++to; if (to == strokes.symbol_annotations.end()) return E_INVALID; }
			la.type = link_type::tail_head;
			la.from = from->strokes;
			la.to = to->strokes;
			strokes.link_annotations.push_back(la);
		}
	}

  return 0;
}


std::ostream &
operator<<(std::ostream &os, const stroke_collection &strokes)
{
	os << '<';
	for (stroke_collection::const_iterator s = strokes.begin(); s != strokes.end(); ++s) {
		os << *s;
		stroke_collection::const_iterator next = s;
		++next;
		if (next != strokes.end()) {
			os << ", ";
		}
	}

	os << '>';
	return os;
}

std::ostream &
operator<<(std::ostream &os, const symbol_annotation &sa)
	{ return os << "SYMBOL " << sa.strokes << ' ' << sa.name; }

static std::ostream &
operator<<(std::ostream &os, const link_annotation &la)
{
	switch (la.type) {
	case link_type::box_box:
		os << "LINKBB ";
		break;
	case link_type::tail_box:
		os << "LINKTB ";
		break;
	case link_type::box_head:
		os << "LINKBH ";
		break;
	case link_type::tail_head:
		os << "LINKTH ";
		break;
	case link_type::default_type:
		os << "LINK ";
		break;
	default:
		THROW_ERROR(E_INVALID, "invalid link type " << la.type);
	}
	return os << la.from << ' ' << la.rel << ' ' << la.to;
}


int
export_annotated_ink(std::ostream &os, const AnnotatedStrokeGroup &strokes) {
	if (!strokes.source.empty()) {
		os << "#source: " << strokes.source << std::endl;
	}
	if (!strokes.writer.empty()) {
		os << "#writer: " << strokes.writer << std::endl;
	}
    os << static_cast<const RawStrokeGroup &>(strokes) << std::endl;
    
	 os << "ANNOTATIONS\n";

	 for (std::list<symbol_annotation>::const_iterator psa = strokes.symbol_annotations.begin(); psa != strokes.symbol_annotations.end(); ++psa) {
	 	os << *psa << std::endl;
	 }

	 for (std::list<link_annotation>::const_iterator pla = strokes.link_annotations.begin(); pla != strokes.link_annotations.end(); ++pla) {
	 	os << *pla << std::endl;
	 }

    return 0;
}


size_t num_strokes(const AnnotatedStrokeGroup &g) { return g.nstrokes; }


}

