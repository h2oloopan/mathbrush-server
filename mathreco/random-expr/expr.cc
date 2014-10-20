#include "expr.h"
#include "group.h"
#include "ink-io.h"
#include "vector-io.h"

#include <algorithm>
#include <string>


namespace scg
{


static std::string centered_names[] = { "a", "c", "e", "m", "n", "o", "r", "s", "u", "v", "w", "x", "alpha", "mu", "pi", "sigma" };
static size_t num_centered_names = sizeof(centered_names) / sizeof(std::string);

static std::string ascender_names[] = { "zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "b", "d", "f", "h", "k", "l", "beta", "gamma", "lambda", "phi", "psi", "theta", "xi", "zeta", "Delta", "Gamma", "Omega" };
static size_t num_ascender_names = sizeof(ascender_names) / sizeof(std::string);

static std::string descender_names[] = { "g", "p", "q", "y", "z", "rho" };
static size_t num_descender_names = sizeof(descender_names) / sizeof(std::string);

static std::string threequarter_names[] = { "i", "t" };
static size_t num_threequarter_names = sizeof(threequarter_names) / sizeof(std::string);

static std::string threequarterdescender_names[] = { "j" };
static size_t num_threequarterdescender_names = sizeof(threequarterdescender_names) / sizeof(std::string);

static std::string fence_names[] = { "lbracket", "rbracket", "lbrace", "rbrace", "lparen", "rparen", "exclaim" };
static size_t num_fence_names = sizeof(fence_names) / sizeof(std::string);

static std::string centeredop_names[] = { "plus", "eq", "neq", "leq", "geq", "lt", "gt", "infin", "equiv", "plusorminus" };
static size_t num_centeredop_names = sizeof(centeredop_names) / sizeof(std::string);

static std::string largeop_names[] = { "Pi", "Sigma", "Integral" };
static size_t num_largeop_names = sizeof(largeop_names) / sizeof(std::string);

static std::string container_names[] = { "sqrt" };
static size_t num_container_names = sizeof(container_names) / sizeof(std::string);

static std::string short_names[] = { "horzline", "arrow" };
static size_t num_short_names = sizeof(short_names) / sizeof(std::string);


static TerminalType
LookupTerminalType(const std::string &name)
{
	if (std::find(centered_names, centered_names + num_centered_names, name) != centered_names + num_centered_names) {
		return Centered;
	}
	else if (std::find(ascender_names, ascender_names + num_ascender_names, name) != ascender_names + num_ascender_names) {
		return Ascender;
	}
	else if (std::find(descender_names, descender_names + num_descender_names, name) != descender_names + num_descender_names) {
		return Descender;
	}
	else if (std::find(threequarter_names, threequarter_names + num_threequarter_names, name) != threequarter_names + num_threequarter_names) {
		return ThreeQuarter;
	}
	else if (std::find(threequarterdescender_names, threequarterdescender_names + num_threequarterdescender_names, name) != threequarterdescender_names + num_threequarterdescender_names) {
		return ThreeQuarterDescender;
	}
	else if (std::find(fence_names, fence_names + num_fence_names, name) != fence_names + num_fence_names) {
		return Fence;
	}
	else if (std::find(centeredop_names, centeredop_names + num_centeredop_names, name) != centeredop_names + num_centeredop_names) {
		return CenteredOperator;
	}
	else if (std::find(largeop_names, largeop_names + num_largeop_names, name) != largeop_names + num_largeop_names) {
		return LargeOperator;
	}
	else if (std::find(container_names, container_names + num_container_names, name) != container_names + num_container_names) {
		return Container;
	}
	else if (std::find(short_names, short_names + num_short_names, name) != short_names + num_short_names) {
		return Short;
	}

	return static_cast<TerminalType>(-1);
}

std::istream &
operator>>(std::istream &is, scg::ExpressionLink &link)
{
	is >> link.box1;
	is >> link.box2;

	std::string type;
	is >> type;

	if (type == "A") {
		link.type = scg::BelowWide;
	}
	else if (type == "AN") {
		link.type = scg::BelowNarrow;
	}
	else if (type == "L") {
		link.type = scg::Right;
	}
	else if (type == "BR") {
		link.type = scg::BelowRight;
	}
	else if (type == "BL") {
		link.type = scg::AboveRight;
	}
	else if (type == "C") {
		link.type = scg::Contains;
	}

	return is;
}


std::istream &
operator>>(std::istream &is, scg::Expression &expr)
{
	scg::RawStrokeGroup strokes;
	is >> strokes;

	expr.boxes.clear();
	expr.links.clear();

	unsigned nboxes;
	is >> nboxes;

	if (is.eof()) {
	    	throw E_INVALID;
	}

	scg::SymbolBox *boxes = new scg::SymbolBox[nboxes];
	if (!boxes) {
		throw E_OUTOFMEM;
	}

	for (scg::SymbolBox *box = boxes; box < boxes + nboxes; box++) {
		unsigned start;
		unsigned size;

		is >> start >> box->name >> size;

		box->strokes = scg::copy(strokes.strokes + start, strokes.strokes + start + size);
		box->bbox = scg::bbox(box->strokes);
		// TODO: is it possible for types to be determined by analyzing the
		// symbol database?  If so, we would possibly need different types
		// for different models of the same symbol; in this case storing the
		// symbol name in the data file is insufficient.

		box->type = scg::LookupTerminalType(box->name);
	}

	expr.boxes.set_array(boxes, nboxes);
	is >> expr.links;

	return is;
}


std::ostream &
operator<<(std::ostream &os, scg::Expression &expr)
{
	os << "Linked boxes: \n";

	bool *linked = new bool[expr.boxes.size()];
	std::fill(linked, linked + expr.boxes.size(), false);

	for (std::vector<scg::ExpressionLink>::const_iterator i = expr.links.begin(); i != expr.links.end(); ++i) {
		const scg::SymbolBox &box1 = expr.boxes[i->box1];
		const scg::SymbolBox &box2 = expr.boxes[i->box2];

		linked[i->box1] = linked[i->box2] = true;

		os << box1.name << " (" << box1.bbox.left << "," << box1.bbox.top << ")->(" << box1.bbox.right << "," << box1.bbox.bottom << ") , " << box2.name << " (" << box2.bbox.left << "," << box2.bbox.top << ")->(" << box2.bbox.right << "," << box2.bbox.bottom << ") : ";
		if (i->type == scg::BelowWide) {
			os << "A\n";
		}
		else if (i->type == scg::BelowNarrow) {
			os << "AN\n";
		}
		else if (i->type == scg::Right) {
			os << "L\n";
		}
		else if (i->type == scg::BelowRight) {
			os << "BR\n";
		}
		else if (i->type == scg::AboveRight) {
			os << "BL\n";
		}
		else if (i->type == scg::Contains) {
			os << "C\n";
		}
	}

	os << "Unlinked boxes: \n";
	for (size_t i = 0; i < expr.boxes.size(); i++) {
		if (!linked[i]) {
			const scg::Rect<long> &box = expr.boxes[i].bbox;
			os << "(" << box.left << "," << box.top << ")->(" << box.right << "," << box.bottom << ")\n";
		}
	}

	delete[] linked;

	return os;
}

}

std::istream &
operator>>(std::istream &is, scg::Expression &expr)
	{ return scg::operator>>(is, expr); }
