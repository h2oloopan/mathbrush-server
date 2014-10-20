#ifndef EXPRESSION_H_
#define EXPRESSION_H_


#include "group.h"
#include "rect.h"

#include <istream>
#include <ostream>
#include <vector>


namespace scg
{



enum TerminalType
{
	Centered,
	Ascender,
	Descender,
	ThreeQuarter,
	ThreeQuarterDescender,
	Fence,
	CenteredOperator,
	LargeOperator,
	Container,
	Short,
	NumSymbolTypes
};

struct SymbolBox
{
	RawStrokeGroup strokes;
	Rect<long> bbox;
	//TerminalType type;
	int type;
	std::string name;
};


enum LinkType
{
	AboveRight,
	Right,
	BelowRight,
	BelowWide,
	BelowNarrow,
	Contains,
	StrokeOrder,
	NumLinkTypes
};

struct ExpressionLink
{
	size_t box1;
	size_t box2;
	LinkType type;
};

struct Expression
{
	OwnedArrayType(SymbolBox) boxes;
	std::vector<ExpressionLink> links;
};


std::istream &
operator>>(std::istream &is, scg::Expression &expr);

std::ostream &
operator<<(std::ostream &os, scg::Expression &expr);


}


#endif
