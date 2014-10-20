#ifndef SCG_MATHRECOTYPES_H_
#define SCG_MATHRECOTYPES_H_

#include "dlldecl.h"

#include "grammar-values.h"
#include "grammar-fwd.h"
#include "stroke.h"
#include "group.h"
#include "symbols.h"

#include <ostream>
#include <istream>
#include <stddef.h>


#ifdef WIN32
#define USING_TABLETPC
#define NOMINMAX
#define WIN32_EXTRA_LEAN
#include <msinkaut.h>
#endif



namespace scg
{

class SurrogateTree;


struct DLLDECL ExpressionBox
// Represents the bounding box of a symbol or expression.
{
	virtual int left() const = 0;
	virtual int top() const = 0;
	virtual int right() const = 0;
	virtual int bottom() const = 0;
};


struct DLLDECL MathSymbol
{
	// fill in these first 2 fields before passing the structure
	// into the recognizer
	unicode_char symbol;
	Rect<long> bounds;

	inline bool operator==(const MathSymbol &rhs) const;
	inline bool operator!=(const MathSymbol &rhs) const;

	MathSymbol();
	MathSymbol(unicode_char symbol_, const Rect<long> &bounds_);
};


struct DLLDECL ExpressionTree;
struct DLLDECL ExpressionIterator {
	virtual ~ExpressionIterator() { }
	virtual void release() = 0;
	virtual const ExpressionTree *next() = 0;
	virtual const ExpressionTree *nth(size_t i) = 0;
	virtual size_t count() const = 0;
};


struct DLLDECL ExpressionTree
// Represents a node in an expression tree.
{
	virtual ~ExpressionTree() { }

	void release() const { delete this; }
	// Deallocate this object.

	virtual SemanticId type() const = 0;

	virtual const char *long_str() const = 0;
	virtual const char *latex_str() const = 0;

	virtual const char *str() const = 0;
	// Returns a string representation of the tree rooted at this node.
	// This function uses textual symbol names for terminals, eg. "sigma".

	const char *getstr(int id) const;

	//std::string portable_str() const;

#ifndef NO_WSTRING
	virtual const wchar_t *wstr() const = 0;
#endif
	// Returns a Unicode-string representation of the tree rooted at this node.
	// This function uses Unicode character codes for terminals.

	virtual double score() const = 0;
	// Returns the recognition score of the tree rooted at this node.

	virtual const ExpressionBox *box() const = 0;
	// Returns the bounding box of the tree rooted at this node.

	virtual size_t nchildren() const = 0;
	// Returns the number of children this node has.

	virtual const ExpressionTree *child(size_t i) const = 0;
	// Returns the i'th (from 0) child of this node, or null if the index
	// is invalid.

	virtual bool HasLongForm() const = 0;
	virtual ExpressionIterator *CreateLongFormIterator() const = 0;

	virtual int lock() const = 0;
	virtual int unlock() const = 0;
	virtual bool is_locked() const = 0;

	int serialize(char *buf, size_t *n) const;
};

DLLDECL ExpressionTree *CreateExpressionTree(const char *buf, size_t n);
DLLDECL ExpressionTree *ConvertSurrogateTreeToExpressionTree(const SurrogateTree *);

bool are_trees_equivalent(const ExpressionTree *lhs, const ExpressionTree *rhs);


extern DLLDECL const int DEFAULT_FLAGS;


struct DLLDECL MathRecognizer {
	virtual ~MathRecognizer() { }
	virtual void release() = 0;

	virtual int Clear() = 0;

	virtual void Translate(long x, long y) = 0;

	virtual int Save(std::ostream &os) = 0;

	virtual const ExpressionTree *GetTopExpression() = 0;

	virtual ExpressionIterator *CreateIteratorForExpression(const ExpressionTree *expr, bool ownexpr, bool wrap_mathml = true) = 0;

	virtual ExpressionIterator *CreateDefaultIterator(bool wrap_mathml = true) = 0;
	virtual ExpressionIterator *CreateSubtreeIterator(const ExpressionTree *base, bool wrap_mathml = true) = 0;    
	virtual ExpressionIterator *CreateIteratorForStrokes(const RawStroke **strokes, size_t nstrokes, bool wrap_mathml = true) = 0;
	virtual ExpressionIterator *CreateIteratorForStrokesByIndex(const size_t *strokes, size_t nstrokes, bool wrap_mathml = true) = 0;
	
	virtual ExpressionIterator *CreateDefaultSemanticIterator(SemanticId sid, bool wrap_mathml = true) = 0;
	virtual ExpressionIterator *CreateSemanticIteratorForStrokes(SemanticId sid, const RawStroke **strokes, size_t nstrokes, bool wrap_mathml = true) = 0;
	virtual ExpressionIterator *CreateSemanticIteratorForStrokesByIndex(SemanticId sid, const size_t *strokes, size_t nstrokes, bool wrap_mathml = true) = 0;

    
    virtual int RemoveStrokesByIndex(const size_t *strokes, size_t nstrokes) = 0;
    // Remove the strokes at the specified indices from the recognizer object.  The
    // strokes are considered to be indexed the same way they were added.  Note that removing
    // strokes will change the indices of the strokes that follow.
    
	virtual int AddKnownSymbols(const MathSymbol *symbols, size_t n) = 0;
	virtual int RemoveKnownSymbols(const MathSymbol *symbols, size_t n) = 0;

    virtual int AddStrokes(const RawStroke *strokes, size_t nstrokes) = 0;    
    virtual int AddStrokes(const RawStrokeGroup &group) = 0;
	virtual int RemoveStroke(const RawStroke *stk) = 0;

	virtual size_t GetStrokes(const RawStroke **strokes, size_t n) = 0;


#ifdef USING_TABLETPC
	virtual ExpressionIterator *CreateIteratorForStrokesById(const long *strokes, size_t nstrokes, bool wrap_mathml = true) = 0;
    virtual int RemoveStrokesById(const long *strokes, size_t nstrokes) = 0;
    virtual int AddStrokes(IInkStrokeDisp **strokes, size_t nstrokes) = 0;
    virtual int AddStrokes(IInkDisp *ink) = 0;
	
	virtual int GetInk(IInkDisp **ink) = 0;
	
	virtual void ReleaseInk() = 0;
	
#endif
};


}


#endif
