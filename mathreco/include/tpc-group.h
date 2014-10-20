/*
 * tpc-group.h
 *
 * This file defines the interface to a stroke group object for TabletPC
 * strokes.  It is primarily used to convert between TabletPC stroke format
 * and the format used internally by libmathrec.
 */
 
#ifndef TPC_GROUP_H_
#define TPC_GROUP_H_

#ifdef WIN32


#include "group.h"
#include "dlldecl.h"

#define NOMINMAX
#include <msinkaut.h>
#include <windows.h>

namespace scg
{


struct TPC_InkWindow;


struct TPC_StrokeGroup
{
    IInkStrokeDisp **strokes;
    
    size_t nstrokes;

    enum Ownership {
        WEAK,
        STRONG
    } own;
    
    TPC_StrokeGroup() : own(WEAK), strokes(0), nstrokes(0) {}
	 TPC_StrokeGroup(IInkDisp *ink);
    TPC_StrokeGroup(IInkStrokeDisp **s, size_t n, Ownership own_ = STRONG) : strokes(s), nstrokes(n), own(own_) {}
    TPC_StrokeGroup(TPC_StrokeGroup &g)
    {
        *this = g;
    }
    
    ~TPC_StrokeGroup() { clear(); }
    
    void clear()
    {
        if (own == STRONG) {
            IInkStrokeDisp **end = strokes + nstrokes;
            for (IInkStrokeDisp **stroke = strokes; stroke != end; stroke++) {
                (*stroke)->Release();
            }
        
            delete[] strokes;
        }
        strokes = 0;
        nstrokes = 0;
    }
    
    void set_strokes(IInkStrokeDisp **s, size_t n, Ownership own_ = STRONG)
    {
        clear();
        strokes = s;
        nstrokes = n;
        own = own_;
    }
    
	 int add_to_ink(IInkDisp *ink);

    TPC_StrokeGroup &operator=(TPC_StrokeGroup &g)
    {
        clear();
        
        nstrokes = g.nstrokes;
        strokes = g.strokes;
        own = g.own;
        
        g.strokes = 0;
        g.nstrokes = 0;

        return *this;
    }
};

extern TPC_StrokeGroup TPC_EmptyStrokeGroup;


IInkDisp *make_empty_ink_object();


size_t num_strokes(const TPC_StrokeGroup &g);


int convert(RawStrokeGroup &dest, TPC_StrokeGroup &src, size_t first, size_t n);
int convert(RawStrokeGroup &dest, TPC_StrokeGroup &src);

int convert(RawStrokeGroup &dest, IInkDisp *ink);

int convert(TPC_StrokeGroup &dest, const RawStroke *first, const RawStroke *last, IInkDisp *ink);
int convert(TPC_StrokeGroup &dest, const RawStroke *src, size_t n, TPC_InkWindow &ink);
int convert(TPC_StrokeGroup &dest, const RawStrokeGroup &src, TPC_InkWindow &ink);

DLLDECL int convert(const RawStroke &stroke, IInkDisp *ink, IInkStrokeDisp **ms_stroke);


}


#endif

#endif
