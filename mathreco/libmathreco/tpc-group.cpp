#ifdef WIN32

#include "tpc-group.h"
#include "tabletpc.h"
#include "error.h"
#include "memory.h"


namespace scg
{


TPC_StrokeGroup TPC_EmptyStrokeGroup(0, 0);

TPC_StrokeGroup::TPC_StrokeGroup(IInkDisp *ink)
{
	IInkStrokes *strokes;
	HRESULT hr = ink->get_Strokes(&strokes);
	if (FAILED(hr)) {
		THROW_ERROR(MAKE_API_ERROR(hr), "cannot retrieve strokes collection");
	}

	long n;
	hr = strokes->get_Count(&n);
	if (FAILED(hr)) {
		THROW_ERROR(MAKE_API_ERROR(hr), "cannot retrieve stroke count");
	}

	IInkStrokeDisp **S = new IInkStrokeDisp *[n];
	IInkStrokeDisp **s = S;
	for (long i = 0; i < n; ++i) {
		hr = strokes->Item(i, s);
		if (FAILED(hr)) {
			THROW_ERROR(MAKE_API_ERROR(hr), "cannot retrieve stroke at index " << i);
		}
		++s;
	}
	
	strokes->Release();

	set_strokes(S, n);
	own = STRONG;
}


int
TPC_StrokeGroup::add_to_ink(IInkDisp *ink)
{
	if (nstrokes == 0) {
		WARNING(E_NOTFOUND, "TPC_StrokeGroup is empty");
		return 0;
	}

	IInkStrokes *msstrokes;
	VARIANT dummy;
	VariantInit(&dummy);

	HRESULT hr = ink->CreateStrokes(dummy, &msstrokes);
	if (FAILED(hr)) {
		ERR(MAKE_API_ERROR(hr), "cannot create IInkStrokes");
		return 0;
	}
	
	for (IInkStrokeDisp **ps = strokes; ps != strokes + nstrokes; ++ps) {
		hr = msstrokes->Add(*ps);
		if (FAILED(hr)) {
			int e = MAKE_API_ERROR(hr);
			ERR(e, "failed adding tablet stroke to stroke collection");
			return e;
		}
	}
	/*
	IInkRectangle *bounds;
	hr = msstrokes->GetBoundingBox(IBBM_Default, &bounds);
	if (FAILED(hr)) {
		int e = MAKE_API_ERROR(hr);
		ERR(e, "failed to obtain MS stroke bounds");
		msstrokes->Release();
		return e;
	}

	hr = ink->AddStrokesAtRectangle(msstrokes, bounds);
	if (FAILED(hr)) {
		int e = MAKE_API_ERROR(hr);
		ERR(e, "could not add strokes to ink object");
		bounds->Release();
		msstrokes->Release();
		return e;
	}
	
	bounds->Release();*/
	msstrokes->Release();
	return 0;
}


unsigned
num_strokes(const TPC_StrokeGroup &g)
{
    return g.nstrokes;
}


static int
convert(RawStroke &dest, IInkStrokeDisp *src)
{
	HRESULT hr;
	VARIANT varPts;
	SAFEARRAY *APIpoints;
	
	hr = src->GetPoints(ISC_FirstElement, ISC_AllElements, &varPts);
	if (FAILED(hr)) {
	    return MAKE_API_ERROR(hr);
	}
	if (varPts.vt != (VT_ARRAY | VT_I4)) {
		return E_INVALID;
	}

	APIpoints = varPts.parray;
	if (SafeArrayGetDim(APIpoints) != 1) {
		return E_INVALID;
	}
	
	long lo, hi;
	long *ptsarr;

	hr = SafeArrayGetLBound(APIpoints, 1, &lo);
	if (FAILED(hr)) {
		return MAKE_API_ERROR(hr);
	}
	
	hr = SafeArrayGetUBound(APIpoints, 1, &hi);
	if (FAILED(hr)) {
		return MAKE_API_ERROR(hr);
	}
	
	hr = SafeArrayAccessData(APIpoints, (void HUGEP **)&ptsarr);
	if (FAILED(hr)) {
		return MAKE_API_ERROR(hr);
	}
		
	size_t sz = hi - lo + 1;

	long *x = new long[sz / 2];
	if (!x) {
	    SafeArrayUnaccessData(APIpoints);
	    return E_OUTOFMEM;
	}
	
	long *y = new long[sz / 2];
	if (!y) {
	    SafeArrayUnaccessData(APIpoints);
	    return E_OUTOFMEM;
	}

	long *px = x;
	long *py = y;

	for (long *pt = ptsarr; pt != ptsarr + sz; pt += 2) {
		*(px++) = *pt;
		*(py++) = *(pt + 1);
	}
	
	SafeArrayUnaccessData(APIpoints);

    dest.set_points(x, y, sz / 2);
    
    return 0;
}


int
convert(RawStrokeGroup &dest, IInkDisp *ink)
{
    int e;
    
	IInkStrokes *strokes;
	HRESULT hr = ink->get_Strokes(&strokes);
	if (FAILED(hr)) {
		e = MAKE_API_ERROR(hr);
		ERR(e, "cannot retrieve strokes collection");
		return e;
	}

	long n;
	hr = strokes->get_Count(&n);
	if (FAILED(hr)) {
		e = MAKE_API_ERROR(hr);
		ERR(e, "cannot retrieve stroke count");
		return e;
	}

    RawStroke *dest_strokes = new RawStroke[n];
    if (!dest_strokes) {
        return E_OUTOFMEM;
    }

    RawStroke *dest_stroke = dest_strokes;
	for (long i = 0; i < n; ++i) {
		IInkStrokeDisp *s;
		hr = strokes->Item(i, &s);
		if (FAILED(hr)) {
			THROW_ERROR(MAKE_API_ERROR(hr), "cannot retrieve stroke at index " << i);
		}
        e = convert(*(dest_stroke++), s);
        if (FAILURE(e)) {
            return e;
        }
	}
	
	strokes->Release();
    dest.set_strokes(dest_strokes, n);
        
    return 0;
}


int
convert(RawStrokeGroup &dest, TPC_StrokeGroup &src, size_t first, size_t n)
{
    int e;
    
    if (first > num_strokes(src)) {
        return E_INVALID;
    }
    if (first + n > num_strokes(src)) {
        n = num_strokes(src) - first;
    }
    
    if (n == 0) {
        return 0;
    }
    
    RawStroke *dest_strokes = new RawStroke[n];
    if (!dest_strokes) {
        return E_OUTOFMEM;
    }
    
    IInkStrokeDisp **end_of_src_strokes = src.strokes + first + n;    
    RawStroke *dest_stroke = dest_strokes;
    
    for (IInkStrokeDisp **src_stroke = src.strokes + first; src_stroke != end_of_src_strokes; src_stroke++) {
        e = convert(*(dest_stroke++), *src_stroke);
        if (FAILURE(e)) {
            return e;
        }
    }

    dest.set_strokes(dest_strokes, n);
        
    return 0;
}


IInkDisp *
make_empty_ink_object()
{
	IInkDisp *ink;
	HRESULT hr = CoCreateInstance(CLSID_InkDisp, NULL, CLSCTX_INPROC_SERVER, IID_IInkDisp, (void **)&ink);
	if (FAILED(hr)) {
		ERR(MAKE_API_ERROR(hr), "cannot create IInkDisp");
		return 0;
	}
	return ink;
}

int
convert(RawStrokeGroup &dest, TPC_StrokeGroup &src)
{
    return convert(dest, src, 0, num_strokes(src));
}

int
convert(TPC_StrokeGroup &dest, const RawStroke *first, const RawStroke *last, TPC_InkWindow &ink)
{
	IInkDisp *APIink = TPC_get_ink(ink);
    if (!APIink) {
        int e = errval;
        errval = 0;
        return e;
    }

	int e = convert(dest, first, last, APIink);
	APIink->Release();
	return e;
}

int
convert(const RawStroke &stroke, IInkDisp *ink, IInkStrokeDisp **ms_stroke)
{
    HRESULT hr;
 	
    SAFEARRAY *pts = SafeArrayCreateVector(VT_I4, 0, 2 * static_cast<ULONG>(stroke.npoints));
    if (!pts) {
	    return MAKE_API_ERROR(GetLastError());
    }
 	
    long *ptsarr;
    hr = SafeArrayAccessData(pts, reinterpret_cast<void HUGEP **>(&ptsarr));
     if (FAILED(hr)) {
         return MAKE_API_ERROR(hr);
     }
     
     long *x = stroke.x, *y = stroke.y;
    while (x != stroke.x + stroke.npoints) {
	    *(ptsarr++) = *(x++);
	    *(ptsarr++) = *(y++);
    }
 	
    hr = SafeArrayUnaccessData(pts);
     if (FAILED(hr)) {
         return MAKE_API_ERROR(hr);
     }
 	
    VARIANT dummy;
    VARIANT stroke_data;

    VariantInit(&dummy);
    VariantInit(&stroke_data);

    stroke_data.vt = VT_ARRAY | VT_I4;
    stroke_data.parray = pts;

	IInkStrokeDisp *stk;
    hr = ink->CreateStroke(stroke_data, dummy, &stk);
     if (FAILED(hr)) {
         return MAKE_API_ERROR(hr);
     }
     
     VariantClear(&dummy);
     VariantClear(&stroke_data);
     
     if (ms_stroke) {
		*ms_stroke = stk;
     }
     
     return 0;
}


int
convert(TPC_StrokeGroup &dest, const RawStroke *first, const RawStroke *last, IInkDisp *ink)
{
    size_t n = last - first;
    
    IInkStrokeDisp **strokes = new IInkStrokeDisp *[n];
    if (!strokes) {
        return E_OUTOFMEM;
    }

    IInkStrokeDisp **stroke = strokes;
        
    while (first != last) {
	    HRESULT hr;
    	
	    SAFEARRAY *pts = SafeArrayCreateVector(VT_I4, 0, 2 * static_cast<ULONG>(first->npoints));
	    if (!pts) {
		    return MAKE_API_ERROR(GetLastError());
	    }
    	
	    long *ptsarr;
	    hr = SafeArrayAccessData(pts, reinterpret_cast<void HUGEP **>(&ptsarr));
        if (FAILED(hr)) {
            return MAKE_API_ERROR(hr);
        }
        
        long *x = first->x, *y = first->y;
	    while (x != first->x + first->npoints) {
		    *(ptsarr++) = *(x++);
		    *(ptsarr++) = *(y++);
	    }
    	
	    hr = SafeArrayUnaccessData(pts);
        if (FAILED(hr)) {
            return MAKE_API_ERROR(hr);
        }
    	
	    VARIANT dummy;
	    VARIANT stroke_data;

	    VariantInit(&dummy);
	    VariantInit(&stroke_data);

	    stroke_data.vt = VT_ARRAY | VT_I4;
	    stroke_data.parray = pts;

	    hr = ink->CreateStroke(stroke_data, dummy, stroke);
        if (FAILED(hr)) {
            return MAKE_API_ERROR(hr);
        }
        
        VariantClear(&dummy);
        VariantClear(&stroke_data);
    
        ++first;
        ++stroke;
    }
    
    dest.set_strokes(strokes, n);
    
    return 0;
}


int
convert(TPC_StrokeGroup &dest, const RawStrokeGroup &src, TPC_InkWindow &ink)
{
    return convert(dest, src.strokes, src.strokes + src.nstrokes, ink);
}


}

#endif
