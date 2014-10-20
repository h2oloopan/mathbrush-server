#ifdef WIN32

#include "tabletpc.h"
#include "error.h"
#include "memory.h"

#include <new>

#include <msinkaut_i.c>


namespace scg
{


int
TPC_create_window(HWND hwnd, TPC_InkWindow &ink)
{
	HRESULT hr;

	ink.collector = 0;
	ink.overlay = 0;
	ink.renderer = 0;

	hr = CoCreateInstance(CLSID_InkCollector, NULL, CLSCTX_INPROC_SERVER, IID_IInkCollector, (void **)&ink.collector);
	if (FAILED(hr)) {
	    return MAKE_API_ERROR(hr);
	}
    
    /*    
	SAFEARRAY *prop_array = SafeArrayCreateVector(VT_BSTR, 0, 1);
	if (!prop_array) {
	    ink.collector->Release();
	    ink.collector = 0;
	    return MAKE_API_ERROR(GetLastError());	
	}
	
	BSTR *prop_strs;
	hr = SafeArrayAccessData(prop_array, reinterpret_cast<void HUGEP **>(&prop_strs));
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    
    prop_strs[0] = SysAllocString(STR_GUID_TIMERTICK);

    SafeArrayUnaccessData(prop_array);
    
	VARIANT properties;
	VariantInit(&properties);
	properties.vt = VT_ARRAY | VT_BSTR;
	properties.parray = prop_array;
	
	hr = ink.collector->put_DesiredPacketDescription(properties);
	VariantClear(&properties);

	if (FAILED(hr)) {
	    ink.collector->Release();
	    ink.collector = 0;
	    return MAKE_API_ERROR(hr);
	}
	*/
	
	hr = ink.collector->put_hWnd(reinterpret_cast<LONG_PTR>(hwnd));
	if (FAILED(hr)) {
	    ink.collector->Release();
	    ink.collector = 0;
	    return MAKE_API_ERROR(hr);
	}

	return 0;
}


int
TPC_release_window(TPC_InkWindow &ink)
{
    if (ink.overlay) {
        ink.overlay->Release();
    }	
    if (ink.renderer) {
        ink.renderer->Release();
    }
    if (ink.collector) {
        ink.collector->Release();
    }
    
    return 0;
}


int
TPC_enable_collector(TPC_InkWindow &ink)
{
    if (!ink.collector) {
        return E_INVALID;
    }
    
    HRESULT hr;
    
    hr = ink.collector->put_Enabled(TRUE);
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    
    return 0;
}


int
TPC_disable_collector(TPC_InkWindow &ink)
{
    if (!ink.collector) {
        return E_INVALID;
    }

    HRESULT hr;
    
    hr = ink.collector->put_Enabled(FALSE);
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    
    return 0;
}


int
TPC_set_ink(TPC_InkWindow &inkw, IInkDisp *ink)
{
	if (!inkw.collector) {
		return E_INVALID;
	}
	
	HRESULT hr;
	hr = inkw.collector->putref_Ink(ink);
	if (FAILED(hr)) {
		return MAKE_API_ERROR(hr);
	}
	return 0;
}


IInkDisp *
TPC_get_ink(TPC_InkWindow &ink)
{
    if (!ink.collector) {
        errval = E_INVALID;
        return 0;
    }
    
    IInkDisp *ret;
    
    HRESULT hr;
    
    hr = ink.collector->get_Ink(&ret);
    if (FAILED(hr)) {
        errval = MAKE_API_ERROR(hr);
        return 0;
    }
    
    return ret;
}


IInkRenderer *
TPC_get_renderer(TPC_InkWindow &ink)
{
    if (!ink.renderer) {
        if (!ink.collector) {
            errval = E_INVALID;
            return 0;
        }
        
        HRESULT hr;
        
        hr = ink.collector->get_Renderer(&ink.renderer);
        if (FAILED(hr)) {
            errval = MAKE_API_ERROR(hr);
            return 0;
        }
    }
    
    return ink.renderer;
}


IInkOverlay *
TPC_get_overlay(TPC_InkWindow &ink)
{
    if (!ink.overlay) {
        if (!ink.collector) {
            errval = E_INVALID;
            return 0;
        }
        
        HRESULT hr;
        
        hr = CoCreateInstance(CLSID_InkOverlay, NULL, CLSCTX_INPROC_SERVER, IID_IInkOverlay, (void **)&ink.overlay);
        if (FAILED(hr)) {
            errval = MAKE_API_ERROR(hr);
            return 0;
        }
                
        HWND host;
        hr = ink.collector->get_hWnd(reinterpret_cast<LONG_PTR *>(&host));
        if (FAILED(hr)) {
            errval = MAKE_API_ERROR(hr);
            return 0;
        }
        
        hr = ink.overlay->put_hWnd(reinterpret_cast<LONG_PTR>(host));
        if (FAILED(hr)) {
            errval = MAKE_API_ERROR(hr);
            return 0;
        }
        
        IInkDisp *APIink = TPC_get_ink(ink);
        if (!APIink) {
            return 0;
        }
                
        hr = ink.overlay->putref_Ink(APIink);
        if (FAILED(hr)) {
            errval = MAKE_API_ERROR(hr);
            return 0;
        }
        
    }
    
    return ink.overlay;
}


int
TPC_pixel_to_inkspace(TPC_InkWindow &ink, long &x, long &y)
{
    HRESULT hr;
    
    if (!ink.collector) {
        return E_INVALID;
    }
    
    
    HWND host;
    hr = ink.collector->get_hWnd(reinterpret_cast<long *>(&host));
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    
	HDC dc = GetDC(host);
	if (!dc) {
		return MAKE_API_ERROR(GetLastError());
	}
	
	
	IInkRenderer *renderer = TPC_get_renderer(ink);
	
	hr = renderer->PixelToInkSpace(reinterpret_cast<LONG_PTR>(dc), &x, &y);
	if (FAILED(hr)) {
        ReleaseDC(host, dc);
    	return MAKE_API_ERROR(hr);
	}
	
	if (!ReleaseDC(host, dc)) {
		return MAKE_API_ERROR(GetLastError());
	}
	
	return 0;
}


int
TPC_inkspace_to_pixel(TPC_InkWindow &ink, long &x, long &y)
{
    HRESULT hr;
    
    if (!ink.collector) {
        return E_INVALID;
    }
    
    
    HWND host;
    hr = ink.collector->get_hWnd(reinterpret_cast<long *>(&host));
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    
	HDC dc = GetDC(host);
	if (!dc) {
		return MAKE_API_ERROR(GetLastError());
	}
	
	
	IInkRenderer *renderer = TPC_get_renderer(ink);
	
	hr = renderer->InkSpaceToPixel(reinterpret_cast<LONG_PTR>(dc), &x, &y);
	if (FAILED(hr)) {
        ReleaseDC(host, dc);
    	return MAKE_API_ERROR(hr);
	}
	
	if (!ReleaseDC(host, dc)) {
		return MAKE_API_ERROR(GetLastError());
	}
	
	return 0;
}


int
TPC_clear_ink(TPC_InkWindow &ink)
{
    if (!ink.collector) {
        return E_INVALID;
    }
    
    IInkDisp *APIink = TPC_get_ink(ink);
    if (!APIink) {
        int e = errval;
        errval = 0;
        return e;
    }
    
    HRESULT hr;
    hr = APIink->DeleteStrokes();
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    
    APIink->Release();
    
    return 0;
}


int
TPC_rescale_window(TPC_InkWindow &ink)
{
	if (!ink.collector) {
        return E_INVALID;
	}
	
	HRESULT hr;
	
	HWND host;
	hr = ink.collector->get_hWnd(reinterpret_cast<long *>(&host));
	if (FAILED(hr)) {
	    return MAKE_API_ERROR(hr);
	}
	
	RECT bbox;
	if (GetClientRect(host, &bbox) == FALSE) {
		return MAKE_API_ERROR(GetLastError());
	}
	
	bbox.left += 5;
	bbox.right -= 5;
	bbox.top += 5;
	bbox.bottom -= 5;

    TPC_pixel_to_inkspace(ink, bbox.left, bbox.top);
    TPC_pixel_to_inkspace(ink, bbox.right, bbox.bottom);
    
	IInkRectangle *ink_bbox;
	hr = CoCreateInstance(CLSID_InkRectangle, NULL, CLSCTX_INPROC_SERVER, IID_IInkRectangle, reinterpret_cast<void **>(&ink_bbox));
	if (FAILED(hr)) {
	    return MAKE_API_ERROR(hr);
	}
		
	IInkDisp *APIink = TPC_get_ink(ink);
	
	IInkStrokes *APIstrokes;
	hr = APIink->get_Strokes(&APIstrokes);
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    	
	IInkRectangle *ink_src_rect;
	hr = APIstrokes->GetBoundingBox(IBBM_PointsOnly, &ink_src_rect);
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    	
	RECT src_rect;
	ink_src_rect->get_Bottom(&src_rect.bottom);
	ink_src_rect->get_Top(&src_rect.top);
	ink_src_rect->get_Left(&src_rect.left);
	ink_src_rect->get_Right(&src_rect.right);
	
	ink_src_rect->Release();
	
    if (src_rect.bottom - src_rect.top > src_rect.right - src_rect.left) {
        double ratio = static_cast<double>(src_rect.right - src_rect.left) / (src_rect.bottom - src_rect.top);
        hr = ink_bbox->SetRectangle(bbox.top, bbox.left, bbox.bottom, static_cast<long>(bbox.left + (bbox.right - bbox.left) * ratio));
        if (FAILED(hr)) {
            return MAKE_API_ERROR(hr);
        }
    }
    else {
        double ratio = static_cast<double>(src_rect.bottom - src_rect.top) / (src_rect.right - src_rect.left);
        hr = ink_bbox->SetRectangle(bbox.top, bbox.left, static_cast<long>(bbox.top + (bbox.bottom - bbox.top) * ratio), bbox.right);
        if (FAILED(hr)) {
            return MAKE_API_ERROR(hr);
        }
    }

	hr = APIstrokes->ScaleToRectangle(ink_bbox);
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
	
    APIstrokes->Release();
    APIink->Release();
    ink_bbox->Release();
    
	return 0;
}


int
TPC_get_strokes(TPC_InkWindow &ink, TPC_StrokeGroup &strokes)
{   
	HRESULT hr;

    if (!ink.collector) {
        return E_INVALID;
    }
    
    
    IInkStrokes *APIstrokes;
    
	IInkDisp *APIink;
	hr = ink.collector->get_Ink(&APIink);
	if (FAILED(hr)) {
	    return MAKE_API_ERROR(hr);
	}
	
	hr = APIink->get_Strokes(&APIstrokes);
	if (FAILED(hr)) {
	    APIink->Release();
	    return MAKE_API_ERROR(hr);
	}
	
	
	long num_strokes;
	hr = APIstrokes->get_Count(&num_strokes);
    if (FAILED(hr)) {
        APIstrokes->Release();
        APIink->Release();
        return MAKE_API_ERROR(hr);
    }

    IInkStrokeDisp **TPCstrokes;
    
    TPCstrokes = new IInkStrokeDisp *[num_strokes];
    if (!TPCstrokes) {
        APIstrokes->Release();
        APIink->Release();
        return E_OUTOFMEM;
    }
        
	for (long i = 0; i < num_strokes; i++) {
		IInkStrokeDisp *stroke;

		hr = APIstrokes->Item(i, &stroke);
		if (FAILED(hr)) {
		    APIstrokes->Release();
		    APIink->Release();
            return MAKE_API_ERROR(hr);
		}
		
		TPCstrokes[i] = stroke;
	}
	
	strokes.set_strokes(TPCstrokes, num_strokes);
	
	return 0;
}


int
TPC_add_strokes(TPC_InkWindow &ink, TPC_StrokeGroup &strokes)
{
    HRESULT hr;
    
    if (!ink.collector) {
        return E_INVALID;
    }
    
    IInkDisp *APIink = TPC_get_ink(ink);
    
    IInkStrokes *APIstrokes;
    
    VARIANT null;
    VariantInit(&null);
    
    hr = APIink->CreateStrokes(null, &APIstrokes);
    if (FAILED(hr)) {
        APIink->Release();
        return MAKE_API_ERROR(hr);
    }
    
    VariantClear(&null);
    
//    IInkRectangle *APIbbox;
    for (IInkStrokeDisp **stroke = strokes.strokes; stroke != strokes.strokes + num_strokes(strokes); stroke++) {
        hr = APIstrokes->Add(*stroke);
        if (FAILED(hr)) {
            APIink->Release();
            APIstrokes->Release();
            return MAKE_API_ERROR(hr);
        }        
    }

/*    
    hr = APIstrokes->GetBoundingBox(IBBM_Default, &APIbbox);
    if (FAILED(hr)) {
        APIstrokes->Release();
        APIink->Release();
        return MAKE_API_ERROR(hr);
    }

    hr = APIink->AddStrokesAtRectangle(APIstrokes, APIbbox);
    if (FAILED(hr)) {
        APIbbox->Release();
        APIstrokes->Release();
        APIink->Release();
        return MAKE_API_ERROR(hr);
    }

    */
    //APIbbox->Release();
    APIstrokes->Release();
    APIink->Release();
    
    return 0;
}


IInkStrokes *
TPC_get_APIstrokes(TPC_InkWindow &ink)
{
    IInkDisp *APIink = TPC_get_ink(ink);
    
    if (!APIink) {
        return 0;
    }
    
    IInkStrokes *APIstrokes;
    
    HRESULT hr;
    
    hr = APIink->get_Strokes(&APIstrokes);
    if (FAILED(hr)) {
        APIink->Release();
        return 0;
    }
    
    APIink->Release();
    
    return APIstrokes;
}


}

#endif
