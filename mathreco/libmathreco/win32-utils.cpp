#ifdef WIN32

#include "utils.h"
#include "strutils.h"

#define NOMINMAX
#include <msinkaut.h>
#include <windows.h>
#include <cstring>

#include "error.h"
#include "tpc-group.h"


namespace scg
{


static std::string profile_path;
static std::wstring profile_wpath;
static std::wstring user_profile_wpath;

int
GetProfilePathW(std::wstring &path) {
	if (!profile_wpath.empty()) {
	  path = profile_wpath;
	  return 0;
	}
	
	std::string npath;
	int e = GetProfilePath(npath);
	if (FAILURE(e)) {
		return e;
	}
	
	profile_wpath = str2wstr(npath);
	path = profile_wpath;
	return 0;
}

DLLDECL void
SetProfilePath(const char *path) {
	profile_path = path;
}

int
GetUserProfilePathW(std::wstring &path) {
	if (!user_profile_wpath.empty()) {
	  path = user_profile_wpath;
	  return 0;
	}
	
	std::string npath;
	int e = GetUserProfilePath(npath);
	if (FAILURE(e)) {
		return e;
	}
	
	user_profile_wpath = str2wstr(npath);
	path = user_profile_wpath;
	return 0;
}

int
GetProfilePath(std::string &path)
{
	if (!profile_path.empty()) {
		path = profile_path;
		return 0;
	}
	
#ifdef _WINDLL // only check DLL directory for DLL version
	HMODULE module = ::GetModuleHandleA("libmathreco-dll");
	if (!module) {
		goto steal_svn_path;
	}
	
	size_t maxsz = 512;
	char *filename = new char[maxsz];
	DWORD sz;
	do {
		sz = ::GetModuleFileNameA(module, filename, static_cast<DWORD>(maxsz));
		if (sz == 0) {
			delete[] filename;
			if (::GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
				goto steal_svn_path;
			}
			maxsz *= 2;
			filename = new char[maxsz];
		}
	} while (sz == 0);
	
	// trim dll filename from path
	char *p = std::strrchr(filename, '\\');
	if (!p) {
		goto steal_svn_path;
	}
	*p = '\0';
	
	profile_path = filename;
	delete[] filename;
	profile_path += "\\data-files";

	DWORD attribs = ::GetFileAttributesA(profile_path.c_str());
	if ((attribs != INVALID_FILE_ATTRIBUTES) && (attribs & FILE_ATTRIBUTE_DIRECTORY)) {
		path = profile_path;
		return 0;
	}
#endif

steal_svn_path:

	char *svn_path = std::getenv("SVNPEN");
	if (!svn_path) {
		ERR(E_NOTFOUND, "The SVNPEN environment variable must be set.");
		return E_NOTFOUND;
	}
	
	profile_path = svn_path;
	profile_path += "\\Code\\C++\\mathreco\\TRUNK\\src\\data-files\\";
	path = profile_path;
	return 0;
}

int
ReadMicrosoftStrokeGroupInk(std::istream &is, RawStrokeGroup &strokes, unsigned ink_sz)
{
    HRESULT hr;
    
    IInkDisp *tpc_ink;
    
	hr = CoCreateInstance(CLSID_InkDisp, NULL, CLSCTX_INPROC_SERVER, IID_IInkDisp, (void **)&tpc_ink);
	if (FAILED(hr)) {
	    return MAKE_API_ERROR(hr);
	}

    SAFEARRAY *ink_buf = SafeArrayCreateVector(VT_UI1, 0, static_cast<ULONG>(ink_sz));
	if (!ink_buf) {
	    return MAKE_API_ERROR(GetLastError());
	}
    
	char *inkdata;
	hr = SafeArrayAccessData(ink_buf, reinterpret_cast<void HUGEP **>(&inkdata));
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    
    is.read(inkdata, ink_sz);
    if (is.bad() || is.fail()) {
        return E_IO;
    }
    

    VARIANT buf;
    VariantInit(&buf);
    buf.vt = VT_ARRAY | VT_UI1;
    buf.parray = ink_buf;
    
    hr = tpc_ink->Load(buf);
    
    if (FAILED(hr)) {
        tpc_ink->Release();   
        return MAKE_API_ERROR(hr);
    }

    // NOTE!!! The order of the following three operations is CRUCIAL! Any
    // other order and the program will abort!  (Why is this? Windows is dumb.)
    VariantClear(&buf);
 	hr = SafeArrayUnaccessData(ink_buf);
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    SafeArrayDestroy(ink_buf);    

    IInkStrokes *tpc_strokes;
    hr = tpc_ink->get_Strokes(&tpc_strokes);
    if (FAILED(hr)) {
        tpc_ink->Release();
        return MAKE_API_ERROR(hr);
    }
    
    LONG nstrokes;
    hr = tpc_strokes->get_Count(&nstrokes);
    if (FAILED(hr)) {
        tpc_strokes->Release();
        tpc_ink->Release();
        return MAKE_API_ERROR(hr);
    }
    
    IInkStrokeDisp **tpc_stroke_ptr = new IInkStrokeDisp *[nstrokes];
    if (!tpc_stroke_ptr) {
        tpc_strokes->Release();
        tpc_ink->Release();
        return E_OUTOFMEM;
    }
    
    for (LONG i = 0; i < nstrokes; i++) {
        hr = tpc_strokes->Item(i, tpc_stroke_ptr + i);
        if (FAILED(hr)) {
            delete[] tpc_stroke_ptr;
            tpc_strokes->Release();
            tpc_ink->Release();
            return MAKE_API_ERROR(hr);
        }
    }
    
    TPC_StrokeGroup tpc_group(tpc_stroke_ptr, static_cast<unsigned>(nstrokes));
    
    int err = convert(strokes, tpc_group);
    if (FAILURE(err)) {
        tpc_strokes->Release();
        tpc_ink->Release();
        return err;
    }
    
    tpc_strokes->Release();
    tpc_ink->Release();
    
    return 0;
}


int
ReadMicrosoftStrokeGroupInk(std::istream &is, RawStrokeGroup &strokes)
{
    HRESULT hr;
    
    IInkDisp *tpc_ink;
    
	hr = CoCreateInstance(CLSID_InkDisp, NULL, CLSCTX_INPROC_SERVER, IID_IInkDisp, (void **)&tpc_ink);
	if (FAILED(hr)) {
	    return MAKE_API_ERROR(hr);
	}

    unsigned ink_sz;
    is.read(reinterpret_cast<char *>(&ink_sz), sizeof(ink_sz));
    if (is.fail()) {
        return E_IO;
    }
    
    SAFEARRAY *ink_buf = SafeArrayCreateVector(VT_UI1, 0, static_cast<ULONG>(ink_sz));
	if (!ink_buf) {
	    return MAKE_API_ERROR(GetLastError());
	}
    
	char *inkdata;
	hr = SafeArrayAccessData(ink_buf, reinterpret_cast<void HUGEP **>(&inkdata));
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    
    is.read(inkdata, ink_sz);
    if (is.bad() || is.fail()) {
        return E_IO;
    }
    

    VARIANT buf;
    VariantInit(&buf);
    buf.vt = VT_ARRAY | VT_UI1;
    buf.parray = ink_buf;
    
    hr = tpc_ink->Load(buf);
    
    if (FAILED(hr)) {
        tpc_ink->Release();   
        return MAKE_API_ERROR(hr);
    }

    // NOTE!!! The order of the following three operations is CRUCIAL! Any
    // other order and the program will abort!  (Why is this? Windows is dumb.)
    VariantClear(&buf);
 	hr = SafeArrayUnaccessData(ink_buf);
    if (FAILED(hr)) {
        return MAKE_API_ERROR(hr);
    }
    SafeArrayDestroy(ink_buf);    

    IInkStrokes *tpc_strokes;
    hr = tpc_ink->get_Strokes(&tpc_strokes);
    if (FAILED(hr)) {
        tpc_ink->Release();
        return MAKE_API_ERROR(hr);
    }
    
    LONG nstrokes;
    hr = tpc_strokes->get_Count(&nstrokes);
    if (FAILED(hr)) {
        tpc_strokes->Release();
        tpc_ink->Release();
        return MAKE_API_ERROR(hr);
    }
    
    IInkStrokeDisp **tpc_stroke_ptr = new IInkStrokeDisp *[nstrokes];
    if (!tpc_stroke_ptr) {
        tpc_strokes->Release();
        tpc_ink->Release();
        return E_OUTOFMEM;
    }
    
    for (LONG i = 0; i < nstrokes; i++) {
        hr = tpc_strokes->Item(i, tpc_stroke_ptr + i);
        if (FAILED(hr)) {
            delete[] tpc_stroke_ptr;
            tpc_strokes->Release();
            tpc_ink->Release();
            return MAKE_API_ERROR(hr);
        }
    }
    
    TPC_StrokeGroup tpc_group(tpc_stroke_ptr, static_cast<unsigned>(nstrokes));
    
    int err = convert(strokes, tpc_group);
    if (FAILURE(err)) {
        tpc_strokes->Release();
        tpc_ink->Release();
        return err;
    }
    
    tpc_strokes->Release();
    tpc_ink->Release();
    
    return 0;
}


}

#endif
