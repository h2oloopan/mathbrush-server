#define _WIN32_WINNT 0x0500
#define NOMINMAX
#include <windows.h>

#include <comdef.h>
#include <msinkaut.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <iostream>
#include <clocale>
#include <sstream>
#include <time.h>

#include "resource.h"

#include "annotate.h"
#include "dist.h"
#include "error.h"
#include "group.h"
#include "ink-io.h"
#include "stroke.h"
#include "tabletpc.h"
#include "utils.h"
#include "segment.h"
#include "strutils.h"
#include "MathRecognizer.h"
#include "interp.h"
#include "normalize.h"
#include "caclient.h"
//#include "MatrixAnalyzer.h"
//#include "RangeExpander.h"
//#include "functions.h"

static scg::MathRecognizer *rec = 0;


enum LinkType {
    Above,
    Left,
    BelowLeft,
    AboveLeft,
    Inside
};

struct LinkAnnotation
{
    unsigned box1;
    unsigned box2;
    LinkType link;
};


std::vector<LinkAnnotation> link_annotations;

scg::AnnotatedStrokeGroup annotated_strokes;

static int have_link_annotations = 0;
static bool have_symbol_annotations = false;
static bool annotation_mode = false;
static bool USE_PARSER_CHECKED = false;


void
show_error(const std::string &msg)
{
	std::wstring ws = scg::str2wstr(msg);
	MessageBox(0, ws.c_str(), L"error", MB_OK | MB_ICONSTOP);
}


static LRESULT
GetNameDlgProc(HWND dlg, UINT msg, WPARAM wparam, LPARAM lparam)
{
	switch (msg) {
	case WM_INITDIALOG:
		return TRUE;

	case WM_COMMAND:
		if (wparam == IDOK) {
			wchar_t *s = new wchar_t[64];
			if (GetDlgItemText(dlg, IDC_NAME, s, 64) == 0) {
				MessageBox(dlg, L"You must enter a name.", L"Goofy user", MB_OK);
			}
			else {
				EndDialog(dlg, reinterpret_cast<int>(s));
			}
		}
		return TRUE;
	}

	return FALSE;
}


int
SaveInk(const std::string &filename, scg::TPC_InkWindow &ink)
{
    int e;
    
	std::ofstream os(filename.c_str());
	if (!os) {
	    return E_IO;
	}

    if (!have_symbol_annotations) {
        scg::TPC_StrokeGroup TPCgroup;
    
        e = scg::TPC_get_strokes(ink, TPCgroup);
        if (FAILURE(e)) {
            return e;
        }
    
        scg::RawStrokeGroup group;
        e = scg::convert(group, TPCgroup);
        if (FAILURE(e)) {
            return e;
        }
    
        os << group;
    }
    else {
        scg::export_annotated_ink(os, annotated_strokes);
        
        /*if (have_link_annotations) {
            os << "\n";
            os << link_annotations.size() << std::endl;
            for (std::vector<LinkAnnotation>::const_iterator i = link_annotations.begin(); i != link_annotations.end(); ++i) {
                unsigned b = 0;
                for (unsigned s = 0; s < i->box1; s++) {
                    if (annotated_strokes.annotations[s].nstrokes > 0) {
                        b++;
                    }
                }
                os << b << " ";
                b = 0;
                for (unsigned s = 0; s < i->box2; s++) {
                    if (annotated_strokes.annotations[s].nstrokes > 0) {
                        b++;
                    }
                }
                os << b << " ";
                //os << i->box1 << " " << i->box2 << " ";
                switch (i->link) {
                case BelowLeft:
                    os << "BL" << std::endl;
                    break;
                case AboveLeft:
                    os << "BR" << std::endl;
                    break;
                case Above:
                    os << "A" << std::endl;
                    break;
                case Left:
                    os << "L" << std::endl;
                    break;
                case Inside:
                    os << "C" << std::endl;
                }
            }
        }*/
    }
    
	return 0;
}

int
LoadIpadInk(const std::string &filename, scg::TPC_InkWindow &ink) {
    int e;
	std::ifstream is(filename.c_str());
	if (!is) {
		return E_IO;
	}

    scg::RawStrokeGroup group;

	is >> group;
	long scale = 2540.0/132.0;//std::max(bounds.width(), bounds.height());
    //group = subdivide(annotated_strokes);//.copy();
	scg::RawStroke *strokes = new scg::RawStroke[group.nstrokes];
	for (size_t i = 0; i < group.nstrokes; ++i) {
		scg::RawStroke &stk = group.strokes[i];
		//const scg::NormalizedStroke &stk = ss.strokes[i];
		long *x = new long[stk.npoints], *y = new long[stk.npoints];
		for (size_t j = 0; j < stk.npoints; ++j) {
			//stk.x[j] = (long)(stk.x[j] * 2540.0/scg::TABLETPC_DPI);
			//stk.y[j] = (long)(stk.y[j] * 2540.0/scg::TABLETPC_DPI);
			x[j] = (long)(stk.x[j] * scale);
			y[j] = (long)(stk.y[j] * scale);
		}
		strokes[i].set_points(x, y, stk.npoints);
	}
	group.set_strokes(strokes, group.nstrokes);
    
    e = scg::TPC_clear_ink(ink);
    if (FAILURE(e)) {
        return e;
    }
    
    scg::TPC_StrokeGroup TPCgroup;
    e = scg::convert(TPCgroup, group, ink);
    if (FAILURE(e)) {
        return e;
    }
    
    return scg::TPC_add_strokes(ink, TPCgroup);
}
int
LoadInk(const std::string &filename, scg::TPC_InkWindow &ink)
{
    int e;
    
    int mode = std::ios_base::in;
    bool is_binary = (filename[filename.length() - 1] == 'f') || (filename.length() > 5 && filename.substr(filename.length() - 5, 5) == "msink");
    if (is_binary) {
        mode |= std::ios_base::binary;
    }
    
	std::ifstream is(filename.c_str(), mode);
	if (!is) {
		return E_IO;
	}

    scg::NormalizedStrokeGroup group;
	scg::RawStrokeGroup raw_group;
    if (is_binary) {
        std::istream::pos_type start_of_data = is.tellg();
        is.seekg(0, std::ios_base::end);
        std::istream::pos_type end_of_data = is.tellg();
        is.seekg(start_of_data);
	    e = scg::ReadMicrosoftStrokeGroupInk(is, raw_group, end_of_data - start_of_data);
	    if (FAILURE(e)) {
	        return e;
	    }
    }
    else {
        /*static const std::string FileHeader("SCG_INK");
        unsigned i = 0;
        while (is.peek() == FileHeader[i]) {
            i++;
            is.get();
        }
        if (i > 0 && i != FileHeader.length()) {
            return E_INVALID;
        }*/
        
        //scg::import_annotated_ink(is, annotated_strokes);
		is >> group;
		//scg::NormalizedStrokeGroup ns = normalize(annotated_strokes);
		//scg::NormalizedStrokeGroup ss = subdivide(ns);
		//scg::Rect<long> bounds = annotated_strokes.bounds();
		double dpi = 2540;
		double scale = 2540.0/dpi;//std::max(bounds.width(), bounds.height());
		double dx = 0;//-20000;
		double dy = 0;//-60000;
        //group = subdivide(annotated_strokes);//.copy();
		scg::RawStroke *strokes = new scg::RawStroke[group.nstrokes];
		for (size_t i = 0; i < group.nstrokes; ++i) {
			scg::NormalizedStroke &stk = group.strokes[i];
			//const scg::NormalizedStroke &stk = ss.strokes[i];
			long *x = new long[stk.npoints], *y = new long[stk.npoints];
			for (size_t j = 0; j < stk.npoints; ++j) {
				//stk.x[j] = (long)(stk.x[j] * 2540.0/scg::TABLETPC_DPI);
				//stk.y[j] = (long)(stk.y[j] * 2540.0/scg::TABLETPC_DPI);
				x[j] = (long)((stk.x[j]+dx) * scale);
				y[j] = (long)((stk.y[j]+dy) * scale);
			}
			strokes[i].set_points(x, y, stk.npoints);
		}
		raw_group.set_strokes(strokes, group.nstrokes);
		
		//group.translate(-1000, -1000);
		//group.translate(10000, 10000);
        /*if (annotated_strokes.annotations[0].nstrokes > 0) {
            have_symbol_annotations = true;
        }
        else {
            have_symbol_annotations = false;
            have_link_annotations = false;
        }

        link_annotations.clear();
                
        if (have_symbol_annotations) {
            unsigned nboxes;
            is >> nboxes;
            
            if (!is.eof()) {
                while (nboxes-- > 0) {
                    LinkAnnotation annot;
                    std::string type;
                    unsigned b;
                    is >> b;
                    b++;
                    for (unsigned i = 0; i < num_strokes(annotated_strokes); i++) {
                        if (annotated_strokes.annotations[i].nstrokes) {
                            b--;
                            if (b == 0) {
                                annot.box1 = i;
                                break;
                            }
                        }
                    }
                    
                    is >> b;
                    b++;
                    for (unsigned i = 0; i < num_strokes(annotated_strokes); i++) {
                        if (annotated_strokes.annotations[i].nstrokes) {
                            b--;
                            if (b == 0) {
                                annot.box2 = i;
                                break;
                            }
                        }
                    }
                    
                    is >> type;
                    
                    if (type == "A") {
                        annot.link = Above;
                    }
                    else if (type == "L") {
                        annot.link = Left;
                    }
                    else if (type == "BR") {
                        annot.link = AboveLeft;
                    }
                    else if (type == "BL") {
                        annot.link = BelowLeft;
                    }
                    else if (type == "C") {
                        annot.link = Inside;
                    }
                    
                    link_annotations.push_back(annot);
                }
                have_link_annotations = true;
            }
            else {
                have_link_annotations = false;
            }
        }*/
    }
    
    e = scg::TPC_clear_ink(ink);
    if (FAILURE(e)) {
        return e;
    }
    
    scg::TPC_StrokeGroup TPCgroup;
    e = scg::convert(TPCgroup, raw_group, ink);
    if (FAILURE(e)) {
        return e;
    }
    
    return scg::TPC_add_strokes(ink, TPCgroup);
}


static int
Recognize(scg::TPC_InkWindow &ink)
{
	return 0;
}


static void
RunMicrosoftRecognizer(scg::TPC_InkWindow &ink)
{
	// run microsoft's tablet sdk recognizer on the input
	HRESULT hr;

	IInkRecognizers *recognizers;
	hr = CoCreateInstance(CLSID_InkRecognizers, NULL, CLSCTX_INPROC_SERVER, IID_IInkRecognizers, (void **)&recognizers);
	if (FAILED(hr)) {
		show_error("couldn't enumerate the installed recognizers");
		return;
	}
	
	long n;
	hr = recognizers->get_Count(&n);
	if (FAILED(hr)) {
		show_error("couldn't enumerate the installed recognizers");
		return;
	}

	if (n == 0) {
		show_error("no Microsoft recognizers are installed");
		return;
	}
	
	
	IInkRecognizerContext *recognizer;
	hr = CoCreateInstance(CLSID_InkRecognizerContext, NULL, CLSCTX_INPROC_SERVER, IID_IInkRecognizerContext, (void **)&recognizer);
	if (FAILED(hr)) {
		show_error("could not instantiate a Microsoft recognizer");
		return;
	}

	IInkStrokes *strokes;
	IInkDisp *APIink = scg::TPC_get_ink(ink);
	hr = APIink->get_Strokes(&strokes);
	if (FAILED(hr)) {
		show_error("couldn't get ink strokes");
		return;
	}

    BSTR tmp = SysAllocString(L"NONE");
    
    hr = recognizer->put_Factoid(tmp);
    if (FAILED(hr)) {
		show_error("couldn't set factoid");
		return;
    }
    
	hr = recognizer->putref_Strokes(strokes);
	if (FAILED(hr)) {
		show_error("couldn't tell recognizer what to recognize");
		return;
	}
	
	InkRecognitionStatus status = IRS_NoError;
	IInkRecognitionResult *result;
	hr = recognizer->Recognize(&status, &result);
	if (FAILED(hr)) {
		show_error("Microsoft recognition failed because the cows are sick");
		return;
	}
	
	BSTR match;
	hr = result->get_TopString(&match);
	if (FAILED(hr)) {
		show_error("couldn't extract matched string from recognizer");
		return;
	}
	
	// just dump the matched string to the screen...it's not really worth the work to figure out
	// which character matches which set of input strokes
	HDC dc = GetDC(GetActiveWindow());
	RECT dim = {10, 10, 200, 20};
	DrawTextW(dc, match, SysStringLen(match), &dim, DT_CALCRECT);
	DrawTextW(dc, match, SysStringLen(match), &dim, 0);
	SysFreeString(match);
	ReleaseDC(GetActiveWindow(), dc);
	
	result->Release();
	recognizer->Release();
}


static IInkRecognizerContext *tablet_reco = 0;


static void
RunTabletRecognizer(scg::TPC_InkWindow &ink)
{
	HRESULT hr;
	
	if (!tablet_reco) {
	    IInkRecognizers *recognizers;
	    hr = CoCreateInstance(CLSID_InkRecognizers, NULL, CLSCTX_INPROC_SERVER, IID_IInkRecognizers, (void **)&recognizers);
	    if (FAILED(hr)) {
		    show_error("couldn't enumerate the installed recognizers");
		    return;
	    }
    	
	    long n;
	    hr = recognizers->get_Count(&n);
	    if (FAILED(hr)) {
		    show_error("couldn't enumerate the installed recognizers");
		    return;
	    }

	    if (n == 0) {
		    show_error("no Microsoft recognizers are installed");
		    return;
	    }


        BSTR desired_recognizer = SysAllocString(OLESTR("UW_SCG_MathRecognizer-Research"));
        
        IInkRecognizer *rec;
	    while (n) {
	        n--;
	        recognizers->Item(n, &rec);
	        BSTR name;
	        rec->get_Name(&name);                                                                               
	        if (VarBstrCmp(name, desired_recognizer, GetUserDefaultLCID(), 0) == VARCMP_EQ) {
	            break;
	        }
	    }
    	
	    SysFreeString(desired_recognizer);
    	
    	
	    hr = rec->CreateRecognizerContext(&tablet_reco);
	    if (FAILED(hr)) {
		    show_error("could not instantiate a MathRecognizer context");
		    return;
	    }
	    
	    rec->Release();
    }
    
	IInkStrokes *strokes;
	IInkDisp *APIink = scg::TPC_get_ink(ink);
	hr = APIink->get_Strokes(&strokes);
	if (FAILED(hr)) {
		show_error("couldn't get ink strokes");
		return;
	}
	
	hr = tablet_reco->putref_Strokes(strokes);
	if (FAILED(hr)) {
		show_error("couldn't tell recognizer what to recognize");
		return;
	}
	
	InkRecognitionStatus status = IRS_NoError;
	IInkRecognitionResult *result;
	hr = tablet_reco->Recognize(&status, &result);
	if (FAILED(hr)) {
		show_error("Tablet recognition failed because the pigs are ill");
		return;
	}
	
	if (result) {
	    BSTR match;
	    hr = result->get_TopString(&match);
	    if (FAILED(hr)) {
		    show_error("couldn't extract matched string from recognizer");
		    return;
	    }
        
        HDC dc = GetDC(GetActiveWindow());
    	
	    // just dump the matched string to the screen...it's not really worth the work to figure out
	    // which character matches which set of input strokes
	    RECT dim = {10, 10, 200, 20};
	    DrawTextW(dc, match, SysStringLen(match), &dim, DT_CALCRECT);
	    DrawTextW(dc, match, SysStringLen(match), &dim, 0);
	    SysFreeString(match);
/*
	    int xmm = GetDeviceCaps(dc, HORZSIZE);
	    int ymm = GetDeviceCaps(dc, VERTSIZE);
	    int xpels = GetDeviceCaps(dc, HORZRES);
	    int ypels = GetDeviceCaps(dc, VERTRES);
*/
        UINT len = SysStringLen(match);
        for (unsigned i = 0; i < len; i++) {
            IInkRecognitionAlternates *alts;
            hr = result->AlternatesFromSelection(i, 1, 5, &alts);
            if (FAILED(alts)) {
                break;
            }

            IInkStrokes *strokes;
            hr = alts->get_Strokes(&strokes);
            if (FAILED(hr)) {
                alts->Release();
                break;
            }
            
            long nstrokes;
            hr = alts->get_Count(&nstrokes);
            if (FAILED(hr)) {
                strokes->Release();
                alts->Release();
                break;
            }
            
            IInkRectangle *bounds;
            hr = strokes->GetBoundingBox(IBBM_Default, &bounds);
            if (FAILED(hr)) {
                strokes->Release();
                alts->Release();
                break;
            }
            
            RECT rc;
            bounds->get_Left(&rc.left);
            bounds->get_Top(&rc.top);
            bounds->get_Right(&rc.right);
            bounds->get_Bottom(&rc.bottom);
            
            scg::TPC_inkspace_to_pixel(ink, rc.left, rc.top);
            scg::TPC_inkspace_to_pixel(ink, rc.right, rc.bottom);
            
            FrameRect(dc, &rc, (HBRUSH)GetStockObject(GRAY_BRUSH));

            bounds->Release();
            strokes->Release();
            
            IInkRecognitionAlternate *top_alt;
            hr = alts->Item(0, &top_alt);
            if (FAILED(hr)) {
                alts->Release();
                break;
            }
            
	        BSTR char_name;
	        hr = top_alt->get_String(&char_name);
	        if (FAILED(hr)) {
	            top_alt->Release();
	            alts->Release();
	            break;
	        }
	        
	        RECT dim = {rc.right, rc.bottom, 200, 20};
	        DrawTextW(dc, char_name, SysStringLen(char_name), &dim, DT_CALCRECT);
	        DrawTextW(dc, char_name, SysStringLen(char_name), &dim, 0);
	        
	        SysFreeString(char_name);
            top_alt->Release();
            alts->Release();
        }
        
	    ReleaseDC(GetActiveWindow(), dc);
    	SysFreeString(match);
	    result->Release();
	}
}


static void
SetMenuRecognizer(HWND hwnd, int item)
{
	// add a check mark to the "item"th entry in menu "hwnd"
	
	HMENU menu = GetMenu(hwnd);
	menu = GetSubMenu(menu, 1);
	int n = GetMenuItemCount(menu);
	while (n--) {
		CheckMenuItem(menu, n - 1, MF_BYPOSITION | MF_UNCHECKED);
	}
	CheckMenuItem(menu, item, MF_BYCOMMAND | MF_CHECKED);

	DrawMenuBar(hwnd);
}


/*
static int
RunFancySubdivision(HWND hwnd, scg::TPC_InkWindow ink)
{
    HRESULT hr;
    IInkStrokes *strokes;
    IInkDisp *APIink = scg::TPC_get_ink(ink);
    hr = APIink->get_Strokes(&strokes);
    if (FAILED(hr)) {
        show_error("couldn't get ink strokes");
        return MAKE_API_ERROR(hr);
    }
    
    IInkStrokeDisp *first_stroke;
    hr = strokes->Item(0, &first_stroke);
    if (FAILED(hr)) {
        show_error("couldn't extract first stroke");
        return MAKE_API_ERROR(hr);
    }
    
    scg::TPC_StrokeGroup tpc_strokes(&first_stroke, 1);
    
    int e;
    scg::RawStrokeGroup group;
    e = scg::convert(group, tpc_strokes);

    if (FAILURE(e)) {
        return e;
    }
    
    
    scg::NormalizedStroke stroke = scg::normalize(group.strokes()[0]);
    
    if (scg::num_points(stroke) <= 12) {
        return E_INVALID;
    }
    
    std::vector<unsigned> selected_points;
    selected_points.push_back(0);
    
    std::vector<double> xdiff;
    std::vector<double> ydiff;
    
    double *x, *y;
    double *endx = stroke.x + scg::num_points(stroke) - 1;
    for (x = stroke.x + 1, y = stroke.y + 1; x < endx; x++, y++) {
        double xd = (*(x + 1) - *x) - (*x - *(x - 1));
        double yd = (*(y + 1) - *y) - (*y - *(y - 1));
        double D = std::sqrt((*(x + 1) - *x) * (*(x + 1) - *x) + (*(y + 1) - *y) * (*(y + 1) - *y))
                 + std::sqrt((*x - *(x - 1)) * (*x - *(x - 1)) + (*y - *(y - 1)) * (*y - *(y - 1)));
        xdiff.push_back(xd / std::max(std::abs(*(x + 1) - *(x - 1)), 1.0));
        ydiff.push_back(yd / std::max(std::abs(*(y + 1) - *(y - 1)), 1.0));
    }
    
    double ceiling = std::numeric_limits<double>::infinity();
    
    double min_dist = 10.0;
    
    while (selected_points.size() != 9) {
        std::vector<double>::const_iterator xi, yi;
        double max = 0.0;
        unsigned max_index = std::numeric_limits<unsigned>::max();
        unsigned index = 0;
        for (xi = xdiff.begin(), yi = ydiff.begin(); xi != xdiff.end(); ++xi, ++yi) {
            std::vector<unsigned>::const_iterator i;
            for (i = selected_points.begin(); i != selected_points.end(); i++) {
                if (scg::dist_sq(stroke.x[*i], stroke.y[*i], stroke.x[index], stroke.y[index]) <= min_dist * min_dist) {
                    break;
                }
            }
            if (i == selected_points.end()) {
                double total = *xi + *yi;
                if (std::abs(total) > std::abs(max) && std::abs(total) < ceiling) {
                    max = total;
                    max_index = index;
                }
            }
            index++;
        }
        
        if (max_index == std::numeric_limits<unsigned>::max()) {
            min_dist *= 0.5;
        }
        else {
            selected_points.push_back(max_index + 1);
            ceiling = std::abs(max);
        }
    }
    
    selected_points.push_back(scg::num_points(stroke) - 1);
    
    HDC dc = GetDC(hwnd);
    
    for (std::vector<unsigned>::const_iterator i = selected_points.begin(); i != selected_points.end(); ++i) {
        long xp, yp;
        xp = group.strokes()[0].x[*i];
        yp = group.strokes()[0].y[*i];
        scg::TPC_inkspace_to_pixel(ink, xp, yp);
        
        Ellipse(dc, xp - 3, yp - 3, xp + 3, yp + 3);
    }
    
    ReleaseDC(hwnd, dc);
    
    tpc_strokes.strokes = 0;
    tpc_strokes.nstrokes = 0;
    
    return 0;
}
*/

static void
write_node(std::ostream &os, const scg::ExpressionTree *expr)
{
	os << "result: " << expr->score() << " with string " << expr->long_str() << std::endl;
	os << "string: " << expr->str() << std::endl;
	os << "latex: " << expr->latex_str() << std::endl;

	size_t nbuf;
	int e = expr->serialize(0, &nbuf);
	assert(e == 0);
	char *buf = new char[nbuf];
	e = expr->serialize(buf, &nbuf);
	assert(e == 0);
	std::ofstream fout("treeser", std::ios::binary);
	fout.write((char *)&nbuf, sizeof(nbuf));
	fout.write(buf, nbuf);
	fout.close();

	memset(buf, 0, nbuf);

	std::ifstream fin("treeser", std::ios::binary);
	fin.read((char *)&nbuf, sizeof(nbuf));
	fin.read(buf, nbuf);
	assert(fin.good());
	expr = scg::CreateExpressionTree(buf, nbuf);
	assert(expr);

	os << "AFTER SERIALIZATION:\n";
	os << "result: " << expr->score() << " with string " << expr->long_str() << std::endl;
	os << "string: " << expr->str() << std::endl;
	os << "latex: " << expr->latex_str() << std::endl;
	/*os << "node: type=" << expr->type() << "; string=" << expr->str() << " has " << expr->nchildren() << " children and scored " << expr->score() << "\n";
	for (unsigned kk = 0; kk < expr->nchildren(); ++kk) {
        write_node(os, expr->child(kk));
	}*/
}

static void
write_wide_node(std::wostream &wos, const scg::ExpressionTree *expr)
{
	wos << L"node: type=" << expr->type() << L"; string=";
	std::wstring ws = expr->wstr();
	wos.write(expr->wstr(), std::wcslen(expr->wstr()));
	wos << L" has " << expr->nchildren() << L" children\n";
}


/*static void
DoMatrixRecognition(scg::TPC_InkWindow &msink)
{
	//scg::SetVerbosity(2);
	//scg::VerboseOutputToFile("C:\\Documents and Settings\\Robert Amlard\\Desktop\\log.txt");

	//scg::Matrix m(NULL);
	//scg::MatrixElement* elem[6];
	//for (int i = 0; i < 6; i++) {
	//	elem[i] = new scg::MatrixElement;
	//	m.insert(*elem[i], i % 3, i % 2);
	//}
	//m.print();

	//m.remove(*elem[0]);
	//m.remove(*elem[3]);
	//m.print();

	scg::TPC_StrokeGroup TPCgroup;
	RECT draw_dim = {10, 10, 200, 200};
	HDC dc = GetDC(GetActiveWindow());

	int e = scg::TPC_get_strokes(msink, TPCgroup);
	if (FAILURE(e)) {
		THROW_LAST();
	}

	scg::RawStrokeGroup group;
	e = scg::convert(group, TPCgroup);
	if (FAILURE(e)) {
		THROW_LAST();
	}

	time_t start, end;

	start = time(NULL);
	scg::MathRecognizer* parser = scg::CreateMathRecognizer();
	parser->AddStrokes(group, true);
	end = time(NULL);

	std::stringstream ss;

	ss << "Parser Added strokes (" << (end - start) << "s)\n";

	if (USE_PARSER_CHECKED) {
		scg::ExpressionIterator* eit = parser->CreateDefaultIterator();
		const scg::ExpressionTree* _exprTree;
		bool isLF = false;
		while (eit != NULL && (_exprTree = eit->next()) != NULL) {

			if (_exprTree->type() == scg::PAREN_EXPR) {
				_exprTree = _exprTree->child(scg::PAREN_CONTENTS);
				if (_exprTree->type() != scg::MATRIX_EXPR) continue;
			} else if (_exprTree->type() != scg::MATRIX_EXPR) {
				continue;
			}

			// We now know this is a matrix
			std::string numS = _exprTree->child(scg::MATRIX_NUMROWS)->str();
			int numRows = atoi(numS.c_str());
			numS = _exprTree->child(scg::MATRIX_NUMCOLS)->str();
			int numCols = atoi(numS.c_str());
			for (int i = 0; i < numRows; i++)
			{
				scg::ExpressionTree* rowExpr =  (scg::ExpressionTree *) _exprTree->child(scg::MATRIX_ROWS)->child(i);
				if (rowExpr->type() != scg::MATRIXROW_EXPR)
				{
					THROW_ERROR(E_INVALID, "nullinatorium");
				}
				for (int j = 0; j < numCols; j++) {
					scg::ExpressionTree* colExpr = (scg::ExpressionTree *) rowExpr->child(j);
					if (colExpr->nchildren() == 0) {
						ss << "-\t";
					} else {
						ss << colExpr->long_str() << "\t";
					}
				}
				ss << std::endl;
			}

			if (!isLF && _exprTree->HasLongForm()) {
				ss << std::endl << std::endl;
				eit->release();
				eit = _exprTree->CreateLongFormIterator();
				isLF = true;
			} else {
				break;
			}
		}

		eit->release();
	} else {
		scg::MatrixAnalyzer* ma = scg::CreateMatrixAnalyzer(parser);
		const scg::context ctx = parser->get_context();
		std::vector<scg::segment*> segs = ctx.segments;

		start = time(NULL);
		for (std::vector<scg::segment*>::iterator it = segs.begin(); it != segs.end(); it++) {
			ma->addStroke(&(*it)->stk);
		}
		end = time(NULL);

		ss << "MA Added strokes (" << (end - start) << "s)" << std::endl;
		
		//for (std::vector<scg::segment*>::const_iterator it = segs.begin(); it != segs.end(); it++) {
		//	scg::RawStroke *stk = new scg::RawStroke;
		//	*stk = scg::copy((*it)->stk);
		//	ma->addStroke(stk);
		//}

		/*ma->addStroke(&group.strokes[2]);
		ma->addStroke(&group.strokes[3]);
		ma->addStroke(&group.strokes[4]);
		ma->addStroke(&group.strokes[5]);
		ma->addStroke(&group.strokes[6]);* /

//		ma->addStrokes(group);
		//ma->removeStroke(&group.strokes[2]);

		ss << "Short form:" << std::endl;
		start = time(NULL);
		scg::MatrixIterator* it = ma->createIterator();
		while (it->hasNext()) {
			scg::Matrix* m = it->next();
			m->rebuild();
			const char* output = m->str();
			ss << output << "\n";

			double conf;
			ma->getConfidence(m, conf);
			ss << "Confidence: " << conf << "\n";

			unsigned rows, cols;
			m->getDimensions(rows, cols);
			for (unsigned int i = 0; i < rows; i++) {
				for (unsigned int j = 0; j < cols; j++) {
					scg::MatrixElement* me =(*m)(i, j);
					if (me == NULL) continue;
					scg::Rect<long> r = me->getBoundingRect();

					RECT rc;
					rc.left = r.left;
					rc.top = r.top;
					rc.bottom = r.bottom;
					rc.right = r.right;
		            
					scg::TPC_inkspace_to_pixel(msink, rc.left, rc.top);
					scg::TPC_inkspace_to_pixel(msink, rc.right, rc.bottom);
		            
					FrameRect(GetDC(GetActiveWindow()), &rc, (HBRUSH)GetStockObject(GRAY_BRUSH));
					
					/*scg::ExpressionIterator* results;
					m->getCell(i, j, results);
					results->count();* /
				}
			}
		}
		end = time(NULL);
		ss << (end - start) << "s" << std::endl;

		delete it;

		if (!ma->hasValidLongForm()) {
			ss << "\n(not valid long form)";
		}

		ss << "\nLong form:" << std::endl;
		start = time(NULL);
		scg::LongFormMatrixIterator* lfmi = ma->createLongFormIterator();
		while (lfmi->hasNext()) {
			scg::UnderspecifiedMatrix* m = dynamic_cast<scg::UnderspecifiedMatrix*>(lfmi->next());
			m->rebuild();
			const char* output = m->str();
			ss << output << "\n";

			double conf;
			ma->getConfidence(m, conf);
			ss << "Confidence: " << conf << "\n";
		}
		end = time(NULL);
		ss << (end - start) << "s" << std::endl;

		delete lfmi;
		DestroyMatrixAnalyzer(ma);
	}

	std::string ssOutput = ss.str();
	const char* output = ssOutput.c_str();
	//const unsigned outputLength = ssOutput.length();
	//WCHAR* wideOutput = new WCHAR[outputLength + 1];
	//MultiByteToWideChar(0, 0, output, ssOutput.length(), wideOutput, outputLength + 1);
	//DrawTextW(GetDC(GetActiveWindow()), wideOutput, -1, &draw_dim, DT_NOCLIP);
	DrawTextA(GetDC(GetActiveWindow()), output, -1, &draw_dim, DT_NOCLIP);
	//delete [] wideOutput;

	std::ofstream ofs("C:\\Documents and Settings\\Robert Amlard\\Desktop\\out.txt");
	ofs << ssOutput;
	ofs.close();
}*/


static LRESULT CALLBACK
WndProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	enum {
		RECOGNIZER_DEFAULT,
		RECOGNIZER_MICROSOFT,
		RECOGNIZER_TABLET
	};
	
	static scg::TPC_InkWindow msink;
	
	static int recognizer_special = RECOGNIZER_TABLET;
	static int selected_box = -1;
	static int selected_box2 = -1;
	
	static int flags = 0;
		
	int prev_selected_box;
	int e;
	int x, y;
    char *label;
    RECT dim;
    	
    static LARGE_INTEGER time;
    static LARGE_INTEGER freq;

	switch (uMsg) {
	case WM_CREATE:
	  {
        QueryPerformanceFrequency(&freq);
        time.QuadPart = 0;
        
		e = scg::TPC_create_window(hwnd, msink);
		if (FAILURE(e)) {
		    show_error(scg::error_message(e));
		}
		    scg::TPC_enable_collector(msink);
						
		SetMenuRecognizer(hwnd, ID_RECOGNIZER_TABLET);

		//scg::SetTrainingPath("C:\\training");
		//scg::SetGrammar("diffeq");
		e = scg::InitializeRecognizer();
		if (FAILURE(e)) {
			std::stringstream ss;
			scg::error e = scg::get_error();
			ss << "Error " << e.code << ": ";
			if (e.msg.empty()) ss << scg::error_message(e.code) << ". No more information is available.";
			else ss << e.msg;
			std::wstring ws = scg::str2wstr(ss.str());
			MessageBox(0, ws.c_str(), L"Error", MB_OK|MB_ICONERROR);
			std::exit(0);
		}
		
		return 0;
	  }
/*
    case WM_PAINT:
      {
        //if (annotation_mode) {
            if (have_symbol_annotations) {
            PAINTSTRUCT ps;
            
            BeginPaint(hwnd, &ps);
            
            LOGBRUSH red_brush_info;
            red_brush_info.lbStyle = BS_SOLID;
            red_brush_info.lbColor = RGB(255, 0, 0);            
            HBRUSH red_brush = CreateBrushIndirect(&red_brush_info);
            LOGBRUSH blue_brush_info;
            blue_brush_info.lbStyle = BS_SOLID;
            blue_brush_info.lbColor = RGB(0, 0, 255);            
            HBRUSH blue_brush = CreateBrushIndirect(&blue_brush_info);
            LOGBRUSH green_brush_info;
            green_brush_info.lbStyle = BS_SOLID;
            green_brush_info.lbColor = RGB(0, 255, 0);
            HBRUSH green_brush = CreateBrushIndirect(&green_brush_info);
            
            for (unsigned i = 0; i < num_strokes(annotated_strokes); i++) {
                if (annotated_strokes.annotations[i].nstrokes) {
                    scg::Rect<long> bounds = scg::bbox(annotated_strokes.strokes() + i, annotated_strokes.annotations[i].nstrokes);
                    
                    RECT rc;
                    rc.left = bounds.left;
                    rc.top = bounds.top;
                    rc.right = bounds.right;
                    rc.bottom = bounds.bottom;
            
                    scg::TPC_inkspace_to_pixel(msink, rc.left, rc.top);
                    scg::TPC_inkspace_to_pixel(msink, rc.right, rc.bottom);
            
                    if (i == selected_box || i == selected_box2) {
                        FrameRect(ps.hdc, &rc, red_brush);                    
                    }
                    else {
                        FrameRect(ps.hdc, &rc, blue_brush);
                    }
                }
            }
           
            if (have_link_annotations) {
                SelectObject(ps.hdc, green_brush);
                for (std::vector<LinkAnnotation>::const_iterator i = link_annotations.begin(); i != link_annotations.end(); ++i) {
                    scg::Rect<long> bounds1 = scg::bbox(annotated_strokes.strokes() + i->box1, annotated_strokes.annotations[i->box1].nstrokes);
                    scg::Rect<long> bounds2 = scg::bbox(annotated_strokes.strokes() + i->box2, annotated_strokes.annotations[i->box2].nstrokes);

                    scg::TPC_inkspace_to_pixel(msink, bounds1.left, bounds1.top);
                    scg::TPC_inkspace_to_pixel(msink, bounds1.right, bounds1.bottom);
                    scg::TPC_inkspace_to_pixel(msink, bounds2.left, bounds2.top);
                    scg::TPC_inkspace_to_pixel(msink, bounds2.right, bounds2.bottom);
                    
                    MoveToEx(ps.hdc, (bounds1.left + bounds1.right) / 2, (bounds1.top + bounds1.bottom) / 2, NULL);
                    LineTo(ps.hdc, (bounds2.left + bounds2.right) / 2, (bounds2.top + bounds2.bottom) / 2);

                    switch (i->link) {
                    case Left:
                        label = "left";
                        break;
                    case Above:
                        label = "above";
                        break;
                    case BelowLeft:
                        label = "below-left";
                        break;
                    case AboveLeft:
                        label = "above-left";
                        break;
                    case Inside:
                        label = "inside";
                        break;
                    }
                    
        	        dim.left = (((bounds1.left + bounds1.right) / 2) + ((bounds2.left + bounds2.right) / 2)) / 2;
        	        dim.bottom = (((bounds1.top + bounds1.bottom) / 2) + ((bounds2.top + bounds2.bottom) / 2)) / 2;
        	        dim.right = dim.left + 200;
        	        dim.top = dim.bottom - 20;
	                DrawText(ps.hdc, label, strlen(label), &dim, DT_CALCRECT);
	                DrawText(ps.hdc, label, strlen(label), &dim, 0);
                }
            }
            
            EndPaint(hwnd, &ps);
        //}
        }
        return 0;
      }*/
        
    case WM_LBUTTONUP:
      {
        x = LOWORD(lParam);
        y = HIWORD(lParam);
        /*if (annotation_mode) {
            prev_selected_box = selected_box;
            for (unsigned i = 0; i < num_strokes(annotated_strokes); i++) {
                if (annotated_strokes.annotations[i].nstrokes) {
                    scg::Rect<long> bounds = scg::bbox(annotated_strokes.strokes + i, annotated_strokes.strokes + i + annotated_strokes.annotations[i].nstrokes);
                    
                    scg::TPC_inkspace_to_pixel(msink, bounds.left, bounds.top);
                    scg::TPC_inkspace_to_pixel(msink, bounds.right, bounds.bottom);

                    if (bounds.left <= x && x <= bounds.right
                     && bounds.top <= y && y <= bounds.bottom) {
                        if (selected_box != -1) {
                            scg::Rect<long> sel_bounds = scg::bbox(annotated_strokes.strokes + selected_box, annotated_strokes.strokes + selected_box + annotated_strokes.annotations[selected_box].nstrokes);
                            scg::TPC_inkspace_to_pixel(msink, sel_bounds.left, sel_bounds.top);
                            scg::TPC_inkspace_to_pixel(msink, sel_bounds.right, sel_bounds.bottom);
                            if (sel_bounds.left <= x && x <= sel_bounds.right
                             && sel_bounds.top <= y && y <= sel_bounds.bottom) {
                                if (sel_bounds.left <= bounds.left && bounds.left <= sel_bounds.right
                                 || sel_bounds.left <= bounds.right && bounds.right <= sel_bounds.right
                                 || sel_bounds.top <= bounds.top && bounds.top <= sel_bounds.bottom
                                 || sel_bounds.top <= bounds.bottom && bounds.bottom <= sel_bounds.bottom) {
                                    selected_box = i;
                                    InvalidateRect(hwnd, NULL, TRUE);
                                }
                            }
                            else {
                                selected_box = i;
                                InvalidateRect(hwnd, NULL, TRUE);
                            }
                        }
                        else {
                            selected_box = i;
                            InvalidateRect(hwnd, NULL, TRUE);
                        }
                    }
                }
            }
            
            if (prev_selected_box >= 0) {
                selected_box2 = selected_box;
                selected_box = prev_selected_box;
            }
        }*/
        return 0;
      }
      
	case WM_COMMAND:
		if (wParam == ID_FILE_CLEAR) {
            e = scg::TPC_clear_ink(msink);
            if (FAILURE(e)) {
                show_error(scg::error_message(e));
            }
            InvalidateRect(hwnd, 0, TRUE);
            have_symbol_annotations = false;
            have_link_annotations = false;
            selected_box = selected_box2 = -1;
            if (rec) rec->Clear();
		}

		else if (wParam == ID_FILE_LOADINK || wParam == ID_FILE_LOADIPADINK) {
            selected_box = selected_box2 = -1;
            e = scg::TPC_clear_ink(msink);
            if (FAILURE(e)) {
                show_error(scg::error_message(e));
            }
            else {
				wchar_t filename[1024] = {0};
				wchar_t filetitle[256] = {0};
			
				OPENFILENAME ofn;
				ofn.lStructSize = sizeof(ofn);
				ofn.hwndOwner = hwnd;
				ofn.lpstrFilter = L"MSInk files\0*.msink\0Ink files\0*.ink\0All files\0*.*\0\0";
				ofn.lpstrCustomFilter = NULL;
				ofn.nFilterIndex = 1;
				ofn.lpstrFile = filename;
				ofn.nMaxFile = sizeof(filename);
				ofn.lpstrFileTitle = filetitle;
				ofn.nMaxFileTitle = sizeof(filetitle);
				ofn.lpstrInitialDir = NULL;
				ofn.lpstrTitle = L"Load Ink";
				ofn.Flags = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;
				ofn.lpstrDefExt = L"ink";
				
				if (GetOpenFileName(&ofn)) {
				    char *buf;
				    size_t n = 1+wcstombs(0, ofn.lpstrFile, 0);
				    buf = new char[n];
				    wcstombs(buf, ofn.lpstrFile, n);
					if (wParam == ID_FILE_LOADINK) {
						e = LoadInk(buf, msink);
					}
					else {
						e = LoadIpadInk(buf, msink);
					}
					delete[] buf;
					if (FAILURE(e)) {
					    show_error(scg::error_message(e));
					}
					else {
					    InvalidateRect(hwnd, 0, TRUE);
					}
				}		
			}
		}

		else if (wParam == ID_FILE_SAVEINK) {
			wchar_t filename[1024] = {0};
			wchar_t filetitle[256] = {0};
		
			OPENFILENAME ofn;
			ofn.lStructSize = sizeof(ofn);
			ofn.hwndOwner = hwnd;
			ofn.lpstrFilter = L"Ink files\0*.ink\0All files\0*.*\0\0";
			ofn.lpstrCustomFilter = NULL;
			ofn.nFilterIndex = 1;
			ofn.lpstrFile = filename;
			ofn.nMaxFile = sizeof(filename);
			ofn.lpstrFileTitle = filetitle;
			ofn.nMaxFileTitle = sizeof(filetitle);
			ofn.lpstrInitialDir = NULL;
			ofn.lpstrTitle = L"Save Ink";
			ofn.Flags = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;
			ofn.lpstrDefExt = L"ink";
			
			if (GetSaveFileName(&ofn)) {
			    char *buf;
			    size_t n = 1+wcstombs(0, ofn.lpstrFile, 0);
			    buf = new char[n];
			    wcstombs(buf, ofn.lpstrFile, n);
				e = SaveInk(buf, msink);
				delete[] buf;
				if (FAILURE(e)) {
				    show_error(scg::error_message(e));
				}
			}
		}

		/*else if (wParam == ID_TEST_POLYNOM) {
			std::ofstream ofs("C:\\Documents and Settings\\Robert Amlard\\Desktop\\out.txt");
			int numPassed = 0;
			int totalNum = 0;

			std::vector<std::string> answers, results;

			// 1 + 1
			answers.push_back("2");
			scg::BoundedPolynomialFunction f1, g1;
			scg::PolyTerm t1(1,0);
			f1.addTerm(t1);
			g1.addTerm(t1);
			std::string poly = (f1+g1).str();
			results.push_back(poly);

			// x^2-1 + x^2+1
			answers.push_back("2x^2");
			scg::BoundedPolynomialFunction f2, g2;
			scg::PolyTerm pos1(1,0);
			scg::PolyTerm neg1(-1,0);
			scg::PolyTerm x_sq(1,2);
			f2.addTerm(pos1);
			f2.addTerm(x_sq);
			g2.addTerm(neg1);
			g2.addTerm(x_sq);
			poly = (f2+g2).str();
			results.push_back(poly);

			// (x^2-1) * (x^2+1)
			answers.push_back("x^4-1");
			poly = (f2*g2).str();
			results.push_back(poly);

			// 3*((x^2-1) * (x^2+1))
			answers.push_back("3x^4-3");
			poly = (3*(f2*g2)).str();
			results.push_back(poly);

			// ((x^2-1) * (x^2+1)) / 2
			answers.push_back("0.5x^4-0.5");
			poly = ((f2*g2) / 2).str();
			results.push_back(poly);

			// 1 * (x-2)/-1*(x-3)/-2 + 4*(x-1)/1*(x-3)/-1 + 9*(x-1)/2*(x-2)/1
			answers.push_back("x^2");
			scg::BoundedPolynomialFunction f3[2], g3[2], h3[2];
			scg::PolyTerm neg2(-2,0);
			scg::PolyTerm neg3(-3,0);
			scg::PolyTerm x(1,1);
			f3[0].addTerm(x);
			f3[0].addTerm(neg2);
			f3[1].addTerm(x);
			f3[1].addTerm(neg3);
			g3[0].addTerm(x);
			g3[0].addTerm(neg1);
			g3[1].addTerm(x);
			g3[1].addTerm(neg3);
			h3[0].addTerm(x);
			h3[0].addTerm(neg1);
			h3[1].addTerm(x);
			h3[1].addTerm(neg2);
			poly = (1 * (f3[0]/-1) * (f3[1]/-2) + 4 * (g3[0]/1) * (g3[1]/-1) + 9 * (h3[0]/2) * (h3[1]/1)).str();
			results.push_back(poly);

			// 1 * (x-2)/-1*(x-3)/-2 + 8*(x-1)/1*(x-3)/-1 + 27*(x-1)/2*(x-2)/1
			answers.push_back("6x^2-11x+6");
			poly = (1 * (f3[0]/-1) * (f3[1]/-2) + 8 * (g3[0]/1) * (g3[1]/-1) + 27 * (h3[0]/2) * (h3[1]/1)).str();
			results.push_back(poly);

			// Compare results to answers
			for (unsigned i = 0; i < answers.size(); i++) {
				bool passed = (results[i] == answers[i]);
				std::string passOrFail = passed ? "PASS" : "FAIL";
				ofs << answers[i] << ": " << results[i] << " --" << passOrFail << std::endl;
				if (passed) numPassed++;
			}
			totalNum += answers.size();

			answers.clear();
			results.clear();
			ofs << "Lagrange Polynomials:" << std::endl;

			// (1,1), (2,4), (3,9)
			answers.push_back("x^2");
			std::vector<std::pair<int,int> > points;
			points.push_back(std::make_pair(1,1));
			points.push_back(std::make_pair(2,4));
			points.push_back(std::make_pair(3,9));
			scg::BoundedPolynomialFunction f = scg::createLagrangePolynomial(points);
			poly = f.str();
			results.push_back(poly);

			// (1,1), (2,8), (3,27)
			answers.push_back("6x^2-11x+6");
			points.clear();
			points.push_back(std::make_pair(1,1));
			points.push_back(std::make_pair(2,8));
			points.push_back(std::make_pair(3,27));
			f = scg::createLagrangePolynomial(points);
			poly = f.str();
			results.push_back(poly);

			// (0,0), (1,1), (2,8), (3,27)
			answers.push_back("x^3");
			points.clear();
			points.push_back(std::make_pair(0,0));
			points.push_back(std::make_pair(1,1));
			points.push_back(std::make_pair(2,8));
			points.push_back(std::make_pair(3,27));
			f = scg::createLagrangePolynomial(points);
			poly = f.str();
			results.push_back(poly);

			// Compare results to answers
			for (unsigned i = 0; i < answers.size(); i++) {
				bool passed = (results[i] == answers[i]);
				std::string passOrFail = passed ? "PASS" : "FAIL";
				ofs << answers[i] << ": " << results[i] << " --" << passOrFail << std::endl;
				if (passed) numPassed++;
			}
			totalNum += answers.size();

			ofs.close();

			std::stringstream ss;
			ss << numPassed << "/" << totalNum << " passed";
			std::string stat_str = ss.str();
			const unsigned statLength = stat_str.length();
			const char* stat = stat_str.c_str();
			WCHAR* wideOutput = new WCHAR[statLength + 1];
			MultiByteToWideChar(0, 0, stat, statLength, wideOutput, statLength + 1);
			HDC dc = GetDC(GetActiveWindow());
			RECT dim = {10, 10, 200, 200};
			DrawTextW(dc, wideOutput, -1, &dim, DT_NOCLIP);
			delete [] wideOutput;
		}

		else if (wParam == ID_RANGE) {
			// To test ranges, separate elements in the range by the character "alpha"
			// For example, to test the range "1 4 9 ... 81" (i.e. a range that is defined by x^2), write:
			// 1 alpha 2 alpha 4 alpha 81
			// The last element will always be treated as the bounding element in the range.

			scg::SetVerbosity(2);
			scg::VerboseOutputToFile("C:\\Documents and Settings\\Robert Amlard\\Desktop\\log.txt");

			scg::TPC_StrokeGroup TPCgroup;
			e = scg::TPC_get_strokes(msink, TPCgroup);
			if (FAILURE(e)) {
				return e;
			}
	    
			scg::RawStrokeGroup group;
			e = scg::convert(group, TPCgroup);
			if (FAILURE(e)) {
				return e;
			}
			
			scg::MathRecognizer* parser = scg::CreateMathRecognizer();
			scg::RangeExpander re(parser);

			parser->AddStrokes(group.strokes, group.nstrokes, true);

			scg::ExpressionIterator* eit;
			unsigned curStroke[1];
			unsigned* indices;
			unsigned curElemStrokeNum = 0;
			for (unsigned i = 0; i < group.nstrokes; i++) {
				curStroke[0] = i;
				eit = parser->CreateIteratorForStrokesByIndex(curStroke, 1);
				const scg::ExpressionTree* tree = eit->next();

				// separate elements by alpha
				if (tree != NULL && strcmp(tree->child(0)->str(), "alpha") == 0) {
					indices = new unsigned[curElemStrokeNum];

					for (unsigned j = 0; j < curElemStrokeNum; j++) {
						indices[j] = (i - 1) - j;
					}
					scg::ExpressionIterator* elemIter = parser->CreateIteratorForStrokesByIndex(indices, curElemStrokeNum);
					re.appendStartExpression(elemIter->next());
					elemIter->release();
					
					delete [] indices;
					curElemStrokeNum = 0;
				} else {
					curElemStrokeNum++;
				}

				eit->release();
				if (tree != NULL) tree->release();
			}
			
			// last element
			indices = new unsigned[curElemStrokeNum];
			for (unsigned j = 0; j < curElemStrokeNum; j++) {
				indices[j] = (group.nstrokes - 1) - j;
			}
			eit = parser->CreateIteratorForStrokesByIndex(indices, curElemStrokeNum);
			re.setEndExpression(eit->next());
			eit->release();
			delete [] indices;

			re.updateTemplateTree();

			unsigned i;
			std::ofstream ofs("C:\\Documents and Settings\\Robert Amlard\\Desktop\\out.txt");
			for (i = 0; i < re.size(); i++) {
				scg::ExpressionTree* tree = re.at(i);
				if (tree != NULL) {
					ofs << tree->long_str() << "\n";
				}
			}
			re.at(i);
			ofs.close();
		}

		else if (wParam == ID_MATRIX) {
			DoMatrixRecognition(msink);
		}*/

		else if (wParam == ID_USE_PARSER) {

			USE_PARSER_CHECKED = !USE_PARSER_CHECKED;

			HMENU menu = GetMenu(hwnd);
			menu = GetSubMenu(menu, 2);
			if (USE_PARSER_CHECKED) {
				CheckMenuItem(menu, ID_USE_PARSER, MF_BYCOMMAND | MF_CHECKED);
			} else {
				CheckMenuItem(menu, ID_USE_PARSER, MF_BYCOMMAND | MF_UNCHECKED);
			}
			DrawMenuBar(hwnd);
		}

		else if (wParam == ID_RECOGNIZE_RECOGNIZE) {
			switch (recognizer_special) {
			case RECOGNIZER_MICROSOFT:
				RunMicrosoftRecognizer(msink);
				break;
			case RECOGNIZER_TABLET:
			    RunTabletRecognizer(msink);
			    break;
			default:
				Recognize(msink);
				break;
			}
		}

		else if (wParam == ID_RECOGNIZE_DEFORMABLETEMPLATES) {
			recognizer_special = RECOGNIZER_DEFAULT;
			SetMenuRecognizer(hwnd, ID_RECOGNIZE_DEFORMABLETEMPLATES);
		}
		else if (wParam == ID_RECOGNIZE_ELASTICMATCHER) {
			recognizer_special = RECOGNIZER_DEFAULT;
			SetMenuRecognizer(hwnd, ID_RECOGNIZE_ELASTICMATCHER);
		}

		else if (wParam == ID_RECOGNIZE_STRUCTURALMATCHER) {
			recognizer_special = RECOGNIZER_DEFAULT;
			SetMenuRecognizer(hwnd, ID_RECOGNIZE_STRUCTURALMATCHER);
		}

		else if (wParam == ID_RECOGNIZE_STRUCTURALOPTIMIZER) {
			recognizer_special = RECOGNIZER_DEFAULT;
			SetMenuRecognizer(hwnd, ID_RECOGNIZE_STRUCTURALOPTIMIZER);
		}

		else if (wParam == ID_RECOGNIZER_FOURIERANALYZER) {
		   	recognizer_special = RECOGNIZER_DEFAULT;
			SetMenuRecognizer(hwnd, ID_RECOGNIZER_FOURIERANALYZER);
		}

		else if (wParam == ID_RECOGNIZER_CITRECOGNIZER) {
			recognizer_special = RECOGNIZER_DEFAULT;
			SetMenuRecognizer(hwnd, ID_RECOGNIZER_CITRECOGNIZER);
		}

		else if (wParam == ID_RECOGNIZER_MICROSOFT) {
			recognizer_special = RECOGNIZER_MICROSOFT;
			SetMenuRecognizer(hwnd, ID_RECOGNIZER_MICROSOFT);
		}
		
		else if (wParam == ID_RECOGNIZER_TABLET) {
		    recognizer_special = RECOGNIZER_TABLET;
		    SetMenuRecognizer(hwnd, ID_RECOGNIZER_TABLET);
		}

		else if (wParam == ID_RECOGNIZE_RESET) {
            e = scg::TPC_clear_ink(msink);
            if (FAILURE(e)) {
                show_error(scg::error_message(e));
            }
		}
		
		else if (wParam == ID_RECOGNIZER_SUBDIVIDE) {
			scg::SetTabletResolution(132);
			//RunFancySubdivision(hwnd, msink);
		}
		
		else if (wParam == ID_RECOGNIZER_ALLOWSUBTREEALTERNATES) {
			//flags = flags ? 0 : scg::ExpressionIterator::Flags::AllowSubtreeAlternates;
		}
		else if (wParam == ID_RECOGNIZER_SEMANTICS) {
			scg::NoVerboseOutput();
			scg::VerboseOutputToFile("expr.txt");
			scg::SetVerbosity(1);
			if (!rec) {
				//scg::ShutdownRecognizer(rh);
				//rh = scg::InitializeRecognizer();
				rec = scg::CreateMathRecognizer(scg::TPC_get_ink(msink));
				//std::ifstream inkin("C:\\pen-math\\bad-files\\Example8reco");
				//rec = scg::CreateMathRecognizer(inkin);
		    //rec->BuildParseTable();
			std::ofstream testout("savefile", std::ofstream::binary);
			/*rec->Save(testout);
			rec->release();
			testout.close();
			std::ifstream testin("savefile", std::ifstream::binary);
			rec = scg::CreateMathRecognizer(testin);*/
			//rec->BuildParseTable();
			}
			else {
				rec->release();
				rec = scg::CreateMathRecognizer(scg::TPC_get_ink(msink));
		    //rec->BuildParseTable();
			}
			
		    if (!rec) {
		        abort();
		    }
/*
		   scg::TPC_clear_ink(msink);

			rec->Translate(-1000,-1000);
			
			IInkDisp *ink;		   
			rec->GetInk(&ink);
		   scg::TPC_disable_collector(msink);
		   scg::TPC_set_ink(msink, ink);
		   scg::TPC_enable_collector(msink);
*/    	
		   std::ofstream ofs("exprout.txt");
    		std::wofstream wofs("wexprout.txt");
			const scg::ExpressionTree *top = rec->GetTopExpression();
			/*if (top) {
				ofs << "top result:\n";
				write_node(ofs, top);
				ofs << "type was " << top->type() << std::endl;
				int e = casclient_connect();
				if (e != 0) {
					ofs << "FAILED TO CONNECT\n";
				}
				else {
					casclient_disconnect();
				}
			}*/
    		scg::ExpressionIterator *it1 = rec->CreateDefaultIterator();
    		/*unsigned n = 12;
    		scg::ExpressionTree *ph = scg::CreateBlankExpression();
			scg::ExpressionIterator *it2 = rec->CreateIteratorForExpression(ph);
			ofs << "FAKE ITERATOR GIVES:\n";
	  		for (const scg::ExpressionTree *expr2 = it2->next(); expr2; expr2 = it2->next()) {
				write_node(ofs, expr2);
	  			expr2->release();
			}
			it2->release();
			ph->release();
			ofs << "DONE FAKE ITERATION\n";*/
			if (it1) {
			for (const scg::ExpressionTree *expr = it1->next(/*scg::ExpressionIterator::MATHML_WRAPPER*/); expr; expr = it1->next(/*scg::ExpressionIterator::MATHML_WRAPPER*/)) {
				write_node(ofs, expr);
				write_wide_node(wofs, expr);
				/*scg::ExpressionIterator *it2 = rec->CreateIteratorForExpression(expr);
				ofs << "FAKE ITERATOR GIVES:\n";
	  			for (const scg::ExpressionTree *expr2 = it2->next(); expr2; expr2 = it2->next()) {
					write_node(ofs, expr2);
	  				expr2->release();
				}
				it2->release();*/
				/*if (expr->HasLongForm()) {
					ofs << "LONG FORM DETECTED; iterating\n";
					scg::ExpressionIterator *liter = expr->CreateLongFormIterator();
					for (const scg::ExpressionTree *lexpr = liter->next(scg::ExpressionIterator::MATHML_WRAPPER); lexpr; lexpr = liter->next(scg::ExpressionIterator::MATHML_WRAPPER)) {
						write_node(ofs, lexpr);
						const scg::ExpressionTree *c = lexpr->child(1);
						if (c->type() == scg::MATRIX_EXPR) {
							ofs << " MATRIX! There are " << c->child(scg::MATRIX_NUMROWS)->str() << " rows and " << c->child(scg::MATRIX_NUMCOLS)->str() << " cols\n";
						}
						lexpr->release();
					}
					liter->release();
					ofs << "DONE LONG-FORM ITERATION\n";
				}*/
				expr->release();
			}
			it1->release();
			}
			//delete it1;
			//rec->release();
		}

		else if (wParam == ID_RECOGNIZER_MULTILINE) {
			scg::NoVerboseOutput();
			scg::VerboseOutputToFile("expr.txt");
			scg::SetVerbosity(1);
			if (rec) {
				rec->release();
			}
			rec = scg::CreateMathRecognizer(scg::TPC_get_ink(msink));
		    if (!rec) {
		        abort();
		    }
			std::ofstream ofs("exprout.txt");
			scg::ExpressionIterator *it = rec->CreateDefaultSemanticIterator(scg::MULTI_EXPR);
			for (const scg::ExpressionTree *expr = it->next(); expr; expr = it->next()) {
				write_node(ofs, expr);
				expr->release();
			}
			delete it;
		}
		
		else if (wParam == ID_MODE_DRAW) {
		    have_link_annotations = false;
		    have_symbol_annotations = false;
		    scg::TPC_enable_collector(msink);
		    annotation_mode = false;
		}
		
		else if (wParam == ID_MODE_ANNOTATE) {
		    if (have_symbol_annotations) {
		        //link_annotations.clear();
    		    scg::TPC_disable_collector(msink);
		        annotation_mode = true;
		        have_link_annotations = true;
		    }
		}
		
		else if (wParam == ID_ANNOTATION_LEFT) {
		    if (selected_box >= 0 && selected_box2 >= 0) {
		        LinkAnnotation link;
		        link.box1 = selected_box;
		        link.box2 = selected_box2;
		        link.link = Left;
		        link_annotations.push_back(link);
		        selected_box = selected_box2 = -1;
		        InvalidateRect(hwnd, NULL, TRUE);		        
		    }
		}
		
		else if (wParam == ID_ANNOTATION_ABOVE) {
		    if (selected_box >= 0 && selected_box2 >= 0) {
		        LinkAnnotation link;
		        link.box1 = selected_box;
		        link.box2 = selected_box2;
		        link.link = Above;
		        link_annotations.push_back(link);
		        selected_box = selected_box2 = -1;
		        InvalidateRect(hwnd, NULL, TRUE);
		    }
		}
		
		else if (wParam == ID_ANNOTATION_BELOWLEFT) {
		    if (selected_box >= 0 && selected_box2 >= 0) {
		        LinkAnnotation link;
		        link.box1 = selected_box;
		        link.box2 = selected_box2;
		        link.link = BelowLeft;
		        link_annotations.push_back(link);
		        selected_box = selected_box2 = -1;
		        InvalidateRect(hwnd, NULL, TRUE);
		    }
		}
		
		else if (wParam == ID_ANNOTATION_BELOWRIGHT) {
		    if (selected_box >= 0 && selected_box2 >= 0) {
		        LinkAnnotation link;
		        link.box1 = selected_box;
		        link.box2 = selected_box2;
		        link.link = AboveLeft;
		        link_annotations.push_back(link);
		        selected_box = selected_box2 = -1;
		        InvalidateRect(hwnd, NULL, TRUE);
		    }
		}
		else if (wParam == ID_ANNOTATION_INSIDE) {
		    if (selected_box >= 0 && selected_box2 >= 0) {
		        LinkAnnotation link;
		        link.box1 = selected_box;
		        link.box2 = selected_box2;
		        link.link = Inside;
		        link_annotations.push_back(link);
		        selected_box = selected_box2 = -1;
		        InvalidateRect(hwnd, NULL, TRUE);
		    }
		}
		
		else if (wParam == ID_FILE_EXIT) {
			DestroyWindow(hwnd);
		}

		return 0;

	case WM_DESTROY:
		scg::ShutdownRecognizer();
        scg::TPC_release_window(msink);
		if (tablet_reco) {
		    tablet_reco->Release();
		    tablet_reco = 0;
		}
		PostQuitMessage(0);
		return 0;
	}

	return DefWindowProc(hwnd, uMsg, wParam, lParam);
}


static bool
RegisterWindowClass(HINSTANCE hInstance, WNDPROC proc, int menu, const wchar_t *window_class)
{
    WNDCLASSEX WndClassEx;

    WndClassEx.cbSize        = sizeof(WndClassEx);
    WndClassEx.style         = CS_HREDRAW | CS_VREDRAW;
    WndClassEx.lpfnWndProc   = proc;
    WndClassEx.cbClsExtra    = 0;
    WndClassEx.cbWndExtra    = 0;
    WndClassEx.hInstance     = hInstance;
    WndClassEx.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    WndClassEx.hIconSm       = WndClassEx.hIcon;
    WndClassEx.hCursor       = LoadCursor(NULL, IDC_ARROW);
    WndClassEx.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
    WndClassEx.lpszMenuName  = MAKEINTRESOURCE(menu);
    WndClassEx.lpszClassName = window_class;

    if (!RegisterClassEx(&WndClassEx)) {
        show_error("Failed to register window class!");
        return false;
    }

    return true;
}


int APIENTRY
WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	if (!RegisterWindowClass(hInstance, WndProc, IDR_MENU, L"TestPad")) {
		show_error("Cannot register window class");
		abort();
	}

	HRESULT hr = CoInitializeEx(NULL, COINIT_MULTITHREADED);
	if (FAILED(hr) && hr != S_FALSE) { // S_FALSE ==> COM already initialized
		show_error("CoInitialize failed");
		abort();
	}

    HWND hWnd = CreateWindowEx(WS_EX_CLIENTEDGE, L"TestPad", L"TestPad",
                               WS_OVERLAPPEDWINDOW,
                               CW_USEDEFAULT, CW_USEDEFAULT, 
                               CW_USEDEFAULT, CW_USEDEFAULT,
                               NULL, NULL, hInstance, NULL);
    if (hWnd == NULL)
    {
		show_error("Failed to create main window");
		abort();		
    }

    ShowWindow(hWnd, nCmdShow);
    UpdateWindow(hWnd);
/*
	char *p;
	unsigned n;
	scg::GetUserProfile(0, &n);
	p = new char[n];
	scg::GetUserProfile(p, &n);
	
	for (scg::UnicodeIterator i = scg::FirstUnicodeSymbol(); i.unicode_value; i = scg::NextUnicodeSymbol(i)) {
		scg::unicode_char a = i.unicode_value;
	}
*/
/*
	int nprofiles = scg::CountAvailableProfiles();
	std::vector<std::string> paths;
	paths.insert(paths.end(), nprofiles, std::string());
	scg::GetAvailableProfiles(&paths[0], paths.size());
	
	for (std::vector<std::string>::const_iterator i = paths.begin(); i != paths.end(); ++i) {
		std::string test = *i;
		int a = 0;
	}*/
	
	 MSG msg; 
    while (GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    CoUninitialize();
    
    return (int)msg.wParam;
}