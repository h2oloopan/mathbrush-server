#ifndef INK_IO_H_
#define INK_IO_H_

#include <string>
#include <istream>
#include <ostream>

#include "error.h"
#include "group.h"
#include "memory.h"
#include "stroke.h"


namespace scg
{


template <typename T>
std::istream &
operator>>(std::istream &is, Stroke<T> &stroke)
{
    unsigned npoints;
    is >> npoints;
	 if (!is || is.eof()) {
	 	throw E_IO;
	 }
    
    T *x = new T[npoints];
    if (!x) {
        throw E_OUTOFMEM;
    }
    
    T *y = new T[npoints];
    if (!y) {
				delete[] x;
        throw E_OUTOFMEM;
    }
    
	if (npoints == 0) {
		stroke.set_points(x, y, 0);
		return is;
	}

		std::string ind;
		is >> ind;

    T *px = x;
    T *py = y;
    unsigned long *pt = 0;

		unsigned long *t = 0;
    
		if (ind == "xyt") {
			t = new unsigned long[npoints];
			if (!t) {
				delete[] x;
				delete[] y;
				throw E_OUTOFMEM;
			}
			pt = t;
		}
		else {
			std::stringstream ss;
			ss << ind;
			ss >> *px;
			is >> *py;
			++px;
			++py;
		}

    while (px != x + npoints) {
        is >> *px >> *py;
				if (pt) {
					is >> *pt;
				}

        /*
         Test for duplicate points; the tablet seems to generate these fairly often,
           possibly due to writer hesitation, etc.  However, we do not want duplicates,
           since they cause arclength between points to go to zero.  (arc length is commonly
           used as a denominator, for example when subdividing strokes.)
        */        
        if (px > x && (*px == *(px - 1) && *py == *(py - 1))) {
            npoints--;
        }
        else {
            px++;
            py++;
						if (pt) {
							pt++;
						}
        }
    }
    
		if (pt) {
    	stroke.set_points(x, y, t, npoints);
		}
		else {
    	stroke.set_points(x, y, npoints);
		}
    
    return is;
}


template <typename T>
std::ostream &
operator<<(std::ostream &os, const Stroke<T> &stroke)
{
    os << stroke.npoints << std::endl;
	if (stroke.time) {
		os << "xyt\n";
	}
    
    const T *px = stroke.x;
    const T *py = stroke.y;
    const unsigned long *pt = stroke.time;
	for (size_t i = 0; i < stroke.npoints; ++i) {
		os << stroke.x[i] << ' ' << stroke.y[i];
		if (stroke.time) {
			os << ' ' << stroke.time[i];
		}
		os << std::endl;
	}
    
    return os;
}


//DLLDECL extern const std::string StrokeGroupFileHeader;


template <typename T>
std::istream &
operator>>(std::istream &is, StrokeGroup<T> &group)
{
	const static std::string StrokeGroupFileHeader("SCG_INK");
    std::ws(is);
    for (size_t i = 0; i < StrokeGroupFileHeader.length(); i++) {
        if (is.peek() == StrokeGroupFileHeader[i]) {
            is.get();
        }
        else {
            if (i == 0) {
                break;
            }
            else {
                throw E_INVALID;
            }
        }
    }
    
    unsigned nstrokes;
    is >> nstrokes;
	 if (!is || is.eof()) {
	 	throw E_IO;
	 }
    
    Stroke<T> *strokes = new Stroke<T>[nstrokes];
    if (!strokes) {
        throw E_OUTOFMEM;
    }
    
    Stroke<T> *stroke = strokes;
    
    while (stroke != strokes + nstrokes) {
        is >> *(stroke++);
    }
    
    group.set_strokes(strokes, nstrokes);
    
    return is;
}


template <typename T>
std::ostream &
operator<<(std::ostream &os, const StrokeGroup<T> &group)
{
    //os << scg::StrokeGroupFileHeader << std::endl;
	os << "SCG_INK" << std::endl;
    os << group.nstrokes << std::endl;
    
	for (size_t i = 0; i < group.nstrokes; ++i) {
        os << group.strokes[i] << std::endl;
    }
    
    return os;
}


}

#endif

