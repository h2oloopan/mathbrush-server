/*
 * tabletpc.h
 *
 * This file defines the interfaces that wrap some TabletPC SDK functionality.
 *
 * It is not used for recognition per se, but rather as a shared tool for test
 * programs, etc., which must interact with the tablet PC.
 */
 
#ifndef TABLET_H_
#define TABLET_H_


#include "tpc-group.h"

#define NOMINMAX
#include <windows.h>
#include <msinkaut.h>


namespace scg
{


struct TPC_InkWindow
{
    IInkCollector *collector;
    IInkRenderer *renderer;
    IInkOverlay *overlay;
};


int TPC_create_window(HWND hwnd, TPC_InkWindow &ink);
int TPC_release_window(TPC_InkWindow &ink);


int TPC_enable_collector(TPC_InkWindow &ink);
int TPC_disable_collector(TPC_InkWindow &ink);


int TPC_set_ink(TPC_InkWindow &inkw, IInkDisp *ink);

IInkDisp *TPC_get_ink(TPC_InkWindow &ink);
IInkRenderer *TPC_get_renderer(TPC_InkWindow &ink);
IInkOverlay *TPC_get_overlay(TPC_InkWindow &ink);


int TPC_pixel_to_inkspace(TPC_InkWindow &ink, long &x, long &y);
int TPC_inkspace_to_pixel(TPC_InkWindow &ink, long &x, long &y);


int TPC_clear_ink(TPC_InkWindow &ink);

int TPC_rescale_window(TPC_InkWindow &ink);


int TPC_get_strokes(TPC_InkWindow &ink, TPC_StrokeGroup &strokes);

IInkStrokes *TPC_get_APIstrokes(TPC_InkWindow &ink);

int TPC_add_strokes(TPC_InkWindow &ink, TPC_StrokeGroup &strokes);


}


#endif

