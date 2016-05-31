//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Utility routines for visualizing and dumping out edge/curve/line/ellipse detection results
//
//  History:
//      2012/02/05-szeliski
//          Split off from main EdgeDetectTest.cpp file, for easier re-use
//
//------------------------------------------------------------------------

#include "vt_edgedetect.h"
#include "vt_linedetector.h"
#include "vt_fitEllipse.h"

namespace vt {

#pragma once

// Visualization Code that dumps out text files and images for edgels, curves, lines, elipses and bitangents
HRESULT VtDumpEdgels(CRGBImg& dst, const CRGBImg& src, float attenuation, const vector<EdgeSegment> &el, bool useZeroCrossings);
HRESULT VtDumpCurves(CRGBImg& dst, const CRGBImg& src, float attenuation, const vector< vector<int> > &curves, const vector<EdgeSegment> &el, bool useZeroCrossings);
HRESULT VtDumpLines(CRGBImg& dst, const CRGBImg& src, float attenuation, const vector<LineSegment> &lines);
HRESULT VtDumpEllipses(CRGBImg& dst, const CRGBImg& src, float attenuation, const vector<EdgeSegment> &edgel_list,  const vector<EllipseSegment> &ellipses, bool useZeroCrossings);
HRESULT VtDumpBitangents(CRGBImg& dst, LineSegment bitangents[2]);
HRESULT VtDumpVanishingPoints(CRGBImg& dst, const CRGBImg& src, int& numValid, float attenuation, const vector<VanishingPoint>& vanishingPoints, const vector<LineSegment>& lineSegments);


// Code to write data to disk
void VtDumpCurves(const wchar_t* filename, const vector<EdgeSegment> &edgelList, const vector< vector<int> > &curves, bool Matlab);
void VtDumpEllipses(const wchar_t* filename, const vector<EllipseSegment> &ellipses, bool Matlab);
void VtDumpEdgels(const wchar_t* filename, const vector<EdgeSegment> &edgelList, bool Matlab, bool extraColumns = false);
void VtDumpBitangents(const wchar_t* filename, const LineSegment bitangents[2], const CVec3f vanishingPoints[2], bool Matlab);

};