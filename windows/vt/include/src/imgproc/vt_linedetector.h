//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routine for extracting curves and lines
//
//  History:
//      2011/11/10-kramnath
//          Created
//
//------------------------------------------------------------------------


#include "vt_edgedetect_common.h"

namespace vt {
#pragma once

/// \ingroup edgedetect
/// <summary> Line Detector Params </summary>
struct LineDetectorParams
{
    /// <summary> Min count of edge segments in a line </summary>
    int     minEdgelCount;				
    /// <summary> Threshold used to prune smaller line segments </summary>
    int		minLineSegmentLength;		
    /// <summary> Distance Threshold (in pixels) used to decide whether two line segments can be merged.  </summary>
    float	lineSegmentMergingDistThreshold;  
    /// <summary> Angle Threshold used to decide if two line segments can be merged </summary>
    float	lineSegmentMergingAngleThreshold; 
    /// <summary> Threshold for fitting line equation to a set of edgelsd </summary>
    float  lineEqnFitThreshold;		
    /// <summary> Minimum threshold to ensure edges are oriented in the same direction between neighbouring
    /// edge segments. Used only when edges are detected without zero crossings </summary>
    float edgeOrientationThreshold;
    /// <summary> Threshold to indicate when the line simplification thing should terminate </summary>
    float lineSimplificationThreshold; 

    LineDetectorParams()
    {
        // default values
        minLineSegmentLength   = 4;
        minEdgelCount		   = 4;
        lineEqnFitThreshold    = 0.75f;
        edgeOrientationThreshold = 0.7f;
        lineSimplificationThreshold = 2.f;
        lineSegmentMergingDistThreshold  = 4.f;		// in pixels
        lineSegmentMergingAngleThreshold = 0.99619f; //COS(5.0 degrees)

    }
    ~LineDetectorParams() {};

};

/// \ingroup edgedetect
/// <summary> Finds connected edge segments and forms curves from them. 
/// The algorithm is greedy - it raster scans the image and for each pixel
/// it tries to walk on both sides of the edge segment to find connected edge segments
/// if two neighboring edge segments are not oriented similar then the curve is ended
/// at that point
/// See \ref edgeextract "Edge Extraction and Related Operations" for sample outputs</summary>
/// </summary>
/// <param name="curves"> output: set of curves are the output (each curve is indices into
/// the edgelList)  </param>
/// <param name="edgelList">  input set of \link vt::EdgeSegment edge segments \endlink</param>
/// <param name="roundToNearestQuadOrPixel"> true if the edge segments in edgelList are detected using zero 
/// crossings (e.g. when using DoG for edge detection) so we have to round to the nearest quad of pixels.
/// false if not (e.g. with steerable filters where we look for max of an energy function to define a contour),
/// in which case we round to the nearest pixel. </param>
/// <param name="width"> width of the input image </param>
/// <param name="height"> height of the input image </param>
/// <param name="ldParams"> \link vt::LineDetectorParams line detector parameters \endlink </param>
HRESULT VtFindCurveSegmentsFromEdges(OUT vector< vector<int> >& curves,  vector<EdgeSegment>& edgelList, bool 
								roundToNearestQuadOrPixel, int width, int height, const LineDetectorParams& ldParams);
	


/// \ingroup edgedetect
/// <summary> Find lines from curves. 
/// This is an implementation of the Ramer Douglas Peucker algorithm for line simplification
/// At each stage a curve is examined and split at the point that is the farthest from 
/// the line joining the end points of the curve. If the distance is less than a threshold
/// then the simplification stops otherwise it recursively operates on the two sub curves
/// split at the farthest point.
/// See \ref edgeextract "Edge Extraction and Related Operations" for sample outputs</summary>
/// </summary>
/// <param name="lineSegments"> output: set of \link vt::LineSegment line segments \endlink (two end points of a line) vector as output </param>
/// <param name="edgelList">  input set of \link vt::EdgeSegment edge segments \endlink</param>
/// <param name="curves"> input curves to detect lines from </param>
/// <param name="ldParams"> \link vt::LineDetectorParams line detector parameters \endlink </param>
HRESULT VtFindLinesFromCurves(OUT vector<LineSegment>& lineSegments, const vector<EdgeSegment>& edgelList,  
							const vector< vector<int> >& curves, const LineDetectorParams& ldParams);


/// \ingroup edgedetect
/// <summary> Break curves into smaller curves by looking for inflection points.
/// This module used the curvature of the curves to split them (based on EdgeSegment.curvature value).
/// </summary>
/// </summary>
/// <param name="outputCurves"> output: set of split curves </param>
/// <param name="edgelList">  input set of \link vt::EdgeSegment edge segments \endlink</param>
/// <param name="curves"> input curves to detect lines from </param>
/// <param name="threshold"> the curvature values should be above this threshold (default = 0.4f) </param>
HRESULT VtBreakCurveBasedOnCurvature(vector< vector<int> > &outputCurves, const vector<EdgeSegment> &edgelList, 
                                  const vector< vector<int> > &curves, float threshold = 0.4f);
};