//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      EllipseSegment fitting from linked edge lists
//
//  History:
//      2011/07/07-szeliski
//          Created from earlier Sho/Python port of AWF's Matlab code
//
//      2011/10/31-szeliski
//          Added bi-tangent code.
//
//------------------------------------------------------------------------
#include "vt_edgedetect_common.h"
#pragma once

#ifdef MAKE_DOXYGEN_WORK
void foo();
#endif
namespace vt
{
/// \ingroup edgedetect
/// <summary> Recovered ellipse parameters </summary>
/// <DL><DT> Remarks: </DT></DL>
///		- The equation of the ellipse is
/// \code
///   u = (x-cx)*cos(theta) + (y-cy)*sin(theta)
///   v = (x-cx)*sin(theta) - (y-cy)*cos(theta)
///   (u/ax)^2 + (y/ay)^2 = 1
/// \endcode
///		- To draw the ellipse, for t=[0..2 pi], use
/// \code
///   u = ax*cos(t)
///   v = ay*sin(t)
///   x = u*cos(theta) - v*sin(theta) + cx
///   y = u*sin(theta) + v*cos(theta) + cy
/// \endcode
struct EllipseSegment
{
	/// <summary> X coordinate of ellipse center </summary>
	float cx;
	/// <summary> Y coordinate of ellipse center </summary>
	float cy;
	/// <summary> Length of "x" (major) axis </summary>
	float ax;
	/// <summary> Length of "y" (minor) axis </summary>
	float ay;
	/// <summary> EllipseSegment tilt angle in radians </summary>
	float theta;
	/// <summary> Minimum arc-length (angle) parameter value </summary>
	float min_t;
	/// <summary> Maximum arc-length (angle) parameter value </summary>
	float max_t;
	/// <summary> Average strength of edgels (local gradient magnitudes) </summary>
	float strength;
	/// <summary> Average distance between successive edgels </summary>
	float density;
	/// <summary> Heuristic score for "goodness" of ellipse fit (length * strength / density) </summary>
	float score;
	/// <summary> Number of (inlier) points used in fit </summary>
	int n_points;
	/// <summary> List of constituent points </summary>
	vector<int> points;

	/// <summary> Generate sample points on the ellipse. </summary>
	/// <param name="nPoints"> Number of generated points </param>
	/// <param name="edgels"> List of computed edgels </param>
	/// <param name="l"> Dummy list of edgel indices [0,n_points-1] </param>
	/// <param name="noiseLevel"> Standard deviation of noise added to ellipse </param>
	HRESULT GeneratePoints(int nPoints, vector<EdgeSegment>& edgels, vector<int>& l, float noiseLevel = 0);

	/// <summary> Estimate distance from ellipse and "angle" along arc of ellipse ('t'). </summary>
	/// <param name="edgel"> Edge element being evaluated </param>
	/// <param name="t"> Estimate of "angle" (parameter) along the ellipse [0,2 PI] </param>
	/// <param name="d2"> Squared distance from ellipse </param>
	/// <param name="a"> Angle between the edgel direction and the ellipse tangent </param>
	void td2aFromXY(const EdgeSegment& edgel, float& t, float& d2, float& a) const;
};

/// \ingroup edgedetect
/// <summary> EllipseSegment fitting parameters </summary>
struct FitEllipseParams
{
	/// <summary> Minimum x and y extents (in pixels) of points being fitted </summary>
	float min_extent;
	/// <summary> Minimum number of edgles needed for a valid fit </summary>
	int min_points;
	/// <summary> Maximum distance (in pixels) between edgel and ellipse fit </summary>
	float max_dist;
	/// <summary> Maximum angle (in radians) between edgel and ellipse fit </summary>
	float max_angle;
	/// <summary> Number of robust inlier re-fitting iterations </summary>
	int n_iter;
	/// <summary> Maximum distance between ellipse endpoints for merging </summary>
	float max_endpoint_dist;
	/// <summary> Maximum number of ellipse fragments to consider merging </summary>
	int max_merge_candidates;
	/// <summary> Minimum fraction of merged points </summary>
	float min_merge_frac;

	FitEllipseParams()
	{
		min_extent = 2;
		min_points = 8;
		max_dist   = 0.5f;
		max_angle  = 0.2f;
		n_iter     = 2;
		max_endpoint_dist = 8;
		max_merge_candidates = 12;
		min_merge_frac = 0.9f;
	}
};

/// \ingroup edgedetect
/// <summary> Fit ellipses to all the curves and sort by length/strength 
/// See \ref edgeextract "Edge Extraction and Related Operations" for sample outputs</summary>
/// <param name="ellipses"> Output list of fitted \link vt::EllipseSegment ellipses \endlink and inlier points </param>
/// <param name="curves"> Input list of curves (indices of edgels) </param>
/// <param name="edgelList"> Input list of \link vt::EdgeSegment edge segments \endlink</param>
/// <param name="feParams"> Input \link vt::FitEllipseParams fit parameters \endlink </param>
HRESULT VtFitEllipsesToCurves(OUT vector<EllipseSegment>& ellipses,
    const vector< vector<int> > &curves, const vector<EdgeSegment> &edgelList,
    const FitEllipseParams &feParams);

/// \ingroup edgedetect
/// <summary> Fit bi-tangents to the pair of dominant ellipses
/// See \ref edgeextract "Edge Extraction and Related Operations" for sample outputs</summary>
/// <param name="bitangents"> Output pair outside bitangent \link vt::LineSegment segments \endlink </param>
/// <param name="vanishingPoints"> Output pair vanishing points </param>
/// <param name="ellipses"> Input list of \link vt::EllipseSegment ellipses \endlink </param>
HRESULT VtFindBiTangentLines(OUT LineSegment bitangents[2], OUT CVec3f vanishingPoints[2],
    IN const vector<EllipseSegment>& ellipses);
};