//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for extracting edges
//		Uses two ways to extract edges
//		1) Use Difference of Gaussians to generate signed image and detect zero crossings
//		2) Use Steerable filters and max of orientation energy to detect lines and step edges
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
/// <summary> Parameters for edge detection using DoG and Steerable filters </summary>
struct EdgeDetectParams
{
	/// <summary> parameter used for steerable filters, provides a min strength threshold 
	/// for the orientation energy for finding edge segments </summary>
	float orientationEnergyThreshold;
	/// <summary> parameter used for DoG, provides a min strength threshold for the
	/// squared gradient magnitude </summary>
	float gradientMagnitudeThreshold;
	/// <summary> parameter used for steerable filters, provides a min strength threshold 
	/// for the derivative image for finding even edges </summary>
	float evenDerivativeThreshold;
	/// <summary> parameter used for steerable filters, provides a min strength threshold 
	/// for the derivative image for finding odd edges </summary>
	float oddDerivativeThreshold;
	/// <summary> scale of the steerable filters, default is 1.f </summary>
	float steerableFilterScale;
	/// <summary> specifies the number of levels of the DoG pyramid </summary>
	int numLevelsForDoG;
	/// <summary> margin from the edge of the image for computations (usually 4 pix to
	/// allow for sampling etc.)</summary>
	int margin;

	EdgeDetectParams()
	{	
		// Default values
		orientationEnergyThreshold = 400.f;
		gradientMagnitudeThreshold = 4.f;
		evenDerivativeThreshold = 50.f;
		oddDerivativeThreshold = 20.f;
		steerableFilterScale = 1.f;
		numLevelsForDoG = 2;
		margin = 4;
	}
};

/// \ingroup edgedetect
/// <summary> Creates a set of edge segments starting from a difference of gaussian 
/// image and detecting zero crossings in the 2x2 quad of pixels in the DoG image:
/// v00 v01
/// v10 v11
/// See \ref edgeextract "Edge Extraction and Related Operations" for sample outputs</summary>
/// <param name="edgelList"> output: list of detected \link vt::EdgeSegment edge segments \endlink </param>
/// <param name="src"> source image </param>
/// <param name="edParams"> \link vt::EdgeDetectParams edge detection parameters \endlink</param>
HRESULT VtCreateEdgeSegmentListUsingDoG(OUT vector<EdgeSegment>& edgelList, const CImg &src, 
											const EdgeDetectParams& edParams);

/// \ingroup edgedetect
/// <summary> Creates a set of edge segments using second order gaussian steerable filters.
/// The steerable filters are useful as they have the same response to step edges and ridges
/// which can be useful in cases where we have lines (as well as edges in the image).
/// We use a second order quadrature pair of filters (G2 and H2 ref. paper - The Design and 
/// Use of Steerable Filters by Freeman and Adelson). These are the steps in computing an edge response:
///		1) We compute the dominant orientation using the analytical solution for the orientation energy (this
///		only finds one dominant direction as we use the second order filters). 
///		2) Then we use the computed orientation angle to steer the filters in that direction and get filter
///		responses G2 and H2. 
///		3) We then find a peak in the magnitude = (G2*G2 + H2*H2) by searching for a maxima in the direction
///		perpendicular to the direction of the dominant orientation.
///		4) If we find a peak in step (3) that indicates the presence of a line or edge.
///		5) We do a quadratic fit on the magnitude to find the x,y position of the edge reponse up to sub-pixel
///		accuracy
/// See \ref edgeextract "Edge Extraction and Related Operations" for sample outputs</summary>
/// </summary>
/// <param name="edgelList"> output: list of detected \link vt::EdgeSegment edge segments \endlink </param>
/// <param name="src"> src image </param>
/// <param name="edParams"> \link vt::EdgeDetectParams edge detection parameters \endlink</param>
HRESULT VtCreateEdgeSegmentListUsingSteerableFilters(OUT vector<EdgeSegment>& edgelList, const CImg& src, 
											const EdgeDetectParams& edParams);
};