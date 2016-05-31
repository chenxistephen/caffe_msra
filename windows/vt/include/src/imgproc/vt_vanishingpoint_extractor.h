//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routine for extracting vanishing points
//
//  History:
//      2011/11/25-kramnath
//          Created
//
//------------------------------------------------------------------------


#include "vt_edgedetect_common.h"

namespace vt {
#pragma once

/// \ingroup edgedetect
/// <summary> Parameters for vanishing point extraction </summary>
struct VPExtractorParams
{
    ///<summary> number of requested vanishing points. In most cases
    /// the number of vanishing points that can be reliably estimated
    /// in an image is about 2-3. If more than these numbers are 
    /// requested then spurious points may result. </summary> 
	int numVPoints;
    ///<summary> threshold in cos of the angle in deg for fitting
    /// lines to vanishing points while using ransac and least sq.
    /// This angle threshold depends on the noise level in the image.
    /// Do not vary this parameter unless you have a better estimate.
    ///</summary> 
    float ransacAndLeastSqCosAngleThreshold;
    /// <summary> threshold in cos of the angle in deg for fitting
    /// lines to vanishing points once a good estimate of VP is found.
    /// This angle threshold depends on the noise level in the image.
    /// Do not vary this parameter unless you have a better estimate.
    /// </summary> 
    float fitLinesToVPCosAngleThreshold;
    /// <summary> threshold in cos of the angle in deg for fitting
    /// lines to prior estimate of the vanishing point. This threshold
    /// can be relaxed based on how confident the user is about the 
    /// initial estimate that is provided to VtExtractVanishingPointUsingPrior.
    /// This angle threshold depends on the noise level in the image.
    /// Do not vary this parameter unless you have a better estimate.
    /// </summary> 
    float fitLinesToPriorVPCosineAngleThreshold;
    /// <summary> minimum length of the line for it to be considered
    /// for vanishing point estimation.
    /// </summary>
    float minLength;
    /// <summary> min number of ransac tries for selecting the 
    /// vanishing point. The number of tries is usually the max of
    /// this value and number of lines</summary>
    int minRansacTries;
    /// <summary> the minimum number of inliers required for a vanishing
    /// point to be determined as valid </summary>
    int inlierThreshold;

	VPExtractorParams()
	{
        Initialize();
	}

    VPExtractorParams(int imgWidth, int imgHeight)
	{
        Initialize();
        int maxDim = VtMax(imgWidth, imgHeight);
        //  calculated as a function of image dimension
        minLength = VtMax(0.01f*maxDim, 4.f); 
        // inlier threshold found by fitting a line to the
        // data points collected by Sudipta in his code
        int y = (int) ((0.061668) * maxDim + 1.7062);
        inlierThreshold = (y < 30) ? 30 : y; 
        inlierThreshold = (y > 300) ? 300 : y;
	}

    void Initialize()
    {
		numVPoints = 3;
		ransacAndLeastSqCosAngleThreshold = 0.99998f; //COS(1 degrees)
		fitLinesToVPCosAngleThreshold = 0.99619f; //COS(5.0 degrees)
        fitLinesToPriorVPCosineAngleThreshold = 0.8660f; //COS(30.0 degrees) 
        minLength = 4.f; 
        minRansacTries = 1000; 
        inlierThreshold = 30; 
    }
};

/// \ingroup edgedetect
/// <summary> Extract a certain number of vanishing points in an image
/// given the set of lines in an image. The estimation is done using ransac
/// followed by a least squares fit to the inliers. Typically, about 2-3 vanishing
/// points exist in an image, any number more than that may result in spurious
/// points.
/// See \ref edgeextract "Edge Extraction and Related Operations" for sample outputs</summary>
/// <param name="vanishingPoints"> Output set of \link vt::VanishingPoint vanishing points \endlink that are extracted</param>
/// <param name="lineSegments">  Input set of \link vt::LineSegment line segments \endlink  in an image </param>
/// <param name="vpParams"> \link vt::VPExtractorParams Vanishing point extractor parameters \endlink </param>
HRESULT 
    VtExtractVanishingPoints(
    vector<VanishingPoint>& vanishingPoints, const vector<LineSegment>& lineSegments, 
    const VPExtractorParams& vpParams );

/// \ingroup edgedetect
/// <summary> Extracts a vanishing point in an image given a prior (approximate)
/// estimate of a vanishing point. The priorVanishingPoint can also be a direction
/// to a vanishing point at infinity. For example, to find the vertical vanishing
/// point the prior can be [1 0 0] and for the horizontal [0 1 0]. The estimation is
/// done using ransac followed by a least squares fit to the inliers. User can also
/// specify a region of interest, if specified, only lines within that region of interest
/// will be used for computation. This can be used for example to specify one facade of
/// a building for which we need a vanishing point estimation.</summary>
/// <param name="vanishingPoint"> Output \link vt::VanishingPoint vanishing point \endlink estimate</param>
/// <param name="priorVanishingPoint">  Input a prior \link vt::VanishingPoint vanishing point \endlink or direction in 
/// homogeneous coordinates. </param>
/// <param name="lineSegments"> Input set of \link vt::LineSegment line segments \endlink </param>
/// <param name="vpParams"> \link vt::VPExtractorParams Vanishing point extractor parameters \endlink </param>
/// <param name="roiRect"> Region of interest for limiting the lines considered for 
/// vanishing point extraction </param>
HRESULT
    VtExtractVanishingPointUsingPrior(
    VanishingPoint& vanishingPoint, const CVec3d& priorVanishingPoint, const vector<LineSegment>& lineSegments, 
    const VPExtractorParams& vpParams, const CRect* roiRect= NULL);

/// \ingroup edgedetect
/// <summary> Extracts a vanishing point in an image on the horizon line, given the line 
/// equation of the horizon line. For example, if the horizon line is considered to be 
/// at the center of the image then the line equation would be [0 1 -h/2]. The estimation 
/// is done using ransac followed by a least squares  fit to the inliers. User can also 
/// specify a region of interest, if specified, only lines within that region of interest 
/// will be used for computation. </summary>
/// <param name="vanishingPoint"> Output \link vt::VanishingPoint vanishing point \endlink estimate</param>
/// <param name="horizonLine">  Input horizon line equation. </param>
/// <param name="lineSegments"> Input set of \link vt::LineSegment line segments \endlink </param>
/// <param name="vpParams"> \link vt::VPExtractorParams Vanishing point extractor parameters \endlink </param>
/// <param name="roiRect"> Region of interest for limiting the lines considered for 
/// vanishing point extraction </param>
HRESULT
    VtExtractVanishingPointOnHorizon(
    VanishingPoint& vanishingPoint, const CVec3d& horizonLine, const vector<LineSegment>& lineSegments, 
    const VPExtractorParams& vpParams, const CRect* roiRect = NULL);

HRESULT
    VtDoMetricRectification(CImg& rectifiedImage, const CImg& src, const vector<VanishingPoint>& vanishingPoints, float focalLengthInPixels);

};