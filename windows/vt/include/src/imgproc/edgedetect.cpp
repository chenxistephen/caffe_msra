//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for extracting edges
//
//  History:
//      2008/7/14-swinder
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_edgedetect.h"

using namespace vt;

//Local function declarations **** //
 
// Initializes the provided edgesegment if a zero crossing exists in the
// 2x2 quad of pixels:
// v00 v01
// v10 v11
//
// return true if a crossing is present and the edgel was intialized
// otherwise returns false
bool 
InitializeEdgeSegmentIfZeroCrossing(OUT EdgeSegment& edgel, int x, int y, 
									float v00, float v01, float v10, float v11,
									float minStrength = 0.f);

// Gets the difference of the base image and the image nLevels up
HRESULT 
GetDifferenceOfGaussianImage(OUT CLumaFloatImg& dog, const CLumaFloatImg& src, 
									int nLevels);

// Steerable filter update by detecting zero crossings to find (even or odd) edge segments
HRESULT 
UpdateInternalWithZeroCrossings(OUT vector<EdgeSegment>& edgelList, CSteerableFilter& steerableFilter,
								int width, int height, eSteerableEdgeType edgeType, const EdgeDetectParams &edParams);

// Steerable filter update to find edge segments as maximas of orientation energy
HRESULT 
UpdateInternal(OUT vector<EdgeSegment>& edgelList, CSteerableFilter& steerableFilter,
				int width, int height, const EdgeDetectParams& edParams);

// ****  local function declarations END//

#define COS_THRESH				0.5


inline UINT fsign(float v) { return ((*(UINT*)&v) & 0x80000000) >> 31; }

//
// Extract an edgel as a (singular) zero crossing pair
//
bool 
	InitializeEdgeSegmentIfZeroCrossing(OUT EdgeSegment& e, int x, int y, 
									float v00, float v01, float v10, float v11,
									float minStrength)
{
	// Test how many sign changes there are - use bitwise operators to
	// extract the sign bits from float values 
	UINT s00 = fsign(v00);
	UINT s01 = fsign(v01);
	UINT s10 = fsign(v10);
	UINT s11 = fsign(v11);
	int ex0 = (s00 ^ s01), ex1 = (s10 ^ s11);
	int ey0 = (s00 ^ s10), ey1 = (s01 ^ s11);
	int n_crossings = ex0 + ex1 + ey0 + ey1;
	if (n_crossings != 2)
		return false;       // ignore double edgel condition...

	// Compute strength from the gradient magnitude
	float gradX, gradY;

	// krishnan - changed to include y condition
	//if (ex0 == ex1)
	if (ex0 == ex1 || ey0 == ey1)
	{
		gradX = 0.5f*((v01 - v00) + (v11 - v10));
		gradY = 0.5f*((v10 - v00) + (v11 - v01));
	}
	else
	{
		gradX = ex0? (v01 - v00): (v11 - v10);
		gradY = ey0? (v10 - v00): (v11 - v01);
	}
	
	e.strength = gradX * gradX + gradY * gradY;

	// 2011/11/10 kramnath - not sure of this criterion
	// float avg = i00+i01+i10+i11;
	// if( (fabs(grad_x)+fabs(gradY)) < minStrength*avg )
	// 	 return false;
	// 2011/11/10 kramnath - added this criterion from Sudipta's code 	
	if(e.strength < minStrength)
		return false;

	// analyzer thinks that k can increment beyond 1, but
	// no actual combination of ex{0,1}, ey{0,1} will allow
	// that due to the n_crossings testt above
#pragma warning( push )
#pragma warning ( disable : 6386 )
	// Compute the two e endpoints
	CVec2f pt[2];
	int k = 0;
	if (ex0)
		pt[k].x = v00 / (v00 - v01), pt[k].y = 0, k++;
	if (ex1)
		pt[k].x = v10 / (v10 - v11), pt[k].y = 1, k++;
	if (ey0)
		pt[k].y = v00 / (v00 - v10), pt[k].x = 0, k++;
	if (ey1)
		pt[k].y = v01 / (v01 - v11), pt[k].x = 1, k++;
#pragma warning( pop )

	// Compute mid-point, and set orientation to difference
	e.x   = float(x) + 0.5f * (pt[0].x + pt[1].x);
	e.y   = float(y) + 0.5f * (pt[0].y + pt[1].y);
	e.n_x = pt[0].x - pt[1].x;
	e.n_y = pt[0].y - pt[1].y;

	// set the end points
	e.endpoint1.x = float(x) + pt[0].x;
	e.endpoint1.y = float(y) + pt[0].y;

	e.endpoint2.x = float(x) + pt[1].x;
	e.endpoint2.y = float(y) + pt[1].y;

	// Reverse the orientation if not compatible with gradient
	if (gradX * e.n_y - gradY * e.n_x < 0.0)
		e.n_x = -e.n_x, e.n_y = -e.n_y;

	// ... the length and theta are set in NormalizeN and ComputeTheta
	e.NormalizeNormalVector();
	e.ComputeThetaFromNormalVector();

	e.phi = e.curvature = e.lineLength = e.sigma = 0.0;  // only used for fitted lines

	return true;
}

HRESULT 
	GetDifferenceOfGaussianImage(CLumaFloatImg& dog, const CLumaFloatImg& src, int nLevels)
{
	VT_HR_BEGIN()

	CLumaFloatPyramid pyramid;

	// the blur filter 
	C14641Kernel k14641;

	// the corresponding up-sampling filter
	float arfKrnlUp0[] = {1.f/8.f, 6.f/8.f, 1.f/8.f};
	float arfKrnlUp1[] = {4.f/8.f, 4.f/8.f};
	C1dKernelSet ksUp;
	VT_HR_EXIT( ksUp.Create(2, 1) );
	VT_HR_EXIT( ksUp.Set(0, -1, 3, arfKrnlUp0) );
	VT_HR_EXIT( ksUp.Set(1, 0, 2,  arfKrnlUp1) );

	// blur the base level
	CLumaFloatImg blur;
	VT_HR_EXIT( blur.Create(src.Width(), src.Height()) );
	VT_HR_EXIT( VtSeparableFilter(blur, src,
				k14641.AsKernel(), k14641.AsKernel()) );

	// create a pyramid
	PYRAMID_PROPERTIES props;
	VT_HR_EXIT( pyramid.Create(blur, &props) );

	// create the output
	VT_HR_EXIT( dog.Create(src.Width(), src.Height()) );

	// upsample successive levels of the pyramid
	for( int i = nLevels-1; i >= 0; i-- )
	{
		CImg& dst = (i==0)? dog : pyramid.GetLevel(i);
		VT_HR_EXIT( VtSeparableFilter(dst, dst.Rect(), pyramid[i+1], 
	                              vt::CPoint(0,0), ksUp) );
	}

	// subtract to form the DoG
	VT_HR_EXIT( VtSubImages( dog, pyramid[0], dog) );

	VT_HR_END()
}

//
// based on some observations the code extracts about 1 edge / 20 pixels
// use this to pre-allocate space in the edgel vector
//
#define PIXELS_PER_EDGE 20

HRESULT 
	vt::VtCreateEdgeSegmentListUsingDoG( OUT vector<EdgeSegment>& edgelList, 
									const CImg& src, const EdgeDetectParams& edParams)
{
	VT_HR_BEGIN()
	if (!src.IsValid())
		VT_HR_EXIT( E_INVALIDSRC );
	
	int x, y;
	int width = src.Width(), height = src.Height();
	
	// Scale/conver the image
	CLumaFloatImg dst;
	VT_HR_EXIT( dst.Create(width, height) );
	VT_HR_EXIT( VtScaleImage(dst, src, 255.0) );

	int margin = edParams.margin;
	VT_HR_EXIT( edgelList.resize(0) );
	VT_HR_EXIT( edgelList.reserve(width*height/PIXELS_PER_EDGE) );
	CLumaFloatImg dogImage;
	VT_HR_EXIT( dogImage.Create(width, height) );

	GetDifferenceOfGaussianImage(dogImage, dst, edParams.numLevelsForDoG);

	for (y = margin; y < height-1-margin; y++)
	{
		const float *pD0 = dogImage.Ptr(y);
		const float *pD1 = dogImage.Ptr(y+1);

		for (x = margin; x < width-1-margin; x++)
		{
			// Find the zero crossing, and if successful, push edgel
			EdgeSegment e;
			if ( InitializeEdgeSegmentIfZeroCrossing(e, x, y, 
				pD0[x], pD0[x+1], pD1[x], pD1[x+1], 
				edParams.gradientMagnitudeThreshold) ) 
			{
				VT_HR_EXIT( edgelList.push_back(e) );
			}
		}
	}

	VT_HR_END()
}

HRESULT 
	vt::VtCreateEdgeSegmentListUsingSteerableFilters(OUT vector<EdgeSegment>& edgelList, const CImg& src, 
											 const EdgeDetectParams& edParams)
{
	VT_HR_BEGIN()

	if (!src.IsValid())
		VT_HR_EXIT( E_INVALIDSRC );

	int width = src.Width();
	int height = src.Height();

	CLumaFloatImg dst;
	VT_HR_EXIT( dst.Create(width, height) );

	// Need single channel scaled float image 
	VT_HR_EXIT( VtScaleImage(dst, src, 255.0) );

	// zero out edgeList
	VT_HR_EXIT( edgelList.resize(0) );

	// Setup second order steerable filter
	CSteerableFilter steerableFilter;
	VT_HR_EXIT( steerableFilter.Create(width, height, FilterBoth, SecondOrder, edParams.steerableFilterScale) );
	VT_HR_EXIT( steerableFilter.Update(dst) );

	VT_HR_EXIT( UpdateInternal(edgelList, steerableFilter, width, height, edParams) );

	/// Note: Another potential way to detect edge segments using steerable filters is to generate a 
	/// signed image to detect zero crossings. These are the steps involved to do that using steerable filters:
	/// 1) As before, get G2 and H2 filter responses for the image after detecting the dominant orientation
	/// 2) The magnitude of the response is given by G2*Sin(phi)+H2*Cos(phi) where phi is the phase difference between the quadrature pair G2 and H2
	/// 3) Since phi is 0 or 90 deg at the step or ridge edge and varies significantly on either ends of this we detect zero crossings in the 
	/// derivative image for either the even (0 deg) edges or the odd (90 deg) edges. So if the user needs even or odd edge in 
	/// the "edgeType" argument then the call (of finding zero crossings) below can be used  to detect edges. If all edges are required, then use
	/// the method described in the Freeman paper (UpdateInternal).
	/// NOTE: the call below is not fully tested and may yield poor results
	// eSteerableEdgeType eType = eSteerableEdgeTypeAllEdges;
	// VT_HR_EXIT( UpdateInternalWithZeroCrossings(edgelList, steerableFilter, width, height, eType, edParams) );
	
	VT_HR_END()
}

#define PI 3.14159265358979323846

HRESULT 
	UpdateInternal(vector<EdgeSegment>& edgelList, CSteerableFilter& steerableFilter, int width, 
					int height, const EdgeDetectParams& edParams)
{
	VT_HR_BEGIN()

	CFloatImg magImg;
	VT_HR_EXIT( magImg.Create(width, height) );

	for(int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			// Get G2 and H2
			float orientMag, th;
			steerableFilter.GetLocalMagnitudeAndOrientation(x, y, orientMag, th);
			float g2 = steerableFilter.GetPixelEvenFilter(x, y, th);
			float h2 = steerableFilter.GetPixelOddFilter(x, y, th);

			// Squared magnitude of the quadrature pair
			magImg(x,y) = g2*g2 + h2*h2;
		}
	}
	int margin = edParams.margin;
	for(int y = margin; y < height-1-margin; y++)
	{
		for (int x = margin; x < width-1-margin; x++)
		{
			float m3, m2, m1;
			// Get the dominant orientation by using the analytical solution
			float orientMag, th;
			steerableFilter.GetLocalMagnitudeAndOrientation(x, y, orientMag, th);
			float fdx = cos(th);
			float fdy = -sin(th);

			m3 = magImg(x,y);
			// check magnitude 
			if(m3 < edParams.orientationEnergyThreshold)
			{
				continue;
			}
			VtSampleBicubic(magImg, (float)(x - fdx), (float)(y - fdy), (float *)NULL, &m1);
			VtSampleBicubic(magImg, (float)(x + fdx), (float)(y + fdy), (float *)NULL, &m2);

			if((m3 > m1) && (m3 > m2))
			{
				float mm;
				float d = VtQuadraticFit1D(m1, m3, m2, mm);

				EdgeSegment e;
				e.x = x + d * fdx;
				e.y = y + d * fdy;
				e.n_x = fdy; // not sure why this needs my normal rotated by 90degrees
				e.n_y = -fdx;
				// For phi - Get G2 and H2
				float g2 = steerableFilter.GetPixelEvenFilter(x, y, th);
				float h2 = steerableFilter.GetPixelOddFilter(x, y, th);
				// Phase not used, but documented
				e.phi = atan2(h2, g2) * (float) (180/PI); // in degrees
				e.curvature = e.lineLength = e.sigma = 0.0;  // only used for fitted lines
				e.strength = m3;
				e.NormalizeNormalVector();
				e.ComputeThetaFromNormalVector();
				// Compute the end points and store them
				e.ComputeEndpoints();
				VT_HR_EXIT( edgelList.push_back(e) );
			}
		}
	}

	VT_HR_END()
}


// For experimental testing only, may yield poor results
HRESULT 
	UpdateInternalWithZeroCrossings(vector<EdgeSegment>& edgelList, CSteerableFilter& steerableFilter, int width, 
									int height, eSteerableEdgeType edgeType, const EdgeDetectParams& edParams)
{
	VT_HR_BEGIN()

	CFloatImg derivativeImg;
	VT_HR_EXIT( derivativeImg.Create(width, height) );
	CFloatImg phiImg;
	VT_HR_EXIT( phiImg.Create(width, height) );
	CFloatImg magImg;
	VT_HR_EXIT( magImg.Create(width, height) );

	for(int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			float mag, th;
			steerableFilter.GetLocalMagnitudeAndOrientation(x,y,mag,th);

			float g2 = steerableFilter.GetPixelEvenFilter(x, y, th);
			float h2 = steerableFilter.GetPixelOddFilter(x, y, th);

			phiImg(x,y) = atan2(h2, g2);
			magImg(x,y) = sqrt(g2*g2 + h2*h2);

			if(edgeType == eSteerableEdgeTypeBarEdges)
				derivativeImg(x,y) = h2;
			else 
				derivativeImg(x,y) = g2; 

		}
	}

	int margin = edParams.margin;
	float threshold = (edgeType == eSteerableEdgeTypeBarEdges) ? edParams.evenDerivativeThreshold : edParams.oddDerivativeThreshold;
	for(int y = margin; y < height-1-margin; ++y)
	{
		for (int x = margin; x < width-1-margin; ++x)
		{
			EdgeSegment e;
			// Threshold is currently on the derivative image
			// as threshold on magnitude gives bad responses,
			// maybe there is a better way to do this
			if (fabs(derivativeImg(x,y)) > threshold &&  InitializeEdgeSegmentIfZeroCrossing(e, x, y, 
				derivativeImg(x,y), derivativeImg(x+1, y), derivativeImg(x, y+1), derivativeImg(x+1, y+1), 
				0)) 
			{
				e.phi = (float)(phiImg(x,y) * 180 / PI);
				VT_HR_EXIT(edgelList.push_back(e) );
			}
		}
	}

	VT_HR_END()

}

