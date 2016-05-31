//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description: 3x3 transform functions used in this library
//
//------------------------------------------------------------------------------
#include "stdafx.h"

#include "transforms.h"

using namespace vt;

void ComputeSimParams(FEAT_SIMILARITY& sim, const CMtx3x3f& mtx, 
                      int iFrameWidth, int iFrameHeight)
{
    PointMatch pm[4];
    TransformCorners(pm[0].p1, pm[1].p1, pm[2].p1, pm[3].p1,
                     mtx, iFrameWidth, iFrameHeight);
    
	float fW2 = float(iFrameWidth)/2.f;
	float fH2 = float(iFrameHeight)/2.f;

    pm[0].p0 = CVec2f(-fW2, fH2);
    pm[1].p0 = CVec2f( fW2, fH2);
    pm[2].p0 = CVec2f(-fW2,-fH2);
    pm[3].p0 = CVec2f( fW2,-fH2);
    
    // find the best matching similarity matrix for this frame to frame
    // transformation
    CMtx3x3d mS;
    VtSimilarityFromPointMatches2D(mS, pm, 4);
    
    // keep only the translation and scale from this similarity
    float s = (float)sqrt(mS(0,0)*mS(0,0)+mS(0,1)*mS(0,1));
    float r = (float)atan2(mS(1,0), mS(0,0));
    sim = FEAT_SIMILARITY((float)mS(0,2), (float)mS(1,2), s, r);
}

CMtx3x3f
ClampSimilarityTransform(const CMtx3x3f& mS, float fCrop, int iW, int iH,
						 const CVec2f* pRSOff, int iRSOffCnt)
{
	// compute the offset due to rolling shutter
	float fLeftEx = 0, fRightEx = 0, fTopEx = 0, fBotEx = 0;
	if( iRSOffCnt && pRSOff!=NULL )
	{
		int crop_bound = F2I(float(iRSOffCnt)*(1.f-fCrop));

		fLeftEx = fRightEx = pRSOff[0].x;
		fTopEx  = pRSOff[0].y;
		fBotEx  = pRSOff[iRSOffCnt-1].y;
		for( int y = 1; y < iRSOffCnt; y++ )
		{
			const CVec2f& v = pRSOff[y];
			if( y < crop_bound )
			{
				fTopEx = VtMin(fTopEx, v.y);
			}
			else if( y >= iRSOffCnt-crop_bound )
			{
				fBotEx = VtMax(fBotEx, v.y);
			}
			fLeftEx  = VtMin(fLeftEx, v.x);
			fRightEx = VtMax(fRightEx, v.x);
		}
	}

	// compute scale
	float s = sqrtf(mS(0,0)*mS(0,0)+mS(0,1)*mS(0,1));

	// scale the crop box
	fCrop *= s;

	// compute max allowed translate
	float fMaxC = (1.f-fCrop)/2.f;
	float fMaxX = float(iW)*fMaxC;
	float fMaxY = float(iH)*fMaxC;

	float tx, ty, r;
	if( fabs(fLeftEx) > fMaxX || fabs(fRightEx) > fMaxX ||
		fabs(fTopEx)  > fMaxY || fabs(fBotEx)   > fMaxY )
	{
		// nothing we can do because rolling-shutter calc pushes us OOB
		tx = ty = r = 0;
	}
	else
	{
		// clamp translate components such that they don't go out-of-bounds 
		// include the affect of rolling shutter
		tx  = VtClamp(mS(0,2), -fMaxX-fLeftEx, fMaxX-fRightEx);
		ty  = VtClamp(mS(1,2), -fMaxY+fBotEx,  fMaxY+fTopEx);

		// compute max translate in normalized image dimensions 
		float fMaxD = VtMax(VtMax(-(tx+fLeftEx), tx+fRightEx) / float(iW),
							VtMax(-(ty-fBotEx),  ty-fTopEx)   / float(iH));

		// compute max allowed angle
		float fMaxA = float(VT_PI/4.)-
			acosf((0.5f-fMaxD)/(fCrop*float(VT_SQRT1_2)));
		VT_ASSERT( fMaxA >= 0 );

		// clamp the rotation
		r = atan2f(mS(1,0), mS(0,0));
		r = VtClamp(r, -fMaxA, fMaxA);
	}

	// reform the matrix
	return MatrixFromSimParams(FEAT_SIMILARITY(tx,ty,s,r));
}

//------------------------------------------------------------------------------
// end