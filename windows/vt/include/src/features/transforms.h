//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description: 3x3 transform functions used in this library
//
//------------------------------------------------------------------------------
#pragma once

#include "vt_features_base.h"

inline void
TransformCorners(vt::CVec2f& tl, vt::CVec2f& tr, vt::CVec2f& bl, vt::CVec2f& br,
                 const vt::CMtx3x3f& mtx, int iFrameWidth, int iFrameHeight)
{  
	float fW2 = float(iFrameWidth)/2.f;
	float fH2 = float(iFrameHeight)/2.f;

	tl = vt::CVec3f(mtx*vt::CVec3f(-fW2, fH2, 1)).Dehom();
    tr = vt::CVec3f(mtx*vt::CVec3f( fW2, fH2, 1)).Dehom();
    bl = vt::CVec3f(mtx*vt::CVec3f(-fW2,-fH2, 1)).Dehom();
    br = vt::CVec3f(mtx*vt::CVec3f( fW2,-fH2, 1)).Dehom();
}

inline vt::CMtx3x3f 
MatrixFromSimParams(const FEAT_SIMILARITY& sim)
{
    float sin_z, cos_z;
    vt::VtSinCos(sim.r, &sin_z, &cos_z);
    cos_z *= sim.s;
    sin_z *= sim.s;
    return vt::CMtx3x3f(cos_z, -sin_z, sim.tx, sin_z, cos_z, sim.ty, 0, 0, 1.f);
}

inline vt::CMtx3x3f
MatrixCenterToTopLeft(int iW, int iH)
{ return vt::CMtx3x3f(1.f, 0, float(iW)/2.f, 0, -1.f, float(iH)/2.f, 0, 0, 1.f); }

inline vt::CMtx3x3f
MatrixTopLeftToCenter(int iW, int iH)
{ return vt::CMtx3x3f(1.f, 0, -float(iW)/2.f, 0, -1.f, float(iH)/2.f, 0, 0, 1.f); }

void 
ComputeSimParams(FEAT_SIMILARITY& sim, const vt::CMtx3x3f& mtx, 
                 int iFrameWidth, int iFrameHeight);

vt::CMtx3x3f
ClampSimilarityTransform(const vt::CMtx3x3f& mS, float fCrop, int iW, int iH,
						 const vt::CVec2f* pRSOff=NULL, int iRSOffCnt = 0);
