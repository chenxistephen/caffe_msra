//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for image resampling
//
//  History:
//      2004/11/08-mattu
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_image.h"

namespace vt {

//+-----------------------------------------------------------------------
//
// Functions to point sample an image and return an interpolated value
//
// These functions all take an image with one or more bands and write the 
// sampled value into the array pRtn which has size equal to the number of bands 
// in the source image. The default value is used when the point (fX, fY) lies
// outside the source image. The default pDefault can be NULL in which case the 
// elements of the return value are set to zero.
// Note that if pDefault is NULL then make sure that you cast NULL to the template
// type T *, otherwise you will have compile errors.
// note that bicubic interpolation can result in negative pixel values
//
//------------------------------------------------------------------------
template <class T>
void VtSampleNearest(const CTypedImg<T> &imgSrc, float fX, float fY, 
                     T *pDefault, T *pRtn);
template <class T>
void VtSampleBilinear(const CTypedImg<T> &imgSrc, float fX, float fY, 
                      T *pDefault, T *pRtn);
template <class T>
void VtSampleBicubic(const CTypedImg<T> &imgSrc, float fX, float fY, 
                     T *pDefault, T *pRtn);


// Same functions as above but faster and without any testing of the range of x and y.
// You must allow a 3 pixel border on the source image for interpolation.
template <class T>
void VtSampleBilinearNoRangeTests(const CTypedImg<T> &imgSrc, float fX, float fY, T *pRtn);
template <class T>
void VtSampleBicubicNoRangeTests(const CTypedImg<T> &imgSrc, float fX, float fY, T *pRtn);


// these functions sample along a row starting from x,y at intervals of dx,dy for count pixels.
// for whichever points are outside the image, the pixel value pointed to by default is used.
// multiple bands are supported, in which case default (if not null) must point to an array and
// prtn must point to an array large enough to hold bands * count values of type T.
template <class T>
void VtSampleRowBilinear(const CTypedImg<T> &imgSrc, float fX, float fY, float fDX, float fDY, int iCount,
					T *pDefault, T *pRtn);

template <class T>
void VtSampleRowBicubic(const CTypedImg<T> &imgSrc, float fX, float fY, float fDX, float fDY, int iCount,
					T *pDefault, T *pRtn);

// these functions sample the source image to return a patch of size given by the size of imgdest using the
// affine transform provided to translate from the patch coordinates with origin in the top left into the source
// image coordinate frame.

template <class T>
void VtSamplePatchBilinear(const CTypedImg<T> &imgSrc, const CMtx3x3f &mAff, CTypedImg<T> &imgDest, T *pDefault);

template <class T>
void VtSamplePatchBicubic(const CTypedImg<T> &imgSrc, const CMtx3x3f &mAff, CTypedImg<T> &imgDest, T *pDefault);


//+-----------------------------------------------------------------------------
// 
// Function: VtResampleNearest
// 
// Synopsis: resample an image using nearest-neighbor resampling.
// 
//------------------------------------------------------------------------------
inline int VtResampleNearestCoord(int c, float s, float t)
{ return F2I( float(c)*s+t ); }

inline vt::CPoint VtResampleNearestAddr(vt::CPoint& pt, float sx, float tx,
										float sy, float ty)
{
	return vt::CPoint( VtResampleNearestCoord(pt.x, sx, tx), 
					   VtResampleNearestCoord(pt.y, sy, ty) );
}

inline vt::CRect
GetRequiredSrcRectResampleNearest(const vt::CRect& rctDst, float sx, float tx,
								 float sy, float ty )
{
	int x0 = VtResampleNearestCoord(rctDst.left,     sx, tx);
	int y0 = VtResampleNearestCoord(rctDst.top,      sy, ty);
	int x1 = VtResampleNearestCoord(rctDst.right-1,  sx, tx) + 1;
	int y1 = VtResampleNearestCoord(rctDst.bottom-1, sy, ty) + 1;
	return vt::CRect(VtMin(x0,x1), VtMin(y0,y1), VtMax(x0,x1), VtMax(y0,y1));
}

inline vt::CRect
GetAffectedDstRectResampleNearest(const vt::CRect& rctSrc, float sx, float tx,
								 float sy, float ty )
{ 
	int x0 = (int)floor((float(rctSrc.left) - tx - 0.5f) / sx);
	int y0 = (int)floor((float(rctSrc.top)  - ty - 0.5f) / sy);
	int x1 = (int)ceil((float(rctSrc.right-1)  - tx + 0.5f) / sx);
	int y1 = (int)ceil((float(rctSrc.bottom-1) - ty + 0.5f) / sy);
	return vt::CRect(VtMin(x0,x1), VtMin(y0,y1), VtMax(x0,x1), VtMax(y0,y1));
}

inline vt::CRect
GetResultingDstRectResampleNearest(const vt::CRect& rctSrc, float sx, float tx,
								 float sy, float ty )
{ 
	return GetAffectedDstRectResampleNearest(rctSrc, sx, tx, sy, ty);
}

void VtResampleNearest(const CImg& src, float sx, float tx, float sy, float ty, 
					   CImg& dst);

#if !defined(VT_GCC)
//+-----------------------------------------------------------------------------
//
// Class: CResampleNearestTransform
// 
// Synposis: IImageTransform implementation of nearest-neighbor resampling
// 
//------------------------------------------------------------------------------
class CResampleNearestTransform: 
    public CImageTransformUnaryGeo<CResampleNearestTransform, false>
{
public:
	CResampleNearestTransform()
	{}

	CResampleNearestTransform(float sx, float tx, float sy, float ty)
	{ Initialize(sx, tx, sy, ty); }

	void Initialize(float sx, float tx, float sy, float ty)  
	{ m_fSx = sx, m_fTx = tx, m_fSy = sy, m_fTy = ty; }

	vt::CRect GetRequiredSrcRect(const vt::CRect& rDst)
	{ return GetRequiredSrcRectResampleNearest(rDst, m_fSx, m_fTx, m_fSy, m_fTy); }

	vt::CRect GetAffectedDstRect(const vt::CRect& rSrc)
	{ return GetAffectedDstRectResampleNearest(rSrc, m_fSx, m_fTx, m_fSy, m_fTy); }

	vt::CRect GetResultingDstRect(const vt::CRect& rSrc)
	{ return GetResultingDstRectResampleNearest(rSrc, m_fSx, m_fTx, m_fSy, m_fTy); }

	HRESULT Transform(CImg* pimgDst, const vt::CRect& rctDst, 
					  const CImg& src, const vt::CPoint& ptSrc)
	{
		float tx = float(rctDst.left)*m_fSx+m_fTx-float(ptSrc.x);
		float ty = float(rctDst.top)*m_fSy+m_fTy-float(ptSrc.y);
		VtResampleNearest(src, m_fSx, tx, m_fSy, ty, *pimgDst);
		return S_OK;
	}

	virtual HRESULT Clone(ITaskState **ppState)
	{ 
		return CloneTaskState(ppState, VT_NOTHROWNEW CResampleNearestTransform(
			m_fSx, m_fTx, m_fSy, m_fTy));
    }

protected:
	float m_fSx, m_fSy, m_fTx, m_fTy;
};
#endif

};