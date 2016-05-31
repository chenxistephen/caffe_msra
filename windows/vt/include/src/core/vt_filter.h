//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for image filtering
//
//  History:
//      2004/11/08-swinder
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_separablefilter.h"

#ifndef VT_NO_XFORMS
#include "vt_transform.h"
#endif

namespace vt {

/// \ingroup filter
/// <summary>Applies a Gaussian filter kernel to the source image.  The rctDst 
/// and ptSrc agruments enable the filtering operation to be applied to a 
/// sub-region of the source. For more control over the Gaussian filter see
/// \link Create1dGaussKernel() Create1dGaussKernel \endlink and
/// \link VtSeparableFilter() VtSeparableFilter \endlink .</summary>
/// <param name="imgDst"> Output filtered image. </param>
/// <param name="rctDst"> Indicates the location and dimensions of the output image. </param>
/// <param name="imgSrc"> Input image. </param>
/// <param name="ptSrc"> The origin of the input image </param>
/// <param name="fSigma">Controls the amount of smoothing. Large sigma yields more
/// smoothing, but a more expensive to apply filter kernel.</param>
/// <param name="ex"> The \link IMAGE_EXTEND extension mode \endlink to apply 
/// when the filter kernel requires support outside of the source image.  The
/// default is Extend.</param>
/// <returns> 
/// - S_OK on success
/// - E_INVALIDSRC if imgSrc is not valid
/// - E_INVALIDDST if imgDst shares memory with imgSrc
/// - E_INVALIDARG if conversion from imgSrc to imgDst type is not supported
/// - E_OUTOFMEMORY if imgDst allocation was necessary and failed 
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///	- Implements all \ref stdconvrules.
/// - The filter can not be applied in place; the provided input and
/// output images must be different.

inline HRESULT
VtGaussianSmooth(CImg &imgDst, const vt::CRect& rctDst, 
				 const CImg &imgSrc, vt::CPoint ptSrc, 
				 float fSigma, const IMAGE_EXTEND& ex = IMAGE_EXTEND(Extend))
{
    C1dKernel k;
    HRESULT hr = Create1dGaussKernel(k, fSigma);
	if( hr == S_OK )
	{
		hr = VtSeparableFilter(imgDst, rctDst, imgSrc, ptSrc, k, k, ex);
	}
	return hr;
}
/// \ingroup filter
/// <summary>Applies a Gaussian filter kernel to the source image. For more 
/// control over the Gaussian filter see 
/// \link Create1dGaussKernel() Create1dGaussKernel \endlink and
/// \link VtSeparableFilter() VtSeparableFilter \endlink .</summary>
/// <param name="imgDst"> Output filtered image. </param>
/// <param name="imgSrc"> Input image. </param>
/// <param name="fSigma">Controls the amount of smoothing. Large sigma yields more
/// smoothing, but a more expensive to apply filter kernel.</param>
/// <param name="ex"> The \link IMAGE_EXTEND extension mode \endlink to apply 
/// when the filter kernel requires support outside of the source image.  The
/// default is Extend.</param>
/// <returns> 
/// - S_OK on success
/// - E_INVALIDSRC if imgSrc is not valid
/// - E_INVALIDDST if imgDst shares memory with imgSrc
/// - E_INVALIDARG if conversion from imgSrc to imgDst type is not supported
/// - E_OUTOFMEMORY if imgDst allocation was necessary and failed 
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///	- Implements all \ref stdconvrules.
/// - The filter can not be applied in place; the provided input and
/// output images must be different.

inline HRESULT
VtGaussianSmooth(CImg &imgDst, const CImg &imgSrc, float fSigma,
			     const IMAGE_EXTEND& ex = IMAGE_EXTEND(Extend))
{
	return VtGaussianSmooth(imgDst, imgSrc.Rect(), imgSrc, vt::CPoint(0,0), 
							fSigma, ex);
}

inline HRESULT 
VtGaussianBesselSmooth(CImg &imgDst, const vt::CRect& rctDst, 
					   const CImg &imgSrc, vt::CPoint ptSrc,
					   float fSigma, const IMAGE_EXTEND& ex = IMAGE_EXTEND(Extend))
{
    C1dKernel k;
    HRESULT hr = Create1dGaussBesselKernel(k, fSigma);
	if( hr == S_OK )
	{
		hr = VtSeparableFilter(imgDst, rctDst, imgSrc, ptSrc, k, k, ex);
	}
	return hr;
}

inline HRESULT 
VtGaussianBesselSmooth(CImg &imgDst, const CImg &imgSrc, float fSigma,
					   const IMAGE_EXTEND& ex = IMAGE_EXTEND(Extend))
{
	return VtGaussianBesselSmooth(imgDst, imgSrc.Rect(), imgSrc, vt::CPoint(0,0), 
							      fSigma, ex);
}

/// \ingroup filter
/// <summary>Applies a Laplacian filter kernel to the source image.</summary>
/// <param name="imgDst"> Output filtered image. </param>
/// <param name="imgSrc"> Input image. </param>
/// <param name="fScale"> Scales the output values by this value.</param>
/// <param name="pSrcRegion"> The source region. </param>
/// <returns> 
/// - S_OK on success
/// - E_INVALIDSRC if imgSrc is not valid
/// - E_INVALIDDST if imgDst shares memory with imgSrc
/// - E_INVALIDARG if conversion from imgSrc to imgDst type is not supported
/// - E_OUTOFMEMORY if imgDst allocation was necessary and failed 
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - The Laplacian is an approximation of the second derivative of an image and is
/// implemented by applying the following 2D kernel:
/// \code
/// [ 0  1  0 ]
/// [ 1 -4  1 ]
/// [ 0  1  0 ]
/// \endcode
/// - If the provided output image does not have the same dimensions as 
/// imgSrc, then it will be recreated as a 1 banded float image with the 
/// dimensions of the source.
/// - The filter can not be done in place; the provided input and
/// output images must be different.
/// - Since the Laplacian filter is extremely sensitive to noise, you should
/// first filter the source with a Gaussian.

HRESULT VtComputeLaplacian(CFloatImg &imgDst, const CFloatImg &imgSrc, 
						   float fScale, const CRect *pSrcRegion = NULL);



/// \ingroup filter
/// <summary>Computes the 2D gradient of an image.</summary>
/// <param name="imgDst"> Output filtered image. </param>
/// <param name="imgSrc"> Input image. </param>
/// <param name="fScale"> The intensity of the output values will be 
/// scaled by this amount.</param>
/// <param name="pSrcRegion"> The source region. </param>
/// <returns> 
/// - S_OK on success
/// - E_INVALIDSRC if imgSrc is not valid
/// - E_INVALIDDST if imgDst shares memory with imgSrc
/// - E_INVALIDARG if conversion from imgSrc to imgDst type is not supported
/// - E_OUTOFMEMORY if imgDst allocation was necessary and failed 
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - The gradient is computed by applying the kernel `[-1 0 1]` twice: horizontally
/// to calculate the partial derivative with respect to x, and vertically to
/// calculate the partial derivative with respect to y.
/// - If the provided output image does not have the same dimensions as 
/// imgSrc, then it will be recreated as a CVec2Img image with the 
/// dimensions of the source.
/// - The filter can not be done in place; the provided input and
/// output images must be different.

HRESULT VtComputeGradient(CVec2Img &imgDst, const CFloatImg &imgSrc, 
						  float fScale, const CRect *pSrcRegion = NULL);

#ifndef VT_NO_XFORMS

/// \ingroup filtertransforms
/// <summary>Implementation of \ref IImageTransform that applies
/// a \ref VtComputeLaplacian() "Laplacian filter".</summary> 

class CLaplacianTransform:
    public CImageTransformUnaryGeo<CLaplacianTransform, false>
{
public:
	void    GetSrcPixFormat(IN OUT int* pfrmtSrcs, 
                            IN UINT  /*uSrcCnt*/,
                            IN int /*frmtDst*/)
	{ pfrmtSrcs[0] = VT_IMG_FIXED_PIXFRMT_MASK|VT_IMG_FIXED_ELFRMT_MASK|OBJ_LUMAFLOATIMG; }

	void    GetDstPixFormat(OUT int& frmtDst,
                            IN  const int* /*pfrmtSrcs*/, 
                            IN  UINT  /*uSrcCnt*/)
	{ frmtDst = VT_IMG_FIXED_PIXFRMT_MASK|VT_IMG_FIXED_ELFRMT_MASK|OBJ_LUMAFLOATIMG; }

    vt::CRect GetRequiredSrcRect(const vt::CRect& rctDst)
	{
        // Calculate size of source region needed.
        return vt::CRect(rctDst.left-1,  rctDst.top-1,
                         rctDst.right+1, rctDst.bottom+1);
	}

	vt::CRect GetAffectedDstRect(const vt::CRect& rctSrc)
	{
		// Calculate size of dest region affected by a change in src region
		return vt::CRect(rctSrc.left-1,  rctSrc.top-1,
						 rctSrc.right+1, rctSrc.bottom+1);
	}

	vt::CRect GetResultingDstRect(const vt::CRect& rctSrc)
	{
		return rctSrc;
	}

	HRESULT Transform(CImg* pimgDst, IN  const CRect& rctDst,
					  const CImg& imgSrc, const CPoint& ptSrc)
	{ 
		const CRect r = vt::CRect(rctDst-ptSrc);
		return VtComputeLaplacian((CFloatImg&)(*pimgDst), (CFloatImg&)imgSrc, 
			                      m_fScale, &r); 
	}

    HRESULT Clone(ITaskState **ppState)
	{ 
		CLaplacianTransform* pClone = VT_NOTHROWNEW CLaplacianTransform();
		HRESULT hr = CloneTaskState(ppState, pClone); 
		if( hr == S_OK )
		{
			ANALYZE_ASSUME( pClone != NULL );
			pClone->Initialize(m_fScale);
		}
		return hr;
	}

public:
    CLaplacianTransform() : m_fScale(1)
    {}

	/// <summary> Initialize the scale parameter of the filter. </summary> 
	/// <param name="fScale">The intensity of the output values will be
	/// scaled by this amount.  If Initialize is not called then the scale
	/// parameter will be 1.0 .</param> 
    void Initialize(float fScale)
    { m_fScale = fScale; }

protected:
    float m_fScale;
};

#endif

// smooths using edge preserving smoothing
// sigma is the standard deviation of the Gaussian kernel
// mu is the standard deviation of the Gaussian image value difference measure
// the difference measure is computed on image luminance
// alpha values are passed through but ignored during filtering
HRESULT VtBilateralFilter(CRGBAFloatImg &imgOut, const CRGBAFloatImg &imgIn, float fSigma, float fMu);

// this is the *single* banded float version
// if you call this with multiple bands it will fail
HRESULT VtBilateralFilter(CFloatImg &imgOut, const CFloatImg &imgIn, float fSigma, float fMu);

// same as above but with byte images. note that these dont have SSE speed up at present - swinder 10/1/07
HRESULT VtBilateralFilter(CRGBAImg &imgOut, const CRGBAImg &imgIn, float fSigma, float fMu);
HRESULT VtBilateralFilter(CByteImg &imgOut, const CByteImg &imgIn, float fSigma, float fMu);

};