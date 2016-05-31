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

#include "stdafx.h"

#include "vt_filter.h"
#include "vt_utils.h"
#include "vt_convert.h"

using namespace vt;

////////////////////////////////////////////////////////////////////////////////
// Filter Functions
////////////////////////////////////////////////////////////////////////////////
void ZeroExtendedBorders(int iWidth, int iHeight, const CRect &rctin, 
	CRect &rctout, CImg &imgdst, int &iDstX, int &iDstY)
{
	iDstX = 0;
	iDstY = 0;

	rctout = rctin;

	if(rctout.left==0)
	{
		rctout.left++;
		iDstX++;
		imgdst.Clear(CRect(0, 0, 1, rctin.Height())); // clear left edge
		if(rctout.IsRectEmpty())
			return;
	}
	if(rctout.right==iWidth)
	{
		// clear right edge
		rctout.right--;
		imgdst.Clear(CRect(rctin.Width()-1, 0, rctin.Width(), rctin.Height())); 
		if(rctout.IsRectEmpty())
			return;
	}
	if(rctout.top==0)
	{
		// clear top row
		rctout.top++;
		iDstY++;
		imgdst.Clear(CRect(iDstX, 0, iDstX + rctout.Width(), 1)); 
		if(rctout.IsRectEmpty())
			return;
	}
	if(rctout.bottom==iHeight)
	{  
		// clear bottom row
		rctout.bottom--;
		imgdst.Clear(CRect(iDstX, rctin.Height()-1, iDstX + rctout.Width(), 
			rctin.Height()));
	}
}


HRESULT vt::VtComputeLaplacian(CFloatImg &imgdst, const CFloatImg &imgsrc, 
	float fScale, const CRect *pRect)
{
	HRESULT hr = NOERROR;

	if(imgsrc.Bands()!=1)
		VT_HR_RET(E_NOTIMPL);

	CRect rctsrc = imgsrc.ClipRect(pRect);
	if(rctsrc.IsRectEmpty())
		return hr;
	if(imgdst.Width() < rctsrc.Width() || imgdst.Height() < rctsrc.Height())
		VT_HR_RET(E_INVALIDARG);

	int iDstX, iDstY;
	CRect rct;

	ZeroExtendedBorders(imgsrc.Width(), imgsrc.Height(), rctsrc, rct, imgdst, 
		iDstX, iDstY);
	if(rct.IsRectEmpty())
		return hr;

	int iSrcStrideBytes = (int)imgsrc.StrideBytes();
#if (defined(_M_IX86) || defined(_M_AMD64))
	int iDstStrideBytes = (int)imgdst.StrideBytes();
#endif

	int iProcessWidth = 256/sizeof(float);
	int iAlignmentWidth = 16/sizeof(float);
	int iProcessWidthMask = iProcessWidth - 1;
	int iAlignmentWidthMask = iAlignmentWidth - 1;

	int iSpan;
	int x;
	for(x = rct.left; x<rct.right; x += iSpan)
	{
		if(x & iAlignmentWidthMask) // not yet sse aligned
			iSpan = iAlignmentWidth - (x & iAlignmentWidthMask); // move to sse boundary
		else if(x & iProcessWidthMask)
			iSpan = iProcessWidth - (x & iProcessWidthMask);
		else
			iSpan = iProcessWidth;
		iSpan = VtMin((long)iSpan, rct.right - x); // don't go too far
		if(iSpan > iAlignmentWidth)
			iSpan = iSpan & (~iAlignmentWidthMask);  // keep aligned if possible (at end of row)

#if (defined(_M_IX86) || defined(_M_AMD64))
		if(g_SupportSSE2() && imgsrc.IsSSEAligned() && iSpan >= iAlignmentWidth && x+iSpan <= imgsrc.Width()-4)
		{
			__m128 mmneg4 = _mm_set1_ps(-4.0);
			__m128 mmscale = _mm_set1_ps(fScale);

			int y;
			const float *pfSrcRow = imgsrc.Ptr(x, rct.top);
			float *pfDstRow = imgdst.Ptr(x - rct.left + iDstX, iDstY);
			for(y = rct.top; y < rct.bottom; y++,
				pfSrcRow = PointerOffset(pfSrcRow, iSrcStrideBytes),
				pfDstRow = PointerOffset(pfDstRow, iDstStrideBytes) )
			{
				const float *pfSrc = pfSrcRow;
				float *pfDst = pfDstRow;
				__m128 mmleft = _mm_load_ss(pfSrc - 1); // initialize with pixel to the left

				int i;
				if(IsAligned16(pfDst))
				{
					for(i=0; i<iSpan; i += 4)
					{
						__m128 mmcenter = _mm_load_ps(pfSrc + i);
						__m128 mmsum = _mm_mul_ps(mmcenter, mmneg4);

						__m128 mmup = _mm_load_ps(PointerOffset(pfSrc + i, -iSrcStrideBytes));
						mmsum = _mm_add_ps(mmsum, mmup);

						__m128 mmdown = _mm_load_ps(PointerOffset(pfSrc + i, iSrcStrideBytes));
						mmsum = _mm_add_ps(mmsum, mmdown);

						__m128 mmnx1 = _mm_unpacklo_ps(mmleft, mmcenter);
						__m128 mmnx2 = _mm_shuffle_ps(mmnx1, mmcenter, _MM_SHUFFLE(2, 1, 1, 0));
						mmsum = _mm_add_ps(mmsum, mmnx2);

						mmleft = _mm_shuffle_ps(mmcenter, mmcenter, _MM_SHUFFLE(3, 3, 3, 3));

						__m128 mmright = _mm_load_ps(pfSrc + i + 4);
						__m128 mmpx1 = _mm_move_ss(mmcenter, mmright);
						__m128 mmpx2 = _mm_shuffle_ps(mmpx1, mmpx1, _MM_SHUFFLE(0, 3, 2, 1));
						mmsum = _mm_add_ps(mmsum, mmpx2);

						mmsum = _mm_mul_ps(mmsum, mmscale);
						_mm_store_ps(pfDst + i, mmsum);
					}
				}
				else
				{
					for(i=0; i<iSpan; i += 4)
					{
						__m128 mmcenter = _mm_load_ps(pfSrc + i);
						__m128 mmsum = _mm_mul_ps(mmcenter, mmneg4);

						__m128 mmup = _mm_load_ps(PointerOffset(pfSrc + i, -iSrcStrideBytes));
						mmsum = _mm_add_ps(mmsum, mmup);

						__m128 mmdown = _mm_load_ps(PointerOffset(pfSrc + i, iSrcStrideBytes));
						mmsum = _mm_add_ps(mmsum, mmdown);

						__m128 mmnx1 = _mm_unpacklo_ps(mmleft, mmcenter);
						__m128 mmnx2 = _mm_shuffle_ps(mmnx1, mmcenter, _MM_SHUFFLE(2, 1, 1, 0));
						mmsum = _mm_add_ps(mmsum, mmnx2);

						mmleft = _mm_shuffle_ps(mmcenter, mmcenter, _MM_SHUFFLE(3, 3, 3, 3));

						__m128 mmright = _mm_load_ps(pfSrc + i + 4);
						__m128 mmpx1 = _mm_move_ss(mmcenter, mmright);
						__m128 mmpx2 = _mm_shuffle_ps(mmpx1, mmpx1, _MM_SHUFFLE(0, 3, 2, 1));
						mmsum = _mm_add_ps(mmsum, mmpx2);

						mmsum = _mm_mul_ps(mmsum, mmscale);
						_mm_storeu_ps(pfDst + i, mmsum);
					}
				}
			}
		}
		else
#endif
		{
			int y;
			for(y = rct.top; y < rct.bottom; y++)
			{
				const float *pfSrc = imgsrc.Ptr(x, y);
				float *pfDst = imgdst.Ptr(x - rct.left + iDstX, y - rct.top + iDstY);
				int i;
				for(i=0; i<iSpan; i++, pfSrc++)
				{
					pfDst[i] = fScale * (pfSrc[-1] + pfSrc[+1] 
					+ *PointerOffset(pfSrc, -iSrcStrideBytes) + *PointerOffset(pfSrc, iSrcStrideBytes)
						- 4 * pfSrc[0]);
				}
			}
		}
	}


	return NOERROR;
}

HRESULT vt::VtComputeGradient(CVec2Img &imgdst, const CFloatImg &imgsrc, 
	float fScale, const CRect *pRect)
{
	HRESULT hr = NOERROR;

	if(imgsrc.Bands()!=1)
		VT_HR_RET(E_NOTIMPL);

	CRect rctsrc = imgsrc.ClipRect(pRect);
	if(rctsrc.IsRectEmpty())
		return hr;
	if(imgdst.Width() < rctsrc.Width() || imgdst.Height() < rctsrc.Height())
		VT_HR_RET(E_INVALIDARG);

	int iDstX, iDstY;
	CRect rct;

	ZeroExtendedBorders(imgsrc.Width(), imgsrc.Height(), rctsrc, rct, imgdst, iDstX, iDstY);
	if(rct.IsRectEmpty())
		return hr;

	int iSrcStrideBytes = (int)imgsrc.StrideBytes();
#if (defined(_M_IX86) || defined(_M_AMD64))
	int iDstStrideBytes = (int)imgdst.StrideBytes();
#endif

	int iProcessWidth = 256/sizeof(float);
	int iAlignmentWidth = 16/sizeof(float);
	int iProcessWidthMask = iProcessWidth - 1;
	int iAlignmentWidthMask = iAlignmentWidth - 1;

	int iSpan;
	int x;
	for(x = rct.left; x<rct.right; x += iSpan)
	{
		if(x & iAlignmentWidthMask) // not yet sse aligned
			iSpan = iAlignmentWidth - (x & iAlignmentWidthMask); // move to sse boundary
		else if(x & iProcessWidthMask)
			iSpan = iProcessWidth - (x & iProcessWidthMask);
		else
			iSpan = iProcessWidth;
		iSpan = VtMin((long)iSpan, rct.right - x); // don't go too far
		if(iSpan > iAlignmentWidth)
			iSpan = iSpan & (~iAlignmentWidthMask);  // keep aligned if possible (at end of row)

#if (defined(_M_IX86) || defined(_M_AMD64))
		if(g_SupportSSE2() && imgsrc.IsSSEAligned() && iSpan >= iAlignmentWidth && x+iSpan <= imgsrc.Width()-4)
		{
			__m128 mmscale = _mm_set1_ps(fScale);

			int y;
			const float *pfSrcRow = imgsrc.Ptr(x, rct.top);
			CVec2f *pvDstRow = imgdst.Ptr(x - rct.left + iDstX, iDstY);
			for(y = rct.top; y < rct.bottom; y++,
				pfSrcRow = PointerOffset(pfSrcRow, iSrcStrideBytes),
				pvDstRow = PointerOffset(pvDstRow, iDstStrideBytes) )
			{
				const float *pfSrc = pfSrcRow;
				CVec2f *pvDst = pvDstRow;
				__m128 mmleft = _mm_load_ss(pfSrc - 1); // initialize with pixel to the left

				int i;
				if(IsAligned16(pvDstRow))
				{
					for(i=0; i<iSpan; i += 4, pvDst += 4)
					{
						__m128 mm1 = _mm_load_ps(PointerOffset(pfSrc + i, -iSrcStrideBytes));
						__m128 mm3 = _mm_load_ps(PointerOffset(pfSrc + i, iSrcStrideBytes));
						__m128 mmdy = _mm_sub_ps(mm3, mm1);

						__m128 mm2 = _mm_load_ps(pfSrc + i);
						__m128 mmright = _mm_load_ps(pfSrc + i + 4);

						__m128 mmnx1 = _mm_unpacklo_ps(mmleft, mm2);
						__m128 mmnx2 = _mm_shuffle_ps(mmnx1, mm2, _MM_SHUFFLE(2, 1, 1, 0));
						__m128 mmpx1 = _mm_move_ss(mm2, mmright);
						__m128 mmpx2 = _mm_shuffle_ps(mmpx1, mmpx1, _MM_SHUFFLE(0, 3, 2, 1));
						__m128 mmdx = _mm_sub_ps(mmpx2, mmnx2);

						mmleft = _mm_shuffle_ps(mm2, mm2, _MM_SHUFFLE(3, 3, 3, 3));

						mmdx = _mm_mul_ps(mmdx, mmscale);
						mmdy = _mm_mul_ps(mmdy, mmscale);

						__m128 mms0 = _mm_unpacklo_ps(mmdx, mmdy);
						_mm_store_ps((float *)pvDst, mms0);
						__m128 mms1 = _mm_unpackhi_ps(mmdx, mmdy);
						_mm_store_ps((float *)(pvDst + 2), mms1);
					}
				}
				else
				{
					for(i=0; i<iSpan; i += 4, pvDst += 4)
					{
						__m128 mm1 = _mm_load_ps(PointerOffset(pfSrc + i, -iSrcStrideBytes));
						__m128 mm3 = _mm_load_ps(PointerOffset(pfSrc + i, iSrcStrideBytes));
						__m128 mmdy = _mm_sub_ps(mm3, mm1);

						__m128 mm2 = _mm_load_ps(pfSrc + i);
						__m128 mmright = _mm_load_ps(pfSrc + i + 4);

						__m128 mmnx1 = _mm_unpacklo_ps(mmleft, mm2);
						__m128 mmnx2 = _mm_shuffle_ps(mmnx1, mm2, _MM_SHUFFLE(2, 1, 1, 0));
						__m128 mmpx1 = _mm_move_ss(mm2, mmright);
						__m128 mmpx2 = _mm_shuffle_ps(mmpx1, mmpx1, _MM_SHUFFLE(0, 3, 2, 1));
						__m128 mmdx = _mm_sub_ps(mmpx2, mmnx2);

						mmleft = _mm_shuffle_ps(mm2, mm2, _MM_SHUFFLE(3, 3, 3, 3));

						mmdx = _mm_mul_ps(mmdx, mmscale);
						mmdy = _mm_mul_ps(mmdy, mmscale);

						__m128 mms0 = _mm_unpacklo_ps(mmdx, mmdy);
						_mm_storeu_ps((float *)pvDst, mms0);
						__m128 mms1 = _mm_unpackhi_ps(mmdx, mmdy);
						_mm_storeu_ps((float *)(pvDst + 2), mms1);
					}
				}
			}
		}
		else
#endif
		{
			int y;
			for(y = rct.top; y < rct.bottom; y++)
			{
				const float *pfSrc = imgsrc.Ptr(x, y);
				CVec2f *pvDst = imgdst.Ptr(x - rct.left + iDstX, y - rct.top + iDstY);
				int i;
				for(i=0; i<iSpan; i++, pfSrc++)
				{
					pvDst[i].x = fScale * (pfSrc[1] - pfSrc[-1]);
					pvDst[i].y = fScale * (*PointerOffset(pfSrc, iSrcStrideBytes) 
						- *PointerOffset(pfSrc, -iSrcStrideBytes));
				}
			}
		}
	}


	return NOERROR;
}

// special SSE four band version
// this version uses luma information as the filter difference
HRESULT vt::VtBilateralFilter(CRGBAFloatImg &imgOut, const CRGBAFloatImg &imgIn, 
	float fSigma, float fMu)
{
	VT_HR_BEGIN();

	int x,y,w,h;
	w = imgIn.Width();
	h = imgIn.Height();
	CFloatImg imgLuma;
	C1dKernel kern;

	if(!imgIn.IsValid() || !imgOut.IsValid() || imgOut.Width()!=w || imgOut.Height()!=h)
		VT_HR_EXIT( E_INVALIDARG );
	if(fSigma<=0 || fMu<0)
		VT_HR_EXIT( E_INVALIDARG );

	VT_HR_EXIT( Create1dGaussKernel(kern, fSigma) );

	int ks = kern.Center();
	int kw = kern.Width();

	float rgfExp[256];
	int i;
	for(i=0;i<256;i++)
		rgfExp[i] = exp(-i*i/(2.0f*75.0f*75.0f));
	float fDiffScale = 75.0f/fMu;

	VT_HR_EXIT( imgLuma.Create(w,h) );

	for(y=0; y<h; y++)
	{
		float *pLuma = imgLuma.Ptr(y);
		const RGBAFloatPix *pIn = imgIn.Ptr(y);
		for(x=0; x<w; x++, pLuma++, pIn++)
			*pLuma = fDiffScale * VtLumaFromRGB_CCIR601YPbPr(pIn);
	}

	for(y=0;y<h; y++)
	{
		for(x=0;x<w;x++)
		{
			RGBAFloatPix *pOutVal = imgOut.Ptr(x,y);
			float fCenterVal = imgLuma(x,y);
			int yy;
			int ky = 0;
			int kxs = VtMax(0, ks-x);
			int kxe = VtMin(kw, ks+w-x);
			int xoffset = x-ks+kxs;

#if (defined(_M_IX86) || defined(_M_AMD64))
			if(g_SupportSSE1())
			{
				__m128 xmSum = _mm_setzero_ps();
				__m128 xmOne = _mm_set_ss(1.0f);

				for(yy=-ks; yy<=ks; yy++, ky++)
				{
					int yc = yy+y;
					if(yc<0 || yc>h-1)
						continue;
					float fKY = kern[ky];
					int kx;
					const RGBAFloatPix *pOffsetVal = imgIn.Ptr(xoffset,yc);
					float *pOffsetLuma = imgLuma.Ptr(xoffset,yc);

					for(kx=kxs; kx<kxe; kx++, pOffsetVal++, pOffsetLuma++)
					{
						float fDiff = fabs(*pOffsetLuma - fCenterVal);
						if(fDiff<255)
						{
							float fWeight = kern[kx] * fKY * rgfExp[F2I(fDiff)];
							__m128 xmWeight = _mm_set1_ps(fWeight);
							__m128 xmARGB = _mm_loadu_ps((float *)pOffsetVal);
							__m128 xmRGBX = _mm_shuffle_ps(xmARGB, xmARGB, _MM_SHUFFLE(2, 1, 0, 0));
							xmRGBX = _mm_move_ss(xmRGBX, xmOne);
							xmWeight = _mm_mul_ps(xmWeight, xmRGBX);
							xmSum = _mm_add_ps(xmSum, xmWeight);
						}
					}
				}

				const float *pfInA = &imgIn(x,y).a;
				__m128 xmWeightSum = _mm_shuffle_ps(xmSum, xmSum, _MM_SHUFFLE(0, 0, 0, 0));
				xmSum = _mm_div_ps(xmSum, xmWeightSum);
				__m128 xmA = _mm_load_ss(pfInA);
				xmSum = _mm_move_ss(xmSum, xmA);
				xmSum = _mm_shuffle_ps(xmSum, xmSum, _MM_SHUFFLE(0, 3, 2, 1));
				_mm_storeu_ps((float *)pOutVal, xmSum);
			}
			else
#endif
			{
				float fSumR = 0;
				float fSumG = 0;
				float fSumB = 0;
				float fWeightSum = 0;

				for(yy=-ks; yy<=ks; yy++, ky++)
				{
					int yc = yy+y;
					if(yc<0 || yc>h-1)
						continue;
					float fKY = kern[ky];
					int kx;
					const RGBAFloatPix *pOffsetVal = imgIn.Ptr(xoffset,yc);
					float *pOffsetLuma = imgLuma.Ptr(xoffset,yc);

					for(kx=kxs; kx<kxe; kx++, pOffsetVal++, pOffsetLuma++)
					{
						float fDiff = fabs(*pOffsetLuma - fCenterVal);
						if(fDiff<255)
						{
							int ival = F2I(fDiff); // does round to nearest
							float fWeight = kern[kx] * fKY * rgfExp[ival];
							fSumR += fWeight * pOffsetVal->r;
							fSumG += fWeight * pOffsetVal->g;
							fSumB += fWeight * pOffsetVal->b;
							fWeightSum += fWeight;
						}
					}
				}

				fWeightSum = 1/fWeightSum;
				pOutVal->r = fSumR * fWeightSum;
				pOutVal->g = fSumG * fWeightSum;
				pOutVal->b = fSumB * fWeightSum;
				pOutVal->a = imgIn(x,y).a;

			}
		}
	}

	VT_HR_END();
}


HRESULT vt::VtBilateralFilter(CRGBAImg &imgOut, const CRGBAImg &imgIn, float fSigma, float fMu)
{
    VT_HR_BEGIN();

    int x,y,w,h;
    w = imgIn.Width();
    h = imgIn.Height();
    CFloatImg imgLuma;
    C1dKernel kern;

    if(!imgIn.IsValid() || !imgOut.IsValid() || imgOut.Width()!=w || imgOut.Height()!=h)
        VT_HR_EXIT( E_INVALIDARG );
    if(fSigma<=0 || fMu<0)
        VT_HR_EXIT( E_INVALIDARG );

    VT_HR_EXIT( Create1dGaussKernel(kern, fSigma) );

    int ks = kern.Center();
    int kw = kern.Width();

    float rgfExp[256];
    int i;
    for(i=0;i<256;i++)
        rgfExp[i] = exp(-i*i/(2.0f*75.0f*75.0f));
    float fDiffScale = 75.0f/fMu;

    VT_HR_EXIT( imgLuma.Create(w,h) );

    for(y=0; y<h; y++)
    {
        float *pLuma = imgLuma.Ptr(y);
        const RGBAPix *pIn = imgIn.Ptr(y);
        for(x=0; x<w; x++, pLuma++, pIn++)
            *pLuma = fDiffScale * VtLumaFromRGB_CCIR601YPbPr(pIn);
    }

    for(y=0;y<h; y++)
	{
        for(x=0;x<w;x++)
        {
            RGBAPix *pOutVal = imgOut.Ptr(x,y);
            float fCenterVal = imgLuma(x,y);
            int yy;
            int ky = 0;
            int kxs = VtMax(0, ks-x);
            int kxe = VtMin(kw, ks+w-x);
            int xoffset = x-ks+kxs;

            if(g_SupportSSE1())
            {
                /*
                __m128 xmSum = _mm_setzero_ps();
                __m128 xmOne = _mm_set_ss(1.0f);

                for(yy=-ks; yy<=ks; yy++, ky++)
                {
                int yc = yy+y;
                if(yc<0 || yc>h-1)
                continue;
                float fKY = kern[ky];
                int kx;
                RGBAFloatPix *pOffsetVal = imgIn.Ptr(xoffset,yc);
                float *pOffsetLuma = imgLuma.Ptr(xoffset,yc);

                for(kx=kxs; kx<kxe; kx++, pOffsetVal++, pOffsetLuma++)
                {
                float fDiff = fabs(*pOffsetLuma - fCenterVal);
                if(fDiff<255)
                {
                float fWeight = kern[kx] * fKY * rgfExp[F2I(fDiff)];
                __m128 xmWeight = _mm_set1_ps(fWeight);
                __m128 xmARGB = _mm_loadu_ps((float *)pOffsetVal);
                __m128 xmRGBX = _mm_shuffle_ps(xmARGB, xmARGB, _MM_SHUFFLE(2, 1, 0, 0));
                xmRGBX = _mm_move_ss(xmRGBX, xmOne);
                xmWeight = _mm_mul_ps(xmWeight, xmRGBX);
                xmSum = _mm_add_ps(xmSum, xmWeight);
                }
                }
                }

                float *pfInA = &imgIn(x,y).a;
                __m128 xmWeightSum = _mm_shuffle_ps(xmSum, xmSum, _MM_SHUFFLE(0, 0, 0, 0));
                xmSum = _mm_div_ps(xmSum, xmWeightSum);
                __m128 xmA = _mm_load_ss(pfInA);
                xmSum = _mm_move_ss(xmSum, xmA);
                xmSum = _mm_shuffle_ps(xmSum, xmSum, _MM_SHUFFLE(0, 3, 2, 1));
                _mm_storeu_ps((float *)pOutVal, xmSum);
                */
            }
            else
            {
                float fSumR = 0;
                float fSumG = 0;
                float fSumB = 0;
                float fWeightSum = 0;

                for(yy=-ks; yy<=ks; yy++, ky++)
                {
                    int yc = yy+y;
                    if(yc<0 || yc>h-1)
                        continue;
                    float fKY = kern[ky];
                    int kx;
                    const RGBAPix *pOffsetVal = imgIn.Ptr(xoffset,yc);
                    float *pOffsetLuma = imgLuma.Ptr(xoffset,yc);

                    for(kx=kxs; kx<kxe; kx++, pOffsetVal++, pOffsetLuma++)
                    {
                        float fDiff = fabs(*pOffsetLuma - fCenterVal);
                        if(fDiff<255)
                        {
                            int ival = F2I(fDiff); // does round to nearest
                            float fWeight = kern[kx] * fKY * rgfExp[ival];
                            fSumR += fWeight * pOffsetVal->r;
                            fSumG += fWeight * pOffsetVal->g;
                            fSumB += fWeight * pOffsetVal->b;
                            fWeightSum += fWeight;
                        }
                    }
                }

                fWeightSum = 1/fWeightSum;
                VtClip(&pOutVal->r, fSumR * fWeightSum);
                VtClip(&pOutVal->g, fSumG * fWeightSum);
                VtClip(&pOutVal->b,  fSumB * fWeightSum);
                pOutVal->a = imgIn(x,y).a;

            }
        }
	}

	VT_HR_END();
}

// one band version

HRESULT vt::VtBilateralFilter(CFloatImg &imgOut, const CFloatImg &imgIn, float fSigma, float fMu)
{
    VT_HR_BEGIN();

    int x,y,w,h;
    w = imgIn.Width();
    h = imgIn.Height();
    C1dKernel kern;

    if( (imgIn.Bands( ) != 1) || (imgOut.Bands( ) != 1) )
    {
        VT_HR_EXIT( E_INVALIDARG);
    }

    if(!imgIn.IsValid() || !imgOut.IsValid() || imgOut.Width()!=w || imgOut.Height()!=h)
        VT_HR_EXIT( E_INVALIDARG );
    if(fSigma<=0 || fMu<0)
        VT_HR_EXIT( E_INVALIDARG );

    VT_HR_EXIT( Create1dGaussKernel(kern, fSigma) );

    int ks = kern.Center();
    int kw = kern.Width();

    float rgfExp[256];
    int i;
    for(i=0;i<256;i++)
        rgfExp[i] = exp(-i*i/(2.0f*75.0f*75.0f));
    float fDiffScale = 75.0f/fMu;
    float fDiffThresh = 255/fDiffScale;

    for(y=0;y<h; y++)
	{
        for(x=0;x<w;x++)
        {
            float fCenterVal = imgIn(x,y);
            int yy;
            int ky = 0;
            int kxs = VtMax(0, ks-x);
            int kxe = VtMin(kw, ks+w-x);
            int xoffset = x-ks+kxs;

            float fSum = 0;
            float fWeightSum = 0;
            for(yy=-ks; yy<=ks; yy++, ky++)
            {
                int yc = yy+y;
                if(yc<0 || yc>h-1)
                    continue;
                float fKY = kern[ky];
                int kx;
                const float *pOffsetVal = imgIn.Ptr(xoffset, yc);

                for(kx=kxs; kx<kxe; kx++, pOffsetVal++)
                {
                    float fDiff = fabs(*pOffsetVal - fCenterVal);
                    if(fDiff<fDiffThresh)
                    {
                        float fWeight = kern[kx] * fKY * rgfExp[F2I(fDiffScale * fDiff)];
                        fSum += fWeight * *pOffsetVal;
                        fWeightSum += fWeight;
                    }
                }
            }

            imgOut(x,y) = fSum / fWeightSum;
        }
	}

	VT_HR_END();
}

HRESULT vt::VtBilateralFilter(CByteImg &imgOut, const CByteImg &imgIn, 
                              float fSigma, float fMu)
{
    VT_HR_BEGIN();

    int x,y,w,h;
    w = imgIn.Width();
    h = imgIn.Height();
    C1dKernel kern;

    if( (imgIn.Bands( ) != 1) || (imgOut.Bands( ) != 1) )
    {
        VT_HR_EXIT( E_INVALIDARG);
    }

    if(!imgIn.IsValid() || !imgOut.IsValid() || imgOut.Width()!=w || imgOut.Height()!=h)
        VT_HR_EXIT( E_INVALIDARG );
    if(fSigma<=0 || fMu<0)
        VT_HR_EXIT( E_INVALIDARG );

    VT_HR_EXIT( Create1dGaussKernel( kern, fSigma) );

    int ks = kern.Center();
    int kw = kern.Width();

    float rgfExp[256];
    int i;
    for(i=0;i<256;i++)
        rgfExp[i] = exp(-i*i/(2.0f*75.0f*75.0f));
    float fDiffScale = 75.0f/fMu;
    float fDiffThresh = 255/fDiffScale;

    for(y=0;y<h; y++)
	{
        for(x=0;x<w;x++)
        {
            int iCenterVal = imgIn(x,y);
            int yy;
            int ky = 0;
            int kxs = VtMax(0, ks-x);
            int kxe = VtMin(kw, ks+w-x);
            int xoffset = x-ks+kxs;

            float fSum = 0;
            float fWeightSum = 0;
            for(yy=-ks; yy<=ks; yy++, ky++)
            {
                int yc = yy+y;
                if(yc<0 || yc>h-1)
                    continue;
                float fKY = kern[ky];
                int kx;
				const Byte *pOffsetVal = imgIn.Ptr(xoffset, yc);

                for(kx=kxs; kx<kxe; kx++, pOffsetVal++)
                {
                    float fDiff = (float)abs(*pOffsetVal - iCenterVal);
                    if(fDiff<fDiffThresh)
                    {
                        float fWeight = kern[kx] * fKY * rgfExp[F2I(fDiffScale * fDiff)];
                        fSum += fWeight * *pOffsetVal;
                        fWeightSum += fWeight;
                    }
                }
            }

            VtClip(imgOut.Ptr(x,y), fSum / fWeightSum);
        }
	}

	VT_HR_END();
}
