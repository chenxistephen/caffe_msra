//+---------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation 2009.  All rights reserved.
//
//  Description:
//     Image Pyramid Filtering
//
//  History:
//      2011/7/26-sbaker
//			Created 
// 
//----------------------------------------------------------------------------

#include "stdafx.h"
#include "pyrfilter.h"

HRESULT FilterHoriz121NoSSE(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	HRESULT hr = S_OK;
	
	int iSrcWidth = imgSrc.Width();
	int iSrcHeight = imgSrc.Height();
	VT_HR_EXIT((iDstWidth < 2) ? E_FAIL : S_OK);
	VT_HR_EXIT((imgDst.Height() < iDstHeight || imgDst.Width() < iDstWidth) ? E_FAIL : S_OK);
	VT_HR_EXIT((iSrcHeight < iDstHeight || iSrcWidth < iDstWidth) ? E_FAIL : S_OK);

	for(int iY=0; iY<iDstHeight; iY++)
	{
		float *pfSrcRow = imgSrc.Ptr(iY);
		float *pfDstRow = imgDst.Ptr(iY);
		float *pfDstRowEnd = pfDstRow+iDstWidth-1;
		float fPrev = *pfSrcRow++;
		float fCurr = *pfSrcRow++;
		*pfDstRow++ = 0.25f*(3.0f*fPrev+fCurr);
		while(pfDstRow<pfDstRowEnd)
		{
			float fNext = *pfSrcRow++;
			*pfDstRow++ = 0.25f*(fPrev+2.0f*fCurr+fNext);
			fPrev = fCurr;
			fCurr = fNext;
		}
		*pfDstRow++ = 0.25f*(fPrev+3.0f*fCurr);
	}

Exit:
	return hr;
}

HRESULT FilterHoriz121SSE2(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	HRESULT hr = S_OK;

#if (defined(_M_IX86) || defined(_M_AMD64))
	int iSrcWidth = imgSrc.Width();
	int iSrcHeight = imgSrc.Height();
	__m128 m128Two = _mm_set_ps1(2.0f);
	__m128 m128Quarter = _mm_set_ps1(0.25f);
    __m128 m128One = _mm_set1_ps(1.0f);
	VT_HR_EXIT((iDstWidth < 2) ? E_FAIL : S_OK);
	VT_HR_EXIT((imgDst.Height() < iDstHeight || imgDst.Width() < iDstWidth) ? E_FAIL : S_OK);
	VT_HR_EXIT((iSrcHeight < iDstHeight || iSrcWidth < iDstWidth) ? E_FAIL : S_OK);

	for(int iY=0; iY<iDstHeight; iY++)
	{
		float *pfSrcRow = imgSrc.Ptr(iY);
		float *pfDstRow = imgDst.Ptr(iY);
		int iWidth1 = VtMin(iDstWidth-1, 4);
		int iWidth2 = 4*((iDstWidth-1)/4);
		int iWidth3 = iDstWidth-1;
		float *pfDstRowEnd1 = pfDstRow+iWidth1;
		float *pfDstRowEnd2 = pfDstRow+iWidth2;
		float *pfDstRowEnd3 = pfDstRow+iWidth3;
		float fPrev = *pfSrcRow++;
		float fCurr = *pfSrcRow++;
		*pfDstRow++ = 0.25f*(3.0f*fPrev+fCurr);
		while(pfDstRow<pfDstRowEnd1)
		{
			float fNext = *pfSrcRow++;
			*pfDstRow++ = 0.25f*(fPrev+2.0f*fCurr+fNext);
			fPrev = fCurr;
			fCurr = fNext;
		}
		if (pfDstRow<pfDstRowEnd2)
		{
			__m128 m128PrevX = _mm_set_ps1(fPrev);
			pfSrcRow--;
			__m128 m128CurrentX = _mm_load_ps(pfSrcRow);
			while(pfDstRow<pfDstRowEnd2)
			{
				// float fNext = *pfSrcRow++;
				// *pfDstRow++ = 0.25f*(fPrev+2.0f*fCurr+fNext);
				// fPrev = fCurr;
				// fCurr = fNext;

				pfSrcRow += 4;
				__m128 m128NextX = _mm_load_ps(pfSrcRow);

				// Shuffle for Left
				__m128 m128LXT1 = _mm_unpacklo_ps(m128PrevX, m128CurrentX);
				__m128 m128LXT2 = _mm_shuffle_ps(m128LXT1, m128CurrentX, _MM_SHUFFLE(2, 1, 1, 0));

				// Shuffle for Right
				__m128 m128RXT1 = _mm_move_ss(m128CurrentX, m128NextX);
				__m128 m128RXT2 = _mm_shuffle_ps(m128RXT1, m128RXT1, _MM_SHUFFLE(0, 3, 2, 1));

				// Do the math
				__m128 m128SumX = _mm_mul_ps(m128Two, m128CurrentX);
				m128SumX = _mm_add_ps(m128SumX, m128LXT2);
				m128SumX = _mm_add_ps(m128SumX, m128RXT2);
				m128SumX = _mm_mul_ps(m128Quarter, m128SumX);

				// Store
				_mm_storeu_ps(pfDstRow, m128SumX);
				pfDstRow += 4;

				// Update Left and Current
				m128PrevX = _mm_shuffle_ps(m128CurrentX, m128CurrentX, _MM_SHUFFLE(3, 3, 3, 3));
				m128CurrentX = _mm_mul_ps(m128NextX, m128One);
			}
			pfSrcRow--;
			fPrev = *pfSrcRow++;
			fCurr = *pfSrcRow++;
		}
		while(pfDstRow<pfDstRowEnd3)
		{
			float fNext = *pfSrcRow++;
			*pfDstRow++ = 0.25f*(fPrev+2.0f*fCurr+fNext);
			fPrev = fCurr;
			fCurr = fNext;
		}
		*pfDstRow++ = 0.25f*(fPrev+3.0f*fCurr);
	}
Exit:
#else
    imgSrc,imgDst,iDstWidth,iDstHeight;
#endif
	return hr;
}

HRESULT FilterHoriz121SubsNoSSE(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	HRESULT hr = S_OK;
	
	int iSrcWidth = imgSrc.Width();
	int iSrcHeight = imgSrc.Height();	
	VT_HR_EXIT((iDstWidth < 1) ? E_FAIL : S_OK);
	VT_HR_EXIT((imgDst.Height() < iDstHeight || imgDst.Width() < iDstWidth) ? E_FAIL : S_OK);
	VT_HR_EXIT((iSrcHeight < iDstHeight || iSrcWidth < 2*iDstWidth) ? E_FAIL : S_OK);

	for(int iY=0; iY<iDstHeight; iY++)
	{
		float *pfSrcRow = imgSrc.Ptr(iY);
		float *pfDstRow = imgDst.Ptr(iY);
		float *pfDstRowEnd = pfDstRow+iDstWidth;
		float fCurr = *pfSrcRow++;
		float fNext = *pfSrcRow++;
		*pfDstRow++ = 0.25f*(3.0f*fCurr+fNext);
		while(pfDstRow<pfDstRowEnd)
		{
			float fPrev = fNext;
			fCurr = *pfSrcRow++;
			fNext = *pfSrcRow++;
			*pfDstRow++ = 0.25f*(fPrev+2.0f*fCurr+fNext);
		}
	}

Exit:
	return hr;
}

HRESULT FilterHoriz121SubsSSE2(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	HRESULT hr = S_OK;
	
#if (defined(_M_IX86) || defined(_M_AMD64))
	int iSrcWidth = imgSrc.Width();
	int iSrcHeight = imgSrc.Height();
	__m128 m128Two = _mm_set_ps1(2.0f);
	__m128 m128Quarter = _mm_set_ps1(0.25f);
    __m128 m128One = _mm_set1_ps(1.0f);
	VT_HR_EXIT((iDstWidth < 1) ? E_FAIL : S_OK);
	VT_HR_EXIT((imgDst.Height() < iDstHeight || imgDst.Width() < iDstWidth) ? E_FAIL : S_OK);
	VT_HR_EXIT((iSrcHeight < iDstHeight || iSrcWidth < 2*iDstWidth) ? E_FAIL : S_OK);

	for(int iY=0; iY<iDstHeight; iY++)
	{
		float *pfSrcRow = imgSrc.Ptr(iY);
		float *pfDstRow = imgDst.Ptr(iY);
		int iWidth1 = VtMin(iDstWidth-1, 4);
		int iWidth2 = 4*((iDstWidth-1)/4);
		int iWidth3 = iDstWidth;
		float *pfDstRowEnd1 = pfDstRow+iWidth1;
		float *pfDstRowEnd2 = pfDstRow+iWidth2;
		float *pfDstRowEnd3 = pfDstRow+iWidth3;
		float fCurr = *pfSrcRow++;
		float fNext = *pfSrcRow++;
		*pfDstRow++ = 0.25f*(3.0f*fCurr+fNext);
		while(pfDstRow<pfDstRowEnd1)
		{
			float fPrev = fNext;
			fCurr = *pfSrcRow++;
			fNext = *pfSrcRow++;
			*pfDstRow++ = 0.25f*(fPrev+2.0f*fCurr+fNext);
		}
		if (pfDstRow<pfDstRowEnd2)
		{
			__m128 m128PrevX = _mm_set_ps1(fNext);
			__m128 m128CurrentX = _mm_load_ps(pfSrcRow);
			while(pfDstRow<pfDstRowEnd2)
			{
				// float fNext = *pfSrcRow++;
				// *pfDstRow++ = 0.25f*(fPrev+2.0f*fCurr+fNext);
				// fPrev = fCurr;
				// fCurr = fNext;

				pfSrcRow += 4;
				__m128 m128NextX = _mm_load_ps(pfSrcRow);

				// Shuffle for Left
				__m128 m128LXT1 = _mm_unpacklo_ps(m128PrevX, m128CurrentX);
				__m128 m128LXT2 = _mm_shuffle_ps(m128LXT1, m128CurrentX, _MM_SHUFFLE(2, 1, 1, 0));

				// Shuffle for Right
				__m128 m128RXT1 = _mm_move_ss(m128CurrentX, m128NextX);
				__m128 m128RXT2 = _mm_shuffle_ps(m128RXT1, m128RXT1, _MM_SHUFFLE(0, 3, 2, 1));

				// Do the math
				__m128 m128SumX = _mm_mul_ps(m128Two, m128CurrentX);
				m128SumX = _mm_add_ps(m128SumX, m128LXT2);
				m128SumX = _mm_add_ps(m128SumX, m128RXT2);
				m128SumX = _mm_mul_ps(m128Quarter, m128SumX);

				// Store
#ifndef VT_GCC
				*pfDstRow++ = m128SumX.m128_f32[0];
				*pfDstRow++ = m128SumX.m128_f32[2];
#else
				*pfDstRow++ = SSE2_mm_extract_epf32(m128SumX, 0);
				*pfDstRow++ = SSE2_mm_extract_epf32(m128SumX, 2);
#endif

				// Update Left and Current
				m128PrevX = _mm_shuffle_ps(m128CurrentX, m128CurrentX, _MM_SHUFFLE(3, 3, 3, 3));
				m128CurrentX = _mm_mul_ps(m128NextX, m128One);
			}
			pfSrcRow--;
			fNext = *pfSrcRow++;
		}
		while(pfDstRow<pfDstRowEnd3)
		{
			float fPrev = fNext;
			fCurr = *pfSrcRow++;
			fNext = *pfSrcRow++;
			*pfDstRow++ = 0.25f*(fPrev+2.0f*fCurr+fNext);
		}
	}

Exit:
#else
    imgSrc,imgDst,iDstWidth,iDstHeight;
#endif
	return hr;
}

HRESULT FilterVert121NoSSE(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	HRESULT hr = S_OK;

	int iSrcWidth = imgSrc.Width();
	int iSrcHeight = imgSrc.Height();
	VT_HR_EXIT((iDstHeight < 2) ? E_FAIL : S_OK);
	VT_HR_EXIT((imgDst.Height() < iDstHeight || imgDst.Width() < iDstWidth) ? E_FAIL : S_OK);
	VT_HR_EXIT((iSrcHeight < iDstHeight || iSrcWidth < iDstWidth) ? E_FAIL : S_OK);

	for(int iY=0; iY<iDstHeight; iY++)
	{
		float *pfSrcRow = imgSrc.Ptr(iY);
		float *pfSrcAboveRow = NULL;
		float *pfSrcBelowRow = NULL;
		if (iY==0)
		{
			pfSrcAboveRow = pfSrcRow;
			pfSrcBelowRow = imgSrc.Ptr(iY+1);
		}
		else if (iY==iDstHeight-1)
		{
			pfSrcAboveRow = imgSrc.Ptr(iY-1);
			pfSrcBelowRow = pfSrcRow;
		}
		else
		{
			pfSrcAboveRow = imgSrc.Ptr(iY-1);
			pfSrcBelowRow = imgSrc.Ptr(iY+1);
		}
		float *pfDstRow = imgDst.Ptr(iY);
		float *pfDstRowEnd = pfDstRow+iDstWidth;
		while(pfDstRow<pfDstRowEnd)
		{
			*pfDstRow++ = 0.25f*(*pfSrcAboveRow++ + (2.0f * (*pfSrcRow++)) + *pfSrcBelowRow++);
		}
	}

Exit:
	return hr;
}

HRESULT FilterVert121SSE2(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	HRESULT hr = S_OK;

#if (defined(_M_IX86) || defined(_M_AMD64))
	int iSrcWidth = imgSrc.Width();
	int iSrcHeight = imgSrc.Height();
	__m128 m128Two = _mm_set_ps1(2.0f);
	__m128 m128Quarter = _mm_set_ps1(0.25f);
	VT_HR_EXIT((iDstHeight < 2) ? E_FAIL : S_OK);
	VT_HR_EXIT((imgDst.Height() < iDstHeight || imgDst.Width() < iDstWidth) ? E_FAIL : S_OK);
	VT_HR_EXIT((iSrcHeight < iDstHeight || iSrcWidth < iDstWidth) ? E_FAIL : S_OK);

	for(int iY=0; iY<iDstHeight; iY++)
	{
		float *pfSrcRow = imgSrc.Ptr(iY);
		float *pfSrcAboveRow = NULL;
		float *pfSrcBelowRow = NULL;
		if (iY==0)
		{
			pfSrcAboveRow = pfSrcRow;
			pfSrcBelowRow = imgSrc.Ptr(iY+1);
		}
		else if (iY==iDstHeight-1)
		{
			pfSrcAboveRow = imgSrc.Ptr(iY-1);
			pfSrcBelowRow = pfSrcRow;
		}
		else
		{
			pfSrcAboveRow = imgSrc.Ptr(iY-1);
			pfSrcBelowRow = imgSrc.Ptr(iY+1);
		}
		float *pfDstRow = imgDst.Ptr(iY);
		int iDstWidthTrunc = 4*(iDstWidth/4);
		float *pfDstRowEnd = pfDstRow+iDstWidthTrunc;
		float *pfDstRowEnd2 = pfDstRow+iDstWidth;
		while(pfDstRow<pfDstRowEnd)
		{
			// *pfDstRow++ = 0.25f*(*pfSrcAboveRow++ + (2.0f * (*pfSrcRow++)) + *pfSrcBelowRow++);
			__m128 m128Above = _mm_load_ps(pfSrcAboveRow);
			pfSrcAboveRow +=4;
			__m128 m128ThisRow = _mm_loadu_ps(pfSrcRow);
			pfSrcRow +=4;
			__m128 m128Below = _mm_loadu_ps(pfSrcBelowRow);
			pfSrcBelowRow +=4;
			m128ThisRow = _mm_mul_ps(m128ThisRow, m128Two);
			m128ThisRow = _mm_add_ps(m128ThisRow, m128Above);
			m128ThisRow = _mm_add_ps(m128ThisRow, m128Below);
			m128ThisRow = _mm_mul_ps(m128ThisRow, m128Quarter);
			_mm_store_ps(pfDstRow, m128ThisRow);
			pfDstRow += 4;
		}
		while(pfDstRow<pfDstRowEnd2)
		{
			*pfDstRow++ = 0.25f*(*pfSrcAboveRow++ + (2.0f * (*pfSrcRow++)) + *pfSrcBelowRow++);
		}
	}

Exit:
#else
    imgSrc,imgDst,iDstWidth,iDstHeight;
#endif
	return hr;
}


HRESULT FilterVert121SubsNoSSE(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	HRESULT hr = S_OK;

	int iSrcWidth = imgSrc.Width();
	int iSrcHeight = imgSrc.Height();
	VT_HR_EXIT((iDstHeight < 2) ? E_FAIL : S_OK);
	VT_HR_EXIT((imgDst.Height() < iDstHeight || imgDst.Width() < iDstWidth) ? E_FAIL : S_OK);
	VT_HR_EXIT((iSrcHeight < 2*iDstHeight || iSrcWidth < iDstWidth) ? E_FAIL : S_OK);

	for(int iY=0; iY<iDstHeight; iY++)
	{
		float *pfSrcRow = imgSrc.Ptr(2*iY);
		float *pfSrcAboveRow = NULL;
		float *pfSrcBelowRow = NULL;
		if (iY==0)
		{
			pfSrcAboveRow = pfSrcRow;
			pfSrcBelowRow = imgSrc.Ptr(iY+1);
		}
		else
		{
			pfSrcAboveRow = imgSrc.Ptr(2*iY-1);
			pfSrcBelowRow = imgSrc.Ptr(2*iY+1);
		}
		float *pfDstRow = imgDst.Ptr(iY);
		float *pfDstRowEnd = pfDstRow+iDstWidth;
		while(pfDstRow<pfDstRowEnd)
		{
			*pfDstRow++ = 0.25f*(*pfSrcAboveRow++ + (2.0f * (*pfSrcRow++)) + *pfSrcBelowRow++);
		}
	}

Exit:
	return hr;
}

HRESULT FilterVert121SubsSSE2(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	HRESULT hr = S_OK;

#if (defined(_M_IX86) || defined(_M_AMD64))
	int iSrcWidth = imgSrc.Width();
	int iSrcHeight = imgSrc.Height();
	__m128 m128Two = _mm_set_ps1(2.0f);
	__m128 m128Quarter = _mm_set_ps1(0.25f);
	VT_HR_EXIT((iDstHeight < 2) ? E_FAIL : S_OK);
	VT_HR_EXIT((imgDst.Height() < iDstHeight || imgDst.Width() < iDstWidth) ? E_FAIL : S_OK);
	VT_HR_EXIT((iSrcHeight < 2*iDstHeight || iSrcWidth < iDstWidth) ? E_FAIL : S_OK);

	for(int iY=0; iY<iDstHeight; iY++)
	{
		float *pfSrcRow = imgSrc.Ptr(2*iY);
		float *pfSrcAboveRow = NULL;
		float *pfSrcBelowRow = NULL;
		if (iY==0)
		{
			pfSrcAboveRow = pfSrcRow;
			pfSrcBelowRow = imgSrc.Ptr(iY+1);
		}
		else
		{
			pfSrcAboveRow = imgSrc.Ptr(2*iY-1);
			pfSrcBelowRow = imgSrc.Ptr(2*iY+1);
		}
		float *pfDstRow = imgDst.Ptr(iY);
		int iDstWidthTrunc = 4*(iDstWidth/4);
		float *pfDstRowEnd = pfDstRow+iDstWidthTrunc;
		float *pfDstRowEnd2 = pfDstRow+iDstWidth;
		while(pfDstRow<pfDstRowEnd)
		{
			// *pfDstRow++ = 0.25f*(*pfSrcAboveRow++ + (2.0f * (*pfSrcRow++)) + *pfSrcBelowRow++);
			__m128 m128Above = _mm_load_ps(pfSrcAboveRow);
			pfSrcAboveRow +=4;
			__m128 m128ThisRow = _mm_loadu_ps(pfSrcRow);
			pfSrcRow +=4;
			__m128 m128Below = _mm_loadu_ps(pfSrcBelowRow);
			pfSrcBelowRow +=4;
			m128ThisRow = _mm_mul_ps(m128ThisRow, m128Two);
			m128ThisRow = _mm_add_ps(m128ThisRow, m128Above);
			m128ThisRow = _mm_add_ps(m128ThisRow, m128Below);
			m128ThisRow = _mm_mul_ps(m128ThisRow, m128Quarter);
			_mm_store_ps(pfDstRow, m128ThisRow);
			pfDstRow += 4;
		}
		while(pfDstRow<pfDstRowEnd2)
		{
			*pfDstRow++ = 0.25f*(*pfSrcAboveRow++ + (2.0f * (*pfSrcRow++)) + *pfSrcBelowRow++);
		}
	}

Exit:
#else
    imgSrc,imgDst,iDstWidth,iDstHeight;
#endif
	return hr;
}
