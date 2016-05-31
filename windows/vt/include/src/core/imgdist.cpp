//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for computing distance between images
//
//  History:
//      2004/11/08-mattu
//			Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_imgdist.h"
#include "vt_pixeliterators.h"
#include "vt_convert.h"

using namespace vt;

// forward declarations
#if (defined(_M_IX86) || defined(_M_AMD64))
int VtSADSpanSSE(const BYTE* pA,  const BYTE* pB, UINT64& dist, int iSpan);
int VtSADSpanSSE(const UInt16* pA,  const UInt16* pB, UINT64& dist, int iSpan);
int VtSADSpanSSE(const float* pA, const float* pB, double& dist, int iSpan);
#endif

// inline utilties
int   sadabs(Byte v1, Byte v2)     { return abs(int(v1)-int(v2)); }
int   sadabs(UInt16 v1, UInt16 v2) { return abs(int(v1)-int(v2)); }
float sadabs(float v1, float v2)   { return abs(v1-v2); }

//+-----------------------------------------------------------------------
//
// Utilities -
// 
// TODO: seems like this type of error check code is used in lots of places
//       i.e. imgmath.cpp would probably be good to have a central source file
//       for these functions i.e. errorcheck.cpp
//------------------------------------------------------------------------
static HRESULT CheckArgs(const CImg &imgsrc1, const CImg &imgsrc2)
{
	if(!imgsrc1.IsValid() || !imgsrc2.IsValid())
		return E_INVALIDARG;
	if(imgsrc1.Width()!=imgsrc2.Width() || imgsrc1.Height()!=imgsrc2.Height())
		return E_INVALIDARG;
	if(imgsrc1.Bands()!=imgsrc2.Bands())
		return E_INVALIDARG;
	if(imgsrc1.ElSize()!=imgsrc2.ElSize())
		return E_INVALIDARG;
	return NOERROR;
}

//+-----------------------------------------------------------------------
//
// Compute the sum of absolute differences (SAD) between two images
//
//------------------------------------------------------------------------
template <class T, class TResult>
HRESULT VtSADImagesInternal( const CTypedImg<T>& imgA, const CTypedImg<T>& imgB,
	TResult& dist)
{
	HRESULT hr;

	VT_HR_EXIT( CheckArgs(imgA, imgB) );

	dist = 0;

	int iSpan = imgA.Width() * imgA.Bands();

	for( int y = 0; y < imgA.Height(); y++ )
	{
		dist += VtSADSpan<T, TResult>(imgA.Ptr(y), imgB.Ptr(y), iSpan);
	}

Exit:
	return hr;
}

template <class T, class TResult>
TResult vt::VtSADSpan(const T* pA, const T* pB, int iSpan)
{
	TResult dist = 0;

	int iX = 0;
#if (defined(_M_IX86) || defined(_M_AMD64))
	if( g_SupportSSE2() )
	{
		iX = VtSADSpanSSE(pA, pB, dist, iSpan);
		pA += iX;
		pB += iX;
	}
#endif

	for( ; iX < iSpan; iX++, pA++, pB++ )
	{
		dist += sadabs(*pA, *pB);
	}

	return dist;
}

template <class T, class TResult>
TResult vt::VtSADSpan(const RGBAType<T>* pA, const RGBAType<T>* pB, int iSpan)
{
	TResult dist = 0;

	int iX = 0;

	for( ; iX < iSpan; iX++, pA++, pB++ )
	{
		dist += sadabs(pA->r, pB->r);
		dist += sadabs(pA->g, pB->g);
		dist += sadabs(pA->b, pB->b);
	}

	return dist;
}

// force instantiation
template UINT64 vt::VtSADSpan(const RGBAPix* pA, const RGBAPix* pB, int iSpan);  
template UINT64 vt::VtSADSpan(const RGBAShortPix* pA, const RGBAShortPix* pB, int iSpan);  
template double vt::VtSADSpan(const RGBAFloatPix* pA, const RGBAFloatPix* pB, int iSpan);  

#if (defined(_M_IX86) || defined(_M_AMD64))
//+-----------------------------------------------------------------------------
//
// Byte specialization
//
//------------------------------------------------------------------------------
int VtSADSpanSSE(const BYTE* pA, const BYTE* pB, UINT64& dist, int iSpan)
{
	int x = 0;

	// test if it's possible to use aligned loads
	if( (INT_PTR(pA) & 0xf) == (INT_PTR(pB) & 0xf) )
	{
		// align inputs to 16byte boundary
		for( ; (INT_PTR(pA) & 0xf) != 0 && x < iSpan; x++, pA++, pB++ )
		{
			dist += abs(int(*pA)-int(*pB));
		}

		// use aligned SSE loads for remainder
		__m128i x0 = _mm_setzero_si128();
		for( ; x <= iSpan-16; x+=16, pA+=16, pB+=16 )
		{
			__m128i x1 = _mm_load_si128 ((__m128i*)pA);
			__m128i x2 = _mm_load_si128 ((__m128i*)pB);
			x1 = _mm_sad_epu8(x1, x2);
			x0 = _mm_add_epi64 (x0, x1);
		}

		dist += x0.m128i_u64[0] + x0.m128i_u64[1];
	}
	else
	{
		// use non-aligned SSE loads for everything
		__m128i x0 = _mm_setzero_si128();
		for( ; x <= iSpan-16; x+=16, pA+=16, pB+=16 )
		{
			__m128i x1 = _mm_loadu_si128 ((__m128i*)pA);
			__m128i x2 = _mm_loadu_si128 ((__m128i*)pB);
			x1 = _mm_sad_epu8(x1, x2);
			x0 = _mm_add_epi64 (x0, x1);
		}

		dist += x0.m128i_u64[0] + x0.m128i_u64[1];
	}

	return x;
}

//+-----------------------------------------------------------------------------
//
// UInt16 specialization
//
//------------------------------------------------------------------------------
int VtSADSpanSSE(const UInt16* pA, const UInt16* pB, UINT64& dist, int iSpan)
{
	int x = 0;

	__m128i x0 = _mm_setzero_si128();
	__m128i x5 = _mm_setzero_si128();
	__m128i x6 = _mm_set1_epi16(short(-0x8000));

    // test if it's possible to use aligned loads
	if( (INT_PTR(pA) & 0xf) == (INT_PTR(pB) & 0xf) )
	{
		// align inputs to 16byte boundary
		for( ; (INT_PTR(pA) & 0xf) != 0 && x < iSpan; x++, pA++, pB++ )
		{
			dist += abs(int(*pA) - int(*pB));
		}

		// use aligned SSE loads for remainder
		for( ; x <= iSpan - 8; x += 8, pA += 8, pB += 8 )
		{
			__m128i x1 = _mm_load_si128 ((__m128i*) pA);
			__m128i x2 = _mm_load_si128 ((__m128i*) pB);

            x1 = _mm_sub_epi16(x1, x6);         // convert to signed
            x2 = _mm_sub_epi16(x2, x6);

            __m128i x3 = _mm_min_epi16(x1, x2); // abs
			__m128i x4 = _mm_max_epi16(x1, x2);
			x1 = _mm_sub_epi16(x4, x3);

            x2 = _mm_unpackhi_epi16(x1, x0);    // convert to UINT32
            x3 = _mm_unpacklo_epi16(x1, x0);
			x4 = _mm_add_epi32(x2, x3);         // sum

            x2 = _mm_unpackhi_epi32(x4, x0);    // convert to UINT64
            x3 = _mm_unpacklo_epi32(x4, x0);
			x4 = _mm_add_epi64(x2, x3);         // sum

			x5 = _mm_add_epi64(x5, x4);         // total
		}

		dist += x5.m128i_u64[0] + x5.m128i_u64[1];
	}
	else
	{
		// use non-aligned SSE loads for everything
		for( ; x <= iSpan - 8; x += 8, pA += 8, pB += 8 )
		{
			__m128i x1 = _mm_loadu_si128 ((__m128i*) pA);
			__m128i x2 = _mm_loadu_si128 ((__m128i*) pB);

            x1 = _mm_sub_epi16(x1, x6);         // convert to signed
            x2 = _mm_sub_epi16(x2, x6);

            __m128i x3 = _mm_min_epi16(x1, x2); // abs
            __m128i x4 = _mm_max_epi16(x1, x2);
			x1 = _mm_sub_epi16(x4, x3);

            x2 = _mm_unpackhi_epi16(x1, x0);    // convert to UINT32
            x3 = _mm_unpacklo_epi16(x1, x0);
			x4 = _mm_add_epi32(x2, x3);         // sum

            x2 = _mm_unpackhi_epi32(x4, x0);    // convert to UINT64
            x3 = _mm_unpacklo_epi32(x4, x0);
			x4 = _mm_add_epi64(x2, x3);         // sum

			x5 = _mm_add_epi64(x5, x4);         // total
		}

		dist += x5.m128i_u64[0] + x5.m128i_u64[1];
	}

	return x;
}

//+-----------------------------------------------------------------------------
//
// float specialization
//
//------------------------------------------------------------------------------
int VtSADSpanSSE(const float* pA, const float* pB, double& dist, int iSpan)
{ 
	int x = 0;

	// test if it's possible to use aligned loads
	if( (INT_PTR(pA) & 0xf) == (INT_PTR(pB) & 0xf) )
	{
		// align inputs to 16byte boundary
		for( ; (INT_PTR(pA) & 0xf) != 0 && x < iSpan; x++, pA++, pB++ )
		{
			dist += abs(*pA-*pB);
		}

		// use aligned SSE loads for remainder
		__m128d x0 = _mm_setzero_pd();
		for( ; x <= iSpan-4; x+=4, pA+=4, pB+=4 )
		{
			__m128 x1 = _mm_load_ps (pA);    // load
			__m128 x2 = _mm_load_ps (pB);

			__m128 x3 = _mm_min_ps(x1,x2);   // abs
			x1 = _mm_max_ps(x1,x2);
			x1 = _mm_sub_ps(x1,x3);

			__m128d x4 = _mm_cvtps_pd (x1); // convert to double
			x1 = _mm_shuffle_ps(x1, x1, _MM_SHUFFLE(1,0,3,2));
			__m128d x5 = _mm_cvtps_pd (x1);

			x4 = _mm_add_pd(x4, x5);         // sum
			x0 = _mm_add_pd(x0, x4);
		}

		dist += x0.m128d_f64[0] + x0.m128d_f64[1];
	}
	else
	{
		// use non-aligned SSE loads for everything
		__m128d x0 = _mm_setzero_pd();
		for( ; x <= iSpan-4; x+=4, pA+=4, pB+=4 )
		{
			__m128 x1 = _mm_loadu_ps (pA);   // load
			__m128 x2 = _mm_loadu_ps (pB);

			__m128 x3 = _mm_min_ps(x1,x2);   // abs
			x1 = _mm_max_ps(x1,x2);
			x1 = _mm_sub_ps(x1,x3);

			__m128d x4 = _mm_cvtps_pd (x1); // convert to double
			x1 = _mm_shuffle_ps(x1, x1, _MM_SHUFFLE(1,0,3,2));
			__m128d x5 = _mm_cvtps_pd (x1);

			x4 = _mm_add_pd(x4, x5);         // sum
			x0 = _mm_add_pd(x0, x4);
		}

		dist += x0.m128d_f64[0] + x0.m128d_f64[1];
	}

	return x;
}
#endif

//+-----------------------------------------------------------------------
//
// Compute the sum of absolute differences (SAD) between two images
//
//------------------------------------------------------------------------
HRESULT vt::VtSADImages(const CImg& imgA, const CImg& imgB, double& dist)
{
	HRESULT hr;

	VT_HR_EXIT(CheckArgs(imgA, imgB));

	switch (EL_FORMAT(imgA.GetType()))
	{
	case EL_FORMAT_BYTE:
		{
			UINT64 int_dist;
			
			VT_HR_EXIT(VtSADImagesInternal((CTypedImg<Byte>&)imgA, 
				(CTypedImg<Byte>&)imgB, int_dist));

			dist = (double) int_dist;			
			break;
		}
	case EL_FORMAT_SHORT:
		{
			UINT64 int_dist;
			
			VT_HR_EXIT(VtSADImagesInternal((CTypedImg<UInt16>&)imgA, 
				(CTypedImg<UInt16>&)imgB, int_dist));

			dist = (double) int_dist;			
			break;
		}
	case EL_FORMAT_FLOAT:
		{
			VT_HR_EXIT(VtSADImagesInternal((CTypedImg<float>&)imgA, 
				(CTypedImg<float>&)imgB, dist));
			break;
		}
	default:
		{
			VT_ASSERT(false);
		}
	}

Exit:
	return hr;
}


//+-----------------------------------------------------------------------------
//
// Function: CompareImages
// 
// Synopsis: Compare the pixel values in 2 images and report differences
//
//------------------------------------------------------------------------------
template <class T>
bool CompareBlock(const CTypedImg<T>& imgBlk1, const CTypedImg<T>& imgBlk2, 
				  const vt::CRect& rctBlk, float fMaxAllowedDiff,
				  vt::vector<vt::wstring>* pvecMessages, int& iMessageCount, 
				  int iMaxMessageCount)
{
	int iBands = imgBlk1.Bands();

	bool bMatch = true;
	for( int y = 0; y < rctBlk.Height(); y++ )
	{
		const int c_bufsize = 256;
		float buff1[c_bufsize];
		float buff2[c_bufsize];
		int   useablebuf = (c_bufsize / iBands) * iBands;

		for( int x = 0; x < rctBlk.Width()*iBands; x+=useablebuf )
		{
			int curbuf = VtMin(useablebuf, rctBlk.Width()*iBands-x);
			HRESULT hr1 = VtConvertSpan(buff1, VT_IMG_MAKE_TYPE(EL_FORMAT_FLOAT, iBands),
										imgBlk1.Ptr(y)+x, imgBlk1.GetType(),
										curbuf);
			HRESULT hr2 = VtConvertSpan(buff2, VT_IMG_MAKE_TYPE(EL_FORMAT_FLOAT, iBands),
				                        imgBlk2.Ptr(y)+x, imgBlk2.GetType(),
						                curbuf);

			UNREFERENCED_PARAMETER(hr1);
			UNREFERENCED_PARAMETER(hr2);
			VT_ASSERT( hr1==S_OK && hr2==S_OK );

			const float* p1 = buff1;
			const float* p2 = buff2;
			for( int c = 0; c < curbuf; c++, p1++, p2++ )
			{
				float fDiff = fabs(*p1-*p2);
				if( fDiff > fMaxAllowedDiff )
				{
					if( pvecMessages && (iMessageCount < iMaxMessageCount) )
					{
						int xpix = (x+c)/iBands;
						int band = ((x+c)-xpix)*iBands;
	
						vt::wstring msg;
						msg.format_with_resize(L"VtCompareImages: pixel mismatch at"\
							L" location(%d,%d,%d), %f != %f",
							rctBlk.left+xpix, rctBlk.top+y, band, *p1, *p2);
						pvecMessages->push_back(msg);
						iMessageCount++;
					}
					else
					{
						return false;
					}
					bMatch = false;
				}
			}

			x += curbuf;
		}
	}
	return bMatch;
}

bool 
vt::VtCompareImages(IImageReader* pImg1, IImageReader* pImg2,
					int iMaxMessageCount, vt::vector<vt::wstring>* pvecMessages,
					float fMaxAllowedDiff, const vt::CRect* pRect)
{
	vt::wstring msg;

	CImgInfo info1 = pImg1->GetImgInfo();
	CImgInfo info2 = pImg2->GetImgInfo();

	if ( info1.width != info2.width || info1.height != info2.height )
	{
		if( pvecMessages && iMaxMessageCount )
		{
			msg.format_with_resize(L"VtCompareImages: image dimensions mismatch"\
				L" (%d,%d) vs. (%d,%d)",
				info1.width, info1.height, 
				info2.width, info2.height);
			pvecMessages->push_back(msg);
		}
		return false;
	}
	if ( info1.PixSize() != info2.PixSize() )
	{
		if( pvecMessages && iMaxMessageCount )
		{
			msg.format_with_resize(L"VtCompareImages: pixel size mismatch");
			pvecMessages->push_back(msg);
		}
		return false;
	}

	bool bMatch = true;

	CImg imgBlk1, imgBlk2;
	int iMessageCount = 0;
	for ( CBlockIterator bi(pRect? *pRect: info1.Rect(), 256); 
		!bi.Done(); bi.Advance() )
	{
		CRect rctBlk = bi.GetCompRect();
		pImg1->ReadRegion(rctBlk, imgBlk1);
		pImg2->ReadRegion(rctBlk, imgBlk2);
		switch( EL_FORMAT(info1.type) )
		{
		case EL_FORMAT_BYTE:
			bMatch = CompareBlock((CByteImg&)imgBlk1, (CByteImg&)imgBlk2, 
								  rctBlk, fMaxAllowedDiff, pvecMessages, 
								  iMessageCount, iMaxMessageCount);
			break;
		case EL_FORMAT_SHORT:
			bMatch = CompareBlock((CShortImg&)imgBlk1, (CShortImg&)imgBlk2, 
								  rctBlk, fMaxAllowedDiff, pvecMessages,
								  iMessageCount, iMaxMessageCount);
			break;
		case EL_FORMAT_HALF_FLOAT:
			bMatch = CompareBlock((CHalfFloatImg&)imgBlk1, (CHalfFloatImg&)imgBlk2,
								  rctBlk, fMaxAllowedDiff, pvecMessages,
								  iMessageCount, iMaxMessageCount);
			break;
		case EL_FORMAT_FLOAT:
			bMatch = CompareBlock((CFloatImg&)imgBlk1, (CFloatImg&)imgBlk2,
								  rctBlk, fMaxAllowedDiff, pvecMessages,
								  iMessageCount, iMaxMessageCount);
			break;
		default:
			bMatch = false;
			break;
		}
		if ( !bMatch && (pvecMessages==NULL || iMessageCount>=iMaxMessageCount) )
		{
			return false;
		}
	}
	return bMatch;
}

