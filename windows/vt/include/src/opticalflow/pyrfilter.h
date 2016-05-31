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

#pragma once

#include "vtcore.h"
using namespace vt;

/// <summary> SSE2 implementation of FilterHoriz121() </summary>
HRESULT FilterHoriz121SSE2(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight);
/// <summary> Straight C++ implementation of FilterHoriz121() </summary>
HRESULT FilterHoriz121NoSSE(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight);
/// <summary> Function to filter an image with a 121 filter in horizontal direction </summary>
inline HRESULT FilterHoriz121(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	if (g_SupportSSE2())
	{
		return FilterHoriz121SSE2(imgSrc, imgDst, iDstWidth, iDstHeight);	
	}
	else 
	{
		return FilterHoriz121NoSSE(imgSrc, imgDst, iDstWidth, iDstHeight);	
	}
}

/// <summary> SSE2 implementation of FilterVert121() </summary>
HRESULT FilterVert121SSE2(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight);
/// <summary> Straight C++ implementation of FilterVert121() </summary>
HRESULT FilterVert121NoSSE(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight);
/// <summary> Function to filter an image with a 121 filter in vertical direction </summary>
inline HRESULT FilterVert121(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	if (g_SupportSSE2())
	{
		return FilterVert121SSE2(imgSrc, imgDst, iDstWidth, iDstHeight);	
	}
	else 
	{
		return FilterVert121NoSSE(imgSrc, imgDst, iDstWidth, iDstHeight);	
	}
}

/// <summary> SSE2 implementation of FilterHoriz121Subs() </summary>
HRESULT FilterHoriz121SubsSSE2(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight);
/// <summary> Straight C++ implementation of FilterHoriz121Subs() </summary>
HRESULT FilterHoriz121SubsNoSSE(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight);
/// <summary> Function to filter an image with a 121 filter in horizontal direction and subsample
/// to the next higher level in a pyramid </summary>
inline HRESULT FilterHoriz121Subs(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	if (g_SupportSSE2())
	{
		return FilterHoriz121SubsSSE2(imgSrc, imgDst, iDstWidth, iDstHeight);	
	}
	else 
	{
		return FilterHoriz121SubsNoSSE(imgSrc, imgDst, iDstWidth, iDstHeight);	
	}
}

/// <summary> SSE2 implementation of FilterVert121Subs() </summary>
HRESULT FilterVert121SubsSSE2(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight);
/// <summary> Straight C++ implementation of FilterVert121Subs() </summary>
HRESULT FilterVert121SubsNoSSE(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight);
/// <summary> Function to filter an image with a 121 filter in vertical direction and subsample
/// to the next higher level in a pyramid </summary>
inline HRESULT FilterVert121Subs(CFloatImg &imgSrc, CFloatImg &imgDst, int iDstWidth, int iDstHeight)
{
	if (g_SupportSSE2())
	{
		return FilterVert121SubsSSE2(imgSrc, imgDst, iDstWidth, iDstHeight);	
	}
	else 
	{
		return FilterVert121SubsNoSSE(imgSrc, imgDst, iDstWidth, iDstHeight);	
	}
}


