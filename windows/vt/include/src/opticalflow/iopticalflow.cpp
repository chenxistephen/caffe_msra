//+---------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation 2010.  All rights reserved.
//
//  Description:
//     Interface for Optical Flow Algorithms
//
//  History:
//      2011/8/11-sbaker
//          Created
//
//----------------------------------------------------------------------------

#include "stdafx.h"
#include "iopticalflow.h"
#include "cvof.h"
#ifndef VT_GCC
#include "cvof2.h"
#endif

using namespace vt;

HRESULT vt::CreateHornSchunck1Channel(IOpticalFlow* &pPyramidFlow, 
                                      int iImageWidth, int iImageHeight,
                                      int iSpeed, float fLambda)
{
    pPyramidFlow = NULL;

    VT_HR_BEGIN() 

    VT_HR_EXIT( (iSpeed < 0 || iSpeed > 5) ? E_INVALIDARG : S_OK );
    VT_HR_EXIT( (fLambda <= 0.0f) ? E_INVALIDARG : S_OK );

    pPyramidFlow = VT_NOTHROWNEW CVOF();
    VT_PTR_OOM_EXIT(pPyramidFlow);

    HORN_SCHUNCK_1CHANNEL_PARAMS.m_iSpeed = iSpeed;
    HORN_SCHUNCK_1CHANNEL_PARAMS.m_pfSmoothnessParams[0] = fLambda;
    VT_HR_EXIT( ((CVOF*)pPyramidFlow)->Allocate(iImageWidth, iImageHeight, 
                                                &HORN_SCHUNCK_1CHANNEL_PARAMS) ); 

    VT_HR_EXIT_LABEL()

    if ( hr != S_OK )
    {
        delete pPyramidFlow;
        pPyramidFlow = NULL;
    }
    return hr;
}

HRESULT vt::CreateOpticalFlow2(IOpticalFlow* &pPyramidFlow, int iImageWidth, int iImageHeight)
{
    pPyramidFlow = NULL;

#ifndef VT_GCC
    VT_HR_BEGIN() 

    pPyramidFlow = VT_NOTHROWNEW CVOF2();
    VT_PTR_OOM_EXIT(pPyramidFlow);

    VT_HR_EXIT( ((CVOF2*)pPyramidFlow)->Allocate(iImageWidth, iImageHeight, 
                                                &OPTICAL_FLOW_V2_PARAMS) ); 

    VT_HR_EXIT_LABEL()

    if ( hr != S_OK )
    {
        delete pPyramidFlow;
        pPyramidFlow = NULL;
    }

	return hr;
#else
	return E_NOTFOUND;
#endif
}

#ifndef VT_GCC

#if (defined(_M_IX86) || defined(_M_AMD64))
inline void BilinearAddressComputationSSSE3(__m128 &m128fSrc, __m128i &m128iSrc, 
     __m128i &m128iMinusOne,  __m128i &m128iZero, __m128i &m128iMax, __m128i &m128iValidity)
{
    m128iSrc = _mm_cvttps_epi32(m128fSrc);
    __m128 m128fSrc2 = _mm_cvtepi32_ps(m128iSrc);
    __m128 m128Mask = _mm_cmplt_ps(m128fSrc, m128fSrc2);
    m128iSrc = _mm_add_epi32(m128iSrc, _mm_castps_si128(m128Mask));
    __m128i m128MaskLower = _mm_cmplt_epi32(m128iMinusOne, m128iSrc);
    __m128i m128MaskUpper = _mm_cmplt_epi32(m128iZero, _mm_sub_epi32(m128iMax, m128iSrc));
    m128iValidity = _mm_and_si128(m128iValidity, m128MaskLower);
    m128iValidity = _mm_and_si128(m128iValidity, m128MaskUpper);
    m128iSrc = _mm_sign_epi32(m128iSrc, m128iValidity);
}
#endif

void vt::OpticalFlow_WarpImageSSE(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, CFloatImg &imgDst, 
                         CFloatImg &imgDataTermWeight, CRect rctDst, CFloatImg& imgSrcExtend)
{
    int iSrcWidth = imgSrc.Width();
    int iSrcHeight = imgSrc.Height();
    int iSrcWidthMax = iSrcWidth-1;
    int iSrcHeightMax = iSrcHeight-1;
    int iDstWidth = rctDst.Width();
    int iDstWidthTrunc = 4*(iDstWidth/4);
    int iDstHeight = rctDst.Height();
    unsigned char *pucSrcData = (unsigned char*) imgSrc.Ptr();
    unsigned int uiSrcStride = imgSrc.StrideBytes();
#if (defined(_M_IX86) || defined(_M_AMD64))
    __m128 m128Four = _mm_set1_ps(4.0f);
    __m128 m128One = _mm_set1_ps(1.0f);
    __m128i m128iMinusOne  = _mm_set1_epi32(-1);
    __m128 m128fZero = _mm_set1_ps(0.0f);
    __m128i m128iZero = _mm_cvttps_epi32(m128fZero);
    __m128 m128fMaxWidth = _mm_set1_ps(float(iSrcWidthMax));
    __m128i m128iMaxWidth = _mm_cvttps_epi32(m128fMaxWidth);
    __m128 m128fMaxHeight = _mm_set1_ps(float(iSrcHeightMax));
    __m128i m128iMaxHeight = _mm_cvttps_epi32(m128fMaxHeight);
#endif

    for(int iY=0; iY<iDstHeight; iY++)
    {
#if (defined(_M_IX86) || defined(_M_AMD64))
        float *fpSrcExtendRow = imgSrcExtend.Ptr(rctDst.left, rctDst.top+iY);
#endif
        float *fpDstRow = imgDst.Ptr(rctDst.left, rctDst.top+iY);
        float *fpFlowRowX = imgFlowX.Ptr(rctDst.left, rctDst.top+iY);
        float *fpFlowRowY = imgFlowY.Ptr(rctDst.left, rctDst.top+iY);
        float *fpDataTermWeight = imgDataTermWeight.Ptr(rctDst.left, rctDst.top+iY);
#if (defined(_M_IX86) || defined(_M_AMD64))
        float *fpDataTermWeightEnd = fpDataTermWeight+iDstWidthTrunc;
#endif

        float fLeft = float(rctDst.left);
#if (defined(_M_IX86) || defined(_M_AMD64))
        __m128 m128iX = _mm_set_ps(fLeft+3.0f, fLeft+2.0f, fLeft+1.0f, fLeft);
        float fY = float(rctDst.top+iY);
        __m128 m128iY = _mm_set1_ps(fY);
        while(fpDataTermWeight < fpDataTermWeightEnd)
        {
            // float fSrcX = float(iX) + *fpFlowRowX++;
            // float fSrcY = float(iY) + *fpFlowRowY++;
            __m128 m128fSrcX = _mm_load_ps(fpFlowRowX);
            fpFlowRowX+=4;
            m128fSrcX = _mm_add_ps(m128fSrcX, m128iX);
            m128iX = _mm_add_ps(m128iX, m128Four);
            __m128 m128fSrcY = _mm_load_ps(fpFlowRowY);
            fpFlowRowY+=4;
            //Adds the y-direction flow strength to each pixel y-location to get the new pixel location
            m128fSrcY = _mm_add_ps(m128fSrcY, m128iY);

            //Start computing the weight based on the newly calculated warped pixel locations.
            //Set the initial weight to 1 to assume that pixels are within boundaries.
            // int iDataTermWeight = 1;
            __m128i m128iDataTermWeight = _mm_cvttps_epi32(m128One);

            //The logic to compute weight is equivalent to the following Non-SSE version:
            // int iSrcX = int(fSrcX);
            // int iSrcY = int(fSrcY);
            // if (iSrcX < 0 || iSrcX > iSrcWidthMax || iSrcY < 0 || iSrcY > iSrcHeightMax) 
            // {
            //	iDataTermWeight = 0;
            // }
            __m128i m128iSrcX;
            __m128i m128iSrcY;
            BilinearAddressComputationSSSE3(m128fSrcY, m128iSrcY, m128iMinusOne, m128iZero,  m128iMaxHeight, m128iDataTermWeight);
            BilinearAddressComputationSSSE3(m128fSrcX, m128iSrcX, m128iMinusOne, m128iZero, m128iMaxWidth, m128iDataTermWeight);

            // *fpDataTermWeight++ = float(iDataTermWeight);
            __m128 m128fDataTermWeight = _mm_cvtepi32_ps(m128iDataTermWeight);
            _mm_store_ps(fpDataTermWeight, m128fDataTermWeight);
            fpDataTermWeight+=4;

            int iSrcX0 = m128iSrcX.m128i_i32[0];
            int iSrcY0 = m128iSrcY.m128i_i32[0];
            int iSrcX1 = m128iSrcX.m128i_i32[1];
            int iSrcY1 = m128iSrcY.m128i_i32[1];
            int iSrcX2 = m128iSrcX.m128i_i32[2];
            int iSrcY2 = m128iSrcY.m128i_i32[2];
            int iSrcX3 = m128iSrcX.m128i_i32[3];
            int iSrcY3 = m128iSrcY.m128i_i32[3];

            // float fT = fSrcX - iSrcX;
            // float fU = fSrcY - iSrcY;
            // float fTT = 1.0f - fT;
            // float fUU = 1.0f - fU;
            __m128 m128fSrcX2 = _mm_cvtepi32_ps(m128iSrcX);
            __m128 m128fSrcY2 = _mm_cvtepi32_ps(m128iSrcY);
            __m128 m128fT = _mm_sub_ps(m128fSrcX, m128fSrcX2);
            __m128 m128fU = _mm_sub_ps(m128fSrcY, m128fSrcY2);
            __m128 m128fTT = _mm_sub_ps(m128One, m128fT);
            __m128 m128fUU = _mm_sub_ps(m128One, m128fU);

            unsigned char *pucSrcData01 = pucSrcData+iSrcY0*uiSrcStride+sizeof(float)*iSrcX0;
            unsigned char *pucSrcData02 = pucSrcData01+uiSrcStride;
            float *pfSrcData01 = ((float*) pucSrcData01);
            float fPix01 = *pfSrcData01;
            pfSrcData01++;
            float fPix02 = *pfSrcData01;
            float *pfSrcData02 = ((float*) pucSrcData02);
            float fPix03 = *pfSrcData02;
            pfSrcData02++;
            float fPix04 = *pfSrcData02;

            unsigned char *pucSrcData11 = pucSrcData+iSrcY1*uiSrcStride+sizeof(float)*iSrcX1;
            unsigned char *pucSrcData12 = pucSrcData11+uiSrcStride;
            float *pfSrcData11 = ((float*) pucSrcData11);
            float fPix11 = *pfSrcData11;
            pfSrcData11++;
            float fPix12 = *pfSrcData11;
            float *pfSrcData12 = ((float*) pucSrcData12);
            float fPix13 = *pfSrcData12;
            pfSrcData12++;
            float fPix14 = *pfSrcData12;

            unsigned char *pucSrcData21 = pucSrcData+iSrcY2*uiSrcStride+sizeof(float)*iSrcX2;
            unsigned char *pucSrcData22 = pucSrcData21+uiSrcStride;
            float *pfSrcData21 = ((float*) pucSrcData21);
            float fPix21 = *pfSrcData21;
            pfSrcData21++;
            float fPix22 = *pfSrcData21;
            float *pfSrcData22 = ((float*) pucSrcData22);
            float fPix23 = *pfSrcData22;
            pfSrcData22++;
            float fPix24 = *pfSrcData22;

            unsigned char *pucSrcData31 = pucSrcData+iSrcY3*uiSrcStride+sizeof(float)*iSrcX3;
            unsigned char *pucSrcData32 = pucSrcData31+uiSrcStride;
            float *pfSrcData31 = ((float*) pucSrcData31);
            float fPix31 = *pfSrcData31;
            pfSrcData31++;
            float fPix32 = *pfSrcData31;
            float *pfSrcData32 = ((float*) pucSrcData32);
            float fPix33 = *pfSrcData32;
            pfSrcData32++;
            float fPix34 = *pfSrcData32;

            // *fpDstRow++ = fUU*(fTT*fPix1+fT*fPix2)+fU*(fTT*fPix3+fT*fPix4);
            __m128 m128Pix1 = _mm_set_ps(fPix31, fPix21, fPix11, fPix01);
            __m128 m128Pix2 = _mm_set_ps(fPix32, fPix22, fPix12, fPix02);
            __m128 m128Pix3 = _mm_set_ps(fPix33, fPix23, fPix13, fPix03);
            __m128 m128Pix4 = _mm_set_ps(fPix34, fPix24, fPix14, fPix04);
            __m128 m128Tmp1 = _mm_mul_ps(m128fTT, m128Pix1);
            __m128 m128Tmp2 = _mm_mul_ps(m128fT, m128Pix2);
            __m128 m128Tmp3 = _mm_mul_ps(m128fTT, m128Pix3);
            __m128 m128Tmp4 = _mm_mul_ps(m128fT, m128Pix4);
            m128Tmp2 = _mm_add_ps(m128Tmp1, m128Tmp2);
            m128Tmp4 = _mm_add_ps(m128Tmp3, m128Tmp4);
            m128Tmp2 = _mm_mul_ps(m128Tmp2, m128fUU);
            m128Tmp4 = _mm_mul_ps(m128Tmp4, m128fU);
            m128Tmp2 = _mm_add_ps(m128Tmp4, m128Tmp2);
            m128Tmp2 = _mm_mul_ps(m128Tmp2, m128fDataTermWeight);
            
            //Load the extend image data
            __m128 m128fSrcExtend = _mm_load_ps(fpSrcExtendRow);

            //If the weight is zero then it means that warped pixel for the current location was 
            //out of range so we should use the extend image data for this location as the warped value.
            __m128 m128Tmp5 = ::_mm_mul_ps(::_mm_xor_ps(m128fDataTermWeight, m128One), m128fSrcExtend);
            m128Tmp2 = ::_mm_add_ps(m128Tmp2, m128Tmp5);

            _mm_store_ps(fpDstRow, m128Tmp2);
            fpDstRow+=4;
            fpSrcExtendRow+=4;
        }
#endif
        for(int iX=iDstWidthTrunc; iX<iDstWidth; iX++)
        {
            float fSrcX = fLeft + float(iX) + *fpFlowRowX++;
            float fSrcY = float(rctDst.top+iY) + *fpFlowRowY++;
            int iSrcX = int(floor(fSrcX));
            int iSrcY = int(floor(fSrcY));

            int iDataTermWeight = 0;
            if (iSrcX >= 0 && iSrcX < iSrcWidthMax && iSrcY >= 0 && iSrcY < iSrcHeightMax) 
            {
                iDataTermWeight = 1;
            }
            else
            {
                float *pfSource = imgSrcExtend.Ptr(rctDst.left + iX, rctDst.top + iY);
                float *pfDestination = imgDst.Ptr(rctDst.left + iX, rctDst.top + iY);
                *pfDestination = *pfSource;
                *fpDataTermWeight++ = 0.0f;
                continue;
            }
            float fDataTermWeight = float(iDataTermWeight);
            *fpDataTermWeight++ = fDataTermWeight;
            iSrcX *= iDataTermWeight;
            iSrcY *= iDataTermWeight;

            float fT = fSrcX - iSrcX;
            float fU = fSrcY - iSrcY;
            float fTT = 1.0f - fT;
            float fUU = 1.0f - fU;
            unsigned char *pucSrcData1 = pucSrcData+iSrcY*uiSrcStride+sizeof(float)*iSrcX;
            unsigned char *pucSrcData2 = pucSrcData1+uiSrcStride;
            float *pfSrcData1 = ((float*) pucSrcData1);
            float fPix1 = *pfSrcData1;
            pfSrcData1++;
            float fPix2 = *pfSrcData1;
            float *pfSrcData2 = ((float*) pucSrcData2);
            float fPix3 = *pfSrcData2;
            pfSrcData2++;
            float fPix4 = *pfSrcData2;
            *fpDstRow++ = fDataTermWeight*(fUU*(fTT*fPix1+fT*fPix2)+fU*(fTT*fPix3+fT*fPix4));
        }
    }
}

#endif
