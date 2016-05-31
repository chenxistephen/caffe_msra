//+---------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation 2010.  All rights reserved.
//
//  Description:
//     Variational Optical Flow
//          Implements a class of optical flow algorithms such as Horn-Schunk
//          Black-Anandan, etc. The caller passes a weighting function that
//          determines the exact algorithm implemented.
//
//  History:
//      2010/4/13-sbkang
//          Prototype Implementation of Black-Anandan Optical Flow
//      2010/4/28-sbaker
//			Created and Modified for Improved Efficiency
//            Black Anandan -> Horn Schunck
//            Bidirectional -> Unidirectional
//            Inverse Compositional Update
//            Added an Inner Loop to the Linear Solver
//            Misc Small Changes in Filtering, Etc
//      2011/4/27-sbaker
//            Modifed computation of linear system and SOR for more efficiency
//            Added ability to weight the smoothness and data term
//      2011/8/11-sbaker
//            Added parameter struct
//            Added option to pass weighting function
//      2011/8/24-sbaker
//            Experimented with multicore: In end, decided to leave out
//            The benefits for 640x480 and below are marginal (20% best on 4-core)
//            Decided best approach for VGA and below is to add multicore at image
//            level or above (ie. call ComputeFlow() for multiple images on multiple
//            threads.)
//
//----------------------------------------------------------------------------

#include "stdafx.h"
#include "cvof.h"

// Timing code
#ifndef VT_GCC
#include "oftimer.h"
#else
#define TIMER_START_RGB2Y
#define TIMER_STOP_RGB2Y
#define TIMER_START_PYRAMID
#define TIMER_STOP_PYRAMID
#define TIMER_START_DFLOW
#define TIMER_STOP_DFLOW
#define TIMER_START_W2
#define TIMER_STOP_W2
#define TIMER_START_ERRORIM
#define TIMER_STOP_ERRORIM
#define TIMER_START_LAP
#define TIMER_STOP_LAP
#define TIMER_START_WEIGHT
#define TIMER_STOP_WEIGHT
#define TIMER_START_LS
#define TIMER_STOP_LS
#define TIMER_START_SOR
#define TIMER_STOP_SOR
#define TIMER_START_EXPAND
#define TIMER_STOP_EXPAND
#define TIMER_START_DERIVS
#define TIMER_STOP_DERIVS
#endif

// Filtering Routines for Pyramid
#include "pyrfilter.h"

// Simple Macros
#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a > b ? b : a)


/// Horn-Schunck 
/// A 1-channel implementation of Horn-Schunck Optical Flow
/// The underlying representation of the image is a single channel. RGB images
/// are first converted to a single Y channel and optical flow esimated on that.
/// The algorithm has two parameters. 
/// 1) "Lambda" (default 500.f) is the smoothing parameter. Increasing lambda results
/// in more and more smoothing. As Horn-Schunck using L2 penalty functions, the
/// flow fields are always somewhat over-smoothed across depth/motion
/// discontinuities. Horn-Schunck flow is sufficient, however, for many video 
/// enhancement applications such as stabilization and rolling shutter correction.
/// 2) "Speed" (default 4). When Speed = 5, the number of interations of the 
/// underlying solvers is reduced to an absolute minimum, giving the fastest
/// possible execution time. When Speed is reduced to 4, 3, 2, 1, or 0, more
/// iterations are performed. The quality of the flow steady improves, but the
/// algorithm takes more and more time.
void HORN_SCHUNCK_1CHANNEL_WEIGHTING_FN(CFloatImg &imgSmoothness, float *pfSmoothnessParams)
{
    float fLambda = pfSmoothnessParams[0];
    int iWidth = imgSmoothness.Width();
    int iHeight = imgSmoothness.Height();
    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfSmoothnessRow = imgSmoothness.Ptr(iY);
        float *pfSmoothnessRowEnd = pfSmoothnessRow + iWidth;
        while(pfSmoothnessRow < pfSmoothnessRowEnd)
        {
            *pfSmoothnessRow++ = fLambda;
        }
    }
}		
static int HORN_SCHUNCK_1CHANNEL_ALGORITHM = 0;
static float HORN_SCHUNCK_1CHANNEL_DEFAULT_LAMBDA = 500.0f;
static const float HORN_SCHUNCK_1CHANNEL_DEFAULT_W = 1.95f;
static const int HORN_SCHUNCK_1CHANNEL_PREBLURCOUNT = 0;
static const int HORN_SCHUNCK_1CHANNEL_TOP_LEVEL_MIN_PIXEL_COUNT = 600;
static const int HORN_SCHUNCK_1CHANNEL_MAX_ITERATION_LEVEL = 3;
static const int HORN_SCHUNCK_1CHANNEL_DEFAULT_SPEED = 4;
static int HORN_SCHUNCK_1CHANNEL_FAST_MAX_INNER_ITERATIONS[HORN_SCHUNCK_1CHANNEL_MAX_ITERATION_LEVEL+1] = { 3, 3, 5, 7 };
static int HORN_SCHUNCK_1CHANNEL_FAST_MAX_OUTER_ITERATIONS[HORN_SCHUNCK_1CHANNEL_MAX_ITERATION_LEVEL+1] = { 1, 2, 3, 4 };
static int HORN_SCHUNCK_1CHANNEL_FULL_MAX_INNER_ITERATIONS[HORN_SCHUNCK_1CHANNEL_MAX_ITERATION_LEVEL+1] = { 8, 8, 10, 12 };
static int HORN_SCHUNCK_1CHANNEL_FULL_MAX_OUTER_ITERATIONS[HORN_SCHUNCK_1CHANNEL_MAX_ITERATION_LEVEL+1] = { 6, 7, 8, 9 };
CVOFParams vt::HORN_SCHUNCK_1CHANNEL_PARAMS(HORN_SCHUNCK_1CHANNEL_PREBLURCOUNT, 
                                            HORN_SCHUNCK_1CHANNEL_TOP_LEVEL_MIN_PIXEL_COUNT,
                                            &HORN_SCHUNCK_1CHANNEL_DEFAULT_LAMBDA, 
                                            HORN_SCHUNCK_1CHANNEL_DEFAULT_W,
                                            HORN_SCHUNCK_1CHANNEL_DEFAULT_SPEED,
                                            HORN_SCHUNCK_1CHANNEL_MAX_ITERATION_LEVEL, 
                                            HORN_SCHUNCK_1CHANNEL_FAST_MAX_INNER_ITERATIONS,
                                            HORN_SCHUNCK_1CHANNEL_FAST_MAX_OUTER_ITERATIONS, 
                                            HORN_SCHUNCK_1CHANNEL_FULL_MAX_INNER_ITERATIONS,
                                            HORN_SCHUNCK_1CHANNEL_FULL_MAX_OUTER_ITERATIONS,
                                            &HORN_SCHUNCK_1CHANNEL_WEIGHTING_FN);
const int HORN_SCHUNCK_1CHANNEL_PUBLIC_PARAM_COUNT = 2;
const int HORN_SCHUNCK_1CHANNEL_PUBLIC_PARAM_MAX_LENGTH = 20;
char HORN_SCHUNCK_1CHANNEL_PUBLIC_PARAM_LAMBDA[HORN_SCHUNCK_1CHANNEL_PUBLIC_PARAM_MAX_LENGTH];
char HORN_SCHUNCK_1CHANNEL_PUBLIC_PARAM_SPEED[HORN_SCHUNCK_1CHANNEL_PUBLIC_PARAM_MAX_LENGTH];

void CVOF::Init()
{
    m_bAllocated = false;
    m_iCurrentImage = 1;
    m_iPreviousImage = 0;
    m_iNumberImages = 0;
    m_iInputWidth = 0;
    m_iInputHeight = 0;
    m_iBaseWidth = 0;
    m_iBaseHeight = 0;
    m_iSubsample = 1;
    m_iLevels = 0;
    m_iPyramidTopLevelPixelCount = 100;
    m_iPreBlurCount = 1;
    m_pfSmoothnessParams = NULL;
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
    m_pfnWeightingFn = NULL;
#endif
    m_fW = 1.775f;
    m_fOneMinusW = 1.0f-m_fW;
    m_bSSE2 = false;
    m_bSSSE3 = false;
    m_bSSE4_1 = false;
    m_iMaxIterationLevel = 0;
    m_piInnerIterations = NULL;
    m_piOuterIterations = NULL;
    if (g_SupportSSE2())
    {
        m_bSSE2 = true;
    }
    if (g_SupportSSSE3())
    {
        m_bSSSE3 = true;
    }
    if (g_SupportSSE4_1())
    {
        m_bSSE4_1 = true;
    }
//	m_ri = NULL;
}

HRESULT CVOF::Allocate(int iInputWidth, int iInputHeight, CVOFParams *cParam)
{
    VT_HR_BEGIN()

    Deallocate();
   
    // Process Params
    m_pfSmoothnessParams = cParam->m_pfSmoothnessParams;
    m_fW = cParam->m_fW;
    m_fOneMinusW = 1.0f-m_fW;
    m_iMaxIterationLevel = cParam->m_iMaxIterationLevel;
    m_iPreBlurCount = cParam->m_iPreBlurCount;
    m_iPyramidTopLevelPixelCount = cParam->m_iPyramidTopLevelPixelCount;
    m_piInnerIterationsFast = cParam->m_piInnerIterationsFast;
    m_piOuterIterationsFast = cParam->m_piOuterIterationsFast;
    m_piInnerIterationsFull = cParam->m_piInnerIterationsFull;
    m_piOuterIterationsFull = cParam->m_piOuterIterationsFull;
    m_iSpeed = cParam->m_iSpeed;
    VT_PTR_EXIT(m_piInnerIterations = VT_NOTHROWNEW int[m_iMaxIterationLevel+1]);
    VT_PTR_EXIT(m_piOuterIterations = VT_NOTHROWNEW int[m_iMaxIterationLevel+1]);
    SetInterationCountsFromSpeed();

#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
    m_pfnWeightingFn = cParam->m_pfnWeightingFn;
#endif

    // Set the Input Width and Height
    m_iInputWidth = iInputWidth;
    m_iInputHeight = iInputHeight;
    
    // Check Base Image Size is at least 2x2
    m_iBaseWidth = int(ceil(float(m_iInputWidth) / float(m_iSubsample)));
    m_iBaseHeight = int(ceil(float(m_iInputHeight) / float(m_iSubsample)));
    VT_HR_EXIT((m_iBaseWidth < 2 || m_iBaseHeight < 2) ? E_INVALIDARG : S_OK);

    // Compute Number of Levels Reasonable
    m_iLevels = 1;
    int iBaseWidth = m_iBaseWidth;
    int iBaseHeight = m_iBaseHeight;
    while(iBaseWidth*iBaseHeight>4*m_iPyramidTopLevelPixelCount)
    {
        m_iLevels++;
        iBaseWidth /= 2;
        iBaseHeight /= 2;
    }

    // none of these pyramids should autogenerate their levels
    PYRAMID_PROPERTIES pyrprops;
    pyrprops.eAutoFilter = ePyramidFilterNone;

    // Allocate the Input Image

    VT_HR_EXIT(m_imgTmpImage.Create(m_iBaseWidth, m_iBaseHeight));
    VT_HR_EXIT(m_imgTmpImage2.Create(m_iBaseWidth, m_iBaseHeight));

    // Set the Base Levels, Force Allocation, and Add Kernel for the Input Image Pyramids
    VT_HR_EXIT(m_ppyrInputImages[0].Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_ppyrInputImages[1].Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));

    // Intermediate Image Pyramids
    VT_HR_EXIT(m_pyrDerivsX.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrDerivsY.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrWarpedInput.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrErrorImage.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));

    // Flow Pyramids
    VT_HR_EXIT(m_pyrFlowX.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrFlowY.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrDFlowX.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrDFlowY.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrFlowLaplacianX.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrFlowLaplacianY.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));

    // Regularization
    VT_HR_EXIT(m_pyrDataTermWeight.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
    VT_HR_EXIT(m_pyrSmoothnessWeight.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
#endif

    // Linear System
    VT_HR_EXIT(m_pyrLinearSystem1.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrLinearSystem2.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrLinearSystem3.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrLinearSystem4.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrLinearSystem5.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrLinearSystem6.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_imgSORX.Create(m_iBaseWidth, 1, 1));
    VT_HR_EXIT(m_imgSORY.Create(m_iBaseWidth, 1, 1));
#ifdef OUTPUT_INNER_ITERATION_ERROR
    VT_HR_EXIT(m_pyrInnerErrorFactorX.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
    VT_HR_EXIT(m_pyrInnerErrorFactorY.Create(m_iBaseWidth, m_iBaseHeight, &pyrprops));
#endif

//	VT_HR_EXIT(((m_ri = VT_NOTHROWNEW CVTResourceManagerInit()) == NULL) ? E_FAIL : S_OK);

    // Finally set allocated flag
    m_bAllocated = true;

    VT_HR_EXIT_LABEL()

    if (hr != S_OK)
    {
        Deallocate();
    }
    return hr;
}

void CVOF::SetInterationCountsFromSpeed()
{
    if (m_piInnerIterations != NULL)
    {
        // Clamp speed to valid range
        m_iSpeed = MAX(0, MIN(5, m_iSpeed));
        for(int i=0; i<=m_iMaxIterationLevel; i++)
        {
            m_piInnerIterations[i] = (m_iSpeed * m_piInnerIterationsFast[i] +
                                      (5 - m_iSpeed) * m_piInnerIterationsFull[i])/5; 
            m_piOuterIterations[i] = (m_iSpeed * m_piOuterIterationsFast[i] +
                                      (5 - m_iSpeed) * m_piOuterIterationsFull[i])/5; 
        }
    }
}

void CVOF::Deallocate()
{
    m_imgTmpImage.Deallocate();
    m_imgTmpImage2.Deallocate();
    m_ppyrInputImages[0].Deallocate();
    m_ppyrInputImages[1].Deallocate();
    m_pyrDerivsX.Deallocate();
    m_pyrDerivsY.Deallocate();
    m_pyrWarpedInput.Deallocate();
    m_pyrWarpedInput.Deallocate();
    m_pyrErrorImage.Deallocate();
    m_pyrFlowX.Deallocate();
    m_pyrFlowY.Deallocate();
    m_pyrDFlowX.Deallocate();
    m_pyrDFlowY.Deallocate();
    m_pyrFlowLaplacianX.Deallocate();
    m_pyrFlowLaplacianY.Deallocate();
    m_pyrDataTermWeight.Deallocate();
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
    m_pyrSmoothnessWeight.Deallocate();
#endif
    m_pyrLinearSystem1.Deallocate();
    m_pyrLinearSystem2.Deallocate();
    m_pyrLinearSystem3.Deallocate();
    m_pyrLinearSystem4.Deallocate();
    m_pyrLinearSystem5.Deallocate();
    m_pyrLinearSystem6.Deallocate();
    m_imgSORX.Deallocate();
    m_imgSORY.Deallocate();

    if (m_piInnerIterations != NULL)
    {
        delete [] m_piInnerIterations;
        m_piInnerIterations = NULL;
    }
    if (m_piOuterIterations != NULL)
    {
        delete [] m_piOuterIterations;
        m_piOuterIterations = NULL;
    }
#ifdef OUTPUT_INNER_ITERATION_ERROR
    m_pyrInnerErrorFactorX.Deallocate();
    m_pyrInnerErrorFactorY.Deallocate();
#endif
//	if (m_ri != NULL)
//	{
//		delete m_ri;
//		m_ri = NULL;
//	}
    Init();
}

HRESULT CVOF::AddBGRAImage(unsigned char* pucBuffer, int iBufferWidth, int iBufferHeight, 
                           size_t uiBufferStrideBytes)
{
    VT_HR_BEGIN()

    // Check Allocated and Input Buffer is the Right Size
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iBufferWidth != m_iInputWidth || iBufferHeight != m_iInputHeight) ? E_INVALIDARG : S_OK);

    // Copy into Base Level, Converting RGB 2 Luminance
    TIMER_START_RGB2Y;
    m_iCurrentImage = m_iPreviousImage;
    m_iPreviousImage = 1-m_iCurrentImage;
    m_iNumberImages++;
    unsigned char *pucDestinationData;
    unsigned int uiDestinationStride;
#if 0
    if (iPreBlurCount%2 == 0)
    {
        pucDestinationData = (unsigned char*) m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Ptr();
        uiDestinationStride = m_ppyrInputImages[m_iCurrentImage].GetLevel(0).StrideBytes();
    }
    else
    {
        pucDestinationData = (unsigned char*) m_imgTmpImage.Ptr();
        uiDestinationStride = m_imgTmpImage.StrideBytes();
    }
#else
    pucDestinationData = (unsigned char*) m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Ptr();
    uiDestinationStride = m_ppyrInputImages[m_iCurrentImage].GetLevel(0).StrideBytes();
#endif

	unsigned int uiSkip = 2+(m_iSubsample-1)*4;
	for(int iY=0; iY<m_iInputHeight/m_iSubsample; iY+=1)
	{
		const unsigned char *pucSourceRow = pucBuffer+m_iSubsample*iY*uiBufferStrideBytes;
		float *pfDestinationRow = ((float*) (pucDestinationData+iY*uiDestinationStride));
		float *pfDestinationRowEnd = pfDestinationRow + m_iBaseWidth;
		while(pfDestinationRow < pfDestinationRowEnd)
		{
//			*pfDestinationRow++ = float(((25 * (*pucSourceRow++) + 129 * (*pucSourceRow++) + 66 * (*pucSourceRow) + 128) >> 8) + 16);
			*pfDestinationRow++ = 0.114f*float(*pucSourceRow++) + 0.587f*float(*pucSourceRow++) + 0.299f*float(*pucSourceRow);
			pucSourceRow += uiSkip;
		}
	}
    TIMER_STOP_RGB2Y;

    // Optionally Pre-Blur the Image
    TIMER_START_PYRAMID;
	{
		int iPreBlurCount = m_iPreBlurCount;
		while(iPreBlurCount>0)
		{
			VT_HR_EXIT(FilterHoriz121(m_ppyrInputImages[m_iCurrentImage].GetLevel(0), m_imgTmpImage,
									 m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Width(),
									 m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Height()));
			VT_HR_EXIT(FilterVert121(m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(0),
									m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Width(),
									m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Height()));
			iPreBlurCount--;
		}
	}

    // Create the Pyramid
    for(int iLevel=0; iLevel<m_iLevels-1; iLevel++)
    {
        VT_HR_EXIT(FilterHoriz121(m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel), m_imgTmpImage,
                                 m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Width(),
                                 m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterVert121(m_imgTmpImage, m_imgTmpImage2,
                                m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Width(),
                                m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterHoriz121Subs(m_imgTmpImage2,
                                     m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Width(),
                                     m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterVert121Subs(m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1),
                                    m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Width(),
                                    m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Height()));
    }
    TIMER_STOP_PYRAMID;

	VT_HR_EXIT_LABEL()

    if (hr != S_OK)
    {
        Deallocate();
    }
    return hr;
}

HRESULT CVOF::AddYImage(float* pfBuffer, int iBufferWidth, int iBufferHeight, 
                         size_t uiBufferStrideBytes)
{
    VT_HR_BEGIN()

    // Check Allocated and Input Buffer is the Right Size
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iBufferWidth != m_iInputWidth || iBufferHeight != m_iInputHeight) ? E_INVALIDARG : S_OK);

    // Copy into Base Level
    TIMER_START_RGB2Y;
    m_iCurrentImage = m_iPreviousImage;
    m_iPreviousImage = 1-m_iCurrentImage;
    m_iNumberImages++;
    unsigned char *pucDestinationData;
    unsigned int uiDestinationStride;
    pucDestinationData = (unsigned char*) m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Ptr();
    uiDestinationStride = m_ppyrInputImages[m_iCurrentImage].GetLevel(0).StrideBytes();

	unsigned int uiSkip = m_iSubsample;
	for(int iY=0; iY<m_iInputHeight/m_iSubsample; iY+=1)
	{
		float *pfSourceRow = (float*) (((unsigned char*) pfBuffer)+m_iSubsample*iY*uiBufferStrideBytes);
		float *pfDestinationRow = ((float*) (pucDestinationData+iY*uiDestinationStride));
		float *pfDestinationRowEnd = pfDestinationRow + m_iBaseWidth;
		while(pfDestinationRow < pfDestinationRowEnd)
		{
			*pfDestinationRow++ = *pfSourceRow;
			pfSourceRow += uiSkip;
		}
	}
    TIMER_STOP_RGB2Y;

    // Optionally Pre-Blur the Image
    TIMER_START_PYRAMID;
	int iPreBlurCount = m_iPreBlurCount;
	while(iPreBlurCount>0)
	{
		VT_HR_EXIT(FilterHoriz121(m_ppyrInputImages[m_iCurrentImage].GetLevel(0), m_imgTmpImage,
									m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Width(),
									m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Height()));
		VT_HR_EXIT(FilterVert121(m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(0),
								m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Width(),
								m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Height()));
		iPreBlurCount--;
	}

    // Create the Pyramid
    for(int iLevel=0; iLevel<m_iLevels-1; iLevel++)
    {
        VT_HR_EXIT(FilterHoriz121(m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel), m_imgTmpImage,
                                 m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Width(),
                                 m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterVert121(m_imgTmpImage, m_imgTmpImage2,
                                m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Width(),
                                m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterHoriz121Subs(m_imgTmpImage2,
                                     m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Width(),
                                     m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterVert121Subs(m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1),
                                    m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Width(),
                                    m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Height()));
    }
    TIMER_STOP_PYRAMID;

	VT_HR_EXIT_LABEL()

    if (hr != S_OK)
    {
        Deallocate();
    }
    return hr;
}

HRESULT CVOF::AddYImage(unsigned char *pucBuffer, int iBufferWidth, int iBufferHeight, 
                         size_t uiBufferStrideBytes)
{
    VT_HR_BEGIN()

    // Check Allocated and Input Buffer is the Right Size
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iBufferWidth != m_iInputWidth || iBufferHeight != m_iInputHeight) ? E_INVALIDARG : S_OK);

    // Copy into Base Level
    TIMER_START_RGB2Y;
    m_iCurrentImage = m_iPreviousImage;
    m_iPreviousImage = 1-m_iCurrentImage;
    m_iNumberImages++;
    unsigned char *pucDestinationData;
    unsigned int uiDestinationStride;
    pucDestinationData = (unsigned char*) m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Ptr();
    uiDestinationStride = m_ppyrInputImages[m_iCurrentImage].GetLevel(0).StrideBytes();

	unsigned int uiSkip = m_iSubsample;
	for(int iY=0; iY<m_iInputHeight/m_iSubsample; iY+=1)
	{
		unsigned char *pucSourceRow = pucBuffer+m_iSubsample*iY*uiBufferStrideBytes;
		float *pfDestinationRow = ((float*) (pucDestinationData+iY*uiDestinationStride));
		float *pfDestinationRowEnd = pfDestinationRow + m_iBaseWidth;
		while(pfDestinationRow < pfDestinationRowEnd)
		{
			*pfDestinationRow++ = *pucSourceRow;
			pucSourceRow += uiSkip;
		}
	}
    TIMER_STOP_RGB2Y;

    // Optionally Pre-Blur the Image
    TIMER_START_PYRAMID;
	int iPreBlurCount = m_iPreBlurCount;
	while(iPreBlurCount>0)
	{
		VT_HR_EXIT(FilterHoriz121(m_ppyrInputImages[m_iCurrentImage].GetLevel(0), m_imgTmpImage,
									m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Width(),
									m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Height()));
		VT_HR_EXIT(FilterVert121(m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(0),
								m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Width(),
								m_ppyrInputImages[m_iCurrentImage].GetLevel(0).Height()));
		iPreBlurCount--;
	}

    // Create the Pyramid
    for(int iLevel=0; iLevel<m_iLevels-1; iLevel++)
    {
        VT_HR_EXIT(FilterHoriz121(m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel), m_imgTmpImage,
                                 m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Width(),
                                 m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterVert121(m_imgTmpImage, m_imgTmpImage2,
                                m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Width(),
                                m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterHoriz121Subs(m_imgTmpImage2,
                                     m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Width(),
                                     m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height()));
        VT_HR_EXIT(FilterVert121Subs(m_imgTmpImage, m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1),
                                    m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Width(),
                                    m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel+1).Height()));
    }
    TIMER_STOP_PYRAMID;

	VT_HR_EXIT_LABEL()

    if (hr != S_OK)
    {
        Deallocate();
    }
    return hr;
}

HRESULT CVOF::ComputeFlow(int iLowestLevel)
{
    VT_HR_BEGIN()

    // Check Inputs
    VT_HR_EXIT((iLowestLevel < 0 || iLowestLevel>=m_iLevels) ? E_INVALIDARG : S_OK);

    // Start at top level and work down to iLowestLevel
    for(int iLevel=m_iLevels-1; iLevel>=iLowestLevel; iLevel--)
    {
        TIMER_START_EXPAND;
        // Initialize the Flow
        if (iLevel==m_iLevels-1)
        {
            // Need to Zero the Flow in the Top Level -> Might have been overwritten
            VT_HR_EXIT(m_pyrFlowX.GetLevel(iLevel).Clear());
            VT_HR_EXIT(m_pyrFlowY.GetLevel(iLevel).Clear());
        }
        else
        {
            // Interpolate Flow from iLevel+1 to iLevel
            ExpandFlow(iLevel);
        }
        TIMER_STOP_EXPAND;


        // Compute flow for this level
        int iIterationLevel = MIN(m_iMaxIterationLevel, iLevel);
		if (m_piOuterIterations[iIterationLevel]!=0)
		{
			VT_HR_EXIT(ComputeFlowForLevel(iLevel, m_piOuterIterations[iIterationLevel],
					                          m_piInnerIterations[iIterationLevel]));
		}
	}

	VT_HR_END()
}

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

void CVOF::WarpImageSSE3(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, CFloatImg &imgDst, 
				         CFloatImg &imgDataTermWeight, CRect rctDst)
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
	//__m128i m128iOne = _mm_cvttps_epi32(m128One);
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
		float *fpDstRow = imgDst.Ptr(rctDst.left, rctDst.top+iY);
		float *fpFlowRowX = imgFlowX.Ptr(rctDst.left, rctDst.top+iY);
		float *fpFlowRowY = imgFlowY.Ptr(rctDst.left, rctDst.top+iY);
		float *fpDataTermWeight = imgDataTermWeight.Ptr(rctDst.left, rctDst.top+iY);

		float fLeft = float(rctDst.left);
#if (defined(_M_IX86) || defined(_M_AMD64))
		float *fpDataTermWeightEnd = fpDataTermWeight+iDstWidthTrunc;
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
			m128fSrcY = _mm_add_ps(m128fSrcY, m128iY);

			// int iDataTermWeight = 1;
			__m128i m128iDataTermWeight = _mm_cvttps_epi32(m128One);

//			__m128i m128iSrcX;
//			BilinearAddressComputationSSE4_1(m128fSrcX, m128iSrcX, m128iOne, m128iZero, m128iMaxWidth, m128iDataTermWeight);
//			__m128i m128iSrcY;
//			BilinearAddressComputationSSE4_1(m128fSrcY, m128iSrcY, m128iOne, m128iZero, m128iMaxHeight, m128iDataTermWeight);

			__m128i m128iSrcX;
			__m128i m128iSrcY;
			BilinearAddressComputationSSSE3(m128fSrcY, m128iSrcY, m128iMinusOne, m128iZero,  m128iMaxHeight, m128iDataTermWeight);
			BilinearAddressComputationSSSE3(m128fSrcX, m128iSrcX, m128iMinusOne, m128iZero, m128iMaxWidth, m128iDataTermWeight);

			// int iSrcX = int(fSrcX);
			// int iSrcY = int(fSrcY);
			// __m128i m128iSrcX = _mm_cvttps_epi32(_mm_floor_ps(m128fSrcX));
			// __m128i m128iSrcY = _mm_cvttps_epi32(_mm_floor_ps(m128fSrcY));

			// if (iSrcX < 0 || iSrcX > iSrcWidthMax || iSrcY < 0 || iSrcY > iSrcHeightMax) 
			// {
			//	iDataTermWeight = 0;
			// }
			// __m128i m128iMaxiSrcX = _mm_add_epi32(m128iSrcX, m128iOne);
			// m128iMaxiSrcX = _mm_max_epi32(m128iMaxiSrcX, m128iZero);
			// m128iDataTermWeight = _mm_sign_epi32(m128iDataTermWeight, m128iMaxiSrcX);
			// m128iMaxiSrcX = _mm_sub_epi32(m128iMaxWidth, m128iSrcX); 
			// m128iMaxiSrcX = _mm_max_epi32(m128iMaxiSrcX, m128iZero);
			// m128iDataTermWeight = _mm_sign_epi32(m128iDataTermWeight, m128iMaxiSrcX);
			// __m128i m128iMaxiSrcY = _mm_add_epi32(m128iSrcY, m128iOne);
			// m128iMaxiSrcY = _mm_max_epi32(m128iMaxiSrcY, m128iZero);
			// m128iDataTermWeight = _mm_sign_epi32(m128iDataTermWeight, m128iMaxiSrcY);
			// m128iMaxiSrcY = _mm_sub_epi32(m128iMaxHeight, m128iSrcY); 
			// m128iMaxiSrcY = _mm_max_epi32(m128iMaxiSrcY, m128iZero);
			// m128iDataTermWeight = _mm_sign_epi32(m128iDataTermWeight, m128iMaxiSrcY);

			// iSrcX *= iDataTermWeight;
			// iSrcY *= iDataTermWeight;
			// m128iSrcX = _mm_sign_epi32(m128iSrcX, m128iDataTermWeight);
			// m128iSrcY = _mm_sign_epi32(m128iSrcY, m128iDataTermWeight);

			// *fpDataTermWeight++ = float(iDataTermWeight);
			__m128 m128fDataTermWeight = _mm_cvtepi32_ps(m128iDataTermWeight);
			_mm_store_ps(fpDataTermWeight, m128fDataTermWeight);
			fpDataTermWeight+=4;

			// _mm_store_si128((__m128i*) ipiSrcX, m128iSrcX);
			// _mm_store_si128((__m128i*) ipiSrcY, m128iSrcY);

#ifndef VT_GCC
			int iSrcX0 = m128iSrcX.m128i_i32[0];
			int iSrcY0 = m128iSrcY.m128i_i32[0];
			int iSrcX1 = m128iSrcX.m128i_i32[1];
			int iSrcY1 = m128iSrcY.m128i_i32[1];
			int iSrcX2 = m128iSrcX.m128i_i32[2];
			int iSrcY2 = m128iSrcY.m128i_i32[2];
			int iSrcX3 = m128iSrcX.m128i_i32[3];
			int iSrcY3 = m128iSrcY.m128i_i32[3];
#else
			int iSrcX0 = SSE2_mm_extract_epi32(m128iSrcX, 0);
			int iSrcY0 = SSE2_mm_extract_epi32(m128iSrcY, 0);
			int iSrcX1 = SSE2_mm_extract_epi32(m128iSrcX, 1);
			int iSrcY1 = SSE2_mm_extract_epi32(m128iSrcY, 1);
			int iSrcX2 = SSE2_mm_extract_epi32(m128iSrcX, 2);
			int iSrcY2 = SSE2_mm_extract_epi32(m128iSrcY, 2);
			int iSrcX3 = SSE2_mm_extract_epi32(m128iSrcX, 3);
			int iSrcY3 = SSE2_mm_extract_epi32(m128iSrcY, 3);
#endif

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
//			_mm_store_ps(fpDstRow, m128Tmp2);
			_mm_store_ps(fpDstRow, m128Tmp2);
			fpDstRow+=4;
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

HRESULT CVOF::ExpandFlow(int iLevel)
{
    VT_HR_BEGIN()

    // Check Inputs
    VT_HR_EXIT((iLevel < 0 || iLevel>m_iLevels-2) ? E_INVALIDARG : S_OK);

#if 0
    VT_HR_EXIT(VtZoom(m_pyrFlow.GetLevel(iLevel), m_pyrFlow.GetLevel(iLevel+1)));
    VT_HR_EXIT(VtScaleImage(m_pyrFlow.GetLevel(iLevel), m_pyrFlow.GetLevel(iLevel), 2.0f));
#else
	int iWidthDown = m_pyrFlowX.GetLevel(iLevel).Width();
	int iHeightDown = m_pyrFlowX.GetLevel(iLevel).Height();
	int iWidthUp = m_pyrFlowX.GetLevel(iLevel+1).Width();
	int iHeightUp = m_pyrFlowX.GetLevel(iLevel+1).Height();
	for(int iYDown=0; iYDown<iHeightDown; iYDown++)
	{
		int iYUp = iYDown/2;
		int iYUp2 = iYDown/2 + iYDown%2;
		if (iYUp > iHeightUp-1)
		{
			iYUp = iHeightUp-1;
		}
		if (iYUp2 > iHeightUp-1)
		{
			iYUp2 = iHeightUp-1;
		}
		float *pfDstRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(iYDown);
		float *pfDstRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(iYDown);
		float *pfDstRowEndX = m_pyrFlowX.GetLevel(iLevel).Ptr(iWidthDown-1, iYDown);
		float *pfSrcRow1X = m_pyrFlowX.GetLevel(iLevel+1).Ptr(iYUp);
		float *pfSrcRow1Y = m_pyrFlowY.GetLevel(iLevel+1).Ptr(iYUp);
		float *pfSrcRow2X = m_pyrFlowX.GetLevel(iLevel+1).Ptr(iYUp2);
		float *pfSrcRow2Y = m_pyrFlowY.GetLevel(iLevel+1).Ptr(iYUp2);
		float *pfSrcRow1EndX = m_pyrFlowX.GetLevel(iLevel+1).Ptr(iWidthUp-1, iYUp);
    
		// First Columns
		float fCurr1X = *pfSrcRow1X++;
		float fCurr1Y = *pfSrcRow1Y++;
		float fCurr2X = *pfSrcRow2X++;
		float fCurr2Y = *pfSrcRow2Y++;
		*pfDstRowX++ = fCurr1X+fCurr2X;
		*pfDstRowY++ = fCurr1Y+fCurr2Y;

		// Middle Columns
		while(pfSrcRow1X<pfSrcRow1EndX)
		{
			float fNext1X = *pfSrcRow1X++;
			float fNext1Y = *pfSrcRow1Y++;
			float fNext2X = *pfSrcRow2X++;
			float fNext2Y = *pfSrcRow2Y++;
			*pfDstRowX++ = 0.5f*(fCurr1X+fCurr2X+fNext1X+fNext2X);
			*pfDstRowY++ = 0.5f*(fCurr1Y+fCurr2Y+fNext1Y+fNext2Y);
			*pfDstRowX++ = fNext1X+fNext2X;
			*pfDstRowY++ = fNext1Y+fNext2Y;
			fCurr1X = fNext1X;
			fCurr1Y = fNext1Y;
			fCurr2X = fNext2X;
			fCurr2Y = fNext2Y;
		}

		// Last Columns
		while(pfDstRowX<=pfDstRowEndX)
		{
			*pfDstRowX++ = fCurr1X+fCurr2X;
			*pfDstRowY++ = fCurr1Y+fCurr2Y;
		}
	}
#endif

	VT_HR_END()
}


HRESULT CVOF::ComputeErrorImageNoSSE(int iLevel, float &fOuterError)
{
    VT_HR_BEGIN()

    // Check Inputs
    VT_HR_EXIT((iLevel < 0 || iLevel>=m_iLevels) ? E_INVALIDARG : S_OK);

    int iWidth = m_pyrErrorImage.GetLevel(iLevel).Width();
    int iHeight = m_pyrErrorImage.GetLevel(iLevel).Height();
#ifdef OUTPUT_OUTER_ITERATION_ERROR
    fOuterError = 0.0f;
#else
    UNREFERENCED_PARAMETER(fOuterError);
#endif
    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfSrc1 = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY);
        float *pfSrc1End = pfSrc1+iWidth;
        float *pfSrc2 = m_pyrWarpedInput.GetLevel(iLevel).Ptr(iY);
        float *pfDst = m_pyrErrorImage.GetLevel(iLevel).Ptr(iY);
        while(pfSrc1<pfSrc1End)
        {
            float fError = *pfSrc1++ - *pfSrc2++;
            *pfDst++ = fError;
#ifdef OUTPUT_OUTER_ITERATION_ERROR
            fOuterError += fError*fError;
#endif
        }
    }

	VT_HR_END()
}

HRESULT CVOF::ComputeErrorImageSSE2(int iLevel, float &fOuterError)
{
    VT_HR_BEGIN()

    // Check Inputs
    VT_HR_EXIT((iLevel < 0 || iLevel>=m_iLevels) ? E_INVALIDARG : S_OK);

    int iWidth = m_pyrErrorImage.GetLevel(iLevel).Width();
#if (defined(_M_IX86) || defined(_M_AMD64))
    int iWidthTrunc = 4*(iWidth/4);
#endif
    int iHeight = m_pyrErrorImage.GetLevel(iLevel).Height();
#ifdef OUTPUT_OUTER_ITERATION_ERROR
    fOuterError = 0.0f;
#else
    UNREFERENCED_PARAMETER(fOuterError);
#endif
    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfSrc1 = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY);
#if (defined(_M_IX86) || defined(_M_AMD64))
        float *pfSrc1EndTrunc = pfSrc1+iWidthTrunc;
#endif
        float *pfSrc1End = pfSrc1+iWidth;
        float *pfSrc2 = m_pyrWarpedInput.GetLevel(iLevel).Ptr(iY);
        float *pfDst = m_pyrErrorImage.GetLevel(iLevel).Ptr(iY);
#if (defined(_M_IX86) || defined(_M_AMD64))
        while(pfSrc1<pfSrc1EndTrunc)
        {
            // float fError = *pfSrc1++ - *pfSrc2++;
            __m128 m128Src1 = _mm_load_ps(pfSrc1);
            pfSrc1 += 4;
            __m128 m128Src2 = _mm_load_ps(pfSrc2);
            pfSrc2 += 4;
            // *pfDst++ = fError;
            __m128 m128Error = _mm_sub_ps(m128Src1, m128Src2);
            _mm_store_ps(pfDst, m128Error);
            pfDst += 4;
#ifdef OUTPUT_OUTER_ITERATION_ERROR
            // fOuterError += fError*fError;
            fOuterError += m128Error.m128_f32[0] + m128Error.m128_f32[1] + 
                           m128Error.m128_f32[2] + m128Error.m128_f32[3];
#endif
        }
#endif
        while(pfSrc1<pfSrc1End)
        {
            float fError = *pfSrc1++ - *pfSrc2++;
            *pfDst++ = fError;
#ifdef OUTPUT_OUTER_ITERATION_ERROR
            fOuterError += fError*fError;
#endif
        }
    }

	VT_HR_END()
}

HRESULT CVOF::ComputeFlowForLevel(int iLevel, int iMaxOuterIterations, int iMaxInnerIterations)
{
    VT_HR_BEGIN()

	float fOuterError;
	int iOuterIterationCount;

    // Check Inputs
    VT_HR_EXIT((iLevel < 0 || iLevel>=m_iLevels) ? E_INVALIDARG : S_OK);

    // Compute X and Y Derivs of Image
    VT_HR_EXIT(ComputeDerivativesForLevel(iLevel));

    // Zero the DFlow
    TIMER_START_DFLOW;
    VT_HR_EXIT(m_pyrDFlowX.GetLevel(iLevel).Clear());
    VT_HR_EXIT(m_pyrDFlowY.GetLevel(iLevel).Clear());
    TIMER_STOP_DFLOW;

    // Outer Loop
    fOuterError = 1e10;
    iOuterIterationCount = 0;
    while(fOuterError >= 0.0f && iOuterIterationCount<iMaxOuterIterations)
    {

        // Compute Error Image
        TIMER_START_W2;
        VT_HR_EXIT(WarpImage(m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel), m_pyrFlowX.GetLevel(iLevel), m_pyrFlowY.GetLevel(iLevel),
                             m_pyrWarpedInput.GetLevel(iLevel), m_pyrDataTermWeight.GetLevel(iLevel), iLevel));
        TIMER_STOP_W2;

        TIMER_START_ERRORIM;
        VT_HR_EXIT(ComputeErrorImage(iLevel, fOuterError));
        TIMER_STOP_ERRORIM;

        TIMER_START_LAP;
        VT_HR_EXIT(ComputeFlowLaplacian(iLevel));
        TIMER_STOP_LAP;
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
        TIMER_START_WEIGHT;
		if (m_pfnWeightingFn != NULL)
		{
			VT_HR_EXIT(ComputeWeights(iLevel));
		}
        TIMER_STOP_WEIGHT;
#endif
        TIMER_START_LS;
        VT_HR_EXIT(ComputeLinearSystem(iLevel));
        TIMER_STOP_LS;

#ifdef OUTPUT_OUTER_ITERATION_ERROR
        float fSmoothnessWeight = m_pfSmoothnessParams[0];
        int iHeight = m_pyrFlowLaplacianX.GetLevel(iLevel).Height();
        int iWidth = m_pyrFlowLaplacianX.GetLevel(iLevel).Width();
        for(int iY=0; iY<iHeight; iY++)
        {
            float *pfLaplacianRowX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(iY);
            float *pfLaplacianRowY = m_pyrFlowLaplacianY.GetLevel(iLevel).Ptr(iY);
            float *pfLaplacianRowEndX = pfLaplacianRowX+iWidth;
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
            float *pfSmoothnessRow = m_pyrSmoothnessWeight.GetLevel(iLevel).Ptr(iY);
#endif
            while(pfLaplacianRowX<pfLaplacianRowEndX)
            {
                float fLaplacianX = *pfLaplacianRowX++;
                float fLaplacianY = *pfLaplacianRowY++;
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
                fSmoothnessWeight = *pfSmoothnessRow++;
#endif
                fOuterError += fSmoothnessWeight*(fLaplacianX*fLaplacianX+fLaplacianY*fLaplacianY);
            }
        }
        fOuterError /= float(iWidth*iHeight);
        fOuterError = sqrt(fOuterError);
        wprintf(L"Level: %d Outer Iteration: %d Error: %lf\n", iLevel, iOuterIterationCount, fOuterError);
#endif

		// Solve the Linear System using GS-SOR
		float fInnerError = 1e10;
		int iInnerIterationCount = 0;
		while(fInnerError >= 0.0f && iInnerIterationCount<iMaxInnerIterations)
		{
			TIMER_START_SOR;
			VT_HR_EXIT(SOR(iLevel, fInnerError));
			TIMER_STOP_SOR;

#ifdef OUTPUT_INNER_ITERATION_ERROR
			wprintf(L"Level: %d Outer Iteration: %d Inner Iteration: %d Error: %lf\n", iLevel,
					iOuterIterationCount, iInnerIterationCount, fInnerError);
#endif

			iInnerIterationCount++;
		}

        // Update the Flow using Inverse Compositional Update
        TIMER_START_DFLOW;
        // VT_HR_EXIT(VtSubImages(m_pyrFlowX.GetLevel(iLevel), m_pyrFlowX.GetLevel(iLevel), m_pyrDFlowX.GetLevel(iLevel)));
        // VT_HR_EXIT(VtSubImages(m_pyrFlowY.GetLevel(iLevel), m_pyrFlowY.GetLevel(iLevel), m_pyrDFlowY.GetLevel(iLevel)));
        UpdateDFlow(m_pyrFlowX.GetLevel(iLevel), m_pyrDFlowX.GetLevel(iLevel));
        UpdateDFlow(m_pyrFlowY.GetLevel(iLevel), m_pyrDFlowY.GetLevel(iLevel));
        TIMER_STOP_DFLOW;

        iOuterIterationCount++;
    }
            
#ifdef OUTPUT_OUTER_ITERATION_ERROR
    VT_HR_EXIT(WarpImage(m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel), m_pyrFlowX.GetLevel(iLevel), m_pyrFlowY.GetLevel(iLevel),
                         m_pyrWarpedInput.GetLevel(iLevel), m_pyrDataTermWeight.GetLevel(iLevel), iLevel));
    int iWidth = m_pyrErrorImage.GetLevel(iLevel).Width();
    int iHeight = m_pyrErrorImage.GetLevel(iLevel).Height();
    fOuterError = 0.0f;
    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfSrc1 = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY);	
        float *pfSrc1End = pfSrc1+iWidth;
        float *pfSrc2 = m_pyrWarpedInput.GetLevel(iLevel).Ptr(iY);
        float *pfDst = m_pyrErrorImage.GetLevel(iLevel).Ptr(iY);
        while(pfSrc1<pfSrc1End)
        {
            float fError = *pfSrc1++ - *pfSrc2++;	
            *pfDst++ = fError;
            fOuterError += fError*fError;
        }
    }
    float fSmoothnessWeight = m_fLambda;
    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfLaplacianRowX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(iY);
        float *pfLaplacianRowY = m_pyrFlowLaplacianY.GetLevel(iLevel).Ptr(iY);
        float *pfLaplacianRowEndX = pfLaplacianRowX+iWidth;
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
        float *pfSmoothnessRow = m_pyrSmoothnessWeight.GetLevel(iLevel).Ptr(iY);
#endif
        while(pfLaplacianRowX<pfLaplacianRowEndX)
        {
            float fLaplacianX = *pfLaplacianRowX++;
            float fLaplacianY = *pfLaplacianRowY++;
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS	
            fSmoothnessWeight = *pfSmoothnessRow++;
#endif
            fOuterError += fSmoothnessWeight*(fLaplacianX*fLaplacianX+fLaplacianY*fLaplacianY);	
        }		
    }
    fOuterError /= float(iWidth*iHeight);
    fOuterError = sqrt(fOuterError);
    wprintf(L"Level: %d Outer Iteration: %d Error: %lf\n", iLevel, iOuterIterationCount, fOuterError);
#endif

	VT_HR_END()
}

HRESULT CVOF::ComputeDerivativesForLevel(int iLevel)
{
    VT_HR_BEGIN()

    // Check Inputs
    VT_HR_EXIT((iLevel < 0 || iLevel>=m_iLevels) ? E_INVALIDARG : S_OK);
    int iWidth = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Width();
    int iHeight = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Height();
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

    TIMER_START_DERIVS;
    for(int iY=0; iY<iHeight; iY++)
    {
        float *fpSrcRow = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY);
        float *fpDstRowX = m_pyrDerivsX.GetLevel(iLevel).Ptr(iY);
        float fLeft = *fpSrcRow++;
        float fCurrent = fLeft; 
        float fRight = *fpSrcRow++;
        *fpDstRowX++ = fRight-fLeft;
        float *fpSrcRowEnd = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY)+iWidth;
        while(fpSrcRow<fpSrcRowEnd)
        {
            fLeft = fCurrent;
            fCurrent = fRight;
            fRight = *fpSrcRow++;
            *fpDstRowX++ = 0.5f*(fRight-fLeft);
        }
        *fpDstRowX = fRight-fCurrent;
    }
    for(int iY=0; iY<iHeight; iY++)
    {
        float fFactor = 0.5f;
        float *fpSrcRowAbove;
        float *fpSrcRowBelow;
        if (iY==0)
        {
            fpSrcRowAbove = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY);
            fpSrcRowBelow = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY+1);
            fFactor = 1.0f;
        }
        else if (iY==iHeight-1)
        {
            fpSrcRowAbove = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY-1);
            fpSrcRowBelow = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY);
            fFactor = 1.0f;
        }
        else
        {
            fpSrcRowAbove = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY-1);
            fpSrcRowBelow = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr(iY+1);
        }
        float *fpDstRowY = m_pyrDerivsY.GetLevel(iLevel).Ptr(iY);
        float *fpDstRowEndY = m_pyrDerivsY.GetLevel(iLevel).Ptr(iY)+ iWidth;
        while(fpDstRowY<fpDstRowEndY)
        {
            *fpDstRowY++ = fFactor*(*fpSrcRowBelow++ - *fpSrcRowAbove++);
        }
    }
    TIMER_STOP_DERIVS;

	VT_HR_END()
}

#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS

HRESULT CVOF::ComputeWeights(int iLevel)
{
    VT_HR_BEGIN()

    // Check Inputs
    VT_HR_EXIT((iLevel < 0 || iLevel>=m_iLevels) ? E_INVALIDARG : S_OK);
    VT_HR_EXIT((m_pfnWeightingFn == NULL) ? E_POINTER : S_OK);
    m_pfnWeightingFn(m_pyrSmoothnessWeight.GetLevel(iLevel), m_pfSmoothnessParams);

	VT_HR_END()
}

#endif

/*

class CTransformGraphNoSrcNode: public CTransformGraphNode
{
public:
    virtual UINT                       GetSourceCount() const
    { return 0; }
    virtual eSourceType                GetSourceType(UINT uIndex) const
    { return eReaderSource; }
    virtual const IMAGEREADER_SOURCE*  GetReaderSource(UINT uIndex) const
    { return NULL; }
    virtual const CTransformGraphNode* GetTransformSource(UINT uIndex) const
    { return NULL; }
    virtual UINT                       GetDestCount() const
    { return 0; }
public:
    CTransformGraphNoSrcNode() 
    {} 

//	CTransformGraphNoSrcNode(IImageTransform* pTransform, IImageWriter* pDst=NULL): 
 //       CTransformGraphNode(pTransform, pDst)
    CTransformGraphNoSrcNode(IImageTransform* pTransform): 
        CTransformGraphNode(pTransform)
    {}
};

class CCustomWarpTransform : public IImageTransform
{
    // IImageTransform implementation
public:
    virtual bool RequiresCloneForConcurrency()
    { return false; }

    virtual void    GetSrcPixFormat(IN OUT TRANSFORM_PIX_FRMT* pfrmtSrcs, 
                                    IN UINT  //uSrcCnt
                                    , IN const TRANSFORM_PIX_FRMT& frmtDst)
    { }

    virtual void    GetDstPixFormat(OUT TRANSFORM_PIX_FRMT& frmtDst,
                                    IN  const TRANSFORM_PIX_FRMT* pfrmtSrcs, 
                                    IN  UINT  //uSrcCnt
                                    )
    { frmtDst = TRANSFORM_PIX_FRMT(OBJ_FLOATIMG, 1); }

    virtual HRESULT GetRequiredSrcRect(OUT TRANSFORM_SOURCE_DESC* pSrcReq,
                                      OUT UINT& uSrcReqCount,
                                      IN  UINT  //uSrcCnt
                                      , 
                                      IN  const CRect& rctLayerDst
                                      )
    { 
        uSrcReqCount = 0;
        return S_OK;
    }

    virtual HRESULT GetAffectedDstRect(OUT CRect& rctDst,
                                      IN  const CRect& rctSrc,
                                      IN  UINT // uSrcIndex
                                      ,
                                      IN  UINT // uSrcCnt
                                      )
    { 
        rctDst = m_imgFlowX.Rect();
        return S_OK;
    }

    virtual HRESULT GetResultingDstRect(OUT CRect& rctDst,
                                        IN  const CRect& rctSrc,
                                        IN  UINT //uSrcIndex
                                        ,
                                        IN  UINT //uSrcCnt
                                        )
    { 
        rctDst = m_imgFlowX.Rect();
        return S_OK;
    }

    virtual HRESULT Transform(OUT CImg* pimgDstRegion, 
                              IN  const CRect& rctLayerDst,
                              IN  CImg *const *ppimgSrcRegions,
                              IN  const TRANSFORM_SOURCE_DESC* pSrcDesc,
                              IN  UINT  //uSrcCnt
                              )
    { 
        WarpImageSSE3(m_imgSrc, m_imgFlowX, m_imgFlowY, m_imgDst, m_imgDataTermWeight, rctLayerDst);
        return S_OK;
    }

public:
    HRESULT Initialize(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, CFloatImg &imgDst, CFloatImg &imgDataTermWeight)
    {
        VT_HR_BEGIN();

        VT_HR_EXIT(imgSrc.Share(m_imgSrc));
        VT_HR_EXIT(imgFlowX.Share(m_imgFlowX));
        VT_HR_EXIT(imgFlowY.Share(m_imgFlowY));
        VT_HR_EXIT(imgDst.Share(m_imgDst));
        VT_HR_EXIT(imgDataTermWeight.Share(m_imgDataTermWeight));

        VT_HR_END();
    }

    virtual HRESULT Clone(ITaskState **ppState)
    {
        return CloneTaskState<CCustomWarpTransform>(ppState, [this](CCustomWarpTransform* pN) -> HRESULT
        { 
            return pN->Initialize(m_imgSrc, m_imgFlowX, m_imgFlowY, m_imgDst, m_imgDataTermWeight);
        });  
    }

protected:
    CFloatImg m_imgSrc;
    CFloatImg m_imgFlowX;
    CFloatImg m_imgFlowY;
    CFloatImg m_imgDst;
    CFloatImg m_imgDataTermWeight;
};

HRESULT CVOF::WarpImageCustomWarp(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, CFloatImg &imgDst, 
                                  CFloatImg &imgDataTermWeight)
{
    VT_HR_BEGIN();

//	CImgReaderWriter<CImg> dst;
//	VT_HR_EXIT(imgDst.Share(dst));

    CCustomWarpTransform x; 
    VT_HR_EXIT(x.Initialize(imgSrc, imgFlowX, imgFlowY, imgDst, imgDataTermWeight));

    CTransformGraphNoSrcNode g(&x);

    // choose block size so that at least 4 blocks of work are done and blklsize 
    // is a multiple of 8
    int iBlkSize = (VtMin(imgDst.Width(), imgDst.Height()) / 2) & ~7;
    iBlkSize = VtMin(iBlkSize, 128);
    VT_ASSERT( iBlkSize != 0 );

    VT_TRANSFORM_TASK_OPTIONS opts;
    opts.blocksize.cx = opts.blocksize.cy = iBlkSize;
    VT_HR_EXIT(PushTransformTaskAndWait(imgDst.Rect(), &g, (CTaskProgress*)NULL, &opts));

    VT_HR_END();
};
*/

void CVOF::MaskDataTerm(CFloatImg &imgDataTermWeight, int iLevel)
{
    int iSrcWidth = imgDataTermWeight.Width();
    int iSrcHeight = imgDataTermWeight.Height();
    int iDstWidth = imgDataTermWeight.Width();
    int iDstHeight = imgDataTermWeight.Height();
    int iMaskWidth = MAX(2,iLevel+1);

    // Mask out Pixels close to border
    // Even if interpolation valid, derivative will be questionable
    for(int iY=0; iY<iMaskWidth; iY++)
    {
        float *fpDataTermWeight = imgDataTermWeight.Ptr(iY);
        float *fpDataTermWeightEnd = fpDataTermWeight+iDstWidth;
        while(fpDataTermWeight < fpDataTermWeightEnd)
        {
            *fpDataTermWeight++ = 0.0f;
        }
        fpDataTermWeight = imgDataTermWeight.Ptr(iSrcHeight-1-iY);
        fpDataTermWeightEnd = fpDataTermWeight+iDstWidth;
        for(int iX=0; iX<iDstWidth; iX++)
        {
            *fpDataTermWeight++ = 0.0f;
        }
    }
    for(int iY=0; iY<iDstHeight; iY++)
    {
        float *fpDataTermWeight = imgDataTermWeight.Ptr(iY);
        float *fpDataTermWeightEnd = fpDataTermWeight+iMaskWidth;
        while(fpDataTermWeight < fpDataTermWeightEnd)
        {
            *fpDataTermWeight++ = 0.0f;
        }
        fpDataTermWeight = imgDataTermWeight.Ptr(iSrcWidth-iMaskWidth, iY);
        fpDataTermWeightEnd = fpDataTermWeight+iMaskWidth;
        while(fpDataTermWeight < fpDataTermWeightEnd)
        {
            *fpDataTermWeight++ = 0.0f;
        }
    }
}

void CVOF::UpdateDFlowNoSSE(CFloatImg &imgFlow, CFloatImg &imgDFlow)
{
    int iWidth = imgFlow.Width();
    int iHeight = imgFlow.Height();

    for(int iY=0; iY<iHeight; iY++)
    {
        float *fpFlowRow = imgFlow.Ptr(iY);
        float *fpFlowRowEnd = fpFlowRow + iWidth;
        float *fpDFlowRow = imgDFlow.Ptr(iY);

        while(fpFlowRow < fpFlowRowEnd)
        {
            float fDFlow = *fpDFlowRow;
            *fpDFlowRow++ = 0.0f;
            float fFlow = *fpFlowRow;
            *fpFlowRow = fFlow - fDFlow;
            fpFlowRow++;
        }
    }
}

void CVOF::UpdateDFlowSSE2(CFloatImg &imgFlow, CFloatImg &imgDFlow)
{
    int iWidth = imgFlow.Width();
    int iHeight = imgFlow.Height();
#if (defined(_M_IX86) || defined(_M_AMD64))
    int iWidthTrunc = 4*(iWidth/4);
    __m128 m128Zero = _mm_set_ps1(0.0f);
#endif

    for(int iY=0; iY<iHeight; iY++)
    {
        float *fpFlowRow = imgFlow.Ptr(iY);
#if (defined(_M_IX86) || defined(_M_AMD64))
        float *fpFlowRowEndTrunc = fpFlowRow + iWidthTrunc;
#endif
        float *fpFlowRowEnd = fpFlowRow + iWidth;
        float *fpDFlowRow = imgDFlow.Ptr(iY);

#if (defined(_M_IX86) || defined(_M_AMD64))
        while(fpFlowRow < fpFlowRowEndTrunc)
        {
            // float fDFlow = *fpDFlowRow;
            __m128 m128DFlow = _mm_load_ps(fpDFlowRow);
            // *fpDFlowRow++ = 0.0f;
            _mm_store_ps(fpDFlowRow, m128Zero);
            fpDFlowRow += 4;
            // float fFlow = *fpFlowRow;
            __m128 m128Flow = _mm_load_ps(fpFlowRow);
            // *fpFlowRow = fFlow - fDFlow;
            m128Flow = _mm_sub_ps(m128Flow, m128DFlow);
            _mm_store_ps(fpFlowRow, m128Flow);
            // fpFlowRow++;
            fpFlowRow += 4;
        }
#endif
        while(fpFlowRow < fpFlowRowEnd)
        {
            float fDFlow = *fpDFlowRow;
            *fpDFlowRow++ = 0.0f;
            float fFlow = *fpFlowRow;
            *fpFlowRow = fFlow - fDFlow;
            fpFlowRow++;
        }
    }
}

void CVOF::WarpImageNoSSE(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, CFloatImg &imgDst, 
                          CFloatImg &imgDataTermWeight)
{
    int iSrcWidth = imgSrc.Width();
    int iSrcHeight = imgSrc.Height();
    int iSrcWidthMax = iSrcWidth-2;
    int iSrcHeightMax = iSrcHeight-2;
    int iDstWidth = imgDst.Width();
    int iDstHeight = imgDst.Height();
    unsigned char *pucSrcData = (unsigned char*) imgSrc.Ptr();
    unsigned int uiSrcStride = imgSrc.StrideBytes();

    for(int iY=0; iY<iDstHeight; iY++)
    {
        float *fpDstRow = imgDst.Ptr(iY);
        float *fpFlowRowX = imgFlowX.Ptr(iY);
        float *fpFlowRowY = imgFlowY.Ptr(iY);
        float *fpDataTermWeight = imgDataTermWeight.Ptr(iY);

        for(int iX=0; iX<iDstWidth; iX++)
        {
            float fSrcX = float(iX) + *fpFlowRowX++;
            float fSrcY = float(iY) + *fpFlowRowY++;
            int iSrcX = int(floor(fSrcX));
            int iSrcY = int(floor(fSrcY));

            if (iSrcX >= 0 && iSrcX <= iSrcWidthMax && iSrcY >= 0 && iSrcY <= iSrcHeightMax) 
            {
                // Inside Image -> Mark as Valid			
                *fpDataTermWeight++ = 1.0f;

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
                float fTotal = fUU*(fTT*fPix1+fT*fPix2)+fU*(fTT*fPix3+fT*fPix4);
                *fpDstRow++ = fTotal;
            }
            else 
            {		
                // Outside Image -> Mark as Invalid			
                *fpDataTermWeight++ = 0.0f;

                // NN Interpolation - Should be infrequent so don't worry about speed
                iSrcX = MIN(iSrcWidth-1, MAX(0, int(fSrcX+0.5f))); 
                iSrcY = MIN(iSrcHeight-1, MAX(0, int(fSrcY+0.5f))); 
                *fpDstRow++ = 0.0f;
            }
        }
    }
}

HRESULT CVOF::ComputeFlowLaplacianSSE2(int iLevel)
{
    VT_HR_BEGIN()
    
    int iWidth = m_pyrFlowX.GetLevel(iLevel).Width();
    int iHeight = m_pyrFlowX.GetLevel(iLevel).Height();
#if (defined(_M_IX86) || defined(_M_AMD64))
    __m128 m128Four = _mm_set1_ps(4.0f);
    __m128 m128One = _mm_set1_ps(1.0f);
#endif
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfFlowRightRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowRightRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowAboveRowX;
        float *pfFlowAboveRowY;
        if (iY>0)
        {
            pfFlowAboveRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY-1);
            pfFlowAboveRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY-1);
        }
        else
        {
            pfFlowAboveRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY+1);
            pfFlowAboveRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY+1);
        }
        float *pfFlowBelowRowX;
        float *pfFlowBelowRowY;
        if (iY<iHeight-1)
        {
            pfFlowBelowRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY+1);
            pfFlowBelowRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY+1);
        }
        else
        {
            pfFlowBelowRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY-1);
            pfFlowBelowRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY-1);
        }
        float *pfFlowLaplacianRowX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowLaplacianRowY = m_pyrFlowLaplacianY.GetLevel(iLevel).Ptr(0, iY);
        int iWidth1 = VtMin(iWidth-1, 4);
        float *pfFlowLaplacianRowEndX1 = pfFlowLaplacianRowX + iWidth1;
#if (defined(_M_IX86) || defined(_M_AMD64))
        int iWidth2 = 4*((iWidth-1)/4);
        float *pfFlowLaplacianRowEndX2 = pfFlowLaplacianRowX + iWidth2;
#endif
        int iWidth3 = iWidth-1;
        float *pfFlowLaplacianRowEndX3 = pfFlowLaplacianRowX + iWidth3;

        // First Column
        float fFlowCurrX = *pfFlowRightRowX++;
        float fFlowCurrY = *pfFlowRightRowY++;
        float fFlowRightX = *pfFlowRightRowX++;
        float fFlowRightY = *pfFlowRightRowY++;
        float fFlowLeftX = fFlowRightX;
        float fFlowLeftY = fFlowRightY;
        float fFlowAboveX = *pfFlowAboveRowX++;
        float fFlowAboveY = *pfFlowAboveRowY++;
        float fFlowBelowX = *pfFlowBelowRowX++;
        float fFlowBelowY = *pfFlowBelowRowY++;
        *pfFlowLaplacianRowX++ = 4.0f*fFlowCurrX - (fFlowLeftX+fFlowRightX+fFlowAboveX+fFlowBelowX);
        *pfFlowLaplacianRowY++ = 4.0f*fFlowCurrY - (fFlowLeftY+fFlowRightY+fFlowAboveY+fFlowBelowY);

        // Columns 1, 2, 3
        while(pfFlowLaplacianRowX < pfFlowLaplacianRowEndX1)
        {
            fFlowLeftX = fFlowCurrX;
            fFlowLeftY = fFlowCurrY;
            fFlowCurrX = fFlowRightX;
            fFlowCurrY = fFlowRightY;
            fFlowRightX = *pfFlowRightRowX++;
            fFlowRightY = *pfFlowRightRowY++;
            fFlowAboveX = *pfFlowAboveRowX++;
            fFlowAboveY = *pfFlowAboveRowY++;
            fFlowBelowX = *pfFlowBelowRowX++;
            fFlowBelowY = *pfFlowBelowRowY++;		
            *pfFlowLaplacianRowX++ = 4.0f*fFlowCurrX - (fFlowLeftX+fFlowRightX+fFlowAboveX+fFlowBelowX);
            *pfFlowLaplacianRowY++ = 4.0f*fFlowCurrY - (fFlowLeftY+fFlowRightY+fFlowAboveY+fFlowBelowY);		
        }

        // Main Columns
#if (defined(_M_IX86) || defined(_M_AMD64))
        __m128 m128LeftX = _mm_set_ps1(fFlowCurrX);
        __m128 m128LeftY = _mm_set_ps1(fFlowCurrY);
#endif
        pfFlowRightRowX--;
        pfFlowRightRowY--;
#if (defined(_M_IX86) || defined(_M_AMD64))
        __m128 m128CurrentX = _mm_load_ps(pfFlowRightRowX);
        __m128 m128CurrentY = _mm_load_ps(pfFlowRightRowY);
        while(pfFlowLaplacianRowX < pfFlowLaplacianRowEndX2)
        {
            // Get the Right
            pfFlowRightRowX+=4;
            pfFlowRightRowY+=4;
            __m128 m128RightX = _mm_load_ps(pfFlowRightRowX);
            __m128 m128RightY = _mm_load_ps(pfFlowRightRowY);

            // Get Above and Below
            __m128 m128AboveX = _mm_load_ps(pfFlowAboveRowX);
            __m128 m128AboveY = _mm_load_ps(pfFlowAboveRowY);
            pfFlowAboveRowX+=4;
            pfFlowAboveRowY+=4;
            __m128 m128BelowX = _mm_load_ps(pfFlowBelowRowX);
            __m128 m128BelowY = _mm_load_ps(pfFlowBelowRowY);
            pfFlowBelowRowX+=4;
            pfFlowBelowRowY+=4;

            // Shuffle for Left
            __m128 m128LXT1 = _mm_unpacklo_ps(m128LeftX, m128CurrentX);
            __m128 m128LYT1 = _mm_unpacklo_ps(m128LeftY, m128CurrentY);
            __m128 m128LXT2 = _mm_shuffle_ps(m128LXT1, m128CurrentX, _MM_SHUFFLE(2, 1, 1, 0));
            __m128 m128LYT2 = _mm_shuffle_ps(m128LYT1, m128CurrentY, _MM_SHUFFLE(2, 1, 1, 0));

            // Shuffle for Right
            __m128 m128RXT1 = _mm_move_ss(m128CurrentX, m128RightX);
            __m128 m128RYT1 = _mm_move_ss(m128CurrentY, m128RightY);
            __m128 m128RXT2 = _mm_shuffle_ps(m128RXT1, m128RXT1, _MM_SHUFFLE(0, 3, 2, 1));
            __m128 m128RYT2 = _mm_shuffle_ps(m128RYT1, m128RYT1, _MM_SHUFFLE(0, 3, 2, 1));

            // Do the Math
            __m128 m128SumX = _mm_mul_ps(m128Four, m128CurrentX);
            __m128 m128SumY = _mm_mul_ps(m128Four, m128CurrentY);
            m128SumX = _mm_sub_ps(m128SumX, m128AboveX);
            m128SumY = _mm_sub_ps(m128SumY, m128AboveY);
            m128SumX = _mm_sub_ps(m128SumX, m128BelowX);
            m128SumY = _mm_sub_ps(m128SumY, m128BelowY);
            m128SumX = _mm_sub_ps(m128SumX, m128LXT2);
            m128SumY = _mm_sub_ps(m128SumY, m128LYT2);
            m128SumX = _mm_sub_ps(m128SumX, m128RXT2);
            m128SumY = _mm_sub_ps(m128SumY, m128RYT2);

            // Store
            _mm_storeu_ps(pfFlowLaplacianRowX, m128SumX);
            _mm_storeu_ps(pfFlowLaplacianRowY, m128SumY);
            pfFlowLaplacianRowX+=4;
            pfFlowLaplacianRowY+=4;

            // Update Left and Current
            m128LeftX = _mm_shuffle_ps(m128CurrentX, m128CurrentX, _MM_SHUFFLE(3, 3, 3, 3));
            m128LeftY = _mm_shuffle_ps(m128CurrentY, m128CurrentY, _MM_SHUFFLE(3, 3, 3, 3));
            m128CurrentX = _mm_mul_ps(m128RightX, m128One);
            m128CurrentY = _mm_mul_ps(m128RightY, m128One);
        }
#endif
        pfFlowRightRowX--;
        pfFlowRightRowY--;
        fFlowCurrX = *pfFlowRightRowX++;
        fFlowCurrY = *pfFlowRightRowY++;
        fFlowRightX = *pfFlowRightRowX++;
        fFlowRightY = *pfFlowRightRowY++;

        // Last Few
        while(pfFlowLaplacianRowX < pfFlowLaplacianRowEndX3)
        {
            fFlowLeftX = fFlowCurrX;
            fFlowLeftY = fFlowCurrY;
            fFlowCurrX = fFlowRightX;
            fFlowCurrY = fFlowRightY;
            fFlowRightX = *pfFlowRightRowX++;
            fFlowRightY = *pfFlowRightRowY++;
            fFlowAboveX = *pfFlowAboveRowX++;
            fFlowAboveY = *pfFlowAboveRowY++;
            fFlowBelowX = *pfFlowBelowRowX++;
            fFlowBelowY = *pfFlowBelowRowY++;		
            *pfFlowLaplacianRowX++ = 4.0f*fFlowCurrX - (fFlowLeftX+fFlowRightX+fFlowAboveX+fFlowBelowX);
            *pfFlowLaplacianRowY++ = 4.0f*fFlowCurrY - (fFlowLeftY+fFlowRightY+fFlowAboveY+fFlowBelowY);		
        }

        // Last Column
        fFlowLeftX = fFlowCurrX;
        fFlowLeftY = fFlowCurrY;
        fFlowCurrX = fFlowRightX;
        fFlowCurrY = fFlowRightY;
        fFlowRightX = fFlowLeftX;
        fFlowRightY = fFlowLeftY;
        fFlowAboveX = *pfFlowAboveRowX++;
        fFlowAboveY = *pfFlowAboveRowY++;
        fFlowBelowX = *pfFlowBelowRowX++;
        fFlowBelowY = *pfFlowBelowRowY++;			
        *pfFlowLaplacianRowX++ = 4.0f*fFlowCurrX - (fFlowLeftX+fFlowRightX+fFlowAboveX+fFlowBelowX);
        *pfFlowLaplacianRowY++ = 4.0f*fFlowCurrY - (fFlowLeftY+fFlowRightY+fFlowAboveY+fFlowBelowY);
    }

	VT_HR_END()
}

HRESULT CVOF::ComputeFlowLaplacianNoSSE(int iLevel)
{
    VT_HR_BEGIN()
    
    int iWidth = m_pyrFlowX.GetLevel(iLevel).Width();
    int iHeight = m_pyrFlowX.GetLevel(iLevel).Height();
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfFlowRightRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowRightRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowAboveRowX;
        float *pfFlowAboveRowY;
        if (iY>0)
        {
            pfFlowAboveRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY-1);
            pfFlowAboveRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY-1);
        }
        else
        {
            pfFlowAboveRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY+1);
            pfFlowAboveRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY+1);
        }
        float *pfFlowBelowRowX;
        float *pfFlowBelowRowY;
        if (iY<iHeight-1)
        {
            pfFlowBelowRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY+1);
            pfFlowBelowRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY+1);
        }
        else
        {
            pfFlowBelowRowX = m_pyrFlowX.GetLevel(iLevel).Ptr(0, iY-1);
            pfFlowBelowRowY = m_pyrFlowY.GetLevel(iLevel).Ptr(0, iY-1);
        }
        float *pfFlowLaplacianRowX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowLaplacianRowY = m_pyrFlowLaplacianY.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowLaplacianRowEndX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(iWidth-1, iY);

        // First Column
        float fFlowCurrX = *pfFlowRightRowX++;
        float fFlowCurrY = *pfFlowRightRowY++;
        float fFlowRightX = *pfFlowRightRowX++;
        float fFlowRightY = *pfFlowRightRowY++;
        float fFlowLeftX = fFlowRightX;
        float fFlowLeftY = fFlowRightY;
        float fFlowAboveX = *pfFlowAboveRowX++;
        float fFlowAboveY = *pfFlowAboveRowY++;
        float fFlowBelowX = *pfFlowBelowRowX++;
        float fFlowBelowY = *pfFlowBelowRowY++;
        *pfFlowLaplacianRowX++ = 4.0f*fFlowCurrX - (fFlowLeftX+fFlowRightX+fFlowAboveX+fFlowBelowX);
        *pfFlowLaplacianRowY++ = 4.0f*fFlowCurrY - (fFlowLeftY+fFlowRightY+fFlowAboveY+fFlowBelowY);

        // Middle Columns
        while(pfFlowLaplacianRowX < pfFlowLaplacianRowEndX)
        {
            fFlowLeftX = fFlowCurrX;
            fFlowLeftY = fFlowCurrY;
            fFlowCurrX = fFlowRightX;
            fFlowCurrY = fFlowRightY;
            fFlowRightX = *pfFlowRightRowX++;
            fFlowRightY = *pfFlowRightRowY++;
            fFlowAboveX = *pfFlowAboveRowX++;
            fFlowAboveY = *pfFlowAboveRowY++;
            fFlowBelowX = *pfFlowBelowRowX++;
            fFlowBelowY = *pfFlowBelowRowY++;		
            *pfFlowLaplacianRowX++ = 4.0f*fFlowCurrX - (fFlowLeftX+fFlowRightX+fFlowAboveX+fFlowBelowX);
            *pfFlowLaplacianRowY++ = 4.0f*fFlowCurrY - (fFlowLeftY+fFlowRightY+fFlowAboveY+fFlowBelowY);		
        }

        // Last Column
        fFlowLeftX = fFlowCurrX;
        fFlowLeftY = fFlowCurrY;
        fFlowCurrX = fFlowRightX;
        fFlowCurrY = fFlowRightY;
        fFlowRightX = fFlowLeftX;
        fFlowRightY = fFlowLeftY;
        fFlowAboveX = *pfFlowAboveRowX++;
        fFlowAboveY = *pfFlowAboveRowY++;
        fFlowBelowX = *pfFlowBelowRowX++;
        fFlowBelowY = *pfFlowBelowRowY++;			
        *pfFlowLaplacianRowX++ = 4.0f*fFlowCurrX - (fFlowLeftX+fFlowRightX+fFlowAboveX+fFlowBelowX);
        *pfFlowLaplacianRowY++ = 4.0f*fFlowCurrY - (fFlowLeftY+fFlowRightY+fFlowAboveY+fFlowBelowY);
    }

	VT_HR_END()
}

HRESULT CVOF::ComputeLinearSystemNoSSE(int iLevel)
{
    VT_HR_BEGIN()

	float fSmoothnessWeight;
	float fTotalSmoothnessWeight;
    
    int iWidth = m_pyrFlowX.GetLevel(iLevel).Width();
    int iHeight = m_pyrFlowX.GetLevel(iLevel).Height();
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

	fSmoothnessWeight = m_pfSmoothnessParams[0];
	fTotalSmoothnessWeight = 4.0f*fSmoothnessWeight;
	for(int iY=0; iY<iHeight; iY++)
	{
		float *pfFlowLaplacianRowX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(0, iY);
		float *pfFlowLaplacianRowY = m_pyrFlowLaplacianY.GetLevel(iLevel).Ptr(0, iY);
		float *pfLinearSystemRow1 = m_pyrLinearSystem1.GetLevel(iLevel).Ptr(0, iY);
		float *pfLinearSystemRow2 = m_pyrLinearSystem2.GetLevel(iLevel).Ptr(0, iY);
		float *pfLinearSystemRow3 = m_pyrLinearSystem3.GetLevel(iLevel).Ptr(0, iY);
		float *pfLinearSystemRow4 = m_pyrLinearSystem4.GetLevel(iLevel).Ptr(0, iY);
		float *pfLinearSystemRow5 = m_pyrLinearSystem5.GetLevel(iLevel).Ptr(0, iY);
		float *pfLinearSystemRow6 = m_pyrLinearSystem6.GetLevel(iLevel).Ptr(0, iY);
		float *pfLinearSystemRowEnd1 = m_pyrLinearSystem1.GetLevel(iLevel).Ptr(iWidth-1, iY);
		float *fpDataTermWeightRow = m_pyrDataTermWeight.GetLevel(iLevel).Ptr(iY);
		float *fpImageDerivsRowX = m_pyrDerivsX.GetLevel(iLevel).Ptr(iY);
		float *fpImageDerivsRowY = m_pyrDerivsY.GetLevel(iLevel).Ptr(iY);
		float *fpErrorImageRow = m_pyrErrorImage.GetLevel(iLevel).Ptr(iY);
#ifdef OUTPUT_INNER_ITERATION_ERROR
		float *fpInnerErrorRowX = m_pyrInnerErrorFactorX.GetLevel(iLevel).Ptr(iY);
		float *fpInnerErrorRowY = m_pyrInnerErrorFactorY.GetLevel(iLevel).Ptr(iY);
#endif
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
		float *fpSmoothnessWeightRow = m_pyrSmoothnessWeight.GetLevel(iLevel).Ptr(iY);
#endif

		while(pfLinearSystemRow1 <= pfLinearSystemRowEnd1)
		{
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
			fSmoothnessWeight = *fpSmoothnessWeightRow++;
			fTotalSmoothnessWeight = 4.0f*fSmoothnessWeight;
#endif			
			float fFlowLaplacianX = *pfFlowLaplacianRowX++;
			float fFlowLaplacianY = *pfFlowLaplacianRowY++;
			float fDataTermWeight = *fpDataTermWeightRow++;
			float fIX = *fpImageDerivsRowX++;
			float fIY = *fpImageDerivsRowY++;
			float fIT = *fpErrorImageRow++;
			float fDWIX = fDataTermWeight*fIX;
			float fDWIY = fDataTermWeight*fIY;
			float fDivisorX = - m_fW/(fDWIX*fIX + fTotalSmoothnessWeight);
			*pfLinearSystemRow1++ = fDivisorX*(fDWIX*fIT - fSmoothnessWeight*fFlowLaplacianX);
			*pfLinearSystemRow2++ = fDivisorX*fDWIX*fIY;
			*pfLinearSystemRow3++ = -fDivisorX*fSmoothnessWeight;
			float fDivisorY = - m_fW/(fDWIY*fIY + fTotalSmoothnessWeight);
			*pfLinearSystemRow4++ = fDivisorY*(fDWIY*fIT - fSmoothnessWeight*fFlowLaplacianY);
			*pfLinearSystemRow5++ = fDivisorY*fDWIY*fIX;
			*pfLinearSystemRow6++ = -fDivisorY*fSmoothnessWeight;
#ifdef OUTPUT_INNER_ITERATION_ERROR
			*fpInnerErrorRowX++ = -1.0f/fDivisorX;
			*fpInnerErrorRowY++ = -1.0f/fDivisorY;
#endif		
		}
	}

	VT_HR_END()
}

HRESULT CVOF::ComputeLinearSystemSSE2(int iLevel)
{
    VT_HR_BEGIN()
    
    int iWidth = m_pyrFlowX.GetLevel(iLevel).Width();
    int iHeight = m_pyrFlowX.GetLevel(iLevel).Height();
    float fSmoothnessWeight = m_pfSmoothnessParams[0];
    float fTotalSmoothnessWeight = 4.0f*fSmoothnessWeight;
#if (defined(_M_IX86) || defined(_M_AMD64))
    float fMinusW = -m_fW;
    float fMinusOne = -1.0f;
    __m128 m128MinusW = _mm_set_ps1(fMinusW);
    __m128 m128MinusOne = _mm_set_ps1(fMinusOne);
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
    float fSmoothnessMultiplier = 4.0f;
    __m128 m128SmoothnessMultiplier = _mm_load_ps1(&fSmoothnessMultiplier);
#endif
#endif
#ifdef OUTPUT_INNER_ITERATION_ERROR
    float fOne = 1.0f;
    __m128 m128One = _mm_load_ps1(&fOne);
#endif
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfFlowLaplacianRowX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowLaplacianRowY = m_pyrFlowLaplacianY.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow1 = m_pyrLinearSystem1.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow2 = m_pyrLinearSystem2.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow3 = m_pyrLinearSystem3.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow4 = m_pyrLinearSystem4.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow5 = m_pyrLinearSystem5.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow6 = m_pyrLinearSystem6.GetLevel(iLevel).Ptr(0, iY);
#if (defined(_M_IX86) || defined(_M_AMD64))
        int iRowLengthTrunc = 4*(iWidth/4);
        float *pfLinearSystemRowEnd1 = pfLinearSystemRow1 + iRowLengthTrunc;
#endif
        float *pfLinearSystemRowEnd2 = pfLinearSystemRow1 + iWidth;
        float *fpDataTermWeightRow = m_pyrDataTermWeight.GetLevel(iLevel).Ptr(iY);
        float *fpImageDerivsRowX = m_pyrDerivsX.GetLevel(iLevel).Ptr(iY);
        float *fpImageDerivsRowY = m_pyrDerivsY.GetLevel(iLevel).Ptr(iY);
        float *fpErrorImageRow = m_pyrErrorImage.GetLevel(iLevel).Ptr(iY);
#ifdef OUTPUT_INNER_ITERATION_ERROR
        float *fpInnerErrorRowX = m_pyrInnerErrorFactorX.GetLevel(iLevel).Ptr(iY);
        float *fpInnerErrorRowY = m_pyrInnerErrorFactorY.GetLevel(iLevel).Ptr(iY);
#endif
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
        float *fpSmoothnessWeightRow = m_pyrSmoothnessWeight.GetLevel(iLevel).Ptr(iY);
#endif

#if (defined(_M_IX86) || defined(_M_AMD64))
        __m128 m128SmoothnessWeight = _mm_load_ps1(&fSmoothnessWeight);
        __m128 m128TotalSmoothnessWeight = _mm_load_ps1(&fTotalSmoothnessWeight);
        while(pfLinearSystemRow1 < pfLinearSystemRowEnd1)
        {
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
            m128SmoothnessWeight = _mm_load_ps(fpSmoothnessWeightRow);
            fpSmoothnessWeightRow += 4;
            m128TotalSmoothnessWeight = _mm_mul_ps(m128SmoothnessWeight, m128SmoothnessMultiplier);
#endif			
            // float fFlowLaplacianX = *pfFlowLaplacianRowX++;
            __m128 m128fFlowLaplacianX = _mm_load_ps(pfFlowLaplacianRowX);
            pfFlowLaplacianRowX += 4;
            // float fDataTermWeight = *fpDataTermWeightRow++;
            __m128 m128fDataTermWeight = _mm_load_ps(fpDataTermWeightRow);
            fpDataTermWeightRow += 4;
            // float fIX = *fpImageDerivsRowX++;
            __m128 m128fIX = _mm_load_ps(fpImageDerivsRowX);
            fpImageDerivsRowX += 4;
            // float fIY = *fpImageDerivsRowY++;
            __m128 m128fIY = _mm_load_ps(fpImageDerivsRowY);
            fpImageDerivsRowY += 4;
            // float fIT = *fpErrorImageRow++;
            __m128 m128fIT = _mm_load_ps(fpErrorImageRow);
            fpErrorImageRow += 4;
            // float fDWIX = fDataTermWeight*fIX;
            __m128 m128DWIX = _mm_mul_ps(m128fDataTermWeight, m128fIX);
            // float fDivisorX = - m_fW/(fDWIX*fIX + fTotalSmoothnessWeight);
            __m128 m128DWIXIX = _mm_mul_ps(m128DWIX, m128fIX);
            __m128 m128DWIXIXPTS = _mm_add_ps(m128DWIXIX, m128TotalSmoothnessWeight);
            __m128 m128DivisorX = _mm_div_ps(m128MinusW, m128DWIXIXPTS);
            // *pfLinearSystemRow1++ = fDivisorX*(fDWIX*fIT - fSmoothnessWeight*fFlowLaplacianX);
            __m128 m128DWIXIT = _mm_mul_ps(m128DWIX, m128fIT);
            __m128 m128SWFLX = _mm_mul_ps(m128SmoothnessWeight, m128fFlowLaplacianX);
            __m128 m128DWIXITSubSWFLX = _mm_sub_ps(m128DWIXIT, m128SWFLX);
            __m128 m128LS1 = _mm_mul_ps(m128DivisorX, m128DWIXITSubSWFLX);
            _mm_store_ps(pfLinearSystemRow1, m128LS1);
            pfLinearSystemRow1+=4;
            // *pfLinearSystemRow2++ = fDivisorX*fDWIX*fIY;
            __m128 m128DWIXIY = _mm_mul_ps(m128DWIX, m128fIY);
            __m128 m128LS2 = _mm_mul_ps(m128DivisorX, m128DWIXIY);
            _mm_store_ps(pfLinearSystemRow2, m128LS2);
            pfLinearSystemRow2+=4;
            // *pfLinearSystemRow3++ = -fDivisorX*fSmoothnessWeight;
            __m128 m128MinusDivisorX = _mm_mul_ps(m128MinusOne, m128DivisorX);
            __m128 m128LS3 = _mm_mul_ps(m128MinusDivisorX, m128SmoothnessWeight);
            _mm_store_ps(pfLinearSystemRow3, m128LS3);
            pfLinearSystemRow3+=4;

            // float fDWIY = fDataTermWeight*fIY;
            __m128 m128DWIY = _mm_mul_ps(m128fDataTermWeight, m128fIY);
            // float fDivisorY = - m_fW/(fDWIY*fIY + fTotalSmoothnessWeight);
            __m128 m128DWIYIY = _mm_mul_ps(m128DWIY, m128fIY);
            __m128 m128DWIYIYPTS = _mm_add_ps(m128DWIYIY, m128TotalSmoothnessWeight);
            __m128 m128DivisorY = _mm_div_ps(m128MinusW, m128DWIYIYPTS);
            // float fFlowLaplacianY = *pfFlowLaplacianRowY++;
            __m128 m128fFlowLaplacianY = _mm_load_ps(pfFlowLaplacianRowY);
            pfFlowLaplacianRowY += 4;
            // *pfLinearSystemRow4++ = fDivisorY*(fDWIY*fIT - fSmoothnessWeight*fFlowLaplacianY);
            __m128 m128DWIYIT = _mm_mul_ps(m128DWIY, m128fIT);
            __m128 m128SWFLY = _mm_mul_ps(m128SmoothnessWeight, m128fFlowLaplacianY);
            __m128 m128DWIYITSubSWFLY = _mm_sub_ps(m128DWIYIT, m128SWFLY);
            __m128 m128LS4 = _mm_mul_ps(m128DivisorY, m128DWIYITSubSWFLY);
            _mm_store_ps(pfLinearSystemRow4, m128LS4);
            pfLinearSystemRow4+=4;
            // *pfLinearSystemRow5++ = fDivisorY*fDWIY*fIX;
            // __m128 m128DWIYIX = _mm_mul_ps(m128DWIY, m128fIX);
            __m128 m128LS5 = _mm_mul_ps(m128DivisorY, m128DWIXIY);
            _mm_store_ps(pfLinearSystemRow5, m128LS5);
            pfLinearSystemRow5+=4;
            // *pfLinearSystemRow6++ = -fDivisorY*fSmoothnessWeight;
            __m128 m128MinusDivisorY = _mm_mul_ps(m128MinusOne, m128DivisorY);
            __m128 m128LS6 = _mm_mul_ps(m128MinusDivisorY, m128SmoothnessWeight);
            _mm_store_ps(pfLinearSystemRow6, m128LS6);
            pfLinearSystemRow6+=4;
#ifdef OUTPUT_INNER_ITERATION_ERROR
            // *fpInnerErrorRowX++ = -1.0f/fDivisorX;
            __m128 m128IERX = _mm_div_ps(m128One, m128DivisorX);
            _mm_store_ps(fpInnerErrorRowX, m128IERX);
            fpInnerErrorRowX += 4;
            // *fpInnerErrorRowY++ = -1.0f/fDivisorY;
            __m128 m128IERY = _mm_div_ps(m128One, m128DivisorY);
            _mm_store_ps(fpInnerErrorRowY, m128IERY);
            fpInnerErrorRowY += 4;
#endif		
        }
#endif

        while(pfLinearSystemRow1 < pfLinearSystemRowEnd2)
        {
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
            fSmoothnessWeight = *fpSmoothnessWeightRow++;
            fTotalSmoothnessWeight = 4.0f*fSmoothnessWeight;
#endif			
            float fFlowLaplacianX = *pfFlowLaplacianRowX++;
            float fDataTermWeight = *fpDataTermWeightRow++;
            float fIX = *fpImageDerivsRowX++;
            float fIY = *fpImageDerivsRowY++;
            float fIT = *fpErrorImageRow++;
            float fDWIX = fDataTermWeight*fIX;
            float fDivisorX = - m_fW/(fDWIX*fIX + fTotalSmoothnessWeight);
            *pfLinearSystemRow1++ = fDivisorX*(fDWIX*fIT - fSmoothnessWeight*fFlowLaplacianX);
            *pfLinearSystemRow2++ = fDivisorX*fDWIX*fIY;
            *pfLinearSystemRow3++ = -fDivisorX*fSmoothnessWeight;
            float fDWIY = fDataTermWeight*fIY;
            float fDivisorY = - m_fW/(fDWIY*fIY + fTotalSmoothnessWeight);
            float fFlowLaplacianY = *pfFlowLaplacianRowY++;
            *pfLinearSystemRow4++ = fDivisorY*(fDWIY*fIT - fSmoothnessWeight*fFlowLaplacianY);
            *pfLinearSystemRow5++ = fDivisorY*fDWIY*fIX;
            *pfLinearSystemRow6++ = -fDivisorY*fSmoothnessWeight;
#ifdef OUTPUT_INNER_ITERATION_ERROR
            *fpInnerErrorRowX++ = -1.0f/fDivisorX;
            *fpInnerErrorRowY++ = -1.0f/fDivisorY;
#endif		
        }
    }

	VT_HR_END()
}

HRESULT CVOF::ComputeLinearSystemSSE2NoWeight(int iLevel)
{
    VT_HR_BEGIN()
    
    int iWidth = m_pyrFlowX.GetLevel(iLevel).Width();
    int iHeight = m_pyrFlowX.GetLevel(iLevel).Height();
    float fSmoothnessWeight = m_pfSmoothnessParams[0];
    float fTotalSmoothnessWeight = 4.0f*fSmoothnessWeight;
#if (defined(_M_IX86) || defined(_M_AMD64))
    float fMinusW = -m_fW;
    float fMinusOne = -1.0f;
    __m128 m128MinusW = _mm_set_ps1(fMinusW);
    __m128 m128MinusOne = _mm_set_ps1(fMinusOne);
#ifdef OUTPUT_INNER_ITERATION_ERROR
    float fOne = 1.0f;
    __m128 m128One = _mm_load_ps1(&fOne);
#endif
#endif
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfFlowLaplacianRowX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowLaplacianRowY = m_pyrFlowLaplacianY.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow1 = m_pyrLinearSystem1.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow2 = m_pyrLinearSystem2.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow3 = m_pyrLinearSystem3.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow4 = m_pyrLinearSystem4.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow5 = m_pyrLinearSystem5.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow6 = m_pyrLinearSystem6.GetLevel(iLevel).Ptr(0, iY);
#if (defined(_M_IX86) || defined(_M_AMD64))
        int iRowLengthTrunc = 4*(iWidth/4);
        float *pfLinearSystemRowEnd1 = pfLinearSystemRow1 + iRowLengthTrunc;
#endif
        float *pfLinearSystemRowEnd2 = pfLinearSystemRow1 + iWidth;
        float *fpDataTermWeightRow = m_pyrDataTermWeight.GetLevel(iLevel).Ptr(iY);
        float *fpImageDerivsRowX = m_pyrDerivsX.GetLevel(iLevel).Ptr(iY);
        float *fpImageDerivsRowY = m_pyrDerivsY.GetLevel(iLevel).Ptr(iY);
        float *fpErrorImageRow = m_pyrErrorImage.GetLevel(iLevel).Ptr(iY);
#ifdef OUTPUT_INNER_ITERATION_ERROR
        float *fpInnerErrorRowX = m_pyrInnerErrorFactorX.GetLevel(iLevel).Ptr(iY);
        float *fpInnerErrorRowY = m_pyrInnerErrorFactorY.GetLevel(iLevel).Ptr(iY);
#endif

#if (defined(_M_IX86) || defined(_M_AMD64))
        __m128 m128SmoothnessWeight = _mm_load_ps1(&fSmoothnessWeight);
        __m128 m128TotalSmoothnessWeight = _mm_load_ps1(&fTotalSmoothnessWeight);
        while(pfLinearSystemRow1 < pfLinearSystemRowEnd1)
        {
			// float fFlowLaplacianX = *pfFlowLaplacianRowX++;
            __m128 m128fFlowLaplacianX = _mm_load_ps(pfFlowLaplacianRowX);
            pfFlowLaplacianRowX += 4;
            // float fDataTermWeight = *fpDataTermWeightRow++;
            __m128 m128fDataTermWeight = _mm_load_ps(fpDataTermWeightRow);
            fpDataTermWeightRow += 4;
            // float fIX = *fpImageDerivsRowX++;
            __m128 m128fIX = _mm_load_ps(fpImageDerivsRowX);
            fpImageDerivsRowX += 4;
            // float fIY = *fpImageDerivsRowY++;
            __m128 m128fIY = _mm_load_ps(fpImageDerivsRowY);
            fpImageDerivsRowY += 4;
            // float fIT = *fpErrorImageRow++;
            __m128 m128fIT = _mm_load_ps(fpErrorImageRow);
            fpErrorImageRow += 4;
            // float fDWIX = fDataTermWeight*fIX;
            __m128 m128DWIX = _mm_mul_ps(m128fDataTermWeight, m128fIX);
            // float fDivisorX = - m_fW/(fDWIX*fIX + fTotalSmoothnessWeight);
            __m128 m128DWIXIX = _mm_mul_ps(m128DWIX, m128fIX);
            __m128 m128DWIXIXPTS = _mm_add_ps(m128DWIXIX, m128TotalSmoothnessWeight);
            __m128 m128DivisorX = _mm_div_ps(m128MinusW, m128DWIXIXPTS);
            // *pfLinearSystemRow1++ = fDivisorX*(fDWIX*fIT - fSmoothnessWeight*fFlowLaplacianX);
            __m128 m128DWIXIT = _mm_mul_ps(m128DWIX, m128fIT);
            __m128 m128SWFLX = _mm_mul_ps(m128SmoothnessWeight, m128fFlowLaplacianX);
            __m128 m128DWIXITSubSWFLX = _mm_sub_ps(m128DWIXIT, m128SWFLX);
            __m128 m128LS1 = _mm_mul_ps(m128DivisorX, m128DWIXITSubSWFLX);
            _mm_store_ps(pfLinearSystemRow1, m128LS1);
            pfLinearSystemRow1+=4;
            // *pfLinearSystemRow2++ = fDivisorX*fDWIX*fIY;
            __m128 m128DWIXIY = _mm_mul_ps(m128DWIX, m128fIY);
            __m128 m128LS2 = _mm_mul_ps(m128DivisorX, m128DWIXIY);
            _mm_store_ps(pfLinearSystemRow2, m128LS2);
            pfLinearSystemRow2+=4;
            // *pfLinearSystemRow3++ = -fDivisorX*fSmoothnessWeight;
            __m128 m128MinusDivisorX = _mm_mul_ps(m128MinusOne, m128DivisorX);
            __m128 m128LS3 = _mm_mul_ps(m128MinusDivisorX, m128SmoothnessWeight);
            _mm_store_ps(pfLinearSystemRow3, m128LS3);
            pfLinearSystemRow3+=4;

            // float fDWIY = fDataTermWeight*fIY;
            __m128 m128DWIY = _mm_mul_ps(m128fDataTermWeight, m128fIY);
            // float fDivisorY = - m_fW/(fDWIY*fIY + fTotalSmoothnessWeight);
            __m128 m128DWIYIY = _mm_mul_ps(m128DWIY, m128fIY);
            __m128 m128DWIYIYPTS = _mm_add_ps(m128DWIYIY, m128TotalSmoothnessWeight);
            __m128 m128DivisorY = _mm_div_ps(m128MinusW, m128DWIYIYPTS);
            // float fFlowLaplacianY = *pfFlowLaplacianRowY++;
            __m128 m128fFlowLaplacianY = _mm_load_ps(pfFlowLaplacianRowY);
            pfFlowLaplacianRowY += 4;
            // *pfLinearSystemRow4++ = fDivisorY*(fDWIY*fIT - fSmoothnessWeight*fFlowLaplacianY);
            __m128 m128DWIYIT = _mm_mul_ps(m128DWIY, m128fIT);
            __m128 m128SWFLY = _mm_mul_ps(m128SmoothnessWeight, m128fFlowLaplacianY);
            __m128 m128DWIYITSubSWFLY = _mm_sub_ps(m128DWIYIT, m128SWFLY);
            __m128 m128LS4 = _mm_mul_ps(m128DivisorY, m128DWIYITSubSWFLY);
            _mm_store_ps(pfLinearSystemRow4, m128LS4);
            pfLinearSystemRow4+=4;
            // *pfLinearSystemRow5++ = fDivisorY*fDWIY*fIX;
            // __m128 m128DWIYIX = _mm_mul_ps(m128DWIY, m128fIX);
            __m128 m128LS5 = _mm_mul_ps(m128DivisorY, m128DWIXIY);
            _mm_store_ps(pfLinearSystemRow5, m128LS5);
            pfLinearSystemRow5+=4;
            // *pfLinearSystemRow6++ = -fDivisorY*fSmoothnessWeight;
            __m128 m128MinusDivisorY = _mm_mul_ps(m128MinusOne, m128DivisorY);
            __m128 m128LS6 = _mm_mul_ps(m128MinusDivisorY, m128SmoothnessWeight);
            _mm_store_ps(pfLinearSystemRow6, m128LS6);
            pfLinearSystemRow6+=4;
#ifdef OUTPUT_INNER_ITERATION_ERROR
            // *fpInnerErrorRowX++ = -1.0f/fDivisorX;
            __m128 m128IERX = _mm_div_ps(m128One, m128DivisorX);
            _mm_store_ps(fpInnerErrorRowX, m128IERX);
            fpInnerErrorRowX += 4;
            // *fpInnerErrorRowY++ = -1.0f/fDivisorY;
            __m128 m128IERY = _mm_div_ps(m128One, m128DivisorY);
            _mm_store_ps(fpInnerErrorRowY, m128IERY);
            fpInnerErrorRowY += 4;
#endif		
        }
#endif

        while(pfLinearSystemRow1 < pfLinearSystemRowEnd2)
        {
            float fFlowLaplacianX = *pfFlowLaplacianRowX++;
            float fDataTermWeight = *fpDataTermWeightRow++;
            float fIX = *fpImageDerivsRowX++;
            float fIY = *fpImageDerivsRowY++;
            float fIT = *fpErrorImageRow++;
            float fDWIX = fDataTermWeight*fIX;
            float fDivisorX = - m_fW/(fDWIX*fIX + fTotalSmoothnessWeight);
            *pfLinearSystemRow1++ = fDivisorX*(fDWIX*fIT - fSmoothnessWeight*fFlowLaplacianX);
            *pfLinearSystemRow2++ = fDivisorX*fDWIX*fIY;
            *pfLinearSystemRow3++ = -fDivisorX*fSmoothnessWeight;
            float fDWIY = fDataTermWeight*fIY;
            float fDivisorY = - m_fW/(fDWIY*fIY + fTotalSmoothnessWeight);
            float fFlowLaplacianY = *pfFlowLaplacianRowY++;
            *pfLinearSystemRow4++ = fDivisorY*(fDWIY*fIT - fSmoothnessWeight*fFlowLaplacianY);
            *pfLinearSystemRow5++ = fDivisorY*fDWIY*fIX;
            *pfLinearSystemRow6++ = -fDivisorY*fSmoothnessWeight;
#ifdef OUTPUT_INNER_ITERATION_ERROR
            *fpInnerErrorRowX++ = -1.0f/fDivisorX;
            *fpInnerErrorRowY++ = -1.0f/fDivisorY;
#endif		
        }
    }

	VT_HR_END()
}

HRESULT CVOF::ComputeLinearSystemNoSSENoWeight(int iLevel)
{
    VT_HR_BEGIN()

	float fSmoothnessWeight;
	float fTotalSmoothnessWeight;
    
    int iWidth = m_pyrFlowX.GetLevel(iLevel).Width();
    int iHeight = m_pyrFlowX.GetLevel(iLevel).Height();
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

    fSmoothnessWeight = m_pfSmoothnessParams[0];
    fTotalSmoothnessWeight = 4.0f*fSmoothnessWeight;
    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfFlowLaplacianRowX = m_pyrFlowLaplacianX.GetLevel(iLevel).Ptr(0, iY);
        float *pfFlowLaplacianRowY = m_pyrFlowLaplacianY.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow1 = m_pyrLinearSystem1.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow2 = m_pyrLinearSystem2.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow3 = m_pyrLinearSystem3.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow4 = m_pyrLinearSystem4.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow5 = m_pyrLinearSystem5.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow6 = m_pyrLinearSystem6.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRowEnd1 = m_pyrLinearSystem1.GetLevel(iLevel).Ptr(iWidth-1, iY);
        float *fpDataTermWeightRow = m_pyrDataTermWeight.GetLevel(iLevel).Ptr(iY);
        float *fpImageDerivsRowX = m_pyrDerivsX.GetLevel(iLevel).Ptr(iY);
        float *fpImageDerivsRowY = m_pyrDerivsY.GetLevel(iLevel).Ptr(iY);
        float *fpErrorImageRow = m_pyrErrorImage.GetLevel(iLevel).Ptr(iY);
#ifdef OUTPUT_INNER_ITERATION_ERROR
        float *fpInnerErrorRowX = m_pyrInnerErrorFactorX.GetLevel(iLevel).Ptr(iY);
        float *fpInnerErrorRowY = m_pyrInnerErrorFactorY.GetLevel(iLevel).Ptr(iY);
#endif

        while(pfLinearSystemRow1 <= pfLinearSystemRowEnd1)
        {
            float fFlowLaplacianX = *pfFlowLaplacianRowX++;
            float fFlowLaplacianY = *pfFlowLaplacianRowY++;
            float fDataTermWeight = *fpDataTermWeightRow++;
            float fIX = *fpImageDerivsRowX++;
            float fIY = *fpImageDerivsRowY++;
            float fIT = *fpErrorImageRow++;
            float fDWIX = fDataTermWeight*fIX;
            float fDWIY = fDataTermWeight*fIY;
            float fDivisorX = - m_fW/(fDWIX*fIX + fTotalSmoothnessWeight);
            *pfLinearSystemRow1++ = fDivisorX*(fDWIX*fIT - fSmoothnessWeight*fFlowLaplacianX);
            *pfLinearSystemRow2++ = fDivisorX*fDWIX*fIY;
            *pfLinearSystemRow3++ = -fDivisorX*fSmoothnessWeight;
            float fDivisorY = - m_fW/(fDWIY*fIY + fTotalSmoothnessWeight);
            *pfLinearSystemRow4++ = fDivisorY*(fDWIY*fIT - fSmoothnessWeight*fFlowLaplacianY);
            *pfLinearSystemRow5++ = fDivisorY*fDWIY*fIX;
            *pfLinearSystemRow6++ = -fDivisorY*fSmoothnessWeight;
#ifdef OUTPUT_INNER_ITERATION_ERROR
            *fpInnerErrorRowX++ = -1.0f/fDivisorX;
            *fpInnerErrorRowY++ = -1.0f/fDivisorY;
#endif		
        }
    }

	VT_HR_END()
}

HRESULT CVOF::SORNoSSE(int iLevel, float &fTotalError)
{
    VT_HR_BEGIN()
    
    int iWidth = m_pyrDFlowX.GetLevel(iLevel).Width();
    int iHeight = m_pyrDFlowX.GetLevel(iLevel).Height();
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

    fTotalError = 0.0f;
    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfDFlowRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY);
        float *pfDFlowRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY);
        float *pfDFlowRowEndX = m_pyrDFlowX.GetLevel(iLevel).Ptr(iWidth-1, iY);
        float *pfDFlowRightRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(1, iY);
        float *pfDFlowRightRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(1, iY);
        float *pfLinearSystemRow1 = m_pyrLinearSystem1.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow2 = m_pyrLinearSystem2.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow3 = m_pyrLinearSystem3.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow4 = m_pyrLinearSystem4.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow5 = m_pyrLinearSystem5.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow6 = m_pyrLinearSystem6.GetLevel(iLevel).Ptr(0, iY);
        float *pfDFlowAboveRowX;
        float *pfDFlowAboveRowY;
        if (iY>0)
        {
            pfDFlowAboveRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY-1);
            pfDFlowAboveRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY-1);
        }
        else
        {
            pfDFlowAboveRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY+1);
            pfDFlowAboveRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY+1);
        }
        float *pfDFlowBelowRowX;
        float *pfDFlowBelowRowY;
        if (iY<iHeight-1)
        {
            pfDFlowBelowRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY+1);
            pfDFlowBelowRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY+1);
        }
        else
        {
            pfDFlowBelowRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY-1);
            pfDFlowBelowRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY-1);
        }
#ifdef OUTPUT_INNER_ITERATION_ERROR
        float *fpInnerErrorRowX = m_pyrInnerErrorFactorX.GetLevel(iLevel).Ptr(iY);
        float *fpInnerErrorRowY = m_pyrInnerErrorFactorY.GetLevel(iLevel).Ptr(iY);
#endif

        // First Column
        float fDFlowCurrX = *pfDFlowRowX;
        float fDFlowCurrY = *pfDFlowRowY;
        float fDFlowRightX = *pfDFlowRightRowX++;
        float fDFlowRightY = *pfDFlowRightRowY++;
        float fDFlowLeftX = fDFlowRightX;
        float fDFlowLeftY = fDFlowRightY;
        float fDFlowAboveX = *pfDFlowAboveRowX++;
        float fDFlowAboveY = *pfDFlowAboveRowY++;
        float fDFlowBelowX = *pfDFlowBelowRowX++;
        float fDFlowBelowY = *pfDFlowBelowRowY++;
        float fLS1X = *pfLinearSystemRow1++;
        float fLS2X = *pfLinearSystemRow2++;
        float fLS3X = *pfLinearSystemRow3++;
        float fLS1Y = *pfLinearSystemRow4++;
        float fLS2Y = *pfLinearSystemRow5++;
        float fLS3Y = *pfLinearSystemRow6++;
        float fTotalX = fLS1X + fLS2X*fDFlowCurrY + fLS3X*(fDFlowLeftX+fDFlowRightX+fDFlowAboveX+fDFlowBelowX);
        *pfDFlowRowX++ = fDFlowLeftX = fTotalX+m_fOneMinusW*fDFlowCurrX;
        // Note: Must use new value of fDFlowLeftX !!
        float fTotalY = fLS1Y + fLS2Y*fDFlowLeftX + fLS3Y*(fDFlowLeftY+fDFlowRightY+fDFlowAboveY+fDFlowBelowY);
        *pfDFlowRowY++ = fDFlowLeftY = fTotalY+m_fOneMinusW*fDFlowCurrY;

#ifdef OUTPUT_INNER_ITERATION_ERROR
        float fErrorX = (*fpInnerErrorRowX++)*(fDFlowCurrX - fDFlowLeftX);
        float fErrorY = (*fpInnerErrorRowY++)*(fDFlowCurrY - fDFlowLeftY);
        fTotalError += fErrorX*fErrorX + fErrorY*fErrorY;
#endif

        // Middle Columns
        while(pfDFlowRowX < pfDFlowRowEndX)
        {
            fDFlowCurrX = fDFlowRightX;
            fDFlowCurrY = fDFlowRightY;
            fDFlowRightX = *pfDFlowRightRowX++;
            fDFlowRightY = *pfDFlowRightRowY++;
            fDFlowAboveX = *pfDFlowAboveRowX++;
            fDFlowAboveY = *pfDFlowAboveRowY++;
            fDFlowBelowX = *pfDFlowBelowRowX++;
            fDFlowBelowY = *pfDFlowBelowRowY++;
            fLS1X = *pfLinearSystemRow1++;
            fLS2X = *pfLinearSystemRow2++;
            fLS3X = *pfLinearSystemRow3++;
            fLS1Y = *pfLinearSystemRow4++;
            fLS2Y = *pfLinearSystemRow5++;
            fLS3Y = *pfLinearSystemRow6++;
            fTotalX = fLS1X + fLS2X*fDFlowCurrY + fLS3X*(fDFlowLeftX+fDFlowRightX+fDFlowAboveX+fDFlowBelowX);
            *pfDFlowRowX++ = fDFlowLeftX = fTotalX+m_fOneMinusW*fDFlowCurrX;
            // Note: Must use new value of fDFlowLeftX !!
            fTotalY = fLS1Y + fLS2Y*fDFlowLeftX + fLS3Y*(fDFlowLeftY+fDFlowRightY+fDFlowAboveY+fDFlowBelowY);
            *pfDFlowRowY++ = fDFlowLeftY = fTotalY+m_fOneMinusW*fDFlowCurrY;

#ifdef OUTPUT_INNER_ITERATION_ERROR
            fErrorX = (*fpInnerErrorRowX++)*(fDFlowCurrX - fDFlowLeftX);
            fErrorY = (*fpInnerErrorRowY++)*(fDFlowCurrY - fDFlowLeftY);
            fTotalError += fErrorX*fErrorX + fErrorY*fErrorY;
#endif
        }

        // Last Column
        fDFlowCurrX = fDFlowRightX;
        fDFlowCurrY = fDFlowRightY;
        fDFlowRightX = fDFlowLeftX;
        fDFlowRightY = fDFlowLeftY;
        fDFlowAboveX = *pfDFlowAboveRowX++;
        fDFlowAboveY = *pfDFlowAboveRowY++;
        fDFlowBelowX = *pfDFlowBelowRowX++;
        fDFlowBelowY = *pfDFlowBelowRowY++;
        fLS1X = *pfLinearSystemRow1++;
        fLS2X = *pfLinearSystemRow2++;
        fLS3X = *pfLinearSystemRow3++;
        fLS1Y = *pfLinearSystemRow4++;
        fLS2Y = *pfLinearSystemRow5++;
        fLS3Y = *pfLinearSystemRow6++;
        fTotalX = fLS1X + fLS2X*fDFlowCurrY + fLS3X*(fDFlowLeftX+fDFlowRightX+fDFlowAboveX+fDFlowBelowX);
        *pfDFlowRowX++ = fDFlowLeftX = fTotalX+m_fOneMinusW*fDFlowCurrX;
        // Note: Must use new value of fDFlowLeftX !!
        fTotalY = fLS1Y + fLS2Y*fDFlowLeftX + fLS3Y*(fDFlowLeftY+fDFlowRightY+fDFlowAboveY+fDFlowBelowY);
        *pfDFlowRowY++ = fTotalY+m_fOneMinusW*fDFlowCurrY;

#ifdef OUTPUT_INNER_ITERATION_ERROR
        fErrorX = (*fpInnerErrorRowX++)*(fDFlowCurrX - fDFlowLeftX);
        fErrorY = (*fpInnerErrorRowY++)*(fDFlowCurrY - fDFlowLeftY);
        fTotalError += fErrorX*fErrorX + fErrorY*fErrorY;
#endif	
    }

#ifdef OUTPUT_INNER_ITERATION_ERROR			
    fTotalError = sqrt(fTotalError / (iWidth*iHeight));
#endif	

	VT_HR_END()
}

HRESULT CVOF::SORSSE2(int iLevel, float &fTotalError)
{
    VT_HR_BEGIN()
    
    int iWidth = m_pyrDFlowX.GetLevel(iLevel).Width();
    int iHeight = m_pyrDFlowX.GetLevel(iLevel).Height();
#if (defined(_M_IX86) || defined(_M_AMD64))
    __m128 m128OneMinusW = _mm_set1_ps(m_fOneMinusW);
    __m128 m128One = _mm_set1_ps(1.0f);
#endif
    VT_HR_EXIT((iWidth<2 || iHeight<2) ? E_INVALIDARG : S_OK);

    fTotalError = 0.0f;
    for(int iY=0; iY<iHeight; iY++)
    {
        float *pfDFlowRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY);
        float *pfDFlowRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY);
        int iWidth1 = VtMin(iWidth-1, 4);
        float *pfDFlowRowEndX1 = pfDFlowRowX + iWidth1;
        int iWidth2 = 4*((iWidth-1)/4);
        int iWidth3 = iWidth-1;
        float *pfDFlowRowEndX3 = pfDFlowRowX + iWidth3;
        float *pfDFlowRightRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(1, iY);
        float *pfDFlowRightRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(1, iY);
        float *pfLinearSystemRow1 = m_pyrLinearSystem1.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow2 = m_pyrLinearSystem2.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow3 = m_pyrLinearSystem3.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow4 = m_pyrLinearSystem4.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow5 = m_pyrLinearSystem5.GetLevel(iLevel).Ptr(0, iY);
        float *pfLinearSystemRow6 = m_pyrLinearSystem6.GetLevel(iLevel).Ptr(0, iY);
        float *pfDFlowAboveRowX;
        float *pfDFlowAboveRowY;
        if (iY>0)
        {
            pfDFlowAboveRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY-1);
            pfDFlowAboveRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY-1);
        }
        else
        {
            pfDFlowAboveRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY+1);
            pfDFlowAboveRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY+1);
        }
        float *pfDFlowBelowRowX;
        float *pfDFlowBelowRowY;
        if (iY<iHeight-1)
        {
            pfDFlowBelowRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY+1);
            pfDFlowBelowRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY+1);
        }
        else
        {
            pfDFlowBelowRowX = m_pyrDFlowX.GetLevel(iLevel).Ptr(0, iY-1);
            pfDFlowBelowRowY = m_pyrDFlowY.GetLevel(iLevel).Ptr(0, iY-1);
        }
#ifdef OUTPUT_INNER_ITERATION_ERROR
        float *fpInnerErrorRowX = m_pyrInnerErrorFactorX.GetLevel(iLevel).Ptr(iY);
        float *fpInnerErrorRowY = m_pyrInnerErrorFactorY.GetLevel(iLevel).Ptr(iY);
#endif

        // First Column
        float fDFlowCurrX = *pfDFlowRowX;
        float fDFlowCurrY = *pfDFlowRowY;
        float fDFlowRightX = *pfDFlowRightRowX++;
        float fDFlowRightY = *pfDFlowRightRowY++;
        float fDFlowLeftX = fDFlowRightX;
        float fDFlowLeftY = fDFlowRightY;
        float fDFlowAboveX = *pfDFlowAboveRowX++;
        float fDFlowAboveY = *pfDFlowAboveRowY++;
        float fDFlowBelowX = *pfDFlowBelowRowX++;
        float fDFlowBelowY = *pfDFlowBelowRowY++;
        float fLS1X = *pfLinearSystemRow1++;
        float fLS2X = *pfLinearSystemRow2++;
        float fLS3X = *pfLinearSystemRow3++;
        float fLS1Y = *pfLinearSystemRow4++;
        float fLS2Y = *pfLinearSystemRow5++;
        float fLS3Y = *pfLinearSystemRow6++;
        float fTotalX = fLS1X + fLS2X*fDFlowCurrY + fLS3X*(fDFlowLeftX+fDFlowRightX+fDFlowAboveX+fDFlowBelowX);
        *pfDFlowRowX++ = fDFlowLeftX = fTotalX+m_fOneMinusW*fDFlowCurrX;
        // Note: Must use new value of fDFlowLeftX !!
        float fTotalY = fLS1Y + fLS2Y*fDFlowLeftX + fLS3Y*(fDFlowLeftY+fDFlowRightY+fDFlowAboveY+fDFlowBelowY);
        *pfDFlowRowY++ = fDFlowLeftY = fTotalY+m_fOneMinusW*fDFlowCurrY;

#ifdef OUTPUT_INNER_ITERATION_ERROR
        float fErrorX = (*fpInnerErrorRowX++)*(fDFlowCurrX - fDFlowLeftX);
        float fErrorY = (*fpInnerErrorRowY++)*(fDFlowCurrY - fDFlowLeftY);
        fTotalError += fErrorX*fErrorX + fErrorY*fErrorY;
#endif

        // Columns 2, 3, 4
        while(pfDFlowRowX < pfDFlowRowEndX1)
        {
            fDFlowCurrX = fDFlowRightX;
            fDFlowCurrY = fDFlowRightY;
            fDFlowRightX = *pfDFlowRightRowX++;
            fDFlowRightY = *pfDFlowRightRowY++;
            fDFlowAboveX = *pfDFlowAboveRowX++;
            fDFlowAboveY = *pfDFlowAboveRowY++;
            fDFlowBelowX = *pfDFlowBelowRowX++;
            fDFlowBelowY = *pfDFlowBelowRowY++;
            fLS1X = *pfLinearSystemRow1++;
            fLS2X = *pfLinearSystemRow2++;
            fLS3X = *pfLinearSystemRow3++;
            fLS1Y = *pfLinearSystemRow4++;
            fLS2Y = *pfLinearSystemRow5++;
            fLS3Y = *pfLinearSystemRow6++;
            fTotalX = fLS1X + fLS2X*fDFlowCurrY + fLS3X*(fDFlowLeftX+fDFlowRightX+fDFlowAboveX+fDFlowBelowX);
            *pfDFlowRowX++ = fDFlowLeftX = fTotalX+m_fOneMinusW*fDFlowCurrX;
            // Note: Must use new value of fDFlowLeftX !!
            fTotalY = fLS1Y + fLS2Y*fDFlowLeftX + fLS3Y*(fDFlowLeftY+fDFlowRightY+fDFlowAboveY+fDFlowBelowY);
            *pfDFlowRowY++ = fDFlowLeftY = fTotalY+m_fOneMinusW*fDFlowCurrY;

#ifdef OUTPUT_INNER_ITERATION_ERROR
            fErrorX = (*fpInnerErrorRowX++)*(fDFlowCurrX - fDFlowLeftX);
            fErrorY = (*fpInnerErrorRowY++)*(fDFlowCurrY - fDFlowLeftY);
            fTotalError += fErrorX*fErrorX + fErrorY*fErrorY;
#endif
        }

        // Main Columns
        float *pfTmpX = m_imgSORX.Ptr() + iWidth1;
        float *pfTmpXEnd = m_imgSORX.Ptr() + iWidth2;
        float *pfTmpY = m_imgSORY.Ptr() + iWidth1;
        float *pfLinearSystemRow3Orig = pfLinearSystemRow3;
        float *pfLinearSystemRow5Orig = pfLinearSystemRow5;
        float *pfLinearSystemRow6Orig = pfLinearSystemRow6;
        pfDFlowRightRowX--;
        pfDFlowRightRowY--;
#if (defined(_M_IX86) || defined(_M_AMD64))
        __m128 m128DFlowCurrX = _mm_load_ps(pfDFlowRightRowX);
        __m128 m128DFlowCurrY = _mm_load_ps(pfDFlowRightRowY);
        while(pfTmpX < pfTmpXEnd)
        {
            // Get the Right
            // fDFlowCurrX = fDFlowRightX;
            // fDFlowCurrY = fDFlowRightY;
            // fDFlowRightX = *pfDFlowRightRowX++;
            // fDFlowRightY = *pfDFlowRightRowY++;
            pfDFlowRightRowX+=4;
            pfDFlowRightRowY+=4;
            __m128 m128DFlowRightX = _mm_load_ps(pfDFlowRightRowX);
            __m128 m128DFlowRightY = _mm_load_ps(pfDFlowRightRowY);

            // Get Above and Below
            // fDFlowAboveX = *pfDFlowAboveRowX++;
            // fDFlowAboveY = *pfDFlowAboveRowY++;
            // fDFlowBelowX = *pfDFlowBelowRowX++;
            // fDFlowBelowY = *pfDFlowBelowRowY++;
            __m128 m128DFlowAboveX = _mm_load_ps(pfDFlowAboveRowX);
            __m128 m128DFlowAboveY = _mm_load_ps(pfDFlowAboveRowY);
            pfDFlowAboveRowX+=4;
            pfDFlowAboveRowY+=4;
            __m128 m128DFlowBelowX = _mm_load_ps(pfDFlowBelowRowX);
            __m128 m128DFlowBelowY = _mm_load_ps(pfDFlowBelowRowY);
            pfDFlowBelowRowX+=4;
            pfDFlowBelowRowY+=4;

            // Shuffle for Right
            __m128 m128RXT1 = _mm_move_ss(m128DFlowCurrX, m128DFlowRightX);
            __m128 m128RYT1 = _mm_move_ss(m128DFlowCurrY, m128DFlowRightY);
            __m128 m128RXT2 = _mm_shuffle_ps(m128RXT1, m128RXT1, _MM_SHUFFLE(0, 3, 2, 1));
            __m128 m128RYT2 = _mm_shuffle_ps(m128RYT1, m128RYT1, _MM_SHUFFLE(0, 3, 2, 1));

            // Get the Linear System
            // fLS1X = *pfLinearSystemRow1++;
            // fLS2X = *pfLinearSystemRow2++;
            // fLS3X = *pfLinearSystemRow3++;
            // fLS1Y = *pfLinearSystemRow4++;
            // fLS3Y = *pfLinearSystemRow6++;
            __m128 m128LS1X = _mm_load_ps(pfLinearSystemRow1);
            pfLinearSystemRow1+=4;
            __m128 m128LS2X = _mm_load_ps(pfLinearSystemRow2);
            pfLinearSystemRow2+=4;
            __m128 m128LS3X = _mm_load_ps(pfLinearSystemRow3);
            pfLinearSystemRow3+=4;
            __m128 m128LS1Y = _mm_load_ps(pfLinearSystemRow4);
            pfLinearSystemRow4+=4;
            __m128 m128LS3Y = _mm_load_ps(pfLinearSystemRow6);
            pfLinearSystemRow6+=4;

            // Do the Math and Store
            // *pfTmpX++ = m_fOneMinusW*fDFlowCurrX + fLS1X + fLS2X*fDFlowCurrY + fLS3X*(fDFlowRightX+fDFlowAboveX+fDFlowBelowX);
            // *pfTmpY++ = m_fOneMinusW*fDFlowCurrY + fLS1Y + fLS3Y*(fDFlowRightY+fDFlowAboveY+fDFlowBelowY);
            __m128 m128SumX = _mm_mul_ps(m128OneMinusW, m128DFlowCurrX);
            __m128 m128SumY = _mm_mul_ps(m128OneMinusW, m128DFlowCurrY);
            m128SumX = _mm_add_ps(m128SumX, m128LS1X);
            m128SumY = _mm_add_ps(m128SumY, m128LS1Y);
            __m128 m128TmpX = _mm_add_ps(m128DFlowAboveX, m128DFlowBelowX);
            __m128 m128TmpY = _mm_add_ps(m128DFlowAboveY, m128DFlowBelowY);
            m128TmpX = _mm_add_ps(m128TmpX, m128RXT2);
            m128TmpY = _mm_add_ps(m128TmpY, m128RYT2);
            m128TmpX = _mm_mul_ps(m128TmpX, m128LS3X);
            m128TmpY = _mm_mul_ps(m128TmpY, m128LS3Y);
            m128SumX = _mm_add_ps(m128SumX, m128TmpX);
            m128SumY = _mm_add_ps(m128SumY, m128TmpY);
            m128TmpX = _mm_mul_ps(m128LS2X, m128DFlowCurrY);
            m128SumX = _mm_add_ps(m128SumX, m128TmpX);

            // Store
            _mm_storeu_ps(pfTmpX, m128SumX);
            _mm_storeu_ps(pfTmpY, m128SumY);
            pfTmpX+=4;
            pfTmpY+=4;

            // Update Current
            m128DFlowCurrX = _mm_mul_ps(m128DFlowRightX, m128One);
            m128DFlowCurrY = _mm_mul_ps(m128DFlowRightY, m128One);
        }
#endif
        pfTmpX = m_imgSORX.Ptr() + iWidth1;
        pfTmpY = m_imgSORY.Ptr() + iWidth1;
        pfLinearSystemRow3 = pfLinearSystemRow3Orig;
        pfLinearSystemRow5 = pfLinearSystemRow5Orig;
        pfLinearSystemRow6 = pfLinearSystemRow6Orig;
        while(pfTmpX < pfTmpXEnd)
        {
#ifdef OUTPUT_INNER_ITERATION_ERROR
            fDFlowCurrX = *pfDFlowRowX;
            fDFlowCurrY = *pfDFlowRowY;
#endif

            fLS3X = *pfLinearSystemRow3++;
            fLS2Y = *pfLinearSystemRow5++;
            fLS3Y = *pfLinearSystemRow6++;
            *pfDFlowRowX++ = fDFlowLeftX = (*pfTmpX++)+fLS3X*fDFlowLeftX;
            *pfDFlowRowY++ = fDFlowLeftY = (*pfTmpY++)+fLS2Y*fDFlowLeftX+fLS3Y*fDFlowLeftY;

#ifdef OUTPUT_INNER_ITERATION_ERROR
            fErrorX = (*fpInnerErrorRowX++)*(fDFlowCurrX - fDFlowLeftX);
            fErrorY = (*fpInnerErrorRowY++)*(fDFlowCurrY - fDFlowLeftY);
            fTotalError += fErrorX*fErrorX + fErrorY*fErrorY;
#endif
        }
        pfDFlowRightRowX++;
        pfDFlowRightRowY++;

        // Last Few
        fDFlowRightX = *pfDFlowRowX;
        fDFlowRightY = *pfDFlowRowY;
        while(pfDFlowRowX < pfDFlowRowEndX3)
        {
            fDFlowCurrX = fDFlowRightX;
            fDFlowCurrY = fDFlowRightY;
            fDFlowRightX = *pfDFlowRightRowX++;
            fDFlowRightY = *pfDFlowRightRowY++;
            fDFlowAboveX = *pfDFlowAboveRowX++;
            fDFlowAboveY = *pfDFlowAboveRowY++;
            fDFlowBelowX = *pfDFlowBelowRowX++;
            fDFlowBelowY = *pfDFlowBelowRowY++;
            fLS1X = *pfLinearSystemRow1++;
            fLS2X = *pfLinearSystemRow2++;
            fLS3X = *pfLinearSystemRow3++;
            fLS1Y = *pfLinearSystemRow4++;
            fLS2Y = *pfLinearSystemRow5++;
            fLS3Y = *pfLinearSystemRow6++;
            fTotalX = fLS1X + fLS2X*fDFlowCurrY + fLS3X*(fDFlowLeftX+fDFlowRightX+fDFlowAboveX+fDFlowBelowX);
            *pfDFlowRowX++ = fDFlowLeftX = fTotalX+m_fOneMinusW*fDFlowCurrX;
            // Note: Must use new value of fDFlowLeftX !!
            fTotalY = fLS1Y + fLS2Y*fDFlowLeftX + fLS3Y*(fDFlowLeftY+fDFlowRightY+fDFlowAboveY+fDFlowBelowY);
            *pfDFlowRowY++ = fDFlowLeftY = fTotalY+m_fOneMinusW*fDFlowCurrY;

#ifdef OUTPUT_INNER_ITERATION_ERROR
            fErrorX = (*fpInnerErrorRowX++)*(fDFlowCurrX - fDFlowLeftX);
            fErrorY = (*fpInnerErrorRowY++)*(fDFlowCurrY - fDFlowLeftY);
            fTotalError += fErrorX*fErrorX + fErrorY*fErrorY;
#endif
        }

        // Last Column
        fDFlowCurrX = fDFlowRightX;
        fDFlowCurrY = fDFlowRightY;
        fDFlowRightX = fDFlowLeftX;
        fDFlowRightY = fDFlowLeftY;
        fDFlowAboveX = *pfDFlowAboveRowX++;
        fDFlowAboveY = *pfDFlowAboveRowY++;
        fDFlowBelowX = *pfDFlowBelowRowX++;
        fDFlowBelowY = *pfDFlowBelowRowY++;
        fLS1X = *pfLinearSystemRow1++;
        fLS2X = *pfLinearSystemRow2++;
        fLS3X = *pfLinearSystemRow3++;
        fLS1Y = *pfLinearSystemRow4++;
        fLS2Y = *pfLinearSystemRow5++;
        fLS3Y = *pfLinearSystemRow6++;
        fTotalX = fLS1X + fLS2X*fDFlowCurrY + fLS3X*(fDFlowLeftX+fDFlowRightX+fDFlowAboveX+fDFlowBelowX);
        *pfDFlowRowX++ = fDFlowLeftX = fTotalX+m_fOneMinusW*fDFlowCurrX;
        // Note: Must use new value of fDFlowLeftX !!
        fTotalY = fLS1Y + fLS2Y*fDFlowLeftX + fLS3Y*(fDFlowLeftY+fDFlowRightY+fDFlowAboveY+fDFlowBelowY);
        *pfDFlowRowY++ = fTotalY+m_fOneMinusW*fDFlowCurrY;

#ifdef OUTPUT_INNER_ITERATION_ERROR
        fErrorX = (*fpInnerErrorRowX++)*(fDFlowCurrX - fDFlowLeftX);
        fErrorY = (*fpInnerErrorRowY++)*(fDFlowCurrY - fDFlowLeftY);
        fTotalError += fErrorX*fErrorX + fErrorY*fErrorY;
#endif	
    }

#ifdef OUTPUT_INNER_ITERATION_ERROR			
    fTotalError = sqrt(fTotalError / (iWidth*iHeight));
#endif	

	VT_HR_END()
}

HRESULT CVOF::GetFlowX(CFloatImg* &pimgFlowX, int iLevel)
{
    VT_HR_BEGIN()
    
    // Check Inputs OK
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iLevel>=m_iLevels || iLevel<0) ? E_INVALIDARG : S_OK);

    // Return the Result
    pimgFlowX = &m_pyrFlowX.GetLevel(iLevel);

	VT_HR_END()
}

HRESULT CVOF::GetFlowY(CFloatImg* &pimgFlowY, int iLevel)
{
    VT_HR_BEGIN()
    
    // Check Inputs OK
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iLevel>=m_iLevels || iLevel<0) ? E_INVALIDARG : S_OK);

    // Return the Result
    pimgFlowY = &m_pyrFlowY.GetLevel(iLevel);

	VT_HR_END()
}

HRESULT CVOF::GetCurrentImage(float* &pfBuffer, int &iBufferWidth, int &iBufferHeight, 
                               size_t &uiBufferStrideBytes, int iLevel)
{
    VT_HR_BEGIN()
    
    // Check Inputs OK
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iLevel>=m_iLevels || iLevel<0) ? E_INVALIDARG : S_OK);

    // Return the Data
    pfBuffer = m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Ptr();
    iBufferWidth = m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Width();
    iBufferHeight = m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).Height();
    uiBufferStrideBytes = m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel).StrideBytes();

	VT_HR_END()
}

HRESULT CVOF::GetPreviousImage(float* &pfBuffer, int &iBufferWidth, int &iBufferHeight, 
                                size_t &uiBufferStrideBytes, int iLevel)
{
    VT_HR_BEGIN()
    
    // Check Inputs OK
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iLevel>=m_iLevels || iLevel<0) ? E_INVALIDARG : S_OK);
    VT_HR_EXIT((m_iNumberImages<2) ? E_INVALIDARG : S_OK);

    // Return the Data
    pfBuffer = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Ptr();
    iBufferWidth = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Width();
    iBufferHeight = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).Height();
    uiBufferStrideBytes = m_ppyrInputImages[m_iPreviousImage].GetLevel(iLevel).StrideBytes();

	VT_HR_END()
}

HRESULT CVOF::GetWarpedImage(float* &pfBuffer, int &iBufferWidth, int &iBufferHeight, 
                              size_t &uiBufferStrideBytes, int iLevel)
{
    VT_HR_BEGIN()
    
    // Check Inputs OK
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iLevel>=m_iLevels || iLevel<0) ? E_INVALIDARG : S_OK);
    VT_HR_EXIT((m_iNumberImages<2) ? E_INVALIDARG : S_OK);

    // Return the Data
    pfBuffer = m_pyrWarpedInput.GetLevel(iLevel).Ptr();
    iBufferWidth = m_pyrWarpedInput.GetLevel(iLevel).Width();
    iBufferHeight = m_pyrWarpedInput.GetLevel(iLevel).Height();
    uiBufferStrideBytes = m_pyrWarpedInput.GetLevel(iLevel).StrideBytes();

	VT_HR_END()
}

HRESULT CVOF::GetErrorImage(float* &pfBuffer, int &iBufferWidth, int &iBufferHeight, 
                             size_t &uiBufferStrideBytes, int iLevel)
{
    VT_HR_BEGIN()
    
    // Check Inputs OK
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iLevel>=m_iLevels || iLevel<0) ? E_INVALIDARG : S_OK);
    VT_HR_EXIT((m_iNumberImages<2) ? E_INVALIDARG : S_OK);

    // Return the Data
    pfBuffer = m_pyrErrorImage.GetLevel(iLevel).Ptr();
    iBufferWidth = m_pyrErrorImage.GetLevel(iLevel).Width();
    iBufferHeight = m_pyrErrorImage.GetLevel(iLevel).Height();
    uiBufferStrideBytes = m_pyrErrorImage.GetLevel(iLevel).StrideBytes();

	VT_HR_END()
}
