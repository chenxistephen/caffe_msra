//+---------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation 2010.  All rights reserved.
//
//  Description:
//     Optical Flow
//
//  History:
//      10/03/2011 lhoang
//          Implementation of Optical Flow using celiu's algorithms
//			and adapted from sbaker's code.
//----------------------------------------------------------------------------

#include "stdafx.h"
#include "cvof2.h"

using namespace vt;

//Specify the minimum pixel count that the input image can be downsampled to.
static const int CVOF2_TOP_LEVEL_MIN_PIXEL_COUNT = 100;

//Use pixels value between 0-1, otherwise 0-255.
static bool CVOF2_BNORMALIZEPIXELVALUES = true;

//Outer base number of iterations. Actual iterations vary with pyramid level.
static int CVOF2_IOUTERITERATIONS = 1;

//Inner base number of iterations. Actual iterations vary with pyramid level.
static int CVOF2_IINNERITERATIONS = 1;

//S.O.R base number of iterations. Actual iterations vary with pyramid level.
static int CVOF2_ISORITERATIONS = 2;

//Omega used in S.O.R calculations.
static float CVOF2_OMEGA = 1.8f;

//Alpha used in S.O.R calculations.
static float CVOF2_ALPHA = 1;

//Rho value to be used in 1-Rho-1 filter. Only applies if using smooth pyramid.
//Only matters when smooth pyramid construction is turned on
static float CVOF2_RHO = 8;

//Method to warp image: 1 - Bilinear, 2 - Bicubic.
static int CVOF2_IIMAGEWARPMETHOD = 1;

//Use smooth gaussian gradients Dx, Dy, Dt for S.O.R.
static bool CVOF2_BGAUSSIANGRADIENTS = true;

//Use universal smooth derivatives, only applicable in 3-band mode.
//Only matters when features are turned on (through #DEFINE in header file)
static bool CVOF2_BSMOOTHDERIVATIVES = false;

//Noise model: 1 - Gaussian Mixture, 2 - Laplacian.
static int CVOF2_INOISEMODEL = 2;

//Use SSE computations if available, otherwise disable SSE completely
static bool CVOF2_BUSESSE = true;

CVOF2Params vt::OPTICAL_FLOW_V2_PARAMS(CVOF2_TOP_LEVEL_MIN_PIXEL_COUNT,
    CVOF2_BNORMALIZEPIXELVALUES,
    CVOF2_IOUTERITERATIONS,
    CVOF2_IINNERITERATIONS,
    CVOF2_ISORITERATIONS,
    CVOF2_OMEGA,
    CVOF2_ALPHA,
    CVOF2_RHO,
    CVOF2_IIMAGEWARPMETHOD,
    CVOF2_BGAUSSIANGRADIENTS,
    CVOF2_BSMOOTHDERIVATIVES,
    CVOF2_INOISEMODEL,
    CVOF2_BUSESSE);

void CVOF2::Init()
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
    m_bSSE1 = false;
    m_bSSE2 = false;
    m_bSSSE3 = false;
    m_bSSE4_1 = false;
    if (g_SupportSSE1())
    {
        m_bSSE1 = true;
    }
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
    m_LaplacianParameters = NULL;

    //Input parameters
    m_iTopLevelMinPixelCount = 0;
    m_Rho = 0;
    m_bNormalizePixelValues = false;
    m_iOuterIterations = 0;
    m_iInnerIterations = 0;
    m_iSORIterations = 0;
    m_Omega = 0;
    m_Alpha = 0;
    m_Rho = 0;
    m_iImageWarpMethod = 0;
    m_bGaussianGradients = false;
    m_bSmoothDerivatives = false;
    m_iNoiseModel = 0;
    m_bUseSSE = false;
}

HRESULT CVOF2::Allocate(int iInputWidth, int iInputHeight, CVOF2Params* params, int iSubsample /* = 1 */)
{
    VT_HR_BEGIN()

#ifdef USE_FEATURES
    int featureBands = 5; //Number of channels for feature images and other variants of them.
    int pyramidBands = 3; //Number of channels for pyramid images.
#else
    int featureBands = 1;
    int pyramidBands = 1;
#endif

    // First Deallocate
    Deallocate();

    // Read input file and set up parameters
    this->ReadInput(params);

    // Set the Input Width and Height
    m_iInputWidth = iInputWidth;
    m_iInputHeight = iInputHeight;

    PYRAMID_PROPERTIES pyramidProps;
    pyramidProps.eAutoFilter = ePyramidFilterNone;
    pyramidProps.iSubOctaveStepSize = 1;
    pyramidProps.iSubOctavesPerLevel = 1;

    // Check iSubSample in Plausible Range
    VT_HR_EXIT((iSubsample < 1) ? E_FAIL : S_OK);
    m_iSubsample = iSubsample;
    
    // Check Base Image Size is at least 2x2
    m_iBaseWidth = int(ceil(float(m_iInputWidth) / float(m_iSubsample)));
    m_iBaseHeight = int(ceil(float(m_iInputHeight) / float(m_iSubsample)));
    VT_HR_EXIT((m_iBaseWidth < 2 || m_iBaseHeight < 2) ? E_FAIL : S_OK);

    // Compute Number of Levels Reasonable
    m_iLevels = 1;
    int iBaseWidth = m_iBaseWidth;
    int iBaseHeight = m_iBaseHeight;
    while (iBaseWidth * iBaseHeight > 4 * m_iTopLevelMinPixelCount)
    {
        m_iLevels++;
        iBaseWidth = (int)(iBaseWidth / 2);
        iBaseHeight = (int)(iBaseHeight / 2);
    }

    // Allocate the Input Image
    VT_HR_EXIT(m_TempImage.Create(m_iBaseWidth, m_iBaseHeight, featureBands));

    VT_HR_EXIT(m_SmoothFlowTempImage1.Create(m_iBaseWidth, m_iBaseHeight, featureBands));
    VT_HR_EXIT(m_SmoothFlowTempImage2.Create(m_iBaseWidth, m_iBaseHeight, featureBands));

    // Set the Base Levels, Force Allocation, and Add Kernel for the Input Image Pyramids
    //Pyramid for the first image
    VT_HR_EXIT(m_ppyrInputImages[0].Create(m_iBaseWidth, m_iBaseHeight, pyramidBands, &pyramidProps));

    //Pyramid for the second image
    VT_HR_EXIT(m_ppyrInputImages[1].Create(m_iBaseWidth, m_iBaseHeight, pyramidBands, &pyramidProps));

#ifdef USE_FEATURES
    //Feature images
    VT_HR_EXIT(m_ppyrFeatureImages[0].Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps));

    VT_HR_EXIT(m_ppyrFeatureImages[1].Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps));
#endif

    VT_HR_EXIT(m_pyrDerivsX.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps)); //x-derivative
    VT_HR_EXIT(m_pyrDerivsY.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps)); //y-derivative
    VT_HR_EXIT(m_pyrWarpedInput.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps)); //warped input
    VT_HR_EXIT(m_pyrErrorImage.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps)); //error image

    VT_HR_EXIT(m_imgiSrcX.Create(4, 1, 1));
    VT_HR_EXIT(m_imgiSrcY.Create(4, 1, 1));

    // Flow Pyramids
    VT_HR_EXIT(m_pyrFlowX.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(m_pyrFlowY.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(m_pyrDFlowX.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(m_pyrDFlowY.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(m_pyrFlowLaplacianX.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(m_pyrFlowLaplacianY.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));

    //LH:Todo:Rename these variables
    //Pyramids used for SOR computation ---------------------------
    VT_HR_EXIT(uu.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(vv.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(ux.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(uy.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(vx.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(vy.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));

    VT_HR_EXIT(Phi_1st.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(Psi_1st.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps));
    
    VT_HR_EXIT(m_imgSORX.Create(m_iBaseWidth, 1, 1));
    VT_HR_EXIT(m_imgSORY.Create(m_iBaseWidth, 1, 1));
    VT_HR_EXIT(m_imgSORPhi.Create(m_iBaseWidth, 1, 1));

#ifdef USE_FEATURES
    VT_HR_EXIT(imdxy.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(imdx2.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(imdy2.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(imdtdx.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
    VT_HR_EXIT(imdtdy.Create(m_iBaseWidth, m_iBaseHeight, 1, &pyramidProps));
#endif
    VT_HR_EXIT(ImDxy.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps));
    VT_HR_EXIT(ImDx2.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps));
    VT_HR_EXIT(ImDy2.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps));
    VT_HR_EXIT(ImDtDx.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps));
    VT_HR_EXIT(ImDtDy.Create(m_iBaseWidth, m_iBaseHeight, featureBands, &pyramidProps));

    // Finally set allocated flag
    m_bAllocated = true;

    VT_HR_EXIT_LABEL()

    if (hr != S_OK)
    {
        Deallocate();
    }
    return hr;
}

void CVOF2::Deallocate()
{
    m_TempImage.Deallocate();
    m_imgTmpImage.Deallocate();
    m_imgTmpImage2.Deallocate();
    m_SmoothFlowTempImage1.Deallocate();
    m_SmoothFlowTempImage2.Deallocate();
    m_ppyrInputImages[0].Deallocate();
    m_ppyrInputImages[1].Deallocate();
#ifdef USE_FEATURES
    m_ppyrFeatureImages[0].Deallocate();
    m_ppyrFeatureImages[1].Deallocate();
#endif
    
    m_pyrDerivsX.Deallocate();
    m_pyrDerivsY.Deallocate();
    m_pyrWarpedInput.Deallocate();
    m_pyrErrorImage.Deallocate();
    m_pyrFlowX.Deallocate();
    m_pyrFlowY.Deallocate();
    m_pyrDFlowX.Deallocate();
    m_pyrDFlowY.Deallocate();
    m_pyrFlowLaplacianX.Deallocate();
    m_pyrFlowLaplacianY.Deallocate();

    //LH:Todo:Rename variables
    //SOR pyramids and images
    uu.Deallocate();
    vv.Deallocate();
    ux.Deallocate();
    uy.Deallocate();
    vx.Deallocate();
    vy.Deallocate();
    Phi_1st.Deallocate();
    Psi_1st.Deallocate();
    m_imgSORX.Deallocate();
    m_imgSORY.Deallocate();
    m_imgSORPhi.Deallocate();
#ifdef USE_FEATURES
    imdxy.Deallocate();
    imdx2.Deallocate();
    imdy2.Deallocate();
    imdtdx.Deallocate();
    imdtdy.Deallocate();
#endif
    ImDxy.Deallocate();
    ImDx2.Deallocate();
    ImDy2.Deallocate();
    ImDtDx.Deallocate();
    ImDtDy.Deallocate();

    Init();
}

HRESULT CVOF2::ReadInput(CVOF2Params* params)
{
    HRESULT hr = S_OK;

    m_iTopLevelMinPixelCount = params->m_iTopLevelMinPixelCount;
    m_bNormalizePixelValues = params->m_bNormalizePixelValues;
    m_iOuterIterations = params->m_iOuterIterations;
    m_iInnerIterations = params->m_iInnerIterations;
    m_iSORIterations = params->m_iSORIterations;
    m_Omega = params->m_Omega;
    m_Alpha = params->m_Alpha;
    m_Rho = params->m_Rho;
    m_iImageWarpMethod = params->m_iImageWarpMethod;
    m_bGaussianGradients = params->m_bGaussianGradients;
    m_bSmoothDerivatives = params->m_bSmoothDerivatives;
    m_iNoiseModel = params->m_iNoiseModel;
    m_bUseSSE = params->m_bUseSSE;

    if (hr != S_OK)
    {
        this->Deallocate();
    }
    return hr;
}

HRESULT CVOF2::AddBGRAImage(CRGBAImg &imgInput)
{
    VT_HR_BEGIN()

    int iBufferWidth = imgInput.Width();
    int iBufferHeight = imgInput.Height();

    // Check Allocated and Input Buffer is the Right Size
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iBufferWidth != m_iInputWidth || iBufferHeight != m_iInputHeight) ? E_FAIL : S_OK);

    // Copy into Base Level, Converting RGB 2 Luminance
    //This is simply to index into the right pyramid list
    //The list contains two pyramids with index 0 to be used for the first image and index 1 for the second image.
    m_iCurrentImage = m_iPreviousImage;
    m_iPreviousImage = 1 - m_iCurrentImage;
    m_iNumberImages++;

    CFloatImg* currentBaseLevelImage = &(m_ppyrInputImages[m_iCurrentImage].GetLevel(0));

#ifdef USE_FEATURES
    size_t uiBufferStrideBytes = imgInput.StrideBytes();
    int nInputBands = currentBaseLevelImage->Bands();
    unsigned char *pucDestinationData = (unsigned char*)currentBaseLevelImage->Ptr(); //Get the original size image
    unsigned int uiDestinationStride = currentBaseLevelImage->StrideBytes(); //Get the length of each row of pixels of the original size image

    Byte* pImageStart = imgInput.BytePtr();
    for (int iY = 0; iY < m_iInputHeight / m_iSubsample; iY += 1)
    {
        //Get a pointer to the current working source row based on level iY
        const unsigned char *pucSourceRow = pImageStart + m_iSubsample * iY * uiBufferStrideBytes;

        //Get a pointer to the current working destination row based on level iY
        float *pfDestinationRow = ((float*) (pucDestinationData + iY * uiDestinationStride));
        float *pfDestinationRowEnd = pfDestinationRow + m_iBaseWidth * nInputBands;

        while (pfDestinationRow < pfDestinationRowEnd)
        {
            if (m_bNormalizePixelValues)
            {
                *pfDestinationRow++ = (float)*pucSourceRow++ / 255; // Blue
                *pfDestinationRow++ = (float)*pucSourceRow++ / 255; // Green
                *pfDestinationRow++ = (float)*pucSourceRow++ / 255; // Red
            }
            else
            {
                *pfDestinationRow++ = (float)*pucSourceRow++; // Blue
                *pfDestinationRow++ = (float)*pucSourceRow++; // Green
                *pfDestinationRow++ = (float)*pucSourceRow++; // Red
            }
            pucSourceRow += 1; //Skip Alpha value
        }
    }
#else
    /*VT*/::VtConvertImage(*currentBaseLevelImage, imgInput);
#endif

    // Create the Pyramid
    for (int iLevel = 0; iLevel < m_iLevels - 1; iLevel++)
    {
        CFloatImg* currentImage = &(m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel));
        CFloatImg* nextImage = &(m_ppyrInputImages[m_iCurrentImage].GetLevel(iLevel + 1));
#ifdef USE_FEATURES
        if (!m_imgTmpImage.IsValid())
        {
            VT_HR_EXIT(m_imgTmpImage.Create(m_iBaseWidth, m_iBaseHeight, 3));
            VT_HR_EXIT(m_imgTmpImage2.Create(m_iBaseWidth, m_iBaseHeight, 3));
        }
        double filter[3] = { 1.0 / (m_Rho + 2), m_Rho / (m_Rho + 2), 1.0 / (m_Rho + 2)};
        this->ApplyFilter_Horizontal(*currentImage, m_imgTmpImage, filter, 1, currentImage->Width(), currentImage->Height());
        this->ApplyFilter_Vertical(m_imgTmpImage, m_imgTmpImage2, filter, 1, currentImage->Width(), currentImage->Height());
        CRect dstRect(0, 0, nextImage->Width(), nextImage->Height());
        /*VT*/VT_HR_EXIT(::VtResizeImage(*nextImage, dstRect, m_imgTmpImage2, 2.0f, 1.0f, 2.0f, 1.0f, eSamplerKernelBilinear));
#else
        CRect dstRect(0, 0, nextImage->Width(), nextImage->Height());
        /*VT*/VT_HR_EXIT(::VtResizeImage(*nextImage, dstRect, *currentImage, 2.0f, 1.0f, 2.0f, 1.0f, eSamplerKernelBilinear));
#endif
    }

    VT_HR_EXIT_LABEL()

    if (hr != S_OK)
    {
        Deallocate();
    }
    return hr;
}

HRESULT CVOF2::ComputeFlow(int iLowestLevel)
{
    VT_HR_BEGIN()

    //Initialize noise with according to specified noise model
    switch (m_iNoiseModel)
    {
        //LH:Todo:Implement gaussian noise model
        case GaussianMixture:
            break;
        case Laplacian:
        {
            int size = m_ppyrInputImages[0].GetLevel(0).Bands() + 2;
            m_LaplacianParameters = new double[size];
            for (int i = 0; i < size; i++)
            {
                m_LaplacianParameters[i] = 0.02;
            }
            break;
        }
    }
    // Start at top level and work down to iLowestLevel
    for (int iLevel = m_iLevels - 1; iLevel >= iLowestLevel; iLevel--)
    {
        CFloatImg* pWarpedImage = &(m_pyrWarpedInput.GetLevel(iLevel));
#ifdef USE_FEATURES
        CFloatImg* pFeatureImage2 = &(m_ppyrFeatureImages[1].GetLevel(iLevel));
        //Create the feature images for the two image
        //LH:Todo: Convert this method to take in input images for easier testing
        this->ImageToFeature(0 /*First image*/, iLevel); 
        this->ImageToFeature(1 /*Second image*/, iLevel);
#endif

        if (iLevel == m_iLevels - 1) //Top level
        {
            //Need to Zero the Flow in the Top Level -> Might have been overwritten
            VT_HR_EXIT(m_pyrFlowX.GetLevel(iLevel).Clear());
            VT_HR_EXIT(m_pyrFlowY.GetLevel(iLevel).Clear());

#ifdef USE_FEATURES
            //Copy the second feature image into the warped input image
            /*VT*/pFeatureImage2->CopyTo(*pWarpedImage);
#else
            /*VT*/m_ppyrInputImages[1].GetLevel(iLevel).CopyTo(*pWarpedImage);
#endif
        }
        else
        {
            VT_HR_EXIT(this->ExpandFlow(iLevel));
            //LH:Todo: Convert this method to take in input images for easier testing
            this->WarpImage(iLevel);
        }

        double alpha = m_Alpha;
        //LH:Todo: Convert this method to take in input images for easier testing
        this->SmoothFlowSOR(iLevel, alpha, 
            m_iOuterIterations + iLevel * 2,
            m_iInnerIterations,
            m_iSORIterations + iLevel * 3);
    }

    VT_HR_EXIT_LABEL()
    return hr;
}

HRESULT CVOF2::ComputeDerivativeX(CFloatImg& sourceImage, CFloatImg& derivativeXImage, bool useSmoothDerivative, int derivativeWidth /* = 0 */, int derivativeHeight /* = 0 */)
{
    HRESULT hr = S_OK;

#ifdef USE_FEATURES
    if (useSmoothDerivative)
    {
        double derivativeFilter[5] = { 1.0/12, -8.0/12, 0, 8.0/12, -1.0/12 };
        this->ApplyFilter_Horizontal(sourceImage, derivativeXImage, derivativeFilter, 2, derivativeWidth, derivativeHeight);
    }
    else
    {
        double derivativeFilter[3] = { -0.5, 0, 0.5};
        this->ApplyFilter_Horizontal(sourceImage, derivativeXImage, derivativeFilter, 1, derivativeWidth, derivativeHeight);
    }
#else
    UNREFERENCED_PARAMETER(useSmoothDerivative);
    int iWidth = (derivativeWidth == 0) ? sourceImage.Width() : derivativeWidth;
    int iHeight = (derivativeHeight == 0) ? sourceImage.Height() : derivativeHeight;
#if (defined(_M_IX86) || defined(_M_AMD64))
    if (m_bUseSSE && m_bSSE2)
    {
        __m128 m128Left;
        __m128 m128Right;
        float fRight1, fRight2, fRight3, fRight4;
        float fLeft1, fLeft2, fLeft3, fLeft4;
        
        __m128 half = ::_mm_set_ps1(0.5f);
        
        int remainder = iWidth % 4;
        //The last pixel is computed a little differently so when width is divisible by 4, 
        //we have to process the last 4 pixels separately. Otherwise SSE will bundle and 
        //do them all the same way which is incorrect.
        remainder = (remainder == 0) ? 4 : remainder;
        
        //The first 4 pixels and the last "remainder" pixels will be computed manually
        int sseLength = (iWidth - remainder - 4) / 4;
        
        for (int iY = 0; iY < iHeight; iY++)
        {
            //Perform Non-SSE calculation for the first 4 pixels
            float *fpSrcRow = sourceImage.Ptr(iY);
            float *fpDstRowX = derivativeXImage.Ptr(iY);
            float *fLeft = fpSrcRow;
            float *fRight = fpSrcRow + 1;
            *fpDstRowX++ = *fRight++ - *fLeft;
            *fpDstRowX++ = 0.5f * (*fRight++ - *fLeft++);
            *fpDstRowX++ = 0.5f * (*fRight++ - *fLeft++);
            *fpDstRowX++ = 0.5f * (*fRight++ - *fLeft++);

            //Run SSE calculations for the middle pixels
            __m128* m128DestinationRow = (__m128*)fpDstRowX;
            for (int i = 0; i < sseLength; i++)
            {
                fLeft1 = *fLeft++;
                fLeft2 = *fLeft++;
                fLeft3 = *fLeft++;
                fLeft4 = *fLeft++;
                m128Left = ::_mm_set_ps(fLeft4, fLeft3, fLeft2, fLeft1);

                fRight1 = *fRight++;
                fRight2 = *fRight++;
                fRight3 = *fRight++;
                fRight4 = *fRight++;
                m128Right = ::_mm_set_ps(fRight4, fRight3, fRight2, fRight1);

                *m128DestinationRow++ = ::_mm_mul_ps(half, ::_mm_sub_ps(m128Right, m128Left));
            }
            
            ////Perform Non-SSE calculation for the remainder pixels
            fpDstRowX = (float*)m128DestinationRow;
            for (int i = 0; i < remainder - 1; i++)
            {
                *fpDstRowX++ = 0.5f * (*fRight++ - *fLeft++);
            }
            *fpDstRowX = *(fRight - 1) - *(fLeft);
        }
    }
    else
#endif
    {
        for (int iY = 0; iY < iHeight; iY++)
        {
            float *fpSrcRow = sourceImage.Ptr(iY);
            float *fpDstRowX = derivativeXImage.Ptr(iY);
            float fLeft = *fpSrcRow++;
            float fCurrent = fLeft; 
            float fRight = *fpSrcRow++;
            *fpDstRowX++ = fRight - fLeft;
            float *fpSrcRowEnd = sourceImage.Ptr(iY) + iWidth;
            while (fpSrcRow < fpSrcRowEnd)
            {
                fLeft = fCurrent;
                fCurrent = fRight;
                fRight = *fpSrcRow++;
                *fpDstRowX++ = 0.5f * (fRight - fLeft);
            }
            *fpDstRowX = fRight - fCurrent;
        }
    }
#endif

    return hr;
}

HRESULT CVOF2::ComputeDerivativeY(CFloatImg& sourceImage, CFloatImg& derivativeYImage, bool useSmoothDerivative, int derivativeWidth /* = 0 */, int derivativeHeight /* = 0 */)
{
    HRESULT hr = S_OK;

#ifdef USE_FEATURES
    if (useSmoothDerivative)
    {
        double derivativeFilter[5] = { 1.0/12, -8.0/12, 0, 8.0/12, -1.0/12 };
        this->ApplyFilter_Vertical(sourceImage, derivativeYImage, derivativeFilter, 2, derivativeWidth, derivativeHeight);
    }
    else
    {
        double derivativeFilter[3] = { -.5, 0, .5};
        this->ApplyFilter_Vertical(sourceImage, derivativeYImage, derivativeFilter, 1, derivativeWidth, derivativeHeight);
    }
#else
    UNREFERENCED_PARAMETER(useSmoothDerivative);
    int iWidth = (derivativeWidth == 0) ? sourceImage.Width() : derivativeWidth;
    int iHeight = (derivativeHeight == 0) ? sourceImage.Height() : derivativeHeight;
    
#if (defined(_M_IX86) || defined(_M_AMD64))
    if (m_bUseSSE && m_bSSE2)
    {
        int remainder = iWidth % 4;
        int sseLength = iWidth / 4;
        __m128 m128OneHalf = ::_mm_set_ps1(0.5);

        //Process first row
        __m128* m128SrcRowAbove = (__m128*)sourceImage.Ptr();
        __m128* m128SrcRowBelow = (__m128*)sourceImage.Ptr(1);
        __m128* m128DstRowY = (__m128*)derivativeYImage.Ptr();
        for (int i = 0; i < sseLength; i++)
        {
            *m128DstRowY++ = ::_mm_sub_ps(*m128SrcRowBelow++, *m128SrcRowAbove++);
        }
        float* fpDstRowY = (float*)m128DstRowY;
        float* fpSrcRowBelow = (float*)m128SrcRowBelow;
        float* fpSrcRowAbove = (float*)m128SrcRowAbove;
        for (int i = 0; i < remainder; i++)
        {
            *fpDstRowY++ = (*fpSrcRowBelow++ - *fpSrcRowAbove++);
        }

        //Process middle rows
        for(int iY = 1; iY < iHeight - 1; iY++)
        {
            m128SrcRowAbove = (__m128*)sourceImage.Ptr(iY - 1);
            m128SrcRowBelow = (__m128*)sourceImage.Ptr(iY + 1);
            m128DstRowY = (__m128*)derivativeYImage.Ptr(iY);
            for (int i = 0; i < sseLength; i++)
            {
                *m128DstRowY++ = ::_mm_mul_ps(::_mm_sub_ps(*m128SrcRowBelow++, *m128SrcRowAbove++), m128OneHalf);
            }
            fpDstRowY = (float*)m128DstRowY;
            fpSrcRowBelow = (float*)m128SrcRowBelow;
            fpSrcRowAbove = (float*)m128SrcRowAbove;
            for (int i = 0; i < remainder; i++)
            {
                *fpDstRowY++ = 0.5f * (*fpSrcRowBelow++ - *fpSrcRowAbove++);
            }
        }

        //Process last row
        m128SrcRowAbove = (__m128*)sourceImage.Ptr(iHeight - 2);
        m128SrcRowBelow = (__m128*)sourceImage.Ptr(iHeight - 1);
        m128DstRowY = (__m128*)derivativeYImage.Ptr(iHeight - 1);
        for (int i = 0; i < sseLength; i++)
        {
            *m128DstRowY++ = ::_mm_sub_ps(*m128SrcRowBelow++, *m128SrcRowAbove++);
        }
        fpDstRowY = (float*)m128DstRowY;
        fpSrcRowBelow = (float*)m128SrcRowBelow;
        fpSrcRowAbove = (float*)m128SrcRowAbove;
        for (int i = 0; i < remainder; i++)
        {
            *fpDstRowY++ = *fpSrcRowBelow++ - *fpSrcRowAbove++;
        }
    }
    else
#endif
    {
        float *fpSrcRowAbove = sourceImage.Ptr();
        float *fpSrcRowBelow = sourceImage.Ptr(1);
        float *fpDstRowY = derivativeYImage.Ptr();
        float *fpDstRowEndY = fpDstRowY + iWidth;
        while (fpDstRowY < fpDstRowEndY)
        {
            *fpDstRowY++ = (*fpSrcRowBelow++ - *fpSrcRowAbove++);
        }
    
        for(int iY = 1; iY < iHeight - 1; iY++)
        {
            fpSrcRowAbove = sourceImage.Ptr(iY - 1);
            fpSrcRowBelow = sourceImage.Ptr(iY + 1);
            fpDstRowY = derivativeYImage.Ptr(iY);
            fpDstRowEndY = fpDstRowY + iWidth;
            while (fpDstRowY < fpDstRowEndY)
            {
                *fpDstRowY++ = 0.5f * (*fpSrcRowBelow++ - *fpSrcRowAbove++);
            }
        }

        fpSrcRowAbove = sourceImage.Ptr(iHeight - 2);
        fpSrcRowBelow = sourceImage.Ptr(iHeight - 1);
        fpDstRowY = derivativeYImage.Ptr(iHeight - 1);
        fpDstRowEndY = fpDstRowY + iWidth;
        while (fpDstRowY < fpDstRowEndY)
        {
            *fpDstRowY++ = (*fpSrcRowBelow++ - *fpSrcRowAbove++);
        }
    }
#endif

    return hr;
}

HRESULT CVOF2::ComputePhi(CFloatImg& derivativeXFlowX, CFloatImg& derivativeYFlowX, CFloatImg& derivativeXFlowY, CFloatImg& derivativeYFlowY, CFloatImg& phi, double& epsilon)
{
    HRESULT hr = S_OK;

    int width = phi.Width();
    int height = phi.Height();
    
#if (defined(_M_IX86) || defined(_M_AMD64))
    if (m_bUseSSE && m_bSSE2)
    {
        int remainder = width % 4;
        int sseLength = width / 4;
        __m128 m128temp; //temporary variable to help with computation
        __m128 eps = ::_mm_set_ps1(float(epsilon));
        __m128 scalar = ::_mm_set_ps1(0.5);
        //Compute the weight of phi
        for (int iY = 0; iY < height; iY++)
        {
            __m128* m128uxRow = (__m128*)derivativeXFlowX.Ptr(iY);
            __m128* m128uyRow = (__m128*)derivativeYFlowX.Ptr(iY);
            __m128* m128vxRow = (__m128*)derivativeXFlowY.Ptr(iY);
            __m128* m128vyRow = (__m128*)derivativeYFlowY.Ptr(iY);

            __m128* m128phiDataRow = (__m128*)phi.Ptr(iY);
            
            for (int i = 0; i < sseLength; i++)
            {
                m128temp = ::_mm_add_ps(::_mm_add_ps(::_mm_mul_ps(*m128uxRow, *m128uxRow), ::_mm_mul_ps(*m128uyRow, *m128uyRow)), 
                                    ::_mm_add_ps(::_mm_mul_ps(*m128vxRow, *m128vxRow), ::_mm_mul_ps(*m128vyRow, *m128vyRow)));
                m128temp = ::_mm_sqrt_ps(::_mm_add_ps(m128temp, eps));
                *m128phiDataRow++ = ::_mm_div_ps(scalar, m128temp);
                m128uxRow++;
                m128uyRow++;
                m128vxRow++;
                m128vyRow++;
            }

            if (remainder > 0)
            {
                float* pfuxRow = (float*)m128uxRow;
                float* pfuyRow = (float*)m128uyRow;
                float* pfvxRow = (float*)m128vxRow;
                float* pfvyRow = (float*)m128vyRow;
                float* phiDataRow = (float*)m128phiDataRow;
                for (int i = 0; i < remainder; i++)
                {
                    float temp = (*pfuxRow) * (*pfuxRow++) + (*pfuyRow) * (*pfuyRow++) +
                                 (*pfvxRow) * (*pfvxRow++) + (*pfvyRow) * (*pfvyRow++);
                    
                    *phiDataRow++ = 0.5F / ::sqrt(temp + (float)epsilon);
                }
            }
        }
    }
    else
#endif
    {
        float temp; //temporary variable to help with computation
        //Compute the weight of phi
        for (int iY = 0; iY < height; iY++)
        {
            const float* pfuxRow = derivativeXFlowX.Ptr(iY);
            const float* pfuyRow = derivativeYFlowX.Ptr(iY);
            const float* pfvxRow = derivativeXFlowY.Ptr(iY);
            const float* pfvyRow = derivativeYFlowY.Ptr(iY);
            float* phiDataRow = phi.Ptr(iY);
            float* phiDataRowEnd = phiDataRow + width;
            while (phiDataRow < phiDataRowEnd)
            {
                temp = (*pfuxRow) * (*pfuxRow++) + (*pfuyRow) * (*pfuyRow++) +
                        (*pfvxRow) * (*pfvxRow++) + (*pfvyRow) * (*pfvyRow++);
                    
                *phiDataRow++ = 0.5F / ::sqrt(temp + (float)epsilon);
            }
        }
    }
    return hr;
}

HRESULT CVOF2::ComputePsi(CFloatImg& derivativeX, CFloatImg& derivativeY, CFloatImg& errorImage, CFloatImg& dflowX, CFloatImg& dflowY, CFloatImg& psi, double*& parameters, double& epsilon)
{
    HRESULT hr = S_OK;

    int width = psi.Width();
    int height = psi.Height();
    float temp;
    //Compute the nonlinear term of psi
#ifdef USE_FEATURES
    int nChannels = psi.Bands();
    for (int iY = 0; iY < height; iY++)
    {
        const float* imdxData = derivativeX.Ptr(iY);
        const float* imdyData = derivativeY.Ptr(iY);
        const float* imdtData = errorImage.Ptr(iY);
        const float* duData   = dflowX.Ptr(iY);
        const float* dvData   = dflowY.Ptr(iY);
        float* psiData  = psi.Ptr(iY);

        for (int iX = 0; iX < width; iX++)
        {
            for (int k = 0; k < nChannels; k++)
            {
                temp = (*imdtData++) + (*imdxData++) * (*duData) +
                                       (*imdyData++) * (*dvData);
                temp *= temp;
                switch (m_iNoiseModel)
                {
                    //LH:Todo:Implement gaussian noise model
                    case GaussianMixture:
                        break;
                    case Laplacian:
                        if (parameters[k] < 1E-20)
                            continue;
                        *psiData++ = 1.0f / (2 * ::sqrt(temp + (float)epsilon) * (float)parameters[k]);
                        break;
                }
            }
            duData++;
            dvData++;
        }
    }
#else
#if (defined(_M_IX86) || defined(_M_AMD64))
    if (m_bUseSSE && m_bSSE2)
    {
        switch (m_iNoiseModel)
        {
            //LH:Todo:Implement gaussian noise model
            case GaussianMixture:
                break;
            case Laplacian:
            {
                float noiseValue = float(parameters[0]);
                if (noiseValue >= 1E-20)
                {
                    int remainder = width % 4;
                    int sseLength = width / 4;
                    __m128 m128temp; //temporary variable to help with computation
                    __m128 one = ::_mm_set_ps1(1.0f);
                    __m128 eps = ::_mm_set_ps1((float)epsilon);
                    __m128 noise = ::_mm_set_ps1(2 * noiseValue);

                    for (int iY = 0; iY < height; iY++)
                    {
                        __m128* imdxData = (__m128*)derivativeX.Ptr(iY);
                        __m128* imdyData = (__m128*)derivativeY.Ptr(iY);
                        __m128* imdtData = (__m128*)errorImage.Ptr(iY);
                        __m128* duData   = (__m128*)dflowX.Ptr(iY);
                        __m128* dvData   = (__m128*)dflowY.Ptr(iY);
                        __m128* psiData  = (__m128*)psi.Ptr(iY);

                        for (int i = 0; i < sseLength; i++)
                        {
                            m128temp = ::_mm_add_ps(::_mm_add_ps(::_mm_mul_ps(*imdxData++, *duData++), 
                                                             ::_mm_mul_ps(*imdyData++, *dvData++)),
                                                *imdtData++);
                            m128temp = ::_mm_mul_ps(m128temp, m128temp);
                            m128temp = ::_mm_sqrt_ps(::_mm_add_ps(m128temp, eps));
                            *psiData++ = ::_mm_div_ps(one, ::_mm_mul_ps(m128temp, noise));	
                        }
                        if (remainder > 0)
                        {
                            float* pfimdxData = (float*)imdxData;
                            float* pfimdyData = (float*)imdyData;
                            float* pfimdtData = (float*)imdtData;
                            float* pfduData   = (float*)duData;
                            float* pfdvData   = (float*)dvData;
                            float* pfpsiData  = (float*)psiData;
                            for (int i = 0; i < remainder; i++)
                            {
                                float ftemp = (*pfimdtData++) + (*pfimdxData++) * (*pfduData++) +
                                                                (*pfimdyData++) * (*pfdvData++);
                                ftemp *= ftemp;
                                *pfpsiData++ = 1 / (2 * ::sqrt(ftemp + (float)epsilon) * (float)parameters[0]);
                            }
                        }
                    }
                }
                break;
            }
        }
    }
    else
#endif
    {
        switch (m_iNoiseModel)
        {
            //LH:Todo:Implement gaussian noise model
            case GaussianMixture:
                break;
            case Laplacian:
            {
                if (parameters[0] >= 1E-20)
                {
                    for (int iY = 0; iY < height; iY++)
                    {
                        const float* imdxData = derivativeX.Ptr(iY);
                        const float* imdyData = derivativeY.Ptr(iY);
                        const float* imdtData = errorImage.Ptr(iY);
                        const float* duData   = dflowX.Ptr(iY);
                        const float* dvData   = dflowY.Ptr(iY);
                        float* psiData  = psi.Ptr(iY);
                        for (int iX = 0; iX < width; iX++)
                        {
                            temp = (*imdtData++) + (*imdxData++) * (*duData++) +
                                                   (*imdyData++) * (*dvData++);
                            temp *= temp;
                            *psiData++ = 1 / (2 * ::sqrt(temp + (float)epsilon) * (float)parameters[0]);
                        }
                    }
                }
                break;
            }
        }
    }
#endif

    return hr;
}

void CVOF2::WarpImage(int iLevel /* = 0 */)
{
    //The SSE implementation of warp image is currently only applicable to 1-band images.
#ifndef USE_FEATURES
    if (m_bUseSSE && m_bSSE2 && m_bSSSE3)
    {
        CFloatImg* pFeatureImage1 = &(m_ppyrInputImages[0].GetLevel(iLevel));
        CFloatImg* pFeatureImage2 = &(m_ppyrInputImages[1].GetLevel(iLevel));
        CFloatImg* pFlowXImage	  = &(m_pyrFlowX.GetLevel(iLevel));
        CFloatImg* pFlowYImage	  = &(m_pyrFlowY.GetLevel(iLevel));
        CFloatImg* pWarpedImage   = &(m_pyrWarpedInput.GetLevel(iLevel));
        vt::OpticalFlow_WarpImageSSE(*pFeatureImage2, *pFlowXImage, *pFlowYImage, *pWarpedImage, m_SmoothFlowTempImage1, pWarpedImage->Rect(), *pFeatureImage1);
    }
    else
#endif
    {
        switch (m_iImageWarpMethod)
        {
            case Bilinear:
            {
                this->WarpImage_Bilinear(iLevel);
                break;
            }
            default:
            {
                break;
            }
        }
    }
}

void CVOF2::WarpImage_Bilinear(int iLevel /* = 0 */)
{
#ifdef USE_FEATURES
    CFloatImg* pFeatureImage1 = &(m_ppyrFeatureImages[0].GetLevel(iLevel));
    CFloatImg* pFeatureImage2 = &(m_ppyrFeatureImages[1].GetLevel(iLevel));
#else
    CFloatImg* pFeatureImage1 = &(m_ppyrInputImages[0].GetLevel(iLevel));
    CFloatImg* pFeatureImage2 = &(m_ppyrInputImages[1].GetLevel(iLevel));
#endif
    CFloatImg* pWarpedImage   = &(m_pyrWarpedInput.GetLevel(iLevel));

    int width = pFeatureImage2->Width();
    int height = pFeatureImage2->Height();
    int channels = pFeatureImage2->Bands();
    size_t strideBytes = pFeatureImage2->StrideBytes();

    for (int iY = 0; iY < height; iY++)
    {
        //Get the feature images
        float* pfFeatureImage1Row = pFeatureImage1->Ptr(iY);

        //Get the warped image
        float* pfWarpedImageRow = pWarpedImage->Ptr(iY);

        //Get the flow images for x-component and y-component
        float* pfFlowX = m_pyrFlowX.GetLevel(iLevel).Ptr(iY);
        float* pfFlowY = m_pyrFlowY.GetLevel(iLevel).Ptr(iY);

        for (int iX = 0; iX < width; iX++)
        {
            double x = (*pfFlowX++) + iX; //Compute x-coordinate of warped pixel
            double y = (*pfFlowY++) + iY; //Compute y-coordinate of warped pixel
            
            //Offset into the correct pixel of the warped image
            //since each flow pixel corresponds to "channels" number of pixels in warped image.
            int offset = iX * channels;

            //If calculated warped pixel is out of image boundary then simply
            //reset it to be same as the corresponding pixel in the first image
            if (x < 0 || x > width - 1 || y < 0 || y > height - 1)
            {
                float* pfWarpedImagePixel = pfWarpedImageRow + offset;
                float* pfFeatureImage1Pixel = pfFeatureImage1Row + offset;
                for (int iC = 0; iC < channels; iC++)
                {
                    *pfWarpedImagePixel++ = *pfFeatureImage1Pixel++;
                }
                continue;
            }
            CVOF2::BilinearInterpolate(pFeatureImage2->Ptr(), width, height, strideBytes, channels, x, y, pfWarpedImageRow + offset);
        }
    }
}

void CVOF2::BilinearInterpolate(float* image, int width, int height, size_t strideBytes, int channels, double x, double y, float* resultImage)
{
    int xx = (int)x;
    int yy = (int)y;

    //Ensure value is in [0,1]
    float dx = ::VtMax<float>((float)x - xx, 0);
    float dy = ::VtMax<float>((float)y - yy, 0);

    float dxx = 1 - dx;
    float dyy = 1 - dy;

    int previousX = ::VtMin<int>(::VtMax<int>(xx, 0), width - 1);
    int nextX = ::VtMin<int>(::VtMax<int>(xx + 1, 0), width - 1);
    
    int previousY = ::VtMin<int>(::VtMax<int>(yy, 0), height - 1);
    int nextY = ::VtMin<int>(::VtMax<int>(yy + 1, 0), height - 1);

    float* previousRow = (float*)((Byte*)image + previousY * strideBytes);
    float* nextRow = (float*)((Byte*)image + nextY * strideBytes);

#ifdef USE_FEATURES
    int previousXOffset = previousX * channels;
    int nextXOffset = nextX * channels;

    float* f1 = previousRow + previousXOffset;
    float* f2 = previousRow + nextXOffset;
    float* f3 = nextRow + previousXOffset;
    float* f4 = nextRow + nextXOffset;

    for (int iC = 0; iC < channels; iC++)
    {
        *resultImage++ = (float)(dxx * dyy * (*f1++) 
            + dx  * dyy * (*f2++) 
            + dxx * dy  * (*f3++) 
            + dx  * dy  * (*f4++));
    }
#else
    UNREFERENCED_PARAMETER(channels);
    float* f1 = previousRow + previousX;
    float* f2 = previousRow + nextX;
    float* f3 = nextRow + previousX;
    float* f4 = nextRow + nextX;
    *resultImage = dxx * dyy * (*f1) + dx  * dyy * (*f2) + dxx * dy  * (*f3) + dx  * dy  * (*f4);
#endif
}

void CVOF2::SmoothFlowSOR(int iLevel, double alpha, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations)
{
    double varepsilon_phi = pow(0.001, 2);
    double varepsilon_psi = pow(0.001, 2);
    
    //LH:Todo:Rename these variables
    CFloatImg* puu = &(uu.GetLevel(iLevel));
    CFloatImg* pvv = &(vv.GetLevel(iLevel));
    CFloatImg* pux = &(ux.GetLevel(iLevel));
    CFloatImg* puy = &(uy.GetLevel(iLevel));
    CFloatImg* pvx = &(vx.GetLevel(iLevel));
    CFloatImg* pvy = &(vy.GetLevel(iLevel));

    CFloatImg* pPhi_1st = &(Phi_1st.GetLevel(iLevel));
    CFloatImg* pPsi_1st = &(Psi_1st.GetLevel(iLevel));
#ifdef USE_FEATURES
    CFloatImg* pimdxy   = &(imdxy.GetLevel(iLevel));
    CFloatImg* pimdx2   = &(imdx2.GetLevel(iLevel));
    CFloatImg* pimdy2   = &(imdy2.GetLevel(iLevel));
    CFloatImg* pimdtdx  = &(imdtdx.GetLevel(iLevel));
    CFloatImg* pimdtdy  = &(imdtdy.GetLevel(iLevel));
#else
    CFloatImg* pimdxy   = 0;
    CFloatImg* pimdx2   = 0;
    CFloatImg* pimdy2   = 0;
    CFloatImg* pimdtdx  = 0;
    CFloatImg* pimdtdy  = 0;
#endif
    CFloatImg* pImDxy   = &(ImDxy.GetLevel(iLevel));
    CFloatImg* pImDx2   = &(ImDx2.GetLevel(iLevel));
    CFloatImg* pImDy2   = &(ImDy2.GetLevel(iLevel));
    CFloatImg* pImDtDx  = &(ImDtDx.GetLevel(iLevel));
    CFloatImg* pImDtDy  = &(ImDtDy.GetLevel(iLevel));

    CFloatImg* pFlowImageX  = &(m_pyrFlowX.GetLevel(iLevel));
    CFloatImg* pFlowImageY  = &(m_pyrFlowY.GetLevel(iLevel));
    CFloatImg* pDFlowImageX = &(m_pyrDFlowX.GetLevel(iLevel));
    CFloatImg* pDFlowImageY = &(m_pyrDFlowY.GetLevel(iLevel));

#ifdef USE_FEATURES
    CFloatImg* pFeatureImage1 = &(m_ppyrFeatureImages[0].GetLevel(iLevel));
#else
    CFloatImg* pFeatureImage1 = &(m_ppyrInputImages[0].GetLevel(iLevel));
#endif
    CFloatImg* pWarpedImage = &(m_pyrWarpedInput.GetLevel(iLevel));

    CFloatImg* pimdx = &(m_pyrDerivsX.GetLevel(iLevel));
    CFloatImg* pimdy = &(m_pyrDerivsY.GetLevel(iLevel));
    CFloatImg* pimdt = &(m_pyrErrorImage.GetLevel(iLevel));

    CFloatImg* foo1  = &(m_pyrFlowLaplacianX.GetLevel(iLevel));
    CFloatImg* foo2  = &(m_pyrFlowLaplacianY.GetLevel(iLevel));

    int width = pFeatureImage1->Width();
    int height = pFeatureImage1->Height();

    float sorAlpha = (float)alpha * 0.05f;

    for (int outerIteration = 0; outerIteration < nOuterFPIterations; outerIteration++)
    {
        // compute the gradient and error image
        this->ComputeGradients(iLevel); //Smooth version

        //Clear out the derivative flow field
        pDFlowImageX->Clear();
        pDFlowImageY->Clear();
        //for (int innerIteration = 0; innerIteration < nInnerFPIterations; innerIteration++)
        for (int innerIteration = 0; innerIteration < nInnerFPIterations; innerIteration++)
        {
            //Compute the next flow values: uu = u + du
            if (innerIteration == 0)
            {
                /*VT*/pFlowImageX->CopyTo(*puu);
                /*VT*/pFlowImageY->CopyTo(*pvv);
            }
            else
            {
                //Add flow image and derivative flow image together and output accordingly
                /*VT*/::VtAddImages(*puu, *pFlowImageX, *pDFlowImageX);
                /*VT*/::VtAddImages(*pvv, *pFlowImageY, *pDFlowImageY);
            }

            //Compute coarse derivatives for the new flow image value
            this->ComputeDerivativeX(*puu, *pux /*output*/, m_bSmoothDerivatives); //Horizontal flow component
            this->ComputeDerivativeY(*puu, *puy /*output*/, m_bSmoothDerivatives);
            this->ComputeDerivativeX(*pvv, *pvx /*output*/, m_bSmoothDerivatives); //Vertical flow component
            this->ComputeDerivativeY(*pvv, *pvy /*output*/, m_bSmoothDerivatives);

            this->ComputePhi(*pux, *puy, *pvx, *pvy, *pPhi_1st, varepsilon_phi);
            this->ComputePsi(*pimdx, *pimdy, *pimdt, *pDFlowImageX, *pDFlowImageY, *pPsi_1st, m_LaplacianParameters, varepsilon_psi);

            // prepare the components of the large linear system
            this->MultiplyImage(*pPsi_1st, *pimdx, *pimdy, *pImDxy);
            this->MultiplyImage(*pPsi_1st, *pimdx, *pimdx, *pImDx2);
            this->MultiplyImage(*pPsi_1st, *pimdy, *pimdy, *pImDy2);
            this->MultiplyImage(*pPsi_1st, *pimdx, *pimdt, *pImDtDx);
            this->MultiplyImage(*pPsi_1st, *pimdy, *pimdt, *pImDtDy);

#ifdef USE_FEATURES
            this->CollapseImage(*pImDxy,  *pimdxy);
            this->CollapseImage(*pImDx2,  *pimdx2);
            this->CollapseImage(*pImDy2,  *pimdy2);
            this->CollapseImage(*pImDtDx, *pimdtdx);
            this->CollapseImage(*pImDtDy, *pimdtdy);
#else
            pimdxy = pImDxy;
            pimdx2 = pImDx2;
            pimdy2 = pImDy2;
            pimdtdx = pImDtDx;
            pimdtdy = pImDtDy;
#endif
            //Laplacian filtering of the current flow field
            this->ApplyFilters_Laplacian(*foo1, *pFlowImageX, *pPhi_1st);
            this->ApplyFilters_Laplacian(*foo2, *pFlowImageY, *pPhi_1st);
            
#if (defined(_M_IX86) || defined(_M_AMD64))
            if (m_bUseSSE && m_bSSE2)
            {
                int remainder = width % 4;
                int sseLength = width / 4;
                __m128 m128Alpha = ::_mm_set_ps1(float(alpha));
                __m128 m128Zero = ::_mm_set_ps1(0.0f);
                __m128 temp;
                for (int iY = 0; iY < height; iY++)
                {
                    __m128* imdtdxData = (__m128*)pimdtdx->Ptr(iY);
                    __m128* imdtdyData = (__m128*)pimdtdy->Ptr(iY);
                    __m128* foo1Data   = (__m128*)foo1->Ptr(iY);
                    __m128* foo2Data   = (__m128*)foo2->Ptr(iY);
                    for (int i = 0; i < sseLength; i++)
                    {
                        temp = ::_mm_mul_ps(*foo1Data++, m128Alpha);
                        *imdtdxData++ = ::_mm_sub_ps(::_mm_sub_ps(m128Zero, *imdtdxData), temp);
                        temp = ::_mm_mul_ps(*foo2Data++, m128Alpha);
                        *imdtdyData++ = ::_mm_sub_ps(::_mm_sub_ps(m128Zero, *imdtdyData), temp);
                    }
                    if (remainder > 0)
                    {
                        float* pfimdtdxData = (float*)imdtdxData;
                        float* pfimdtdyData = (float*)imdtdyData;
                        float* pffoo1Data   = (float*)foo1Data;
                        float* pffoo2Data   = (float*)foo2Data;
                        for (int i = 0; i < remainder; i++)
                        {
                            *pfimdtdxData++ = - (*pfimdtdxData) - (*pffoo1Data++) * (float)alpha;
                            *pfimdtdyData++ = - (*pfimdtdyData) - (*pffoo2Data++) * (float)alpha;
                        }
                    }
                }
            }
            else
#endif
            {
                for (int iY = 0; iY < height; iY++)
                {
                    float* imdtdxData = pimdtdx->Ptr(iY);
                    float* imdtdyData = pimdtdy->Ptr(iY);
                    float* foo1Data   = foo1->Ptr(iY);
                    float* foo2Data   = foo2->Ptr(iY);
                    float* imdtdxDataEnd = imdtdxData + width;
                    while (imdtdxData < imdtdxDataEnd)
                    {
                        *imdtdxData++ = - (*imdtdxData) - (*foo1Data++) * (float)alpha;
                        *imdtdyData++ = - (*imdtdyData) - (*foo2Data++) * (float)alpha;
                    }
                }
            }

            //here we start SOR
            this->SolveLinearSystems_SOR(iLevel, *pimdxy, *pimdx2, *pimdy2, *pimdtdx, *pimdtdy, alpha, sorAlpha, nSORIterations);
        } //End inner iteration
        
        //Compute new flow field
        /*VT*/::VtAddImages(*pFlowImageX, *pFlowImageX, *pDFlowImageX);
        /*VT*/::VtAddImages(*pFlowImageY, *pFlowImageY, *pDFlowImageY);

        //Now warp image again
        this->WarpImage(iLevel);

        switch (m_iNoiseModel)
        {
            case GaussianMixture:
                break;
            case Laplacian:
                this->EstimateLaplacianNoise(*pFeatureImage1, *pWarpedImage, m_LaplacianParameters);
                break;
        }
    }//End outer iteration
}

void CVOF2::SolveLinearSystems_SOR(int iLevel, CFloatImg& pimdxy, CFloatImg& pimdx2, CFloatImg& pimdy2, CFloatImg& pimdtdx, CFloatImg& pimdtdy, double alpha, double sorAlpha, int nSORIterations)
{
    if (m_bUseSSE && m_bSSE2)
    {
        this->SolveLinearSystems_SOR_SSE(iLevel, pimdxy, pimdx2, pimdy2, pimdtdx, pimdtdy, alpha, sorAlpha, nSORIterations);
    }
    else
    {
        this->SolveLinearSystems_SOR_NoSSE(iLevel, pimdxy, pimdx2, pimdy2, pimdtdx, pimdtdy, alpha, sorAlpha, nSORIterations);
    }
}

void CVOF2::SolveLinearSystems_SOR_SSE(int iLevel, CFloatImg& pimdxy, CFloatImg& pimdx2, CFloatImg& pimdy2, CFloatImg& pimdtdx, CFloatImg& pimdtdy, double alpha, double sorAlpha, int nSORIterations)
{
#if (defined(_M_IX86) || defined(_M_AMD64))
   CFloatImg* pPhi_1st = &(Phi_1st.GetLevel(iLevel));
    CFloatImg* pDFlowImageX = &(m_pyrDFlowX.GetLevel(iLevel));
    CFloatImg* pDFlowImageY = &(m_pyrDFlowY.GetLevel(iLevel));
    int height = pPhi_1st->Height();
    int width = pPhi_1st->Width();

    float flippedOmega = 1 - m_Omega;
    int remainder = width % 4;
    int sseLength = width / 4 - 1;
    double sigma1 = 0, sigma2 = 0, coeff = 0;

    __m128  m128DFlowXNext;
    __m128  m128DFlowYNext;
            
    __m128  m128Alpha = ::_mm_set_ps1((float)alpha);
    __m128  m128NegativeAlpha = ::_mm_set_ps1((float)-alpha);
    __m128  m128Omega = ::_mm_set_ps1(m_Omega);
    __m128  m128FlippedOmega = ::_mm_set_ps1(flippedOmega);
    __m128  m128SorAlpha = ::_mm_set_ps1((float)sorAlpha);

    __m128  m128Phi, m128PhiAbove, m128PhiPrevious;
    __m128  m128Sigma1, m128Sigma2, m128Coeff;
    __m128  m128Temp1, m128Temp2, m128Temp3, m128Temp4;

    for (int sorIteration = 0; sorIteration < nSORIterations; sorIteration++)
    {
        float* pfPhiDataRow = pPhi_1st->Ptr();
        float* pfDFlowXRow  = pDFlowImageX->Ptr();
        float* pfDFlowYRow  = pDFlowImageY->Ptr();

        float* pfPhiDataRowAbove = 0;
        float* pfDFlowXRowAbove  = 0;
        float* pfDFlowYRowAbove  = 0;

        float* pfDFlowXRowBelow  = pDFlowImageX->Ptr(1);
        float* pfDFlowYRowBelow  = pDFlowImageY->Ptr(1);

        float* pfimdxyRow = pimdxy.Ptr();
        float* pfimdx2Row = pimdx2.Ptr();
        float* pfimdy2Row = pimdy2.Ptr();
        float* pfimdtdxRow = pimdtdx.Ptr();
        float* pfimdtdyRow = pimdtdy.Ptr();

        float fPhiPrevious = 0.0f;
        float fPhi = *pfPhiDataRow;
        float fPhiAbove = 0.0f;

        ///-------------------------------------------------------------////
        ///--------------Process first row------------------------------////
        ///-------------------------------------------------------------////
        for (int iX = 0; iX < width; iX++)
        {
            sigma1 = 0; sigma2 = 0; coeff = 0;
            fPhi = *pfPhiDataRow;
            if (iX > 0)
            {
                fPhiPrevious = *(pfPhiDataRow - 1);
                sigma1 += fPhiPrevious * (*(pfDFlowXRow - 1));
                sigma2 += fPhiPrevious * (*(pfDFlowYRow - 1));
                coeff  += fPhiPrevious;
            }
            if (iX < width - 1)
            {
                sigma1 += fPhi * (*(pfDFlowXRow + 1));
                sigma2 += fPhi * (*(pfDFlowYRow + 1));
                coeff  += fPhi;
            }
            sigma1 += fPhi * (*pfDFlowXRowBelow++);
            sigma2 += fPhi * (*pfDFlowYRowBelow++);
            coeff  += fPhi;
            sigma1 *= (-alpha);
            sigma2 *= (-alpha);
            coeff  *= alpha;
            //Compute horizontal differential flow change
            sigma1 += *pfimdxyRow * (*pfDFlowYRow);
            *pfDFlowXRow = (float)(flippedOmega * (*pfDFlowXRow)
                + m_Omega / ((*pfimdx2Row++) + sorAlpha + coeff) * (*pfimdtdxRow++ - sigma1));
            //Compute vertical differential flow change
            sigma2 += (*pfimdxyRow++) * (*pfDFlowXRow++);
            *pfDFlowYRow++ = (float)(flippedOmega * (*pfDFlowYRow)
                + m_Omega / (*pfimdy2Row++ + sorAlpha + coeff) * (*pfimdtdyRow++ - sigma2));
            pfPhiDataRow++;
        }

        ///-------------------------------------------------------------////
        ///--------------Process middle rows----------------------------////
        ///-------------------------------------------------------------////
        for (int iY = 1; iY < height - 1; iY++)
        {
            pfPhiDataRow = pPhi_1st->Ptr(iY);
            pfDFlowXRow  = pDFlowImageX->Ptr(iY);
            pfDFlowYRow  = pDFlowImageY->Ptr(iY);

            pfPhiDataRowAbove = pPhi_1st->Ptr(iY - 1);
            pfDFlowXRowAbove  = pDFlowImageX->Ptr(iY - 1);
            pfDFlowYRowAbove  = pDFlowImageY->Ptr(iY - 1);

            pfDFlowXRowBelow  = pDFlowImageX->Ptr(iY + 1);
            pfDFlowYRowBelow  = pDFlowImageY->Ptr(iY + 1);

            pfimdxyRow = pimdxy.Ptr(iY);
            pfimdx2Row = pimdx2.Ptr(iY);
            pfimdy2Row = pimdy2.Ptr(iY);
            pfimdtdxRow = pimdtdx.Ptr(iY);
            pfimdtdyRow = pimdtdy.Ptr(iY);

            fPhi = *pfPhiDataRow;
            fPhiAbove = *pfPhiDataRowAbove++;
            //-----------Process first 4 pixels-----------------------//
            sigma1 = fPhi * (*(pfDFlowXRow + 1)) +
                     fPhiAbove * (*(pfDFlowXRowAbove++)) +
                     fPhi * (*pfDFlowXRowBelow++);
            sigma2 = fPhi * (*(pfDFlowYRow + 1)) +
                     fPhiAbove * (*pfDFlowYRowAbove++) +
                     fPhi * (*pfDFlowYRowBelow++);
            coeff  = fPhi + fPhiAbove + fPhi;
            sigma1 *= (-alpha);
            sigma2 *= (-alpha);
            coeff  *= alpha;
            //Compute horizontal differential flow change
            sigma1 += *pfimdxyRow * (*pfDFlowYRow);
            *pfDFlowXRow = (float)(flippedOmega * (*pfDFlowXRow)
                + m_Omega / ((*pfimdx2Row++) + sorAlpha + coeff) * (*pfimdtdxRow++ - sigma1));
            //Compute vertical differential flow change
            sigma2 += (*pfimdxyRow++) * (*pfDFlowXRow++);
            *pfDFlowYRow++ = (float)(flippedOmega * (*pfDFlowYRow)
                + m_Omega / (*pfimdy2Row++ + sorAlpha + coeff) * (*pfimdtdyRow++ - sigma2));
            pfPhiDataRow++;
            for (int iX = 1; iX < 4; iX++)
            {
                fPhiPrevious = *(pfPhiDataRow - 1);
                fPhi = *pfPhiDataRow;
                fPhiAbove = *pfPhiDataRowAbove++;

                sigma1 = fPhiPrevious * (*(pfDFlowXRow - 1)) +
                         fPhi * (*(pfDFlowXRow + 1)) +
                         fPhiAbove * (*pfDFlowXRowAbove++) +
                         fPhi * (*pfDFlowXRowBelow++);
                sigma2 = fPhiPrevious * (*(pfDFlowYRow - 1)) +
                         fPhi * (*(pfDFlowYRow + 1)) +
                         fPhiAbove * (*pfDFlowYRowAbove++) +
                         fPhi * (*pfDFlowYRowBelow++);
                coeff  = fPhiPrevious + fPhi + fPhi + fPhiAbove;
                sigma1 *= (-alpha);
                sigma2 *= (-alpha);
                coeff  *= alpha;
                //Compute horizontal differential flow change
                sigma1 += *pfimdxyRow * (*pfDFlowYRow);
                *pfDFlowXRow = (float)(flippedOmega * (*pfDFlowXRow)
                    + m_Omega / ((*pfimdx2Row++) + sorAlpha + coeff) * (*pfimdtdxRow++ - sigma1));
                //Compute vertical differential flow change
                sigma2 += (*pfimdxyRow++) * (*pfDFlowXRow++);
                *pfDFlowYRow++ = (float)(flippedOmega * (*pfDFlowYRow)
                    + m_Omega / (*pfimdy2Row++ + sorAlpha + coeff) * (*pfimdtdyRow++ - sigma2));
                pfPhiDataRow++;
            } //End one pixel

            //-----------Process main pixels--------------------------//
            __m128* pm128Phi = (__m128*)pfPhiDataRow;
            __m128* pm128PhiAbove = (__m128*)pfPhiDataRowAbove;
            
            __m128* pm128DFlowX      = (__m128*)pfDFlowXRow;
            __m128* pm128DFlowXAbove = (__m128*)pfDFlowXRowAbove;
            __m128* pm128DFlowXBelow = (__m128*)pfDFlowXRowBelow;
            
            __m128* pm128DFlowY      = (__m128*)pfDFlowYRow;
            __m128* pm128DFlowYAbove = (__m128*)pfDFlowYRowAbove;
            __m128* pm128DFlowYBelow = (__m128*)pfDFlowYRowBelow;

            __m128* pm128Imdxy  = (__m128*)pfimdxyRow;
            __m128* pm128Imdx2  = (__m128*)pfimdx2Row;
            __m128* pm128Imdtdx = (__m128*)pfimdtdxRow;
            __m128* pm128Imdy2  = (__m128*)pfimdy2Row;
            __m128* pm128Imdtdy = (__m128*)pfimdtdyRow;

            float *pfTmpX = m_imgSORX.Ptr() + 4;
            float *pfTmpY = m_imgSORY.Ptr() + 4;
            float *pfTmpPhi = m_imgSORPhi.Ptr() + 4;
            __m128* pm128TmpX = (__m128*)pfTmpX;
            __m128* pm128TmpY = (__m128*)pfTmpY;
            __m128* pm128TmpPhi = (__m128*)pfTmpPhi;

            for (int i = 1; i < sseLength ; i++)
            {
                //fPhiPrevious = *(pfPhiDataRow - 1);
                //fPhi = *pfPhiDataRow;
                //fPhiAbove = *(pfPhiDataRowAbove + iX);
                m128Phi = *pm128Phi++;
                m128PhiAbove = *pm128PhiAbove++;
                m128PhiPrevious = ::_mm_set_ps(*(pfPhiDataRow + 2), *(pfPhiDataRow + 1), *pfPhiDataRow, *(pfPhiDataRow - 1));

                *pm128TmpPhi++ = ::_mm_mul_ps(m128PhiPrevious, m128NegativeAlpha);

                m128DFlowXNext = ::_mm_set_ps(*(pfDFlowXRow + 4), *(pfDFlowXRow + 3), *(pfDFlowXRow + 2), *(pfDFlowXRow + 1));
                //sigma1 = fPhi * (*(pfDFlowXRow + 1)) +
                //		 fPhiAbove * (*(pfDFlowXRowAbove + iX)) +
                //		 fPhi * (*(pfDFlowXRowBelow + iX));
                //m128Temp1 = ::_mm_mul_ps(m128PhiPrevious, m128DFlowXPrevious);
                m128Temp2 = ::_mm_mul_ps(m128Phi, m128DFlowXNext);
                m128Temp3 = ::_mm_mul_ps(m128PhiAbove, *pm128DFlowXAbove++);
                m128Temp4 = ::_mm_mul_ps(m128Phi, *pm128DFlowXBelow++);

                m128Temp1 = ::_mm_add_ps(m128Temp3, m128Temp4);
                //m128Temp2 = ::_mm_add_ps(m128Temp3, m128Temp4);

                m128Sigma1 = ::_mm_add_ps(m128Temp1, m128Temp2);

                //sigma2 = fPhi * (*(pfDFlowYRow + 1)) +
                //		 fPhiAbove * (*(pfDFlowYRowAbove + iX)) +
                //		 fPhi * (*(pfDFlowYRowBelow + iX));
                m128DFlowYNext = ::_mm_set_ps(*(pfDFlowYRow + 4), *(pfDFlowYRow + 3), *(pfDFlowYRow + 2), *(pfDFlowYRow + 1));
                m128Temp2 = ::_mm_mul_ps(m128Phi, m128DFlowYNext);
                m128Temp3 = ::_mm_mul_ps(m128PhiAbove, *pm128DFlowYAbove++);
                m128Temp4 = ::_mm_mul_ps(m128Phi, *pm128DFlowYBelow++);
                m128Temp1 = ::_mm_add_ps(m128Temp3, m128Temp4);
                m128Sigma2 = ::_mm_add_ps(m128Temp1, m128Temp2);

                //coeff  = fPhiPrevious + fPhi + fPhi + fPhiAbove;
                m128Temp1 = ::_mm_add_ps(m128PhiPrevious, m128Phi);
                m128Temp2 = ::_mm_add_ps(m128PhiAbove, m128Phi);
                m128Coeff  = ::_mm_add_ps(m128Temp1, m128Temp2);
                
                //sigma1 *= (-alpha);
                //sigma2 *= (-alpha);
                //coeff  *= alpha;
                m128Sigma1 = ::_mm_mul_ps(m128Sigma1, m128NegativeAlpha);
                m128Sigma2 = ::_mm_mul_ps(m128Sigma2, m128NegativeAlpha);
                m128Coeff  = ::_mm_mul_ps(m128Coeff,  m128Alpha);
                
                //Compute horizontal differential flow change
                //sigma1 += *pfimdxyRow * (*pfDFlowYRow);
                m128Temp1 = ::_mm_mul_ps(*pm128Imdxy, *pm128DFlowY);
                m128Sigma1 = ::_mm_add_ps(m128Sigma1, m128Temp1);

                //*pfDFlowXRow = flippedOmega * (*pfDFlowXRow)
                //	+ m_Omega / ((*pfimdx2Row++) + sorAlpha + coeff) * (*pfimdtdxRow++ - sigma1);
                m128Temp1 = ::_mm_mul_ps(m128FlippedOmega, *pm128DFlowX);
                m128Temp2 = ::_mm_add_ps(m128Coeff, ::_mm_add_ps(*pm128Imdx2++, m128SorAlpha));
                *pm128TmpX = ::_mm_div_ps(m128Omega, m128Temp2);
                m128Temp4 = ::_mm_mul_ps(*pm128TmpX++, ::_mm_sub_ps(*pm128Imdtdx++, m128Sigma1));
                *pm128DFlowX = ::_mm_add_ps(m128Temp1, m128Temp4);
                
                //Compute vertical differential flow change
                //sigma2 += (*pfimdxyRow++) * (*pfDFlowXRow++);
                m128Temp1 = ::_mm_mul_ps(*pm128Imdxy++, *pm128DFlowX++);
                m128Sigma2 = ::_mm_add_ps(m128Sigma2, m128Temp1);

                //*pfDFlowYRow++ = flippedOmega * (*pfDFlowYRow)
                //	+ m_Omega / (*pfimdy2Row++ + sorAlpha + coeff) * (*pfimdtdyRow++ - sigma2);
                m128Temp1 = ::_mm_mul_ps(m128FlippedOmega, *pm128DFlowY);
                m128Temp2 = ::_mm_add_ps(m128Coeff, ::_mm_add_ps(*pm128Imdy2++, m128SorAlpha));
                *pm128TmpY = ::_mm_div_ps(m128Omega, m128Temp2);
                m128Temp4 = ::_mm_mul_ps(*pm128TmpY++, ::_mm_sub_ps(*pm128Imdtdy++, m128Sigma2));
                *pm128DFlowY++ = ::_mm_add_ps(m128Temp1, m128Temp4);

                pfPhiDataRow += 4;
                pfDFlowXRow  += 4;
                pfDFlowYRow  += 4;
            }
            pfPhiDataRow = pPhi_1st->Ptr(iY) + 4;
            pfDFlowXRow  = pDFlowImageX->Ptr(iY) + 4;
            pfDFlowYRow  = pDFlowImageY->Ptr(iY) + 4;
            float dFlowXIncrement;
            float* pfDFlowXPrevious = pfDFlowXRow - 1;
            float* pfDFlowYPrevious = pfDFlowYRow - 1;
            for (int iX = 4; iX < width - remainder - 4; iX++)
            {
                sigma1 = *pfTmpPhi * (*pfDFlowXPrevious++);
                sigma2 = *pfTmpPhi * (*pfDFlowYPrevious++);
                dFlowXIncrement = (float)(*pfTmpX++ * (-sigma1));
                *pfDFlowXRow++ += dFlowXIncrement;
                sigma2 += (*pfimdxyRow++) * dFlowXIncrement;
                *pfDFlowYRow++ += (float)(*pfTmpY++ * (-sigma2));
                
                pfPhiDataRow++;
                pfimdtdxRow++;
                pfimdtdyRow++;
                pfimdx2Row++;
                pfimdy2Row++;
                pfTmpPhi++;
            }
            pfDFlowXRowAbove  = (float*)pm128DFlowXAbove;
            pfDFlowXRowBelow  = (float*)pm128DFlowXBelow;
            pfDFlowYRowAbove  = (float*)pm128DFlowYAbove;
            pfDFlowYRowBelow  = (float*)pm128DFlowYBelow;
            pfPhiDataRowAbove = (float*)pm128PhiAbove;
            //-----------Process last remainder pixels------------------------//
            for (int iX = width - remainder - 4; iX < width - 1; iX++)
            {
                fPhiPrevious = *(pfPhiDataRow - 1);
                fPhi = *pfPhiDataRow;
                fPhiAbove = *pfPhiDataRowAbove++;
                sigma1 = fPhiPrevious * (*(pfDFlowXRow - 1)) +
                         fPhi * (*(pfDFlowXRow + 1)) +
                         fPhiAbove * (*pfDFlowXRowAbove++) +
                         fPhi * (*pfDFlowXRowBelow++);
                sigma2 = fPhiPrevious * (*(pfDFlowYRow - 1)) +
                         fPhi * (*(pfDFlowYRow + 1)) +
                         fPhiAbove * (*pfDFlowYRowAbove++) +
                         fPhi * (*pfDFlowYRowBelow++);
                coeff  = fPhiPrevious + fPhi + fPhi + fPhiAbove;
                sigma1 *= (-alpha);
                sigma2 *= (-alpha);
                coeff  *= alpha;
                //Compute horizontal differential flow change
                sigma1 += *pfimdxyRow * (*pfDFlowYRow);
                *pfDFlowXRow = (float)(flippedOmega * (*pfDFlowXRow)
                    + m_Omega / ((*pfimdx2Row++) + sorAlpha + coeff) * (*pfimdtdxRow++ - sigma1));
                //Compute vertical differential flow change
                sigma2 += (*pfimdxyRow++) * (*pfDFlowXRow++);
                *pfDFlowYRow++ = (float)(flippedOmega * (*pfDFlowYRow)
                    + m_Omega / (*pfimdy2Row++ + sorAlpha + coeff) * (*pfimdtdyRow++ - sigma2));
                pfPhiDataRow++;
            }
            fPhiPrevious = *(pfPhiDataRow - 1);
            fPhi = *pfPhiDataRow;
            fPhiAbove = *pfPhiDataRowAbove;
            sigma1 = fPhiPrevious * (*(pfDFlowXRow - 1)) +
                     fPhiAbove * (*pfDFlowXRowAbove) +
                     fPhi * (*pfDFlowXRowBelow);
            sigma2 = fPhiPrevious * (*(pfDFlowYRow - 1)) +
                     fPhiAbove * (*pfDFlowYRowAbove) +
                     fPhi * (*pfDFlowYRowBelow);
            coeff  = fPhiPrevious + fPhiAbove + fPhi;
            sigma1 *= (-alpha);
            sigma2 *= (-alpha);
            coeff  *= alpha;
            //Compute horizontal differential flow change
            sigma1 += *pfimdxyRow * (*pfDFlowYRow);
            *pfDFlowXRow = (float)(flippedOmega * (*pfDFlowXRow)
                + m_Omega / ((*pfimdx2Row++) + sorAlpha + coeff) * (*pfimdtdxRow++ - sigma1));
            //Compute vertical differential flow change
            sigma2 += (*pfimdxyRow++) * (*pfDFlowXRow++);
            *pfDFlowYRow++ = (float)(flippedOmega * (*pfDFlowYRow)
                + m_Omega / (*pfimdy2Row++ + sorAlpha + coeff) * (*pfimdtdyRow++ - sigma2));
            pfPhiDataRow++;
        } //End one row

        ///-------------------------------------------------------------////
        ///--------------Process last row-------------------------------////
        ///-------------------------------------------------------------////
        int iY = height - 1;
        pfPhiDataRow = pPhi_1st->Ptr(iY);
        pfDFlowXRow  = pDFlowImageX->Ptr(iY);
        pfDFlowYRow  = pDFlowImageY->Ptr(iY);

        pfPhiDataRowAbove = pPhi_1st->Ptr(iY - 1);
        pfDFlowXRowAbove  = pDFlowImageX->Ptr(iY - 1);
        pfDFlowYRowAbove  = pDFlowImageY->Ptr(iY - 1);

        pfimdxyRow = pimdxy.Ptr(iY);
        pfimdx2Row = pimdx2.Ptr(iY);
        pfimdy2Row = pimdy2.Ptr(iY);
        pfimdtdxRow = pimdtdx.Ptr(iY);
        pfimdtdyRow = pimdtdy.Ptr(iY);

        for (int iX = 0; iX < width; iX++)
        {
            fPhi = *pfPhiDataRow;
            fPhiAbove = *pfPhiDataRowAbove++;
            sigma1 = 0; sigma2 = 0; coeff = 0;
            if (iX > 0)
            {
                fPhiPrevious = *(pfPhiDataRow - 1);
                sigma1 += fPhiPrevious * (*(pfDFlowXRow - 1));
                sigma2 += fPhiPrevious * (*(pfDFlowYRow - 1));
                coeff  += fPhiPrevious;
            }
            if (iX < width - 1)
            {
                sigma1 += fPhi * (*(pfDFlowXRow + 1));
                sigma2 += fPhi * (*(pfDFlowYRow + 1));
                coeff  += fPhi;
            }
            sigma1 += fPhiAbove * (*(pfDFlowXRowAbove + iX));
            sigma2 += fPhiAbove * (*(pfDFlowYRowAbove + iX));
            coeff  += fPhiAbove;
            sigma1 *= (-alpha);
            sigma2 *= (-alpha);
            coeff  *= alpha;
            //Compute horizontal differential flow change
            sigma1 += *pfimdxyRow * (*pfDFlowYRow);
            *pfDFlowXRow = (float)(flippedOmega * (*pfDFlowXRow)
                + m_Omega / ((*pfimdx2Row++) + sorAlpha + coeff) * (*pfimdtdxRow++ - sigma1));
            //Compute vertical differential flow change
            sigma2 += (*pfimdxyRow++) * (*pfDFlowXRow++);
            *pfDFlowYRow++ = (float)(flippedOmega * (*pfDFlowYRow)
                + m_Omega / (*pfimdy2Row++ + sorAlpha + coeff) * (*pfimdtdyRow++ - sigma2));

            pfPhiDataRow++;
        } //End one pixel

    } //End SOR iteration
#else
    iLevel,pimdxy,pimdx2,pimdy2,pimdtdx,pimdtdy,alpha,sorAlpha,nSORIterations;
#endif
}

void CVOF2::SolveLinearSystems_SOR_NoSSE(int iLevel, CFloatImg& pimdxy, CFloatImg& pimdx2, CFloatImg& pimdy2, CFloatImg& pimdtdx, CFloatImg& pimdtdy, double alpha, double sorAlpha, int nSORIterations)
{
    CFloatImg* pPhi_1st = &(Phi_1st.GetLevel(iLevel));
    CFloatImg* pDFlowImageX = &(m_pyrDFlowX.GetLevel(iLevel));
    CFloatImg* pDFlowImageY = &(m_pyrDFlowY.GetLevel(iLevel));
    int height = pPhi_1st->Height();
    int width = pPhi_1st->Width();
    double sigma1 = 0, sigma2 = 0, coeff = 0;
    float oneMinusOmega = 1 - m_Omega;
    for (int sorIteration = 0; sorIteration < nSORIterations; sorIteration++)
    {
        for (int iY = 0; iY < height; iY++)
        {
            float* pfPhiDataRow = pPhi_1st->Ptr(iY);
            float* pfDFlowXRow  = pDFlowImageX->Ptr(iY);
            float* pfDFlowYRow  = pDFlowImageY->Ptr(iY);

            float* pfPhiDataRowAbove = (iY > 0) ? pPhi_1st->Ptr(iY - 1)     : 0;
            float* pfDFlowXRowAbove  = (iY > 0) ? pDFlowImageX->Ptr(iY - 1) : 0;
            float* pfDFlowYRowAbove  = (iY > 0) ? pDFlowImageY->Ptr(iY - 1) : 0;

            float* pfDFlowXRowBelow  = (iY < height - 1) ? pDFlowImageX->Ptr(iY + 1) : 0;
            float* pfDFlowYRowBelow  = (iY < height - 1) ? pDFlowImageY->Ptr(iY + 1) : 0;

            float* pfimdxyRow = pimdxy.Ptr(iY);
            float* pfimdx2Row = pimdx2.Ptr(iY);
            float* pfimdy2Row = pimdy2.Ptr(iY);
            float* pfimdtdxRow = pimdtdx.Ptr(iY);
            float* pfimdtdyRow = pimdtdy.Ptr(iY);

            for (int iX = 0; iX < width; iX++)
            {
                sigma1 = 0; sigma2 = 0; coeff = 0;
                if (iY == 0 || iY == height - 1 || iX == 0 || iX == width - 1)
                {
                    if (iX > 0)
                    {
                        //LH:Todo:Consider optimizing calculations here
                        sigma1 += *(pfPhiDataRow - 1) * (*(pfDFlowXRow - 1));
                        sigma2 += *(pfPhiDataRow - 1) * (*(pfDFlowYRow - 1));
                        coeff  += *(pfPhiDataRow - 1);
                    }

                    if (iX < width - 1)
                    {
                        //LH:Todo:Consider optimizing calculations here
                        sigma1 += *(pfPhiDataRow) * (*(pfDFlowXRow + 1));
                        sigma2 += *(pfPhiDataRow) * (*(pfDFlowYRow + 1));
                        coeff  += *(pfPhiDataRow);
                    }

                    if (iY > 0)
                    {
                        //LH:Todo:Consider optimizing calculations here
                        sigma1 += *(pfPhiDataRowAbove + iX) * (*(pfDFlowXRowAbove + iX));
                        sigma2 += *(pfPhiDataRowAbove + iX) * (*(pfDFlowYRowAbove + iX));
                        coeff  += *(pfPhiDataRowAbove + iX);
                    }
                            
                    if (iY < height - 1)
                    {
						ANALYZE_ASSUME( pfDFlowXRowBelow != NULL );
						ANALYZE_ASSUME( pfDFlowYRowBelow != NULL );

                        //LH:Todo:Consider optimizing calculations here
                        sigma1 += *(pfPhiDataRow) * (*(pfDFlowXRowBelow + iX));
                        sigma2 += *(pfPhiDataRow) * (*(pfDFlowYRowBelow + iX));
                        coeff  += *(pfPhiDataRow);
                    }
                }
                else
                {
                    float fPhiPrevious = *(pfPhiDataRow - 1);
                    float fPhi = *pfPhiDataRow;
                    float fPhiAbove = *(pfPhiDataRowAbove + iX);
                    
                    float fDFlowXPrevious = *(pfDFlowXRow - 1);
                    float fDFlowXNext = *(pfDFlowXRow + 1);
                    float fDFlowXAbove = *(pfDFlowXRowAbove + iX);
                    float fDFlowXBelow = *(pfDFlowXRowBelow + iX);

                    float fDFlowYNext = *(pfDFlowYRow + 1);
                    float fDFlowYPrevious = *(pfDFlowYRow - 1);
                    float fDFlowYAbove = *(pfDFlowYRowAbove + iX);
                    float fDFlowYBelow = *(pfDFlowYRowBelow + iX);

                    sigma1 += fPhiPrevious * fDFlowXPrevious + 
                              fPhi * fDFlowXNext +
                              fPhiAbove * fDFlowXAbove +
                              fPhi * fDFlowXBelow;
                    sigma2 += fPhiPrevious * fDFlowYPrevious +
                              fPhi * fDFlowYNext +
                              fPhiAbove * fDFlowYAbove +
                              fPhi * fDFlowYBelow;
                    coeff +=  fPhiPrevious + fPhi + fPhiAbove + fPhi;
                }

                sigma1 *= (-alpha);
                sigma2 *= (-alpha);
                coeff  *= alpha;
                //Compute horizontal differential flow change
                sigma1 += *pfimdxyRow * (*pfDFlowYRow);
                *pfDFlowXRow = (float)(oneMinusOmega * (*pfDFlowXRow)
                    + m_Omega / ((*pfimdx2Row++) + sorAlpha + coeff) * (*pfimdtdxRow++ - sigma1));

                //Compute vertical differential flow change
                sigma2 += (*pfimdxyRow++) * (*pfDFlowXRow++);
                *pfDFlowYRow++ = (float)(oneMinusOmega * (*pfDFlowYRow)
                    + m_Omega / (*pfimdy2Row++ + sorAlpha + coeff) * (*pfimdtdyRow++ - sigma2));

                pfPhiDataRow++;
            } //End one pixel
        } //End one row
    } //End SOR iteration
}

void CVOF2::ComputeGradients(int iLevel /* = 0 */)
{
#ifdef USE_FEATURES
    CFloatImg* pFeatureImage1 = &(m_ppyrFeatureImages[0].GetLevel(iLevel));
#else
    CFloatImg* pFeatureImage1 = &(m_ppyrInputImages[0].GetLevel(iLevel));
#endif
    CFloatImg* pWarpedImage = &(m_pyrWarpedInput.GetLevel(iLevel));
    CFloatImg* pErrorImage = &(m_pyrErrorImage.GetLevel(iLevel));
    CFloatImg* pTempImage1 = &m_SmoothFlowTempImage1;
    CFloatImg* pTempImage2 = &m_SmoothFlowTempImage2;
    CFloatImg* pDerivativeImageX = &(m_pyrDerivsX.GetLevel(iLevel));
    CFloatImg* pDerivativeImageY = &(m_pyrDerivsY.GetLevel(iLevel));

    int width = pFeatureImage1->Width();
    int height = pFeatureImage1->Height();
    CRect tempShareRect(0, 0, width, height);
    if (m_bGaussianGradients)
    {
        //float gaussianFilter[5] = { 0.02, 0.11, 0.74, 0.11, 0.02 };
        float gaussianFilter[3] = { 0.1f, 0.8f, 0.1f};
        C1dKernel gaussianKernel;
        gaussianKernel.Create(3, 1, gaussianFilter);

        CFloatImg tempShareImage1;
        pTempImage1->Share(tempShareImage1, &tempShareRect);
        
        //Apply the smooth filter to the first feature image and output to the first temp image
        /*VT*/::VtSeparableFilter(tempShareImage1, tempShareRect, *pFeatureImage1, CPoint(0,0), gaussianKernel, gaussianKernel);

        CFloatImg tempShareImage2;
        pTempImage2->Share(tempShareImage2, &tempShareRect);
        //Apply the smooth filter to the warped image and output to the second temp image
        /*VT*/::VtSeparableFilter(tempShareImage2, tempShareRect, *pWarpedImage, CPoint(0,0), gaussianKernel, gaussianKernel);

        //Compute the error image between filtered warped image and filtered feature image 1
        ::VtSubImages(*pErrorImage, tempShareImage2, tempShareImage1);

        //Multiply m_smoothFlowTempImage1 by 0.4 and multiply m_smoothFlowTempImage2 by 0.6 then
        //add them together and output back into the first temp image.
        ::VtBlendImages(tempShareImage1, tempShareImage1, tempShareImage2, 0.4f, 0.6f);
    }
    else
    {
        /*VT*/::VtSubImages(*pErrorImage, *pWarpedImage, *pFeatureImage1);

        CFloatImg tempShareImage1;
        pTempImage1->Share(tempShareImage1, &tempShareRect);
        /*VT*/::VtBlendImages(tempShareImage1, *pFeatureImage1, *pWarpedImage, 0.4f, 0.6f);
    }
    //Now compute the derivative images for the gaussian smoothed image.
    this->ComputeDerivativeX(*pTempImage1, *pDerivativeImageX, m_bSmoothDerivatives, width, height);
    this->ComputeDerivativeY(*pTempImage1, *pDerivativeImageY, m_bSmoothDerivatives, width, height);
}

void CVOF2::ApplyFilters_Laplacian(CFloatImg& output, CFloatImg& input, CFloatImg& weight)
{
    if (m_bUseSSE && m_bSSE2)
    {
        this->ApplyFilters_LaplacianSSE(output, input, weight);
    }
    else
    {
        this->ApplyFilters_LaplacianNoSSE(output, input, weight);
    }
}

void CVOF2::ApplyFilters_LaplacianSSE(CFloatImg& output, CFloatImg& input, CFloatImg& weight)
{
#if (defined(_M_IX86) || defined(_M_AMD64))
    int width = input.Width();
    int height = input.Height();
    
    CFloatImg* tempImage = &m_TempImage;
    //Horizontal filtering
    int remainder = width % 4;
    remainder = (remainder == 0) ? 4 : remainder;
    int sseLength = (width - remainder) / 4;
    for (int iY = 0; iY < height; iY++)
    {
        __m128* m128TempRow = (__m128*)tempImage->Ptr(iY);
        __m128* m128InputRow = (__m128*)input.Ptr(iY);
        __m128* m128WeightRow = (__m128*)weight.Ptr(iY);
        __m128 pfInputRowNext;
        float* pfInputRow = input.Ptr(iY);
        for (int i = 0; i < sseLength; i++)
        {
            pfInputRowNext = ::_mm_set_ps(*(pfInputRow + 4), *(pfInputRow + 3), *(pfInputRow + 2), *(pfInputRow + 1));
            *m128TempRow++ = ::_mm_mul_ps(*m128WeightRow++, ::_mm_sub_ps(pfInputRowNext, *m128InputRow++));
            pfInputRow += 4;
        }
        float* pfTempRow = (float*)m128TempRow;
        float* pfWeightRow = (float*)m128WeightRow;
        for (int i = 0; i < remainder; i++)
        {
            *pfTempRow++ = (*(pfInputRow + 1) - (*pfInputRow++)) * (*pfWeightRow++);
        }
    }
    for (int iY = 0; iY < height; iY++)
    {
        float* pfTempRow = tempImage->Ptr(iY);
        float* pfOutputRow = output.Ptr(iY);
        
        *pfOutputRow++ = -(*pfTempRow++); //Set first 4 pixels
        *pfOutputRow++ = *(pfTempRow - 1) - (*pfTempRow++);
        *pfOutputRow++ = *(pfTempRow - 1) - (*pfTempRow++);
        *pfOutputRow++ = *(pfTempRow - 1) - (*pfTempRow++);

        __m128* m128OutputRow = (__m128*)pfOutputRow;
        __m128* m128TempRow = (__m128*)pfTempRow;
        __m128  m128TempRowPrev;
        for (int i = 0; i < sseLength - 1; i++) //Set middle pixels
        {
            m128TempRowPrev = ::_mm_set_ps(*(pfTempRow + 2), *(pfTempRow + 1), *pfTempRow, *(pfTempRow - 1));
            *m128OutputRow++ = ::_mm_sub_ps(m128TempRowPrev, *m128TempRow++);
            pfTempRow += 4;
        }
        pfOutputRow = (float*)m128OutputRow;
        for (int i = 0; i < remainder - 1; i++)
        {
            *pfOutputRow++ = *(pfTempRow - 1) - (*pfTempRow++);
        }
        *pfOutputRow = *(pfTempRow - 1); //Set last pixel
    }

    //Vertical filtering
    for (int iY = 0; iY < height - 1; iY++)
    {
        __m128* m128TempRow = (__m128*)tempImage->Ptr(iY);
        __m128* m128InputRow = (__m128*)input.Ptr(iY);
        __m128* m128InputRowBelow = (__m128*)input.Ptr(iY + 1);
        __m128* m128WeightRow = (__m128*)weight.Ptr(iY);
        for (int i = 0; i < sseLength; i++)
        {
            *m128TempRow++ = ::_mm_mul_ps(*m128WeightRow++, ::_mm_sub_ps(*m128InputRowBelow++, *m128InputRow++));
        }
        float* pfTempRow = (float*)m128TempRow;
        float* pfInputRow = (float*)m128InputRow;
        float* pfInputRowBelow = (float*)m128InputRowBelow;
        float* pfWeightRow = (float*)m128WeightRow;
        for (int i = 0; i < remainder; i++)
        {
            *pfTempRow++ = (*pfInputRowBelow++ - *pfInputRow++) * (*pfWeightRow++);
        }
    }
    //Process the first row
    __m128* m128TempRow = (__m128*)tempImage->Ptr();
    __m128* m128OutputRow = (__m128*)output.Ptr();
    for (int i = 0; i < sseLength; i++)
    {
        *m128OutputRow++ = ::_mm_sub_ps(*m128OutputRow, *m128TempRow++);
    }
    float* pfTempRow = (float*)m128TempRow;
    float* pfOutputRow = (float*)m128OutputRow;
    for (int i = 0; i < remainder; i++)
    {
        *pfOutputRow++ -= *pfTempRow++;
    }

    //Process the middle rows
    for (int iY = 1; iY < height - 1; iY++)
    {
        m128TempRow = (__m128*)tempImage->Ptr(iY);
        m128OutputRow = (__m128*)output.Ptr(iY);
        __m128* m128TempRowAbove = (__m128*)tempImage->Ptr(iY - 1);
        for (int i = 0; i <= sseLength; i++)
        {
            *m128OutputRow++ = ::_mm_sub_ps(::_mm_add_ps(*m128TempRowAbove++, *m128OutputRow), *m128TempRow++);
        }
        pfTempRow = (float*)m128TempRow;
        pfOutputRow = (float*)m128OutputRow;
        float* pfTempRowAbove = (float*)m128TempRowAbove;
        for (int i = 0; i < remainder; i++)
        {
            *pfOutputRow++ += (*pfTempRowAbove++ - *pfTempRow++);
        }
    }

    //Process the last row
    m128OutputRow = (__m128*)output.Ptr(height - 1);
    __m128* m128TempRowAbove = (__m128*)tempImage->Ptr(height - 2);
    for (int i = 0; i <= sseLength; i++)
    {
        *m128OutputRow++ = ::_mm_add_ps(*m128OutputRow, *m128TempRowAbove++);
    }
    pfOutputRow = (float*)m128OutputRow;
    float* pfTempRowAbove = (float*)m128TempRowAbove;
    for (int i = 0; i < remainder; i++)
    {
        *pfOutputRow++ += *pfTempRowAbove++;
    }
#else
    output,input,weight;
#endif
}

void CVOF2::ApplyFilters_LaplacianNoSSE(CFloatImg& output, CFloatImg& input, CFloatImg& weight)
{
    int width = input.Width();
    int height = input.Height();
    
    CFloatImg* tempImage = &m_TempImage;
    //Horizontal filtering
    for (int iY = 0; iY < height; iY++)
    {
        float* pfTempRow = tempImage->Ptr(iY);
        const float* pfInputRow = input.Ptr(iY);
        const float* pfWeightRow = weight.Ptr(iY);
        for (int iX = 0; iX < width - 1; iX++)
        {
            *pfTempRow++ = (*(pfInputRow + 1) - (*pfInputRow++)) * (*pfWeightRow++);
        }
    }
    for (int iY = 0; iY < height; iY++)
    {
        float* pfTempRow = tempImage->Ptr(iY);
        float* pfOutputRow = output.Ptr(iY);
        
        *pfOutputRow++ = -(*pfTempRow++); //Set first pixel
        for (int iX = 1; iX < width - 1; iX++) //Set middle pixels
        {
            *pfOutputRow++ = *(pfTempRow - 1) - (*pfTempRow++);
        }
        *pfOutputRow = *(pfTempRow - 1); //Set last pixel
    }

    //Vertical filtering
    for (int iY = 0; iY < height - 1; iY++)
    {
        float* pfTempRow = tempImage->Ptr(iY);
        const float* pfInputRow = input.Ptr(iY);
        const float* pfInputRowBelow = input.Ptr(iY + 1);
        const float* pfWeightRow = weight.Ptr(iY);
        for (int iX = 0; iX < width; iX++)
        {
            *pfTempRow++ = (*pfInputRowBelow++ - *pfInputRow++) * (*pfWeightRow++);
        }
    }
    //Process the first row
    float* pfTempRow = tempImage->Ptr();
    float* pfOutputRow = output.Ptr();
    for (int iX = 0; iX < width; iX++)
    {
        *pfOutputRow++ -= *pfTempRow++;
    }
    //Process the middle rows
    for (int iY = 1; iY < height - 1; iY++)
    {
        pfTempRow = tempImage->Ptr(iY);
        pfOutputRow = output.Ptr(iY);
        float* pfTempRowAbove = tempImage->Ptr(iY - 1);
        for (int iX = 0; iX < width; iX++)
        {
            *pfOutputRow++ += (*pfTempRowAbove++ - *pfTempRow++);
        }
    }
    //Process the last row
    pfOutputRow = output.Ptr(height - 1);
    float* pfTempRowAbove = tempImage->Ptr(height - 2);
    for (int iX = 0; iX < width; iX++)
    {
        *pfOutputRow++ += *pfTempRowAbove++;
    }
}

void CVOF2::ApplyFilter_Horizontal(CFloatImg& sourceImage, CFloatImg& destinationImage, double* filter, int filterSize, int iWidth /* = 0 */, int iHeight /* = 0 */)
{
    int width = (iWidth == 0) ? sourceImage.Width() : iWidth;
    int height = (iHeight == 0) ? sourceImage.Height() : iHeight;
#ifdef USE_FEATURES
    int nChannels = sourceImage.Bands();
    for (int iY = 0; iY < height; iY++)
    {
        float* pfSourceRow = sourceImage.Ptr(iY);
        float* pfDestinationRow = destinationImage.Ptr(iY);
        
        //First loop through the first few columns that when applied by filters are out of range
        for (int iX = 0; iX < filterSize; iX++)
        {
            for (int iC = 0; iC < nChannels; iC++)
            {
                *(pfDestinationRow + iX * nChannels + iC) = (float)(*(pfSourceRow + iC) * (*filter));
            }

            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                double weight = *(filter + iF + filterSize);
                int underPixelX = ::VtMax<int>(iX + iF, 0);
                for (int iC = 0; iC < nChannels; iC++)
                {
                    *(pfDestinationRow + iX * nChannels + iC) += (float)(*(pfSourceRow + underPixelX * nChannels + iC) * weight);
                }
            }
        }

        //Loop through the middle columns
        for (int iX = filterSize; iX < width - filterSize; iX++)
        {
            int underPixelX = iX - filterSize;
            for (int iC = 0; iC < nChannels; iC++)
            {
                *(pfDestinationRow + iX * nChannels + iC) = (float)(*(pfSourceRow + underPixelX * nChannels + iC) * (*filter));
            }

            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                double weight = *(filter + iF + filterSize);
                underPixelX = iX + iF;
                for (int iC = 0; iC < nChannels; iC++)
                {
                    *(pfDestinationRow + iX * nChannels + iC) += (float)(*(pfSourceRow + underPixelX * nChannels + iC) * weight);
                }
            }
        }

        //Loop through the last few columns that when applied by filters are out of range
        for (int iX = width - filterSize; iX < width; iX++)
        {
            int underPixelX = iX - filterSize;
            for (int iC = 0; iC < nChannels; iC++)
            {
                *(pfDestinationRow + iX * nChannels + iC) = (float)(*(pfSourceRow + underPixelX * nChannels + iC) * (*filter));
            }

            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                double weight = *(filter + iF + filterSize);
                underPixelX = ::VtMin<int>(iX + iF, width - 1);
                for (int iC = 0; iC < nChannels; iC++)
                {
                    *(pfDestinationRow + iX * nChannels + iC) += (float)(*(pfSourceRow + underPixelX * nChannels + iC) * weight);
                }
            }
        }
    }
#else
    double* weight = 0;
    float* pfSourcePixel = 0;
    for (int iY = 0; iY < height; iY++)
    {
        float* pfSourceRow = sourceImage.Ptr(iY);
        float* pfDestinationRow = destinationImage.Ptr(iY);
        
        //First loop through the first few columns that when applied by filters are out of range
        for (int iX = 0; iX < filterSize; iX++)
        {
            *pfDestinationRow = (*pfSourceRow) * (float)(*filter);
            weight = filter + 1;
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                int underPixelX = ::VtMax<int>(iX + iF, 0);
                *pfDestinationRow += *(pfSourceRow + underPixelX) * (float)(*weight++);
            }
            pfDestinationRow++;
        }

        //Loop through the middle columns
        for (int iX = filterSize; iX < width - filterSize; iX++)
        {
            pfSourcePixel = pfSourceRow + iX - filterSize;
            *pfDestinationRow = (*pfSourcePixel++) * (float)(*filter);
            weight = filter + 1;
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                *pfDestinationRow += (*pfSourcePixel++) * (float)(*weight++);
            }
            pfDestinationRow++;
        }

        //Loop through the last few columns that when applied by filters are out of range
        for (int iX = width - filterSize; iX < width; iX++)
        {
            int underPixelX = iX - filterSize;
            *pfDestinationRow = *(pfSourceRow + underPixelX) * (float)(*filter);
            weight = filter + 1;
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                underPixelX = ::VtMin<int>(iX + iF, width - 1);
                *pfDestinationRow += *(pfSourceRow + underPixelX) * (float)(*weight++);
            }
            pfDestinationRow++;
        }
    }
#endif
}

void CVOF2::ApplyFilter_Vertical(CFloatImg& sourceImage, CFloatImg& destinationImage, double* filter, int filterSize, int iWidth /* = 0 */, int iHeight /* = 0 */)
{
    int width = (iWidth == 0) ? sourceImage.Width() : iWidth;
    int height = (iHeight == 0) ? sourceImage.Height() : iHeight;
#ifdef USE_FEATURES
    int nChannels = sourceImage.Bands();
    for (int iY = 0; iY < filterSize; iY++) //Loop through first set of rows that are out of range
    {
        float* pfSourcePixel = sourceImage.Ptr(iY);
        float* pfDestinationPixel = destinationImage.Ptr(iY);
        for (int iX = 0; iX < width; iX++)
        {
            for (int iC = 0; iC < nChannels; iC++)
            {
                pfSourcePixel = sourceImage.Ptr() + iX * nChannels + iC;
                *(pfDestinationPixel + iX * nChannels + iC) = (float)((*pfSourcePixel) * (*filter));
            }
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                double weight = *(filter + iF + filterSize);
                int underPixelY = ::VtMax<int>(iY + iF, 0);
                for (int iC = 0; iC < nChannels; iC++)
                {
                    pfSourcePixel = sourceImage.Ptr(underPixelY) + iX * nChannels + iC;
                    *(pfDestinationPixel + iX * nChannels + iC) += (float)((*pfSourcePixel) * weight);
                }
            }
        }
    }

    for (int iY = filterSize; iY < height - filterSize; iY++) //Loop through middle rows
    {
        float* pfSourcePixel = sourceImage.Ptr(iY);
        float* pfDestinationPixel = destinationImage.Ptr(iY); 
        for (int iX = 0; iX < width; iX++)
        {
            int underPixelY = iY - filterSize;
            for (int iC = 0; iC < nChannels; iC++)
            {
                pfSourcePixel = sourceImage.Ptr(underPixelY) + iX * nChannels + iC;
                *(pfDestinationPixel + iX * nChannels + iC) = (float)((*pfSourcePixel) * (*filter));
            }
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                double weight = *(filter + iF + filterSize);
                underPixelY = iY + iF;
                for (int iC = 0; iC < nChannels; iC++)
                {
                    pfSourcePixel = sourceImage.Ptr(underPixelY) + iX * nChannels + iC;
                    *(pfDestinationPixel + iX * nChannels + iC) += (float)((*pfSourcePixel) * weight);
                }
            }
        }
    }

    for (int iY = height - filterSize; iY < height; iY++) //Loop through last rows
    {
        float* pfSourcePixel = sourceImage.Ptr(iY);
        float* pfDestinationPixel = destinationImage.Ptr(iY); 
        for (int iX = 0; iX < width; iX++)
        {
            int underPixelY = iY - filterSize;
            for (int iC = 0; iC < nChannels; iC++)
            {
                pfSourcePixel = sourceImage.Ptr(underPixelY) + iX * nChannels + iC;
                *(pfDestinationPixel + iX * nChannels + iC) = (float)((*pfSourcePixel) * (*filter));
            }
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                double weight = *(filter + iF + filterSize);
                underPixelY = ::VtMin<int>(iY + iF, height - 1);
                for (int iC = 0; iC < nChannels; iC++)
                {
                    pfSourcePixel = sourceImage.Ptr(underPixelY) + iX * nChannels + iC;
                    *(pfDestinationPixel + iX * nChannels + iC) += (float)((*pfSourcePixel) * weight);
                }
            }
        }
    }
#else
    double* weight = 0;
    for (int iY = 0; iY < filterSize; iY++) //Loop through first set of rows that are out of range
    {
        float* pfSourcePixel = sourceImage.Ptr(iY);
        float* pfDestinationPixel = destinationImage.Ptr(iY);
        for (int iX = 0; iX < width; iX++)
        {
            weight = filter + 1;
            pfSourcePixel = sourceImage.Ptr() + iX;
            *pfDestinationPixel = (*pfSourcePixel) * (float)(*filter);
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                int underPixelY = ::VtMax<int>(iY + iF, 0);
                pfSourcePixel = sourceImage.Ptr(underPixelY) + iX;
                *pfDestinationPixel += (*pfSourcePixel) * (float)(*weight++);
            }
            pfDestinationPixel++;
        }
    }

    for (int iY = filterSize; iY < height - filterSize; iY++) //Loop through middle rows
    {
        float* pfSourcePixel = sourceImage.Ptr(iY);
        float* pfDestinationPixel = destinationImage.Ptr(iY); 
        for (int iX = 0; iX < width; iX++)
        {
            weight = filter + 1;
            int underPixelY = iY - filterSize;
            pfSourcePixel = sourceImage.Ptr(underPixelY) + iX;
            *pfDestinationPixel = (*pfSourcePixel) * (float)(*filter);
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                underPixelY = iY + iF;
                pfSourcePixel = sourceImage.Ptr(underPixelY) + iX;
                *pfDestinationPixel += (*pfSourcePixel) * (float)(*weight++);
            }
            pfDestinationPixel++;
        }
    }

    for (int iY = height - filterSize; iY < height; iY++) //Loop through last rows
    {
        float* pfSourcePixel = sourceImage.Ptr(iY);
        float* pfDestinationPixel = destinationImage.Ptr(iY); 
        for (int iX = 0; iX < width; iX++)
        {
            weight = filter + 1;
            int underPixelY = iY - filterSize;
            pfSourcePixel = sourceImage.Ptr(underPixelY) + iX;
            *pfDestinationPixel = (*pfSourcePixel) * (float)(*filter);
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                underPixelY = ::VtMin<int>(iY + iF, height - 1);
                pfSourcePixel = sourceImage.Ptr(underPixelY) + iX;
                *pfDestinationPixel += (*pfSourcePixel) * (float)(*weight++);
            }
            pfDestinationPixel++;
        }
    }
#endif
}

void CVOF2::ApplyFilter_Horizontal(Byte* pbSourceImage, size_t strideBytes, Byte* pbDestinationImage, 
    int width, int height, int nChannels, double* filter, int filterSize)
{
    for (int iY = 0; iY < height; iY++)
    {
        int offset = iY * (int)strideBytes; //Number of bytes offset to the beginning of row iY
        float* pfSourceRow = (float*)(pbSourceImage + offset); //Index to the pixel
        float* pfDestinationRow = (float*)(pbDestinationImage + offset); 
        
        for (int iX = 0; iX < width; iX++)
        {
            //Find the underlying pixel x-coordinate for this location and if it's out of bound then
            //reset it to the boundary, either 0 or width - 1.
            int underPixelX = ::VtMin<int>(::VtMax<int>(iX - filterSize, 0), width - 1);
            for (int iC = 0; iC < nChannels; iC++)
            {
                *(pfDestinationRow + iX * nChannels + iC) = (float)(*(pfSourceRow + underPixelX * nChannels + iC) * (*filter));
            }

            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                double weight = *(filter + iF + filterSize);
                
                //Find the underlying pixel x-coordinate for this location and if it's out of bound then
                //reset it to the boundary, either 0 or width - 1.
                underPixelX = ::VtMin<int>(::VtMax<int>(iX + iF, 0), width - 1);
                for (int iC = 0; iC < nChannels; iC++)
                {
                    *(pfDestinationRow + iX * nChannels + iC) += (float)(*(pfSourceRow + underPixelX * nChannels + iC) * weight);
                }
            }
        }
    }
}

void CVOF2::ApplyFilter_Vertical(Byte* pbSourceImage, size_t strideBytes, Byte* pbDestinationImage, 
    int width, int height, int nChannels, double* filter, int filterSize)
{
    float* pfSourcePixel = (float*)pbSourceImage;
    for (int iY = 0; iY < height; iY++)
    {
        //Index to the beginning of row iY
        float* pfDestinationPixel = (float*)(pbDestinationImage + iY * strideBytes); 

        for (int iX = 0; iX < width; iX++)
        {
            //Find the underlying pixel y-coordinate for this location and if it's out of bound then
            //reset it to the boundary, either 0 or height - 1.
            int underPixelY = ::VtMin<int>(::VtMax<int>(iY - filterSize, 0), height - 1);

            for (int iC = 0; iC < nChannels; iC++)
            {
                pfSourcePixel = (float*)(pbSourceImage + underPixelY * strideBytes) + iX * nChannels + iC;
                *(pfDestinationPixel + iX * nChannels + iC) = (float)((*pfSourcePixel) * (*filter));
            }
            for (int iF = -filterSize + 1; iF <= filterSize; iF++)
            {
                double weight = *(filter + iF + filterSize);
                //Find the underlying pixel y-coordinate for this location and if it's out of bound then
                //reset it to the boundary, either 0 or height - 1.
                underPixelY = ::VtMin<int>(::VtMax<int>(iY + iF, 0), height - 1);

                for (int iC = 0; iC < nChannels; iC++)
                {
                    pfSourcePixel = (float*)(pbSourceImage + underPixelY * strideBytes) + iX * nChannels + iC;
                    *(pfDestinationPixel + iX * nChannels + iC) += (float)((*pfSourcePixel) * weight);
                }
            }
        }
    }
}

void CVOF2::EstimateLaplacianNoise(const CFloatImg& image1, const CFloatImg& image2, double*& parameter)
{
    int height = image1.Height();
    int width = image1.Width();

    int paramSize = image1.Bands() + 2;
    for (int i = 0; i < paramSize; i++)
    {
        parameter[i] = 0.0;
    }
    double* total = new double[paramSize];

#ifdef USE_FEATURES
    int nChannels = image1.Bands();
    for (int k = 0; k < nChannels; k++)
    {
        total[k] = 0;
    }
    double temp;
    for (int iY = 0; iY < height; iY++)
    {
        const float* pfImageRow1 = image1.Ptr(iY);
        const float* pfImageRow2 = image2.Ptr(iY);

        for (int iX = 0; iX < width; iX++)
        {
            int offset = iX * nChannels;
            for (int iC = 0; iC < nChannels; iC++)
            {
                temp = ::abs((*(pfImageRow1 + offset + iC)) - (*(pfImageRow2 + offset + iC)));
                if (temp > 0 && temp < 1000000)
                {
                    parameter[iC] += temp;
                    total[iC]++;
                }
            }
        }
    }

    for (int iC = 0; iC < nChannels; iC++)
    {
        if (total[iC] == 0)
        {
            parameter[iC] = 0.001;
        }
        else
        {
            parameter[iC] /= total[iC];
        }
    }
#else
    total[0] = 0.0;
    double temp;
    for (int iY = 0; iY < height; iY++)
    {
        const float* pfImageRow1 = image1.Ptr(iY);
        const float* pfImageRow2 = image2.Ptr(iY);
        for (int iX = 0; iX < width; iX++)
        {
            temp = ::abs((*pfImageRow1++) - (*pfImageRow2++));
            if (temp > 0 && temp < 1000000)
            {
                parameter[0] += temp;
                total[0]++;
            }
        }
    }
    if (total[0] == 0)
    {
        parameter[0] = 0.001;
    }
    else
    {
        parameter[0] /= total[0];
    }
#endif
}

void CVOF2::MultiplyImage(CFloatImg& image1, CFloatImg& image2, CFloatImg& image3, CFloatImg& output)
{
    int width = image1.Width();
    int height = image1.Height();
    int nChannels = image1.Bands();
    int nRowElements = width * nChannels;
    
#if (defined(_M_IX86) || defined(_M_AMD64))
    if (m_bUseSSE && m_bSSE2)
    {
        int remainder = nRowElements % 4;
        int sseLength = nRowElements / 4;
        for (int iY = 0; iY < height; iY++)
        {
            __m128* m128Image1Row = (__m128*)image1.Ptr(iY);
            __m128* m128Image2Row = (__m128*)image2.Ptr(iY);
            __m128* m128Image3Row = (__m128*)image3.Ptr(iY);
            __m128* m128OutputRow = (__m128*)output.Ptr(iY);

            for (int i = 0; i < sseLength; i++)
            {
                *m128OutputRow++ = ::_mm_mul_ps(::_mm_mul_ps(*m128Image1Row++, *m128Image2Row++), *m128Image3Row++);
            }
            for (int i = 0; i < remainder; i++)
            {
                float* pfImage1Row = (float*)m128Image1Row + i;
                float* pfImage2Row = (float*)m128Image2Row + i;
                float* pfImage3Row = (float*)m128Image3Row + i;
                float* pfOutputRow = (float*)m128OutputRow + i;
                *pfOutputRow = (*pfImage1Row) * (*pfImage2Row) * (*pfImage3Row);
            }
        }
    }
    else
#endif
    {
        for (int iY = 0; iY < height; iY++)
        {
            float* pfImage1Row = image1.Ptr(iY);
            float* pfImage2Row = image2.Ptr(iY);
            float* pfImage3Row = image3.Ptr(iY);
            float* pfOutputRow = output.Ptr(iY);
            float* pfImage1RowEnd = pfImage1Row + nRowElements;

            while (pfImage1Row < pfImage1RowEnd)
            {
                *pfOutputRow++ = (*pfImage1Row++) * (*pfImage2Row++) * (*pfImage3Row++);
            }
        }
    }
}

void CVOF2::CollapseImage(CFloatImg& sourceImage, CFloatImg& destinationImage)
{
    int height = destinationImage.Height();
    int width = destinationImage.Width();
    int nChannels = sourceImage.Bands();

    for (int iY = 0; iY < height; iY++)
    {
        float* pfSourceRow = sourceImage.Ptr(iY);
        float* pfDestinationRow = destinationImage.Ptr(iY);

        for (int iX = 0; iX < width; iX++)
        {
            float averagePixelValue = 0;
            for (int iC = 0; iC < nChannels; iC++)
            {
                averagePixelValue += *(pfSourceRow + iC);
            }

            *pfDestinationRow++ = averagePixelValue / nChannels;
            pfSourceRow += nChannels;
        }
    }
}

HRESULT CVOF2::ExpandFlow(int iLevel)
{
    VT_HR_BEGIN()

    float scalar = 2.0f;

    CFloatImg* pFlowImageXUp = &(m_pyrFlowX.GetLevel(iLevel + 1));
    CFloatImg* pCurrentFlowImageX = &(m_pyrFlowX.GetLevel(iLevel));
    pCurrentFlowImageX->Clear(); //LH:Todo:Consider removing this clear call
    CRect dstXRect(0, 0, pCurrentFlowImageX->Width(), pCurrentFlowImageX->Height());
    /*VT*/VT_HR_EXIT(::VtResizeImage(*pCurrentFlowImageX, dstXRect, *pFlowImageXUp, eSamplerKernelBilinear));
    /*VT*/VT_HR_EXIT(::VtScaleImage(*pCurrentFlowImageX, *pCurrentFlowImageX, scalar));

    CFloatImg* pFlowImageYUp = &(m_pyrFlowY.GetLevel(iLevel + 1));
    CFloatImg* pCurrentFlowImageY = &(m_pyrFlowY.GetLevel(iLevel));
    pCurrentFlowImageY->Clear(); //LH:Todo:Consider removing this clear call
    CRect dstYRect(0, 0, pCurrentFlowImageY->Width(), pCurrentFlowImageY->Height());
    /*VT*/VT_HR_EXIT(::VtResizeImage(*pCurrentFlowImageY, dstYRect, *pFlowImageYUp, eSamplerKernelBilinear));
    /*VT*/VT_HR_EXIT(::VtScaleImage(*pCurrentFlowImageY, *pCurrentFlowImageY, scalar));
    
    VT_HR_END()
}

#ifdef USE_FEATURES
HRESULT CVOF2::ImageToFeature(int imageIndex, int iLevel)
{
    VT_HR_BEGIN()
    
    CFloatImg* pInputImage			= &(m_ppyrInputImages[imageIndex].GetLevel(iLevel));
    CFloatImg* pFeatureImage		= &(m_ppyrFeatureImages[imageIndex].GetLevel(iLevel));
    CFloatImg* pDerivativeXImage	= &m_SmoothFlowTempImage1;
    CFloatImg* pDerivativeYImage	= &m_SmoothFlowTempImage2;

    int width = pInputImage->Width();
    int height = pInputImage->Height();
    int numberOfChannels = pInputImage->Bands();

    switch (numberOfChannels)
    {
        case 1:
        {
            VT_HR_EXIT(this->ComputeDerivativeX(*pInputImage, *pDerivativeXImage, m_bSmoothDerivatives, width, height));
            VT_HR_EXIT(this->ComputeDerivativeY(*pInputImage, *pDerivativeYImage, m_bSmoothDerivatives, width, height));

            for (int iY = 0; iY < height; iY++)
            {
                //Get a pointer to the beginning of feature image's row iY
                float* fpFeatureImageRow = pFeatureImage->Ptr(iY);
                //Since each pixel in the original image corresponds to a group of "n" pixels
                //in the feature image, the length of one feature row is "n" times more.
                float* fpFeatureImageRowEnd = fpFeatureImageRow + pFeatureImage->Bands() * width;

                //Get a pointer to the beginning of pyramid image's row iY
                float* fpPyramidImageRow = pInputImage->Ptr(iY);
                
                //Get pointers to the derivate rows
                float* fpXDerivativeRow  = pDerivativeXImage->Ptr(iY);
                float* fpYDerivativeRow  = pDerivativeYImage->Ptr(iY);
                
                //Embed the image pixel value and its x and y-derivative.
                while (fpFeatureImageRow < fpFeatureImageRowEnd)
                {
                    *fpFeatureImageRow++ = *fpPyramidImageRow++;
                    *fpFeatureImageRow++ = *fpXDerivativeRow++;
                    *fpFeatureImageRow++ = *fpYDerivativeRow++;
                }
            }
            break;
        }
        case 3:
        {
            //Only use the DFlow image here for temporary calculation. This shouldn't affect
            //the true values of DFlow later done in S.O.R. iterations.
            CFloatImg* pTempGrayImage = &(m_pyrDFlowX.GetLevel(iLevel));
            for (int iY = 0; iY < height; iY++) //Desaturate the original image
            {
                float* pfInputRow = pInputImage->Ptr(iY);
                float* pfInputRowEnd = pfInputRow + width * numberOfChannels;
                float* pfGrayRow = pTempGrayImage->Ptr(iY);
                while (pfInputRow < pfInputRowEnd)
                {
                    /*VT*/*pfGrayRow++ = ::VtLumaFromRGB_CCIR601YPbPr(*(pfInputRow + 2), *(pfInputRow + 1), *pfInputRow);
                    pfInputRow += 3;
                }
            }
            //Now compute the derivatives of the gray image
            VT_HR_EXIT(this->ComputeDerivativeX(*pTempGrayImage, *pDerivativeXImage, m_bSmoothDerivatives, width, height));
            VT_HR_EXIT(this->ComputeDerivativeY(*pTempGrayImage, *pDerivativeYImage, m_bSmoothDerivatives, width, height));
            
            //Now create the features
            for (int iY = 0; iY < height; iY++)
            {
                float* pfInputRow = pInputImage->Ptr(iY);
                float* pfInputRowEnd = pfInputRow + width * numberOfChannels;
                float* pfGrayRow = pTempGrayImage->Ptr(iY);
                float* pfFeatureImageRow = pFeatureImage->Ptr(iY);
                float* pfDerivativeXRow = pDerivativeXImage->Ptr(iY);
                float* pfDerivativeYRow = pDerivativeYImage->Ptr(iY);
                while (pfInputRow < pfInputRowEnd)
                {
                    *pfFeatureImageRow++ = *pfGrayRow++;
                    *pfFeatureImageRow++ = *pfDerivativeXRow++;
                    *pfFeatureImageRow++ = *pfDerivativeYRow++;
                    *pfFeatureImageRow++ = *(pfInputRow + 1) - *pfInputRow;
                    *pfFeatureImageRow++ = *(pfInputRow + 1) - *(pfInputRow + 2);
                    pfInputRow += 3;
                }
            }
            break;
        }
        default:
        {
            break;
        }
    }
    VT_HR_EXIT_LABEL()
    return hr;
}
#endif

HRESULT CVOF2::GetFlowX(CFloatImg* &pimgFlowX, int iLevel)
{
    VT_HR_BEGIN()
    
    // Check Inputs OK
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iLevel>=m_iLevels || iLevel<0) ? E_INVALIDARG : S_OK);

    // Return the Result
    pimgFlowX = &m_pyrFlowX.GetLevel(iLevel);

    VT_HR_EXIT_LABEL()
    return hr;
}

HRESULT CVOF2::GetFlowY(CFloatImg* &pimgFlowY, int iLevel)
{
    VT_HR_BEGIN()
    
    // Check Inputs OK
    VT_HR_EXIT((!Allocated()) ? E_FAIL : S_OK);
    VT_HR_EXIT((iLevel>=m_iLevels || iLevel<0) ? E_INVALIDARG : S_OK);

    // Return the Result
    pimgFlowY = &m_pyrFlowY.GetLevel(iLevel);

    VT_HR_EXIT_LABEL()
    return hr;
}