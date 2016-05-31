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

/// <summary>
/// When features are turned on, a denser feature image is created and is used to compute flow
/// instead of the original input image. This feature image contains both the original RGB values as well
/// as the input's derivatives information. Turn this on to get higher quality flow result but sacrifice execution speed.
/// </summary>
//#define USE_FEATURES

#pragma once

#include "iopticalflow.h"
#include "vtcore.h"
#include "vtfileio.h"

namespace vt
{
/// <summary> The internal parameters of a variational optical flow algorithm. Not all these
/// parameters are exposed in the IOpticalFlow interface. </summary>
struct CVOF2Params
{
    int		m_iTopLevelMinPixelCount; //Specify the minimum pixel count that the input image can be downsampled to
    bool	m_bNormalizePixelValues; //Normalize pixel values to [0,1]
    int		m_iOuterIterations; //Outer number of iterations
    int		m_iInnerIterations; //Inner number of iterations
    int		m_iSORIterations; //S.O.R number of iterations
    float	m_Omega; //Omega used in S.O.R calculations
    float	m_Alpha; //Alpha used in S.O.R calculations
    float	m_Rho; //Rho value to be used in 1-Rho-1 filter, only applies if using smooth pyramid
    int		m_iImageWarpMethod; //Method to warp image, see enum ImageWarpMethod definition
    bool	m_bGaussianGradients; //Use smooth gaussian gradients Dx, Dy, Dt for S.O.R
    bool	m_bSmoothDerivatives; //Use universal smooth derivatives
    int		m_iNoiseModel; //Noise model, see enum NoiseModel definition
    bool	m_bUseSSE; //Use SSE computations if available, otherwise disable SSE completely

    CVOF2Params(int     iTopLevelMinPixelCount,
                bool    bNormalizePixelValues,
                int     iOuterIterations,
                int     iInnerIterations,
                int     iSORIterations,
                float   Omega,
                float   Alpha,
                float   Rho,
                int     iImageWarpMethod,
                bool    bGaussianGradients,
                bool    bSmoothDerivatives,
                int     iNoiseModel,
                bool    bUseSSE) :
        m_iTopLevelMinPixelCount(iTopLevelMinPixelCount),
        m_bNormalizePixelValues(bNormalizePixelValues),
        m_iOuterIterations(iOuterIterations),
        m_iInnerIterations(iInnerIterations),
        m_iSORIterations(iSORIterations),
        m_Omega(Omega),
        m_Alpha(Alpha),
        m_Rho(Rho),
        m_iImageWarpMethod(iImageWarpMethod),
        m_bGaussianGradients(bGaussianGradients),
        m_bSmoothDerivatives(bSmoothDerivatives),
        m_iNoiseModel(iNoiseModel),
        m_bUseSSE(bUseSSE){}
};

/// <summary>
/// Parameters class to use for CVOF2 flow computation.
/// </summary>
extern CVOF2Params OPTICAL_FLOW_V2_PARAMS;

class CVOF2 : public IOpticalFlow
{
public:

    /// <summary> Constructor. </summary>
    CVOF2() { Init(); }

protected:

    /// <summary> Initializes member variables. </summary>
    void Init();

public:

    /// <summary> Destructor. </summary>
    ~CVOF2() { Deallocate(); }

    /// <summary>
    /// Allocates memory needed for images, pyramids, etc...
    /// </summary>
    /// <param name = "iInputWidth"> The width of the input image. </param>
    /// <param name = "iInputHeight"> The height of the input image. </param>
    /// <param name = "params"> A pointer to a CVOF2Params object that contains parameter settings. </param>
    /// <param name = "iSubsample"> Specify the subsample level. </param>
    /// <returns>HRESULT</returns>
    HRESULT Allocate(int iInputWidth, int iInputHeight, CVOF2Params* params, int iSubsample = 1);
    
    /// <summary>
    /// Deallocates any memory previously reserved.
    /// </summary>
    void Deallocate();

    /// <summary>
    /// Returns whether appropriate memory has been allocated in order to compute optical flow.
    /// </summary>
    /// <returns>bool</returns>
    bool Allocated() { return m_bAllocated; }
    
    /// <summary>
    /// Reads and sets up appropriate input parameters to the flow computation.
    /// </summary>
    /// <param name = "params"> A pointer to a CVOF2Params object that contains parameter settings. </param>
    /// <returns>HRESULT</returns>
    HRESULT ReadInput(CVOF2Params* params);

    /// <summary> Enumeration of different noise models used for estimation while computing flow. </summary>
    enum NoiseModel { GaussianMixture = 1, Laplacian };

    /// <summary> Enumeration of different warping method. </summary>
    enum ImageWarpMethod { Bilinear = 1, Bicubic };

    /// <summary> Adds a 32Bit CRGBAImg to the optical flow algorithm. </summary>
    virtual HRESULT AddImage(CRGBAImg &imgInput)
    {
        return AddBGRAImage(imgInput);
    }
    
    /// <summary> Adds an 8bit CLumaImg to the optical flow algorithm. </summary>
    /// <remarks> Not currently supported. </remarks>
    virtual HRESULT AddImage(CLumaImg &imgInput)
    {
        UNREFERENCED_PARAMETER(imgInput);
        return E_NOTIMPL;
    }

public:
    
    /// <summary> Gets a pointer to a CFloatImg containing the X component of the flow. </summary>
    /// <param name="pimgFlowX"> A pointer that receives the flow. </param>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
    /// <returns>HRESULT</returns>
    virtual HRESULT GetFlowX(CFloatImg* &pimgFlowX, int iLevel = 0);
    
    /// <summary> Gets a pointer to a CFloatImg containing the Y component of the flow. </summary>
    /// <param name="pimgFlowY"> A pointer that receives the flow. </param>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
    /// <returns>HRESULT</returns>
    virtual HRESULT GetFlowY(CFloatImg* &pimgFlowY, int iLevel = 0);

private:

    /// <summary>
    /// Adds an RGBA image as input.
    /// </summary>
    /// <param name = "imgInput"> An RGBA image to add as input to the optical flow algorithm. </param>
    /// <returns>HRESULT</returns>
    HRESULT AddBGRAImage(CRGBAImg &imgInput);

    /// <summary>
    /// Compute flow using Ce's algorithm.
    /// </summary>
    /// <param name = "iLowestLevel"> The lowest level of the image pyramid. </param>
    /// <returns>HRESULT</returns>
    virtual HRESULT ComputeFlow(int iLowestLevel = 0);

    /// <summary> Returns the number of levels in the pyramid. </summary>
    virtual int GetNumberOfLevels() { return m_iLevels; }

    /// <summary>
    /// Computes the x-derivative image for the specified image.
    /// </summary>
    /// <param name = "sourceImage"> The source image to compute derivatives for. </param>
    /// <param name = "derivativeXImage"> The output derivative image. </param>
    /// <param name = "useSmoothDerivative"> Whether to compute smooth derivatives. </param>
    /// <param name = "derivativeWidth"> Optional; the width of the derivative region to compute within. </param>
    /// <param name = "derivativeHeight"> Optional; the height of the derivative region to compute within. </param>
    /// <remarks>
    /// The number of color channels of both images should match.
    /// </remarks>
    HRESULT ComputeDerivativeX(CFloatImg& sourceImage, CFloatImg& derivativeXImage, bool useSmoothDerivative, int derivativeWidth = 0, int derivativeHeight = 0);

    /// <summary>
    /// Computes the y-derivative image for the specified image.
    /// </summary>
    /// <param name = "sourceImage"> The source image to compute derivatives for. </param>
    /// <param name = "derivativeYImage"> The output derivative image. </param>
    /// <param name = "useSmoothDerivative"> Whether to compute smooth derivatives. </param>
    /// <param name = "derivativeWidth"> Optional; the width of the derivative region to compute within. </param>
    /// <param name = "derivativeHeight"> Optional; the height of the derivative region to compute within. </param>
    /// <remarks>
    /// The number of color channels of both images should match.
    /// </remarks>
    HRESULT ComputeDerivativeY(CFloatImg& sourceImage, CFloatImg& derivativeYImage, bool useSmoothDerivative, int derivativeWidth = 0, int derivativeHeight = 0);
    
    /// <summary>
    /// Computes the phi component of the linear system.
    /// </summary>
    /// <param name = "derivativeXFlowX"> The X derivative of the X-Flow image. </param>
    /// <param name = "derivativeYFlowX"> The Y derivative of the X-Flow image. </param>
    /// <param name = "derivativeXFlowY"> The X derivative of the Y-Flow image. </param>
    /// <param name = "derivativeYFlowY"> The Y derivative of the Y-Flow image. </param>
    /// <param name = "phi"> The phi to be computed. </param>
    /// <param name = "epsilon"> The constant to pad in the computation. </param>
    HRESULT ComputePhi(CFloatImg& derivativeXFlowX, CFloatImg& derivativeYFlowX, CFloatImg& derivativeXFlowY, CFloatImg& derivativeYFlowY, CFloatImg& phi, double& epsilon);
    
    /// <summary>
    /// Computes the psi component of the linear system.
    /// </summary>
    /// <param name = "derivativeX"> The X derivative of the feature image. </param>
    /// <param name = "derivativeY"> The Y derivative of the feature image. </param>
    /// <param name = "errorImage"> The error image. </param>
    /// <param name = "dflowX"> The differential change image of X-Flow. </param>
    /// <param name = "dflowY"> The differential change image of Y-Flow. </param>
    /// <param name = "psi"> The psi to be computed. </param>
    /// <param name = "parameters"> The noise parameters. </param>
    /// <param name = "epsilon"> The constant to pad in the computation. </param>
    HRESULT ComputePsi(CFloatImg& derivativeX, CFloatImg& derivativeY, CFloatImg& errorImage, CFloatImg& dflowX, CFloatImg& dflowY, CFloatImg& psi, double*& parameters, double& epsilon);

    /// <summary>
    /// Warps image based on the interpolation method.
    /// </summary>
    /// <param name = "iLevel"> The current level of the pyramid to work with. </param>
    void WarpImage(int iLevel = 0);

    /// <summary>
    /// Warps image using the bilinear method.
    /// </summary>
    /// <param name = "iLevel"> The current level of the pyramid to work with. </param>
    void WarpImage_Bilinear(int iLevel = 0);

    /// <summary>
    /// Interpolate the specified image into the result image.
    /// </summary>
    /// <param name = "image"> The image to be interpolated from. </param>
    /// <param name = "width"> The width of the image. </param>
    /// <param name = "height"> The height of the image. </param>
    /// <param name = "strideBytes"> The offset between rows as number of bytes. </param>
    /// <param name = "channels"> Number of color channels of the image. </param>
    /// <param name = "x"> The x-location of the calculated warped pixel. </param>
    /// <param name = "y"> The y-location of the calculated warped pixel. </param>
    /// <param name = "resultImage"> The resulting image. </param>
    /// <remarks>
    /// Inline method. Examines the location of the calculated warped pixel (need not be 
    /// integer value pairs) and determines the 4 pixels that form a 1-1 square which contains it.
    /// These 4 pixel values are interpolated to form the final value of the warped pixel.
    /// </remarks>
    static inline void BilinearInterpolate(float* image, int width, int height, size_t strideBytes, int channels, double x, double y, float* resultImage);
    
    /// <summary>
    /// Generates and solves the linear system.
    /// </summary>
    /// <param name = "alpha"> </param>
    /// <param name = "nOuterFPIterations"> Number of iterations to determine new flow through differential steps. </param>
    /// <param name = "nInnerFPIterations"> Number of inner iterations to calculate each differential flow step.</param>
    /// <param name = "nSORIterations"> Number of S.O.R iterations to solve the linear system. </param>
    void SmoothFlowSOR(int iLevel, double alpha, int nOuterFPIterations, int nInnerFPIterations, int nSORIterations);
    
    /// <summary>
    /// Solves the linear system using S.O.R method.
    /// </summary>
    /// <param name = "iLevel"> The current level on the pyramid. </param>
    /// <param name = "pimdxy"> Part of the linear system: Ix * Iy. </param>
    /// <param name = "pimdx2"> Part of the linear system: Ix * Ix. </param>
    /// <param name = "pimdy2"> Part of the linear system: Iy * Iy. </param>
    /// <param name = "pimdtdx"> Part of the linear system: Iz * Ix. </param>
    /// <param name = "pimdtdy"> Part of the linear system: Iz * Iy. </param>
    /// <param name = "alpha">Constant that defines how much neighboring pixels should move together. </param>
    /// <param name = "sorAlpha"> Precomputed scaled alpha constant to ease calculations.</param>
    /// <param name = "nSORIterations"> Number of S.O.R iterations to solve the linear system. </param>
    void SolveLinearSystems_SOR(int iLevel, CFloatImg& pimdxy, CFloatImg& pimdx2, CFloatImg& pimdy2, CFloatImg& pimdtdx, CFloatImg& pimdtdy, double alpha, double sorAlpha, int nSORIterations);
    
    //See documentation from CVOF2::SolveLinearSystems_SOR
    void SolveLinearSystems_SOR_NoSSE(int iLevel, CFloatImg& pimdxy, CFloatImg& pimdx2, CFloatImg& pimdy2, CFloatImg& pimdtdx, CFloatImg& pimdtdy, double alpha, double sorAlpha, int nSORIterations);
    
    //See documentation from CVOF2::SolveLinearSystems_SOR
    void SolveLinearSystems_SOR_SSE(int iLevel, CFloatImg& pimdxy, CFloatImg& pimdx2, CFloatImg& pimdy2, CFloatImg& pimdtdx, CFloatImg& pimdtdy, double alpha, double sorAlpha, int nSORIterations);
    
    /// <summary>
    /// Computes gradients for the combined image between the input image and warped image.
    /// </summary>
    /// <param name = "iLevel"> The current level of the pyramid. </param>
    void ComputeGradients(int iLevel = 0);

    /// <summary>
    /// Applies the laplacian filter combined with the <see cref="weight" /> values to the <see cref="input"/> image and output to the <see cref="output"/> image.
    /// </summary>
    /// <param name = "iLevel"> The level of the pyramid. </param>
    /// <param name = "output"> The output image. </param>
    /// <param name = "input"> The input image. </param>
    /// <param name = "weight"> The weight image. </param>
    void ApplyFilters_Laplacian(CFloatImg& output, CFloatImg& input, CFloatImg& weight);
    
    //See documentation from CVOF2::ApplyFilters_Laplacian
    void ApplyFilters_LaplacianNoSSE(CFloatImg& output, CFloatImg& input, CFloatImg& weight);

    //See documentation from CVOF2::ApplyFilters_Laplacian
    void ApplyFilters_LaplacianSSE(CFloatImg& output, CFloatImg& input, CFloatImg& weight);
    
    /// <summary>
    /// Applies horizontal filter to the source image and output to the destination image.
    /// </summary>
    /// <param name = "pbSourceImage"> The source image byte pointer. </param>
    /// <param name = "pbDestinationImage"> The destination image byte pointer. </param>
    /// <param name = "width"> The width of the image. </param>
    /// <param name = "height"> The height of the image. </param>
    /// <param name = "nChannels"> Number of color channels. </param>
    /// <param name = "filter"> The filter to apply. </param>
    /// <param name = "filterSize"> The size of the filter. </param>
    void ApplyFilter_Horizontal(Byte* pbSourceImage, size_t strideBytes, Byte* pbDestinationImage, 
        int width, int height, int nChannels, double* filter, int filterSize);
    
    /// <summary>
    /// Applies vertical filter to the source image and output to the destination image.
    /// </summary>
    /// <param name = "pbSourceImage"> The source image byte pointer. </param>
    /// <param name = "pbDestinationImage"> The destination image byte pointer. </param>
    /// <param name = "width"> The width of the image. </param>
    /// <param name = "height"> The height of the image. </param>
    /// <param name = "nChannels"> Number of color channels. </param>
    /// <param name = "filter"> The filter to apply. </param>
    /// <param name = "filterSize"> The size of the filter. </param>
    void ApplyFilter_Vertical(Byte* pbSourceImage, size_t strideBytes, Byte* pbDestinationImage, 
        int width, int height, int nChannels, double* filter, int filterSize);
    
    /// <summary>
    /// Applies horizontal filter to the source image and output to the destination image.
    /// </summary>
    /// <param name = "sourceImage"> The source image. </param>
    /// <param name = "destinationImage"> The destination image. </param>
    /// <param name = "filter"> The filter to apply. </param>
    /// <param name = "filterSize"> The size of the filter. </param>
    /// <param name = "iWidth"> Optional; The width of the region to apply filters to. </param>
    /// <param name = "iHeight"> Optional; The height of the region to apply filters to. </param>
    void ApplyFilter_Horizontal(CFloatImg& sourceImage, CFloatImg& destinationImage, double* filter, int filterSize, int width = 0, int height = 0);
    
    /// <summary>
    /// Applies vertical filter to the source image and output to the destination image.
    /// </summary>
    /// <param name = "sourceImage"> The source image. </param>
    /// <param name = "destinationImage"> The destination image. </param>
    /// <param name = "filter"> The filter to apply. </param>
    /// <param name = "filterSize"> The size of the filter. </param>
    /// <param name = "iWidth"> Optional; The width of the region to apply filters to. </param>
    /// <param name = "iHeight"> Optional; The height of the region to apply filters to. </param>
    void ApplyFilter_Vertical(CFloatImg& sourceImage, CFloatImg& destinationImage, double* filter, int filterSize, int width = 0, int height = 0);

    /// <summary>
    /// Multiplies images together and output to the specified image.
    /// </summary>
    /// <param name = "image1"> The first image. </param>
    /// <param name = "image2"> The second image. </param>
    /// <param name = "image3"> The third image. </param>
    /// <param name = "output"> The resulting image. </param>
    void MultiplyImage(CFloatImg& image1, CFloatImg& image2, CFloatImg& image3, CFloatImg& output);

    /// <summary>
    /// Collapses the source image and output to the destination image.
    /// </summary>
    /// <param name = "sourceImage"> The source image. </param>
    /// <param name = "destinationImage"> The destination image which should be in grayscale format. </param>
    /// <remarks>
    /// If the images are 1-band then this is equivalent to a copy.
    /// </remarks>
    void CollapseImage(CFloatImg& sourceImage, CFloatImg& destinationImage);
    
    /// <summary>
    /// Estimates noise using the Laplacian model and the two specified input images.
    /// </summary>
    /// <param name = "image1"> The first image. </param>
    /// <param name = "image2"> The second image. </param>
    /// <param name = "parameter"> The laplacian parameter list. </param>
    void EstimateLaplacianNoise(const CFloatImg& image1, const CFloatImg& image2, double*& parameter);
    
    /// <summary>
    /// Expands flow from higher level down to the current level using bilinear interpolation. 
    /// </summary>
    /// <param name = "iLevel"> The current level on the pyramid. </param>
    /// <returns>HRESULT</returns>
    HRESULT ExpandFlow(int iLevel);
    
    /// <summary>
    /// Creates feature images for the original images. Feature image contains information about both 
    /// the original image as well as its gradient. 
    /// </summary>
    /// <param name = "imageIndex"> Specify the index of one of the original images to be converted to feature image. </param>
    /// <param name = "iLevel"> The level on the pyramid that the feature image corresponds to. </param>
    /// <returns>HRESULT</returns>
    HRESULT ImageToFeature(int imageIndex, int iLevel);

private:
    // Allocated Flag
    bool m_bAllocated;
    bool m_bSSE1;
    bool m_bSSE2;
    bool m_bSSSE3;
    bool m_bSSE4_1;

    // Input Image and Pyramids
    int m_iCurrentImage;
    int m_iPreviousImage;
    int m_iNumberImages;
    int m_iInputWidth;
    int m_iInputHeight;
    int m_iBaseWidth;
    int m_iBaseHeight;
    int m_iSubsample;
    int m_iLevels;
    int m_iTopLevelMinPixelCount;
    CFloatImg m_imgTmpImage;
    CFloatImg m_imgTmpImage2;
    CFloatImg m_TempImage;
    CFloatImg m_SmoothFlowTempImage1;
    CFloatImg m_SmoothFlowTempImage2;
    CFloatPyramid m_ppyrInputImages[2];
    CFloatPyramid 
        uu, vv, 
        ux, uy, 
        vx, vy, 
        Phi_1st, Psi_1st, 
        imdxy, imdx2, imdy2, imdtdx, imdtdy,
        ImDxy, ImDx2, ImDy2, ImDtDx, ImDtDy;
    CFloatImg m_imgSORX;
    CFloatImg m_imgSORY;
    CFloatImg m_imgSORPhi;

#ifdef USE_FEATURES
    //Feature images
    CFloatPyramid m_ppyrFeatureImages[2];
#endif

    // Intermediate Images
    CFloatPyramid m_pyrDerivsX;
    CFloatPyramid m_pyrDerivsY;
    CFloatPyramid m_pyrWarpedInput;
    CFloatPyramid m_pyrErrorImage;
    CIntImg m_imgiSrcX;
    CIntImg m_imgiSrcY;

    // Flow
    CFloatPyramid m_pyrFlowX;
    CFloatPyramid m_pyrFlowY;
    CFloatPyramid m_pyrDFlowX;
    CFloatPyramid m_pyrDFlowY;
    CFloatPyramid m_pyrFlowLaplacianX;
    CFloatPyramid m_pyrFlowLaplacianY;

    //Noise model
    double* m_LaplacianParameters;

    //Constants
    bool	m_bNormalizePixelValues; //Normalize pixel values to 0-1
    int		m_iOuterIterations; //Outer number of iterations
    int		m_iInnerIterations; //Inner number of iterations
    int		m_iSORIterations; //S.O.R number of iterations
    float	m_Omega; //Omega used in S.O.R calculations
    float   m_Alpha; //Alpha used in S.O.R calculations
    float	m_Rho; //Rho value to be used in 1-Rho-1 filter, only applies if using smooth pyramid
    int		m_iImageWarpMethod; //Method to warp image, see enum ImageWarpMethod definition
    bool	m_bGaussianGradients; //Use smooth gaussian gradients Dx, Dy, Dt for S.O.R
    bool	m_bSmoothDerivatives; //Use universal smooth derivatives
    int		m_iNoiseModel; //Noise model, see enum NoiseModel definition
    bool	m_bUseSSE; //Use SSE computations if available, otherwise disable SSE completely
};

};