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
//
//  Enable code to compute weighting functions
#define USE_WEIGHTED_DATA_AND_SMOOTHNESS
//  Turn on debugging code
// #define OUTPUT_OUTER_ITERATION_ERROR
// #define OUTPUT_INNER_ITERATION_ERROR
//
//----------------------------------------------------------------------------

#pragma once
#include "iopticalflow.h"
#include "vtcore.h"
namespace vt
{

#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
typedef void VOFWeightingFn(CFloatImg &imgSmoothness, float *pfSmoothnessParams);
#endif

/// <summary> The internal parameters of a variational optical flow algorithm. Not all these
/// parameters are exposed in the IOpticalFlow interface. </summary>	
struct CVOFParams
{
    int m_iPreBlurCount;
    int m_iPyramidTopLevelPixelCount;
    float *m_pfSmoothnessParams;
    float m_fW;
    int m_iSpeed;
    int m_iMaxIterationLevel;
    int *m_piInnerIterationsFast;
    int *m_piOuterIterationsFast;
    int *m_piInnerIterationsFull;
    int *m_piOuterIterationsFull;
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
    VOFWeightingFn *m_pfnWeightingFn;
    CVOFParams(int iPreBlurCount, int iPyramidTopLevelPixelCount, 
        float *pfSmoothnessParams, float fW, int iSpeed, int iMaxIterationLevel, 
        int *piInnerIterationsFast, int *piOuterIterationsFast, int *piInnerIterationsFull, 
        int *piOuterIterationsFull, VOFWeightingFn *pfnWeightingFn) :
        m_iPreBlurCount(iPreBlurCount),
        m_iPyramidTopLevelPixelCount(iPyramidTopLevelPixelCount),
        m_pfSmoothnessParams(pfSmoothnessParams),
        m_fW(fW),
        m_iSpeed(iSpeed),
        m_iMaxIterationLevel(iMaxIterationLevel),
        m_piInnerIterationsFast(piInnerIterationsFast),
        m_piOuterIterationsFast(piOuterIterationsFast),
        m_piInnerIterationsFull(piInnerIterationsFull),
        m_piOuterIterationsFull(piOuterIterationsFull),
        m_pfnWeightingFn(pfnWeightingFn) {}
#else
    CVOFParams(int iPreBlurCount, int iPyramidTopLevelPixelCount, 
        float *pfSmoothnessParams, float fW, int iSpeed, int iMaxIterationLevel, 
        int *piInnerIterationsFast, int *piOuterIterationsFast, int *piInnerIterationsFull, 
        int *piOuterIterationsFull) :
        m_iPreBlurCount(iPreBlurCount),
        m_iPyramidTopLevelPixelCount(iPyramidTopLevelPixelCount),
        m_pfSmoothnessParams(pfSmoothnessParams),
        m_fW(fW),
        m_iSpeed(iSpeed),
        m_iMaxIterationLevel(iMaxIterationLevel),
        m_piInnerIterationsFast(piInnerIterationsFast),
        m_piOuterIterationsFast(piOuterIterationsFast),
        m_piInnerIterationsFull(piInnerIterationsFull),
        m_piOuterIterationsFull(piOuterIterationsFull) {}
#endif
};

/// <summary> Horn-Schunck 1 Channel
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
/// </summary>
extern CVOFParams HORN_SCHUNCK_1CHANNEL_PARAMS;

/// <summary> An implementation of a subset of variational optical flow algorithms
/// </summary>
class CVOF : public IOpticalFlow
{
public:
    /// <summary> Constructor </summary>
    CVOF() { Init(); }
    /// <summary> Destructor </summary>
    ~CVOF() { Deallocate(); }
    /// <summary> Allocates space for the object </summary>
    /// <param name="iImageWidth"> Width of images to compute flow on </param>
    /// <param name="iImageHeight"> Height of images to compute flow on </param>
    /// <param name="cParam"> A pointer to a CVOFParams object containing the variational
    /// optical flow parameters </param>
    HRESULT Allocate(int iInputWidth, int iInputHeight, CVOFParams *cParam);
    /// <summary> Deallocates all memory </summary>
    void Deallocate();
    /// <summary> Returns whether the optical flow algorithm has been allocated
    /// successfully </summary>
    bool Allocated() { return m_bAllocated; }	

    /// <summary> Adds a 32bit color images to the optical flow algorithm. The pixels
    /// must be laid out BGRA. VT CRGBAImgs and RGB32 Video Data is laid out this way. </summary>
    /// <param name="pucBuffer"> Pointer to Pixel Buffer </param>
    /// <param name="iBufferWidth"> Width of Image in Buffer </param>
    /// <param name="iBufferHeight"> Height of Image in Buffer </param>
    /// <param name="uiBufferStrideBytes"> Number of bytes per stride</param>
    HRESULT AddBGRAImage(unsigned char* pucBuffer, int iBufferWidth, int iBufferHeight, 
                         size_t uiBufferStrideBytes);
    /// <summary> Adds a 32Bit CRGBAImg to the optical flow algorithm.</summary>
    HRESULT AddImage(CRGBAImg &imgInput)
    {
        return AddBGRAImage(imgInput.BytePtr(), imgInput.Width(), imgInput.Height(), imgInput.StrideBytes());
    }
    /// <summary> Adds an 8bit image to the optical flow algorithm. </summary>
    /// <param name="pucBuffer"> Pointer to Pixel Buffer </param>
    /// <param name="iBufferWidth"> Width of Image in Buffer </param>
    /// <param name="iBufferHeight"> Height of Image in Buffer </param>
    /// <param name="uiBufferStrideBytes"> Number of bytes per stride</param>
    HRESULT AddYImage(unsigned char* pucBuffer, int iBufferWidth, int iBufferHeight, 
                      size_t uiBufferStrideBytes);
    /// <summary> Adds an 8bit CLumaImg to the optical flow algorithm. </summary>
    HRESULT AddImage(CLumaImg &imgInput)
    {
        return AddYImage(imgInput.BytePtr(), imgInput.Width(), imgInput.Height(), imgInput.StrideBytes());
    }
    /// <summary> Adds an 8bit float image to the optical flow algorithm. </summary>
    /// <param name="pfBuffer"> Pointer to Pixel Buffer </param>
    /// <param name="iBufferWidth"> Width of Image in Buffer </param>
    /// <param name="iBufferHeight"> Height of Image in Buffer </param>
    /// <param name="uiBufferStrideBytes"> Number of bytes per stride</param>
    HRESULT AddYImage(float* pfBuffer, int iBufferWidth, int iBufferHeight, 
                      size_t uiBufferStrideBytes);

    /// <summary> Computes the flow at a specific level in the pyramid </summary>
    /// <param name="iLowestLevel"> The level in the pyramid down to which the flow should be
    /// computed. iLowestLevel must be less than the number of levels in the pyramidn returned by
    /// GetNumberOfLevels().</param>
    HRESULT ComputeFlow(int iLowestLevel = 0);

    /// <summary> Returns the number of levels in the pyramid. </summary>
    int GetNumberOfLevels() { return m_iLevels; }
    /// <summary> Gets a pointer to a CFloatImg containing the Y component of the flow. </summary>
    /// <param name="pimgFlowX"> A pointer that receives the flow. </param>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
    HRESULT GetFlowX(CFloatImg* &pimgFlowX, int iLevel);
    /// <summary> Gets a pointer to a CFloatImg containing the X component of the flow. </summary>
    /// <param name="pimgFlowY"> A pointer that receives the flow. </param>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
    HRESULT GetFlowY(CFloatImg* &pimgFlowY, int iLevel);
    /// <summary> Gets a buffer containing the internal representation of the last image
    /// added to the optical flow object. </summary>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
    /// <returns> 
    ///  - pfBuffer: A pointer to the buffer
    ///  - iBufferWidth: The width of the buffer
    ///  - iBufferHeight: The height of the buffer
    ///  - uiBufferStrideBytes: The stride of the buffer in bytes
    /// </returns>
    HRESULT GetCurrentImage(float* &pfBuffer, int &iBufferWidth, int &iBufferHeight, 
                            size_t &uiBufferStrideBytes, int iLevel);
    /// <summary> Gets a buffer containing the internal representation of the previous last image
    /// added to the optical flow object. </summary>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
    /// <returns> 
    ///  - pfBuffer: A pointer to the buffer
    ///  - iBufferWidth: The width of the buffer
    ///  - iBufferHeight: The height of the buffer
    ///  - uiBufferStrideBytes: The stride of the buffer in bytes
    /// </returns>
    HRESULT GetPreviousImage(float* &pfBuffer, int &iBufferWidth, int &iBufferHeight, 
                             size_t &uiBufferStrideBytes, int iLevel);
    /// <summary> Gets a buffer containing the internal representation of the current image
    /// warped back onto the previous image. </summary>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
    /// <returns> 
    ///  - pfBuffer: A pointer to the buffer
    ///  - iBufferWidth: The width of the buffer
    ///  - iBufferHeight: The height of the buffer
    ///  - uiBufferStrideBytes: The stride of the buffer in bytes
    /// </returns>
    HRESULT GetWarpedImage(float* &pfBuffer, int &iBufferWidth, int &iBufferHeight, 
                           size_t &uiBufferStrideBytes, int iLevel);
    /// <summary> Gets a buffer containing the internal representation of the error image; ie. the
    /// differnce between the previous images and the current image warped back onto the previous 
    /// image. </summary>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
    /// <returns> 
    ///  - pfBuffer: A pointer to the buffer
    ///  - iBufferWidth: The width of the buffer
    ///  - iBufferHeight: The height of the buffer
    ///  - uiBufferStrideBytes: The stride of the buffer in bytes
    /// </returns>
    HRESULT GetErrorImage(float* &pfBuffer, int &iBufferWidth, int &iBufferHeight, 
                          size_t &uiBufferStrideBytes, int iLevel);

protected:
    /// <summary> Initializes the variable in the flow object. Performs no allocaton. </summary>
    void Init();
    /// <summary> Computes flow at a specific level in the pyramid, using a fixed maxiumum
    /// number of inner and outer iterations. Inner iterations just perform an iteration of the
    /// linear solver. Outer iterations warp the current image onto the previous, compute the
    /// error, linearize the flow constraint equations, and solve. </summary>
    HRESULT ComputeFlowForLevel(int iLevel, int iMaxOuterIterations, int iMaxInnerIterations);
    /// <summary> Upsamples the flow from level iLevel to iLevel-1 in the pyramid. Used to
    /// initialize the flow at the next level in the coarse to fine framework. </summary>
    HRESULT ExpandFlow(int iLevel);
    /// <summary> Computes the difference between the previous image and the current image warped
    /// backwards onto the previous image using the current estimate of the flow. ComputeErrorImage()
    /// returns the current estimate of the error in the outer loop if the appropriate code is
    /// turned on with an appropriate #define option. </summary>
    HRESULT ComputeErrorImage(int iLevel, float &fOuterError)
    {
        if (m_bSSE2)
        {
            return ComputeErrorImageSSE2(iLevel, fOuterError);
        }
        else
        {
            return ComputeErrorImageNoSSE(iLevel, fOuterError);
        }
    }
    /// <summary> SSE2 Implementation of ComputeErrorImage() </summary>
    HRESULT ComputeErrorImageSSE2(int iLevel, float &fOuterError);
    /// <summary> Straight C++ Implementation of ComputeErrorImage() </summary>
    HRESULT ComputeErrorImageNoSSE(int iLevel, float &fOuterError);

    /// <summary> Computes the X and Y derivatives of the previous image </summary>
    HRESULT ComputeDerivativesForLevel(int iLevel);
    /// <summary> Computes the Laplacian of the flow; ie. the error in the smoothness
    /// constraint at each pixel. </summary>
    HRESULT ComputeFlowLaplacian(int iLevel)
    {
        if (m_bSSE2 && m_pyrFlowX.GetLevel(iLevel).Width() > 4)
        {
            return ComputeFlowLaplacianSSE2(iLevel);
        }
        else
        {
            return ComputeFlowLaplacianNoSSE(iLevel);
        }
    }
    /// <summary> SSE2 Implementation of ComputeFlowLaplacian() </summary>
    HRESULT ComputeFlowLaplacianSSE2(int iLevel);
    /// <summary> Straight C++ Implementation of ComputeFlowLaplacian() </summary>
    HRESULT ComputeFlowLaplacianNoSSE(int iLevel);
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
    /// <summary> Computes a weighting function for none L2 penalty functions </summary>
    HRESULT ComputeWeights(int iLevel);
#endif
    /// <summary> Computes and simplifies the linear system defined by the linearized
    /// optical flow constraints. </summary>
    HRESULT ComputeLinearSystem(int iLevel)
    {
		if (m_pfnWeightingFn != NULL)
		{
			if (m_bSSE2)
			{
				return ComputeLinearSystemSSE2(iLevel);
			}
			else
			{
				return ComputeLinearSystemNoSSE(iLevel);
			}
		}
		else
		{
			if (m_bSSE2)
			{
				return ComputeLinearSystemSSE2NoWeight(iLevel);
			}
			else
			{
				return ComputeLinearSystemNoSSENoWeight(iLevel);
			}
		}
    }
    /// <summary> SSE2 Implementation of ComputeLinearSystem() </summary>
    HRESULT ComputeLinearSystemNoSSE(int iLevel);
    /// <summary> Straight C++ Implementation of ComputeLinearSystem() </summary>
    HRESULT ComputeLinearSystemSSE2(int iLevel);
    /// <summary> SSE2 Implementation of ComputeLinearSystem() NoWeight </summary>
    HRESULT ComputeLinearSystemNoSSENoWeight(int iLevel);
    /// <summary> Straight C++ Implementation of ComputeLinearSystem() NoWeight </summary>
    HRESULT ComputeLinearSystemSSE2NoWeight(int iLevel);
    /// <summary> Performs 1 iteration of Successfive Over-Relaxation to solve the linear
    /// system. See http://en.wikipedia.org/wiki/Successive_over-relaxation. </summary>
    HRESULT SOR(int iLevel, float &fTotalError)
    {
        if (m_bSSE2 && m_pyrFlowX.GetLevel(iLevel).Width() > 4)
        {
            return SORSSE2(iLevel, fTotalError);
        }
        else
        {
            return SORNoSSE(iLevel, fTotalError);
        }
    }
    /// <summary> SSE2 Implementation of SOR() </summary>
    HRESULT SORSSE2(int iLevel, float &fTotalError);
    /// <summary> Straight C++ Implementation of SOR() </summary>
    HRESULT SORNoSSE(int iLevel, float &fTotalError);
    /// <summary> Updates the flow using the inverse compositional update. The DFLOW is
    /// subtracted from the current estimate of the flow. DFLOW is also set to be zero,
    /// read for the next outer iteration. </summary>
    void UpdateDFlow(CFloatImg &imgFlow, CFloatImg &imgDFlow)
    {
        if (m_bSSE2)
        {
            UpdateDFlowSSE2(imgFlow, imgDFlow);
        }
        else
        {
            UpdateDFlowNoSSE(imgFlow, imgDFlow);
        }
    }
    /// <summary> SSE2 Implementation of UpdateDFlow() </summary>
    void UpdateDFlowSSE2(CFloatImg &imgFlow, CFloatImg &imgDFlow);
    /// <summary> Straight C++ Implementation of UpdateDFlow() </summary>
    void UpdateDFlowNoSSE(CFloatImg &imgFlow, CFloatImg &imgDFlow);
    /// <summary> Used to warp the current image back onto the previous image using the
    /// current estimate of the flow. The data term weight is also set by setting all pixels
    /// to be 0 for which the destination is outside the image during the warp. A border of
    /// pixels around the edge of the image are also set to have zero data term weight. </summary>
    HRESULT WarpImage(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, CFloatImg &imgDst, 
                      CFloatImg &imgDataTermWeight, int iLevel = 0)
    {
//		VT_HR_BEGIN();

        if (m_bSSE2 && m_bSSSE3)
        {
//			if(imgDst.Width() >= 96 && imgDst.Height() >= 96 )
//			{
//				VT_HR_EXIT(WarpImageCustomWarp(imgSrc, imgFlowX, imgFlowY, imgDst, imgDataTermWeight));
//			}
//			else
            {
//                m_imgTmpImage.Fill(0.0f);
//                vt::OpticalFlow_WarpImageSSE(imgSrc, imgFlowX, imgFlowY, imgDst, imgDataTermWeight, imgDst.Rect(), m_imgTmpImage);
				WarpImageSSE3(imgSrc, imgFlowX, imgFlowY, imgDst, imgDataTermWeight, imgDst.Rect());
            }
        }
        else
        {
            WarpImageNoSSE(imgSrc, imgFlowX, imgFlowY, imgDst, imgDataTermWeight);
        }
        MaskDataTerm(imgDataTermWeight, iLevel);

//		VT_HR_END();
        return S_OK;
    }

public:
	void WarpImageSSE3(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, CFloatImg &imgDst, 
				       CFloatImg &imgDataTermWeight, CRect rctDst);

	/// <summary> A straight C++ implementation of WarpImage() </summary>
    void WarpImageNoSSE(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, CFloatImg &imgDst, 
                       CFloatImg &imgDataTermWeight);

private:
	/// <summary> A function that sets the data term weight to be zero around the border of the image.
    /// </summary>
    void MaskDataTerm(CFloatImg &imgDataTermWeight, int iLevel);

    /// <summary> Computes the iteration counts from the speed </summary>
    void SetInterationCountsFromSpeed();

private:
//	CVTResourceManagerInit *m_ri;

    /// <summary> Miscellaneous Flags. </summary>
    bool m_bAllocated;
    bool m_bSSE2;
    bool m_bSSSE3;
    bool m_bSSE4_1;

    /// <summary> Member variables related to the input images and pyramids. </summary>
    int m_iCurrentImage;
    int m_iPreviousImage;
    int m_iNumberImages;
    int m_iInputWidth;
    int m_iInputHeight;
    int m_iBaseWidth;
    int m_iBaseHeight;
    int m_iSubsample;
    int m_iLevels;
    int m_iPyramidTopLevelPixelCount;
    int m_iPreBlurCount;
    CFloatImg m_imgTmpImage;
    CFloatImg m_imgTmpImage2;
    CLumaFloatPyramid m_ppyrInputImages[2];

    /// <summary> Miscellaneous intermediate image pyramids. </summary>
    CLumaFloatPyramid m_pyrDerivsX;
    CLumaFloatPyramid m_pyrDerivsY;
    CLumaFloatPyramid m_pyrWarpedInput;
    CLumaFloatPyramid m_pyrErrorImage;

    /// <summary> The flow image pyramids. </summary>
    CLumaFloatPyramid m_pyrFlowX;
    CLumaFloatPyramid m_pyrFlowY;
    CLumaFloatPyramid m_pyrDFlowX;
    CLumaFloatPyramid m_pyrDFlowY;
    CLumaFloatPyramid m_pyrFlowLaplacianX;
    CLumaFloatPyramid m_pyrFlowLaplacianY;

    /// <summary> Member variables related to objective function and regularization/smoothness. </summary>
    float *m_pfSmoothnessParams;
    CLumaFloatPyramid m_pyrDataTermWeight;
#ifdef USE_WEIGHTED_DATA_AND_SMOOTHNESS
    CLumaFloatPyramid m_pyrSmoothnessWeight;
    VOFWeightingFn *m_pfnWeightingFn;
#endif

    /// <summary> Member variables related to the linear system and its solution. </summary>
    float m_fW;
    float m_fOneMinusW;
    int m_iMaxIterationLevel;
    int m_iSpeed;
    int *m_piInnerIterations;
    int *m_piOuterIterations;
    int *m_piInnerIterationsFast;
    int *m_piOuterIterationsFast;
    int *m_piInnerIterationsFull;
    int *m_piOuterIterationsFull;
    CLumaFloatPyramid m_pyrLinearSystem1;
    CLumaFloatPyramid m_pyrLinearSystem2;
    CLumaFloatPyramid m_pyrLinearSystem3;
    CLumaFloatPyramid m_pyrLinearSystem4;
    CLumaFloatPyramid m_pyrLinearSystem5;
    CLumaFloatPyramid m_pyrLinearSystem6;
    CFloatImg m_imgSORX;
    CFloatImg m_imgSORY;
#ifdef OUTPUT_INNER_ITERATION_ERROR
    CLumaFloatPyramid m_pyrInnerErrorFactorX;
    CLumaFloatPyramid m_pyrInnerErrorFactorY;
#endif
};

};