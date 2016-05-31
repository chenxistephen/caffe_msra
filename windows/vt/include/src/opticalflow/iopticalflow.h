//+---------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation 2010.  All rights reserved.
//
//  Description:
//     Interface for Optical Flow Algorithms
//
//  History:
//     2011/8/12-wkienzle/sbaker
//          Created
//
//  Example Usage (Error Handling Omitted)
//
//      // Create OF Object
//      IOpticalFlow *pHSOpticalFlow;
//      CreateHornSchunck1Channel(pHSOpticalFlow, iWidth, iHeight);
//
//      // Load Images
//      CRGBAImg imgImage1;
//      CRGBAImg imgImage2;
//      VtLoadImage(argv[1], imgImage1);
//      VtLoadImage(argv[2], imgImage2);
//
//      // Add Images to Optical Flow Object
//      pHSOpticalFlow->AddImage(imgImage1);
//      pHSOpticalFlow->AddImage(imgImage2);
//
//      // Compute Flow
//      pHSOpticalFlow->ComputeFlow();
//
//      // Get the Flow
//		CFloatImg *pimgFlowX, *pimgFlowY;
//      pHSOpticalFlow->GetFlowX(pimgFlowX);
//      pHSOpticalFlow->GetFlowY(pimgFlowY);
//
//
//+----------------------------------------------------------------------------------------------
#pragma once
#include "vtcore.h"

namespace vt {

/// \ingroup opticalflow
/// <summary>Interface for optical flow computations</summary>
class IOpticalFlow
{
public:
    /// <summary> Destructor </summary>
	virtual ~IOpticalFlow() {}

    /// <summary> Adds a 32Bit CRGBAImg to the optical flow algorithm.</summary>
	virtual HRESULT AddImage(CRGBAImg &imgInput) = 0;
    /// <summary> Adds an 8bit CLumaImg to the optical flow algorithm. </summary>
	virtual HRESULT AddImage(CLumaImg &imgInput) = 0;

    /// <summary> Computes the flow at a specific level in the pyramid </summary>
    /// <param name="iLowestLevel"> The level in the pyramid down to which the flow should be
	/// computed. iLowestLevel must be less than the number of levels in the pyramidn returned by
	/// GetNumberOfLevels(). Algorithms that do not use a pyramid return 1 for GetNumberOfLevels().
	/// For such algorithms, the only option for iLowestLevel is 0.</param>
	virtual HRESULT ComputeFlow(int iLowestLevel = 0) = 0;

    /// <summary> Gets a pointer to a CFloatImg containing the X component of the flow. </summary>
    /// <param name="pimgFlowX"> A pointer that receives the flow. </param>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
	virtual HRESULT GetFlowX(CFloatImg* &pimgFlowX, int iLevel = 0) = 0;
    /// <summary> Gets a buffer containing the X component of the flow. </summary>
    /// <param name="pimgFlowY"> A pointer that receives the flow. </param>
    /// <param name="iLevel"> The level in the pyramid for which the flow is desired. </param>
	virtual HRESULT GetFlowY(CFloatImg* &pimgFlowY, int iLevel = 0) = 0;

	/// <summary> Returns the number of frames that must be added to the buffer before normal
	/// flow computation can begin. The number of times AddImage() must be called before normal 
	/// operation can begin. For exmaple. If GetBufferPreloadCount() returns 1, then  AddImage()
	/// must be first called to load an image. Then images can be added and flow computed after 
	/// each one. E.g. the user can call AddImage(), ComputeFlow(), AddImage(), ComputeFlow(), etc. 
	/// After the final images is added, ComputeFlow() can be called GetBufferPreloadCount() more 
	/// times to drain the buffer. The main purpose of GetBufferPreloadCount() is to support 
	/// algorithms that use more than 2 frames to estimate flow. In the common 2-frame case,  
	/// GetBufferPreloadCount() returns 1, and the calls should be interleaved: AddImage(), 
	/// AddImage(), ComputeFlow(), AddImage(), ComputeFlow(), ..., AddImage(), ComputeFlow().
    /// </summary>
 	virtual int GetBufferPreloadCount() { return 1; }
    /// <summary> Returns the number of levels in the pyramid. If the algorithm does not use a
	/// pyramid representation, GetNumberOfLevels() returns 1. This allows the user to then
	/// compute flow using ComputeFlow() at a level higher than the bottom level.
	/// </summary>
	virtual int GetNumberOfLevels() { return 1; };
};

/// \ingroup opticalflow
/// <summary> Creates a 1-channel implementation of Horn-Schunck Optical Flow. The 
/// underlying representation of the image is a single channel. RGB images
/// are first converted to a single Y channel and optical flow esimated on that.
/// </summary>
/// <param name="pPyramidFlow"> Pointer that receives the created object </param>
/// <param name="iImageWidth"> Width of images to compute flow on </param>
/// <param name="iImageHeight"> Height of images to compute flow on </param>
/// <param name="iSpeed"> (default 4). When Speed = 5, the number of interations of the 
/// underlying solvers is reduced to an absolute minimum, giving the fastest
/// possible execution time. When Speed is reduced to 4, 3, 2, 1, or 0, more
/// iterations are performed. The quality of the flow steady improves, but the
/// algorithm takes more and more time. </param>
/// <param name="fLambda"> (default 500.f) is the smoothing parameter. Increasing 
/// lambda results in more and more smoothing. As Horn-Schunck usese L2 penalty 
/// functions, the flow fields are always somewhat over-smoothed across depth/motion 
/// discontinuities. However, Horn-Schunck flow is sufficient, however, for many video 
/// enhancement applications such as stabilization and rolling shutter correction. </param>
HRESULT CreateHornSchunck1Channel(IOpticalFlow* &pPyramidFlow, int iImageWidth, int iImageHeight,
							      int iSpeed = 4, float fLambda = 500.0f);

/// \ingroup opticalflow
/// <summary> Creates an implementation of Optical Flow based on the algorithm
/// described in Ce Liu's PhD thesis (Appendix A): 
/// http://people.csail.mit.edu/celiu/Thesis/CePhDThesis.pdf. 
/// The input image should be in RGBA format, although the underlying flow computation 
/// can also convert it to a 1-channel image.
/// </summary>
/// <param name="pPyramidFlow"> Pointer that receives the created object </param>
/// <param name="iImageWidth"> Width of images to compute flow on </param>
/// <param name="iImageHeight"> Height of images to compute flow on </param>
HRESULT CreateOpticalFlow2(IOpticalFlow* &pPyramidFlow, int iImageWidth, int iImageHeight);

/// \ingroup opticalflow
/// <summary>
/// SSE implementation of warp image for OpticalFlow-specific algorithm.
/// </summary>
/// <param name="imgSrc"> The source image to be warped from. </param>
/// <param name="imgFlowX"> The x-flow image. </param>
/// <param name="imgFlowY"> The y-flow image. </param>
/// <param name="imgDst"> The destination image to warp to. </param>
/// <param name="imgDataTermWeight">
/// The weight image that will contain pixel-specific weighing values after warping.
/// If a pixel's warped location is out of bound, its weight is zero.
/// </param>
/// <param name="rctDst"> A CRect structure describing the region of warping. </param>
/// <param name="imgSrcExtend">
/// An image containing data that will be used when warping needs to determine values 
/// for out-of-bound locations. For example, once warped, if pixel (x,y) becomes out 
/// of bound, the destination pixel at (x,y) will be copied from the corresponding 
/// location of this extend image.
/// </param>
void OpticalFlow_WarpImageSSE(CFloatImg &imgSrc, CFloatImg &imgFlowX, CFloatImg &imgFlowY, 
    CFloatImg &imgDst, CFloatImg &imgDataTermWeight, CRect rctDst, CFloatImg& imgSrcExtend);
};
