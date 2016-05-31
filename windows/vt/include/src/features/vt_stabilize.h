//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      stabilization API
//
//------------------------------------------------------------------------------
#pragma once

#include "vtcore.h"
#include "vt_features_base.h"
#include "vt_harrisdetect.h"

//------------------------------------------------------------------------------
// parameters for feature tracker part of stabilizer
//------------------------------------------------------------------------------
struct FEATURE_TRACKER_PARAMS
{
    // maximum feature distance, expressed as a fraction of the max image
    // dimension, used during tracking as a maximum radius for searching
    // for feature matches from adjacent frames
    float fMaxDistance;

    // size of video frame buffer for tracker part of video stabilizer;
    // can be set to zero for 'bufferless' mode, which works well for 
    // most cases but does exhibit some visual jarring artifacts occasionally;
    // when non-zero, the tracker runs RANSAC every 'iBufferSize' frames
    // so the number should be set as large as is reasonable (useful range
    // ~8 to 15)
    int iBufferSize;

    // number of iterations for RANSAC used to differentiate between
    // foreground and background features and to disregard poorly
    // tracked or detected features
    int iRANSACIterations;

    // pixel distance threshold used for culling features that have
    // been matched across one or more frames; expressed in pixel coordinates
    // at the detection dimension
    float fRANSACOutlierDist;

    // fraction of tracked features that can be lost before creating a new ref
    float fTrackLossFraction;

    // we establish the baseline of # tracked features for fTrackLossFraction
    // after iMinTrackLength.  iMinTrackLength will be ignored and a new 
    // track will be started regardless if # tracked features goes below
    // iMinTrackedFeatures  
    int iMinTrackLength;

    // min number of allowed feature tracks
    int iMinTrackedFeatures;

    // max ratio between the best match and the second best match distances
    // if the actual ratio is greater than fMatchTestRatio, the match pair is discarded in BriefFindMatch. 
    float fMatchTestRatio;

    FEATURE_TRACKER_PARAMS() :
        fMaxDistance(0.08f),
        iBufferSize(10),
        fTrackLossFraction(0.5f),
        iMinTrackLength(3),
        iMinTrackedFeatures(32),
        iRANSACIterations(250),
        fRANSACOutlierDist(2.0f),
        fMatchTestRatio(0.85f)
    {}
};

//------------------------------------------------------------------------------
// parameters for feature detector part of stabilizer
//------------------------------------------------------------------------------
struct FEATURE_DETECTOR_PARAMS
{
    // increase or decrease the resolution at which feature detection is 
    // done; size is scaled by 2**(-detectDimensionBias), so is made
    // smaller by positive values and larger by negative values
    int detectDimensionBias;

    // pre-blur image at the detection dimension prior to feature detect;
    // no blur of 0.f; used as sigma for Gaussian blur if non-zero
    float preBlur;

    // detector selection: 0: Harris; 1: FAST10
    int detector;

    // FAST10 feature detection controls
    float FAST10Threshold;

    // Harris feature detection controls
    HARRIS_DETECTOR_PARAMS harrisParams;

    // default parameters
    FEATURE_DETECTOR_PARAMS() :
        detectDimensionBias(0),
        preBlur(0.f),
        detector(0),
        FAST10Threshold(25.f)
    {}
};

//------------------------------------------------------------------------------
// parameters and interface for video stabilizer
//------------------------------------------------------------------------------
struct VIDEO_STABILIZER_PARAMS
{
    // internal parameters for the feature tracker
    FEATURE_TRACKER_PARAMS trackParams;

    // internal parameters for the feature detector
    FEATURE_DETECTOR_PARAMS detectParams;

    // stabilization mode (0-fixed camera, 1-smoothed moving camera, 2-DP alg, 3-no stabilize, 4-weighted least square (wls) alg.)
    unsigned int stabilizeMode;

    // gaussian smooth sigma for mode 1
    float fSmoothSigma;

    // depth of stabilization buffer; range is 5 to ~20
    //
    // total buffering for video stabilizer is the sum of the
    // stabilization buffer and the tracker buffer set in
    // trackParams.iBufferSize
    int stabilizeBuffer;

    // apply rolling shutter correction or not
    bool  bRollingShutter;

    // if true, do automatic cropping of transforms
    bool  AutoCrop;
    float fCropBox;

    // these integers correspond to the eMotionModel enum; supported values are:
    // 2 - eRigidScale
    // 3 - eAffine
    // 4 - e3DRot
    // 5 - eHomography
    // eAffine and eHomography are treated the same by CStabilize
    int   motionModel;

    // default parameters
    VIDEO_STABILIZER_PARAMS() :
        motionModel(2),  // similarity by default
        stabilizeMode(0),
        fSmoothSigma(4),
        stabilizeBuffer(20),
        bRollingShutter(false),
        AutoCrop(false),
        fCropBox(0.9f)
    {}
};

class CStabilizerFrameData
{
public:
    CStabilizerFrameData();
    ~CStabilizerFrameData();

    float GetBlurLevel();
    int GetDetectedFeatureCount();
private:
    class CStabilizerFrameDataInternal;
    CStabilizerFrameDataInternal* m_pint;
public:
    CStabilizerFrameDataInternal* Impl();
};

class CVideoStabilizer: public IDelayResult
{
public:
    int GetMaxDelay();

    BUFFER_RANGE GetResultsRange();

public:
    CVideoStabilizer();
    ~CVideoStabilizer();

    // get current values for parameters
    HRESULT GetParams(VIDEO_STABILIZER_PARAMS& params);

    // call prior to first AddFrame
    // TODO: rename this Begin(), add a SetParameters method to match IFeatureWarpCompute API
    HRESULT Start(const VIDEO_STABILIZER_PARAMS& params);

    // call after last AddFrame - can still call GetTransforms for current buffer contents
    HRESULT End(void);

    // adds a frame; downsamples image as needed to get to the computed detection
    // dimension, runs detector, and computes descriptors prior to returning;
    HRESULT PushFrame(vt::CLumaByteImg& img);

    // variant of PushFrame where the image data is preprocessed to an opaque
    // object that can later be passed to stabilizer; useful when stabilizing
    // a small number of frames for which the tracker might be run multiple
    // times, or if the results of PreProcessFrame (which include the blur level)
    // are to be used to remove frames from a stabilization sequence
    HRESULT PreProcessFrame(vt::CLumaByteImg& img, CStabilizerFrameData& data);
    HRESULT PushFrame(CStabilizerFrameData& data);

    // returns transform for the oldest frame in the buffer (i.e. iBufferSize frames
    // back from the most recent AddFrame) between Start() and End(); after End(),
    // GetTransform call returns remaining transforms in buffer continuing in
    // age-sequential order
    //
    // frameNumber is the zero-based frame count for the returned transform
    HRESULT GetResult(vt::CMtx3x3f& transform, int frameNumber);
    // variant returning the rolling shutter correction
    HRESULT GetResult(vt::CMtx3x3f& transform, vt::vector<vt::CVec2f>& vecRScorrect,
                      int frameNumber);

    // optional callback to return selection of various points and event notifications
#define VSCB_DETECTEDFPT    (1<< 0) // detected feature points (FeaturePoint)
#define VSCB_TRACKEDPT      (1<< 1) // tracked points (CVec2f)
#define VSCB_TRACKEDLEN     (1<< 2) // tracked points length (int)
#define VSCB_INLIERPT       (1<< 3) // inlier points (CVec2f)
#define VSCB_TRSOFTRESETT   (1<< 4) // tracker soft reset due to track loss (no data)
#define VSCB_TRSOFTRESETI   (1<< 5) // tracker soft reset due to inlier loss (no data)
#define VSCB_TRHARDRESET    (1<< 6) // tracker hard reset (no data)
#define VSCB_DETECTEDCNT    (1<< 7) // count of detected feature points (int)
#define VSCB_BLURLEVEL      (1<< 8) // blur level (float)
#define VSCB_TRENDBATCH     (1<< 9) // tracker end of batch (no data)
void SetCallback(unsigned int mask, 
                 void (__stdcall *proc)(unsigned int type, const void* data, int frame));

private:
    class CVideoStabilizerInternal;
    CVideoStabilizerInternal* m_pint;
};

//------------------------------------------------------------------------------
// parameters and interface for image stack stabilizer
//------------------------------------------------------------------------------
struct STACK_STABILIZER_PARAMS
{
    // internal parameters for the feature tracker
    FEATURE_TRACKER_PARAMS trackParams;

    // internal parameters for the feature detector
    FEATURE_DETECTOR_PARAMS detectParams;

    // if true, do automatic cropping of transforms
    bool AutoCrop;
    float fCropbox;

    // these integers correspond to the eMotionModel enum; supported values are:
    // 2 - eRigidScale
    // 3 - eAffine
    // 4 - e3DRot
    // 5 - eHomography
    // eAffine and eHomography are treated the same by CStabilize
    int   motionModel;

    // debugging support
    // TODO: remove when done with development
    unsigned int debugFlags; int debugI0; float debugF0,debugF1,debugF2;

    STACK_STABILIZER_PARAMS() :
        motionModel(2),  // similarity by default
        AutoCrop(false),
        debugFlags(0),debugI0(0),debugF0(0.f),debugF1(0.f),debugF2(0.f)
    {}
};

class CStackStabilizer
{
public:
    CStackStabilizer();
    ~CStackStabilizer();

    // get current values for parameters
    HRESULT GetParams(STACK_STABILIZER_PARAMS& params);

    // set current values for parameters - call only prior to AddReferenceFrame
    HRESULT SetParams(const STACK_STABILIZER_PARAMS& params);

    // sets reference frame against which subsequent AlignFrame calls are paired
    //
    // detectedFeatures returns the number of features detected in that frame
    HRESULT SetReferenceFrame(vt::CLumaByteImg& img);

    // computes alignment transform between this frame image and the current
    // reference frame image; returns E_FAIL if images cannot be aligned with
    // current parameter (which does not necessarily mean that they could be
    // aligned with different settings)
    HRESULT AlignFrame(vt::CLumaByteImg& img, vt::CMtx3x3f& transform);

    // optional callback to return selection of various points and event notifications
#define SSCB_DETECTEDFPT  VSCB_DETECTEDFPT 
#define SSCB_INLIERPT     VSCB_INLIERPT     
#define SSCB_TRSOFTRESETT VSCB_TRSOFTRESETT 
#define SSCB_TRSOFTRESETI VSCB_TRSOFTRESETI 
#define SSCB_TRHARDRESET  VSCB_TRHARDRESET  
#define SSCB_DETECTEDCNT  VSCB_DETECTEDCNT
#define SSCB_BLURLEVEL    VSCB_BLURLEVEL
    void SetCallback(unsigned int mask, void (__stdcall *proc)(unsigned int type, const void* data, int frame));

private:
    class CStackStabilizerInternal;
    CStackStabilizerInternal* m_pint;
};

