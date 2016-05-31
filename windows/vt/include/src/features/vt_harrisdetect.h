//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      specialized implementation of Harris feature detector
//
//------------------------------------------------------------------------
#pragma once

struct HARRIS_DETECTOR_PARAMS
{
    int   border;
    float scoreThreshold;
    float cornerEdgeThreshold;
    float detectThreshold;

    // threshold for luminance value in source image at detection dimension 
    // (after any downsampling): features are not generated when the source
    // pixel at the feature position is at or below the threshold
    unsigned int lowSrcThreshold;

    int   blockSize;
    int   featureCountTarget;

    // number of passes of 2:1 decimation (with a 121 filter kernel) to
    // apply prior to detection; default is 0 (no decimation)
    int   decimate2to1Count;

    HARRIS_DETECTOR_PARAMS()
        :border(5)
        ,scoreThreshold(0.f)
        ,cornerEdgeThreshold(5.f)
        ,detectThreshold(0.f)
        ,lowSrcThreshold(0)
        ,blockSize(64)
        ,featureCountTarget(500)
        ,decimate2to1Count(0)
    {}
};

struct HARRIS_DETECTOR_RESULTDATA
{
    float imageBlurLevel;
    HARRIS_DETECTOR_RESULTDATA()
        :imageBlurLevel(0.f)
    {}
};

struct HARRIS_FEATURE_POINT
{
    float x,y;
    float score;
};

class CHarrisDetector
{
public:
    CHarrisDetector() {};
    ~CHarrisDetector() {};
    HRESULT Detect(vt::vector<vt::vector<HARRIS_FEATURE_POINT>>& pts, 
        const vt::CByteImg& src, const vt::CRect& srcRect,
        const HARRIS_DETECTOR_PARAMS& params,
        HARRIS_DETECTOR_RESULTDATA& resultData);
};

