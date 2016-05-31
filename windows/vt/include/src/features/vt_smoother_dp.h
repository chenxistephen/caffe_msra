//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description: implementation of IFeatureWarpCompute that uses a 
//               dynamic-programming style smoothing algo on similarity parameters
//
//------------------------------------------------------------------------------
#include "vt_features_base.h"

#include "transforms.h"

using namespace vt;

//#define PATH_CSV

struct PATH_STATE
{
    float fDistance;
    float fVelocity;
    float fAcceleration;
    short iStartupCount;
    
    PATH_STATE() : fDistance(0), fVelocity(0), fAcceleration(0), iStartupCount(0)
    {}
};

struct DPSMOOTH_PARAMS
{
	int   frameWidth;
	int   frameHeight;
    int   motionModel;
    int   iBufferSize;
    float fCrop;

    DPSMOOTH_PARAMS() : iBufferSize(30), fCrop(0.9f)
	{}
};

class CDPPathSmoother: public IFeatureWarpCompute
{
public:
	CDPPathSmoother()
    {}

    // TODO: add assert that SetParameters is not called between Begin/End
    virtual HRESULT SetParameters(const DPSMOOTH_PARAMS& p)
    { 
		m_params = p;
		return S_OK;
	}

    virtual int GetMaxDelay()
    { return m_params.iBufferSize; }

    virtual BUFFER_RANGE GetResultsRange();

    virtual HRESULT Begin();

    virtual HRESULT PushFeatures(const vt::PointMatch* pMatches,
                                 int iMatchCount);

    virtual HRESULT GetResult(vt::CMtx3x3f& xfrm, int frameNumber);

    virtual HRESULT End();

protected:
	HRESULT DPFilterTransform(CMtx3x3f& result, int iDst);

protected:
    DPSMOOTH_PARAMS m_params;

    // the current cummulative state for DP alg.
    PATH_STATE m_dpstateCur[4];

	// buffers to hold the similarity parameters used in the smoothing case,
    // and to hold the resulting stabilized transforms
    CRollingBuffer<FEAT_SIMILARITY>  m_bufTrackerSim;
    CRollingBuffer<CMtx3x3f>         m_bufResults;
};

