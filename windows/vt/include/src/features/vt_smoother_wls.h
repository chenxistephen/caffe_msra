//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description: implementation of IFeatureWarpCompute that uses a 
//               weighted-least-square style smoothing algo on homography parameters
//
//------------------------------------------------------------------------------
#pragma once

#include <vtcore.h>
#include "vt_features_base.h"

#define SIGMA_MODIFIED

struct WLSSMOOTH_PARAMS
{
    int   frameWidth;
    int   frameHeight;
    int   motionModel;
    int   iBufferSize;
    float cropRatio;

    // Kernel size of gaussian filter
    int   kernelSize;

    // number of output frames per each path planning operation.
    // We use a single-optimzation-with-batch-output approach to improve the performance, that is 
    // after each optimization, we'll output several neighboring frames instead of 1.
    // To remove the gaps (or jitters) between pathes from each optimization, we keep recent pathes and
    // smooth them before output. The cache size parameter specifies how many pathes to be kept and to be used by 
    // smoothing.
    // An imperical value for output frames and cache size is 5.
    int   outFramesPerOptimization;

    // number of warping sets in cache.
    int   cacheSize;

    WLSSMOOTH_PARAMS() : frameWidth(0), frameHeight(0), motionModel(0), iBufferSize(0), cropRatio(0.0f), kernelSize(0), outFramesPerOptimization(0), cacheSize(0)
    {}
};

// forward declaration
class WarpsetWrapper;

class CWLSSmoother : public IFeatureWarpCompute
{
public:
    CWLSSmoother() {}

    virtual HRESULT SetParameters(const WLSSMOOTH_PARAMS& p)
    { 
        m_params = p;
        return S_OK;
    }

    virtual int GetMaxDelay()
    { return m_params.iBufferSize; }

    virtual BUFFER_RANGE GetResultsRange();

    virtual HRESULT Begin();

    virtual HRESULT PushFeatures(const vt::PointMatch* pMatches, int iMatchCount);

    virtual HRESULT GetResult(vt::CMtx3x3f& xfrm, int frameNumber);

    virtual HRESULT End();

private:

    HRESULT AdaptiveSmoothTransform(vt::CMtx3x3f& result, int iDst);

    ///</summary> 
    /// Calculate warping matrixes for the given motion collectino.
    ///</summary>
    ///<param name="motionList">initial transform between frame i-1 to i, i.e., F_i.</param>
    ///<param name="frameIndexInVideo">Starting frame index in video.</param>
    ///<param name="outputFrameIndexInMotionList">Index of output frame in the motion collection.</param>
    HRESULT Smooth(const vt::vector<vt::CMtx3x3f>& motionList, size_t frameIndexInVideo, size_t outputFrameIndexInMotionList);

    // remove residual jitters after enforcing hard crop constraints
    HRESULT PostSmoothing(vt::vector<vt::CMtx3x3f>& updateSet, const vt::vector<vt::CMtx3x3f>& motionList, const vt::vector<float>& lamdaList);

    ///</summary> 
    /// Get optimized warping matrix collection.
    ///</summary>
    ///<param name="frameIndex">Retrieve warping matrix by frame index.</param>
    const vt::CMtx3x3f& GetUpdate(size_t frameIndex) const;

private:
    WLSSMOOTH_PARAMS m_params;

    // buffers to hold the similarity parameters used in the smoothing case,
    // and to hold the resulting stabilized transforms
    CRollingBuffer<vt::CMtx3x3f>  m_bufTrackerData;
    CRollingBuffer<vt::CMtx3x3f>  m_bufResults;

    vt::vector<vt::CMtx3x3f> m_updateSet;  // warping matrix by mapping input frame to output, i.e., B_i

#ifdef SIGMA_MODIFIED		// use an improved and accelerated parameter setting
    static const int	m_nIterations = 10;		// iteration number
#else
    static const int	m_nIterations = 20;
#endif
    
    // number of optimizations that have been performed.
    // this value is used to retrieve the warping set from cache.
    size_t m_optimizationCount;

    // cache of motion sets in order to generate a smoothed motion frames.
    vt::vector<WarpsetWrapper> m_warpsetCache;
};

///<summary> 
/// A wrapper of warping matrix collection generated from path planning.
/// This wrapper keeps the frame index in input video where the warping set starts from and
/// provides interfaces to retrieve warp matrix via frame index.
///</summary>
class WarpsetWrapper
{
public:
    ///</summary> 
    /// Constructor.
    ///</summary>
    WarpsetWrapper() {}

    ///<summary>
    /// Initialize a warping set.
    ///</summary>
    ///<param name="updates">Warping matrix collection</param>
    ///<param name="frameIndex">Index of starting frame in video</param>
    void initialize(const vt::vector<vt::CMtx3x3f>& updates, size_t frameIndex)
    {
        m_updates = updates;
        m_index = frameIndex;
    }

    ///<summary>
    /// Test whether the given frame is in the range of warping set.
    ///</summary>
    ///<param name="frameIndex">Frame index</param>
    bool is_in_range(size_t frameIndex) const
    {
        return (frameIndex >= m_index) && (frameIndex < m_index + m_updates.size());
    }

    ///<summary>
    /// Get warping matrix by frame index.
    ///</summary>
    ///<param name="frameIndex">Frame index</param>
    const vt::CMtx3x3f& operator[](size_t frameIndex) const
    {
        VT_ASSERT(is_in_range(frameIndex));
        return m_updates[frameIndex - m_index];
    }

    ///<summary>
    /// Get the size of warping matrix collection
    ///</summary>
    size_t size() const { return m_updates.size(); }

    ///<summary>
    /// Get the index of starting frame
    ///</summary>
    size_t start_index() const { return m_index; }

    ///<summary>
    /// Swap the content of input warping set and current set.
    ///</summary>
    ///<param name="other">Another warping set.</param>
    void swap(WarpsetWrapper& other)
    {
        m_updates.swap(other.m_updates);

        size_t tmp = m_index;
        m_index = other.m_index;
        other.m_index = tmp;
    }

private:
    // Disable copy constructor and assignment operation to avoid overheads in coyping matrix collection.
    WarpsetWrapper(const WarpsetWrapper&);
    WarpsetWrapper& operator=(const WarpsetWrapper&);

private:
    size_t m_index;
    vt::vector<vt::CMtx3x3f> m_updates;
};