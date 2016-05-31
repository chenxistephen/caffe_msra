//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      feature extraction and matching
//
//------------------------------------------------------------------------------
#pragma once

#include "vt_features_base.h"

#include "brief.h"
#include "overlapblock.h"

#define BRIEF_DESC_SIZE 128
#define BRIEF_PATCH_SIZE 24

// number of frames of track history used for inlier culling in RANSAC run
// after soft reset
#define SOFTRESET_RANSAC_TRACKHISTORY 4

//+-----------------------------------------------------------------------------
// 
//  Class to manage the feature tracks
// 
//------------------------------------------------------------------------------
struct FEATURE_POINT_MATCH
{ 
    CVec2f   vCur;
    CVec2f   vTrackBegin;
    CVec2f   vPrev;

    uint16_t uCurId;    // index of detected feature point for that frame
    uint16_t uTrackId;  // track id for this entry
    uint32_t uTrackLen;
    uint32_t uInlierCnt;
};

struct FRAME_MATCH_TABLE
{ 
    // matches are generated for each frame as detected feature point match
    // an existing track
    vt::vector<FEATURE_POINT_MATCH> vecMatches;

    // track map used to locate this frame's FPM (if any, else set to 0xffff)
    // for a given track id in the current track id space; the track id space
    // is the number of matched features when tracks are established at start
    // or reset and remains that size until the next reset
    vt::vector<uint16_t> vecTrackMap;

    FEATURE_POINT_MATCH& TrackIdToMatch(uint16_t trackid)
    {
        VT_ASSERT( trackid < vecTrackMap.size() &&
                   vecTrackMap[trackid] < vecMatches.size() );
        return vecMatches[vecTrackMap[trackid]];
    }
    
    const FEATURE_POINT_MATCH& TrackIdToMatch(uint16_t trackid) const
    {
        VT_ASSERT( trackid < vecTrackMap.size() &&
                   vecTrackMap[trackid] < vecMatches.size() );
        return vecMatches[vecTrackMap[trackid]];
    }
};

struct FRAME_MATCHES
{
    int iFrameId;
    int iRefFrameId;

    // number of inliers for this frame
    int iInliers;

    FRAME_MATCH_TABLE matchTables;
};

class CBriefFeatureTracker
{
public:
    CBriefFeatureTracker();

	HRESULT Begin(int iFrameWidth, int iFrameHeight, int iDetectDimPow2, 
				  int iBlockDimPow2, FEATURE_TRACKER_PARAMS* pPrms=NULL);

    HRESULT PushFrame(
        const vt::vector<HARRIS_FEATURE_POINT>& vecFP,
        const vt::vector<BriefDesc<BRIEF_DESC_SIZE>>& vecDesc);

	BUFFER_RANGE GetResultsRange();

    HRESULT GetResult(int frameNumber, vt::vector<PointMatch>& vecMatches, 
					  int& iRefFrameResult);

    HRESULT End(void);

    void SetCallback(unsigned int mask, 
        void (__stdcall *proc)(unsigned int type, const void* data, int frame));

protected:
    struct BEST_FEATURE
    {
        uint8_t  uDescDist;
        uint32_t uTrackLen;
        uint32_t uInlierCnt;
        uint16_t uPrevId;
    };

protected:
    HRESULT MatchTrackedFeatures(FRAME_MATCH_TABLE& fmt, int &iInlierCount, float& fAvgDist,
								 const vt::vector<HARRIS_FEATURE_POINT>& vecFP,
								 const vt::vector<BriefDesc<BRIEF_DESC_SIZE>>& vecDesc);
    HRESULT MatchAllFeatures(FRAME_MATCH_TABLE& fmt, 
							 const vt::vector<HARRIS_FEATURE_POINT>& vecFP,
							 const vt::vector<BriefDesc<BRIEF_DESC_SIZE>>& vecDesc);
    HRESULT MultiFrameRansac(bool& bLostFrames, int iOffset);
    HRESULT PostProcessInliers(int& iInlierCount);

protected:
    // current and previous frame of feature positions and descriptors
    // index 0 is current frame; index 1 is previous frame
    CRollingBuffer<vt::vector<CVec2f>> m_rbFrameFP;
    CRollingBuffer<vt::vector<BriefDesc<BRIEF_DESC_SIZE>>> m_rbFrameDesc;

    // sliding window of frame matches
    // index 0 is current frame; indices 1..n are previous frames
    CRollingBuffer<FRAME_MATCHES> m_rbFrameMatchTables;

    // averaged descriptors for the tracked features: there will be two
    // of these to enable propagating averages to new reference frame sequences;
    // index 0 is current track id sequence; index 1 is previous
    CRollingBuffer<vector<BriefAvgDesc<BRIEF_DESC_SIZE>>> m_rbTrackAvgDesc;

    // working vars 
    int m_iCurrentFrameId;
    // current (most recent) reference frame
    int m_iRefFrameId;
    // previous m_iRefFrameId
    int m_iPrevRefFrameId;
    // m_iRefFrameId at last hard reset
    int m_iLastHardResetFrameId;

    // number of tracked features at m_params.iMinTrackLength frames
    // after reset; used to force reset due to loss of tracked features
    int m_iBaseTrackedFeatures;
    // true if 'bufferless', in which case there is no buffer of video frames 
    // for the tracker, so the tracker works from a buffer of old frame info
    bool m_bBufferlessMode;
    // count of number of frames in buffer batch
    int m_iBufferedBatchCount;

    // params set at Begin
    FEATURE_TRACKER_PARAMS m_params;

    // datastructure to speed feature matching
    COverlappedBlocks<int> m_ovlblk;

    // results of feature matching for each new frame: these arrays
    // are the same size (the number of features in the current frame)
    // and hold data about the best match for the corresponding feature;
    // the Id is initially set to the track id of that feature during
    // MatchTrackedFeatures, then is set to the previous frame feature
    // id during MatchAllFeatures
    vt::vector<uint16_t>     m_vecBestFeatureId;
    // bestFeature is a cache of information for the best-matched track
    vt::vector<BEST_FEATURE> m_vecBestFeature;

    // map to enable access to track data across one reset boundary: indexed
    // by the feature id (uCurId) of the first frame after m_iRefFrameId
    // (i.e. the first frame in the current trackId space), and holds the
    // corresponding trackId (if it exists, else 0xffff) of the previous
    // trackId sequence
    vt::vector<uint16_t>    m_vecResetTrackIdMap;

    // current frame dimensions
    int m_iFrameWidth;
    int m_iFrameHeight;
    int m_iDetectDimensionPow2;

    // debug callback
    unsigned int m_callbackmask;
    void (__stdcall *m_callbackproc)(unsigned int type, const void* data, 
									 int frame);

    void FrameIdCallback(unsigned int type)
    {
        if (m_callbackmask & type)
        {
            (m_callbackproc)(type, NULL, m_iCurrentFrameId);
        }
    }
};

