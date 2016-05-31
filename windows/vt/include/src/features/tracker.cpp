//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      feature extraction and matching
//
//------------------------------------------------------------------------------
#include "stdafx.h"

using namespace vt;

#include "vt_stabilize.h"
#include "vt_harrisdetect.h"
#include "tracker.h"
#include "brief.h"
#include "overlapblock.h"

// define this for supporting similarity and affine motion models only
#define MOTION_MODEL_AFFINE 1 

#if !defined(MOTION_MODEL_AFFINE)

// include stitching to support more complex motion models
#define STITCHING_FEATURES_SUBSET

namespace vt
{
    #include "..\..\stitching\common\src\patchmatch.h"
    #include "..\..\stitching\common\src\vt_featurealign.h"
}
HRESULT AlignImagePair(CMtx3x3f& xform, SMatchTableV& matches, 
                       int iFrameWidth, int iFrameHeight, int motionModel);

#endif

// min inlier count for forcing resets
#define MIN_INLIER_THRESHOLD 8 

// size of buffer in tracker for 'bufferless' (i.e. no buffered video frames) 
// operation
#define BUFFERLESS_TRACKER_TRKBUF_SIZE 30

//+-----------------------------------------------------------------------------
// 
//  Class to manage the feature tracks
// 
//------------------------------------------------------------------------------
CBriefFeatureTracker::CBriefFeatureTracker() : m_iBaseTrackedFeatures(-1),
    m_callbackmask(0x0), m_callbackproc(NULL)
{}

HRESULT 
CBriefFeatureTracker::Begin(int iWidth, int iHeight, int iDetectDimPow2, 
                            int iBlockDimPow2, FEATURE_TRACKER_PARAMS* pPrms)
{
    VT_HR_BEGIN()

    m_params = pPrms? *pPrms: FEATURE_TRACKER_PARAMS();

    m_iCurrentFrameId = 0;
    m_iRefFrameId     = -1;
    m_iPrevRefFrameId = -1;
    m_iLastHardResetFrameId = 0;
    m_iBufferedBatchCount  = 0;    // initialize to force MatchAll on second frame

    m_bBufferlessMode = (m_params.iBufferSize == 0);

    int iTrkBufSize = m_bBufferlessMode? (BUFFERLESS_TRACKER_TRKBUF_SIZE): 
        (m_params.iBufferSize + SOFTRESET_RANSAC_TRACKHISTORY);

    // allocate for current and previous frame
    VT_HR_EXIT( m_rbFrameFP.resize(2) );
    VT_HR_EXIT( m_rbFrameDesc.resize(2) );

    // allocate fmt's and initialize frame id (since it is searched
    // on before the full array is set with frame data)
    VT_HR_EXIT( m_rbFrameMatchTables.resize(iTrkBufSize) );

    // allocate for current and previous track id space
    VT_HR_EXIT( m_rbTrackAvgDesc.resize(2) );

    float fMaxDist = float(VtMax(iWidth, iHeight)) * m_params.fMaxDistance;
    int maxDistPow2 = (int)ceil(logf(fMaxDist)/logf(2.f));
    VT_HR_EXIT( m_ovlblk.Create(CRect(0,0,iWidth,iHeight),
                                iBlockDimPow2, maxDistPow2) );

    m_iFrameWidth  = iWidth;
    m_iFrameHeight = iHeight;
    m_iDetectDimensionPow2 = iDetectDimPow2;

    VT_HR_END()
}

// do matching against currently tracked features (if any), and create
// new vecMatches array for this frame
HRESULT CBriefFeatureTracker::MatchTrackedFeatures(
    FRAME_MATCH_TABLE& fmt,
    int& iInlierCount,
    float& fAvgDist,
    const vt::vector<HARRIS_FEATURE_POINT>& vecFP,
    const vt::vector<BriefDesc<BRIEF_DESC_SIZE>>& vecDesc)
{
    VT_HR_BEGIN();

    iInlierCount = 0;  
    CVec2f avgOff(0,0);

    const FRAME_MATCHES& prevfmp      = m_rbFrameMatchTables[m_iCurrentFrameId-1];
    const FRAME_MATCH_TABLE& prevfmtp = prevfmp.matchTables;
    int iSearchCnt = int(prevfmtp.vecMatches.size());
    if (iSearchCnt <= 0) 
    { 
        return hr; 
    }

    // make sure id's line up as expected
    VT_ASSERT( prevfmp.iFrameId == m_iCurrentFrameId-1 );

    // reserve for every track matching
    VT_HR_EXIT( fmt.vecMatches.reserve(iSearchCnt) );
    fmt.vecMatches.resize(0);

    // allocate track id map of same size as previous and prefill with 'no match'
    VT_HR_EXIT( fmt.vecTrackMap.resize(prevfmtp.vecTrackMap.size()) );
    const uint16_t c_0xffff = 0xffff;
    VtFillSpan( fmt.vecTrackMap.begin(), &c_0xffff, sizeof(fmt.vecTrackMap[0]),
                int(prevfmtp.vecTrackMap.size()) );

    // for each track search for a match in the current set of features;
    // this is done in two passes because we want to cull out cases where 
    // multiple tracks match the same feature in the current frame and only
    // keep the closest match in those cases

    for( int i = 0; i < iSearchCnt; i++ )
    {
        // info for this track
        const FEATURE_POINT_MATCH& prevfpm = prevfmtp.vecMatches[i];
        const CVec2f& fpTrack = prevfpm.vCur;
        const BriefDesc<BRIEF_DESC_SIZE>& descTrack = 
            m_rbTrackAvgDesc.buffer(0)[prevfpm.uTrackId];

        int iNearest; // index to current frame's feature point array of nearest
        int iDistance; // distance of nearest
        BriefFindMatch(iNearest, iDistance, vecFP, m_ovlblk, vecDesc, fpTrack, 
                       descTrack, m_iFrameWidth, m_iFrameHeight, 
                       m_params.fMaxDistance, m_params.fMatchTestRatio);
        VT_ASSERT( iDistance < 256 );

        if( iNearest != -1 )
        {
            if( (m_vecBestFeatureId[iNearest] == 0xffff) ||
                (m_vecBestFeature[iNearest].uDescDist > (uint8_t)iDistance) )
            {
                // this feature point is closest (best) so far for this track
                m_vecBestFeatureId[iNearest] = prevfpm.uTrackId;
                BEST_FEATURE& bf = m_vecBestFeature[iNearest];
                bf.uDescDist = (uint8_t)iDistance;
                // use 'best feature' vector to persist track length and inlier
                // count, and make the previous frame position of this track
                // easily available for MatchAllFeatures
                // just copy the inlier count here it will be incremented or
                // zeroed out in one of the RANSAC routines
                bf.uTrackLen  = prevfmtp.vecMatches[i].uTrackLen+1;
                bf.uInlierCnt = prevfmtp.vecMatches[i].uInlierCnt;
                bf.uPrevId    = prevfpm.uCurId;
            }
        }
    }

    // second pass: for each feature, keep the closest match; count
    // matches that were inliers in the previous frame
    for( int i = 0; i < (int)m_vecBestFeatureId.size(); i++ )
    {
        uint16_t uTrackId = m_vecBestFeatureId[i];
        if( uTrackId != 0xffff )
        {
            const BEST_FEATURE& bf = m_vecBestFeature[i];

            // fpm for this track in the previous frame
            const FEATURE_POINT_MATCH& prevfpm = prevfmtp.TrackIdToMatch(uTrackId);

            // create new feature point match struct;
            FEATURE_POINT_MATCH fpm;
            fpm.uCurId      = (uint16_t)i; // index to feature point table
            fpm.vCur        = m_rbFrameFP.buffer(0)[fpm.uCurId]; // position in current frame
            fpm.vPrev       = prevfpm.vCur;
            fpm.uTrackLen   = bf.uTrackLen; // values already incremented
            fpm.uInlierCnt  = bf.uInlierCnt;
            fpm.uTrackId    = uTrackId;
            fpm.vTrackBegin = prevfpm.vTrackBegin;

            // entry in track map is index of this fpm about to be added to the array
            fmt.vecTrackMap[fpm.uTrackId] = (uint16_t)fmt.vecMatches.size();
            fmt.vecMatches.push_back(fpm);

            iInlierCount += (fpm.uInlierCnt!=0)? 1: 0;

            avgOff += fpm.vCur-fpm.vTrackBegin;

            // update the average descriptors for the tracks
            BriefAvgDesc<BRIEF_DESC_SIZE>& descTrack = 
                m_rbTrackAvgDesc.buffer(0)[uTrackId];
            descTrack.Acc(m_rbFrameDesc.buffer(0)[i]);
        }
    }

    fAvgDist = CVec2f(avgOff/float(fmt.vecMatches.size())).Magnitude();

    VT_HR_END();
}

// do matching against all features, and create new vecMatches for this frame;
HRESULT CBriefFeatureTracker::MatchAllFeatures(
    FRAME_MATCH_TABLE& fmt,
    const vt::vector<HARRIS_FEATURE_POINT>& vecFP,
    const vt::vector<BriefDesc<BRIEF_DESC_SIZE>>& vecDesc)
{
    VT_HR_BEGIN();

    int iNumFeatures = (int)m_vecBestFeatureId.size();

    // create an index of previous frame features that have already 
    // been matched, these will not be considered in the search below
    // since they already have a match using the avgDesc
    VT_HR_EXIT( m_vecResetTrackIdMap.resize(m_rbFrameFP.buffer(1).size()) );
  
    const uint16_t c_0xffff = 0xffff;
    VtFillSpan( m_vecResetTrackIdMap.begin(), &c_0xffff, 
                sizeof(m_vecResetTrackIdMap[0]), 
                int(m_vecResetTrackIdMap.size()) );

    if (fmt.vecTrackMap.size() > 0)
    {
        for( int i = 0; i < (int)m_vecBestFeatureId.size(); i++ )
        {
            uint16_t uPrevTrackId = m_vecBestFeatureId[i];
            if( uPrevTrackId != 0xffff )
            {
                const BEST_FEATURE& bf = m_vecBestFeature[i];
                m_vecResetTrackIdMap[bf.uPrevId] = uPrevTrackId;
            }
        }
    }

    // for each feature in the previous frame, search for a match in the
    // current frame's set of features;
    // this is done in two passes because we want to cull out cases where 
    // multiple features match the same feature in the current frame and only
    // keep the closest match in those cases;
    for( int i = 0; i < (int)m_rbFrameFP.buffer(1).size(); i++ )
    {        
        if( m_vecResetTrackIdMap[i] != 0xffff )
        {
            // don't search for previous matched feature points
            continue;
        }

        const CVec2f& fpTrack = m_rbFrameFP.buffer(1)[i];
        const BriefDesc<BRIEF_DESC_SIZE> descTrack = m_rbFrameDesc.buffer(1)[i];

        int iNearest;  // index to current frame feature point array of nearest
        int iDistance; // distance of nearest
        BriefFindMatch(iNearest, iDistance, vecFP, m_ovlblk, vecDesc, fpTrack, 
                       descTrack, m_iFrameWidth, m_iFrameHeight,
                       m_params.fMaxDistance, m_params.fMatchTestRatio);
        VT_ASSERT( iDistance < 256 );

        if( iNearest != -1 )
        {
            // for this case (matching all), MatchTrackedFeatures has always
            // already been called for this frame, so only update the best
            // feature entry if it has not already been set when matching tracks; 
            // the MatchTrackedFeatures matches are slighly better since they
            // compare to the average descriptors, so those matches are always
            // retained
            if( (m_vecBestFeatureId[iNearest] == 0xffff) ||
                ( (m_vecBestFeatureId[iNearest] == 0xfffe) &&
                  (m_vecBestFeature[iNearest].uDescDist > (uint8_t)iDistance) ) )
            {
                // this is a new closest match then update the match info
                m_vecBestFeatureId[iNearest] = 0xfffe; // indicates new track
                BEST_FEATURE& bf = m_vecBestFeature[iNearest];
                bf.uDescDist  = (uint8_t)iDistance;
                bf.uTrackLen  = 1;
                bf.uInlierCnt = 0;
                bf.uPrevId    = (uint16_t)i;
            }
        }
    }

    // swap the average descriptor index; will free the old one below
    m_rbTrackAvgDesc.advance();

    // reset the descriptor list (note: this generally over-reserves...)
    VT_HR_EXIT( fmt.vecMatches.reserve(iNumFeatures) );
    VT_HR_EXIT( fmt.vecTrackMap.reserve(iNumFeatures) );
    VT_HR_EXIT( m_rbTrackAvgDesc.buffer(0).reserve(iNumFeatures) );
    fmt.vecMatches.resize(0);
    fmt.vecTrackMap.resize(0);
    m_rbTrackAvgDesc.buffer(0).resize(0);

    // second pass: for each best feature match, generate a new track match entry
    for( int i = 0; i < iNumFeatures; i++ )
    {
        uint16_t uPrevTrackId = m_vecBestFeatureId[i];
        if( uPrevTrackId != 0xffff )
        {
            const BEST_FEATURE& bf = m_vecBestFeature[i];
            // for each track get the closest match out of the new features
            FEATURE_POINT_MATCH fpm;
            fpm.uCurId    = (uint16_t)i;
            fpm.vCur      = m_rbFrameFP.buffer(0)[i];
            fpm.vPrev     = m_rbFrameFP.buffer(1)[bf.uPrevId];
            // these have to come from BestFeature so they are persisted
            // from the previous reference frame sequence
            fpm.uTrackLen  = bf.uTrackLen;
            fpm.uInlierCnt = bf.uInlierCnt;

            // this is a new reference frame so the track id space is new,
            // so trackid simply increments for each new match
            fpm.uTrackId    = (uint16_t)fmt.vecMatches.size();
            fpm.vTrackBegin = m_rbFrameFP.buffer(1)[bf.uPrevId];
            fmt.vecTrackMap.push_back((uint16_t)fmt.vecMatches.size());
            fmt.vecMatches.push_back(fpm);

            // create and set the new average descriptor array entry
            m_rbTrackAvgDesc.buffer(0).resize(fpm.uTrackId+1);
            BriefAvgDesc<BRIEF_DESC_SIZE>& descTrack = 
                m_rbTrackAvgDesc.buffer(0)[fpm.uTrackId];
            if ( (m_rbTrackAvgDesc.buffer(1).size() > 0) && (uPrevTrackId != 0xfffe) )
            {
                // get the current average descriptor from the previous array
                descTrack = m_rbTrackAvgDesc.buffer(1)[uPrevTrackId];
                // NOTE: track average descriptors already include current frame
                //       because MatchTrackedFeatures was called
            }
            else
            {
                // new match (or second frame), so clear and set with current
                descTrack.Clear();
                descTrack.Acc( m_rbFrameDesc.buffer(0)[i] );
            }
        }
    }

    // update the reference frame-ids
    m_iPrevRefFrameId = m_iRefFrameId;
    m_iRefFrameId     = m_iCurrentFrameId-1;

    VT_HR_END();
}

HRESULT 
CBriefFeatureTracker::PushFrame(
    const vt::vector<HARRIS_FEATURE_POINT>& vecFP,
    const vt::vector<BriefDesc<BRIEF_DESC_SIZE>>& vecDesc)
{
    VT_HR_BEGIN()

    int iNumFeatures = (int)VtMin(vecFP.size(), size_t(0xffff));

    // copy the feature points and descriptors to internal storage
    VT_HR_EXIT( m_rbFrameFP.buffer(0).resize(iNumFeatures) );
    VT_HR_EXIT( m_rbFrameDesc.buffer(0).resize(iNumFeatures) );
    for(  int i = 0; i < iNumFeatures; i++ )
    {
        m_rbFrameFP.buffer(0)[i] = CVec2f(vecFP[i].x, vecFP[i].y);
    }
    VtMemcpy( m_rbFrameDesc.buffer(0).begin(), vecDesc.begin(), 
              sizeof(BriefDesc<BRIEF_DESC_SIZE>)*iNumFeatures );

    // debug info
    if (m_callbackmask & VSCB_DETECTEDFPT)
    {
        for( int i = 0; i < (int)vecFP.size(); i++ )
        {
            (m_callbackproc)(VSCB_DETECTEDFPT, (void*)&vecFP[i], m_iCurrentFrameId);
        }
    }

    // compute frames since last ref frame
    int iFramesSinceRef = m_iCurrentFrameId - m_iRefFrameId;

    // create a match table for this frame
    m_rbFrameMatchTables.advance();
    FRAME_MATCHES& fm  = m_rbFrameMatchTables.last();
    fm.iFrameId        = m_iCurrentFrameId;
    fm.iRefFrameId     = m_iRefFrameId;
    fm.iInliers        = -1;

    FRAME_MATCH_TABLE& fmt = fm.matchTables;
    fmt.vecMatches.resize(0);
    fmt.vecTrackMap.resize(0);

    // on first frame there is nothing to match against, so stop here
    if( m_iCurrentFrameId == 0 )
    {
        m_iCurrentFrameId++;
        m_rbFrameFP.advance();
        m_rbFrameDesc.advance();
        return S_OK;
    }

    // bin the current feature points to make the search faster
    m_ovlblk.Reset();
    for( int i = 0; i < iNumFeatures; i++ )
    { 
        const HARRIS_FEATURE_POINT& fp = vecFP[i];
        VT_HR_EXIT( m_ovlblk.AddElement(fp.x, fp.y, i) );
    }

    // set initial state of bestTrackId structure for this frame - 
    // 0xffff means no track
    VT_HR_EXIT( m_vecBestFeatureId.resize(iNumFeatures) );
    const uint16_t c_0xffff = 0xffff;
    VtFillSpan( m_vecBestFeatureId.begin(), &c_0xffff, 
                sizeof(m_vecBestFeatureId[0]), iNumFeatures);

    // reset the data-structure that holds the running match info for a feature
    VT_HR_EXIT( m_vecBestFeature.resize(iNumFeatures) );

    // unless this is the very first frame always try matching the tracked 
    // features
    bool bStartNewTracks = true;
    if( m_iCurrentFrameId != 1 )
    {
        int iFramesSinceNewTracks = m_iCurrentFrameId - m_iRefFrameId;

        int iInliersAfterMatch;
        float fAvgDist;
        VT_HR_EXIT( MatchTrackedFeatures(fmt, iInliersAfterMatch, fAvgDist, vecFP, 
                                         vecDesc) );

        // true if in bufferless mode and the tracked inliers fall below a
        // threshold
        bool bMinInliers = 
            m_bBufferlessMode && (iInliersAfterMatch < MIN_INLIER_THRESHOLD);

        // true if the number of tracked features has dropped below the 
        // threshold parameter
        bool bMinTracks = 
            ((int)fmt.vecMatches.size() < m_params.iMinTrackedFeatures);

        // true if we've moved too far
        bool bMaxMove = fabs(fAvgDist) > float(m_iFrameWidth/6);
        
        // set the baseline track count used for soft reset due to track loss
        if( iFramesSinceNewTracks == m_params.iMinTrackLength )
        {
            m_iBaseTrackedFeatures = int(fmt.vecMatches.size());
        }

        // true if tracker has reached the minimum track length since latest 
        // new tracks and the number of matches since hitting the min track 
        // length has dropped below the threshold parameter
        bool bLostTracks = (iFramesSinceNewTracks > m_params.iMinTrackLength) && 
            ( (float)fmt.vecMatches.size() < 
              ((float)m_iBaseTrackedFeatures * (1.f-m_params.fTrackLossFraction)) );

        // detect cases where new tracks need to to injected
        bStartNewTracks = bMinInliers || bMinTracks || bLostTracks || bMaxMove;

        if (bMinTracks || bLostTracks) { FrameIdCallback(VSCB_TRSOFTRESETT); }
        if (bMinInliers)               { FrameIdCallback(VSCB_TRSOFTRESETI); }
    }

    // run RANSAC for buffered mode
    if( !m_bBufferlessMode )
    {
        // increment number of frames in buffer note the 0th frame does not count
        if( !bStartNewTracks )
        {
            m_iBufferedBatchCount++;
        }

        // in buffered mode run RANSAC if (1) buffer is full, (2) new tracks
        // are about to be generated 
        if( (m_iBufferedBatchCount == m_params.iBufferSize) || bStartNewTracks )
        {
            // if MultiFrameRansac has too few inliers on any frame then the
            // first argument will return true causing the tracker to restart
            bool bLostFrames;
            VT_HR_EXIT( MultiFrameRansac(bLostFrames, bStartNewTracks? 1: 0) );
            m_iBufferedBatchCount = bStartNewTracks? 1: 0;
            bStartNewTracks = bLostFrames || bStartNewTracks;
        }
    }

    // add new tracks if necessary
restart_tracker:

    bool bTrackerReset = false;
    if( bStartNewTracks )
    {
        VT_HR_EXIT( MatchAllFeatures(fmt, vecFP, vecDesc) );

        if( (int)fmt.vecMatches.size() < m_params.iMinTrackedFeatures )
        {
            FrameIdCallback(VSCB_TRHARDRESET);
            m_iLastHardResetFrameId = m_iCurrentFrameId; 
            fm.iRefFrameId = -1; // mark the fm entry as having a hard reset

            bTrackerReset = true;
        }
        else
        {
            fm.iRefFrameId = m_iRefFrameId;
        }
    }

    // in bufferless mode we run some form of RANSAC (MFR/PostProc) 
    // every frame time
    if (m_bBufferlessMode && !bTrackerReset)
    {
        int iTrkBufSize = (int)m_rbFrameMatchTables.size();

        // determine if MultiFrameRANSAC needs to be run, either for
        // new reference frame, or for startup after hard reset; else
        // just run lighter-weight inlier post processing
        bool bMinInliers = false;
        if ( (iFramesSinceRef == 1) || 
             ((m_iCurrentFrameId - m_iLastHardResetFrameId) < iTrkBufSize) )
        {
            VT_HR_EXIT( MultiFrameRansac(bMinInliers, 0) );
        }
        else
        {
            int iInliersAfterRansac;
            VT_HR_EXIT( PostProcessInliers(iInliersAfterRansac) );
            bMinInliers = (iInliersAfterRansac < MIN_INLIER_THRESHOLD);
        }

        if( bMinInliers && !bStartNewTracks )
        {
            bStartNewTracks = true;
            goto restart_tracker;
        }

        // inliers are set for this frame for bufferless, so callback now
        if (m_callbackmask & VSCB_INLIERPT)
        {
            for( int i = 0; i < (int)fmt.vecMatches.size(); i++ )
            {
                const FEATURE_POINT_MATCH& fpm = fmt.vecMatches[i];
                if( fpm.uInlierCnt )
                {
                    (m_callbackproc)(VSCB_INLIERPT, (void*)&(fpm.vCur), m_iCurrentFrameId);
                }
            }
        }
    }

    // debug callbacks for tracked points
    if (m_callbackmask & VSCB_TRACKEDPT)
    {
        if( m_iCurrentFrameId == 1 )
        {
            // wait until frame 1 to dump tracks from frame 0
            for( int i = 0; i < (int)fmt.vecMatches.size(); i++ )
            {
                (m_callbackproc)(VSCB_TRACKEDPT, &(fmt.vecMatches[i].vTrackBegin), 0);
                (m_callbackproc)(VSCB_TRACKEDLEN, &(fmt.vecMatches[i].uTrackLen), 0);
            }
        }
        for( int i = 0; i < (int)fmt.vecMatches.size(); i++ )
        {
            (m_callbackproc)(VSCB_TRACKEDPT, &(fmt.vecMatches[i].vCur), m_iCurrentFrameId);
            (m_callbackproc)(VSCB_TRACKEDLEN, &(fmt.vecMatches[i].uTrackLen), m_iCurrentFrameId);
        }
    }

    // at the end bump the necessary pointers etc.
    m_iCurrentFrameId++;
    m_rbFrameFP.advance();
    m_rbFrameDesc.advance();

    VT_HR_END()
}

BUFFER_RANGE CBriefFeatureTracker::GetResultsRange()
{
    BUFFER_RANGE r;
    r.frame_count = m_rbFrameMatchTables.get_available_count();
    r.first_frame = m_rbFrameMatchTables.get_first_id();

    if( !m_bBufferlessMode )
    {
        r.frame_count -= m_iBufferedBatchCount;
        r.frame_count = VtMax(0, r.frame_count);
    }

    return r;
}

HRESULT CBriefFeatureTracker::GetResult(int frameNumber,
                                        vt::vector<PointMatch>& vecMatches, 
                                        int& iRefFrameResult)
{
    VT_HR_BEGIN();

    vecMatches.resize(0);

    // find rolling buffer entry that holds this frame number
    BUFFER_RANGE r = GetResultsRange();
    if( frameNumber < r.first_frame || frameNumber >= (r.first_frame+r.frame_count) )
    {
        VT_HR_EXIT( E_INVALIDARG );
    }

    VT_ASSERT( m_rbFrameMatchTables[frameNumber].iFrameId == frameNumber );

    FRAME_MATCHES& fm = m_rbFrameMatchTables[frameNumber];
    if( (frameNumber == 0) || (fm.iRefFrameId < 0) || (fm.iInliers <= 0) )
    {
        iRefFrameResult = -1;
    }
    else
    {
        iRefFrameResult = fm.iFrameId-1;

        // return the inliers
        const FRAME_MATCH_TABLE& fmt = fm.matchTables;
        VT_HR_EXIT( vecMatches.reserve(fmt.vecMatches.size()) );
        for ( int j = 0; j < (int)fmt.vecMatches.size(); j++ )
        {
            const FEATURE_POINT_MATCH& fpm = fmt.vecMatches[j];
#if _DEBUG
            // ensure no duplicate matches
            for( int k = 0; k < j; k++ )
            {
                const FEATURE_POINT_MATCH& fpmk = fmt.vecMatches[k];
                VT_ASSERT( fabs(fpmk.vPrev.x-fpm.vPrev.x) > .0001f ||
                           fabs(fpmk.vPrev.y-fpm.vPrev.y) > .0001f );
            }
#endif
            if ( fpm.uInlierCnt )
            {
                // use reference frame position for fixed camera, previous frame
                // otherwise
                PointMatch mp;
                mp.p0 = fpm.vCur;
                mp.p1 = fpm.vPrev;
                vecMatches.push_back(mp);
            }
        }
        VT_ASSERT( vecMatches.size() == (size_t)fm.iInliers );
        VT_ASSERT( vecMatches.size() >= 4 );
    }

    VT_HR_END();
}

HRESULT CBriefFeatureTracker::End(void)
{
    VT_HR_BEGIN();

    // called when user is not providing any more frames; is a nop
    // for bufferless mode; for buffered mode, need to generate
    // inliers for any buffered frames remaining
    if (!m_bBufferlessMode)
    {
        bool bLostFrames;
        VT_HR_EXIT( MultiFrameRansac(bLostFrames, 1) );
        m_iBufferedBatchCount = 0;
    }

    VT_HR_END();
}

bool MakeSimilarityFromMatchPairs(CMtx3x3d& mS, const CVec2f& p11, const CVec2f& p12,
                                  const CVec2f& p21, const CVec2f& p22)
{
    if( CVec2f(p11-p21).MagnitudeSq() < 4.f || 
        CVec2f(p12-p22).MagnitudeSq() < 4.f )
    {
        return false;
    }

    // extract the point matches and check that they aren't too close
    PointMatch pm[2];
    pm[0].p0 = p11;
    pm[0].p1 = p12;
    pm[1].p0 = p21;
    pm[1].p1 = p22;
    VtSimilarityFromPointMatches2D(mS, pm, 2);

    return true;
}

inline float ComputeMatchDist2(const CMtx3x3d& mS, const CVec2f& pt1, const CVec2f& pt2 )
{
    // TODO: could optimize for similarity since that is all we compute in RANSAC
    CVec3d vX = mS * CVec3d(pt1.x, pt1.y, 1.f);
    return float( CVec2d(pt2.x-vX.x, pt2.y-vX.y).MagnitudeSq() ); 
}

inline float WeightedInlierCount(int iCount, float fWeightedCount, int iTotalTracks)
{
    float b = float(iCount - 2*MIN_INLIER_THRESHOLD) / float(2*MIN_INLIER_THRESHOLD);
    b = VtMax(0.f, VtMin(1.f, b));
    return fWeightedCount*b + (1-b)*float(iCount);
}


template <typename T>
inline bool test_dist(const CVec2<T>& v1, const CVec2<T>& v2, T maxdist2)
{ return CVec2<T>(v1-v2).MagnitudeSq() > maxdist2; }

template <typename T>
bool IsMatrixIdent(const CMtx3x3<T>& m, const CRect& r)
{
    if( m[2][2]==0 )
    {
        return false;
    }

    T maxdist  = T(VtMax(r.Width(), r.Height())) * T(0.005);
    T maxdist2 = maxdist*maxdist;

    CVec3<T> a = m*CVec3<T>(T(r.left), T(r.top), T(1));
    if( a.z==0 || test_dist(a.Dehom(), CVec2<T>(T(r.left), T(r.top)), maxdist2) )
        return false;

    CVec3<T> b = m*CVec3<T>(T(r.right), T(r.top), T(1));
    if( b.z==0 || test_dist(b.Dehom(), CVec2<T>(T(r.right), T(r.top)), maxdist2) )
        return false;

    CVec3<T> c = m*CVec3<T>(T(r.left), T(r.bottom), T(1));
    if( c.z==0 || test_dist(c.Dehom(), CVec2<T>(T(r.left), T(r.bottom)), maxdist2) ) 
        return false;

    CVec3<T> d = m*CVec3<T>(T(r.right), T(r.bottom), T(1));
    if( d.z==0 || test_dist(d.Dehom(), CVec2<T>(T(r.right), T(r.bottom)), maxdist2) )
        return false;

    return true;
}

struct INL_INFO
{
    uint16_t uTrackId;
    uint32_t uInlierCnt;
};

HRESULT 
CBriefFeatureTracker::MultiFrameRansac(bool& bLostFrames, int iOffset)
{
    bLostFrames = false;

    // length of track buffer 
    int iTrkBufSize = (int)m_rbFrameMatchTables.size();
    VT_ASSERT( iTrkBufSize > 1 );

    // index into frame match table rolling buffer
    int iFrameId = m_iCurrentFrameId - iOffset;
 
    int iFramesSinceRef = iFrameId - m_iRefFrameId;

    // scale outlier distance parameter by the scale factor between the (always
    // decimated) coordinate space where detection is done and the coordinate
    // space where RANSAC processing is done (the latter is the same as the image
    // source pixel coordinates)
    float fODScale = (float)(1<<m_iDetectDimensionPow2);
    float fOutlierDist2 = (m_params.fRANSACOutlierDist*fODScale)*
                          (m_params.fRANSACOutlierDist*fODScale);

    VT_HR_BEGIN()

    // set number of (non-reset) frames on which to run RANSAC
    int iNumVerifyFrames, iNumResultFrames;
    if( m_bBufferlessMode )
    {
        iNumVerifyFrames = iNumResultFrames = VtMin(iFramesSinceRef, iTrkBufSize);
    }
    else
    {
        iNumVerifyFrames = VtMin(m_params.iBufferSize, iFramesSinceRef);
        iNumResultFrames = VtMin(m_iBufferedBatchCount, iFramesSinceRef);
    }

    if( iNumResultFrames <= 0 )
    {
        return S_OK;
    }

    int iNumOrigFrames = iNumResultFrames;

    vt::vector<uint16_t> vecInliers1, vecInliers2;
    vt::vector<uint16_t> *pCurInliers = &vecInliers1;
    vt::vector<uint16_t> *pMaxInliers = &vecInliers2;
    //float fMaxIdent = 1.f;

    // inlier count for winning RANSAC run
    int iMaxInliersCount = -1;

    bool bDone = false;
    while ( !bDone )
    {
        FRAME_MATCH_TABLE& matchEnd = m_rbFrameMatchTables[iFrameId].matchTables;
        int iFinalTrackCount = int(matchEnd.vecMatches.size());
    
        VT_HR_EXIT( vecInliers1.reserve(iFinalTrackCount) );
        VT_HR_EXIT( vecInliers2.reserve(iFinalTrackCount) );

        // weight for winning RANSAC run
        float fMaxWeight = 0;

        CRand rnd;
        rnd.Seed(9037);
        // run RANSAC iterations to select best set of inliers to transform
        // from the reference frame to the current frame
        bool bFirstRun = true; // flag to control clearing of incoming inlier flags
        for(int i = 0; i < m_params.iRANSACIterations; i++)
        {
            // pick a random set of matches from begin and end of track
            int iMatch1 = rnd.IRand(iFinalTrackCount);
            int iMatch2 = rnd.IRand(iFinalTrackCount);
            if( iMatch1 == iMatch2 )
            {
                continue;
            }

            // compute a similarity transform
            const FEATURE_POINT_MATCH& fpmEnd1 = matchEnd.vecMatches[iMatch1];
            const FEATURE_POINT_MATCH& fpmEnd2 = matchEnd.vecMatches[iMatch2];

            CMtx3x3d mS;
            if( !MakeSimilarityFromMatchPairs(mS, fpmEnd1.vTrackBegin, fpmEnd1.vCur,
                                              fpmEnd2.vTrackBegin, fpmEnd2.vCur) )
            {
                continue;
            }

            float fIdentBonus = 
                IsMatrixIdent(mS, CRect(0,0,m_iFrameWidth,m_iFrameHeight))? 1.4f: 1.f;

            // compute the inliers for the current frame
            int iCurInliers = 0;
            float fWeightedInliers = 0;

            pCurInliers->resize(0);
            for( int j = 0; j < iFinalTrackCount; j++ )
            {
                FEATURE_POINT_MATCH& fpmEnd = matchEnd.vecMatches[j];

                if( ComputeMatchDist2(mS, fpmEnd.vTrackBegin, fpmEnd.vCur) < fOutlierDist2 )
                {   
                    iCurInliers++;
                    VT_ASSERT( fpmEnd.uTrackLen != 0 );

                    // weight inliers higher when they have been inliers in 
                    // previous frames
                    fWeightedInliers += fIdentBonus*float(1+fpmEnd.uInlierCnt);
                    pCurInliers->push_back((uint16_t)j);
                }
            }

            float fCurW = WeightedInlierCount(iCurInliers, fWeightedInliers, 
                                              iFinalTrackCount);
            if( fCurW <= fMaxWeight )
            {
                // test for early out
                goto EndRansacIter;
            }

            // test the proposed set of inliers against previous frames: 
            // compute the previous-frame similarity for the random match pair 
            // and compare the previous-frame inliers to them, and remove the 
            // inliers and their weighting contribution if outside the distance
            for( int j = 1; j < iNumVerifyFrames; j++ )
            {
                // get FM and FMT for current frame
                const FRAME_MATCHES& fm = m_rbFrameMatchTables[iFrameId-j];
                if ( fm.matchTables.vecMatches.size() == 0 )
                {
                    continue;
                }
                VT_ASSERT( fm.matchTables.vecTrackMap.size() == matchEnd.vecTrackMap.size() );
                const FRAME_MATCH_TABLE& matchCur = fm.matchTables;
                const FEATURE_POINT_MATCH& fpmCur1 = 
                    matchCur.TrackIdToMatch(fpmEnd1.uTrackId);
                const FEATURE_POINT_MATCH& fpmCur2 = 
                    matchCur.TrackIdToMatch(fpmEnd2.uTrackId);

                // compute a similarity for each of the other sets
                if( !MakeSimilarityFromMatchPairs(mS, fpmCur1.vTrackBegin, fpmCur1.vCur,
                                                  fpmCur2.vTrackBegin, fpmCur2.vCur) )
                {
                    goto EndRansacIter;
                }

                for ( int k = 0; k < (int)pCurInliers->size(); k++ )
                {
                    uint16_t inlierid = (*pCurInliers)[k];
                    if (inlierid == 0xffff) 
                    {
                        continue; // track already removed
                    }
                    uint16_t trackId = matchEnd.vecMatches[inlierid].uTrackId;

                    const FEATURE_POINT_MATCH& fpmCur = 
                        matchCur.TrackIdToMatch(trackId);
        
                    if( ComputeMatchDist2(mS, fpmCur.vTrackBegin, fpmCur.vCur) >= fOutlierDist2 )
                    {   
                        // RANSAC-selected inlier for current frame is not 
                        // inlier in this previous frame

                        // 'remove' track from inlier list
                        iCurInliers--;
                        (*pCurInliers)[k] = 0xffff; 

                        // remove the weighting contribution this inlier had made
                        fWeightedInliers -= fIdentBonus * 
                            float(1+matchEnd.TrackIdToMatch(trackId).uInlierCnt);
                        fCurW = WeightedInlierCount(
                            iCurInliers, fWeightedInliers, iFinalTrackCount);
                        if( fCurW <= fMaxWeight )
                        {
                            // early out test
                            goto EndRansacIter;
                        }
                    }
                }
            }

            // due to early out above we should never get here without an new max
            VT_ASSERT( fCurW > fMaxWeight );

            fMaxWeight = fCurW;
            iMaxInliersCount = iCurInliers;
            vt::vector<uint16_t> *pTmp = pMaxInliers;
            pMaxInliers = pCurInliers;
            pCurInliers = pTmp;
            //fMaxIdent = fIdentBonus;

            // if every tracked point is an inlier then quit because it
            // can't get any better
            if( iCurInliers == iFinalTrackCount )
            {
                break;
            }
        EndRansacIter:;
            bFirstRun = false;
        }

        // if the number of inliers is below the min threshold, then give up on 
        // this frame and try to restart RANSAC with earlier frames
        bDone = true;
        if( (!m_bBufferlessMode) && (iMaxInliersCount < MIN_INLIER_THRESHOLD) )
        {
            // mark frame as "lost" by setting iRefFrameId to -1
            m_rbFrameMatchTables[iFrameId].iRefFrameId = -1;
            m_rbFrameMatchTables[iFrameId].iInliers    = 0;

            bLostFrames = true;

            // dumped this frame, now try starting one frame earlier
            if( iNumVerifyFrames <= iNumResultFrames )
            {
                iNumResultFrames--;
                iFrameId--;
            }
            iNumVerifyFrames--;

            // TODO: would be good to flag these as HARDRESET frames in the 
            //       debug output as well

            bDone = (iNumResultFrames == 0);
        }
    }


    // about to reset the inlier counts so zero them out in all the images
    // and record the inliercnts for the last tracked image 
    vt::vector<INL_INFO> vecInputInlierTrackIds;
    for( int i=0; i<iNumOrigFrames; i++ )
    {
        int k = m_iCurrentFrameId - iOffset - i;
        FRAME_MATCH_TABLE& fmt = m_rbFrameMatchTables[k].matchTables;
        if( k==iFrameId )
        {
            VT_HR_EXIT( vecInputInlierTrackIds.resize(fmt.vecMatches.size()) );
            for (int j=0; j<(int)vecInputInlierTrackIds.size(); j++)
            {
                FEATURE_POINT_MATCH& fpmEnd = fmt.vecMatches[j];
                
                INL_INFO& ii = vecInputInlierTrackIds[j];
                ii.uTrackId   = fpmEnd.uTrackId;
                ii.uInlierCnt = fpmEnd.uInlierCnt;
                fpmEnd.uInlierCnt = 0; // zero out after caching
            }
        }
        else
        {
            for (int j=0; j<(int)fmt.vecMatches.size(); j++)
            {
                fmt.vecMatches[j].uInlierCnt = 0;
            }
        }
    }

    // mark the new inliers, and set the number of inliers in the frame match 
    // table entries
    int framesToMark = m_bBufferlessMode? 1: iNumResultFrames; 
    for (int j=0; j<framesToMark; j++)
    {
        FRAME_MATCHES& fm = m_rbFrameMatchTables[iFrameId - j];

        VT_ASSERT( fm.matchTables.vecMatches.size() != 0 );

        int avg = 0;
        int avg_cnt = 0;
        for ( int i = 0; i < (int)pMaxInliers->size(); i++ )
        {
            uint16_t inlierid = (*pMaxInliers)[i];
            if ( inlierid == 0xffff) { continue; }

            const INL_INFO& ii = vecInputInlierTrackIds[inlierid];
            FEATURE_POINT_MATCH& fpm = fm.matchTables.TrackIdToMatch(ii.uTrackId);   
            fpm.uInlierCnt = ii.uInlierCnt + (framesToMark-j);

            if (!m_bBufferlessMode && (m_callbackmask & VSCB_INLIERPT))
            {
                (m_callbackproc)(VSCB_INLIERPT, (void*)&(fpm.vCur), fm.iFrameId);
                // this path is not entered for the first frame, so always return
                // a value for the first while processing the second
                if (1==fm.iFrameId)
                {
                    (m_callbackproc)(VSCB_INLIERPT, (void*)&(fpm.vCur), 0);
                }
            }
            avg += ii.uInlierCnt;
            avg_cnt++;
        }
        
        // debug only stuff
        /*
        if( j == 0)
        {
            printf("  max ident bonus: %f\n",fMaxIdent);
            printf("  avg inlier count: RANSAC tracks %f\n",float(avg)/float(avg_cnt));

            avg = avg_cnt = 0;
            for( int i = 0; i < (int)vecInputInlierTrackIds.size(); i++ )
            {
                const INL_INFO& ii = vecInputInlierTrackIds[i];
                avg+=ii.uInlierCnt;
                avg_cnt++;
            }
            printf("                      incoming tracks %f\n",float(avg)/float(avg_cnt));
        }
        */

        fm.iInliers = iMaxInliersCount;
    }

    VT_HR_EXIT_LABEL()

    return hr;
}

HRESULT 
CBriefFeatureTracker::PostProcessInliers(int& iInlierCount)
{
    VT_HR_BEGIN()

    // run RANSAC on the current set of inliers to cull away
    // any that are not tracking well
    iInlierCount = 0;

    // scale outlier distance parameter by the scale factor between the (always
    // decimated) coordinate space where detection is done and the coordinate
    // space where RANSAC processing is done (the latter is the same as the image
    // source pixel coordinates)
    float fODScale = (float)(1<<m_iDetectDimensionPow2);
    float fOutlierDist2 = (m_params.fRANSACOutlierDist*fODScale)*
                          (m_params.fRANSACOutlierDist*fODScale);

    FRAME_MATCHES& fm = m_rbFrameMatchTables.last();
    FRAME_MATCH_TABLE& fmt = fm.matchTables;
    int iMatchCount = int(fmt.vecMatches.size());

    // form vector of track id's of incoming inliers, and clear the inlier flags
    // at the same time
    // TODO: this may be generally useful and could be part of the tracker 
    // state; can clear the incoming inlier flags during RANSAC just like MFR
    vt::vector<INL_INFO> vecInputInlierTrackIds;
    VT_HR_EXIT( vecInputInlierTrackIds.reserve(iMatchCount) );
    for (int j=0; j<iMatchCount; j++)
    {
        FEATURE_POINT_MATCH& fpm = fmt.vecMatches[j];
        if (fpm.uInlierCnt) 
        { 
            INL_INFO ii;
            ii.uTrackId   = fpm.uTrackId;
            ii.uInlierCnt = fpm.uInlierCnt;
            vecInputInlierTrackIds.push_back(ii);
            fpm.uInlierCnt = 0;
        }
    }

    int iInputInlierCount = (int)vecInputInlierTrackIds.size();

    // two vectors of track ids; one to hold current proposed
    // set of inliers and one to hold current highest weight set
    vt::vector<uint16_t> vecInliers1, vecInliers2;
    VT_HR_EXIT( vecInliers1.reserve(iInputInlierCount) );
    VT_HR_EXIT( vecInliers2.reserve(iInputInlierCount) );
    vt::vector<uint16_t> *pCurInliers = &vecInliers1;
    vt::vector<uint16_t> *pMaxInliers = &vecInliers2;
    float fMaxWeight = 0;

    CRand rnd;
    rnd.Seed(5372);
    // run RANSAC iterations to select best set of inliers; run no more than
    // 2*inlier count iterations
    for (int i = 0; i < VtMin(m_params.iRANSACIterations,2*iInputInlierCount); i++)
    {
        // pick a random set of matches from begin and end of track
        uint16_t iMatch1TrackId = vecInputInlierTrackIds[rnd.IRand(iInputInlierCount)].uTrackId;
        uint16_t iMatch2TrackId = vecInputInlierTrackIds[rnd.IRand(iInputInlierCount)].uTrackId;
        if( iMatch1TrackId == iMatch2TrackId )
        {
            continue;
        }

        // compute a similarity transform
        CMtx3x3d mS;
        const FEATURE_POINT_MATCH& fpmEnd1 = fmt.TrackIdToMatch(iMatch1TrackId);
        const FEATURE_POINT_MATCH& fpmEnd2 = fmt.TrackIdToMatch(iMatch2TrackId);
        if( !MakeSimilarityFromMatchPairs(mS, fpmEnd1.vTrackBegin, fpmEnd1.vCur,
                                          fpmEnd2.vTrackBegin, fpmEnd2.vCur) )
        {
            continue;
        }

        // test the current set of inliers
        int iCurInlierCount = 0;
        int iWeightedInliers = 0;
        pCurInliers->resize(0);
        for( int j = 0; j < iInputInlierCount; j++ )
        {
            FEATURE_POINT_MATCH& fpmEnd = fmt.TrackIdToMatch(vecInputInlierTrackIds[j].uTrackId);
            if( ComputeMatchDist2(mS, fpmEnd.vTrackBegin, fpmEnd.vCur) < fOutlierDist2 )
            {   
                iCurInlierCount++;
                // weight inliers higher when they have been inliers in previous frames
                iWeightedInliers += (1+fpmEnd.uInlierCnt);
                pCurInliers->push_back((uint16_t)j);
            }
        }
        float fCurW = (float)iWeightedInliers;
        if( fCurW <= fMaxWeight )
        {
            // test for early out
            goto EndRansacIter;
        }
        // due to early out above we should never get here without an new max
        VT_ASSERT( fCurW > fMaxWeight );

        fMaxWeight = fCurW;
        vt::vector<uint16_t> *pTmp = pMaxInliers;
        pMaxInliers = pCurInliers;
        pCurInliers = pTmp;

        // if every input inlier is an output inlier then quit because it
        // can't get any better
        if( iCurInlierCount == iInputInlierCount )
        {
            break;
        }
    EndRansacIter:;
    }

    // set the output count and mark the new set of inliers
    iInlierCount = (int)pMaxInliers->size();
    for ( int j = 0; j < iInlierCount; j++ )
    {
        const INL_INFO& ii = vecInputInlierTrackIds[(*pMaxInliers)[j]];
        fmt.TrackIdToMatch(ii.uTrackId).uInlierCnt = ii.uInlierCnt+1;
    }

    fm.iInliers = iInlierCount;
       
    VT_HR_END()
}

void
CBriefFeatureTracker::SetCallback(
    unsigned int mask, void (__stdcall *proc)(unsigned int type, const void* data, int frame))
{
    m_callbackmask = mask;
    m_callbackproc = proc;
}

//+-----------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
#if !defined(MOTION_MODEL_AFFINE)
HRESULT AlignImagePair(CMtx3x3f& xform, SMatchTableV& matches, 
                       int iFrameWidth, int iFrameHeight, int motionModel)
{
    VT_HR_BEGIN()

    ALIGN_IMAGE_INFO_VECTOR vecAlignerInfo;
    ALIGN_IMAGE_DATA_VECTOR vecAlignerData;
    VT_HR_EXIT(vecAlignerInfo.resize(2));
    VT_HR_EXIT(vecAlignerData.resize(2));

    for( int i = 0; i < 2; i++ )
    {
        vecAlignerInfo[i].fExifFocalLen = -1;
        vecAlignerInfo[i].eOrientExif   = eUnknown;
        vecAlignerInfo[i].iWidth        = iFrameWidth;
        vecAlignerInfo[i].iHeight       = iFrameHeight;

        // switch on motion type
        switch(motionModel)
        {
        case eRigidScale:
            vecAlignerData[i].pose.InitRigidScale(0,1,0,0);
            break;
        case e3DRot:
            // note focal length hard-coded to 2.4 or FOV of 45degrees
            // TODO: use known camera focal
            vecAlignerData[i].pose.Init3D(CQuaterniond().MakeI(), CVec3d(0,0,0), 2.4);
            break;
        case eAffine:
        case eHomography:
            {
                CMtx3x3f ident;
                ident.MakeI();
                vecAlignerData[i].pose.InitTransform3((eMotionModel)motionModel, ident);
            }
            break;
        default:
            VT_ASSERT(0);
            break;
        } 

        vecAlignerData[i].pInfo = &vecAlignerInfo[i];
    }

    // params no radial, fix motion, no exif
    SGlobalAlignParam alignerParams;
    alignerParams.radial_distort_correct = 0;
    alignerParams.eMotion = (eMotionModel)motionModel;
    alignerParams.use_exif_focal_if_available = false;

    // bundle adjustment
    double rms_error;
    IterateGlobalAlign(vecAlignerData, matches, alignerParams, 0, false, true, 
                       rms_error, false, true, NULL);
    IterateGlobalAlign(vecAlignerData, matches, alignerParams, 0, false, true, 
                       rms_error, true, true, NULL);

    xform = TransScaleMat3(iFrameWidth, iFrameHeight).Inv() *
        vecAlignerData[0].pose.M3().Inv() * vecAlignerData[1].pose.M3() * 
        TransScaleMat3(iFrameWidth, iFrameHeight);

    VT_HR_END()
}

#endif