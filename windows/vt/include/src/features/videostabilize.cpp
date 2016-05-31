//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      video stabilization API
//
//------------------------------------------------------------------------------
#include "stdafx.h"

using namespace vt;

#include "vt_stabilize.h"
#include "stabilizeint.h"
#include "brief.h"
#include "tracker.h"
#include "vt_harrisdetect.h"
#include "FAST10.h"
#include "rollingshutter.h"
#include "vt_smoother_dp.h"
#include "vt_smoother_wls.h"
#include "transforms.h"

//#define INLIERPOINT_FILE

//------------------------------------------------------------------------------
// debug utilities
//------------------------------------------------------------------------------
#ifdef INLIERPOINT_FILE
FILE *g_fpFeaturesFile = NULL;
#endif

// from stackstabilize.cpp
extern HRESULT 
RunDetectorAndExtractDescriptors(vt::vector<HARRIS_FEATURE_POINT>& fpnts, 
             vt::vector<BriefDesc<BRIEF_DESC_SIZE>>& desc,
             const vt::CLumaByteImg& img,
             const FEATURE_DETECTOR_PARAMS& detectParams, int iExtractLevel,
             const BriefTable<BRIEF_DESC_SIZE>& briefTable,
             float *pBlurLevel,
             CHarrisDetector& hdetector);

//------------------------------------------------------------------------------
// warp matches
//------------------------------------------------------------------------------
inline CVec2f InterpFP(const CVec2f& fp, const CVec2f* pC)
{
    float f0 = floor(fp.y);
    float c0 = f0-fp.y;
    int   i0 = (int)f0;
    return fp-(pC[i0]*(1.f-c0)+pC[i0+1]*c0);
}

HRESULT WarpMatches(vt::vector<PointMatch>& result, 
                    const vt::vector<PointMatch>& src, 
                    const CVec2f* pfCorrection0, const CVec2f* pfCorrection1, 
                    int iHeight)
{
    VT_HR_BEGIN()

    VT_HR_EXIT( result.resize(src.size()) );

    for( int i = 0 ; i < (int)src.size(); i++ )
    {
#if 1
        const PointMatch& s = src[i];
        PointMatch& r = result[i];

        VT_ASSERT( s.p0.y >= float(0) && s.p0.y <= float(iHeight-2) );
        VT_ASSERT( s.p1.y >= float(0) && s.p1.y <= float(iHeight-2) );

        r.p0 = InterpFP(s.p0, pfCorrection0);
        r.p1 = InterpFP(s.p1, pfCorrection1);
#else
        result[i] = src[i];
#endif
    }

    VT_HR_END()
}

//------------------------------------------------------------------------------
// internal class for video stabilizer
//------------------------------------------------------------------------------
class CVideoStabilizer::CVideoStabilizerInternal
{
public:
    CVideoStabilizerInternal()
        :m_bStarted(false)
        ,m_pushedFrameCount(0),m_trackerNextFrame(0)
        ,m_frameWidth(0), m_frameHeight(0)
        ,m_pSmoother(NULL)
        ,m_detectDimensionSteps(0)
        ,m_inBlurReset(false)
        ,m_callbackmask(0x0),m_callbackproc(NULL)
    {
        // create the Brief table 
        m_bt.Initialize(BRIEF_PATCH_SIZE);
    };
    ~CVideoStabilizerInternal();

    // internal methods forwarded from CStabilize methods of the same name
    int GetMaxDelay();
    BUFFER_RANGE GetResultsRange() const;

    HRESULT Start(const VIDEO_STABILIZER_PARAMS& params);
    HRESULT GetParams(VIDEO_STABILIZER_PARAMS& params) const;
    HRESULT PreProcessFrame(vt::CLumaByteImg& img, CStabilizerFrameData& data);
    HRESULT PushFrame(CStabilizerFrameData& data);
    HRESULT PushFrame(vt::CLumaByteImg& img);
    HRESULT GetResult(vt::CMtx3x3f& transform, vt::vector<CVec2f>* pvecRScorrect, 
                      int frameNumber) const;
    HRESULT End(void);
    void SetCallback(unsigned int mask, 
                     void (__stdcall *proc)(unsigned int type, 
                                            const void* data, int frame));

protected:
    HRESULT AdvancePipeline();

    bool m_bStarted;
    int  m_pushedFrameCount;
    int  m_trackerNextFrame;
    int  m_frameWidth;
    int  m_frameHeight;
    int  m_detectDimensionSteps;

    VIDEO_STABILIZER_PARAMS m_params;

    BriefTable<BRIEF_DESC_SIZE> m_bt; 
    CHarrisDetector      m_hdetector;

    CBriefFeatureTracker m_tracker;

    RSC m_rollingShutterCorrect;

    // temp buffer to hold the tracker feature points
    vt::vector<PointMatch> m_vecMatches;

    // interface to the smoother for modes other than '0'
    // TODO: implement mode '0' smoother
    IFeatureWarpCompute  *m_pSmoother;

    // TODO: move this to a simple stabmode==0 impl of IFeatureWarpCompute
    // current chained reference transform 
    // m_mtxPrev is the transformation (reported by the tracker) of the 
    // previous frame 
    CMtx3x3f m_mtxChained;
    CMtx3x3f m_mtxPrev;

    CRollingBuffer<vt::vector<CVec2f>> m_bufRollingShutterCorrect;

    CRollingBuffer<CMtx3x3f> m_bufResultTransforms;

    // true if stabilizer is not tracking due to blur level
    bool    m_inBlurReset;

    // debug members
    unsigned int m_callbackmask;
    void (__stdcall *m_callbackproc)(unsigned int type, const void* data, int frame);
};

CVideoStabilizer::CVideoStabilizer(): m_pint(new CVideoStabilizerInternal()) 
{}

CVideoStabilizer::~CVideoStabilizer() 
{ delete m_pint; }

int CVideoStabilizer::GetMaxDelay()
{ return (m_pint==NULL)? 0: m_pint->GetMaxDelay(); }

BUFFER_RANGE CVideoStabilizer::GetResultsRange()
{ return (m_pint==NULL)? BUFFER_RANGE(): m_pint->GetResultsRange(); }

HRESULT CVideoStabilizer::GetParams(VIDEO_STABILIZER_PARAMS& params)
{ return (m_pint==NULL)? E_OUTOFMEMORY: m_pint->GetParams(params); }

HRESULT CVideoStabilizer::Start(const VIDEO_STABILIZER_PARAMS& params)
{ return (m_pint==NULL)? E_OUTOFMEMORY: m_pint->Start(params); }

HRESULT CVideoStabilizer::End(void)
{ return (m_pint==NULL)? E_OUTOFMEMORY: m_pint->End(); }

HRESULT
CVideoStabilizer::PreProcessFrame(vt::CLumaByteImg& img, CStabilizerFrameData& data)
{ return (m_pint==NULL)? E_OUTOFMEMORY: m_pint->PreProcessFrame(img,data); }

HRESULT
CVideoStabilizer::PushFrame(CStabilizerFrameData& data)
{ return (m_pint==NULL)? E_OUTOFMEMORY: m_pint->PushFrame(data); }

HRESULT
CVideoStabilizer::PushFrame(vt::CLumaByteImg& img)
{ return (m_pint==NULL)? E_OUTOFMEMORY: m_pint->PushFrame(img); }

HRESULT 
CVideoStabilizer::GetResult(vt::CMtx3x3f& transform, int frameNumber)
{
    return (m_pint==NULL)? 
        E_OUTOFMEMORY: m_pint->GetResult(transform, NULL, frameNumber);
}

HRESULT 
CVideoStabilizer::GetResult(vt::CMtx3x3f& transform,
                            vt::vector<CVec2f>& vecRScorrect, int frameNumber)
{
    return (m_pint==NULL)? 
        E_OUTOFMEMORY: m_pint->GetResult(transform, &vecRScorrect, frameNumber);
}

void 
CVideoStabilizer::SetCallback(
    unsigned int mask, void (__stdcall *proc)(unsigned int type, const void* data, int frame))
{
    if( m_pint )
    {
        m_pint->SetCallback(mask, proc);
    }
}

//+-----------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
CVideoStabilizer::CVideoStabilizerInternal::~CVideoStabilizerInternal()
{
    if( m_pSmoother )
    {
        delete m_pSmoother;
        m_pSmoother = NULL;
    }
}

int CVideoStabilizer::CVideoStabilizerInternal::GetMaxDelay()
{ 
    if( m_bStarted )
    {
        // TODO: call the base objects instead of reading their params
        // TODO: add optional rolling shutter correction
        // 3 for rolling shuttuer
        return m_params.stabilizeBuffer + 4 + m_params.trackParams.iBufferSize;
        //return m_params.stabilizeBuffer + m_params.trackParams.iBufferSize;
    }
    else
    {
        return 0;
    }
}

BUFFER_RANGE CVideoStabilizer::CVideoStabilizerInternal::GetResultsRange() const
{
    BUFFER_RANGE r;
    r.frame_count = m_bufResultTransforms.get_available_count();
    r.first_frame = m_bufResultTransforms.get_first_id();
    return r;
}

HRESULT
CVideoStabilizer::CVideoStabilizerInternal::GetParams(VIDEO_STABILIZER_PARAMS& params) const
{
    params = m_params;
    return S_OK;
}

HRESULT
CVideoStabilizer::CVideoStabilizerInternal::Start(
    const VIDEO_STABILIZER_PARAMS& params)
{
    VT_HR_BEGIN();
        
    m_bStarted = true;

    m_params = params;

    m_pushedFrameCount = 0;
    m_trackerNextFrame = 0;

    // reset the chained transform
    m_mtxChained.MakeI();
    m_mtxPrev.MakeI();

    // need to store the rolling buffer correction through the stabilize pipeline
    // since it is returned to the caller
    VT_HR_EXIT( m_bufRollingShutterCorrect.resize(GetMaxDelay()) );
    VT_HR_EXIT( m_bufResultTransforms.resize(GetMaxDelay()) );

    VT_HR_END();
}

HRESULT 
CVideoStabilizer::CVideoStabilizerInternal::PreProcessFrame(
    vt::CLumaByteImg& img, CStabilizerFrameData& data)
{
    VT_HR_BEGIN();

    // TODO: allow different image sizes; for now just use size of first frame
    m_frameWidth  = img.Width();
    m_frameHeight = img.Height();

    // compute scale for detection - based on larger of the two incoming
    // dimensions, start by going to 1/2 size then clamp to no larger
    // than the range 256..511
    //
    // (use max to make sure enough decimation occurs for wide video
    // formats, since things run a lot slower with larger detection 
    // dimensions)
    //
    int iFrameDim     = VtMax(m_frameWidth, m_frameHeight);
    int frameDimSteps = F2I(logf((float)iFrameDim)/logf(2.f));
    m_detectDimensionSteps = VtMax(frameDimSteps-8, 1);

    // then bias but don't allow it to get too small
    m_detectDimensionSteps += m_params.detectParams.detectDimensionBias;
    if ((frameDimSteps-m_detectDimensionSteps) <= 5) 
    {
        // result smaller than 32..63 range, so clamp detection decimation
        m_detectDimensionSteps = VtMax(0,frameDimSteps-5);
    }

    // run the feature detector
    vt::vector<HARRIS_FEATURE_POINT> curFP;
    vt::vector<BriefDesc<BRIEF_DESC_SIZE>> curDesc;
    float blurLevel;
    VT_HR_EXIT( RunDetectorAndExtractDescriptors(
        data.Impl()->curFP,data.Impl()->curDesc,
        img, m_params.detectParams, m_detectDimensionSteps, 
        m_bt, &blurLevel, m_hdetector) );
    data.Impl()->m_blurLevel = blurLevel;
    data.Impl()->m_imgWidth = img.Width();
    data.Impl()->m_imgHeight = img.Height();

    if (m_callbackmask & VSCB_DETECTEDCNT)
    {
        int cnt = (int)curFP.size();
        (m_callbackproc)(VSCB_DETECTEDCNT, (const void*)(&cnt), 0);
    }
    if (m_callbackmask & VSCB_BLURLEVEL)
    {
        (m_callbackproc)(VSCB_BLURLEVEL, (const void*)(&blurLevel), 0);
    }

    VT_HR_END();
}

HRESULT 
CVideoStabilizer::CVideoStabilizerInternal::PushFrame(CStabilizerFrameData& data)
{
    VT_HR_BEGIN();
    VT_ASSERT( m_bStarted == true ); // must call Start first

    if (m_pushedFrameCount == 0)
    {
        // TODO: allow different image sizes; for now just use size of first frame
        m_frameWidth  = data.Impl()->m_imgWidth;
        m_frameHeight = data.Impl()->m_imgHeight;

        // compute scale for detection - based on larger of the two incoming
        // dimensions, start by going to 1/2 size then clamp to no larger
        // than the range 256..511
        //
        // (use max to make sure enough decimation occurs for wide video
        // formats, since things run a lot slower with larger detection 
        // dimensions)
        //
        int iFrameDim     = VtMax(m_frameWidth, m_frameHeight);
        int frameDimSteps = F2I(logf((float)iFrameDim)/logf(2.f));
        m_detectDimensionSteps = VtMax(frameDimSteps-8, 1);

        // then bias but don't allow it to get too small
        m_detectDimensionSteps += m_params.detectParams.detectDimensionBias;
        if ((frameDimSteps-m_detectDimensionSteps) <= 5) 
        {
            // result smaller than 32..63 range, so clamp detection decimation
            m_detectDimensionSteps = VtMax(0,frameDimSteps-5);
        }

        // start-up the tracker
        VT_HR_EXIT( m_tracker.Begin(m_frameWidth, m_frameHeight, 
                                    m_detectDimensionSteps,
                                    m_detectDimensionSteps+5, 
                                    &m_params.trackParams) );

        for( int i = 0; i < (int)m_bufRollingShutterCorrect.size(); i++ )
        {
            VT_HR_EXIT( m_bufRollingShutterCorrect.buffer(i).resize(m_frameHeight) );
        }

        // start-up rolling shutter correct
        const float fSmoothness = 10.0f;
        VT_HR_EXIT( m_rollingShutterCorrect.Begin(
            m_frameWidth, m_frameHeight, fSmoothness) );

        // start-up the smoother
        if( m_params.stabilizeMode == 2 )
        {
            m_pSmoother = new CDPPathSmoother();
            VT_PTR_OOM_EXIT( m_pSmoother );

            // TODO: perhaps rename the m_params members to indicate generic 
            // smoother params, e.g. smootherMotionModel
            DPSMOOTH_PARAMS dpparams;
            dpparams.frameWidth  = m_frameWidth;
            dpparams.frameHeight = m_frameHeight;
            dpparams.fCrop       = m_params.fCropBox;
            dpparams.iBufferSize = m_params.stabilizeBuffer;
            dpparams.motionModel = m_params.motionModel;
            VT_HR_EXIT( ((CDPPathSmoother*)m_pSmoother)->SetParameters(dpparams) );

            VT_HR_EXIT( m_pSmoother->Begin() );
        }
        else if( m_params.stabilizeMode == 4 )
        {
            m_pSmoother = new CWLSSmoother();
            VT_PTR_OOM_EXIT( m_pSmoother );

            WLSSMOOTH_PARAMS wlsparams;
            wlsparams.frameWidth				= m_frameWidth;
            wlsparams.frameHeight				= m_frameHeight;
            wlsparams.kernelSize				= 15;
            wlsparams.outFramesPerOptimization	= 5;
            wlsparams.cacheSize					= 5;
            wlsparams.cropRatio                 = m_params.fCropBox;
            wlsparams.iBufferSize				= m_params.stabilizeBuffer;
            wlsparams.motionModel				= m_params.motionModel;
            VT_HR_EXIT( ((CWLSSmoother*)m_pSmoother)->SetParameters(wlsparams) );

            VT_HR_EXIT( m_pSmoother->Begin() );
        }
    }
    else
    {
        VT_ASSERT( data.Impl()->m_imgWidth == m_frameWidth );
        VT_ASSERT( data.Impl()->m_imgHeight == m_frameHeight );
    }

#if 0
    if (m_inBlurReset && (blurLevel > 20.f)) { m_inBlurReset = false; }
    else if (!m_inBlurReset && (blurLevel < 10.f)) { m_inBlurReset = true; }
    // if in reset, then clear out detected features
    // PERF: could avoid generating descriptors in this case
    if (m_inBlurReset) { curFP.clear(); curDesc.clear(); }
#endif

    if (m_callbackmask & VSCB_DETECTEDCNT)
    {
        int cnt = (int)data.Impl()->curFP.size();
        (m_callbackproc)(VSCB_DETECTEDCNT, (const void*)(&cnt), m_pushedFrameCount);
    }
    if (m_callbackmask & VSCB_BLURLEVEL)
    {
        (m_callbackproc)(VSCB_BLURLEVEL, (const void*)(&data.Impl()->m_blurLevel), m_pushedFrameCount);
    }

    // add frame to tracker
    VT_HR_EXIT( m_tracker.PushFrame(data.Impl()->curFP, data.Impl()->curDesc) );
    m_pushedFrameCount++;

    VT_HR_EXIT( AdvancePipeline() );

    VT_HR_END();
}

HRESULT 
CVideoStabilizer::CVideoStabilizerInternal::PushFrame(vt::CLumaByteImg& img)
{
    VT_HR_BEGIN();

    CStabilizerFrameData frameData;
    PreProcessFrame(img, frameData);
    VT_HR_EXIT( PushFrame(frameData) );

    VT_HR_END();
}

HRESULT
CVideoStabilizer::CVideoStabilizerInternal::GetResult(
    vt::CMtx3x3f& transform, vt::vector<CVec2f>* pvecRScorrect, int frameNumber) const
{
    VT_HR_BEGIN();

    if( m_params.stabilizeMode == 3 )
    {
        transform.MakeI();
    }
    else
    {
        VT_ASSERT( m_params.stabilizeMode < 5 );

        BUFFER_RANGE r = GetResultsRange();
        if( !r.InRange(frameNumber) )
        {
            VT_HR_EXIT( E_INVALIDARG );
        }

        // copy rolling shutter solution
        const vt::vector<CVec2f>& srcrs = m_bufRollingShutterCorrect[frameNumber];

        if (pvecRScorrect)
        {
            VT_HR_EXIT( pvecRScorrect->resize(srcrs.size()) );
            memcpy(pvecRScorrect->begin(), srcrs.begin(),
                   sizeof(CVec2f)*srcrs.size());
        }

        // get smoothed transform
        CMtx3x3f result = m_bufResultTransforms[frameNumber];
        
        // auto crop and (conditionally) clamp the result
        if( m_params.AutoCrop && m_params.stabilizeMode != 4 )
        {
            // don't clamp for fixed camera
            if (m_params.stabilizeMode != 0)
            {
                result = ClampSimilarityTransform(result, m_params.fCropBox, 
                                                  m_frameWidth, m_frameHeight, 
                                                  srcrs.begin(), (int)srcrs.size());
            }

            // apply the crop
            result = CMtx3x3f().MakeScale(m_params.fCropBox, m_params.fCropBox) * 
                result;
        }

        // modify the matrix to return value coordinate system which has 
        // top-left as origin
        if (m_params.stabilizeMode == 4)
        {
            transform = result;
        }
        else
        {
            transform = MatrixCenterToTopLeft(m_frameWidth, m_frameHeight) *
                result * MatrixTopLeftToCenter(m_frameWidth, m_frameHeight);
        }

    }

    VT_HR_END();
}


HRESULT
CVideoStabilizer::CVideoStabilizerInternal::End(void)
{
    // ignore multiple End() calls
    if (!m_bStarted) 
    { 
        return S_OK; 
    }

    VT_HR_BEGIN();

    // call end and flush the various stages of the pipeline in turn
    VT_HR_EXIT( m_tracker.End() );
    while (m_trackerNextFrame < m_pushedFrameCount)
    {
        VT_HR_EXIT( AdvancePipeline() );
    }

    VT_HR_EXIT( m_rollingShutterCorrect.End() );
    while (m_bufRollingShutterCorrect.get_total_count() < m_pushedFrameCount)
    {
        VT_HR_EXIT( AdvancePipeline() );
    }
 
    if( m_params.stabilizeMode == 2 || m_params.stabilizeMode == 4 )
    {		
        VT_ASSERT( m_pSmoother );
        ANALYZE_ASSUME( m_pSmoother != NULL );
        m_pSmoother->End();	 
        while (m_bufResultTransforms.get_total_count() < m_pushedFrameCount )
        {  
            VT_HR_EXIT( AdvancePipeline() );
        }
    }
    else
    {
        // ensure that for stab mode 0 there is no buffering
        // TODO: move mode 0 into a IFeatureWarpCompute
        VT_ASSERT( m_bufResultTransforms.get_total_count() == m_pushedFrameCount );
    }

    VT_HR_EXIT_LABEL()

    m_bStarted = false;
    if( m_pSmoother )
    {
        delete m_pSmoother;
        m_pSmoother = NULL;
    }

#ifdef INLIERPOINT_FILE
    if( g_fpFeaturesFile ) fclose(g_fpFeaturesFile);
#endif

    return hr;
}

void 
CVideoStabilizer::CVideoStabilizerInternal::SetCallback(
    unsigned int mask, void (__stdcall *proc)(unsigned int type, const void* data, int frame))
{
    m_callbackmask = mask;
    m_callbackproc = proc;
    m_tracker.SetCallback(mask, proc);
}

HRESULT
CVideoStabilizer::CVideoStabilizerInternal::AdvancePipeline()
{
    VT_HR_BEGIN()

    // read the tracker if a new result is available
    if ( m_tracker.GetResultsRange().InRange(m_trackerNextFrame) )
    {
        // get the tracked feature points from the tracker
        int iCurRefFrame;
        VT_HR_EXIT( m_tracker.GetResult(
            m_trackerNextFrame, m_vecMatches, iCurRefFrame) );

#if defined(INLIERPOINT_FILE)
        if (g_fpFeaturesFile == NULL)
        {
            _wfopen_s(&g_fpFeaturesFile, L"inlierpoints.txt", L"w");
        }
        if (g_fpFeaturesFile)
        {
            fprintf(g_fpFeaturesFile,"%03d %d\n", 
                    m_trackerNextFrame, m_vecMatches.size());
            for ( int j = 0; j < (int)m_vecMatches.size(); j++ )
            {
                const PointMatch& fpm = m_vecMatches[j];
                fprintf(g_fpFeaturesFile,"  %f %f    %f %f\n",
                        fpm.p1.x, fpm.p1.y, fpm.p0.x, fpm.p0.y);
            }
            fflush(g_fpFeaturesFile);
        }
#endif

        m_trackerNextFrame++;

        // TODO: add ability to skip rolling shutter
        // TODO: handle lost frames case for rolling shutter  
        // TODO: store the rolling buffer of matches here instead of in
        //       the rolling shutter component
        
        // add result to rolling shutter pipeline
        VT_HR_EXIT( m_rollingShutterCorrect.PushFrame(m_vecMatches) );
    }

    // read the rolling shutter component if a new result is available
    int iRSFrame = m_bufRollingShutterCorrect.get_total_count();
    if( m_rollingShutterCorrect.GetResultsRange().InRange(iRSFrame) )
    {
        // get the rolling shutter adjustment
        const float* pfCorrection;
        const vt::vector<PointMatch>* pMatches;
        VT_HR_EXIT( m_rollingShutterCorrect.GetResult(
            pfCorrection, pMatches, iRSFrame) );

        m_bufRollingShutterCorrect.advance();
		memcpy( m_bufRollingShutterCorrect.last().begin(), pfCorrection,
                m_frameHeight*sizeof(CVec2f) );

        // store the pfCorrection,  warp the features
        // TODO: Can most of the time reuse the warped match tracks from
        //       the previous frame.  For now it warps both this frame and
        //       previous.
        if( m_params.bRollingShutter && iRSFrame != 0 )
        {
            VT_HR_EXIT( WarpMatches(m_vecMatches, *pMatches, 
                                    m_bufRollingShutterCorrect[iRSFrame].begin(),
                                    m_bufRollingShutterCorrect[iRSFrame-1].begin(),
                                    m_frameHeight) );
            pMatches = &m_vecMatches;
        }

        if ( m_params.stabilizeMode == 2 || m_params.stabilizeMode == 4 )
        {
            VT_ASSERT( m_pSmoother );
            ANALYZE_ASSUME( m_pSmoother != NULL );
            VT_HR_EXIT( m_pSmoother->PushFeatures(
                pMatches->begin(), (int)pMatches->size()) );
        }
        else if ( m_params.stabilizeMode == 0 )
        {
            // TODO: move mode 0 into a 'smoother_' source file

            // compute the transformation
            m_bufResultTransforms.advance();
            if ( pMatches->size() == 0 )
            {
                // tracker is lost so reset the chained transforms
                m_mtxChained.MakeI();
                m_mtxPrev.MakeI();
                m_bufResultTransforms.last().MakeI();
            }
            else
            {
                CMtx3x3d mtxCur_d;
                if (m_params.motionModel == 2)
                {
                    VT_HR_EXIT( VtSimilarityFromPointMatches2D(
                        mtxCur_d, pMatches->begin(), (int)pMatches->size()) );
                }
                else
                {
                    VT_HR_EXIT( VtAffineFromPointMatches2D(
                        mtxCur_d, pMatches->begin(), (int)pMatches->size()) );
                }

                // all internal processing is done with coordinate system (0,0) at 
                // center of the image so adjust the returned matrix here
                CMtx3x3f mtxCur;
                VtConvertMtx(mtxCur_d, mtxCur);
                mtxCur = MatrixTopLeftToCenter(m_frameWidth, m_frameHeight) *
                         mtxCur * MatrixCenterToTopLeft(m_frameWidth, m_frameHeight);

                // for mode 0 store the transform matrix relative to frame 0
                m_mtxChained = m_mtxChained * m_mtxPrev;
                m_mtxPrev    = mtxCur;
                m_bufResultTransforms.last() = CMtx3x3f(m_mtxChained * mtxCur).Inv();
            }
        }
    }

    if( m_params.stabilizeMode == 2 || m_params.stabilizeMode == 4 )
    {
        VT_ASSERT( m_pSmoother );
        ANALYZE_ASSUME( m_pSmoother != NULL );
		if( m_pSmoother->GetResultsRange().InRange(m_bufResultTransforms.get_total_count()) )
        {
            m_bufResultTransforms.advance();
            VT_HR_EXIT( m_pSmoother->GetResult(
                m_bufResultTransforms.last(), m_bufResultTransforms.get_last_id()) );
        }
    }

    VT_HR_END()
}

//------------------------------------------------------------------------------
// end