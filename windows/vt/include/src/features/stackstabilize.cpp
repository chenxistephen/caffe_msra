//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      image stack stabilization API
//
//------------------------------------------------------------------------------
#include "stdafx.h"

using namespace vt;

#include "vt_stabilize.h"
#include "stabilizeint.h"
#include "brief.h"
#include "tracker.h"
#include "vt_harrisdetect.h"
//#include "FAST10.h"

//------------------------------------------------------------------------------
// CStabilizerFrameData implementation
//------------------------------------------------------------------------------
CStabilizerFrameData::CStabilizerFrameData(): m_pint(new CStabilizerFrameDataInternal) 
{
}
CStabilizerFrameData::~CStabilizerFrameData() 
{ delete m_pint; }

float CStabilizerFrameData::GetBlurLevel()
{
    return (m_pint==NULL)?(0.f):(m_pint->m_blurLevel);
}
int CStabilizerFrameData::GetDetectedFeatureCount()
{
    return (m_pint==NULL)?(0):((int)(m_pint->curFP.size()));
}
CStabilizerFrameData::CStabilizerFrameDataInternal* CStabilizerFrameData::Impl(void)
{
    return m_pint;
}

//------------------------------------------------------------------------------
// internal class for stack stabilizer
//------------------------------------------------------------------------------
class CStackStabilizer::CStackStabilizerInternal
{
public:
    CStackStabilizerInternal()
        :m_frameWidth(0),m_frameHeight(0)
        ,m_detectDimensionSteps(0)
        ,m_framesSinceRef(0)
        ,m_callbackmask(0x0),m_callbackproc(NULL)
    {
        // create the Brief table 
        m_bt.Initialize(BRIEF_PATCH_SIZE);
    };
    ~CStackStabilizerInternal()
    { };

    // internal methods forwarded from CStabilize methods of the same name
    HRESULT GetParams(STACK_STABILIZER_PARAMS& params);
    HRESULT SetParams(const STACK_STABILIZER_PARAMS& params);
    HRESULT SetReferenceFrame(vt::CLumaByteImg& img);
    HRESULT AlignFrame(vt::CLumaByteImg& img, vt::CMtx3x3f& transform);
    void SetCallback(unsigned int mask, void (__stdcall *proc)(unsigned int type, const void* data, int frame));

protected:
    int m_frameWidth, m_frameHeight;
    int m_detectDimensionSteps;

    STACK_STABILIZER_PARAMS m_params;

    BriefTable<BRIEF_DESC_SIZE> m_bt;

    CBriefFeatureTracker m_tracker;

    vt::vector<HARRIS_FEATURE_POINT> m_refFrameFP;
    vt::vector<BriefDesc<BRIEF_DESC_SIZE>> m_refFrameDesc;

    int m_framesSinceRef;

    CHarrisDetector m_hdetector;

	vt::vector<PointMatch> m_vecMatches;

    // debug callback
    unsigned int m_callbackmask;
    void (__stdcall *m_callbackproc)(unsigned int type, const void* data, int frame);
};

CStackStabilizer::CStackStabilizer(): m_pint(new CStackStabilizerInternal) 
{}

CStackStabilizer::~CStackStabilizer() 
{ delete m_pint; }

HRESULT CStackStabilizer::GetParams(STACK_STABILIZER_PARAMS& params)
{
    return (m_pint==NULL)? 
        E_OUTOFMEMORY: m_pint->GetParams(params);
}

HRESULT CStackStabilizer::SetParams(const STACK_STABILIZER_PARAMS& params)
{
    return (m_pint==NULL)? 
        E_OUTOFMEMORY: m_pint->SetParams(params);
}

HRESULT CStackStabilizer::SetReferenceFrame(vt::CLumaByteImg& img)
{
    return (m_pint==NULL)? 
        E_OUTOFMEMORY: m_pint->SetReferenceFrame(img);
}

HRESULT CStackStabilizer::AlignFrame(vt::CLumaByteImg& img, vt::CMtx3x3f& transform)
{
    return (m_pint==NULL)? 
        E_OUTOFMEMORY: m_pint->AlignFrame(img, transform);
}

void 
CStackStabilizer::SetCallback(unsigned int mask, 
    void (__stdcall *proc)(unsigned int type, const void* data, int frame))
{
    if( m_pint )
    {
        m_pint->SetCallback(mask, proc);
    }
}

//+-----------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
#if WINAPI_FAMILY!=WINAPI_FAMILY_PHONE_APP // can't use transform framework on windows phone built via SDK
#define MULTICORE_DOWNSAMPLE
#endif
static HRESULT
CreateFeatureExtractImage(CByteImg& imgDst, const CByteImg& imgSrc, 
                          int iEndLevel, float fBlurSigma)
{
    VT_HR_BEGIN()
#if !defined(MULTICORE_DOWNSAMPLE)
    // build a temporary pyramid - uses single-core Decimate2to1
    CLumaPyramid pyr;
    PYRAMID_PROPERTIES pyrProps;
    pyrProps.eAutoFilter = ePyramidFilter121Primal;
    VT_HR_EXIT( pyr.Create(imgSrc, &pyrProps) );

    // pre-blur the down-sampled image if requested
    if (fBlurSigma > 0.f)
    {
        VT_HR_EXIT( VtGaussianSmooth(imgDst, pyr[iEndLevel], fBlurSigma) );
    }
    else
    {
        pyr[iEndLevel].Share(imgDst);
    }
#else
    // use transform version of 121 decimate
    int iW = imgSrc.Width();
    int iH = imgSrc.Height();
    CImg imgDstI;
    for (int i=0; i<iEndLevel; i++, iW >>= 1, iH >>= 1)
    {
        // source is function input or result of previous pass; transform
        // function takes cimg, so need to wrap anyway
        CImg cimgSrc; (i==0)?(imgSrc.Share(cimgSrc)):(imgDstI.Share(cimgSrc));

        // create image for dest
        VT_HR_EXIT( imgDstI.Create(iW>>1, iH>>1, OBJ_BYTEIMG) );

        // apply 121 via transform framework using no-copy method (no reader/writers)
        CSeparableFilter121Transform x;
        VT_HR_EXIT( x.InitializeDecimate2to1(imgDstI, cimgSrc) );
        CTransformGraphNoSrcNode g(&x);
        // TODO: tune this for larger image sizes; do two blocks (top/bottom halves) for 720p 
		CRasterBlockMap bm(imgDstI.Rect(), VtMin(640,iW/2), VtMin(512,iH/4));
		g.SetDest(NODE_DEST_PARAMS(NULL, &bm));
		VT_HR_EXIT( PushTransformTaskAndWait(&g, (CTaskStatus*)NULL) ); 

        // deallocate temporary source image for middle passes
        if ((i != 0) && (i != (iEndLevel-1))) { cimgSrc.Deallocate(); }
    }
    VT_HR_EXIT( imgDstI.Share(imgDst) );
#endif
    VT_HR_END()
}

HRESULT 
RunDetectorAndExtractDescriptors(vt::vector<HARRIS_FEATURE_POINT>& fpnts, 
                                 vt::vector<BriefDesc<BRIEF_DESC_SIZE>>& desc,
                                 const vt::CLumaByteImg& img, 
                                 const FEATURE_DETECTOR_PARAMS& detectParams, int iExtractLevel,
                                 const BriefTable<BRIEF_DESC_SIZE>& briefTable,
                                 float* pBlurLevel,
                                 CHarrisDetector& hdetector)
{
    VT_HR_BEGIN();

    const float fExtractLevelScale = float(1<<iExtractLevel);

    CByteImg imgFeature;

     // downsample the the requested feature extract level
    if (iExtractLevel > 0)
    {
        VT_HR_EXIT( CreateFeatureExtractImage(imgFeature, img, 
            iExtractLevel, detectParams.preBlur) );
    }
    else
    {
        // no downsample, so just share
        img.Share(imgFeature);
    }

    // extract features
    vt::vector<vt::vector<HARRIS_FEATURE_POINT>> fpa;
//    if (detectParams.detector == 0)
    {
        // note that the border must be set such that no features are returned
        // that would result in the patch not being wholly within the image
        HARRIS_DETECTOR_PARAMS hparams = detectParams.harrisParams;
        hparams.border = (BRIEF_PATCH_SIZE/2)+1;
        HARRIS_DETECTOR_RESULTDATA hresdata;
        VT_HR_EXIT( hdetector.Detect(fpa, imgFeature, imgFeature.Rect(), hparams, hresdata) );
        if (pBlurLevel) { *pBlurLevel = hresdata.imageBlurLevel; }
    }
//    else
//    {
//        VT_HR_EXIT( FAST10Detect(fpnts, imgFeature, imgFeature.Rect(), 
//                                 (BRIEF_PATCH_SIZE/2)+1, 
//                                 detectParams.FAST10Threshold) );
//    }

    // compute the total number of features
    int iNumFeatures = 0;
    for (int i=0; i<(int)fpa.size(); i++)
    {
        iNumFeatures += (int)fpa[i].size();
    }

    // allocate and compute the descriptors for the current frame; allocate and compute
    // feature point array for the current frame - returned feature points are scaled to
    // the original (pre-decimation) image coordinates
    VT_HR_EXIT( fpnts.resize(iNumFeatures) );
    VT_HR_EXIT( desc.resize(iNumFeatures) );
    int fpi = 0; // index for copying content of 2D feature point array to 1D arrays
    for ( int i = 0; i < (int)fpa.size(); i++ )
    {
        for( int j = 0; j < (int)fpa[i].size(); j++ )
        {
            HARRIS_FEATURE_POINT& fp = fpa[i][j];

#if ( (BRIEF_DESC_SIZE == 128) && (BRIEF_PATCH_SIZE == 24) )
            if (imgFeature.StrideBytes() == 640)
            {
                ComputeBriefDescriptor_d128_p24_s640(desc[fpi],fp,imgFeature);
            }
            else if (imgFeature.StrideBytes() == 320)
            {
                ComputeBriefDescriptor_d128_p24_s320(desc[fpi],fp,imgFeature);
            }
            else if (imgFeature.StrideBytes() == 180)
            {
                ComputeBriefDescriptor_d128_p24_s180(desc[fpi],fp,imgFeature);
            }
            else
#endif
            {
                BriefEvaluatePatch<BRIEF_DESC_SIZE, BRIEF_PATCH_SIZE>(desc[fpi], 
                    briefTable, fp, imgFeature);
            }

            fpnts[fpi] = fp;
            fpnts[fpi].x *= fExtractLevelScale;
            fpnts[fpi].y *= fExtractLevelScale;

            fpi++;
        }
    }

    VT_HR_END();
}

//+-----------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
HRESULT
CStackStabilizer::CStackStabilizerInternal::GetParams(STACK_STABILIZER_PARAMS& params)
{
    params = m_params;
    return S_OK;
}

HRESULT
CStackStabilizer::CStackStabilizerInternal::SetParams(const STACK_STABILIZER_PARAMS& params)
{
    m_params = params;
    return S_OK;
}

HRESULT 
CStackStabilizer::CStackStabilizerInternal::SetReferenceFrame(vt::CLumaByteImg& img)
{
    VT_HR_BEGIN();

    // TODO: allow different image sizes; for now just use size of reference frame
    m_frameWidth = img.Width();
    m_frameHeight = img.Height();

    // compute scale for detection - based on larger of the two incoming
    // dimensions, start by going to 1/2 size then clamp to no larger
    // than the range 256..511
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

    VT_HR_EXIT( RunDetectorAndExtractDescriptors(m_refFrameFP, m_refFrameDesc, 
        img, m_params.detectParams, m_detectDimensionSteps, m_bt, NULL, m_hdetector) );

    if (m_callbackmask & VSCB_DETECTEDCNT)
    {
        int cnt = (int)m_refFrameFP.size();
        (m_callbackproc)(VSCB_DETECTEDCNT, (const void*)(&cnt), 0);
    }

    m_framesSinceRef = 0;

    VT_HR_END();
}


HRESULT
CStackStabilizer::CStackStabilizerInternal::AlignFrame(vt::CLumaByteImg& img, 
    vt::CMtx3x3f& transform)
{
    VT_HR_BEGIN();
    
    int refId;
    m_framesSinceRef++;

    // start-up the tracker
    VT_HR_EXIT( m_tracker.Begin(m_frameWidth, m_frameHeight, 
        m_detectDimensionSteps, m_detectDimensionSteps+5, &m_params.trackParams) );

    // add reference frame
    VT_HR_EXIT( m_tracker.PushFrame(m_refFrameFP, m_refFrameDesc) );

    // run detector on alignment frame and add to tracker
    vt::vector<HARRIS_FEATURE_POINT> frameFP;
    vt::vector<BriefDesc<BRIEF_DESC_SIZE>> frameDesc;
    VT_HR_EXIT( RunDetectorAndExtractDescriptors(frameFP, frameDesc, img,
        m_params.detectParams, m_detectDimensionSteps, m_bt, NULL, m_hdetector) );
    if (m_callbackmask & VSCB_DETECTEDCNT)
    {
        int cnt = (int)frameFP.size();
        (m_callbackproc)(VSCB_DETECTEDCNT, (const void*)(&cnt), m_framesSinceRef);
    }
    VT_HR_EXIT( m_tracker.PushFrame(frameFP, frameDesc) );
    VT_HR_EXIT( m_tracker.End() );
    VT_HR_EXIT( m_tracker.GetResult(1, m_vecMatches, refId) );

    // for stack stabilizer, it makes sense to return tracking failure via HRESULT
    if (refId != 0) 
	{
		VT_HR_EXIT( E_FAIL ); 
	}

	// compute the transform from the point matches
	CMtx3x3d transform_d;
	if (m_params.motionModel == 2)
	{
		VT_HR_EXIT( VtSimilarityFromPointMatches2D(
			transform_d, m_vecMatches.begin(), (int)m_vecMatches.size()) );
	}
	else
	{
		VT_HR_EXIT( VtAffineFromPointMatches2D(
			transform_d, m_vecMatches.begin(), (int)m_vecMatches.size()) );
	}
  
	// all internal processing is done with coordinate system (0,0) at 
	// center of the image so adjust the returned matrix here
	VtConvertMtx(transform_d, transform);

    // invert for return so transform can be passed directly to VtWarpImage
    transform = transform.Inv();

    if ( m_params.AutoCrop )
    {
        // for now just a simple centered 94% crop
        float fScale  = 1.f/0.90f;
        float fTransX = float(m_frameWidth)*0.05f;
        float fTransY = float(m_frameHeight)*0.05f;
        CMtx3x3f cropMtx = CMtx3x3f().MakeScale(fScale, fScale)*
                           CMtx3x3f().MakeTranslate(-fTransX,-fTransY);
        transform = transform * cropMtx.Inv();
    }

    VT_HR_END();
}

void 
CStackStabilizer::CStackStabilizerInternal::SetCallback(
    unsigned int mask, void (__stdcall *proc)(unsigned int type, const void* data, int frame))
{
    m_callbackmask = mask;
    m_callbackproc = proc;
    m_tracker.SetCallback(mask, proc);
}

//------------------------------------------------------------------------------
// end