//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      specialized implementation of Harris feature detector
//
//------------------------------------------------------------------------------
#include "stdafx.h"

using namespace vt;

// when defined, uses the transform framework for multi-proc per-block detection
#if WINAPI_FAMILY!=WINAPI_FAMILY_PHONE_APP // can't use transform framework on windows phone built via SDK
#define MULTICORE_HARRISDETECT
#endif

#include "vt_harrisdetect.h"

//------------------------------------------------------------------------------
// from harris.cpp

extern HRESULT HarrisDetectBlock(vt::vector<HARRIS_FEATURE_POINT>& pts, 
    const CImg& srcImg, const CRect& srcRect, const CPoint& ptSrc,
    int iWidth, int iHeight,
    const HARRIS_DETECTOR_PARAMS& prm, 
    bool useHFBlockCache, int ftrPerBlk);

//+-----------------------------------------------------------------------------
//
// multi-proc support for per-block detection via the transform framework;
// source image is passed in initialization routine (to avoid the copy
// currently done by the image reader), so intented to be used with
// CTransformGraphNoSrcNode
//
//------------------------------------------------------------------------------
#if defined(MULTICORE_HARRISDETECT)

class CHarrisDetectBlockTransform: public IImageTransform
{
    // IImageTransform implementation
public:
    virtual bool RequiresCloneForConcurrency()
    { return true; }

    virtual void    GetSrcPixFormat(IN OUT int* pfrmtSrcs, 
                                    IN UINT  /*uSrcCnt*/,
                                    IN int frmtDst)
    {
    }

    virtual void    GetDstPixFormat(OUT int& frmtDst,
                                    IN  const int* pfrmtSrcs, 
                                    IN  UINT  /*uSrcCnt*/)
    {
    }

    virtual HRESULT GetRequiredSrcRect(OUT TRANSFORM_SOURCE_DESC* pSrcReq,
                                      OUT UINT& uSrcReqCount,
                                      IN  UINT  /*uSrcCnt*/,
                                      IN  const CRect& rctLayerDst
                                      )
    { 
        uSrcReqCount = 0;
        return S_OK;
    }

    virtual HRESULT GetAffectedDstRect(OUT CRect& rctDst,
                                      IN  const CRect& rctSrc,
                                      IN  UINT /*uSrcIndex*/,
                                      IN  UINT /*uSrcCnt*/)
    { 
        rctDst = m_imgSrc.Rect();
        return S_OK;
    }

    virtual HRESULT GetResultingDstRect(OUT CRect& rctDst,
                                        IN  const CRect& rctSrc,
                                        IN  UINT /*uSrcIndex*/,
                                        IN  UINT /*uSrcCnt*/)
    { 
        rctDst = m_imgSrc.Rect();
        return S_OK;
    }

    virtual HRESULT Transform(OUT CImg* pimgDstRegion, 
                              IN  const CRect& rctLayerDst,
                              IN  CImg *const *ppimgSrcRegions,
                              IN  const TRANSFORM_SOURCE_DESC* pSrcDesc,
                              IN  UINT  uSrcCnt
                              )
    { 
        // compute index into block array of feature points
        int bx = rctLayerDst.left/m_iBlkW;
        int by = rctLayerDst.top/m_iBlkH;
        int idx = (by*m_iBlkCW)+bx; 
        VT_ASSERT( (*m_pfpa)[idx].size() == 0 );
        VT_ASSERT( idx < (int)m_pfpa->size() );

        HRESULT hr = HarrisDetectBlock((*m_pfpa)[idx], m_imgSrc, rctLayerDst, CPoint::CPoint(0,0),
            m_iWidth, m_iHeight,
            m_prms, true, m_iFtrPerBlk); 
        return hr;
    }

    virtual HRESULT Clone(ITaskState **ppState)
    {
        return CloneTaskState<CHarrisDetectBlockTransform>(ppState, 
            [this](CHarrisDetectBlockTransform* pN)
            { return pN->Initialize(m_pfpa, m_prms, m_iFtrPerBlk, m_imgSrc,m_iBlkW,m_iBlkH); });
    }

    HRESULT Merge(const ITaskState *pClone)
    {
        const CHarrisDetectBlockTransform* pCloneT =
            (CHarrisDetectBlockTransform*)pClone;
        return S_OK;
    }

    HRESULT Initialize(
        vt::vector<vt::vector<HARRIS_FEATURE_POINT>>* ppfa,
        const HARRIS_DETECTOR_PARAMS& prms,
        int iFtrPerBlk,
        CByteImg& imgSrc,
        int blkWidth, int blkHeight)
    {  
        VT_HR_BEGIN()
        m_pfpa = ppfa;
        VT_HR_EXIT( imgSrc.Share(m_imgSrc) );
        m_iWidth = m_imgSrc.Width();
        m_iHeight = m_imgSrc.Height();
        m_iBlkW = blkWidth;
        m_iBlkH = blkHeight;
        m_iBlkCW = (int)ceilf((float)m_iWidth/(float)m_iBlkW);

        m_prms = prms;
        m_iFtrPerBlk = iFtrPerBlk;

        VT_HR_END()
    };

protected:
    HARRIS_DETECTOR_PARAMS m_prms;
    int m_iFtrPerBlk;
    int m_iWidth, m_iHeight;
    CByteImg m_imgSrc;
    int m_iBlkW, m_iBlkH, m_iBlkCW;
public:
    vt::vector<vt::vector<HARRIS_FEATURE_POINT>>* m_pfpa;
};

#endif

//+-----------------------------------------------------------------------------
//
// main entrypoint for Harris detector
// 
//------------------------------------------------------------------------------
HRESULT CHarrisDetector::Detect(vt::vector<vt::vector<HARRIS_FEATURE_POINT>>& pts, 
    const vt::CByteImg& srcImgInput, const vt::CRect& srcRect,
    const HARRIS_DETECTOR_PARAMS& prm,
    HARRIS_DETECTOR_RESULTDATA& resultData)
{
    VT_HR_BEGIN()

    // iterate in ~equal sized blocks; ideally the perimeter blocks would be 
    // larger than the interior blocks to accommodate the border but this
    // regular (and simpler) blocking pushes for a few more features around
    // the perimeter which is beneficial in many cases anyway
    float fblksw = ceilf((float)(srcImgInput.Width()>>prm.decimate2to1Count)/(float)prm.blockSize);
    float fblksh = ceilf((float)(srcImgInput.Height()>>prm.decimate2to1Count)/(float)prm.blockSize);
    UInt32 blksw = (UInt32)ceilf((float)(srcImgInput.Width()>>prm.decimate2to1Count)/fblksw);
    UInt32 blksh = (UInt32)ceilf((float)(srcImgInput.Height()>>prm.decimate2to1Count)/fblksh);

    // compute feature budget per block
    int nblks = (int)(fblksw*fblksh);
    int ftrPerBlk = prm.featureCountTarget/nblks;
    // don't allow less than 4 features per block
    ftrPerBlk = VtMax(4,ftrPerBlk);

    // allocate feature point array entry for each block
    VT_HR_EXIT( pts.resize(nblks) );

    CByteImg srcImg;
    if (prm.decimate2to1Count == 0)
    {
        VT_HR_EXIT( srcImgInput.Share(srcImg) );
    }
    else
    {
        CLumaPyramid pyr;
        PYRAMID_PROPERTIES pyrProps;
        pyrProps.eAutoFilter = ePyramidFilter121Primal;
        VT_HR_EXIT( pyr.Create(srcImgInput, &pyrProps) );
        VT_HR_EXIT( pyr[prm.decimate2to1Count].Share(srcImg) );
    }

#if !defined(MULTICORE_HARRISDETECT)
    CRect srcRectDetectDimension = CRect(
        srcImgInput.Rect().left>>prm.decimate2to1Count,
        srcImgInput.Rect().top>>prm.decimate2to1Count,
        srcImgInput.Rect().right>>prm.decimate2to1Count,
        srcImgInput.Rect().bottom>>prm.decimate2to1Count);
    CBlockIterator bi(BLOCKITER_INIT(srcRectDetectDimension, blksw, blksh));
    for (int idx = 0; !bi.Done(); bi.Advance(), idx++ )
    {
        VT_HR_EXIT( HarrisDetectBlock(pts[idx], srcImg, bi.GetRect(), vt::CPoint(0,0),
            srcImg.Width(), srcImg.Height(),
            prm, true, ftrPerBlk) );
    }
#else
    // use transform framework to use multi-core for block processing
    VT_HR_EXIT( pts.resize(nblks) );
    CHarrisDetectBlockTransform x;
    VT_HR_EXIT( x.Initialize(&pts,prm,ftrPerBlk,srcImg,blksw,blksh) );
    CTransformGraphNoSrcNode g(&x);
	CRasterBlockMap bm(srcImg.Rect(), blksw, blksh);
	g.SetDest(NODE_DEST_PARAMS(NULL, &bm));
    VT_HR_EXIT( PushTransformTaskAndWait(&g, (CTaskStatus*)NULL) ); 
#endif

    {
        int count = 0;
        double sum = 0.;
        for (int i=0; i<(int)pts.size(); i++)
        {
            for (int j=0; j<(int)pts[i].size(); j++)
            {
                sum += pts[i][j].score;
                count++;
            }
        }
        float avgScore = (float)((1.*sum/(double)(count))/(double)0xffff);
        resultData.imageBlurLevel = avgScore;
    }

    VT_HR_END();
}

// end