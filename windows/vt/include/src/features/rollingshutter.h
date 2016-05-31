#pragma once

#ifndef VT_GCC
#include "vt_stabilize.h"
#endif

using namespace vt;

#ifndef VT_GCC
//+-----------------------------------------------------------------------------
//
// Class: RSC
//
//------------------------------------------------------------------------------
class RSC: public IDelayResult
{
public:
    int GetMaxDelay()
    { return m_iRSBufferSize/2+1; }

    BUFFER_RANGE GetResultsRange();

public:
    RSC();

    virtual ~RSC()
    {}

    HRESULT Begin(int iWidth, int iHeight, float fSmoothness);

    HRESULT PushFrame(const vt::vector<PointMatch>& vecFP);

    HRESULT GetResult(const float*& pfCorrection, 
                      const vt::vector<PointMatch>*& pMatches, 
                      int frameNumber);

    HRESULT End(void)
    { 
        m_bStarted = false;
        return S_OK;
    }

protected:
    void    Deallocate();

    HRESULT SolveForMotion();

    void ConjugateGradient(CVec<float> &vecX, CVec<float> &vecB);

    HRESULT ComputeCorrectionFromMotion(int iBufIdx);

private:
    bool m_bStarted;
    int m_iInputWidth;
    int m_iInputHeight;
    int m_iRSBufferSize;
    int m_iFrameCnt;

    // source features
    CRollingBuffer<vt::vector<PointMatch>> m_correspondences;

    // Solver Datastructures
    float m_fTemporalSmoothness;
    int   m_iMotionSamples;
    int   m_iMaxCGInterations;
    float m_fCGThreshold1;
    float m_fCGThreshold2;
    float m_fAlphaMultiplier;
    CMtx<float> m_matA;
    CVec<float> m_vecBX;
    CVec<float> m_vecBY;
    CVec<float> m_vecX;
    CVec<float> m_vecY;
    CVec<float> m_vecR;
    CVec<float> m_vecP;
    CVec<float> m_vecAP;

    vt::vector<int> m_vecStartNonZero;
    vt::vector<int> m_vecEndNonZero;

    // Correction Datastructures
    vt::vector<float> m_vecMotionX;
    vt::vector<float> m_vecMotionY;
    vt::vector<float> m_vecTranslationX;
    vt::vector<float> m_vecTranslationY;
    vt::vector<float> m_vecCorrection;
};
#endif // VT_GCC

//+-----------------------------------------------------------------------------
//
// Class: CRollingShutterAddressGen
//
//------------------------------------------------------------------------------
class CRollingShutterAddressGen: public IAddressGenerator
{
public: 
    CRollingShutterAddressGen() :
        m_iDstW(0), m_iDstH(0), m_bApplyMatrix(false), m_bAffineMatrix(true)
    {}
    CRollingShutterAddressGen(const CVec2f* pAdj, int iDstW, int iDstH, 
                              float fScale=1.f, const CMtx3x3f* pMtx = NULL) :
        m_iDstW(0), m_iDstH(0), m_bApplyMatrix(false), m_bAffineMatrix(true)
    { Initialize(pAdj, iDstW, iDstH, fScale, pMtx); }

public:
    virtual vt::CRect MapDstRectToSrc(IN const CRect& rctDst);
    virtual vt::CRect MapSrcRectToDst(IN const CRect& rctSrc)
    { return vt::CRect(0,0,m_iDstW,m_iDstH); }
    virtual HRESULT   MapDstSpanToSrc(OUT CVec2f* pSpan, 
                                      IN const vt::CPoint &ptDst, int iSpan);
    virtual HRESULT   MapDstAddrToSrc(IN OUT CVec2f* pSpan, int iSpan);
    virtual HRESULT   Clone(IAddressGenerator **ppClone)
    {
        HRESULT hr = CloneAddressGenerator<CRollingShutterAddressGen>(ppClone);
        if( hr == S_OK )
        {
            hr = ((CRollingShutterAddressGen*)*ppClone)->Initialize(
                m_pAdj, m_iDstW, m_iDstH, m_fScale, m_bApplyMatrix?(&m_xfrm):(NULL));
        }
        return hr;
    }

public:
    HRESULT Initialize(const CVec2f* pAdj, int iDstW, int iDstH, float fScale=1.f, const CMtx3x3f* pMtx = NULL)
    { 
        m_iDstW  = iDstW;
        m_iDstH  = iDstH;
        m_fScale = fScale;
        m_fScaleInv = 1.f/fScale;
        m_pAdj = pAdj;

        if (pMtx)
        {
            m_bApplyMatrix = true;
            m_xfrm = *pMtx;
            m_bAffineMatrix = IsMatrixAffine(m_xfrm, vt::CRect(0,0,iDstW,iDstH)); 
        }
        return S_OK;
    }

protected:
    const CVec2f*   m_pAdj;
    int             m_iDstW;
    int             m_iDstH;
    float           m_fScale, m_fScaleInv;
    bool            m_bApplyMatrix;
    bool            m_bAffineMatrix;
    CMtx3x3f        m_xfrm;
};
