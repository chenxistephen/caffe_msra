#pragma once

#include "completion_tiles.h"

namespace vt {

// debug options
//#define DEBUG_OUTPUT
#define DEBUG_MATCH
//#define MCORE_MATCH

// optimizations
#define FAST_RAND

// multi-core completion
#define MULTI_CORE
#ifdef MULTI_CORE
#define MULTI_CORE_INIT
#define MULTI_CORE_ITERATION
#define MULTI_CORE_AVG
#endif MULTI_CORE

wstring wformat(const wchar_t * str, ...);

#ifdef FAST_RAND
struct FastRand
{
    int seed;

    FastRand(int s = 0)
    {
        Seed(s);
    }

    void Seed(int s)
    {
        seed = s;
    }

    int IRand(int mod)
    { 
       seed = (214013*seed+2531011);
       int val = (seed>>16)&0x7FFF;
       return val % mod;
    } 
};
#endif

class Completion
{
private:
    CVTResourceManagerInit m_resourceManager;

public:
    Completion();

	HRESULT init(IImageReaderWriter * pPyr, const wchar_t * pwszTimingFile,
        float discontinuousPenalty, float windowEnergyThresh);
    HRESULT finish();

#ifdef DEBUG_OUTPUT
    bool m_dbgSaveIterationData, m_dbgSaveAverageData;
    int m_dbgGen;
#endif DEBUG_OUTPUT

    HRESULT pass(const int numRepetitions, const bool lockImage);

    template <typename TUV>
    HRESULT upsample();

    const static int s_prad = 3, s_psize = 7;
    const static int s_bdDistWeightLength = 1024;
    const static Byte s_knownBit = 128, s_validBit = 64;
    const static int s_numIters = 2;    // down-up

    IImageReaderWriter * m_imgPyr;
    CImgInCache m_uvPyr;
    CByteImgInCache m_maskPyr;
    int m_level, m_tileLevel;
    int m_blockSize;

    template <typename TRGBA, typename TUV>
    HRESULT iterateBlock(int iter, int seed,
        const CRect& rctCur, const CRect& rctUv, CTypedImg<TUV>& uv,
        const CRect& rctDst, const CCompositeImg<RGBAType<TRGBA>>& imgDst,
        const CRect& rctSrc, const CCompositeImg<RGBAType<TRGBA>>& imgSrc
#ifdef DEBUG_OUTPUT
        , FILE * flog
#endif DEBUG_OUTPUT
        );

    template <typename TRGBA, typename TUV>
    HRESULT averageBlock(
        const CRect& rctDst, CCompositeImg<RGBAType<TRGBA>>& imgDst,
        const CRect& rctSrc, const CCompositeImg<RGBAType<TRGBA>>& imgSrc,
        const CRect& rctUv, const CTypedImg<TUV>& uv,
        const CRect& rctMrk, const CByteImg& marked);

    vector<CRect> m_bboxes[3];  // down iteration, up iteration, average
    vector<CRect> m_tiles;
#ifdef USE_TILE_BBOXES
    vector<CRect> m_unions;
#endif USE_TILE_BBOXES

private:
    template <typename TUV>
    HRESULT initLevel(vector<CVec2<TUV>>* validPixels = NULL, vector<CVec2<TUV>>* holePixels = NULL, vector<CVec2<TUV>>* uvPixels = NULL);
    template <typename TUV>
    HRESULT randomize(vector<CVec2<TUV>>& validPixels, vector<CVec2<TUV>>& holePixels, vector<CVec2<TUV>>& uvPixels);

#ifdef DEBUG_OUTPUT
    HRESULT saveUv(const wstring & fileName);
#endif DEBUG_OUTPUT

    HRESULT initBoundaryDistance();

    HRESULT average();

    HRESULT iteration(int iter);

    bool m_changed;

    Byte m_mark;

    vector<CFloatImg> m_restriction;

#ifdef FAST_RAND
    FastRand m_rnd;
#else FAST_RAND
    CRand m_rnd;
#endif FAST_RAND

    float m_weightGamma;
    float m_discontinuousPenalty;
    float m_windowEnergyThresh;

    CFloatImg m_bdDist;
    int m_bdDistLevel;
    vector<float> m_bdDistWeight;

    wstring m_csvFileName;

    float m_timeCompletion, m_timeInit, m_timeInitTiles,
        m_timePass, m_timePassIteration, m_timePassAverage,
        m_timeUpsample, m_timeUpsampleAverage, m_timeUpsampleInvalid,
        m_timeFusion;

    CTimer m_timer;
};

//+-----------------------------------------------------------------------------
//
// Class: CUnknownBlockMap
// 
// Synposis: Like CRasterBlockMap but returns empry rects for known regions.
// 
//------------------------------------------------------------------------------
class CUnknownBlockMap: public CRasterBlockMap
{
public:
	CUnknownBlockMap()
	{ m_iTotalBlkCnt = 0; }

	CUnknownBlockMap(const BLOCKITER_INIT& bii, 
					 int iBlkCntXInMaj=4, int iBlkCntYInMaj=4, int level=0)
	{ Initialize(bii, iBlkCntXInMaj, iBlkCntYInMaj, level); }

	CUnknownBlockMap(const vt::CRect& region, int iBlockWidth=256, int iBlockHeight=256,
					 int iBlkCntXInMaj=4, int iBlkCntYInMaj=4, int level=0)
	{ 
		Initialize(BLOCKITER_INIT(region,iBlockWidth,iBlockHeight),
	               iBlkCntXInMaj, iBlkCntYInMaj, level);
	}

    HRESULT SetKnown(vector<CRect>& bboxes)
    {
        HRESULT hr = S_OK;
        if ((int) bboxes.size() != m_iTotalBlkCnt)
            return E_INVALIDARG;
        if ((hr = m_vecKnown.resize(m_iTotalBlkCnt)) != S_OK)
            return hr;
        for (int i = 0; i < m_iTotalBlkCnt; i++)
            m_vecKnown[i] = bboxes[i].IsRectEmpty();
        return S_OK;
    }

	virtual HRESULT GetBlockItems(CRect& rct, int& iPrefetchCount, int iBlockIndex)
    {
        HRESULT hr = CRasterBlockMap::GetBlockItems(rct, iPrefetchCount, iBlockIndex);
        if (hr != S_OK)
            return hr;
        int block = rct.left / m_bii.iBlockWidth +
                    (rct.top / m_bii.iBlockHeight) * m_iBlkCntX;
        if (!m_vecKnown.empty() && m_vecKnown[block])
            rct = CRect(0,0,0,0);
        return S_OK;
    }
	  
	virtual CRect GetPrefetchRect(int iBlockIndex, int iPrefetchIndex)
    {
        CRect rct = CRasterBlockMap::GetPrefetchRect(iBlockIndex, iPrefetchIndex);
        int block = rct.left / m_bii.iBlockWidth +
                    (rct.top / m_bii.iBlockHeight) * m_iBlkCntX;
        if (!m_vecKnown.empty() && m_vecKnown[block])
            return CRect(0,0,0,0);
        else
            return rct;
    }

protected:
    vector<bool> m_vecKnown;
};
}
