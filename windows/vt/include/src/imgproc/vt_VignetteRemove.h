#pragma once

#include "vtcommon.h"

namespace vt {

enum eVignetteRemove
{
	eVRNone = 0,
    eVRLum  = 1,
	eVRRGB  = 2,
    eVRAuto = 3
};

struct VIGNETTE_PARAMS
{
    VIGNETTE_PARAMS() : fFactorR(0.f), fFactorG(0.f), fFactorB(0.f),
        fGlobalIntensityScale(1.f) {};

    float fFactorR, fFactorG, fFactorB;
    float fGlobalIntensityScale;
};

inline bool VignetteParamsActive(const VIGNETTE_PARAMS& vp)
{ 
	return vp.fFactorR > 0.f || vp.fFactorG > 0.f || vp.fFactorB > 0.f ||
		   vp.fGlobalIntensityScale < 1.f;
}

HRESULT ComputeVignetteParams(IImageReader* pReader, eVignetteRemove eVR,
                              OUT VIGNETTE_PARAMS& params,
                              CTaskProgress* pProgress = NULL);

class CVignetteAnalyze
{
public:
    CVignetteAnalyze()
    {
        Reset();
    }
    ~CVignetteAnalyze(void){}

	void Reset()
    {
        m_info = CImgInfo();
        m_uSourceId = 0;
        m_fEntropyRGB = m_fEntropyY = 0;
    }

    HRESULT AnalyzeImage(IImageReader* pReader,
                         IN UINT uSourceId);

    HRESULT GetBestImage(OUT UINT& uSourceId, bool bPerChannel)
    { 
        if (m_fEntropyRGB == 0 && (bPerChannel || m_fEntropyY == 0))
            return S_FALSE;
        else
            uSourceId = m_uSourceId;
        return S_OK;
    }

private:
    CImgInfo m_info;
    UINT m_uSourceId;
    float m_fEntropyRGB, m_fEntropyY;
};

//+-----------------------------------------------------------------------------
//
// Class: CVignetteRemoveTransform
// 
// Synposis: IImageTransform implementation of CVignetteRemove::CorrectImage()
// 
//------------------------------------------------------------------------------
class CVignetteRemoveTransform :
    public CImageTransformUnaryPoint<CVignetteRemoveTransform, false>
{
public:
	CVignetteRemoveTransform()
	{}
	 
	CVignetteRemoveTransform(int w, int h, VIGNETTE_PARAMS& params)
	{ Initialize(w, h, params); }

	~CVignetteRemoveTransform()
	{}

	void Initialize(int w, int h, VIGNETTE_PARAMS& params)
	{
		m_w = w;
		m_h = h;
		m_params = params;

	    // scale is the scale of src and dst compared with the original image
	    float scale = 0.5f * (float) w; // assumes aspect ratios are the same

        // note: modify channel because the order is actually BGRA
        m_f_inv[0] = params.fFactorB / scale;
        m_f_inv[1] = params.fFactorG / scale;
        m_f_inv[2] = params.fFactorR / scale;
        m_f_inv[3] = 0.f; // no-op for alpha

        VT_DEBUG_LOG( "IP_LOG_VIGNETTE: CVignetteRemoveTransform: %f %f %f (%f)\n",
            m_f_inv[2], m_f_inv[1], m_f_inv[0], params.fGlobalIntensityScale );
	}

public:
	virtual HRESULT Transform(OUT CImg* pimgDst, 
							  IN  const CRect& rctLayerDst,
							  IN  const CImg& imgSrc);

	virtual HRESULT Clone(ITaskState **ppState)
	{
		return CloneTaskState<CVignetteRemoveTransform>(ppState,
            VT_NOTHROWNEW CVignetteRemoveTransform(m_w, m_h, m_params));
    }

protected:
    int m_w, m_h;
    VIGNETTE_PARAMS m_params;
    float m_f_inv[4]; // vignetting parameters (Kang-Weiss'00)
};

//+-----------------------------------------------------------------------------
//
// function: PushVignetteTransformAndWait
// 
// Synposis: sets up a vignette correcting transform and runs it via the
//           multi-core transform framework 
// 
//------------------------------------------------------------------------------
inline
HRESULT PushVignetteTransformAndWait(IImageReaderWriter* pImage, 
                                     VIGNETTE_PARAMS& params,
									 CTaskProgress* pProgress = NULL, 
									 VT_TRANSFORM_TASK_OPTIONS* pOpts = NULL)
{
    HRESULT hr;

	vt::CRect rctDst = pImage->GetImgInfo().Rect();
    CVignetteRemoveTransform vrt(rctDst.Width(), rctDst.Height(), params);
    CTransformGraphUnaryNode graph(&vrt);
    IMAGE_EXTEND ex(Extend);
	if ((hr = graph.BindSourceToReader(pImage, &ex)) == S_OK)
    {
	    graph.SetDest(pImage);
	    hr = PushTransformTaskAndWait(&graph, pProgress, pOpts);
    }
    return hr;
}

};
