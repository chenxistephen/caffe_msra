//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Steerable filter implementation
//
//  History:
//      2008/4/4-swinder
//			Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_image.h"
#include "vt_kernel.h"
#include "vt_separablefilter.h"
#include "vt_steerablefilter.h"

using namespace vt;

HRESULT CSteerableFilter::Create(int iW, int iH, FilterType eft, FilterOrder efo, 
                                 float fScale, int iOrientationCount)
{
    m_fo = efo;
    m_ft = eft;

    m_iWidth = iW;
    m_iHeight = iH;

    if(fScale<1)
        return E_INVALIDARG;

    HRESULT hr;
    VT_HR_RET( m_v.Create(11) );
    if(iOrientationCount>0)
    {
        VT_HR_RET( m_vvCoefs.resize(iOrientationCount) );
        int i;
        for(i=0; i<iOrientationCount; i++)
        {
            VT_HR_RET( m_vvCoefs[i].Create(11) );
            m_vvCoefs[i] = *GetCoeffsByAngle(i * (float)VT_PI/iOrientationCount);
        }
    }

    if(m_fo==SecondOrder)
    {
        // second order

#if 0
        static float rgf1[] = { 0.0094f, 0.1148f, 0.3964f, -0.0601f, -0.9213f, -0.0601f, 0.3964f, 0.1148f, 0.0094f };
        static float rgf2[] = { 0.0008f, 0.0176f, 0.1660f, 0.6383f, 1.0f, 0.6383f, 0.1660f, 0.0176f, 0.0008f };
        static float rgf3[] = { 0.0028f, 0.0480f, 0.3020f, 0.5806f, 0, -0.5806f, -0.3020f, -0.0480f, -0.0028f };
        static float rgf4[] = { 0.0098f, 0.0618f, -0.0998f, -0.7551f, 0, 0.7551f, 0.0998f, -0.0618f, -0.0098f };
        static float rgf5[] = { 0.0020f, 0.0354f, 0.2225f, 0.4277f, 0, -0.4277f, -0.2225f, -0.0354f, -0.0020f };
        static float rgf6[] = { 0.0048f, 0.0566f, 0.1695f, -0.1889f, -0.7349f, -0.1889f, 0.1695f, 0.0566f, 0.0048f };

        VT_HR_RET( m_kern[0].Create(9, 4, rgf1) );
        VT_HR_RET( m_kern[1].Create(9, 4, rgf2) );
        VT_HR_RET( m_kern[2].Create(9, 4, rgf3) );
        VT_HR_RET( m_kern[3].Create(9, 4, rgf4) );
        VT_HR_RET( m_kern[4].Create(9, 4, rgf5) );
        VT_HR_RET( m_kern[5].Create(9, 4, rgf6) );
#else
        // the following functions recreate the filter kernels above for fScale=1

        int iTaps2 = (int)ceil(fScale * 4);
        int iTaps = iTaps2 * 2 + 1;

        int i;

        for(i=0; i<6; i++)
            VT_HR_RET( m_kern[i].Create(iTaps, iTaps2) );

        for(i = -iTaps2; i<=iTaps2; i++)
        {
            int j = i + iTaps2;
            float x = i * 0.67f / fScale;
            float xs = x * x;
            float fexp = exp(-xs) / fScale; 
            // since its kind of arbitrary how to scale these, I am scaling 
            // down with kernel area

            m_kern[0][j] = 0.9213f * (2*xs -1) * fexp;
            m_kern[1][j] = fexp;
            m_kern[2][j] = -sqrt(1.843f) * x * fexp;
            m_kern[3][j] = -0.978f * x * (-2.254f + xs) * fexp;
            m_kern[4][j] = -x * fexp;
            m_kern[5][j] = 0.978f * (-0.7515f + xs) * fexp;
        }

#endif

        if(m_ft==FilterEven || m_ft==FilterBoth)
        {
            VT_HR_RET( m_img[0].Create(iW, iH) );
            VT_HR_RET( m_img[1].Create(iW, iH) );
            VT_HR_RET( m_img[2].Create(iW, iH) );
        }

        if(m_ft==FilterOdd || m_ft==FilterBoth)
        {
            VT_HR_RET( m_img[5].Create(iW, iH) );
            VT_HR_RET( m_img[6].Create(iW, iH) );
            VT_HR_RET( m_img[7].Create(iW, iH) );
            VT_HR_RET( m_img[8].Create(iW, iH) );
        }

        m_iLastEven = 2;
        m_iLastOdd = 8;
    }
    else
    {
        // fourth order

#if 0
        static float rgf1[] = { 0.0084f, 0.0507f, 0.1084f, -0.1231f, -0.5729f, 0.0606f, 0.9344f,
            0.0606f, -0.5729f, -0.1231f, 0.1084f, 0.0507f, 0.0084f };
        static float rgf2[] = { 0.0001f, 0.0019f, 0.0183f, 0.1054f, 0.3679f, 0.7788f, 1.0f,
            0.7788f, 0.3679f, 0.1054f, 0.0183f, 0.0019f, 0.0001f };
        static float rgf3[] = { 0.0028f, 0.0229f, 0.0916f, 0.1186f, -0.1839f, -0.4867f, 0.0f,
            0.4867f, 0.1839f, -0.1186f, -0.0916f, -0.0229f, -0.0028f };
        static float rgf4[] = { 0.0005f, 0.0060f, 0.0456f, 0.1970f, 0.4583f, 0.4851f, 0.0f,
            -0.4851f, -0.4583f, -0.1970f, -0.0456f, -0.0060f, -0.0005f };
        static float rgf5[] = { 0.0012f, 0.0124f, 0.0715f, 0.2059f, 0.2053f, -0.2173f, -0.5581f,
            -0.2173f, 0.2053f, 0.2059f, 0.0715f, 0.0124f, 0.0012f };

        static float rgf6[] = { 0.0030f, -0.0012f, -0.0993f, -0.2908f, 0.1006f, 0.8322f, 0.0f,
            -0.8322f, -0.1006f, 0.2908f, 0.0993f, 0.0012f, -0.0030f };
        static float rgf7[] = { 0.0021f, 0.0095f, -0.0041f, -0.1520f, -0.3017f, 0.1161f, 0.5715f,
            0.1161f, -0.3017f, -0.1520f, -0.0041f, 0.0095f, 0.0021f };
        static float rgf8[] = { 0.0004f, 0.0048f, 0.0366f, 0.1581f, 0.3679f, 0.3894f, 0.0f,
            -0.3894f, -0.3679f, -0.1581f, -0.0366f, -0.0048f, -0.0004f };
        static float rgf9[] = { 0.0010f, 0.0077f, 0.0258f, 0.0016f, -0.1791f, -0.3057f, 0.0f,
            0.3057f, 0.1791f, -0.0016f, -0.0258f, -0.0077f, -0.0010f };
        static float rgf10[] = { 0.0010f, 0.0108f, 0.0611f, 0.1672f, 0.1237f, -0.3223f, -0.6638f,
            -0.3223f, 0.1237f, 0.1672f, 0.0611f, 0.0108f, 0.0010f };

        VT_HR_RET( m_kern[0].Create(13, 6, rgf1) );
        VT_HR_RET( m_kern[1].Create(13, 6, rgf2) );
        VT_HR_RET( m_kern[2].Create(13, 6, rgf3) );
        VT_HR_RET( m_kern[3].Create(13, 6, rgf4) );
        VT_HR_RET( m_kern[4].Create(13, 6, rgf5) );
        VT_HR_RET( m_kern[5].Create(13, 6, rgf6) );
        VT_HR_RET( m_kern[6].Create(13, 6, rgf7) );
        VT_HR_RET( m_kern[7].Create(13, 6, rgf8) );
        VT_HR_RET( m_kern[8].Create(13, 6, rgf9) );
        VT_HR_RET( m_kern[9].Create(13, 6, rgf10) );

#else
        // the following functions recreate the filter kernels above for fScale=1
        int iTaps2 = (int)ceil(fScale * 6);
        int iTaps = iTaps2 * 2 + 1;

        int i;
        for(i=0; i<10; i++)
            VT_HR_RET( m_kern[i].Create(iTaps, iTaps2) );

        for(i = -iTaps2; i<=iTaps2; i++)
        {
            int j = i + iTaps2;
            float x = i * 0.5f / fScale;
            float xs = x * x;
            float fexp = exp(-xs) / fScale;
            // since its kind of arbitrary how to scale these, I am scaling down 
            // with kernel area

            m_kern[0][j] = 1.246f * (0.75f + xs * (xs - 3)) * fexp;
            m_kern[1][j] = fexp;
            m_kern[2][j] = (1.5f - xs) * x * fexp;
            m_kern[3][j] = -1.246f * x * fexp;
            m_kern[4][j] = sqrt(1.246f) * (xs - 0.5f) * fexp;
            m_kern[5][j] = -0.3975f * (7.189f + xs * (xs - 7.501f)) * x * fexp;
            m_kern[6][j] = 0.3975f * (1.438f + xs * (xs - 4.501f)) * fexp;
            m_kern[7][j] = -x * fexp;
            m_kern[8][j] = -0.3975f * (xs - 2.225f) * x * fexp;
            m_kern[9][j] = (xs - 0.6638f) * fexp;
        }
#endif

        if(m_ft==FilterEven || m_ft==FilterBoth)
        {
            VT_HR_RET( m_img[0].Create(iW, iH) );
            VT_HR_RET( m_img[1].Create(iW, iH) );
            VT_HR_RET( m_img[2].Create(iW, iH) );
            VT_HR_RET( m_img[3].Create(iW, iH) );
            VT_HR_RET( m_img[4].Create(iW, iH) );
        }

        if(m_ft==FilterOdd || m_ft==FilterBoth)
        {
            VT_HR_RET( m_img[5].Create(iW, iH) );
            VT_HR_RET( m_img[6].Create(iW, iH) );
            VT_HR_RET( m_img[7].Create(iW, iH) );
            VT_HR_RET( m_img[8].Create(iW, iH) );
            VT_HR_RET( m_img[9].Create(iW, iH) );
            VT_HR_RET( m_img[10].Create(iW, iH) );
        }

        m_iLastEven = 4;
        m_iLastOdd = 10;
    }

    return NOERROR;
}

void CSteerableFilter::Deallocate()
{
    int i;
    for(i=0; i<11; i++)
        m_img[i].Deallocate();
    m_vvCoefs.clear();
}

HRESULT CSteerableFilter::Update(const CFloatImg &img)
{
    HRESULT hr;
    if(m_fo==SecondOrder)
    {
        if(m_ft==FilterEven || m_ft==FilterBoth)
        {
            VT_HR_RET( VtSeparableFilter(m_img[0], img, m_kern[0], m_kern[1]) );
            VT_HR_RET( VtSeparableFilter(m_img[1], img, m_kern[2], m_kern[2]) );
            VT_HR_RET( VtSeparableFilter(m_img[2], img, m_kern[1], m_kern[0]) );
        }

        if(m_ft==FilterOdd || m_ft==FilterBoth)
        {
            VT_HR_RET( VtSeparableFilter(m_img[5], img, m_kern[3], m_kern[1]) );
            VT_HR_RET( VtSeparableFilter(m_img[6], img, m_kern[5], m_kern[4]) );
            VT_HR_RET( VtSeparableFilter(m_img[7], img, m_kern[4], m_kern[5]) );
            VT_HR_RET( VtSeparableFilter(m_img[8], img, m_kern[1], m_kern[3]) );
        }
    }
    else
    {
        if(m_ft==FilterEven || m_ft==FilterBoth)
        {
            VT_HR_RET( VtSeparableFilter(m_img[0], img, m_kern[0], m_kern[1]) );
            VT_HR_RET( VtSeparableFilter(m_img[1], img, m_kern[2], m_kern[3]) );
            VT_HR_RET( VtSeparableFilter(m_img[2], img, m_kern[4], m_kern[4]) );
            VT_HR_RET( VtSeparableFilter(m_img[3], img, m_kern[3], m_kern[2]) );
            VT_HR_RET( VtSeparableFilter(m_img[4], img, m_kern[1], m_kern[0]) );
        }

        if(m_ft==FilterOdd || m_ft==FilterBoth)
        {
            VT_HR_RET( VtSeparableFilter(m_img[5],  img, m_kern[5], m_kern[1]) );
            VT_HR_RET( VtSeparableFilter(m_img[6],  img, m_kern[6], m_kern[7]) );
            VT_HR_RET( VtSeparableFilter(m_img[7],  img, m_kern[8], m_kern[9]) );
            VT_HR_RET( VtSeparableFilter(m_img[8],  img, m_kern[9], m_kern[8]) );
            VT_HR_RET( VtSeparableFilter(m_img[9],  img, m_kern[7], m_kern[6]) );
            VT_HR_RET( VtSeparableFilter(m_img[10], img, m_kern[1], m_kern[5]) );
        }
    }
    return NOERROR;
}

HRESULT CSteerableFilter::GetImage(int iOrientation, FilterType eft, CFloatImg &imgRtn)
{
    if(iOrientation<0 || iOrientation>=(int)m_vvCoefs.size())
        return E_INVALIDARG;
    if(imgRtn.Width() != m_iWidth || imgRtn.Height() != m_iHeight)
        return E_INVALIDARG;

    if(iOrientation==0) // zero degrees
    {
        if(eft==FilterEven)
            m_img[0].CopyTo(imgRtn);
        else if(eft==FilterOdd)
            m_img[5].CopyTo(imgRtn);
        else
            return E_INVALIDARG;
    } 
    else if(iOrientation*2 == (int)m_vvCoefs.size()) // 90 degrees
    {
        if(eft==FilterEven)
            m_img[m_iLastEven].CopyTo(imgRtn);
        else if(eft==FilterOdd)
        {
            int x,y;
            for(y=0; y<imgRtn.Height(); y++)
            {
                float *pSrc = m_img[m_iLastOdd].Ptr(y);
                float *pDst = imgRtn.Ptr(y);
                for(x=0; x<imgRtn.Width(); x++)
                    pDst[x] = -pSrc[x];
            }
        }
        else
            return E_INVALIDARG;
    }

    return GetImage(&(m_vvCoefs[iOrientation]), eft, imgRtn); 
}

HRESULT CSteerableFilter::GetImage(CVecf *pCoefs, FilterType eft, CFloatImg &imgRtn)
{
    if(imgRtn.Width() != m_iWidth || imgRtn.Height() != m_iHeight)
        return E_INVALIDARG;

    int x,y;
    if(eft==FilterEven)
    {
        if(m_fo==SecondOrder)
        {
            for(y=0; y<imgRtn.Height(); y++)
            {
                float *pRtn = imgRtn.Ptr(y);
                float *pSrc0 = m_img[0].Ptr(y);
                float *pSrc1 = m_img[1].Ptr(y);
                float *pSrc2 = m_img[2].Ptr(y);
                for(x=0; x<imgRtn.Width(); x++)
                    pRtn[x] = pCoefs->El(0) * pSrc0[x] + pCoefs->El(1) * pSrc1[x] + pCoefs->El(2) * pSrc2[x];
            }
        }
        else
        {
            for(y=0; y<imgRtn.Height(); y++)
            {
                float *pRtn = imgRtn.Ptr(y);
                float *pSrc0 = m_img[0].Ptr(y);
                float *pSrc1 = m_img[1].Ptr(y);
                float *pSrc2 = m_img[2].Ptr(y);
                float *pSrc3 = m_img[3].Ptr(y);
                float *pSrc4 = m_img[4].Ptr(y);
                for(x=0; x<imgRtn.Width(); x++)
                    pRtn[x] = pCoefs->El(0) * pSrc0[x] + pCoefs->El(1) * pSrc1[x] + pCoefs->El(2) * pSrc2[x]
                        + pCoefs->El(3) * pSrc3[x] + pCoefs->El(4) * pSrc4[x];
            }
        }
    }
    else if(eft==FilterOdd)
    {
        if(m_fo==SecondOrder)
        {
            for(y=0; y<imgRtn.Height(); y++)
            {
                float *pRtn = imgRtn.Ptr(y);
                float *pSrc5 = m_img[5].Ptr(y);
                float *pSrc6 = m_img[6].Ptr(y);
                float *pSrc7 = m_img[7].Ptr(y);
                float *pSrc8 = m_img[8].Ptr(y);
                for(x=0; x<imgRtn.Width(); x++)
                    pRtn[x] = pCoefs->El(5) * pSrc5[x] + pCoefs->El(6) * pSrc6[x] + pCoefs->El(7) * pSrc7[x]
                        + pCoefs->El(8) * pSrc8[x];
            }
        }
        else
        {
            for(y=0; y<imgRtn.Height(); y++)
            {
                float *pRtn = imgRtn.Ptr(y);
                float *pSrc5 = m_img[5].Ptr(y);
                float *pSrc6 = m_img[6].Ptr(y);
                float *pSrc7 = m_img[7].Ptr(y);
                float *pSrc8 = m_img[8].Ptr(y);
                float *pSrc9 = m_img[9].Ptr(y);
                float *pSrc10 = m_img[10].Ptr(y);
                for(x=0; x<imgRtn.Width(); x++)
                    pRtn[x] = pCoefs->El(5) * pSrc5[x] + pCoefs->El(6) * pSrc6[x] + pCoefs->El(7) * pSrc7[x]
                        + pCoefs->El(8) * pSrc8[x] + pCoefs->El(9) * pSrc9[x] + pCoefs->El(10) * pSrc10[x];
            }
        }
    }
    else
        return E_INVALIDARG;
    return NOERROR;
}

// assumes second order and that both even and odd have been computed
void CSteerableFilter::GetLocalMagnitudeAndOrientation(int x, int y, float &fMag, float &fArg)
{
    float g2a = m_img[0].Pix(x,y);
    float g2b = m_img[1].Pix(x,y);
    float g2c = m_img[2].Pix(x,y);

    float h2a = m_img[5].Pix(x,y);
    float h2b = m_img[6].Pix(x,y);
    float h2c = m_img[7].Pix(x,y);
    float h2d = m_img[8].Pix(x,y);

    float c2 = 0.5f * (g2a * g2a - g2c * g2c) + 0.46875f * (h2a * h2a - h2d * h2d) +
        0.28125f * (h2b * h2b - h2c * h2c) + 0.1875f * (h2a * h2c - h2b * h2d);
    float c3 = -g2a * g2b - g2b * g2c - 0.9375f * (h2c * h2d + h2a * h2b) 
        - 1.6875f * h2b * h2c - 0.1875f * h2a * h2d;

    fMag = sqrt(c2*c2 + c3*c3);
    fArg = atan2(c3, c2)/2;
}


CVecf *CSteerableFilter::GetCoeffsByAngle(float fAngle)
{
    float c = cos(fAngle);
    float s = sin(fAngle);

    if(m_fo==SecondOrder)
    {
        m_v[0] = c*c;
        m_v[1] = -2*c*s;
        m_v[2] = s*s;
        m_v[5] = c*c*c;
        m_v[6] = -3*c*c*s;
        m_v[7] = 3*c*s*s;
        m_v[8] = -s*s*s;
    }
    else
    {
        float cc = c*c;
        float ss = s*s;
        m_v[0] = cc*cc;
        m_v[1] = -4*cc*c*s;
        m_v[2] = 6*cc*ss;
        m_v[3] = -4*c*s*ss;
        m_v[4] = ss*ss;
        m_v[5] = cc*cc*c;
        m_v[6] = -5*cc*cc*s;
        m_v[7] = 10*cc*c*ss;
        m_v[8] = -10*cc*s*ss;
        m_v[9] = 5*c*ss*ss;
        m_v[10] = -ss*ss*s;
    }
    return &m_v;
}
