//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for steerable filters, based on:
//      W. T. Freeman and E. H. Adelson
//      The design and use of steerable filters
//      IEEE Trans. Pattern Analysis and Machine Intelligence volume 13, 891--906, 1991
//
//  History:
//      2008/4/4-swinder
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"

namespace vt {

enum FilterOrder { SecondOrder, FourthOrder };
enum FilterType { FilterEven, FilterOdd, FilterBoth };

class CSteerableFilter
{
public:
    // If iOrientationCount is specified this lets the class pre-compute coefficients
    // for the specific orientations that you need, otherwise it will have to
    // compute them on the fly.
    // Scale is relative to the example filters in the PAMI paper
    // so that a scale of 1 gives the coefficients in the paper. This is the smallest filter
    // that can be generated, so scale cannot be less than 1.
    HRESULT Create(int iWidth, int iHeight, FilterType eft, FilterOrder efo, float fScale = 1, int iOrientationCount = 0);
    HRESULT Update(const CFloatImg &img);
    void Deallocate(); // free up space
    
    // Get the filter response for a particular orientation channel.
    // Filter type can be FilterEven or FilterOdd
    HRESULT GetImage(int iOrientation, FilterType eft, CFloatImg &imgRtn);

    // Get the filter response for a particular angle
    HRESULT GetImage(float fAngle, FilterType eft, CFloatImg &imgRtn)
        { return GetImage(GetCoeffsByAngle(fAngle), eft, imgRtn); }

    // get pixels on the fly (don't compute a whole image)
    float GetPixelOddFilter(int x, int y, int iOrientation)
        { return GetPixelOddFilter(x, y, &(m_vvCoefs[iOrientation])); }

    float GetPixelEvenFilter(int x, int y, int iOrientation)
        { return GetPixelEvenFilter(x, y, &(m_vvCoefs[iOrientation])); }

    float GetPixelOddFilter(int x, int y, float fAngle)
        { return GetPixelOddFilter(x, y, GetCoeffsByAngle(fAngle)); }

    float GetPixelEvenFilter(int x, int y, float fAngle)
        { return GetPixelEvenFilter(x, y, GetCoeffsByAngle(fAngle)); }

    // the following assumes second order and that both even and odd have been computed
    // assumes only one local edge orientation
    // returns the magnitude of the filter responses and the optimal angle (in radians)
    void GetLocalMagnitudeAndOrientation(int x, int y, float &fMag, float &fArg);

    // faster special cases
    float GetPixelOddFilter_0_Degrees(int x, int y) 
        { return m_img[5].Pix(x,y); }
    float GetPixelEvenFilter_0_Degrees(int x, int y)
        { return m_img[0].Pix(x,y); }
    float GetPixelOddFilter_90_Degrees(int x, int y)
        { return -m_img[m_iLastOdd].Pix(x,y); }
    float GetPixelEvenFilter_90_Degrees(int x, int y)
        { return m_img[m_iLastEven].Pix(x,y); }

private:
    CVecf *GetCoeffsByAngle(float fAngle);

    float GetPixelOddFilter(int x, int y, CVecf *pCoefs)
    {
        float f = pCoefs->El(5) * m_img[5](x,y) + pCoefs->El(6) * m_img[6](x,y)
            + pCoefs->El(7) * m_img[7](x,y) + pCoefs->El(8) * m_img[8](x,y);
        if(m_fo==FourthOrder)
            f += pCoefs->El(9) * m_img[9](x,y) + pCoefs->El(10) * m_img[10](x,y);
        return f;
    }

    float GetPixelEvenFilter(int x, int y, CVecf *pCoefs)
    {
        float f = pCoefs->El(0) * m_img[0](x,y) + pCoefs->El(1) * m_img[1](x,y)
            + pCoefs->El(2) * m_img[2](x,y);
        if(m_fo==FourthOrder)
            f += pCoefs->El(3) * m_img[3](x,y) + pCoefs->El(4) * m_img[4](x,y);
        return f;
    }

    HRESULT GetImage(CVecf *pCoefs, FilterType eft, CFloatImg &imgRtn);

    FilterType  m_ft;
    FilterOrder m_fo;
    int m_iWidth;
    int m_iHeight;
    C1dKernel m_kern[10];
    CFloatImg m_img[11];
    vector<CVecf> m_vvCoefs; // one for each orientation
    CVecf m_v;
    int m_iLastEven, m_iLastOdd;
};

};
