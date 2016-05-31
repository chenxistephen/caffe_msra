//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      FFT routines
//
//  History:
//      2004/11/08-swinder
//			Created
//
//------------------------------------------------------------------------

// SPECIAL NOTE: (swinder 7/9/2010)
//
// ************************* FOR INTERNAL RESEARCH ONLY *************************
// 
// ************ DO NOT CALL THE ROUTINES BELOW IN ANY SHIPPING CODE *************
//
//

#pragma once

#include "vtcommon.h"
#include "vt_image.h"

namespace vt {

#define FFFT (1)
#define IFFT (-1)

void VtFindBestSizeForFFT(int &iW, int &iH, int iMaxExpand = 5);
void VtFindBestSizeForFFT1d(int &iSize, int iMaxExpand = 5);

// NOTE! it is importand to create the complex image using the following code:
// cOut.Create(width, height, alignAny)
// otherwise the FFT routine will not work with certain image widths
HRESULT VtForwardFFT(CComplexImg &cOut, const CFloatImg &cIn);

// NOTE! it is importand to create the complex image using the following code:
// cIn.Create(width, height, alignAny)
// otherwise the FFT routine will not work with certain image widths
// if retainsrc==false then input array is overwritten
HRESULT VtInverseFFT(CFloatImg &cOut, CComplexImg &cIn, bool bRetainSrc = false);

// for 1d fft total = size, pass = size, span = size
// for 2d fft (1) total = w*h, pass = w, span = w, (2) total = w*h, pass = h, span = w*h
HRESULT VtFFTFloat(Complexf *pSrc, int iTotal, int iPass, int iSpan, int iSign);

// dedicated 1d fft filters
HRESULT VtInverseFFT1d(CFloatImg &cOut, int iRowOut, CComplexImg &cIn, int iRowIn, bool bRetainSrc = false);

HRESULT VtFFTFilter1d(CFloatImg &cOut, int iRowOut, const CFloatImg &cIn, int iRowIn, Complexf (*callback)(float fW, void *p), void *pParam);
HRESULT VtFFTFilter1dGaussian(CFloatImg &cOut, int iRowOut, const CFloatImg &cIn, int iRowIn, double fS);
HRESULT VtFFTFilter1dDOG(CFloatImg &cOut, int iRowOut, const CFloatImg &cIn, int iRowIn, double fC, double fS, double k = 1.0);
HRESULT VtFFTFilter1dHilbert(CFloatImg &cOut, int iRowOut, const CFloatImg &cIn, int iRowIn);

// dedicated 2d fft filters
HRESULT VtFFTFilter2d(CFloatImg &cOut, const CFloatImg &cIn, Complexf (*callback)(float fWx, float fWy, void *p), void *pParam, bool bRetainSrc = false);

// Gaussian smooth - kernel is k exp(-|r/sigma|^2)
HRESULT VtFFTFilter2dGaussian(CFloatImg &cOut, const CFloatImg &cIn, double fS);

// exponential smooth - kernel is k exp(-|r/sigma|). sigma is radius at y=e^(-1) (=0.3679)
HRESULT VtFFTFilter2dExponential(CFloatImg &cOut, const CFloatImg &cIn, double fS);

// filter is a Gaussian along the orientation axis and a difference-of-Gaussians orthogonal to it
// returns the result (even-phase) and the same result passed through a hilbert transform (odd-phase)
// fCenter must be less than fSurround, the frequency peak is given unit magnitude
HRESULT VtFFTFilter2dOrientedDOG(CFloatImg &cOutEven, CFloatImg &cOutOdd, const CFloatImg &cIn, 
                                 double angle, double fCenter, double fSurround, double fLongAxis);

};
