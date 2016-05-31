//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      FFT routines
//
//  the main mixed radix fft routing was adapted from fortran code written
//  by r. c. singleton, stanford research institute, sept. 1968
//
// R.C. Singleton, An Algorithm for Computing the Mixed Radix F.F.T., IEEE Trans. Audio Electroacoust.,
// AU-1(1969) 93-107. 
// Reprinted in: L.R. Rabiner and C.M. Rader: Digital Signal Processing, IEEE Press New York (1972) 294. 
// also:
// Singleton, R. C. (1979) Mixed Radix Fast Fourier Transforms, in Programs for Digital Signal Processing, 
// IEEE Digital Signal Processing Committee eds. IEEE Press. 
// 
//  History:
//      2004/11/08-swinder
//			Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_image.h"
#include "vt_fft.h"

using namespace vt;

bool AssertStrideLengthForFFT(CComplexImg& img)
{
	const bool strideOK = img.StrideBytes() == img.Width() * img.PixSize();
	
	VT_ASSERT(strideOK);

	return strideOK;
}

static Int64 FFTSizeCost(int iS1, int iS2)
{
    static int rgiPrimes[] = {2,3,5,7,11,13,17,19,23};

    Int64 iCost = 0;
    int i;
    int iS = iS2;
    for(i=0;i<9;i++)
    {
        int iPrime = rgiPrimes[i];
        while((iS/iPrime)*iPrime==iS)
        {
            iS/=iPrime;
            iCost += iPrime;
        }
    }
    if(iS!=1)
    {
        iCost += iS;
    }

    return (Int64)iS1 * (Int64)iS2 * iCost;
}

void vt::VtFindBestSizeForFFT(int &iW, int &iH, int iMaxExpand)
{
    if(iW<1 || iH<1)
        return;

    int iX, iY;
    Int64 iMinCost = -1;
    Int64 iStartCost = -1;
    int iMinX = 0;
    int iMinY = 0;
    for(iY=iH; iY<iH + iMaxExpand; iY++)
        for(iX=iW; iX<iW + iMaxExpand; iX++)
        {
            Int64 iCost = FFTSizeCost(iX, iY) + FFTSizeCost(iY, iX);
            if(iStartCost<0)
               iStartCost = iCost;
                 
            if(iMinCost<0 || iCost<iMinCost)
            {
                iMinX = iX;
                iMinY = iY;
                iMinCost = iCost;
            }
        }

    if(iStartCost - iMinCost > iStartCost / 20) // gain of 5%
    {
        iW = iMinX;
        iH = iMinY;
    }
}

void vt::VtFindBestSizeForFFT1d(int &iSize, int iMaxExpand)
{
    if(iSize<1)
        return;

    int iX;
    Int64 iMinCost = -1;
    Int64 iStartCost = -1;
    int iMinX = 0;
    for(iX=iSize; iX<iSize + iMaxExpand; iX++)
    {
        Int64 iCost = FFTSizeCost(1, iX);
        if(iStartCost<0)
            iStartCost = iCost;
        if(iMinCost<0 || iCost<iMinCost)
        {
            iMinX = iX;
            iMinCost = iCost;
        }
    }

    if(iStartCost - iMinCost > iStartCost / 20) // gain of 5%
        iSize = iMinX;
}

#define SIN_60	0.86602540378443864676372317075294
#define COS_72	0.30901699437494742410229341718282
#define SIN_72	0.95105651629515357211643933337938

#define MAXNFACTORS		30

HRESULT vt::VtForwardFFT(CComplexImg &cOut, const CFloatImg &cIn)
{
	if (!AssertStrideLengthForFFT(cOut))
		return E_INVALIDARG;

	if(!cIn.IsValid() || !cOut.IsValid() || cIn.Width()!=cOut.Width() || cIn.Height()!=cOut.Height())
		return E_INVALIDARG;

	int x, y;
	for(y=0; y<cIn.Height(); y++)
	{
		const float *pf = cIn.Ptr(y);
		Complexf *pc = cOut.Ptr(y);
		
		x=0;

#if (defined(_M_IX86) || defined(_M_AMD64))
		if(g_SupportSSE1() && cIn.Width()>=4)
		{
			__m128 x0 = _mm_setzero_ps();

			for(; x<cIn.Width()-3; x+=4)
			{
				__m128 x1 = _mm_loadu_ps(pf + x);
				__m128 x2 = _mm_unpacklo_ps(x1, x0);
				__m128 x3 = _mm_unpackhi_ps(x1, x0);
				_mm_storeu_ps((float *)(pc + x), x2);
				_mm_storeu_ps((float *)(pc + x + 2), x3);
			}
		}
#endif
		for(; x<cIn.Width(); x++)
		{
			pc[x].Re = pf[x];
			pc[x].Im = 0;
		}
	}

	int iTotal = cIn.Width() * cIn.Height();
	HRESULT hr = NOERROR;

	if(cIn.Width()==1 || cIn.Height()==1)
	{
		hr = VtFFTFloat(cOut.Ptr(), iTotal, iTotal, iTotal, FFFT);
	}
	else
	{
		hr = VtFFTFloat(cOut.Ptr(), iTotal, cIn.Width(), cIn.Width(), FFFT);
		if(FAILED(hr))
			return hr;

		hr = VtFFTFloat(cOut.Ptr(), iTotal, cIn.Height(), iTotal, FFFT);
	}

	return hr;
}

HRESULT vt::VtInverseFFT(CFloatImg &cOut, CComplexImg &cIn, bool bRetainSrc)
{
	if (!AssertStrideLengthForFFT(cIn))
		return E_INVALIDARG;

	if(!cIn.IsValid() || !cOut.IsValid() || cIn.Width()!=cOut.Width() || cIn.Height()!=cOut.Height())
		return E_INVALIDARG;

	HRESULT hr = NOERROR;
	CComplexImg *pcIn = &cIn;
	CComplexImg *pcNew = NULL;

	if(bRetainSrc)
	{
		pcNew = VT_NOTHROWNEW CComplexImg(cIn.Width(), cIn.Height(), alignAny);
		if(pcNew==NULL)
		{
			hr = E_OUTOFMEMORY;
			goto exit_label;
		}

		cIn.CopyTo(*pcNew);
		pcIn = pcNew;
	}

	int iTotal = cIn.Width() * cIn.Height();

	if(cIn.Width()==1 || cIn.Height()==1)
	{
		hr = VtFFTFloat(pcIn->Ptr(), iTotal, iTotal, iTotal, IFFT);
	}
	else
	{
		hr = VtFFTFloat(pcIn->Ptr(), iTotal, cIn.Width(), cIn.Width(), IFFT);
		if(FAILED(hr))
			goto exit_label;

		hr = VtFFTFloat(pcIn->Ptr(), iTotal, cIn.Height(), iTotal, IFFT);
	}

	if(FAILED(hr))
		goto exit_label;

	int x, y;
	float fNorm = 1.0f/iTotal;
	for(y=0; y<cIn.Height(); y++)
	{
		float *pf = cOut.Ptr(y);
		Complexf *pc = pcIn->Ptr(y);

		x = 0;
#if (defined(_M_IX86) || defined(_M_AMD64))
		if(g_SupportSSE1() && cIn.Width()>=4)
		{
			__m128 x0 = _mm_set1_ps(fNorm);

			for(; x<cIn.Width()-3; x+=4)
			{
				__m128 x1 = _mm_loadu_ps((float *)(pc + x));
				__m128 x2 = _mm_loadu_ps((float *)(pc + x + 2));
				__m128 x3 = _mm_shuffle_ps(x1, x2, _MM_SHUFFLE(2, 0, 2, 0));
				__m128 x4 = _mm_mul_ps(x3, x0);
				_mm_storeu_ps(pf + x, x4);
			}
		}
#endif

		for(; x<cIn.Width(); x++)
		{
			pf[x] = pc[x].Re * fNorm;
		}
	}

exit_label:
	if(pcNew!=NULL)
		delete pcNew;
	return hr;
}

HRESULT vt::VtInverseFFT1d(CFloatImg &cOut, int iRowOut, CComplexImg &cIn, int iRowIn, bool bRetainSrc)
{
	if (!AssertStrideLengthForFFT(cIn))
		return E_INVALIDARG;

	if(cIn.Width() != cOut.Width() || iRowIn<0 || iRowIn>=cIn.Height() || iRowOut<0 || iRowOut>=cOut.Height())
		return E_INVALIDARG;

	int iW = cIn.Width();
	HRESULT hr = NOERROR;
	float fNorm = 1.0f/iW;
	int iX;

	if(bRetainSrc)
	{
		CComplexImg cCpx;
		hr = cCpx.Create(iW, 1);
		if(FAILED(hr))
			return hr;
		memcpy(cCpx.Ptr(), cIn.Ptr(iRowIn), iW * sizeof(Complexf));
		hr = VtFFTFloat(cCpx.Ptr(), iW, iW, iW, IFFT);
		if(FAILED(hr))
			return hr;
		for(iX=0; iX<iW; iX++)
			cOut.Pix(iX, iRowOut) = fNorm * cCpx.Pix(iX, 0).Re;
	}
	else
	{
		hr = VtFFTFloat(cIn.Ptr(iRowIn), iW, iW, iW, IFFT);
		if(FAILED(hr))
			return hr;
		for(iX=0; iX<iW; iX++)
			cOut.Pix(iX, iRowOut) = fNorm * cIn.Pix(iX, iRowIn).Re;
	}

	return hr;
}

HRESULT vt::VtFFTFloat(Complexf *pSrc, int iTotal, int iPass, int iSpan, int iSign)
{
	HRESULT hr = NOERROR;

	pSrc--; // original code was fortran
	if(iPass<2)
		return E_INVALIDARG;
	
	int rgiFactor[MAXNFACTORS];
	
	// determine the factors of iPass
	int iFacCount = 0;
	int k = iPass;
	while (k % 16 == 0) {
		rgiFactor [iFacCount++] = 4;
		if(iFacCount>=MAXNFACTORS)
			return E_TOOCOMPLEX;
		k /= 16;
	}
	int j = 3;
	int jj = 9;
	do {
		while (k % jj == 0) {
			rgiFactor [iFacCount++] = j;
			if(iFacCount>=MAXNFACTORS)
				return E_TOOCOMPLEX;
			k /= jj;
		}
		j += 2;
		jj = j * j;
	} while (jj <= k);
	int kt;
	if (k <= 4) {
		kt = iFacCount;
		rgiFactor [iFacCount] = k;
		if (k != 1)
		{
			iFacCount++;
			if(iFacCount>=MAXNFACTORS)
				return E_TOOCOMPLEX;
		}
	} else {
		if (k - (k / 4 << 2) == 0) {
			rgiFactor [iFacCount++] = 2;
			if(iFacCount>=MAXNFACTORS)
				return E_TOOCOMPLEX;
			k /= 4;
		}
		kt = iFacCount;
		j = 2;
		do {
			if (k % j == 0) {
				rgiFactor [iFacCount++] = j;
				if(iFacCount>=MAXNFACTORS)
					return E_TOOCOMPLEX;
				k /= j;
			}
			j = ((j + 1) / 2 << 1) + 1;
		} while (j <= k);
	}
	if (kt) {
		j = kt;
		do {
			rgiFactor [iFacCount++] = rgiFactor [--j];
			if(iFacCount>=MAXNFACTORS)
				return E_TOOCOMPLEX;
		} while (j);
	}
	
	float fS60 = (float)SIN_60;
	float fC72 = (float)COS_72;
	float fS72 = (float)SIN_72;
	float fPi2 = (float)VT_PI;
	
	int inc = iSign;
	if (iSign < 0) {
		fS72 = -fS72;
		fS60 = -fS60;
		fPi2 = -fPi2;
		inc = -inc;		// absolute value
	}
	int nt = inc * iTotal;
	int ns = inc * iSpan;
	int kspan = ns;
	int nn = nt - inc;
	int jc = ns / iPass;
	float radf = fPi2 * jc;
	fPi2 *= 2.0;	// use 2 PI from here on
	int ii = 0;
	int jf = 0;
	
	int k1, k2, k3 = 0, k4, kk;
	float aa, aj, ak, ajm, ajp, akm, akp;
	float bb, bj, bk, bjm, bjp, bkm, bkp;

	int iMaxFactors = iPass;
	int iMaxPerm = iPass;
	float *pfReTmp = VT_NOTHROWNEW float [iPass];
	float *pfImTmp = VT_NOTHROWNEW float [iPass];
	float *pfCos = VT_NOTHROWNEW float [iPass];
	float *pfSin = VT_NOTHROWNEW float [iPass];
	int *piPerm = VT_NOTHROWNEW int [iPass];
	if(piPerm==NULL || pfSin==NULL || pfCos==NULL || pfImTmp==NULL || pfReTmp==NULL)
		goto Memory_Error_Label;
	

	/* compute fourier transform */
	for (;;)
	{
		float c1;
		float c2 = 0;
		float c3 = 0;
		float s1;
		float s2 = 0;
		float s3 = 0;
		
		float sd = radf / kspan;
		float cd = (float)sin(sd);
		cd = (float)(2.0 * cd * cd);
		sd = (float)(sin(sd + sd));
		kk = 1;
		ii++;
		
		switch (rgiFactor [ii - 1])
		{
		case 2:
			/* transform for factor of 2 (including rotation factor) */
			kspan /= 2;
			k1 = kspan + 2;
			do {
				do {
					k2 = kk + kspan;
					ak = pSrc[k2].Re;
					bk = pSrc[k2].Im;
					pSrc[k2].Re = pSrc[kk].Re - ak;
					pSrc[k2].Im = pSrc[kk].Im - bk;
					pSrc[kk].Re += ak;
					pSrc[kk].Im += bk;
					kk = k2 + kspan;
				} while (kk <= nn);
				kk -= nn;
			} while (kk <= jc);
			if (kk > kspan)
				goto Permute_Results_Label;		/* exit infinite loop */
			do {
				c1 = 1.0f - cd;
				s1 = sd;
				do {
					do {
						do {
							k2 = kk + kspan;
							ak = pSrc[kk].Re - pSrc[k2].Re;
							bk = pSrc[kk].Im - pSrc[k2].Im;
							pSrc[kk].Re += pSrc[k2].Re;
							pSrc[kk].Im += pSrc[k2].Im;
							pSrc[k2].Re = c1 * ak - s1 * bk;
							pSrc[k2].Im = s1 * ak + c1 * bk;
							kk = k2 + kspan;
						} while (kk < nt);
						k2 = kk - nt;
						c1 = -c1;
						kk = k1 - k2;
					} while (kk > k2);
					ak = c1 - (cd * c1 + sd * s1);
					s1 = sd * c1 - cd * s1 + s1;
					c1 = 2.0f - (ak * ak + s1 * s1);
					s1 *= c1;
					c1 *= ak;
					kk += jc;
				} while (kk < k2);
				k1 += inc + inc;
				kk = (k1 - kspan) / 2 + jc;
			} while (kk <= jc + jc);
			break;
			
		case 4:			/* transform for factor of 4 */
			iSpan = kspan;
			kspan /= 4;
			
			do {
				c1 = 1.0;
				s1 = 0.0;
				do {
					do {
						k1 = kk + kspan;
						k2 = k1 + kspan;
						k3 = k2 + kspan;
						akp = pSrc[kk].Re + pSrc[k2].Re;
						akm = pSrc[kk].Re - pSrc[k2].Re;
						ajp = pSrc[k1].Re + pSrc[k3].Re;
						ajm = pSrc[k1].Re - pSrc[k3].Re;
						bkp = pSrc[kk].Im + pSrc[k2].Im;
						bkm = pSrc[kk].Im - pSrc[k2].Im;
						bjp = pSrc[k1].Im + pSrc[k3].Im;
						bjm = pSrc[k1].Im - pSrc[k3].Im;
						pSrc[kk].Re = akp + ajp;
						pSrc[kk].Im = bkp + bjp;
						ajp = akp - ajp;
						bjp = bkp - bjp;
						if (iSign < 0) 
						{
							akp = akm + bjm;
							bkp = bkm - ajm;
							akm -= bjm;
							bkm += ajm;
						}
						else 
						{
							akp = akm - bjm;
							bkp = bkm + ajm;
							akm += bjm;
							bkm -= ajm;
						}
						/* avoid useless multiplies */
						if (s1 == 0.0) 
						{
							pSrc[k1].Re = akp;
							pSrc[k2].Re = ajp;
							pSrc[k3].Re = akm;
							pSrc[k1].Im = bkp;
							pSrc[k2].Im = bjp;
							pSrc[k3].Im = bkm;
						} 
						else 
						{
							pSrc[k1].Re = akp * c1 - bkp * s1;
							pSrc[k2].Re = ajp * c2 - bjp * s2;
							pSrc[k3].Re = akm * c3 - bkm * s3;
							pSrc[k1].Im = akp * s1 + bkp * c1;
							pSrc[k2].Im = ajp * s2 + bjp * c2;
							pSrc[k3].Im = akm * s3 + bkm * c3;
						}
						kk = k3 + kspan;
					} while (kk <= nt);
					
					c2 = c1 - (cd * c1 + sd * s1);
					s1 = sd * c1 - cd * s1 + s1;
					c1 = 2.0f - (c2 * c2 + s1 * s1);
					s1 *= c1;
					c1 *= c2;
					/* values of c2, c3, s2, s3 that will get used next time */
					c2 = c1 * c1 - s1 * s1;
					s2 = 2.0f * c1 * s1;
					c3 = c2 * c1 - s2 * s1;
					s3 = c2 * s1 + s2 * c1;
					kk = kk - nt + jc;
				} while (kk <= kspan);
				kk = kk - kspan + inc;
			} while (kk <= jc);
			if (kspan == jc)
				goto Permute_Results_Label;		/* exit infinite loop */
			break;
			
		default:
			/*  transform for odd factors */
			k = rgiFactor [ii - 1];
			iSpan = kspan;
			kspan /= k;
			
			switch (k) 
			{
			case 3:	/* transform for factor of 3 (optional code) */
				do {
					do {
						k1 = kk + kspan;
						k2 = k1 + kspan;
						ak = pSrc[kk].Re;
						bk = pSrc[kk].Im;
						aj = pSrc[k1].Re + pSrc[k2].Re;
						bj = pSrc[k1].Im + pSrc[k2].Im;
						pSrc[kk].Re = ak + aj;
						pSrc[kk].Im = bk + bj;
						ak -= 0.5f * aj;
						bk -= 0.5f * bj;
						aj = (pSrc[k1].Re - pSrc[k2].Re) * fS60;
						bj = (pSrc[k1].Im - pSrc[k2].Im) * fS60;
						pSrc[k1].Re = ak - bj;
						pSrc[k2].Re = ak + bj;
						pSrc[k1].Im = bk + aj;
						pSrc[k2].Im = bk - aj;
						kk = k2 + kspan;
					} while (kk < nn);
					kk -= nn;
				} while (kk <= kspan);
				break;
				
			case 5:	/*  transform for factor of 5 (optional code) */
				c2 = fC72 * fC72 - fS72 * fS72;
				s2 = 2.0f * fC72 * fS72;
				do {
					do {
						k1 = kk + kspan;
						k2 = k1 + kspan;
						k3 = k2 + kspan;
						k4 = k3 + kspan;
						akp = pSrc[k1].Re + pSrc[k4].Re;
						akm = pSrc[k1].Re - pSrc[k4].Re;
						bkp = pSrc[k1].Im + pSrc[k4].Im;
						bkm = pSrc[k1].Im - pSrc[k4].Im;
						ajp = pSrc[k2].Re + pSrc[k3].Re;
						ajm = pSrc[k2].Re - pSrc[k3].Re;
						bjp = pSrc[k2].Im + pSrc[k3].Im;
						bjm = pSrc[k2].Im - pSrc[k3].Im;
						aa = pSrc[kk].Re;
						bb = pSrc[kk].Im;
						pSrc[kk].Re = aa + akp + ajp;
						pSrc[kk].Im = bb + bkp + bjp;
						ak = akp * fC72 + ajp * c2 + aa;
						bk = bkp * fC72 + bjp * c2 + bb;
						aj = akm * fS72 + ajm * s2;
						bj = bkm * fS72 + bjm * s2;
						pSrc[k1].Re = ak - bj;
						pSrc[k4].Re = ak + bj;
						pSrc[k1].Im = bk + aj;
						pSrc[k4].Im = bk - aj;
						ak = akp * c2 + ajp * fC72 + aa;
						bk = bkp * c2 + bjp * fC72 + bb;
						aj = akm * s2 - ajm * fS72;
						bj = bkm * s2 - bjm * fS72;
						pSrc[k2].Re = ak - bj;
						pSrc[k3].Re = ak + bj;
						pSrc[k2].Im = bk + aj;
						pSrc[k3].Im = bk - aj;
						kk = k4 + kspan;
					} while (kk < nn);
					kk -= nn;
				} while (kk <= kspan);
				break;
				
			default:
				if (k != jf) 
				{
					jf = k;
					s1 = fPi2 / k;
					c1 = (float)cos(s1);
					s1 = (float)sin(s1);
					if (jf > iMaxFactors)
						goto Error_Label;
					pfCos [jf - 1] = 1.0;
					pfSin [jf - 1] = 0.0;
					j = 1;
					do {
						pfCos [j - 1] = pfCos [k - 1] * c1 + pfSin [k - 1] * s1;
						pfSin [j - 1] = pfCos [k - 1] * s1 - pfSin [k - 1] * c1;
						k--;
						pfCos [k - 1] = pfCos [j - 1];
						pfSin [k - 1] = -pfSin [j - 1];
						j++;
					} while (j < k);
				}
				do {
					do {
						k1 = kk;
						k2 = kk + iSpan;
						ak = aa = pSrc[kk].Re;
						bk = bb = pSrc[kk].Im;
						j = 1;
						k1 += kspan;
						do {
							k2 -= kspan;
							j++;
							pfReTmp [j - 1] = pSrc[k1].Re + pSrc[k2].Re;
							ak += pfReTmp [j - 1];
							pfImTmp [j - 1] = pSrc[k1].Im+ pSrc[k2].Im;
							bk += pfImTmp [j - 1];
							j++;
							pfReTmp [j - 1] = pSrc[k1].Re - pSrc[k2].Re;
							pfImTmp [j - 1] = pSrc[k1].Im - pSrc[k2].Im;
							k1 += kspan;
						} while (k1 < k2);
						pSrc[kk].Re = ak;
						pSrc[kk].Im = bk;
						k1 = kk;
						k2 = kk + iSpan;
						j = 1;
						do {
							k1 += kspan;
							k2 -= kspan;
							jj = j;
							ak = aa;
							bk = bb;
							aj = 0.0;
							bj = 0.0;
							k = 1;
							do {
								k++;
								ak += pfReTmp [k - 1] * pfCos [jj - 1];
								bk += pfImTmp [k - 1] * pfCos [jj - 1];
								k++;
								aj += pfReTmp [k - 1] * pfSin [jj - 1];
								bj += pfImTmp [k - 1] * pfSin [jj - 1];
								jj += j;
								if (jj > jf) 
								{
									jj -= jf;
								}
							} while (k < jf);
							k = jf - j;
							pSrc[k1].Re = ak - bj;
							pSrc[k1].Im = bk + aj;
							pSrc[k2].Re = ak + bj;
							pSrc[k2].Im = bk - aj;
							j++;
						} while (j < k);
						kk += iSpan;
					} while (kk <= nn);
					kk -= nn;
				} while (kk <= kspan);
				break;
			}
			/*  multiply by rotation factor (except for factors of 2 and 4) */
			if (ii == iFacCount)
				goto Permute_Results_Label;		/* exit infinite loop */
			kk = jc + 1;
			do {
				c2 = 1.0f - cd;
				s1 = sd;
				do {
					c1 = c2;
					s2 = s1;
					kk += kspan;
					do {
						do {
							ak = pSrc[kk].Re;
							pSrc[kk].Re = c2 * ak - s2 * pSrc[kk].Im;
							pSrc[kk].Im = s2 * ak + c2 * pSrc[kk].Im;
							kk += iSpan;
						} while (kk <= nt);
						ak = s1 * s2;
						s2 = s1 * c2 + c1 * s2;
						c2 = c1 * c2 - ak;
						kk = kk - nt + kspan;
					} while (kk <= iSpan);
					c2 = c1 - (cd * c1 + sd * s1);
					s1 += sd * c1 - cd * s1;
					c1 = 2.0f - (c2 * c2 + s1 * s1);
					s1 *= c1;
					c2 *= c1;
					kk = kk - iSpan + jc;
				} while (kk <= kspan);
				kk = kk - kspan + jc + inc;
			} while (kk <= jc + jc);
			break;
		}
	}
   
   /*  permute the results to normal order---done in two stages */
   /*  permutation for square factors of n */
Permute_Results_Label:
	piPerm [0] = ns;
	if (kt) {
		k = kt + kt + 1;
		if (iFacCount < k)
			k--;
		j = 1;
		piPerm [k] = jc;
		do {
			piPerm [j] = piPerm [j - 1] / rgiFactor [j - 1];
			piPerm [k - 1] = piPerm [k] * rgiFactor [j - 1];
			j++;
			k--;
		} while (j < k);
		k3 = piPerm [k];
		kspan = piPerm [1];
		kk = jc + 1;
		k2 = kspan + 1;
		j = 1;
		if (iPass != iTotal) {
			/*  permutation for multivariate transform */
Permute_Multi_Label:
			do {
				do {
					k = kk + jc;
					do {
						/* swap pSrc[kk].Re <> pSrc[k2].Re, pSrc[kk].Im <> pSrc[k2].Im */
						ak = pSrc[kk].Re; pSrc[kk].Re = pSrc[k2].Re; pSrc[k2].Re = ak;
						bk = pSrc[kk].Im; pSrc[kk].Im = pSrc[k2].Im; pSrc[k2].Im = bk;
						kk += inc;
						k2 += inc;
					} while (kk < k);
					kk += ns - jc;
					k2 += ns - jc;
				} while (kk < nt);
				k2 = k2 - nt + kspan;
				kk = kk - nt + jc;
			} while (k2 < ns);
			do {
				do {
					k2 -= piPerm [j - 1];
					j++;
					k2 = piPerm [j] + k2;
				} while (k2 > piPerm [j - 1]);
				j = 1;
				do {
					if (kk < k2)
						goto Permute_Multi_Label;
					kk += jc;
					k2 += kspan;
				} while (k2 < ns);
			} while (kk < ns);
		} else {
			//  permutation for single-variate transform (optional code)
Permute_Single_Label:
			do {
				/* swap pSrc[kk].Re <> pSrc[k2].Re, pSrc[kk].Im <> pSrc[k2].Im */
				ak = pSrc[kk].Re; pSrc[kk].Re = pSrc[k2].Re; pSrc[k2].Re = ak;
				bk = pSrc[kk].Im; pSrc[kk].Im = pSrc[k2].Im; pSrc[k2].Im = bk;
				kk += inc;
				k2 += kspan;
			} while (k2 < ns);
			do {
				do {
					k2 -= piPerm [j - 1];
					j++;
					k2 = piPerm [j] + k2;
				} while (k2 > piPerm [j - 1]);
				j = 1;
				do {
					if (kk < k2)
						goto Permute_Single_Label;
					kk += inc;
					k2 += kspan;
				} while (k2 < ns);
			} while (kk < ns);
		}
		jc = k3;
	}
	
	if ((kt << 1) + 1 >= iFacCount)
	{
		//printf("error?\n");
		goto OK_Label;
	}
	iSpan = piPerm [kt];
	/* permutation for square-free factors of n */
	j = iFacCount - kt;
	rgiFactor [j] = 1;
	do {
		rgiFactor [j - 1] *= rgiFactor [j];
		j--;
	} while (j != kt);
	kt++;
	nn = rgiFactor [kt - 1] - 1;
	if (nn > iMaxPerm)
		goto Error_Label;
	j = jj = 0;

	for (;;)
	{
		k = kt + 1;
		k2 = rgiFactor [kt - 1];
		kk = rgiFactor [k - 1];
		j++;
		if (j > nn)
			break;				/* exit infinite loop */
		jj += kk;
		while (jj >= k2) {
			jj -= k2;
			k2 = kk;
			k++;
			kk = rgiFactor [k - 1];
			jj += kk;
		}
		piPerm [j - 1] = jj;
	}
	/*  determine the permutation cycles of length greater than 1 */
	j = 0;
	
	for (;;) 
	{
		do {
			j++;
			kk = piPerm [j - 1];
		} while (kk < 0);
		if (kk != j) {
			do {
				k = kk;
				kk = piPerm [k - 1];
				piPerm [k - 1] = -kk;
			} while (kk != j);
			k3 = kk;
		} else {
			piPerm [j - 1] = -j;
			if (j == nn)
				break;		/* exit infinite loop */
		}
	}
	iMaxFactors *= inc;
	/*  reorder a and b, following the permutation cycles */
	
	for (;;)
	{
		j = k3 + 1;
		nt -= iSpan;
		ii = nt - inc + 1;
		if (nt < 0)
			break;			/* exit infinite loop */
		do {
			do {
				j--;
			} while (piPerm [j - 1] < 0);
			jj = jc;
			do {
				kspan = jj;
				if (jj > iMaxFactors) {
					kspan = iMaxFactors;
				}
				jj -= kspan;
				k = piPerm [j - 1];
				kk = jc * k + ii + jj;
				k1 = kk + kspan;
				k2 = 0;
				do {
					k2++;
					pfReTmp [k2 - 1] = pSrc[k1].Re;
					pfImTmp [k2 - 1] = pSrc[k1].Im;
					k1 -= inc;
				} while (k1 != kk);
				do {
					k1 = kk + kspan;
					k2 = k1 - jc * (k + piPerm [k - 1]);
					k = -piPerm [k - 1];
					do {
						pSrc[k1].Re = pSrc[k2].Re;
						pSrc[k1].Im = pSrc[k2].Im;
						k1 -= inc;
						k2 -= inc;
					} while (k1 != kk);
					kk = k2;
				} while (k != j);
				k1 = kk + kspan;
				k2 = 0;
				do {
					k2++;
					pSrc[k1].Re = pfReTmp [k2 - 1];
					pSrc[k1].Im = pfImTmp [k2 - 1];
					k1 -= inc;
				} while (k1 != kk);
			} while (jj);
		} while (j != 1);
	}

	goto OK_Label;
Error_Label:
	hr = E_TOOCOMPLEX;
	goto OK_Label;
Memory_Error_Label:
	hr = E_OUTOFMEMORY;
OK_Label:
	if(pfReTmp)
		delete [] pfReTmp;
	if(pfImTmp)
		delete [] pfImTmp;
	if(pfCos)
		delete [] pfCos;
	if(pfSin)
		delete [] pfSin;
	if(piPerm)
		delete [] piPerm;
	return hr;
}

HRESULT vt::VtFFTFilter1d(CFloatImg &cOut, int iRowOut, const CFloatImg &cIn, int iRowIn, Complexf (*pcallback)(float fW, void *p), void *pParam)
{
	if(!cIn.IsValid() || !cOut.IsValid() ||
		iRowIn<0 || iRowIn>=cIn.Height() || iRowOut<0 || iRowOut>=cOut.Height() ||
		cIn.Width() != cOut.Width())
		return E_INVALIDARG;
	int iW = cIn.Width();
	float f2PW = (float)( VT_PI * 2.0 / iW );
	int iW2 = iW / 2 + 1;
	float fNorm = 1.0f / iW;
	HRESULT hr = NOERROR;

	if(cIn.Bands()==1 && cOut.Bands()==1)
	{
		// float to float
		CComplexImg cCpx;
		hr = cCpx.Create(iW, 1);
		if(FAILED(hr))
			return hr;

		int iX;
		for(iX=0; iX<iW; iX++)
		{
			cCpx(iX, 0).Re = cIn(iX, iRowIn);
			cCpx(iX, 0).Im = 0;
		}

		hr = VtFFTFloat(cCpx.Ptr(), iW, iW, iW, FFFT);
		if(FAILED(hr))
			return hr;

		for(iX=0; iX<iW; iX++)
		{
			float fW = f2PW * (iX<iW2 ? iX : iX - iW);
			cCpx(iX, 0) *= (*pcallback)(fW, pParam);
		}

		hr = VtFFTFloat(cCpx.Ptr(), iW, iW, iW, IFFT);
		if(FAILED(hr))
			return hr;

		for(iX=0; iX<iW; iX++)
			cOut(iX, iRowOut) = fNorm * cCpx(iX, 0).Re;

		return NOERROR;
	}
	else if(cIn.Bands()==1 && cOut.Bands()==2)
	{
		// float to complex (in place)
		CComplexImg& cCpx = (CComplexImg &)cOut;

		int iX;
		for(iX=0; iX<iW; iX++)
		{
			cCpx(iX, iRowOut).Re = cIn(iX, iRowIn);
			cCpx(iX, iRowOut).Im = 0;
		}

		hr = VtFFTFloat(cCpx.Ptr(iRowOut), iW, iW, iW, FFFT);
		if(FAILED(hr))
			return hr;

		for(iX=0; iX<iW; iX++)
		{
			float fW = f2PW * (iX<iW2 ? iX : iX - iW);
			cCpx(iX, iRowOut) *= (*pcallback)(fW, pParam);
		}

		return NOERROR;
	}
	else if(cIn.Bands()==2 && cOut.Bands()==1)
	{
		// complex to float
		CComplexImg cCpx;
		hr = cCpx.Create(iW, 1);
		if(FAILED(hr))
			return hr;

		int iX;
		for(iX=0; iX<iW; iX++)
		{
			cCpx(iX, 0).Re = cIn.Pix(iX, iRowIn, 0);
			cCpx(iX, 0).Im = cIn.Pix(iX, iRowIn, 1);
		}

		for(iX=0; iX<iW; iX++)
		{
			float fW = f2PW * (iX<iW2 ? iX : iX - iW);
			cCpx(iX, 0) *= (*pcallback)(fW, pParam);
		}

		hr = VtFFTFloat(cCpx.Ptr(), iW, iW, iW, IFFT);
		if(FAILED(hr))
			return hr;

		for(iX=0; iX<iW; iX++)
			cOut(iX, iRowOut) = fNorm * cCpx(iX, 0).Re;

		return NOERROR;
	}
	else if(cIn.Bands()==2 && cOut.Bands()==2)
	{
		// complex to complex
		CComplexImg& cCpx = (CComplexImg &)cOut;

		int iX;
		for(iX=0; iX<iW; iX++)
		{
			cCpx(iX, iRowOut).Re = cIn.Pix(iX, iRowIn, 0);
			cCpx(iX, iRowOut).Im = cIn.Pix(iX, iRowIn, 1);
		}

		for(iX=0; iX<iW; iX++)
		{
			float fW = f2PW * (iX<iW2 ? iX : iX - iW);
			cCpx(iX, iRowOut) *= (*pcallback)(fW, pParam);
		}

		return NOERROR;
	}
	else
		return E_INVALIDARG;
}

typedef struct
{
	double fConst;
	double fConst2;
	double fK;
} GaussParams;

Complexf fft1dgaussiancb(float fW, void *p)
{
	GaussParams *ppar = (GaussParams *)p;
	Complexf cx;
	cx.Re = (float)exp(ppar->fConst * fW * fW);
	cx.Im = 0;
	return cx;
}

HRESULT vt::VtFFTFilter1dGaussian(CFloatImg &cOut, int iRowOut, const CFloatImg &cIn, int iRowIn, double fS)
{
	GaussParams par;
	par.fConst = - 0.5 * fS * fS;
	
	return VtFFTFilter1d(cOut, iRowOut, cIn, iRowIn, fft1dgaussiancb, (void *)&par);
}

Complexf fft1ddogcb(float fW, void *p)
{
	GaussParams *ppar = (GaussParams *)p;
	Complexf cx;
	double fWS = (double)fW * (double)fW;
	cx.Re = (float)(exp(ppar->fConst * fWS) - ppar->fK * exp(ppar->fConst2 * fWS));
	cx.Im = 0;
	return cx;
}

HRESULT vt::VtFFTFilter1dDOG(CFloatImg &cOut, int iRowOut, const CFloatImg &cIn, int iRowIn, double fC, double fS, double k)
{
	GaussParams par;
	par.fConst = - 0.5 * fC * fC;
	par.fConst2 = - 0.5 * fS * fS;
	par.fK = k;
	return VtFFTFilter1d(cOut, iRowOut, cIn, iRowIn, fft1ddogcb, (void *)&par);
}

Complexf fft1hilbertcb(float fW, void* /*p*/)
{
	Complexf cx;
	cx.Re = 0;
	cx.Im = fW > 0 ? -1.0f : (fW < 0 ? 1.0f : 0);
	return cx;
}

HRESULT vt::VtFFTFilter1dHilbert(CFloatImg &cOut, int iRowOut, const CFloatImg &cIn, int iRowIn)
{
	return VtFFTFilter1d(cOut, iRowOut, cIn, iRowIn, fft1hilbertcb, NULL);
}

HRESULT vt::VtFFTFilter2d(CFloatImg &cOut, const CFloatImg &cIn, Complexf (*pcallback)(float fWx, float fWy, void *p), void *pParam, bool /*bRetainSrc*/)
{
    HRESULT hr = NOERROR;

	if(!cIn.IsValid() || !cOut.IsValid() ||
		cIn.Width() != cOut.Width() || cIn.Height() != cOut.Height())
		VT_HR_EXIT( E_INVALIDARG );
    if(pcallback==NULL)
        VT_HR_EXIT( E_POINTER );

	int iW = cIn.Width();
    int iH = cIn.Height();

	float f2PW = (float)( VT_PI * 2.0 / iW );
    float f2PH = (float)( VT_PI * 2.0 / iH );
	int iW2 = iW / 2 + 1;
    int iH2 = iH / 2 + 1;
    int iX, iY;

    if(cIn.Bands()==1 && cOut.Bands()==1)
	{
		// float to float
		CComplexImg cCpx;
		VT_HR_EXIT( cCpx.Create(iW, iH, alignAny) );

        VT_HR_EXIT( VtForwardFFT(cCpx, cIn) );

        for(iY=0; iY<iH; iY++)
        {
            float fWy = f2PH * (iY<iH2 ? iY : iY - iH);
		    for(iX=0; iX<iW; iX++)
		    {
			    float fWx = f2PW * (iX<iW2 ? iX : iX - iW);
			    cCpx(iX, iY) *= (*pcallback)(fWx, fWy, pParam);
		    }
        }

		VT_HR_EXIT( VtInverseFFT(cOut, cCpx) );
	}
	else if(cIn.Bands()==1 && cOut.Bands()==2)
	{
		// float to complex (in place)
		CComplexImg& cCpx = (CComplexImg &)cOut;

		VT_HR_EXIT( VtForwardFFT(cCpx, cIn) );

		for(iY=0; iY<iH; iY++)
        {
            float fWy = f2PH * (iY<iH2 ? iY : iY - iH);
		    for(iX=0; iX<iW; iX++)
		    {
			    float fWx = f2PW * (iX<iW2 ? iX : iX - iW);
			    cCpx(iX, iY) *= (*pcallback)(fWx, fWy, pParam);
		    }
        }
	}
	else if(cIn.Bands()==2 && cOut.Bands()==1)
	{
		// complex to float
		CComplexImg cCpx;
		VT_HR_EXIT( cCpx.Create(iW, iH, alignAny) );

		VT_HR_EXIT( cIn.CopyTo(cCpx) ); // to avoid trashing input image

		for(iY=0; iY<iH; iY++)
        {
            float fWy = f2PH * (iY<iH2 ? iY : iY - iH);
		    for(iX=0; iX<iW; iX++)
		    {
			    float fWx = f2PW * (iX<iW2 ? iX : iX - iW);
			    cCpx(iX, iY) *= (*pcallback)(fWx, fWy, pParam);
		    }
        }

        VT_HR_EXIT( VtInverseFFT(cOut, cCpx) );
	}
	else if(cIn.Bands()==2 && cOut.Bands()==2)
	{
		// complex to complex
		CComplexImg& cCpx = (CComplexImg &)cOut;

        VT_HR_EXIT( cIn.CopyTo(cCpx) );

		for(iY=0; iY<iH; iY++)
        {
            float fWy = f2PH * (iY<iH2 ? iY : iY - iH);
		    for(iX=0; iX<iW; iX++)
		    {
			    float fWx = f2PW * (iX<iW2 ? iX : iX - iW);
			    cCpx(iX, iY) *= (*pcallback)(fWx, fWy, pParam);
		    }
        }
	}
	else
		VT_HR_EXIT( E_INVALIDARG );

Exit:
    return hr;
}

// exp(-|X/s|^2) which is exp(-|sW|^2) in freq domain
Complexf fft2dgaussiancb(float fWx, float fWy, void *p)
{
	GaussParams *ppar = (GaussParams *)p;
	Complexf cx;
	cx.Re = (float)exp(ppar->fConst * (fWx * fWx + fWy * fWy));
	cx.Im = 0;
	return cx;
}

HRESULT vt::VtFFTFilter2dGaussian(CFloatImg &cOut, const CFloatImg &cIn, double fS)
{
	GaussParams par;
	par.fConst = - 0.5 * fS * fS;
	
	return VtFFTFilter2d(cOut, cIn, fft2dgaussiancb, (void *)&par);
}

// exp (-|X/s|) which is exp() in freq domain
Complexf fft2dexpcb(float fWx, float fWy, void *p)
{
	GaussParams *ppar = (GaussParams *)p;
	Complexf cx;
	cx.Re = (float)(ppar->fConst * pow(ppar->fConst2 + fWx * fWx + fWy * fWy, -1.5));
	cx.Im = 0;
	return cx;
}

HRESULT vt::VtFFTFilter2dExponential(CFloatImg &cOut, const CFloatImg &cIn, double dE)
{
	GaussParams par;
    par.fConst2 = 1/(dE*dE);
	par.fConst = par.fConst2/dE;    
	
	return VtFFTFilter2d(cOut, cIn, fft2dexpcb, (void *)&par);
}

typedef struct {
	int iOdd;
	float dNorm;
	float dRcs, dRss;
	float dRls;
	float dSin, dCos;
} FFTOrientedDOGParams;

Complexf fft2dorienteddog(float fWx, float fWy, void *pData)
{
	Complexf cx;
	FFTOrientedDOGParams *pPar = (FFTOrientedDOGParams *)pData;

	float dX = pPar->dCos * fWx + pPar->dSin * fWy;
	float dY = -pPar->dSin * fWx + pPar->dCos * fWy;
	float dWs = 0.5f * dX * dX;
	float dM = exp(-0.5f * pPar->dRls * dY * dY)
		* (exp(-pPar->dRcs * dWs) - exp(-pPar->dRss * dWs))
		* pPar->dNorm;

	if(dX>0)
	{
		if(pPar->iOdd)
		{
			cx.Re = 0;
			cx.Im = -dM;
		}
		else
		{
			cx.Re = dM;
			cx.Im = 0;
		}
	}
	else
	{
		if(pPar->iOdd)
		{
			cx.Re = 0;
			cx.Im = dM;
		}
		else
		{
			cx.Re = dM;
			cx.Im = 0;
		}
	}

	return cx;
}

// filter is a Gaussian along the orientation axis and a difference-of-Gaussians orthogonal to it
// returns the result (even-phase) and the same result passed through a hilbert transform (odd-phase)
// fCenter must be less than fSurround, the frequency peak is given unit magnitude
HRESULT vt::VtFFTFilter2dOrientedDOG(CFloatImg &cOutEven, CFloatImg &cOutOdd, const CFloatImg &cIn, 
                                 double dT, double dRc, double dRs, double dRl)
{
    HRESULT hr = NOERROR;
    if(dRs<=dRc)
        VT_HR_EXIT(E_INVALIDARG);

	FFTOrientedDOGParams par;

	double dRcs = dRc * dRc;
	double dRss = dRs * dRs;
	double dWcs = 2 * log(dRcs/dRss)/(dRcs-dRss);
	double dPk = exp(-dRcs * dWcs / 2.0) - exp(-dRss * dWcs / 2.0);
	
	par.dNorm = (float)(1.0 / dPk);
	dT = VtWrapAngle(-dT + 0.5 * VT_PI);
	par.dSin = (float)sin(dT);
	par.dCos = (float)cos(dT);
	par.dRls = (float)(dRl * dRl);
	par.dRcs = (float)dRcs;
	par.dRss = (float)dRss;

    if(cOutEven.IsValid())
    {
	    par.iOdd = 0;
        VT_HR_EXIT( VtFFTFilter2d(cOutEven, cIn, fft2dorienteddog, (void *)&par) );
    }

    if(cOutOdd.IsValid())
    {
	    par.iOdd = 1;
	    VT_HR_EXIT( VtFFTFilter2d(cOutOdd, cIn, fft2dorienteddog, (void *)&par) );
    }

Exit:
    return hr;
}


