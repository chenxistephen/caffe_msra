//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Class for generating uniform and Gaussian deviates
//
//  History:
//      2003/11/12-swinder
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_rand.h"
#include "vt_mathutils.h"
#include "vt_matrix.h"
#include "vt_mathutils.inl"

namespace vt {

#define RAND_Q1 16807
#define RAND_Q2 2147483647
#define RAND_Q3 (1.0/RAND_Q2)
#define RAND_Q4 127773
#define RAND_Q5 2836
#define RAND_NORM (1 + (RAND_Q2-1)/RAND_NTAB)
#define RAND_TINY 1.2e-7
#define RAND_RMAX (1.0 - RAND_TINY)


////////////////////////////////////////////////////////////////////////
// CRand - Good Random Number Generator
////////////////////////////////////////////////////////////////////////
CRand::CRand(int iSeed)
{
    Seed(iSeed);
}

// Note: /analyze warning here is a false positive so disabling it.  When
//       a new version of the tools comes out we should check if this false
//       positive is fixed.
#pragma warning( push )
#pragma warning ( disable : 6201 )

void CRand::Seed(int iSeed)
{
	if(iSeed==0)
        iSeed = 1;
    else if(iSeed<0)
        iSeed = -iSeed;

    int i;
    m_iState = iSeed;

    for(i = RAND_NTAB+7; i>=0; i--)
    {
        int iK = m_iState / RAND_Q4;
        m_iState = RAND_Q1 * (m_iState - iK * RAND_Q4) - RAND_Q5 * iK;
        if(m_iState<0)
            m_iState += RAND_Q2;
        if( i < RAND_NTAB)
        {
            m_rgiShuffle[i] = m_iState;
        }
    }

    m_iLast = m_rgiShuffle[0];
    m_bGaussISet = 0;
}

#pragma warning( pop )

double CRand::DRand()
{
    int iK = m_iState / RAND_Q4;
    m_iState = RAND_Q1 * (m_iState - iK * RAND_Q4) - RAND_Q5 * iK;
    if(m_iState<0)
        m_iState += RAND_Q2;
    
    int i = m_iLast / RAND_NORM;
    m_iLast = m_rgiShuffle[i];
    m_rgiShuffle[i] = m_iState;

    double dResult = RAND_Q3 * m_iLast;
    if(dResult>RAND_RMAX)
        dResult = RAND_RMAX;

    return dResult;
}

double CRand::Gauss()
{
    if(!m_bGaussISet)
    {
        double dV1, dV2, dRsq;
        do {
            dV1 = 2.0 * DRand() - 1.0;
            dV2 = 2.0 * DRand() - 1.0;
            dRsq = dV1 * dV1 + dV2 * dV2;
        } while(dRsq>=1.0 || dRsq==0.0);
        
        double dFac = sqrt(-2.0 * log(dRsq) / dRsq);
        m_dGaussGSet = dV1 * dFac;
        m_bGaussISet = 1;

        return dV2 * dFac;
    }
    else
    {
        m_bGaussISet = 0;
        
        return m_dGaussGSet;
    }
}


////////////////////////////////////////////////////////////////////////
// CRC4 - RC4 encryption/decryption class
////////////////////////////////////////////////////////////////////////

HRESULT CRC4::Init(Byte *pchKey, int iKeyLen)
{
    int rgiKey[256];
    int i, j;
    if(iKeyLen<4)
        return E_INVALIDARG;
    if(pchKey==NULL)
        return E_POINTER;

    for(i=0; i<256; i++)
    {
        m_rgiS[i] = i;
        rgiKey[i] = pchKey[i%iKeyLen];
    }
    j = 0;
    for(i=0; i<256; i++)
    {
        j = (j + m_rgiS[i] + rgiKey[i]) & 255;
        int iSTmp = m_rgiS[i];
        m_rgiS[i] = m_rgiS[j];
        m_rgiS[j] = iSTmp;
    }

    m_j = 0;
    m_i = 0;
    VtMemset(rgiKey, 0, 256*sizeof(int));

    return NOERROR;
}

HRESULT CRC4::Process(Byte *pchMsg, int iMsgLen)
{
    if(pchMsg==NULL)
        return E_POINTER;
    for(; iMsgLen>0; pchMsg++, iMsgLen--)
    {
        m_i = (m_i + 1) & 255;
        m_j = (m_j + m_rgiS[m_i]) & 255;
        int iSTmp = m_rgiS[m_i];
        m_rgiS[m_i] = m_rgiS[m_j];
        m_rgiS[m_j] = iSTmp;
        *pchMsg ^= m_rgiS[(m_rgiS[m_i] + m_rgiS[m_j]) & 255];
    }

    return NOERROR;
}

void CRC4::Skip(int iLen)
{
    for(; iLen>0; iLen--)
    {
        m_i = (m_i + 1) & 255;
        m_j = (m_j + m_rgiS[m_i]) & 255;
        int iSTmp = m_rgiS[m_i];
        m_rgiS[m_i] = m_rgiS[m_j];
        m_rgiS[m_j] = iSTmp;
    }
}

CRC4::~CRC4()
{
    VtMemset(m_rgiS, 0, 256*sizeof(int));
}

};