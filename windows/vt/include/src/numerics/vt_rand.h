//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Class for generating uniform and Gaussian deviates
//
//  History:
//      2003/11/12-swinder
//			Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"

namespace vt {

#define RAND_NTAB 32

////////////////////////////////////////////////////////////////////////
// CRand - Good Random Number Generator
////////////////////////////////////////////////////////////////////////

class CRand
{
public:
	CRand(int iSeed=1);

	void Seed(int iSeed);
	double DRand();							// [0,1)
	double URand(double fmin, double fmax)	// [fmin,fmax)
		{ return fmin + (fmax-fmin) * DRand(); }
	int IRand(int iMod)						// 0 to iMod - 1
		{ return (int)(iMod * DRand()); }
	double Gauss();

private:
	int m_iLast;
	int m_iState;
	int m_rgiShuffle[RAND_NTAB];
	bool m_bGaussISet;
	double m_dGaussGSet;
};


////////////////////////////////////////////////////////////////////////
// CRC4 - RC4 encryption/decryption class
////////////////////////////////////////////////////////////////////////

class CRC4
{
public:
	~CRC4();
	HRESULT Init(Byte *pchKey, int iKeyLen);
	HRESULT Process(Byte *pchMsg, int iMsgLen);
	void Skip(int iLen);

private:
	int m_i, m_j;
	int m_rgiS[256];
};

};
