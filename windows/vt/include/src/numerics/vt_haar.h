//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Class to do 2d haar wavelet decomp
//
//  History:
//      2007/10/4-swinder
//			Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"

namespace vt {

class CHaar2D
{
public:
	CHaar2D() {}
	
	// The wavelet is computed on a size x size 2D array.
	// Size *must* be a power of two.
	// Compute the wavelet in place.
	HRESULT Process(float *p, int iSize);

private:
	void Haar1D(float *data, float *tmp, int size);
	void Haar1(float *data, int size, float *tmp);
	void Transpose(float *data, int size, float *tmp);

	vector<float> m_vHaar;
};

};
