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

#include "stdafx.h"

#include "vt_mathutils.h"
#include "vt_haar.h"

using namespace vt;

// Haar transform of 1D vector
// data and tmp both have size elements
void CHaar2D::Haar1D(float *data, float *tmp, int size)
{
    while (size > 1)
    {
        size /= 2;
        int i;
        for (i = 0; i < size; i++)
        {
            tmp[i]        = (data[2 * i] + data[2 * i + 1]) * (float)VT_SQRT1_2;
            tmp[size + i] = (data[2 * i] - data[2 * i + 1]) * (float)VT_SQRT1_2;
        }
        memcpy(data, tmp, 2 * size * sizeof(float));
	}
}

// 1D Haar transform of 2D matrix
void CHaar2D::Haar1(float *data, int size, float *tmp)
{
    float *rowptr = data;
    int i;
    for (i = 0; i < size; i++)
    {
        Haar1D(rowptr, tmp, size);
        rowptr += size;
    }
}

// transpose data
void CHaar2D::Transpose(float *data, int size, float *tmp)
{
    int i,j;
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            tmp[i + j * size] = data[j + i * size];

    memcpy(data, tmp, size * size * sizeof(float));
}

// 2D Haar transform of 2D matrix
HRESULT CHaar2D::Process(float *data, int size)
{
    VT_HR_BEGIN();

    int size2 = size*size;
    if(m_vHaar.size()!=size_t(size2))
        VT_HR_EXIT( m_vHaar.resize(size2) );

    Haar1(data, size, &(m_vHaar[0]));
    Transpose(data, size, &(m_vHaar[0]));
    Haar1(data, size, &(m_vHaar[0]));
    Transpose(data, size, &(m_vHaar[0]));

    VT_HR_END();
}
