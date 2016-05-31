//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Numerical matrix routines
//
//  History:
//      2004/11/08-swinder
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_matrix.h"

using namespace vt;

namespace vt {

template<>
Complex<float> CMtx<Complex<float> >::Det() const
{
    if(IsError() || Rows()!=Cols() || Rows()==0)
        return (Complex<float>)0.0f;
    if(Rows()==1)
        return El(0,0);
    if(Rows()==2)
        return El(0,0) * El(1,1) - El(0,1) * El(1,0);

    // simple but uselessly slow recursive way for any size >6
    // we need to have a complex decomposition class to solve this
    // quickly
    Complex<float> f = 0;
    int i, s = 1;
    for(i=0; i<Cols(); i++, s = -s)
        f += ((float)s) * El(0, i) * DeleteRowCol(0, i).Det();

    return f;
}

template<>
Complex<double> CMtx<Complex<double> >::Det() const
{
    if(IsError() || Rows()!=Cols() || Rows()==0)
        return (Complex<double>)0;
    if(Rows()==1)
        return El(0,0);
    if(Rows()==2)
        return El(0,0) * El(1,1) - El(0,1) * El(1,0);

    // simple but uselessly slow recursive way for any size >6
    // we need to have a complex decomposition class to solve this
    // quickly
    Complex<double> f = 0;
    int i, s = 1;
    for(i=0; i<Cols(); i++, s = -s)
        f += ((double)s) * El(0, i) * DeleteRowCol(0, i).Det();

    return f;
}

}

// multiplies a matrix with itself to form a symmetric matrix
// mdst = msrc.T() * msrc
HRESULT vt::VtFastMatrixMul(CMtx<float> &msrc, CMtx<float> &mdst)
{
	HRESULT hr = NOERROR;

	hr = mdst.Create(msrc.Cols(), msrc.Cols());
	if(FAILED(hr))
		return hr;

	int i,j,k;
			
	if(g_SupportSSE2())
	{
#if (defined(_M_IX86) || defined(_M_AMD64))
		for(i=0; i<msrc.Cols(); i++)
		{
			for(j=i; j<msrc.Cols()-3; j+=4)
			{
				// sum over k
				// msrc(k,i) * [msrc(k,j)msrc(k,j+1)msrc(k,j+2)msrc(k,j+3)]
				// save in [mdst(i,j)mdst(i,j+1)mdst(i,j+2)mdst(i,j+3)]

				__m128 mmsum = _mm_setzero_ps();
				float *prow = msrc.Ptr();
				for(k=0; k<msrc.Rows(); k++, prow += msrc.Cols())
				{
					__m128 mm1 = _mm_set1_ps(prow[i]);
					__m128 mm2 = _mm_loadu_ps(prow + j);
					__m128 mm3 = _mm_mul_ps(mm1, mm2);
					mmsum = _mm_add_ps(mm3, mmsum);
				}
				_mm_storeu_ps(&(mdst(i,j)), mmsum);
			}
			for(; j<msrc.Cols(); j++)
			{
				float fsum = 0;
				for(k=0; k<msrc.Rows(); k++)
					fsum += msrc(k,i) * msrc(k,j);
				mdst(i,j) = fsum;
			}
		}
#endif
	}
	else
	{
		for(i=0; i<msrc.Cols(); i++)
			for(j=i; j<msrc.Cols(); j++)
			{
				float fsum = 0;
				for(k=0; k<msrc.Rows(); k++)
					fsum += msrc(k,i) * msrc(k,j);
				mdst(i,j) = fsum;
			}
	}

	for(i=1; i<msrc.Cols(); i++)
		for(j=0; j<i; j++)
			mdst(i,j) = mdst(j,i);

	return hr;
}
