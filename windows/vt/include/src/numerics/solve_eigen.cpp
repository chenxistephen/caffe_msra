//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      These routines were developed by geoffrey cross MSR Cambridge and
//        adapted for vision tools by swinder.
//
//  History:
//      2004/11/17-swinder
//			Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_matrix.h"
#include "vt_solve_schur.h"
#include "vt_solve_eigen.h"

using namespace vt;

/// stand alone eigenroutines

HRESULT vt::EigHouseholderReduction(CMtxd &m, CVecd &vDiag, CVecd &vOffDiag)
{
    int n = m.Rows();
    HRESULT hr = NOERROR;
	VT_HR_EXIT(vDiag.Create(n));
    VT_HR_EXIT(vOffDiag.Create(n));

    int i, j, k;
    for(i = n-1; i>0; i--)
    {
        double sig = 0.0;
        double sum = 0.0;
        if(i>1)
        {
            for(k=0; k<i; k++)
                sum += fabs(m(i,k));
            if(sum == 0.0)
                vOffDiag[i] = m(i,i-1);
            else
            {
                for(k=0; k<i; k++)
                {
                    m(i,k) /= sum;
                    sig += m(i,k) * m(i,k);
                }
                double f = m(i,i-1);
                double g = f>=0.0 ? -sqrt(sig) : sqrt(sig);
                vOffDiag[i] = sum * g;
                sig -= f * g;
                m(i,i-1) = f - g;
                f = 0.0;
                for(j=0; j<i; j++)
                {
                    m(j,i) = m(i,j) / sig;
                    g = 0.0;
                    for(k=0; k<=j; k++)
                        g += m(j,k) * m(i,k);
                    for(k=j+1; k<i; k++)
                        g += m(k,j) * m(i,k);
                    vOffDiag[j] = g/sig;
                    f += vOffDiag[j] * m(i,j);
                }
                double fk = 0.5 * f / sig;
                for(j=0; j<i; j++)
                {
                    f = m(i,j);
                    g = vOffDiag[j] - fk * f;
                    vOffDiag[j] = g;
                    for(k=0; k<=j; k++)
                        m(j,k) -= f * vOffDiag[k] + g * m(i,k);
                }
            }
        }
        else
            vOffDiag[i] = m(i,i-1);
        vDiag[i] = sig;
    }
    
    vDiag[0] = 0.0;
    vOffDiag[0] = 0.0;

    for(i=0; i<n; i++)
    {
        if(vDiag[i] != 0.0)
        {
            for(j=0; j<i; j++)
            {
                double g = 0.0;
                for(k=0; k<i; k++)
                    g += m(i,k) * m(k,j);
                for(k=0; k<i; k++)
                    m(k,j) -= g * m(k,i);
            }
        }
        vDiag[i] = m(i,i);
        m(i,i) = 1.0;
        for(j=0; j<i; j++)
        {
            m(j,i) = 0.0;
            m(i,j) = 0.0;
        }
    }
Exit:
	return hr;
}

HRESULT vt::EigTridiagonalQLImplicit(CVecd &vDiag, CVecd &vOffDiag, CMtxd &mRot)
{
    int n = vDiag.Size();
    int i, k, j, d;
    for(i=1; i<n; i++)
        vOffDiag[i-1] = vOffDiag[i];
    vOffDiag[n-1] = 0.0;

    for(d=0; d<n; d++)
    {
        int iter = 0;
        do {
            for(j=d; j<n-1; j++)
            {
                double tmp = fabs(vDiag[j]) + fabs(vDiag[j+1]);
                if(fabs(vOffDiag[j]) + tmp == tmp)
                    break;
            }
            if(j != d)
            {
                iter++;
                if(iter>1000)
                    return E_FAIL;
                double u = (vDiag[d+1] - vDiag[d]) / (2*vOffDiag[d]);
                double r = VtHypot(u, 1.0);
                u = vDiag[j] - vDiag[d] + vOffDiag[d] / (u + (u>=0 ? fabs(r) : -fabs(r)));
                double s = 1.0;
                double c = 1.0;
                double q = 0.0;
                for(i=j-1; i>=d; i--)
                {
                    double a = s * vOffDiag[i];
                    double b = c * vOffDiag[i];
                    r = VtHypot(a, u);
                    vOffDiag[i+1] = r;
                    if(r==0.0)
                    {
                        vDiag[i+1] -= q;
                        vOffDiag[j] = 0.0;
                        break;
                    }
                    s = a/r;
                    c = u/r;
                    u = vDiag[i+1] - q;
                    r = (vDiag[i] - u) * s + 2 * c * b;
                    q = s * r;
                    vDiag[i+1] = u + q;
                    u = c * r - b;
                    for(k=0; k<n; k++)
                    {
                        a = mRot(k,i+1);
                        mRot(k,i+1) = s * mRot(k,i) + c * a;
                        mRot(k,i) = c * mRot(k,i) - s * a;
                    }
                }
                if(r==0.0 && i >= d)
                    continue;
                vDiag[d] -= q;
                vOffDiag[d] = u;
                vOffDiag[j] = 0.0;
            }
        } while(j != d);
    }

    return NOERROR;
}

HRESULT vt::VtEigenDecomposition(const CMtxd &mA, CMtxd &mV, CVecd &vD)
{
	HRESULT hr = NOERROR;

    int n = mA.Rows();
    CVecd vOffD;
    mV = mA;
	if(mV.IsError())
		VT_HR_EXIT(E_OUTOFMEMORY);

    VT_HR_EXIT(EigHouseholderReduction(mV, vD, vOffD));

    VT_HR_EXIT(EigTridiagonalQLImplicit(vD, vOffD, mV));

	// sort eigenvectors into descending eigenvalue order
    int i, j;
    for(i=0; i<n-1; i++)
    {
        int big = i;
        double fbig = fabs(vD[i]);
        for(j=i+1; j<n; j++)
        {
            double f = fabs(vD[j]);
            if(f>fbig)
            {
                fbig = f;
                big = j;
            }
        }
        double tmp = vD[i];
        vD[i] = vD[big];
        vD[big] = tmp;
        for(j=0; j<n; j++)
        {
            tmp = mV(j, i);
            mV(j, i) = mV(j, big);
            mV(j, big) = tmp;
        }
    }

Exit:
    return hr;
}
