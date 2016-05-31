//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Powell's minimization
//
//  History:
//      2006/11/20-swinder
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_matrix.h"
#include "vt_powellmin.h"

using namespace vt;

#define RATIO       1.618034f
#define INVRATIO    0.381966f
#define INITXTOL    2.0e-4f     // minimum tolerance for x in 1d line min
#define EPSIL       1.0e-20f    // avoid divide by zero
#define BIGJUMP     100.0f      // big step for bracket function
#define TOL_LIMIT   1.0e-10f    // prevents tol from getting too small in linemin

#define MAXITERS 100

#define Shift(a, b, c, d) (a)=(b); (b)=(c); (c)=(d)

inline float Sign(float a, float b)
{
    return b>=0 ? fabs(a) : -fabs(a);
}

HRESULT vt::VtLineMinimize1D(float xmin, float xmid, float xmax, float &xrtn, float &fxrtn, 
                         HRESULT (*pFunc)(float x, float &fx, void *p), void *pUser)
{
    if(pFunc==NULL)
        return E_INVALIDARG;

    HRESULT hr = NOERROR;
    float e = 0;
    float a = xmin < xmax ? xmin : xmax;
    float b = xmin > xmax ? xmin : xmax;

    float vsimon[4];
    int iCycle = 0;
    int iCount = 0;

    CVec4f vx, vfx;
    vx[0] = vx[1] = vx[2] = xmid;
    VT_HR_RET( (*pFunc)(vx[2], vfx[2], pUser) );

    vsimon[iCycle] = vfx[2];
    iCycle = (iCycle + 1)%4;
    iCount++;

    vfx[0] = vfx[1] = vfx[2];

    float d=0;

    int i;
    for(i=0; i<MAXITERS; i++)
    {
        float xm = 0.5f * (a + b);
        float tol1 = INITXTOL * fabs(vx[0]) + TOL_LIMIT;
        float tol2 = 2*tol1;
        if(fabs(vx[0] - xm) <= (tol2 - 0.5f * (b - a))
            // prevents the situation where we keep evaluating an function unchanging with x
            || iCount==4 && vsimon[0] == vsimon[1] && vsimon[0] == vsimon[2] && vsimon[0] == vsimon[3])
        {
            xrtn = vx[0];
            fxrtn = vfx[0];
            return hr;
        }
        if(fabs(e) > tol1)
        {
            float r = (vx[0] - vx[1]) * (vfx[0] - vfx[2]);
            float q = (vx[0] - vx[2]) * (vfx[0] - vfx[1]);
            float p = (vx[0] - vx[2]) * q - (vx[0] - vx[1]) * r;
            q = 2*(q-r);
            if(q>0)
                p  = -p;
            q = fabs(q);
            float esave = e;
            e = d;
            if(fabs(p) >= fabs(0.5f * q * esave) || p <= q*(a-vx[0]) || p >= q*(b-vx[0]))
            {
                e = vx[0] >= xm ? a-vx[0] : b-vx[0];
                d = INVRATIO * e;
            }
            else
            {
                d = p/q;
                vx[3] = vx[0] + d;
                if(vx[3]-a < tol2 || b-vx[3] < tol2)
                    d = Sign(tol1, xm - vx[0]);
            }
        }
        else
        {
            e = vx[0] >= xm ? a-vx[0] : b-vx[0];
            d = INVRATIO * e;
        }
        vx[3] = fabs(d) >= tol1 ? vx[0]+d : vx[0] + Sign(tol1, d);
        VT_HR_RET( (*pFunc)(vx[3], vfx[3], pUser) );

        vsimon[iCycle] = vfx[3];
        iCycle = (iCycle + 1)%4;
        iCount++;
        if(iCount>4)
            iCount = 4;

        if(vfx[3] <= vfx[0])
        {
            if(vx[3] >= vx[0])
                a = vx[0];
            else
                b = vx[0];
            Shift(vx[2], vx[1], vx[0], vx[3]);
            Shift(vfx[2], vfx[1], vfx[0], vfx[3]);
        }
        else
        {
            if(vx[3] < vx[0])
                a = vx[3];
            else
                b = vx[3];
            if(vfx[3] <= vfx[1] || vx[1] == vx[0])
            {
                vx[2] = vx[1];
                vx[1] = vx[3];
                vfx[2] = vfx[1];
                vfx[1] = vfx[3];
            }
            else if(vfx[3] <= vfx[2] || vx[2] == vx[0] || vx[2] == vx[1])
            {
                vx[2] = vx[3];
                vfx[2] = vfx[3];
            }
        }
    }

    return E_FAIL;
}

class CPowellMin
{
public:
    CPowellMin() : m_fTol(INITXTOL) {}
    void SetLineMinXTolerance(float fTol) { m_fTol = fTol; }
    HRESULT Minimize(CVecf &vP, const CVecf &vDeltas, float &fRtn,
                             HRESULT (*pFunc)(const CVecf &vParams, float &fValueRtn, void *pUserData),
                             void *pUser, int iMaxIters, float fTol);
private:
    HRESULT LineMinimize(CVecf &vPoint, CVecf &vDirection, float &frtn);
    HRESULT LineMinimize1D(const CVec3f &vx, float &xrtn, float &fxrtn);
    HRESULT Bracket1D(CVec3f &vx, CVec3f &vfx);
    HRESULT Func1D(float fVal, float &fRtn)
    {
        return (*m_pFunc)(m_vPoint + fVal * m_vDirection, fRtn, m_pUserData);
    }

    HRESULT (*m_pFunc)(const CVecf &vParams, float &fValueRtn, void *pUserData);
    void *m_pUserData;
    CVecf m_vPoint, m_vDirection;
    float m_fTol;
};

HRESULT vt::VtPowellSearch(CVecf &vStart, const CVecf &vDeltas, float &fRtn,
                         HRESULT (*pFunc)(const CVecf &vParams, float &fValueRtn, void *pUserData),
                         void *pUser, int iMaxIters, float fFuncTolFrac, float fXTolFrac)
{
    CPowellMin powell;
    powell.SetLineMinXTolerance(fXTolFrac);
    return powell.Minimize(vStart, vDeltas, fRtn, pFunc, pUser, iMaxIters, fFuncTolFrac);
}

HRESULT CPowellMin::Minimize(CVecf &vP, const CVecf &vDeltas, float &fRtn,
                             HRESULT (*pFunc)(const CVecf &vParams, float &fValueRtn, void *pUserData),
                             void *pUser, int iMaxIters, float fTol)
{
    HRESULT hr = NOERROR;

    m_pFunc = pFunc;
    m_pUserData = pUser;
    int n = vP.Size();
    CVecf vPTT(n);
    CVecf vPT = vP;
    CVecf vDir(n);
    CMtxf mDirSet(n, n);
    float fptt;

    mDirSet.MakeDiag(vDeltas);

    VT_HR_RET( (*pFunc)(vP, fRtn, pUser) );

    int iters, i, j;    
    for(iters=0;; iters++)
    {   
        float fp = fRtn;
        int iBig = 0;
        float delta = 0;
        for(i=0; i<n; i++)
        {
            for(j=0; j<n; j++)
                vDir[j] = mDirSet[j][i]; // extract column
            fptt = fRtn;
            VT_HR_RET( LineMinimize(vP, vDir, fRtn) );
            if(fabs(fptt-fRtn) > delta)
            {
                delta = fabs(fptt-fRtn);
                iBig = i;
            }
        }

        if(2*fabs(fp - fRtn) < fTol * (fabs(fp)+fabs(fRtn)))
            break;
        if(iters>=iMaxIters)
            return E_NOTFOUND;

        vPTT = 2.0f * vP - vPT;
        vDir = vP - vPT;
        vPT = vP;
        
        VT_HR_RET( (*pFunc)(vPTT, fptt, pUser) );

        if(fptt < fp)
        {
            float t1 = fp - fRtn - delta;
            float t2 = fp - fptt;
            float t = 2*(fp - 2*fRtn + fptt) * t1 * t1 - delta * t2 * t2;
            if(t<0)
            {
                VT_HR_RET( LineMinimize(vP, vDir, fRtn) );
                for(j=0; j<n; j++)
                {
                    mDirSet[j][iBig] = mDirSet[j][n-1];
                    mDirSet[j][n-1] = vDir[j];
                }
            }
        }
    }

    return hr;
}

HRESULT CPowellMin::LineMinimize(CVecf &vPoint, CVecf &vDirection, float &frtn)
{
    HRESULT hr = NOERROR;

    m_vPoint = vPoint;
    m_vDirection = vDirection;

    CVec3f vBracket, vFuncEval;
    vBracket[0] = 0;
    vBracket[1] = 1;
    VT_HR_RET( Bracket1D(vBracket, vFuncEval) );
//#undef printf
    //printf("bracket %g %g %g (%g %g %g)\n", vBracket[0], vBracket[1], vBracket[2], vFuncEval[0], vFuncEval[1], vFuncEval[2]);

    float xmin;
    VT_HR_RET( LineMinimize1D(vBracket, xmin, frtn) );
    //printf("line min %g (%g)\n", xmin, frtn);

    vDirection *= xmin;
    vPoint += vDirection;

    return hr;
}

HRESULT CPowellMin::Bracket1D(CVec3f &vx, CVec3f &vfx)
{
    HRESULT hr = NOERROR;
    float tmp = 0;

    VT_HR_RET( Func1D(vx[0], vfx[0]) );
    VT_HR_RET( Func1D(vx[1], vfx[1]) );
    if(vfx[1] > vfx[0])
    {
        // swap 0,1
        Shift(tmp, vx[0], vx[1], tmp);
        Shift(tmp, vfx[0], vfx[1], tmp);
    }
    // fx[1] is minimum
    vx[2] = vx[1] + RATIO * (vx[1] - vx[0]);
    VT_HR_RET( Func1D(vx[2], vfx[2]) );
    while(vfx[1] > vfx[2])
    {
        float a = (vx[1] - vx[0]) * (vfx[1] - vfx[2]);
        float b = (vx[1] - vx[2]) * (vfx[1] - vfx[0]);
        tmp = fabs(b - a);
        tmp = VtMax(tmp, EPSIL);
        float x = vx[1] - ((vx[1] - vx[2]) * b - (vx[1] - vx[0]) * a)
            / (2 * Sign(tmp, b - a));
        float xupper = vx[1] + BIGJUMP * (vx[2] - vx[1]);
        float fx;
        if((vx[1] - x) * (x - vx[2]) > 0)
        {
            VT_HR_RET( Func1D(x, fx) );
            if(fx < vfx[2])
            {
                vx[0] = vx[1];
                vx[1] = x;
                vfx[0] = vfx[1];
                vfx[1] = fx;
                return hr;
            }
            else if(fx > vfx[1])
            {
                vx[2] = x;
                vfx[2] = fx;
                return hr;
            }
            x = vx[2] + RATIO * (vx[2] - vx[1]);
            VT_HR_RET( Func1D(x, fx) );
        }
        else if((vx[2] - x) * (x - xupper) > 0)
        {
            VT_HR_RET( Func1D(x, fx) );
            if(fx < vfx[2])
            {
                vx[1] = vx[2];
                vx[2] = x;
                x = vx[2] + RATIO * (vx[2] - vx[1]);
                vfx[1] = vfx[2];
                vfx[2] = fx;
                VT_HR_RET( Func1D(x, fx) ); 
            }
        }
        else if((x - xupper) * (xupper - vx[2]) >= 0)
        {
            x = xupper;
            VT_HR_RET( Func1D(x, fx) );
        }
        else
        {
            x = vx[2] + RATIO * (vx[2] - vx[1]);
            VT_HR_RET( Func1D(x, fx) );
        }
        Shift(vx[0], vx[1], vx[2], x);
        Shift(vfx[0], vfx[1], vfx[2], fx);
    }

    return hr;
}


HRESULT CPowellMin::LineMinimize1D(const CVec3f &vu, float &xrtn, float &fxrtn)
{
    HRESULT hr = NOERROR;
    float e = 0;
    float a = vu[0] < vu[2] ? vu[0] : vu[2];
    float b = vu[0] > vu[2] ? vu[0] : vu[2];

    float vsimon[4];
    int iCycle = 0;
    int iCount = 0;

    CVec4f vx, vfx;
    vx[0] = vx[1] = vx[2] = vu[1];
    VT_HR_RET( Func1D(vx[2], vfx[2]) );
    
    vsimon[iCycle] = vfx[2];
    iCycle = (iCycle + 1)%4;
    iCount++;

    vfx[0] = vfx[1] = vfx[2];

    float d=0;

    int i;
    for(i=0; i<MAXITERS; i++)
    {
        float xm = 0.5f * (a + b);
        float tol1 = m_fTol * fabs(vx[0]) + TOL_LIMIT;
        float tol2 = 2*tol1;
        if(fabs(vx[0] - xm) <= (tol2 - 0.5f * (b - a))
            // prevents the situation where we keep evaluating an function unchanging with x
            || iCount==4 && vsimon[0] == vsimon[1] && vsimon[0] == vsimon[2] && vsimon[0] == vsimon[3])
        {
            xrtn = vx[0];
            fxrtn = vfx[0];
            return hr;
        }
        if(fabs(e) > tol1)
        {
            float r = (vx[0] - vx[1]) * (vfx[0] - vfx[2]);
            float q = (vx[0] - vx[2]) * (vfx[0] - vfx[1]);
            float p = (vx[0] - vx[2]) * q - (vx[0] - vx[1]) * r;
            q = 2*(q-r);
            if(q>0)
                p  = -p;
            q = fabs(q);
            float esave = e;
            e = d;
            if(fabs(p) >= fabs(0.5f * q * esave) || p <= q*(a-vx[0]) || p >= q*(b-vx[0]))
            {
                e = vx[0] >= xm ? a-vx[0] : b-vx[0];
                d = INVRATIO * e;
            }
            else
            {
                d = p/q;
                vx[3] = vx[0] + d;
                if(vx[3]-a < tol2 || b-vx[3] < tol2)
                    d = Sign(tol1, xm - vx[0]);
            }
        }
        else
        {
            e = vx[0] >= xm ? a-vx[0] : b-vx[0];
            d = INVRATIO * e;
        }
        vx[3] = fabs(d) >= tol1 ? vx[0]+d : vx[0] + Sign(tol1, d);
        VT_HR_RET( Func1D(vx[3], vfx[3]) );

        vsimon[iCycle] = vfx[3];
        iCycle = (iCycle + 1)%4;
        iCount++;
        if(iCount>4)
            iCount = 4;

        if(vfx[3] <= vfx[0])
        {
            if(vx[3] >= vx[0])
                a = vx[0];
            else
                b = vx[0];
            Shift(vx[2], vx[1], vx[0], vx[3]);
            Shift(vfx[2], vfx[1], vfx[0], vfx[3]);
        }
        else
        {
            if(vx[3] < vx[0])
                a = vx[3];
            else
                b = vx[3];
            if(vfx[3] <= vfx[1] || vx[1] == vx[0])
            {
                vx[2] = vx[1];
                vx[1] = vx[3];
                vfx[2] = vfx[1];
                vfx[1] = vfx[3];
            }
            else if(vfx[3] <= vfx[2] || vx[2] == vx[0] || vx[2] == vx[1])
            {
                vx[2] = vx[3];
                vfx[2] = vfx[3];
            }
        }
    }

    return E_FAIL;
}
