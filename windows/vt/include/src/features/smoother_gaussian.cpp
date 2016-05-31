//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description: 2 implementations of IFeatureWarpCompute that use a simple
//               gaussian smooth on the transform parameters to compute
//               a smooth camera path
//
//------------------------------------------------------------------------------
#include "stdafx.h"

using namespace vt;

//+-----------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
bool ComputeHomography4Points(CMtx3x3f& mtx, 
                              const vt::CVec2f* p1, const vt::CVec2f* p2)
{
    double Astore[64];
    CMtxd A;
    A.Wrap(Astore, 8, 8);
    ZeroMemory(A.Ptr(), A.Size()*sizeof(double));

    double RHSstore[8], LHSstore[8];
    CVecd RHS, LHS;
    RHS.Wrap(RHSstore, 8);
    LHS.Wrap(LHSstore, 8);

    // create the system of equations to solve
    for( int i = 0; i < 4; i++, p1++, p2++ )
    {
        A(i,0)=p1->x, A(i,1)=p1->y, A(i,2)=1, A(i,6)=-p1->x*p2->x, A(i,7)=-p1->y*p2->x;
        RHS(i) = p2->x;

        int j = i+4;
        A(j,3)=p1->x, A(j,4)=p1->y, A(j,5)=1, A(j,6)=-p1->x*p2->y, A(j,7)=-p1->y*p2->y;
        RHS(j) = p2->y;
    }

    // solve
    CMtxd AInv = A.Inv();
    if( AInv.IsError() )
	{
		return false;
	}

	LHS = AInv*RHS;
	mtx = CMtx3x3f(float(LHS[0]), float(LHS[1]), float(LHS[2]),
				   float(LHS[3]), float(LHS[4]), float(LHS[5]),
				   float(LHS[6]), float(LHS[7]), 1.f);
    return true;
}

struct STAB_SIMILARITY
{
    float s;
    float r;
    float tx;
    float ty;

};

HRESULT SmoothTransformsSimilarity(vt::vector<CMtx3x3f>& vecResult)
{
    VT_HR_BEGIN()

    vt::vector<STAB_SIMILARITY> vecSimilarity, vecSmoothedSimilarity;

    VT_HR_EXIT( vecSimilarity.resize(vecResult.size()) );
    VT_HR_EXIT( vecSmoothedSimilarity.resize(vecResult.size()) );

    // extract the similarity parameters
    for( UINT i = 0; i < vecResult.size(); i++ )
    {
        const CMtx3x3f& x    = vecResult[i];
        STAB_SIMILARITY& sim = vecSimilarity[i]; 
        sim.s = sqrtf(x(0,0)*x(0,0)+x(0,1)*x(0,1));
        sim.r = acosf(x(0,0)/sim.s);
        sim.tx = x(0,2);
        sim.ty = x(1,2);
    }

    C1dKernel k;
    Create1dGaussKernel(k, 1.0);
    VtFilter1d(vecSmoothedSimilarity.begin(), EL_FORMAT_FLOAT, 4,
               vecSimilarity.begin(), EL_FORMAT_FLOAT, 4, 
               int(vecSmoothedSimilarity.size()), k );

    for( UINT i = 0; i < vecResult.size(); i++ )
    {
        const STAB_SIMILARITY& sim    = vecSimilarity[i]; 
        const STAB_SIMILARITY& smooth = vecSmoothedSimilarity[i]; 

        printf("%u, %f, %f, %f, %f, %f, %f, %f, %f\n",
               i, sim.s, sim.r*180.f/3.14159f, sim.tx, sim.ty,
               smooth.s, smooth.r*180.f/3.14159f, smooth.tx, smooth.ty);

        float cos_z = smooth.s*cos(smooth.r);
        float sin_z = smooth.s*sin(smooth.r);
        vecResult[i] = CMtx3x3f(cos_z, -sin_z, smooth.tx,
                                sin_z, cos_z,  smooth.ty,
                                    0,     0,  1);
    }

    VT_HR_END()
}


HRESULT SmoothTransformsGeneral(vt::vector<CMtx3x3f>& vecResult, 
                                int iFrameWidth, int iFrameHeight)
{
    VT_HR_BEGIN()

	const int iResultSize = int(vecResult.size());
    vt::vector<CVec2f> vecCorners;
    VT_HR_EXIT( vecCorners.resize(5*iResultSize) );

    // for each transform project the corners into points
    CVec2f *pTL = vecCorners.begin();
    CVec2f *pTR = vecCorners.begin()+1*iResultSize;
    CVec2f *pBR = vecCorners.begin()+2*iResultSize;
    CVec2f *pBL = vecCorners.begin()+3*iResultSize;

    for( int i = 0; i < iResultSize; i++, pTL++, pTR++, pBR++, pBL++ )
    {
        const CMtx3x3f& x = vecResult[i];
        CVec3f tl = x*CVec3f(0,0,1.f);
        CVec3f tr = x*CVec3f(float(iFrameWidth),0,1.f);
        CVec3f br = x*CVec3f(float(iFrameWidth),float(iFrameHeight),1.f);
        CVec3f bl = x*CVec3f(0,float(iFrameHeight),1.f);

        *pTL = tl.Dehom();
        *pTR = tr.Dehom();
        *pBR = br.Dehom();
        *pBL = bl.Dehom();
    }

    // smooth the projected corners
    pTL = vecCorners.begin();
    pTR = vecCorners.begin()+1*iResultSize;
    pBR = vecCorners.begin()+2*iResultSize;
    pBL = vecCorners.begin()+3*iResultSize;
    CVec2f *pTmp = vecCorners.begin()+4*iResultSize;

    //C14641Kernel k;
    C1dKernel k;
    Create1dGaussKernel(k, 2.0);
    //float testK[3] = {0.1f, 0.8f, 0.1f};
    //k.Create(3,1,testK);

    VtFilter1d(pTmp, EL_FORMAT_FLOAT, 2, pTL, EL_FORMAT_FLOAT, 2, iResultSize, k );  //tl
    VtFilter1d(pTL,  EL_FORMAT_FLOAT, 2, pTR, EL_FORMAT_FLOAT, 2, iResultSize, k );  //tr
    VtFilter1d(pTR,  EL_FORMAT_FLOAT, 2, pBR, EL_FORMAT_FLOAT, 2, iResultSize, k );  //br
    VtFilter1d(pBR,  EL_FORMAT_FLOAT, 2, pBL, EL_FORMAT_FLOAT, 2, iResultSize, k );  //bl

    // aliasing of names is a bit confusing here - used rolling buffer to filter into
    pBL = pBR;
    pBR = pTR;
    pTR = pTL;
    pTL = pTmp;

    // for each set of source of smoothed points compute a new transform
    CVec2f pts1[4], pts2[4];
    pts1[0] = CVec2f(0,0);
    pts1[1] = CVec2f(float(iFrameWidth),0);
    pts1[2] = CVec2f(float(iFrameWidth),float(iFrameHeight));
    pts1[3] = CVec2f(0,float(iFrameHeight));

    for( int i = 0; i < iResultSize; i++, pTL++, pTR++, pBR++, pBL++ )
    {
        pts2[0] = *pTL;
        pts2[1] = *pTR;
        pts2[2] = *pBR;
        pts2[3] = *pBL;
        ComputeHomography4Points(vecResult[i], pts1, pts2);
    }

    VT_HR_END()
}

