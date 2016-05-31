//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Model fitting from point matches
//
//  History:
//      2005/5/18-swinder
//          Created
//      2012/12/6-luyuan
//          Modified, Add Homography Model Fitting
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_modelfitting.h"
#include "vt_matrix.h"
#include "vt_solve_svd.h"

using namespace vt;

#define MIN_SQ_LENGTH      1e-6     // square of minimum distance between two source points for computing similarity

// Calculate the appropriate transform from the list of point matches provided.
// The minimum number of matches is 2 for Similarity, 3 for Affine and 4 for Homography.
// If the number of matches provided exceeds the minimum requirement then the solution
// returned is correct in the least squares sense. The minimization assumes a circular
// Gaussian error on each point p1 and exact measurements for p0. The error equation is
// given by E = sum_i ( ||p1_i - transformed(p0_i)||^2 )
HRESULT vt::VtSimilarityFromPointMatches2D(OUT CMtx3x3d &mS, 
                                           IN  const PointMatch* pMatch,
                                           unsigned int uMatchCount)
{
    HRESULT hr = NOERROR;
    if(uMatchCount<2)
        VT_HR_EXIT( E_INVALIDARG );
    if(uMatchCount==2)
    {
        // special case for exact number of matches
        CVec2f vdp0 = pMatch[1].p0 - pMatch[0].p0;
        double fDen = vdp0.x * (double)vdp0.x + vdp0.y * (double)vdp0.y;

        if(fDen<MIN_SQ_LENGTH)
        {
            mS = 0; // case when two initial points very closely coincide
            mS(2,2) = 1;
            goto Exit;
        }
        
        CVec2f vdp1 = pMatch[1].p1 - pMatch[0].p1;
        mS(0,0) = (vdp1.x * (double)vdp0.x + vdp1.y * (double)vdp0.y) / fDen;
        mS(0,1) = (vdp1.x * (double)vdp0.y - vdp1.y * (double)vdp0.x) / fDen;
        mS(1,0) = -mS(0,1);
        mS(1,1) = mS(0,0);
        mS(0,2) = pMatch[0].p1.x - (mS(0,0) * pMatch[0].p0.x + mS(0,1) * pMatch[0].p0.y);
        mS(1,2) = pMatch[0].p1.y - (mS(1,0) * pMatch[0].p0.x + mS(1,1) * pMatch[0].p0.y);
        mS(2,0) = 0;
        mS(2,1) = 0;
        mS(2,2) = 1;
    }
    else
    {
        // least squares estimate of similarity
        // matrix parameterised by:
        // (a b c; -b a d; 0 0 1)
        CVec2d vMean0(0);
        CVec2d vMean1(0);
        for(unsigned int i=0; i<uMatchCount; i++)
        {
            const CVec2f *pv0 = &pMatch[i].p0;
            const CVec2f *pv1 = &pMatch[i].p1;
            vMean0 += CVec2d(pv0->x, pv0->y);
            vMean1 += CVec2d(pv1->x, pv1->y);
        }
        // calculate the mean location of the point cloud
        // we could do this together with the second loop below
        // but there would be more limitation of machine precision
        // in accumulating x*x rather than dx*dx.
        vMean0 /= (double)uMatchCount;
        vMean1 /= (double)uMatchCount;

        CVec2d vU(0);
        double fSumSqXY = 0;

        for(unsigned int i=0; i<uMatchCount; i++)
        {
            // get information about rotation/translation
            // around the mean points
            CVec2d vDiff0 = CVec2d(pMatch[i].p0.x, pMatch[i].p0.y) - vMean0;
            CVec2d vDiff1 = CVec2d(pMatch[i].p1.x, pMatch[i].p1.y) - vMean1;
            fSumSqXY += vDiff0 * vDiff0;
            vU(0) += vDiff0 * vDiff1;
            vU(1) += vDiff0.y * vDiff1.x - vDiff0.x *vDiff1.y;
        }

        if(fSumSqXY < MIN_SQ_LENGTH * uMatchCount)
        {
            mS = 0; // case when all image 0 points very closely coincide
            mS(2,2) = 1;
            goto Exit;
        }

        mS(0,0) = vU(0)/fSumSqXY;
        mS(0,1) = vU(1)/fSumSqXY;
        mS(1,0) = -mS(0,1);
        mS(1,1) = mS(0,0);
        mS(0,2) = vMean1.x - mS(0,0) * vMean0.x - mS(0,1) * vMean0.y;
        mS(1,2) = vMean1.y - mS(1,0) * vMean0.x - mS(1,1) * vMean0.y;
        mS(2,0) = 0;
        mS(2,1) = 0;
        mS(2,2) = 1;
    }

Exit:
    return hr;
}

HRESULT vt::VtAffineFromPointMatches2D(OUT CMtx3x3d &mAff, 
                                       IN  const PointMatch* pMatch,
                                       unsigned int uMatchCount)
{
    HRESULT hr = NOERROR;
    if(uMatchCount<3)
        VT_HR_EXIT( E_INVALIDARG );
    if(uMatchCount==3)
    {
        CVec2d vOrigin0 = CVec2d(pMatch[0].p0.x, pMatch[0].p0.y);
        CVec2d vA0 = CVec2d(pMatch[1].p0.x, pMatch[1].p0.y) - vOrigin0;
        CVec2d vB0 = CVec2d(pMatch[2].p0.x, pMatch[2].p0.y) - vOrigin0;
        CVec2d vOrigin1 = CVec2d(pMatch[0].p1.x, pMatch[0].p1.y);
        CVec2d vA1 = CVec2d(pMatch[1].p1.x, pMatch[1].p1.y) - vOrigin1;
        CVec2d vB1 = CVec2d(pMatch[2].p1.x, pMatch[2].p1.y) - vOrigin1;

        double fDet = vA0.x * vB0.y - vA0.y * vB0.x;
        if(fDet==0)
        {
            mAff = 0; // case when all points are co-linear
            mAff(2,2) = 1;
            goto Exit;
        }

        double fS = 1.0/fDet;

        mAff(0,0) = fS * (vB0.y * vA1.x - vA0.y * vB1.x);
        mAff(0,1) = fS * (-vB0.x * vA1.x + vA0.x * vB1.x);
        mAff(1,0) = fS * (vB0.y * vA1.y - vA0.y * vB1.y);
        mAff(1,1) = fS * (-vB0.x * vA1.y + vA0.x * vB1.y);
        mAff(0,2) = vOrigin1.x - mAff(0,0) * vOrigin0.x - mAff(0,1) * vOrigin0.y;
        mAff(1,2) = vOrigin1.y - mAff(1,0) * vOrigin0.x - mAff(1,1) * vOrigin0.y;
        mAff(2,0) = 0;
        mAff(2,1) = 0;
        mAff(2,2) = 1;
    }
    else
    {
        // least squares estimate of affine
        // matrix parameterised by:
        // (a b c; d e f; 0 0 1)
        CVec2d vMean0(0);
        CVec2d vMean1(0);
        for(unsigned int i=0; i<uMatchCount; i++)
        {
            vMean0 += CVec2d(pMatch[i].p0.x, pMatch[i].p0.y);
            vMean1 += CVec2d(pMatch[i].p1.x, pMatch[i].p1.y);
        }
        // calculate the mean location of the point cloud
        // we could do this together with the second loop below
        // but there would be more limitation of machine precision
        // in accumulating x*x rather than dx*dx.
        vMean0 /= (double)uMatchCount;
        vMean1 /= (double)uMatchCount;

        CVec3d vSum0(0);
        CVec4d vSum1(0);

        for(unsigned int i=0; i<uMatchCount; i++)
        {
            CVec2d vDiff0 = CVec2d(pMatch[i].p0.x, pMatch[i].p0.y) - vMean0;
            CVec2d vDiff1 = CVec2d(pMatch[i].p1.x, pMatch[i].p1.y) - vMean1;
            vSum0(0) += vDiff0.x * vDiff0.x;
            vSum0(1) += vDiff0.x * vDiff0.y;
            vSum0(2) += vDiff0.y * vDiff0.y;
            vSum1(0) += vDiff0.x * vDiff1.x;
            vSum1(1) += vDiff0.y * vDiff1.x;
            vSum1(2) += vDiff0.x * vDiff1.y;
            vSum1(3) += vDiff0.y * vDiff1.y;
        }

        double fDet = vSum0(0) * vSum0(2) - vSum0(1) * vSum0(1);
        if(fDet==0)
        {
            mAff = 0; // case when all points are co-linear
            mAff(2,2) = 1;
            goto Exit;
        }

        double fS = 1.0/fDet;

        mAff(0,0) = fS * (vSum0(2) * vSum1(0) - vSum0(1) * vSum1(1));
        mAff(0,1) = fS * (-vSum0(1) * vSum1(0) + vSum0(0) * vSum1(1));
        mAff(1,0) = fS * (vSum0(2) * vSum1(2) - vSum0(1) * vSum1(3));
        mAff(1,1) = fS * (-vSum0(1) * vSum1(2) + vSum0(0) * vSum1(3));
        mAff(0,2) = vMean1.x - mAff(0,0) * vMean0.x - mAff(0,1) * vMean0.y;
        mAff(1,2) = vMean1.y - mAff(1,0) * vMean0.x - mAff(1,1) * vMean0.y;
        mAff(2,0) = 0;
        mAff(2,1) = 0;
        mAff(2,2) = 1;
    }

Exit:
    return hr;
}


// Calculate the appropriate transform from the list of point matches provided. The minimum number of 
// matches is 4 for Homography. This code follows the normalised direct linear transformation algorithm 
// (i.e., H_i = argmin sum_i ( ||p1_i - H_i * p0_i)||^2 )) which is given by Hartley and Zisserman 
// "Multiple View Geometry in Computer Vision" p91. 

HRESULT vt::VtHomographyFromPointMatches2D(
	OUT CMtx3x3d &mHom, 
	IN  const PointMatch* pMatch,
	unsigned int uMatchCount
	)
{
	HRESULT hr = NOERROR;
	if(uMatchCount < 40)
	{
		if(uMatchCount <= 4)
			mHom = vt::CMtx3x3d().MakeI();
		else
			hr = VtAffineFromPointMatches2D(mHom, pMatch, uMatchCount);	
	}
	else
	{
		// least squares estimate of affine matrix parameterized by:
		// (a b c; d e f; h i 1)
		CVec2f vMean0(0);
		CVec2f vMean1(0);
		for(unsigned int i=0; i<uMatchCount; i++)
		{
			vMean0 += pMatch[i].p0;
			vMean1 += pMatch[i].p1;
		}
		// calculate the mean location of the point cloud
		// we could do this together with the second loop below
		// but there would be more limitation of machine precision
		// in accumulating x*x rather than dx*dx.
		vMean0 /= (float)uMatchCount;
		vMean1 /= (float)uMatchCount;

		float dist_p0 = 0.0f, dist_p1 = 0.0f;
		for(unsigned int i=0; i<uMatchCount; i++)
		{
			CVec2f vDiff0 = pMatch[i].p0 - vMean0;
			CVec2f vDiff1 = pMatch[i].p1 - vMean1;
			dist_p0 += sqrt(vDiff0.x * vDiff0.x + vDiff0.y * vDiff0.y);
			dist_p1 += sqrt(vDiff1.x * vDiff1.x + vDiff1.y * vDiff1.y);
		}
		dist_p0 /= (float)uMatchCount;
		dist_p1 /= (float)uMatchCount;

		float scale_p0 = sqrtf(2.0f) / dist_p0;
		float scale_p1 = sqrtf(2.0f) / dist_p1;
		vt::CMtx3x3f T_p0, T_p1;
		T_p0(0, 0) = scale_p0;	T_p0(0, 1) = 0.0f;		T_p0(0, 2) = -scale_p0 * vMean0.x;	
		T_p0(1, 0) = 0.0f;		T_p0(1, 1) = scale_p0;	T_p0(1, 2) = -scale_p0 * vMean0.y;	
		T_p0(2, 0) = 0.0f;		T_p0(2, 1) = 0.0f;		T_p0(2, 2) = 1.0f;	

		T_p1(0, 0) = scale_p1;	T_p1(0, 1) = 0.0f;		T_p1(0, 2) = -scale_p1 * vMean1.x;	
		T_p1(1, 0) = 0.0f;		T_p1(1, 1) = scale_p1;	T_p1(1, 2) = -scale_p1 * vMean1.y;	
		T_p1(2, 0) = 0.0f;		T_p1(2, 1) = 0.0f;		T_p1(2, 2) = 1.0f;	

		vector<vt::PointMatch> new_pvm;
		new_pvm.resize(uMatchCount);
		for(unsigned int i = 0; i < uMatchCount; i ++)
		{
			new_pvm[i].p0.x = T_p0(0, 0) * pMatch[i].p0.x + T_p0(0, 2);
			new_pvm[i].p0.y = T_p0(1, 1) * pMatch[i].p0.y + T_p0(1, 2);
			new_pvm[i].p1.x = T_p1(0, 0) * pMatch[i].p1.x + T_p1(0, 2);
			new_pvm[i].p1.y = T_p1(1, 1) * pMatch[i].p1.y + T_p1(1, 2);
		}

		vt::CMtxd mAA;
		mAA.Create(3*uMatchCount, 9);
		for(unsigned int i = 0; i < uMatchCount; i ++)
		{
			int j3 = 3*i;

			double fxx = new_pvm[i].p1.x * new_pvm[i].p0.x;
			double fxy = new_pvm[i].p1.x * new_pvm[i].p0.y;
			double fyx = new_pvm[i].p1.y * new_pvm[i].p0.x;
			double fyy = new_pvm[i].p1.y * new_pvm[i].p0.y;

			mAA(j3,0) = 0;
			mAA(j3,1) = 0;
			mAA(j3,2) = 0;
			mAA(j3,3) = -new_pvm[i].p0.x;
			mAA(j3,4) = -new_pvm[i].p0.y;
			mAA(j3,5) = -1;
			mAA(j3,6) = fyx;
			mAA(j3,7) = fyy;
			mAA(j3,8) = new_pvm[i].p1.y;

			mAA(j3+1,0) = new_pvm[i].p0.x;
			mAA(j3+1,1) = new_pvm[i].p0.y;
			mAA(j3+1,2) = 1;
			mAA(j3+1,3) = 0;
			mAA(j3+1,4) = 0;
			mAA(j3+1,5) = 0;
			mAA(j3+1,6) = -fxx;
			mAA(j3+1,7) = -fxy;
			mAA(j3+1,8) = -new_pvm[i].p1.x;

			mAA(j3+2,0) = -fyx;
			mAA(j3+2,1) = -fyy;
			mAA(j3+2,2) = -new_pvm[i].p1.y;
			mAA(j3+2,3) = fxx;
			mAA(j3+2,4) = fxy;
			mAA(j3+2,5) = new_pvm[i].p1.x;
			mAA(j3+2,6) = 0;
			mAA(j3+2,7) = 0;
			mAA(j3+2,8) = 0;
		}
		vt::CSolveSVDd svd;
		// If failed, then the matrix is singular.
		VT_HR_EXIT( svd.Solve(mAA) );

		vt::CVecd vh = svd.GetBestNullSpaceVector();

		vt::CMtx3x3d mQ;
		for(int j=0; j<3; j++){
			for(int k=0; k<3; k++){
				mQ(j,k) = vh[j*3+k];
			}
		}

		mHom = T_p1.Inv() * mQ * T_p0;

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				mHom(i,j) /= mHom(2, 2);
	}
Exit:
	return hr;
}

// Calculate the appropriate transform from the list of point matches provided.
// The minimum number of matches is 2 for Similarity, 3 for Affine and 4 for Homography.
// If the number of matches provided exceeds the minimum requirement then the solution
// returned is the one that minimizes a geometric distance. The minimization assumes a circular
// Gaussian error on both the points p0 and p1. The error equation is
// E = sum_i ( ||p1_i - p1exact_i||^2 + ||p0_i - p0exact_i||^2 )
// If the pointer pvExact is not NULL, the function returns the points p0exact and p1exact
// which are related by p1exact = H p0exact.
HRESULT vt::VtSimilarityFromPointMatchesSymmetric2D(OUT CMtx3x3d &mS,
                                                    IN  const PointMatch* pMatch,
                                                    IN  unsigned int uMatchCount,
                                                    OUT PointMatch* pMatchExact)
{
    HRESULT hr = NOERROR;
    if(uMatchCount<2)
        VT_HR_EXIT( E_INVALIDARG );
    if(uMatchCount==2)
    {
        VT_HR_EXIT( VtSimilarityFromPointMatches2D(mS, pMatch, uMatchCount) );
        if(pMatchExact)
        {
            pMatchExact[0] = pMatch[0];
            pMatchExact[1] = pMatch[1];
        }
    }
    else
    {
        CVec2d vMean0(0.0);
        CVec2d vMean1(0.0);
        for(unsigned int i=0; i<uMatchCount; i++)
        {
            vMean0 += CVec2d(pMatch[i].p0.x, pMatch[i].p0.y);
            vMean1 += CVec2d(pMatch[i].p1.x, pMatch[i].p1.y);
        }
        vMean0 /= (double)uMatchCount;
        vMean1 /= (double)uMatchCount;

        CMtx4x4d mCov(0.0);
        for(unsigned int i=0; i<uMatchCount; i++)
        {
            CVec4d vX(pMatch[i].p0.x - vMean0.x, pMatch[i].p0.y - vMean0.y,
                pMatch[i].p1.x - vMean1.x, pMatch[i].p1.y - vMean1.y);
            mCov(0,0) += vX(0)*vX(0);
            mCov(0,1) += vX(0)*vX(1);
            mCov(0,2) += vX(0)*vX(2);
            mCov(0,3) += vX(0)*vX(3);
            mCov(1,1) += vX(1)*vX(1);
            mCov(1,2) += vX(1)*vX(2);
            mCov(1,3) += vX(1)*vX(3);
            mCov(2,2) += vX(2)*vX(2);
            mCov(2,3) += vX(2)*vX(3);
            mCov(3,3) += vX(3)*vX(3);
        }
        mCov.MakeSymmetric();

        CVec4d vU(1,0,1,0); // initialise to vector that is not likely to be in the null space of mCov.
        for(unsigned int i=0;i<20;i++)
        {
            // Iterate...
            // we expect that mCov has 2 large and 2 small eigenvalues.
            // iteration should either converge to the eigenvector with largest
            // eigenvalue or else cycle around vectors within the 2d space
            // spanned by the largest pair of eigenvectors which is ok since
            // any vector in that space gives the correct solution.
            vU = (mCov * vU).Unit();
//           vU.Dump();
//#undef printf
//          printf("%g\n", vU * (mCov * vU));
        }

        // could solve with svd but this is 20 times slower
        //CSolveSVDd svd(mCov);
        //svd.V().Dump();
        //vU = svd.V().GetCol(0);

        double fScale = 1.0/(vU(0) * vU(0) + vU(1) * vU(1));
        mS(0,0) = fScale * (vU(0) * vU(2) + vU(1) * vU(3));
        mS(0,1) = fScale * (vU(1) * vU(2) - vU(0) * vU(3));
        mS(1,0) = -mS(0,1);
        mS(1,1) = mS(0,0);
        mS(0,2) = vMean1.x - mS(0,0) * vMean0.x - mS(0,1) * vMean0.y;
        mS(1,2) = vMean1.y - mS(1,0) * vMean0.x - mS(1,1) * vMean0.y;
        mS(2,0) = 0;
        mS(2,1) = 0;
        mS(2,2) = 1;

        if(pMatchExact)
        {
            // use sampson error dX to adjust point matches
            // onto the variety.
            double fSq = mS(0,0)*mS(0,0) + mS(0,1)*mS(0,1);
            fScale = 1.0/(1.0 + fSq);
            CMtx4x4d mMap;
            mMap(0,0) = fScale;
            mMap(0,1) = 0;
            mMap(0,2) = fScale*mS(0,0);
            mMap(0,3) = -fScale*mS(0,1);
            mMap(1,1) = fScale;
            mMap(1,2) = -mMap(0,3);
            mMap(1,3) = mMap(0,2);
            mMap(2,2) = fScale*fSq;
            mMap(2,3) = 0;
            mMap(3,3) = mMap(2,2);
            mMap.MakeSymmetric();
            CVec4d vOff(vMean0.x, vMean0.y, vMean1.x, vMean1.y);
            vOff = vOff - mMap * vOff;
           
            for(unsigned int i=0; i<uMatchCount; i++)
            {
                CVec4d vX(pMatch[i].p0.x, pMatch[i].p0.y, pMatch[i].p1.x, pMatch[i].p1.y);
                vX = mMap * vX + vOff;
                pMatchExact[i].p0 = CVec2f((float)vX(0), (float)vX(1));
                pMatchExact[i].p1 = CVec2f((float)vX(2), (float)vX(3));
            }
        }
    }

Exit:
    return hr;
}

HRESULT vt::VtAffineFromPointMatchesSymmetric2D(OUT CMtx3x3d &mAff,
                                                IN  const PointMatch* pMatch,
                                                IN  unsigned int uMatchCount,
                                                OUT PointMatch* pMatchExact)
{
    HRESULT hr = NOERROR;
    if(uMatchCount<3)
        VT_HR_EXIT( E_INVALIDARG );
    if(uMatchCount==3)
    {
        VT_HR_EXIT( VtAffineFromPointMatches2D(mAff, pMatch, uMatchCount) );
        if(pMatchExact)
        {
            pMatchExact[0] = pMatch[0];
            pMatchExact[1] = pMatch[1];
            pMatchExact[2] = pMatch[2];
        }
    }
    else
    {
        CVec2d vMean0(0);
        CVec2d vMean1(0);
        for(unsigned int i=0; i<uMatchCount; i++)
        {
            vMean0 += CVec2d(pMatch[i].p0.x, pMatch[i].p0.y);
            vMean1 += CVec2d(pMatch[i].p1.x, pMatch[i].p1.y);
        }
        vMean0 /= (double)uMatchCount;
        vMean1 /= (double)uMatchCount;

        CMtx4x4d mCov(0);
        for(unsigned int i=0; i<uMatchCount; i++)
        {
            CVec4d vX(pMatch[i].p0.x - vMean0.x, pMatch[i].p0.y - vMean0.y, 
                      pMatch[i].p1.x - vMean1.x, pMatch[i].p1.y - vMean1.y);
            mCov(0,0) += vX(0)*vX(0);
            mCov(0,1) += vX(0)*vX(1);
            mCov(0,2) += vX(0)*vX(2);
            mCov(0,3) += vX(0)*vX(3);
            mCov(1,1) += vX(1)*vX(1);
            mCov(1,2) += vX(1)*vX(2);
            mCov(1,3) += vX(1)*vX(3);
            mCov(2,2) += vX(2)*vX(2);
            mCov(2,3) += vX(2)*vX(3);
            mCov(3,3) += vX(3)*vX(3);
        }
        mCov.MakeSymmetric();

        CSolveSVDd svd(mCov);
        if(svd.IsError())
            VT_HR_EXIT( svd.GetError() );

        //svd.V().Dump();
        // get largest two singular vectors
        // these span the data space
        // if the data is degenerate then the singular values will be close to 
        // zero for one or both of these vectors
        CVec4d v1(svd.V().GetCol(0));
        CVec4d v2(svd.V().GetCol(1));

        CMtx2x2d mA(v1(2), v2(2), v1(3), v2(3));
        CMtx2x2d mB(v1(0), v2(0), v1(1), v2(1));
        CMtx2x2d mH = mA * mB.Inv();

        mAff(0,0) = mH(0,0);
        mAff(0,1) = mH(0,1);
        mAff(1,0) = mH(1,0);
        mAff(1,1) = mH(1,1);
        mAff(0,2) = vMean1.x - mAff(0,0) * vMean0.x - mAff(0,1) * vMean0.y;
        mAff(1,2) = vMean1.y - mAff(1,0) * vMean0.x - mAff(1,1) * vMean0.y;
        mAff(2,0) = 0;
        mAff(2,1) = 0;
        mAff(2,2) = 1;

        if(pMatchExact)
        {
            // use sampson error dX to adjust point matches
            // onto the variety.
            CMtx4x4d mMap;
            CMtxd mJ(2,4);
            VT_HR_EXIT( mJ.GetError() );

            mJ.Update(0,0,mH);
            mJ(0,2) = -1;
            mJ(0,3) = 0;
            mJ(1,2) = 0;
            mJ(1,3) = -1;
            mMap = -mJ.T() * (mJ * mJ.T()).Inv() * mJ;
            mMap(0,0) += 1;
            mMap(1,1) += 1;
            mMap(2,2) += 1;
            mMap(3,3) += 1;

            CVec4d vOff(vMean0.x, vMean0.y, vMean1.x, vMean1.y);
            vOff = vOff - mMap * vOff;
           
            for(unsigned int i=0; i<uMatchCount; i++)
            {
                CVec4d vX(pMatch[i].p0.x, pMatch[i].p0.y, pMatch[i].p1.x, pMatch[i].p1.y);
                vX = mMap * vX + vOff;
                pMatchExact[i].p0 = CVec2f((float)vX(0), (float)vX(1));
                pMatchExact[i].p1 = CVec2f((float)vX(2), (float)vX(3));
            }
        }
    }

Exit:
    return hr;
}

//
// ifdefed out of the build because the implementation of these is not finished
//
#if 0
// Calculate the appropriate transform from the list of point matches provided.
// The minimum number of matches is 2 for Similarity, 3 for Affine and 4 for Homography.
// If the number of matches provided exceeds the minimum requirement then the solution
// returned is the one that minimizes a geometric distance. The inverse covariance matrices must
// be provided for both p0 and p1 for each match. The error
// equation is given by
// E = sum_i ( ||p1_i - p1exact_i)||^2_ic1 ) + sum_i ( ||p0_i - p0exact_i||^2_ic0 )
// where ||x||^2_icov = x' icov x.
// If the pointer pMatchExact is not NULL, the function returns the points p0exact and p1exact
// which are related by p1exact = H p0exact.

// point match structure with inverse covariances
typedef struct {
    PointMatch match;
    CMtx2x2f ic0, ic1;
} PointMatchCov;

HRESULT vt::VtSimilarityFromPointMatchesSymmetricCov2D(OUT CMtx3x3d &mS,
                                                       IN  const PointMatchCov* pMatch,
                                                       IN  unsigned int uMatchCount,
                                                       OUT PointMatch* pMatchExact)
{
    HRESULT hr = NOERROR;
    if(uMatchCount<2)
        VT_HR_EXIT( E_INVALIDARG );
    if(uMatchCount==2)
    {
        PointMatch vpm[2];
        vpm[0] = pMatch[0].match;
        vpm[1] = pMatch[1].match;
        VT_HR_EXIT( VtSimilarityFromPointMatches2D(mS, vpm, 2) );
        if(pMatchExact)
        {
            pMatchExact[0] = vpm[0];
            pMatchExact[1] = vpm[1];
        }
    }
    else
    {
        CVec2d vSum0 = 0;
        CMtx2x2d mSumIC0 = 0;
        CVec2d vSum1 = 0;
        CMtx2x2d mSumIC1 = 0;
        for(unsigned int i=0; i<uMatchCount; i++)
        {
            vSum0 += CVec2d(pMatch[i].ic0 * pMatch[i].match.p0);
            mSumIC0 += CMtx2x2d(pMatch[i].ic0);
            vSum1 += CVec2d(pMatch[i].ic1 * pMatch[i].match.p1);
            mSumIC1 += CMtx2x2d(pMatch[i].ic1);
        }
        CVec2d vMean0 = mSumIC0.Inv() * vSum0;
        CVec2d vMean1 = mSumIC1.Inv() * vSum1;


    }

    hr = E_NOTIMPL;
Exit:
    return hr;
}

HRESULT vt::VtAffineFromPointMatchesSymmetricCov2D(OUT CMtx3x3d &mA,
                                                   IN  const PointMatchCov* pMatch,
                                                   IN  unsigned int uMatchCount,
                                                   OUT PointMatch* pMatchExact)
{
    HRESULT hr = NOERROR;
    if(uMatchCount<3)
        VT_HR_EXIT( E_INVALIDARG );
    if(uMatchCount==3)
    {
        PointMatch vpm[3];
        vpm[0] = pMatch[0].match;
        vpm[1] = pMatch[1].match;
        vpm[2] = pMatch[2].match;
        VT_HR_EXIT( VtAffineFromPointMatches2D(mA, vpm, 3) );
        if(pMatchExact)
        {
            pMatchExact[0] = vpm[0];
            pMatchExact[1] = vpm[1];
            pMatchExact[2] = vpm[1];
        }
    }
    else
    {
    }

    hr = E_NOTIMPL;
Exit:
    return hr;
}
#endif
