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
#pragma once

#include "vtcommon.h"

#include "vt_matrix2x2.h"
#include "vt_matrix3x3.h"

namespace vt {

	// point match structure
	typedef struct {
		CVec2f p0, p1;
	} PointMatch;

	// Implementation Notes
	//
	// fully implemented:
	// VtSimilarityFromPointMatches2D(const vector<PointMatch> &vm, CMtx3x3f &mS);
	// VtAffineFromPointMatches2D(const vector<PointMatch> &vm, CMtx3x3d &mA);
	// VtSimilarityFromPointMatchesSymmetric2D(const vector<PointMatch> &vm, CMtx3x3d &mS, vector<PointMatch> *pvExact);
	// VtAffineFromPointMatchesSymmetric2D(const vector<PointMatch> &vm, CMtx3x3d &mA, vector<PointMatch> *pvExact);
	// all the rest are not implented

	// Calculate the appropriate transform from the list of point matches provided.
	// The minimum number of matches is 2 for Similarity, 3 for Affine and 4 for Homography.
	// If the number of matches provided exceeds the minimum requirement then the solution
	// returned is correct in the least squares sense. The minimization assumes a circular
	// Gaussian error on each point p1 and exact measurements for p0. The error equation is
	// given by E = sum_i ( ||p1_i - transformed(p0_i)||^2 )
	HRESULT VtSimilarityFromPointMatches2D(OUT CMtx3x3d &mS, 
		IN  const PointMatch* pvm, 
		IN  unsigned int uMatchCount);

	HRESULT VtAffineFromPointMatches2D(OUT CMtx3x3d &mA,
		IN  const PointMatch* pvm, 
		IN  unsigned int uMatchCount);

	// Calculate the appropriate transform from the list of point matches provided. The minimum number of 
	// matches is 4 for Homography. This code follows the normalised direct linear transformation algorithm 
	// (i.e., H_i = argmin sum_i ( ||p1_i - H_i * p0_i)||^2 )) which is given by Hartley and Zisserman 
	// "Multiple View Geometry in Computer Vision" p91. 
	HRESULT VtHomographyFromPointMatches2D(OUT CMtx3x3d &mS, 
		IN  const PointMatch* pvm, 
		IN  unsigned int uMatchCount);


	// Calculate the appropriate transform from the list of point matches provided.
	// The minimum number of matches is 2 for Similarity, 3 for Affine and 4 for Homography.
	// If the number of matches provided exceeds the minimum requirement then the solution
	// returned is the one that minimizes a geometric distance. The minimization assumes a circular
	// Gaussian error on both the points p0 and p1. The error equation is
	// E = sum_i ( ||p1_i - p1exact_i||^2 + ||p0_i - p0exact_i||^2 )
	// If the pointer pvExact is not NULL, the function returns the points p0exact and p1exact
	// which are related by p1exact = H p0exact.
	HRESULT VtSimilarityFromPointMatchesSymmetric2D(OUT CMtx3x3d &mS, 
		IN  const PointMatch* pvm, 
		IN  unsigned int uMatchCount, 
		OUT PointMatch* pvExact=NULL);

	HRESULT VtAffineFromPointMatchesSymmetric2D(OUT CMtx3x3d &mA, 
		IN  const PointMatch* pvm, 
		IN  unsigned int uMatchCount, 
		OUT PointMatch* pvExact=NULL);

	// transfers a point and inverse covariance about that point through a homography
	// the inverse covariance is transferred according to the jacobian of the homography around the point
	//
	// pdst ~ H [psrc.x; psrc.y; 1]
	// invCovDst = Jinv' invCovSrc Jinv
	template <class T>
	void VtTransferPoint(const CVec2<T> &vPointSrc, const CMtx2x2<T> &mInvCovSrc, const CMtx3x3<T> &mH, 
		CVec2<T> &vPointDst, CMtx2x2<T> &mInvCovDst)
	{
		CVec3<T> p = mH * CVec3<T>(vPointSrc, (T)1);
		vPointDst = p.Dehom();
		CMtx2x2<T> mJ((mH(0,0) - p.x * mH(2,0)/p.z)/p.z,
			(mH(0,1) - p.x * mH(2,1)/p.z)/p.z,
			(mH(1,0) - p.y * mH(2,0)/p.z)/p.z,
			(mH(1,1) - p.y * mH(2,1)/p.z)/p.z);
		mJ = mJ.Inv();
		mInvCovDst = mJ.T() * mInvCovSrc * mJ;
	}


};
