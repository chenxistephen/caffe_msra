//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Class to do k-means clustering
//
//  History:
//      2005/8/22-swinder
//			Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_rand.h"
#include "vt_solve_svd.h"

namespace vt {

typedef struct {
    CMtxf mCov;
    CMtxf mInvCov;
    CMtxd mCovNew;
    CVecf vMean;
    CVecd vMeanNew;
    int iCount;
    bool bDropped;
} ClusterInfo;

class CKMeans
{
public:
    CKMeans() : m_bUseCov(false) {}
    ~CKMeans() {}

    HRESULT Init(int iDim, int k, int iRandSeed, bool m_bPCAInit, int iMaxIters = 1000000);
    HRESULT Compute(float **pData, int iCount);
    HRESULT ComputeWithVariances(float **pData, int iCount, float fMinVar = 1e-4);
    float *GetClusterCenter(int i) { return m_vecCI[i].vMean.Ptr(); }
    void SetClusterCenter(int i, float *p) { memcpy(m_vecCI[i].vMean.Ptr(), p, m_iDim * sizeof(float)); }
    void GetCovariance(int i, CMtxf &mRtn) { mRtn = m_vecCI[i].mCov; }
    int FindClosestCenter(float *pData, float *pfDistRtn = NULL);
    float GetDistanceToCenter(float *pData, int i);
    int GetNumberOfPoints(int i) { return m_vecCI[i].iCount; }
    double GetTotalVariance() { return m_fSumSq; } // total data sum squared distances to cluster centers

protected:
    HRESULT InitializeClusterCenters(float **pData, int iCount);
    HRESULT ConditionCovariance(int i, float fMinVar);

    vector<ClusterInfo> m_vecCI;
    int m_iDim;
    int m_iK;
    CVecf m_vTmp;
    CRand m_rnd;
    bool m_bPCAInit;
    bool m_bUseCov;
    CSolveSVDf m_svd;
    int m_iMaxIterations;
    double m_fSumSq;
};

};
