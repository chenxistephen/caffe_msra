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


#include "stdafx.h"

#include <stdlib.h>
#include "vt_kmeans.h"

#include "vt_mathutils.h"
#include "vt_matrix.h"
#include "vt_mathutils.inl"

using namespace vt;

#define KMEANS_MOVE_THRESH (1e-12)

static int 
#if defined(MSRVT_WINDOWS_BUILD)
__cdecl
#endif
qsort_compare_increasing(const void *p1, const void *p2)
{
    float *pf1 = (float *)p1;
    float *pf2 = (float *)p2;
    return *pf1 > *pf2 ? 1 : -1;
}

HRESULT CKMeans::Init(int iDim, int k, int iRandSeed, bool bPCAInit, int iMaxIters)
{
    HRESULT hr = NOERROR;
    m_iK = k;
    m_iDim = iDim;
    m_bPCAInit = bPCAInit;

    if(iDim<1 || k<1)
        VT_HR_EXIT( E_INVALIDARG );

    VT_HR_EXIT( m_vTmp.Create(m_iDim) );
    VT_HR_EXIT( m_vecCI.resize(k) );

    int i;
    for(i=0; i<k; i++)
    {
        VT_HR_EXIT( m_vecCI[i].vMean.Create(m_iDim) );
        VT_HR_EXIT( m_vecCI[i].vMeanNew.Create(m_iDim) );
        VT_HR_EXIT( m_vecCI[i].mCov.Create(m_iDim, m_iDim) );
        VT_HR_EXIT( m_vecCI[i].mCovNew.Create(m_iDim, m_iDim) );
        VT_HR_EXIT( m_vecCI[i].mInvCov.Create(m_iDim, m_iDim) );
        m_vecCI[i].mCov.MakeI();
        m_vecCI[i].mInvCov.MakeI();
    }

    m_rnd.Seed(iRandSeed);

    m_iMaxIterations = iMaxIters;

Exit:
    return hr;
}

#undef printf

HRESULT CKMeans::Compute(float **pData, int iCount)
{ 
    HRESULT hr = NOERROR;
    m_bUseCov = false;
	// BSW:
	bool	bPrevNewCluster=false;
	int		iPrevNewClusterNum=0;

    int iIter = 0;

    VT_HR_EXIT( InitializeClusterCenters(pData, iCount) );

	while(iIter++ < m_iMaxIterations)
    {
        // determine which data are closest to each cluster
        int i, j;
        for(i=0; i<m_iK; i++)
        {
            m_vecCI[i].vMeanNew = 0;
            m_vecCI[i].iCount = 0;
        }
#if 0
        FILE *fp = fopen("cl.txt", "a");
        for(i=0; i<m_iK; i++)
        {
            int j;
            for(j=0; j<m_vecCI[i].vMean.Size(); j++)
              fprintf(fp, "%02g ", m_vecCI[i].vMean[j]);
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n\n");
        fclose(fp);
#endif
        m_fSumSq = 0;
        for(i=0; i<iCount; i++)
        {
            float fDist;
            int iCent = FindClosestCenter(pData[i], &fDist);
            m_fSumSq += fDist * fDist;
            // accumulate
            m_vecCI[iCent].iCount++;
            for(j=0; j<m_iDim; j++)
                m_vecCI[iCent].vMeanNew[j] += pData[i][j];
        }
        m_fSumSq /= iCount;

        // determine if we moved
        bool bMoved = false;
        bool bZero = false;
        for(i=0; i<m_iK; i++)
        {
            if(m_vecCI[i].iCount==0)
            {
                if(!bZero)
                {
                    // this center has dropped out, so reinitialize.
                    // we reinit to a location which is the centroid of the
                    // data weighted by the distance to the closest cluster centers
                    double fWeight = 0;
                    int k;
                    for(k=0; k<iCount; k++)
                    {
                        float fDist;
                        FindClosestCenter(pData[k], &fDist);
                        fWeight += fDist;
                        for(j=0; j<m_iDim; j++)
                            m_vecCI[i].vMeanNew[j] += fDist * pData[k][j];
                    }
                    for(j=0; j<m_iDim; j++)
                        m_vecCI[i].vMean[j] = (float)(m_vecCI[i].vMeanNew[j] / fWeight);

					// BSW:
					if(!bPrevNewCluster || iPrevNewClusterNum != i){
	                    bMoved = true; // prevent early exit
					}
					bPrevNewCluster = true;
					iPrevNewClusterNum = i;
                    bZero = true; // if some other one needs reinit, do it on the next iteration
                }
				// BSW: ignore 2nd zero member cluster. 
                //else
                //    bMoved = true;
            }
            else
            {
                double fscale = 1.0/m_vecCI[i].iCount;
                if(!bMoved)
                {
                    // test and update
                    double fsum = 0;
                    for(j=0; j<m_iDim; j++)
                    {
                        float fnew = (float)(fscale * m_vecCI[i].vMeanNew[j]);
                        float fd = fnew - m_vecCI[i].vMean[j];
                        fsum += (double)fd * (double)fd;
                        m_vecCI[i].vMean[j] = fnew;
                    }
                    if(fsum > KMEANS_MOVE_THRESH)
                        bMoved = true;
                }
                else
                {
                    // just update
                    for(j=0; j<m_iDim; j++)
                        m_vecCI[i].vMean[j] = (float)(fscale * m_vecCI[i].vMeanNew[j]);
                }
            }
        }
        if(!bMoved)
            break;
		if(!bZero){
			bPrevNewCluster = false;
		}
    }

Exit:
    return hr;
}

HRESULT CKMeans::ComputeWithVariances(float **pData, int iCount, float fMinVar)
{
    HRESULT hr = NOERROR;
    
	int iIter = 0;
    
	VT_HR_EXIT( InitializeClusterCenters(pData, iCount) );
    m_bUseCov = true;

    int i;
    for(i=0; i<m_iK; i++)
        m_vecCI[i].bDropped = false;

    while(iIter++ < m_iMaxIterations)
    {
        //printf("------------------------------\n");
        // determine which data are closest to each cluster
        int j;
        for(i=0; i<m_iK; i++)
        {
            m_vecCI[i].vMeanNew = 0;
            m_vecCI[i].iCount = 0;
            m_vecCI[i].mCovNew = 0;
        }

        m_fSumSq = 0;
        for(i=0; i<iCount; i++)
        {
            float fDist;
            int iCent = FindClosestCenter(pData[i], &fDist);
            m_fSumSq += fDist * fDist;
            // accumulate mean and covariance of data points
            m_vecCI[iCent].iCount++;
            for(j=0; j<m_iDim; j++)
            {
                double f = pData[i][j];
                m_vecCI[iCent].vMeanNew[j] += f;
                int k;
                for(k=j; k<m_iDim; k++)
                    m_vecCI[iCent].mCovNew(j,k) += f * (double)pData[i][k];
            }
        }
        m_fSumSq /= iCount;

        bool bMoved = false;
        bool bZero = false;
        for(i=0; i<m_iK; i++)
        {
            if(m_vecCI[i].iCount==0)
            {
                if(!m_vecCI[i].bDropped) // if just previously dropped, give up
                {
                    if(!bZero)
                    {
                        m_vecCI[i].bDropped = true;
                        // this center has dropped out, so reinitialize.
                        // we reinit to a location which is the centroid of the
                        // data weighted by the distance to the closest cluster centers
                        double fWeight = 0;
                        int k;
                        for(k=0; k<iCount; k++)
                        {
                            float fDist;
                            FindClosestCenter(pData[k], &fDist);
                            fWeight += fDist;
                            for(j=0; j<m_iDim; j++)
                                m_vecCI[i].vMeanNew[j] += fDist * pData[k][j];
                        }
                        for(j=0; j<m_iDim; j++)
                            m_vecCI[i].vMean[j] = (float)(m_vecCI[i].vMeanNew[j] / fWeight);

                        // re-init covariance to minimum variance
                        m_vecCI[i].mCov = 0;
                        m_vecCI[i].mInvCov = 0;
                        for(j=0; j<m_iDim; j++)
                        {
                            m_vecCI[i].mCov(j,j) = fMinVar;
                            m_vecCI[i].mInvCov(j,j) = 1/fMinVar;
                        }

                        //bMoved = true; // prevent early exit
                        bZero = true; // if some other one needs reinit, do it on the next iteration
                    }
                    else
                        bMoved = true;
                }
            }
            else
            {
                m_vecCI[i].bDropped = false;

                double fscale = 1.0/m_vecCI[i].iCount;
                if(!bMoved)
                {
                    // test and update
                    double fsum = 0;
                    for(j=0; j<m_iDim; j++)
                    {
                        float fnew = (float)(fscale * m_vecCI[i].vMeanNew[j]);
                        float fd = fnew - m_vecCI[i].vMean[j];
                        fsum += (double)fd * (double)fd;
                        m_vecCI[i].vMean[j] = fnew;
                    }
                    if(fsum > KMEANS_MOVE_THRESH)
                        bMoved = true;
                }
                else
                {
                    // just update
                    for(j=0; j<m_iDim; j++)
                        m_vecCI[i].vMean[j] = (float)(fscale * m_vecCI[i].vMeanNew[j]);
                }
                
                int k;
                // recalculate variances
                for(j=0; j<m_iDim; j++)
                    for(k=j; k<m_iDim; k++)
                        m_vecCI[i].mCov(j,k) = (float)(m_vecCI[i].mCovNew(j,k) * fscale 
                            - (double)m_vecCI[i].vMean[j] * (double)m_vecCI[i].vMean[k]);

                m_vecCI[i].mCov.MakeSymmetric();

                //printf("mean %d - %d\n", i, m_vecCI[i].iCount);m_vecCI[i].vMean.Dump();
                //printf("Cov:\n");m_vecCI[i].mCov.Dump();

                VT_HR_EXIT( ConditionCovariance(i, fMinVar) );

                //printf("Covnew:\n");m_vecCI[i].mCov.Dump();
                //printf("InvCovnew:\n");m_vecCI[i].mInvCov.Dump();
            }
        }

        if(!bMoved)
            break;
    }

Exit:
    return hr;
}

HRESULT CKMeans::ConditionCovariance(int i, float fMinVar)
{
    HRESULT hr = NOERROR;
    // prevent covariance from collapsing

    VT_HR_EXIT( m_svd.Solve(m_vecCI[i].mCov) );

    //printf("D:\n");m_svd.D().Dump();
    //printf("V:\n");m_svd.V().Dump();

    int j, k;
    for(j=0; j<m_iDim; j++)
    {
        float fd = m_svd.D()(j,j);
        fd = VtMax(fd, fMinVar);
        for(k=0; k<m_iDim; k++)
        {
            // reconstitute and compute inverse
            m_vecCI[i].mCov(j,k) = fd * m_svd.V()(k,j);
            m_vecCI[i].mInvCov(j,k) = 1/fd * m_svd.V()(k,j);
        }
    }
    m_vecCI[i].mCov = m_svd.U() * m_vecCI[i].mCov;
    m_vecCI[i].mInvCov = m_svd.U() * m_vecCI[i].mInvCov;
    VT_HR_EXIT( m_vecCI[i].mCov.GetError() );
    VT_HR_EXIT( m_vecCI[i].mInvCov.GetError() );
Exit:
    return hr;
}

HRESULT CKMeans::InitializeClusterCenters(float **pData, int iCount)
{
    HRESULT hr = NOERROR;
    if(iCount<m_iK)
        VT_HR_EXIT( E_INVALIDARG );
    if(pData==NULL)
        VT_HR_EXIT( E_POINTER );

    if(!m_bPCAInit)
    {
        CVecf vMin, vMax;
        VT_HR_EXIT( vMin.Create(m_iDim) );
        VT_HR_EXIT( vMax.Create(m_iDim) );
        vMin = VT_MAXFLOAT;
        vMax = -VT_MAXFLOAT;

        int i, j;
        int iC = VtMin(10000, iCount);

        for(i=0; i<iC; i++)
        {
            int ind;
            if(iC<iCount)
                ind = m_rnd.IRand(iCount); // pick a random subset
            else
                ind = i;
            for(j=0; j<m_iDim; j++)
            {
                float f = pData[ind][j];
                if(f>vMax[j])
                    vMax[j] = f;
                if(f<vMin[j])
                    vMin[j] = f;
            }
        }

        for(i=0; i<m_iK; i++)
        {
            for(j=0; j<m_iDim; j++)
                m_vecCI[i].vMean[j] = (float)m_rnd.URand(vMin[j], vMax[j]);
        }
    }
    else
    {
        CVecf vMean, vU, vD;
        CMtxf mCov;
        // find approx covariance of data
        VT_HR_EXIT( VtMeanAndCovariance(pData, m_iDim, iCount, 10000, vMean, mCov) );
        VT_HR_EXIT( vU.Create(m_iDim) );
        VT_HR_EXIT( vD.Create(m_iDim) );

        // quick method to find dominant eigenvector
        int i, k;
        for(i=0; i<m_iDim; i++)
            vU[i] = (float)m_rnd.URand(-1, 1);
        for(i=0; i<20; i++)
            vU = (mCov * vU).Unit();
        vector<float> vecProj;
        int iC = VtMin(VtMax(10000, 2*m_iK), iCount);
        VT_HR_EXIT( vecProj.resize(iC) );

        // project subset of data onto the unit vector
        for(i=0; i<iC; i++)
        {
            int ind;
            if(iC<iCount)
                ind = m_rnd.IRand(iCount); // pick a random subset
            else
                ind = i;
            for(k=0; k<m_iDim; k++)
                vD[k] = pData[ind][k] - vMean[k];
            vecProj[i] = vD * vU;
        }

        qsort(&(vecProj[0]), iC, sizeof(float), qsort_compare_increasing);
    
        float d = (iC-1)/(float)m_iK;
        for(i=0; i<m_iK; i++)
        {
            int ind = (int)((0.5f + i) * d + 0.5f);
            ind = VtMin(iC-1, ind);
            float fProj = vecProj[ind];
            m_vecCI[i].vMean = vU * fProj + vMean;
        }
    }

Exit:
    return hr;
}

int CKMeans::FindClosestCenter(float *pData, float *pDist)
{
    int i;
    float fMinDist = VT_MAXFLOAT;
    int iClosest = 0;
    for(i=0; i<m_iK; i++)
    {
        float fDist = GetDistanceToCenter(pData, i);
        if(fDist < fMinDist)
        {
            fMinDist = fDist;
            iClosest = i;
        }
    }
    if(pDist!=NULL)
        *pDist = fMinDist;

    return iClosest;
}

float CKMeans::GetDistanceToCenter(float *pData, int i)
{
    float fDist = 0;
    if(m_bUseCov)
    {
        int r, c;
        for(r=0; r<m_iDim; r++)
            m_vTmp[r] = pData[r] - m_vecCI[i].vMean[r];

        for(r=0; r<m_iDim; r++)
        {
            float f = 0;
            for(c=0; c<m_iDim; c++)
                f += m_vecCI[i].mInvCov(r, c) * m_vTmp[c];
            fDist += f * m_vTmp[r];
        }
    }
    else
    {
        int r;
        for(r=0; r<m_iDim; r++)
        {
            float fdiff = pData[r] - m_vecCI[i].vMean[r];
            fDist += fdiff * fdiff;
        }
    }
    return sqrt(fDist);
}

