//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description: implementation of IFeatureWarpCompute that uses a 
//               dynamic-programming style smoothing algo on similarity parameters
//
//------------------------------------------------------------------------------
#include "stdafx.h"

#include "vt_smoother_dp.h"

HRESULT CDPPathSmoother::Begin()
{
    VT_HR_BEGIN()

    // allocate the rolling buffers
    VT_HR_EXIT( m_bufTrackerSim.resize(m_params.iBufferSize) );
    VT_HR_EXIT( m_bufResults.resize(m_params.iBufferSize) );

    VT_HR_END()
}

BUFFER_RANGE CDPPathSmoother::GetResultsRange()
{
    BUFFER_RANGE r;
    r.frame_count = m_bufResults.get_available_count();
    r.first_frame = m_bufResults.get_first_id();
    return r;
}

HRESULT CDPPathSmoother::PushFeatures(const vt::PointMatch* pMatches,
                                      int iMatchCount)
{
    VT_HR_BEGIN()

    // compute the transformation
    // TODO: add a set of transform functions for these and those left in the stabilizer
    CMtx3x3d mtxCur_d;
    if ( iMatchCount == 0 )
    {
        mtxCur_d.MakeI();
    }
    else
    {
        if (m_params.motionModel == 2)
        {
            VT_HR_EXIT( VtSimilarityFromPointMatches2D(mtxCur_d, pMatches, iMatchCount) );
        }
        else
        {
            VT_HR_EXIT( VtAffineFromPointMatches2D(mtxCur_d, pMatches, iMatchCount) );
        }
    }

    // all internal processing is done with coordinate system (0,0) at 
    // center of the image so adjust the returned matrix here
    CMtx3x3f mtxCur;
    VtConvertMtx(mtxCur_d, mtxCur);
    mtxCur = MatrixTopLeftToCenter(m_params.frameWidth, m_params.frameHeight) *
        mtxCur * MatrixCenterToTopLeft(m_params.frameWidth, m_params.frameHeight);

    // compute the similarity parameters for this frame and add them to a 
    // rolling buffer
    CMtx3x3f mtxToPrev = (iMatchCount==0 || m_bufTrackerSim.get_total_count()==0)? 
        CMtx3x3f().MakeI(): mtxCur;

    // compute the similarity params
    // TODO: could compute sim parameters directly when tracker in
    //       similarity mode
	m_bufTrackerSim.advance();
    ComputeSimParams(m_bufTrackerSim.last(), mtxToPrev, 
                     m_params.frameWidth, m_params.frameHeight);

#ifdef PATH_CSV
    if ( iSrcFrame == 0 )
    {
        _wfopen_s(&g_fpDbg,L"framedbg.csv",L"w");
    }

    if ( g_fpDbg )
    {
        DBG_PATH_INFO dd;
        dd.sp   = sim;
        dd.iRef = (iMatchCount==0 || m_iSrcCount==0)? -1: m_iSrcCount - 1;
        g_vecSrcDbg.push_back(dd);
    }
#endif

    if ( (int)m_bufTrackerSim.get_total_count() >= m_params.iBufferSize )
    {
        m_bufResults.advance();
        VT_HR_EXIT( DPFilterTransform(
            m_bufResults.last(), m_bufResults.get_last_id()) );
    }

    VT_HR_END()
}

HRESULT
CDPPathSmoother::GetResult(vt::CMtx3x3f& xfrm, int frameNumber)
{
    // results are -1 -> m_iDstCount-1
    if( !GetResultsRange().InRange(frameNumber) )
    {
        return E_INVALIDARG;
    }
    xfrm = m_bufResults[frameNumber];
    return S_OK;
}

HRESULT CDPPathSmoother::End()
{
	VT_HR_BEGIN()
	
	while( m_bufResults.get_total_count() < m_bufTrackerSim.get_total_count() )
	{
        m_bufResults.advance();
        VT_HR_EXIT( DPFilterTransform(
            m_bufResults.last(), m_bufResults.get_last_id()) );
	}

#ifdef PATH_CSV
	PrintCSV();
    if( g_fpDbg ) fclose(g_fpDbg);
#endif

	VT_HR_END()
}

//+-----------------------------------------------------------------------------
//
// Implementation of path smoother that finds piece-wise linear camera motions
//
//------------------------------------------------------------------------------

float g_fMaxDist = 0.667f;
short g_iMinFrameCountForPhase = 4;

typedef float (*T_pfnUpdatePathState)(const FEAT_SIMILARITY&, int, int);

struct OOB_VAL
{
    int   index; 
	float diff;
	float distB;
	float acc;
	short acc_cnt;
	short skip;

	OOB_VAL() : index(-1), acc(0), acc_cnt(0), skip(0)
	{}
};

#ifdef PATH_CSV
inline OOB_VAL testOOBPrint(float fD0, float fV0, float fA, int iACnt, int iCnt,
							int iSkip, const CRollingBuffer<FEAT_SIMILARITY>& srcpath,
							int idx, float fMax, int iW, int iH,
							T_pfnUpdatePathState pfnUpdate)
{
    OOB_VAL r;

    float fDB = 0, fDV = 0, fV = fV0;
    for( int i = 0; i < iCnt; i++ )
    {
        if( i < iACnt )
        {
            fV += fA;
        }
		fDB += pfnUpdate(srcpath[idx+i], iW, iH); 
		fDV += fV;

		float diff = fD0 + fDB - fDV;
		printf("[%d] dist %f\n", i, diff);
        if( fabs(diff) > fMax && i >= iSkip )
        {
            r.index = i;
            r.diff  = diff;
			r.distB = fDB;
            return r;
        }
    }

    r.index = -1;
	r.distB = fDB;
	r.diff  = fD0 + fDB - fDV;
    return r;
}

inline void PrintTests(const vt::vector<OOB_VAL>& vecTests,
					   float fD0, float fV0, int iACnt, int iCnt,
					   const CRollingBuffer<FEAT_SIMILARITY>& srcpath,
					   int idx, float fMax, int iW, int iH,
					   T_pfnUpdatePathState pfnUpdate)
{
	for( int k = 0; k < iCnt; k++ )
	{
		printf("%d", k);
		for( int j = 0; j < (int)vecTests.size(); j++ )
		{
			float fA = vecTests[j].acc;
		
			float fDB = 0, fDV = 0, fV = fV0;
			for( int i = 0; i < k+1; i++ )
			{
				if( i < iACnt )
				{
					fV += fA;
				}
				fDB += pfnUpdate(srcpath[idx+i], iW, iH); 
				fDV += fV;
			}

			if( j==0 )
			{
				printf(", %f", fD0+fDB);
			}
			printf(", %f", fDV);
		}
		printf("\n");
	}
}
#endif // PATH_CSV

//
// Update velocity and distance traveled for the various similarity components
//
float DistTransX(const FEAT_SIMILARITY& curp, int iW, int iH)
{ return curp.tx/float(iW); }

float DistTransY(const FEAT_SIMILARITY& curp, int iW, int iH)
{ return curp.ty/float(iH); }

float DistRotation(const FEAT_SIMILARITY& curp, int iW, int iH)
{ return curp.r; }

float DistScale(const FEAT_SIMILARITY& curp, int iW, int iH)
{ return log(curp.s); }

//
// Update the schedule for the similarity components
//

// return the sum of numbers from 1 to 'n'
inline int sum1_n(int n)
{ return (n&1)? sum1_n(n-1)+n: (1+n)*(n/2); }

// compute an new velocity 'v1' that results in traveling distance 'd' after a
// change in velocity from 'v0' with a constant acceleration lasting 'as' steps
// followed by traveling 'v1' over 'ts-as' steps 
inline float newV(float v0, int as, int ts, float d)
{
    if( ts <= as )
    {
        float sc = sum1_n(ts) / float(as);
        return (d-(ts-sc)*v0) / sc;
    }
    else
    {
        float sc = sum1_n(as) / float(as);
        return (d-(as-sc)*v0) / (float(ts-as)+sc);
    }
}

inline bool OOB1_Better(const OOB_VAL& o1, int i, const OOB_VAL& o2)
{
	int i1 = (o1.index<0)? -1: VtMax(0, o1.index-(o1.acc_cnt-g_iMinFrameCountForPhase));
	int i2 = (o2.index<0)? -1: VtMax(0, o2.index-(o2.acc_cnt-g_iMinFrameCountForPhase));

	return (i2 != -1 && i1 == -1)  || (i2 != -1 && i2 < i1+i) ||
		   (fabs(o2.diff) > fabs(o1.diff) &&((i2 == -1 && i1 == -1) || (i2==i1+i))) ||
		   (o1.skip < o2.skip);
}

inline OOB_VAL testOOB(float fD0, float fV0, float fA, int iACnt, int iCnt,
					   int iSkip, const CRollingBuffer<FEAT_SIMILARITY>& srcpath,
					   int frame, float fMax, int iW, int iH,
					   T_pfnUpdatePathState pfnUpdate)
{
	OOB_VAL r;

	int iAvgEndCnt = VtMin(5, iCnt);
	float fAvgDiff = 0;

    float fDB = 0, fDV = 0, fV = fV0;
    for( int i = 0; i < iCnt; i++ )
    {
        if( i < iACnt )
        {
            fV += fA;
        }
		fDB += pfnUpdate(srcpath[frame+i], iW, iH); 
		fDV += fV;

		float diff = fD0 + fDB - fDV;
        if( fabs(diff) > fMax && i >= iSkip )
        {
            r.index = i;
            r.diff  = diff;
			r.distB = fDB;
            return r;
        }

		if( i >= iCnt-iAvgEndCnt )
		{
			fAvgDiff += diff;
		}
    }

    r.index = -1;
	r.distB = fDB;
	r.diff  = fabs(fAvgDiff / float(iAvgEndCnt));
    return r;
}

PATH_STATE ChangeScheduleInternal(const PATH_STATE& init,
								  const CRollingBuffer<FEAT_SIMILARITY>& srcpath, 
                                  int frame, int iFrameCnt, int iW, int iH, float fMax,
								  T_pfnUpdatePathState pfnUpdate, bool bPrintTests = false)
{
	PATH_STATE result = init;

	// count down on the remaining steps in this phase
	if ( init.iStartupCount && init.fAcceleration != 0 )
	{
		result.fVelocity += result.fAcceleration;
		result.fDistance += pfnUpdate(srcpath[frame], iW, iH) - result.fVelocity;
		if ( --result.iStartupCount == 0 )
		{
			// 
			// TODO: if in acceleration mode and current schedule under 
			// shoots, then modify the acceleration schedule
			//

			// clamp velocity to 0
			if( fabs(result.fVelocity) < 0.0001f )
			{
				result.fVelocity = 0.f;
			}

			// switch from acceleration to constant velocity
			result.fAcceleration = 0;
		}
		return result;
	}

	// test what happens if we stay on the current path
	OOB_VAL oobNoChange = testOOB(init.fDistance, init.fVelocity, 0, 0, 
								  iFrameCnt, 0, srcpath, frame, fMax, iW, iH, 
								  pfnUpdate);

	// compute the total distance traveled to the final sample in the buffer
	float fBaseDist = 0;
	for( int i = 0; i < iFrameCnt; i++ )
	{
		fBaseDist += pfnUpdate(srcpath[frame+i], iW, iH);
	}

	// test what happens if we plot a path to zero the distance at the last
	// frame
	float fTestVelocity = newV(init.fVelocity, g_iMinFrameCountForPhase,
							   iFrameCnt, init.fDistance+fBaseDist);
	float fTestAccel = 
		(fTestVelocity-init.fVelocity) / float(g_iMinFrameCountForPhase);
	OOB_VAL oobChange = testOOB(init.fDistance, init.fVelocity, fTestAccel, 
								g_iMinFrameCountForPhase, iFrameCnt, 0, 
								srcpath, frame, fMax, iW, iH, pfnUpdate);

	// if staying on the current path goes OOB or changing the path is
	// advantageous  
	if( oobNoChange.index != -1 || 
		(result.fVelocity!=0 && oobNoChange.diff > fMax*0.1f && oobChange.diff < fMax*0.1f) )
	{
		OOB_VAL oob0;

#ifdef PATH_CSV
		vt::vector<OOB_VAL> vecTests;
#endif
		// test that waiting one or more cycles doesn't yield a better
		// result
		float fD0 = init.fDistance;
		int test_cnt = oobNoChange.index==-1? 
			g_iMinFrameCountForPhase: VtMax(1,oobNoChange.index);
		test_cnt = VtMin(test_cnt, iFrameCnt);
		for( int i = 0; i < test_cnt; i++ )
		{
			// test accelerating for g_iMinFrameCountForPhase to
			// 2*g_iMinFrameCountForPhase cycles
			//short acc_cnt = g_iMinFrameCountForPhase;
			for( short acc_cnt = g_iMinFrameCountForPhase; 
				 acc_cnt <= 2*g_iMinFrameCountForPhase; acc_cnt++)
			{
				int iSkip = 0;
restart:
				fTestVelocity = newV(init.fVelocity, acc_cnt, iFrameCnt-i, 
									 init.fDistance+fBaseDist - i*init.fVelocity);
				
				// test the new velocity - if it goes out of bounds somewhere
				// in the visible next frames then adjust to avoid the OOB 
				// condition 
				OOB_VAL oob, lastoob;
				while( 1 )
				{
					// if velocity is near zero then clamp
					if( fabs(fTestVelocity) < 0.0005f )
					{
						fTestVelocity = 0.0f;
					}

					fTestAccel = (fTestVelocity-init.fVelocity) / float(acc_cnt);
					oob = testOOB(fD0, init.fVelocity, fTestAccel, acc_cnt,
								  iFrameCnt-i, iSkip, srcpath, frame+i, fMax,
								  iW, iH, pfnUpdate); 
					oob.acc = fTestAccel;
					oob.acc_cnt = acc_cnt;
					oob.skip    = (short)iSkip;

#ifdef PATH_CSV
					if( i==0 )
					{
						vecTests.push_back(oob);
					}
#endif

					// line hits the end frame then quit
					if( oob.index == -1 )
					{
						break;
					}
					
					// should never see the same index as the last iteration go 
					// OOB because we just set up a specific path to make sure 
					// it wasn't OOB
					VT_ASSERT( oob.index != lastoob.index || fTestVelocity == 0 );

					// if new test did worse at staying in-bounds than previous
					// then we keep the previous
					if( oob.index <= lastoob.index )
					{						
						// use last V   
						oob = lastoob;

						// if best case still causes OOB during acceleration
						// then just punt on keeping the nearby values inbounds
						// TODO: need to clip these to inbounds
						if( i==0 &&
							oob.index+1 <= g_iMinFrameCountForPhase &&
							iSkip < g_iMinFrameCountForPhase )
						{
							iSkip++;
							goto restart;
						}

						break;
					}

					VT_ASSERT( fabs(oob.diff) > fMax );

					// compute a new test velocity that just misses going
					// out-of-bounds at the current OOB point
					fTestVelocity = 
						newV(init.fVelocity, acc_cnt, oob.index+1, 
							 fD0+oob.distB+0.999f*((oob.diff < 0)? fMax: -fMax));

					lastoob = oob;
				}

				if( i == 0 )
				{
					if( acc_cnt == g_iMinFrameCountForPhase || 
						OOB1_Better(oob, 0, oob0) )
					{
						oob0 = oob;
					}
				}
				else if( OOB1_Better(oob, i, oob0) )
				{
					// we are better off waiting at least a cycle
					goto done;
				}
			}

			fD0 += pfnUpdate(srcpath[frame+i], iW, iH) - init.fVelocity;
		}

#ifdef PATH_CSV
		/*testOOBPrint(init.fDistance, init.fVelocity, oob0.acc,
						 oob0.acc_cnt, iFrameCnt,
						 iOOBFrameV? 0: 1, srcpath, frame, fMax,
						 iW, iH, pfnUpdate);*/
		if( bPrintTests )
		{
			PrintTests(vecTests, init.fDistance, init.fVelocity, oob0.acc_cnt,
					   iFrameCnt, srcpath, frame, fMax, iW, iH, pfnUpdate);
		}
#endif
		result.iStartupCount = oob0.acc_cnt - 1;
		result.fAcceleration = oob0.acc;
	}

done:
	// update the current state
	result.fVelocity += result.fAcceleration;
	result.fDistance += pfnUpdate(srcpath[frame], iW, iH) - result.fVelocity;
	return result;
}

PATH_STATE ChangeScheduleTransX(const PATH_STATE& init,
                                const CRollingBuffer<FEAT_SIMILARITY>& srcpath,
                                int frame, int iFrameCnt, int iW, int iH, float fCrop)
{
    float fBorder = (1.f - fCrop) / 2.f;
    float fMax    = fBorder * g_fMaxDist;

	return ChangeScheduleInternal(init, srcpath, frame, iFrameCnt, iW, iH, fMax, 
                                  DistTransX);
}

PATH_STATE ChangeScheduleTransY(const PATH_STATE& init,
                                const CRollingBuffer<FEAT_SIMILARITY>& srcpath,
                                int frame, int iFrameCnt, int iW, int iH, float fCrop)
{
    float fBorder = (1.f - fCrop) / 2.f;
    float fMax    = fBorder * g_fMaxDist;

	return ChangeScheduleInternal(init, srcpath, frame, iFrameCnt, iW, iH, fMax,
                                  DistTransY, true);
}

PATH_STATE ChangeScheduleRot(const PATH_STATE& init,
                             const CRollingBuffer<FEAT_SIMILARITY>& srcpath,
                             int frame, int iFrameCnt, int iW, int iH, float fCrop)
{
	float fBorder =  float(VT_PI)/4.f - acos(float(VT_SQRT1_2)/fCrop);
	float fMax    =  VtMin(fBorder * 0.5f, float(VT_PI)/120.f); 

	return ChangeScheduleInternal(init, srcpath, frame, iFrameCnt, iW, iH, fMax,
                                  DistRotation);
}

PATH_STATE ChangeScheduleScale(const PATH_STATE& init,
							   const CRollingBuffer<FEAT_SIMILARITY>& srcpath,
                               int frame, int iFrameCnt, int iW, int iH, float fCrop)
{
	const float fBorder =  2.f * (0.975f - 0.5f);
	float fMax    =  fabs(log(fBorder) * g_fMaxDist); 

	return ChangeScheduleInternal(init, srcpath, frame, iFrameCnt, iW, iH, fMax,
                                  DistScale);
}

HRESULT CDPPathSmoother::DPFilterTransform(CMtx3x3f& result, int iDst)
{
    if( iDst == 0 )
    {  
        m_dpstateCur[0] = m_dpstateCur[1] = m_dpstateCur[2] = m_dpstateCur[3] = 
            PATH_STATE();
    }

    int iCount = m_bufTrackerSim.get_total_count()-iDst;

    // solve for each parameter independently
    m_dpstateCur[0] = ChangeScheduleTransX(
        m_dpstateCur[0], m_bufTrackerSim, iDst, iCount, 
        m_params.frameWidth, m_params.frameHeight, m_params.fCrop);

	m_dpstateCur[1] = ChangeScheduleTransY(
        m_dpstateCur[1], m_bufTrackerSim, iDst, iCount, 
        m_params.frameWidth, m_params.frameHeight, m_params.fCrop);

	m_dpstateCur[2] = ChangeScheduleRot(
		m_dpstateCur[2], m_bufTrackerSim, iDst, iCount, 
        m_params.frameWidth, m_params.frameHeight, m_params.fCrop);

    m_dpstateCur[3] = ChangeScheduleScale(
        m_dpstateCur[3], m_bufTrackerSim, iDst, iCount, 
        m_params.frameWidth, m_params.frameHeight, m_params.fCrop);

    FEAT_SIMILARITY sp(m_dpstateCur[0].fVelocity*float(m_params.frameWidth), 
                       m_dpstateCur[1].fVelocity*float(m_params.frameHeight),
				       exp(m_dpstateCur[3].fVelocity), 
				       m_dpstateCur[2].fVelocity);

#ifdef PATH_CSV
    if( g_fpDbg )
    {
		g_vecSmoothDbg.push_back(sp);
    }
#endif

	CMtx3x3f mtxSmooth = MatrixFromSimParams(sp);
    if( iDst == 0 )
    {
        result = mtxSmooth;
    }
    else
    { 
		CMtx3x3f mtxSimToPrev = MatrixFromSimParams(m_bufTrackerSim[iDst]);
		result = mtxSimToPrev.Inv() * mtxSmooth * m_bufResults[iDst-1];
    }

    return S_OK;
}

