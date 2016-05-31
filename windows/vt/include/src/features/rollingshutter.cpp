#include "stdafx.h"

#include "rollingshutter.h"

#ifndef VT_GCC

#if defined(MSRVT_WINDOWS_BUILD)
static const int   MOTION_SAMPLES = 16;
static const int   MAX_CG_ITERATIONS = 15;
static const int   DEFAULT_BUFFER = 5;
#else
static const int   MOTION_SAMPLES = 16;
static const int   MAX_CG_ITERATIONS = 25;
static const int   DEFAULT_BUFFER = 7;
#endif
static const float TEMPORAL_SMOOTHNESS = 10.0f;
static const float CG_THRESHOLD1 = 1.0f;
static const float CG_THRESHOLD2 = 0.001f;
static const float DEFAULT_ALPHA = 0.75f;

RSC::RSC()
{
    m_iInputWidth  = 0;
    m_iInputHeight = 0;
    m_iRSBufferSize = DEFAULT_BUFFER;
    m_fTemporalSmoothness = TEMPORAL_SMOOTHNESS;
    m_iMotionSamples = MOTION_SAMPLES;
    m_iMaxCGInterations = MAX_CG_ITERATIONS;
    m_fCGThreshold1 = CG_THRESHOLD1;
    m_fCGThreshold2 = CG_THRESHOLD2;
    m_fAlphaMultiplier = DEFAULT_ALPHA;
    m_bStarted = false;
}

BUFFER_RANGE RSC::GetResultsRange()
{
    BUFFER_RANGE r;
    int iFrameCnt = m_correspondences.get_total_count();
    if( m_bStarted )
    {
        r.first_frame = iFrameCnt-GetMaxDelay();
        r.frame_count = (r.first_frame>=0)? 1: 0;
    }
    else
    {  
        r.frame_count = VtMin(GetMaxDelay(), iFrameCnt);
        r.first_frame = iFrameCnt - r.frame_count;
    }
    return r;
}

HRESULT RSC::Begin(int iWidth, int iHeight, float fSmoothness)
{
    VT_HR_BEGIN();

    Deallocate();

    m_iInputWidth  = iWidth;
    m_iInputHeight = iHeight;
    m_fTemporalSmoothness = fSmoothness;

    VT_HR_EXIT(m_vecStartNonZero.resize(m_iMotionSamples*m_iRSBufferSize));
    VT_HR_EXIT(m_vecEndNonZero.resize(m_iMotionSamples*m_iRSBufferSize));

    VT_HR_EXIT(m_vecMotionX.resize(m_iMotionSamples));
    VT_HR_EXIT(m_vecMotionY.resize(m_iMotionSamples));
    VT_HR_EXIT(m_vecTranslationX.resize(m_iMotionSamples+1));
    VT_HR_EXIT(m_vecTranslationY.resize(m_iMotionSamples+1));
    VT_HR_EXIT(m_vecCorrection.resize(2*m_iInputHeight));
    VT_HR_EXIT(m_correspondences.resize(m_iRSBufferSize));

    m_bStarted = true;

    VT_HR_END();
}

HRESULT RSC::PushFrame(const vt::vector<PointMatch>& vecFP)
{
    VT_HR_BEGIN();
    
    // copy the correspondences to the internal rolling buffer
    m_correspondences.advance();
    vt::vector<PointMatch>& matchesBuffer = m_correspondences.last();
    VT_HR_EXIT( matchesBuffer.resize(vecFP.size()) );
    memcpy(matchesBuffer.begin(), vecFP.begin(), 
           vecFP.size()*sizeof(PointMatch));
     
    // Solve for the motion once we have a large enough buffer 
    if( m_correspondences.get_total_count() >= GetMaxDelay() )
    {
        VT_HR_EXIT(SolveForMotion()); 
    }

    VT_HR_END();
}

HRESULT RSC::GetResult(const float*& pfCorrection, 
                       const vt::vector<PointMatch>*& pMatches, int frameNumber)
{
    pfCorrection = NULL;

    VT_HR_BEGIN()

    if( !GetResultsRange().InRange(frameNumber) )
    {
        VT_HR_EXIT(E_INVALIDARG);
    }
 
    // return the matches
    // TODO: rolling shutter component shouldn't even store the matches
    pMatches = &m_correspondences[frameNumber];

    // Compute the Correction
    int iFrameCnt     = m_correspondences.get_total_count();
    int iNumberFrames = VtMin(m_iRSBufferSize, iFrameCnt);
    int iCenter = m_bStarted? VtMin(frameNumber, GetMaxDelay()-1): 
        iNumberFrames-(iFrameCnt-frameNumber);
    VT_HR_EXIT( ComputeCorrectionFromMotion(iCenter) );
    pfCorrection = m_vecCorrection.begin();

    VT_HR_END()
}

void RSC::Deallocate()
{
    m_iInputWidth  = 0;
    m_iInputHeight = 0;

    m_matA.Free();
    m_vecBX.Free();
    m_vecBY.Free();
    m_vecX.Free();
    m_vecY.Free();
    m_vecR.Free();
    m_vecP.Free();
    m_vecAP.Free();
    m_vecMotionX.clear();
    m_vecMotionY.clear();
    m_vecTranslationX.clear();
    m_vecTranslationY.clear();
    m_vecCorrection.clear();
    m_vecStartNonZero.clear();
    m_vecEndNonZero.clear();
    m_correspondences.clear();
}

HRESULT RSC::SolveForMotion()
{
    VT_HR_BEGIN();

    // Compute the Number of Unknowns
    int iNumberFrames = VtMin(m_iRSBufferSize, m_correspondences.get_total_count());
    int iLastFrameId  = m_correspondences.get_last_id();

    VT_ASSERT( iNumberFrames > 0 ); 

    int iUnknowns = iNumberFrames*m_iMotionSamples;

    // Allocate Space
    VT_HR_EXIT(m_matA.Create(iUnknowns, iUnknowns));
    VT_HR_EXIT(m_vecBX.Create(iUnknowns));
    VT_HR_EXIT(m_vecBY.Create(iUnknowns));
    VT_HR_EXIT(m_vecX.Create(iUnknowns));
    VT_HR_EXIT(m_vecY.Create(iUnknowns));
    VT_HR_EXIT(m_vecR.Create(iUnknowns));
    VT_HR_EXIT(m_vecP.Create(iUnknowns));
    VT_HR_EXIT(m_vecAP.Create(iUnknowns));

    // Set Up the Problem Step 1: Zero A and B
    float *fpA = m_matA.Ptr();
    float *fpBX = m_vecBX.Ptr();
    float *fpBY = m_vecBY.Ptr();
    float *fpX = m_vecX.Ptr();
    float *fpY = m_vecY.Ptr();
    for(int iU=0; iU<iUnknowns; iU++)
    {    
        for(int iU2=0; iU2<iUnknowns; iU2++)
        {
            *fpA++ = 0.0f;
        }
        *fpBX++ = 0.0f;
        *fpBY++ = 0.0f;
        *fpX++ = 0.0f;
        *fpY++ = 0.0f;
    }

    // Set Up the Problem Step 2: Data Term
    // Only need iNumberFrames-1 for example for frame 3 and a 7 frame buffer
    // we use the following feature matches:
    //     | 0->1 | 1->2 | 2->3 | 3->4 | 4->5 | 5->6 |
    int iNumberCorrespondences = 0;
    float fMotionSamples = float(m_iMotionSamples);
    float fAlpha = m_fAlphaMultiplier / float(m_iInputHeight);
    for(int iImage=0; iImage<iNumberFrames-1; iImage++)
    {
        float fImage = float(iImage);
        // correspondances are ordered newest to oldest
        const vt::vector<PointMatch>& vecFrameMatches = 
            m_correspondences[iLastFrameId-(iNumberFrames-2)+iImage];
        int iCorrespondenceCount = (int)vecFrameMatches.size();
        for(int iCount = 0; iCount<iCorrespondenceCount; iCount++)
        {
            iNumberCorrespondences++;
            const PointMatch& cCorrespondence = vecFrameMatches[iCount];
            float fY  = cCorrespondence.p1.y;
            float fDX = cCorrespondence.p0.x - cCorrespondence.p1.x;
            float fDY = cCorrespondence.p0.y - cCorrespondence.p1.y;

            // Compute the Start and End Times
            float fTime1 = fMotionSamples * (fImage + fAlpha * fY);
            float fTime1Floor = floor(fTime1);
            int iTime1Floor = int(fTime1Floor);
            float fTime2 = fMotionSamples * (fImage + 1.0f + fAlpha * (fY + fDY));
            float fTime2Ceil = ceil(fTime2);
            int iTime2Ceil = int(fTime2Ceil);

            // Add the Constraint
            // Following code not completely optimized -> just optimized inner loop
            // But it isn't too bad
            if (iTime1Floor >= 0 && iTime2Ceil < iNumberFrames * m_iMotionSamples)
            { 
                float fWeight = 1.0f - (fTime1 - fTime1Floor);
                for (int iIndex = iTime1Floor; iIndex <= iTime2Ceil; iIndex++)
                {
                    float fWeight2 = 1.0f - (fTime1 - fTime1Floor);
                    float *pfData = m_matA.Ptr() + iIndex*iUnknowns + iTime1Floor;
                    int iWidth = iIndex*iUnknowns + iTime2Ceil-1;
                    int iWidthTrunc = (iWidth/4)*4;
                    float *pfDataEnd = m_matA.Ptr() + iWidthTrunc;
                    float *pfDataEnd2 = m_matA.Ptr() + iWidth;
                    // m_matA.El(iIndex, iTime1Floor) += fWeight * fWeight2;
                    *pfData++ += fWeight * fWeight2;
                    // for (int iIndex2 = iTime1Floor+1; iIndex2 < iTime2Ceil-1; iIndex2++)
#if (defined(_M_IX86) || defined(_M_AMD64))
                    if ( g_SupportSSE2() )
                    {
                        __m128 m128Weight = _mm_set_ps1(fWeight);
                        while(pfData<pfDataEnd)
                        {
                            // m_matA.El(iIndex, iIndex2) += fWeight;
                            // *pfData++ += fWeight;
                            __m128 m128Data = _mm_loadu_ps(pfData);
                            m128Data = _mm_add_ps(m128Data, m128Weight);
                            _mm_storeu_ps(pfData, m128Data);
                            pfData+=4;
                        }
                    }
#elif defined(_M_ARM)
#endif
                    while(pfData<pfDataEnd2)
                    {
                        // m_matA.El(iIndex, iIndex2) += fWeight;
                        *pfData++ += fWeight;
                    }
                    fWeight2 = 1.0f - (fTime2Ceil - fTime2);
                    // m_matA.El(iIndex, iTime2Ceil-1) += fWeight * fWeight2;
                    *pfData += fWeight * fWeight2;

                    m_vecBX.El(iIndex) += fWeight * fDX;
                    m_vecBY.El(iIndex) += fWeight * fDY;

                    fWeight = 1.0;
                    if (iIndex == iTime2Ceil-1)
                    {
                        fWeight = 1.0f - (fTime2Ceil - fTime2);
                    }        
                }
            }
        }

/*
        for(int iY=m_iFlowBorder; iY<m_iFlowHeight-m_iFlowBorder; iY++)
        {
            iNumberCorrespondences++;

            // Get the Correspondence
            float fY = float(iY);
            int iCorrespondenceStart = 2*m_iFlowHeight*(iImage%m_iFlowBufferSize);
            float fDX = m_pfCorrespondenceBuffer[iCorrespondenceStart+2*iY];
            float fDY = m_pfCorrespondenceBuffer[iCorrespondenceStart+2*iY+1];

            // Compute the Start and End Times
            float fTime1 = fMotionSamples * (fImage + fAlpha * fY);
            float fTime1Floor = floor(fTime1);
            int iTime1Floor = int(fTime1Floor);
            float fTime2 = fMotionSamples * (fImage + 1.0f + fAlpha * (fY + fDY));
            float fTime2Ceil = ceil(fTime2);
            int iTime2Ceil = int(fTime2Ceil);

            // Add the Constraint
            // Following code not completely optimized -> just optimized inner loop
            // But it isn't too bad
            if (iTime1Floor >= 0 && iTime2Ceil < iNumberFrames * m_iMotionSamples)
            { 
                float fWeight = 1.0f - (fTime1 - fTime1Floor);
                for (int iIndex = iTime1Floor; iIndex <= iTime2Ceil; iIndex++)
                {
                    float fWeight2 = 1.0f - (fTime1 - fTime1Floor);
                    float *pfData = m_matA.Ptr() + iIndex*iUnknowns + iTime1Floor;
                    int iWidth = iIndex*iUnknowns + iTime2Ceil-1;
                    int iWidthTrunc = (iWidth/4)*4;
                    float *pfDataEnd = m_matA.Ptr() + iWidthTrunc;
                    float *pfDataEnd2 = m_matA.Ptr() + iWidth;
                    // m_matA.El(iIndex, iTime1Floor) += fWeight * fWeight2;
                    *pfData++ += fWeight * fWeight2;
                    // for (int iIndex2 = iTime1Floor+1; iIndex2 < iTime2Ceil-1; iIndex2++)
                    while(pfData<pfDataEnd)
                    {
                        // m_matA.El(iIndex, iIndex2) += fWeight;
                        *pfData++ += fWeight;
                    }
                    while(pfData<pfDataEnd2)
                    {
                        // m_matA.El(iIndex, iIndex2) += fWeight;
                        *pfData++ += fWeight;
                    }
                    fWeight2 = 1.0f - (fTime2Ceil - fTime2);
                    // m_matA.El(iIndex, iTime2Ceil-1) += fWeight * fWeight2;
                    *pfData += fWeight * fWeight2;

                    m_vecBX.El(iIndex) += fWeight * fDX;
                    m_vecBY.El(iIndex) += fWeight * fDY;

                    fWeight = 1.0;
                    if (iIndex == iTime2Ceil-1)
                    {
                        fWeight = 1.0f - (fTime2Ceil - fTime2);
                    }        
                }
            }
        }
        */
    }

    float fCorrPerImage = float(iNumberCorrespondences) / float(iNumberFrames);
    float fTemporalSmoothness = m_fTemporalSmoothness * fCorrPerImage;
    
    // Add the Regularization
    for(int iIndex = 0; iIndex < iNumberFrames*m_iMotionSamples - 1; iIndex++)
    {
        if (iIndex%m_iMotionSamples != m_iMotionSamples - 1)
        {
            int iU = iIndex;
            int iU2 = iIndex+1;
            float *fpRow1 = m_matA[iU];
            fpRow1[iU] += fTemporalSmoothness;
            fpRow1[iU2] -= fTemporalSmoothness;
            float *fpRow2 = m_matA[iU2];
            fpRow2[iU] -= fTemporalSmoothness;
            fpRow2[iU2] += fTemporalSmoothness;
        }
    }

    // Compute the Start and End of the Non-Zero Elements
    for(int iUnknown=0; iUnknown<iUnknowns; iUnknown++)
    {
        float *pfMatrixPointer = m_matA.Ptr() + iUnknown*iUnknowns;
        float *pfMatrixPointerEnd = m_matA.Ptr() + iUnknown*iUnknowns + iUnknowns;
        int *piStart = m_vecStartNonZero.begin()+iUnknown;
        *piStart = 0;
        while(*pfMatrixPointer++ == 0.0f && pfMatrixPointer<pfMatrixPointerEnd)
        {
            (*piStart)++;
        }
        pfMatrixPointer = m_matA.Ptr() + iUnknown*iUnknowns + iUnknowns;
        pfMatrixPointerEnd = m_matA.Ptr() + iUnknown*iUnknowns;
        int *piEnd = m_vecEndNonZero.begin()+iUnknown;
        *piEnd = iUnknowns;
        while(*(--pfMatrixPointer) == 0.0f && pfMatrixPointer>pfMatrixPointerEnd)
        {
            (*piEnd)--;
        }
    }

    // Run Conjugate Gradient
    ConjugateGradient(m_vecX, m_vecBX);
    ConjugateGradient(m_vecY, m_vecBY);

    VT_HR_END();
}

void RSC::ConjugateGradient(CVec<float> &vecX, CVec<float> &vecB)
{
    m_vecR = vecB; // - m_matA * vecX; 
    m_vecP = m_vecR;
    float fRROld = m_vecR * m_vecR;
    float fRRNew = fRROld;
    int iUnknowns = m_vecR.Size();
    for(int iIteration=0; iIteration<m_iMaxCGInterations; iIteration++)
    {
        // m_vecAP = m_matA * m_vecP;
        for(int iUnknown = 0; iUnknown<iUnknowns; iUnknown++)
        {
            float *pfResultPointer = m_vecAP.Ptr() + iUnknown;
            float *pfMatrixPointer = m_matA.Ptr() + iUnknowns*iUnknown + m_vecStartNonZero[iUnknown];
            float *pfSourcePointer = m_vecP.Ptr() + m_vecStartNonZero[iUnknown];
            float *pfSourcePointerEnd = m_vecP.Ptr() + m_vecEndNonZero[iUnknown];
            *pfResultPointer = 0.0f;
#if (defined(_M_IX86) || defined(_M_AMD64))
            if ( g_SupportSSE2() )
            {
                int iWidth = m_vecEndNonZero[iUnknown] - m_vecStartNonZero[iUnknown];
                int iWidthTrunc = 4 * (iWidth / 4);
                float *pfSourcePointerEnd2 = pfSourcePointer + iWidthTrunc;
                float fResult = 0.0f;
                while(pfSourcePointer<pfSourcePointerEnd2)
                {
                    // *pfResultPointer += (*pfMatrixPointer++) * (*pfSourcePointer++);
                    __m128 m128Matrix = _mm_loadu_ps(pfMatrixPointer);
                    pfMatrixPointer+=4;
                    __m128 m128Source = _mm_loadu_ps(pfSourcePointer);
                    pfSourcePointer+=4;
                    __m128 m128DP = _mm_mul_ps(m128Matrix, m128Source);
                    fResult +=  m128DP.m128_f32[0] + m128DP.m128_f32[1] +
                        m128DP.m128_f32[2] + m128DP.m128_f32[3];
                }
                *pfResultPointer = fResult;
            }
#elif defined(_M_ARM)
#endif
            while(pfSourcePointer<pfSourcePointerEnd)
            {
                *pfResultPointer += (*pfMatrixPointer++) * (*pfSourcePointer++);
            }
        }
        float fPAP = m_vecP * m_vecAP;
        if (fPAP < m_fCGThreshold1)
        {
            break;
        }
        float fAlpha = fRRNew / fPAP;
        vecX += fAlpha * m_vecP;
        m_vecR -= fAlpha * m_vecAP;
        fRROld = fRRNew;
        fRRNew = m_vecR * m_vecR;
        if (fRRNew < m_fCGThreshold2)
        {
            break;
        }
        float fBeta = fRRNew / fRROld;
        m_vecP *= fBeta;
        m_vecP += m_vecR;
    }
}

HRESULT RSC::ComputeCorrectionFromMotion(int bufidx)
{
    VT_HR_BEGIN();

    // Extract the Motion
    int iCorrectionStart = m_iMotionSamples*bufidx;
    for(int iSample=0; iSample<m_iMotionSamples; iSample++)
    {
        m_vecMotionX[iSample] = m_vecX.El(iCorrectionStart+iSample);
        m_vecMotionY[iSample] = m_vecY.El(iCorrectionStart+iSample);
    }

    // compute the correction
    float *pfCorrection = m_vecCorrection.begin();

    // Compute Render Time
    float fRenderTime = 0.5f*m_fAlphaMultiplier*float(m_iMotionSamples);
    int iTimeBeforeRT = int(floor(fRenderTime));
    int iTimeAfterRT = int(ceil(fRenderTime));
    VT_HR_EXIT((iTimeBeforeRT < 0 || iTimeAfterRT >= m_iMotionSamples) ? E_FAIL : S_OK);

    // Compute Alpha
    float fAlphaInvDivMotionSamples = float(m_iInputHeight) / m_fAlphaMultiplier / float(m_iMotionSamples);

    // Integrate Motion to Give Translation - Using Linear Interpolation
    m_vecTranslationX[iTimeBeforeRT] = -(fRenderTime-float(iTimeBeforeRT))*m_vecMotionX[iTimeBeforeRT];
    m_vecTranslationX[iTimeAfterRT]  = (float(iTimeAfterRT)-fRenderTime)*m_vecMotionX[iTimeBeforeRT];
    m_vecTranslationY[iTimeBeforeRT] = -(fRenderTime-float(iTimeBeforeRT))*m_vecMotionY[iTimeBeforeRT];
    m_vecTranslationY[iTimeAfterRT]  = (float(iTimeAfterRT)-fRenderTime)*m_vecMotionY[iTimeBeforeRT];
    for(int iTime=iTimeAfterRT; iTime<m_iMotionSamples; iTime++)
    {
        m_vecTranslationX[iTime+1] = m_vecTranslationX[iTime] + m_vecMotionX[iTime];
        m_vecTranslationY[iTime+1] = m_vecTranslationY[iTime] + m_vecMotionY[iTime];
    }
    for(int iTime=iTimeBeforeRT-1; iTime>=0; iTime--)
    {
        m_vecTranslationX[iTime] = m_vecTranslationX[iTime+1] - m_vecMotionX[iTime];
        m_vecTranslationY[iTime] = m_vecTranslationY[iTime+1] - m_vecMotionY[iTime];
    }

    float fx = m_vecTranslationX[0];
    float fy = m_vecTranslationY[0];

    // Loop Over the Row in the Image
    for(int iY=0; iY<m_iInputHeight; iY++)
    {
        // Loop Over the Motion Samples
        int iTime = 0;
        bool bFoundSolution = false;
        while(bFoundSolution == false && iTime<m_iMotionSamples)
        {
            // Check for Intersection
            float fStartDiff = fAlphaInvDivMotionSamples * float(iTime) - m_vecTranslationY[iTime] - iY;
            float fEndDiff =  fAlphaInvDivMotionSamples * float(iTime+1) - m_vecTranslationY[iTime+1] - iY;
            if (fStartDiff * fEndDiff <= 0.0)
            {
                if (fStartDiff - fEndDiff != 0)
                {
                    float fSolutionTime = fStartDiff / (fStartDiff - fEndDiff);
                    if (fSolutionTime >= 0 && fSolutionTime <= 1.0f)
                    {
                        pfCorrection[2*iY] = fSolutionTime * m_vecTranslationX[iTime+1] + 
                            (1.0f - fSolutionTime) * m_vecTranslationX[iTime];
                        pfCorrection[2*iY+1] = fSolutionTime * m_vecTranslationY[iTime+1] + 
                            (1.0f - fSolutionTime) * m_vecTranslationY[iTime];
                        bFoundSolution = true;
                    }
                }
            }
            iTime++;
        }
        if (bFoundSolution == false)
        {
            pfCorrection[2*iY] = -99999.9f;
            pfCorrection[2*iY+1] = -99999.9f;
        }
    }

    // Attempt to Fill In Missing Values
    float fLastX = -99999.9f;
    float fLastY = -99999.9f;
    for(int iY=m_iInputHeight/2; iY>=0; iY--)
    {
        if (pfCorrection[2*iY] == -99999.9f)
        {
            pfCorrection[2*iY] = fLastX;
            pfCorrection[2*iY+1] = fLastY;
        }
        else
        {
            fLastX = pfCorrection[2*iY];
            fLastY = pfCorrection[2*iY+1];
        }
    }
    fLastX = -99999.9f;
    fLastY = -99999.9f;
    for(int iY=m_iInputHeight/2; iY<m_iInputHeight; iY++)
    {
        if (pfCorrection[2*iY] == -99999.9f)
        {
            pfCorrection[2*iY] = fLastX;
            pfCorrection[2*iY+1] = fLastY;
        }
        else
        {
            fLastX = pfCorrection[2*iY];
            fLastY = pfCorrection[2*iY+1];
        }
    }

    VT_HR_END();
}

#endif // VT_GCC

//+-----------------------------------------------------------------------------
//
// Class: CRollingShutterAddressGen
//
//------------------------------------------------------------------------------
#if (defined(_M_IX86) || defined(_M_AMD64))
#define VT_RSAG_SIMD_SUPPORTED (g_SupportSSSE3())
#elif defined(_M_ARM)
#define VT_RSAG_SIMD_SUPPORTED true
#endif

vt::CRect CRollingShutterAddressGen::MapDstRectToSrc(IN const CRect& rctDstPar)
{
    CRect rctDst = rctDstPar;

    // apply matrix rect first, if present
    if (m_bApplyMatrix)
    {
        rctDst = MapRegion3x3(m_xfrm, rctDst, NULL);
    }

    // TODO: when we move to a arb. length LU table change all the clamps to
    //       clamp to that
    int iTS = int(VtClamp(m_fScaleInv*float(rctDst.top), 
                          0.f, m_fScaleInv*float(m_iDstH)-1.f));
    int iBS = int(VtClamp(m_fScaleInv*float(rctDst.bottom), 
                          0.f, m_fScaleInv*float(m_iDstH)));

    
    // compute min/max X offsets for this rect    
    float minx = FLT_MAX, maxx = -FLT_MAX;
    for( int y = iTS; y < iBS; y++ )
    {
        const CVec2f& adj = m_pAdj[y];
        minx = VtMin(adj.x, minx);
        maxx = VtMax(adj.x, maxx);
    }

    vt::CRect rctSrc;
    // need floor prior to cast since the int cast truncate towards zero
    // and the address computations use the floor
    rctSrc.left   = int(floorf(float(rctDst.left)   + minx*m_fScale));
    rctSrc.right  = int(floorf(float(rctDst.right)  + maxx*m_fScale))+1;
    rctSrc.top    = int(floorf(float(rctDst.top)    + m_pAdj[iTS].y*m_fScale));
    rctSrc.bottom = int(floorf(float(rctDst.bottom) + m_pAdj[iBS-1].y*m_fScale))+1;

    return rctSrc;
}

HRESULT CRollingShutterAddressGen::MapDstSpanToSrc(
    OUT CVec2f* pOut, IN const vt::CPoint &ptDst, int iSpan)
{
    VT_ASSERT( ptDst.x >= 0 && ptDst.x < m_iDstW );
    VT_ASSERT( ptDst.y >= 0 && ptDst.y < m_iDstH );

    if (m_bApplyMatrix && !m_bAffineMatrix)
    {
        // the Y for RS adjustment lookup is not linear, so can't use the
        // same method as is used for affine; instead do two passes with the
        // first doing SIMD generation of non-RS adjusted addresses (just like
        // C3x3TransformAddressGen) and then in the second pass add the RS
        // adjustment factors, avoiding the lookup when possible via index checks

        float X = float(ptDst.x) * m_xfrm(0,0) + float(ptDst.y) * m_xfrm(0,1) + m_xfrm(0,2);
        float Y = float(ptDst.x) * m_xfrm(1,0) + float(ptDst.y) * m_xfrm(1,1) + m_xfrm(1,2);
        float Z = float(ptDst.x) * m_xfrm(2,0) + float(ptDst.y) * m_xfrm(2,1) + m_xfrm(2,2);
        float Ymax = (m_fScaleInv*float(m_iDstH))-1.f;

        int x = 0;
#if (defined(_M_IX86) || defined(_M_AMD64) || defined(_M_ARM))
        float Zend = Z + (m_xfrm(2,0)*float(iSpan-1));
        if( (VT_RSAG_SIMD_SUPPORTED) && (Z>0.f) && (Zend>0.f) )
        {
            // vector regs for iterating two xyz's
            __m128 x0y0x1y1 = _mm_setr_ps(X,Y,X+m_xfrm(0,0),Y+m_xfrm(1,0));
            __m128 dxy = _mm_setr_ps(2*m_xfrm(0,0),2*m_xfrm(1,0),2*m_xfrm(0,0),2*m_xfrm(1,0));
            __m128 z0z0z1z1 = _mm_setr_ps(Z,Z,Z+m_xfrm(2,0),Z+m_xfrm(2,0));
            __m128 dz = _mm_set1_ps(2*m_xfrm(2,0));
            __m128 adjScale = _mm_set1_ps(m_fScale);

            if (IsAligned16(pOut))
            {
                for( ; x < iSpan - 1; x += 2, pOut += 2)
                {
#if 1
                    __m128 xyOut = _simd_recip_lp_ps(z0z0z1z1);
                    xyOut = _mm_mul_ps(xyOut,x0y0x1y1);
#else
                    __m128 xyOut = _mm_div_ps(x0y0x1y1,z0z0z1z1);
#endif
                    _mm_store_ps((float*)pOut,xyOut);
                    x0y0x1y1 = _mm_add_ps(x0y0x1y1,dxy);
                    z0z0z1z1 = _mm_add_ps(z0z0z1z1,dz);
                }
            }
            else
            {
                for( ; x < iSpan - 1; x += 2, pOut += 2)
                {
#if 1
                    __m128 xyOut = _simd_recip_lp_ps(z0z0z1z1);
                    xyOut = _mm_mul_ps(xyOut,x0y0x1y1);
#else
                    __m128 xyOut = _mm_div_ps(x0y0x1y1,z0z0z1z1);
#endif
                    _mm_storeu_ps((float*)pOut,xyOut);
                    x0y0x1y1 = _mm_add_ps(x0y0x1y1,dxy);
                    z0z0z1z1 = _mm_add_ps(z0z0z1z1,dz);
                }
            }
            // update scalars for remainder
            X = SSE2_mm_extract_epf32(x0y0x1y1, 0);
            Y = SSE2_mm_extract_epf32(x0y0x1y1, 1);
            Z = SSE2_mm_extract_epf32(x0y0x1y1, 0);
        }
#endif
        // remainder or no SIMD support or zero or negative Z's
        for( ; x < iSpan; x ++, pOut++)
        {
            // check for negative Z and set OOB address if so
            if (Z <= 0.f)
            {
                SetAddressOOB(*pOut);
            }
            else
            {
                float recipZ = 1.f/Z;
                pOut->x = X*recipZ;
                pOut->y = Y*recipZ;
            }
            X += m_xfrm(0,0); Y += m_xfrm(1,0); Z += m_xfrm(2,0);
        }

        // back up output pointer and apply RS adjustment
        pOut -= iSpan;
        x = 0;
        int Yadj = 0;
        CVec2f Adj = m_pAdj[Yadj]*m_fScale;
        for( ; x < iSpan; x ++, pOut++)
        {
            int Yadjt = (int)(pOut->y);
            if (Yadjt != Yadj)
            {
                // load RS adjustment
                int idx = int(VtClamp(m_fScaleInv*float(Yadjt), 0.f, Ymax));    
                Adj = m_pAdj[idx]*m_fScale;
                Yadj = Yadjt;
            }
            pOut->x += Adj.x;
            pOut->y += Adj.y;
        }
#if 0 // check results against 'reference implementation'
        X = float(ptDst.x) * m_xfrm(0,0) + float(ptDst.y) * m_xfrm(0,1) + m_xfrm(0,2);
        Y = float(ptDst.x) * m_xfrm(1,0) + float(ptDst.y) * m_xfrm(1,1) + m_xfrm(1,2);
        Z = float(ptDst.x) * m_xfrm(2,0) + float(ptDst.y) * m_xfrm(2,1) + m_xfrm(2,2);
        pOut -= iSpan;
        float tolerance = 1.f/8.f;
#define CMPVAL(_val) (fabs(_val) > tolerance)
        for( x = 0; x < iSpan; x ++, pOut++)
        {
            if (Z <= 0.f)
            {
            }
            else
            {
                float recipZ = 1.f/Z;
                float Xp = X*recipZ; 
                float Yp = Y*recipZ;
                int iYp = (int)Yp;
                int idx = int(VtClamp(m_fScaleInv*float(iYp), 0.f, Ymax));    
                CVec2f Adj = m_pAdj[idx]*m_fScale;
                if ( CMPVAL(pOut->x -(Xp + Adj.x)) || CMPVAL(pOut->y - (Yp + Adj.y)) )
                {
                    VT_ASSERT(false);
                }
            }
            X += m_xfrm(0,0); Y += m_xfrm(1,0); Z += m_xfrm(2,0);
        }
#endif
    }
    else if (m_bApplyMatrix)
    {
        // compute X&Y for this row
        float X = float(ptDst.x) * m_xfrm(0,0) + float(ptDst.y) * m_xfrm(0,1) + 
            m_xfrm(0,2);
        float Y = float(ptDst.x) * m_xfrm(1,0) + float(ptDst.y) * m_xfrm(1,1) + 
            m_xfrm(1,2);
        float Ymax = (m_fScaleInv*float(m_iDstH))-1.f;

        // Y is linear and generally changing slowly, thus the same entry in the
        // RS adjustment table is reused for many steps in X; compute the number
        // of X steps between integral changes in Y (first and subsequent) and
        // move the RS adjustment lookup outside of the inner loops
        int xs = iSpan;             // end of current interval
        int iSubSpanSize = iSpan;   // steps in each interval
        {
            if ( fabs(m_xfrm(1,0)) > .000001f)
            {
                bool stepPos = m_xfrm(1,0) > 0.f;
                iSubSpanSize = (int)(fabs(1.f/m_xfrm(1,0))+.5f);
                float Yfrac = (Y - floorf(Y));
                if (stepPos) { Yfrac = 1.f-Yfrac; }
                xs = (int)((Yfrac*(float)iSubSpanSize)+.5f);
                xs = VtClamp(xs,0,iSpan);
            }
        }

#if (defined(_M_IX86) || defined(_M_AMD64) || defined(_M_ARM))
        __m128 xmm2 = _mm_setr_ps(2.f*m_xfrm(0,0), 2.f*m_xfrm(1,0),
                                  2.f*m_xfrm(0,0), 2.f*m_xfrm(1,0));
#endif        
        int x = 0;
        while ( x < iSpan )
        {
            // load RS adjustment
            int iYAdj = int(VtClamp(m_fScaleInv*Y, 0.f, Ymax));    
            CVec2f Adj = m_pAdj[iYAdj]*m_fScale;

#if (defined(_M_IX86) || defined(_M_AMD64) || defined(_M_ARM))
            if( (VT_RSAG_SIMD_SUPPORTED) && ((xs-x)>=16) )
            {
                // xmm0 = x0 y0 x1 y1
                __m128 xmm0 = _mm_setr_ps(X, Y, X+m_xfrm(0,0), Y+m_xfrm(1,0));

                // rolling shutter adjustment
                __m128 xmm1 = _mm_setr_ps(Adj.x, Adj.y, Adj.x, Adj.y);;
        
                for( ; x<(xs-1) ; x+=2, pOut+=2 )
                {
                    __m128 xmm3 = _mm_add_ps(xmm0, xmm1);
                    _mm_storeu_ps((float*)pOut, xmm3);  // store the result                
                    xmm0 = _mm_add_ps(xmm0, xmm2);      // xmm0 = xmm0 + xmm2
                }
        
                // need to re-init these for vector remainder
                X = SSE2_mm_extract_epf32(xmm0, 0);
                Y = SSE2_mm_extract_epf32(xmm0, 1);
            }
#endif
            for( ;x<xs ; x++, pOut++, X+=m_xfrm(0,0), Y+=m_xfrm(1,0) )
            {
                pOut->x = X + Adj.x;
                pOut->y = Y + Adj.y;
            }
            xs = VtMin(iSpan,x+iSubSpanSize);
            VT_ASSERT(xs>=0);
        }
    }
    else
    {
        int iY = int(VtClamp(m_fScaleInv * ptDst.y, 0.f, m_fScaleInv * m_iDstH - 1.f));
        CVec2f adj = CVec2f((float)ptDst.x, (float)ptDst.y)+ m_pAdj[iY];
        for( int i = 0; i < iSpan; i++ )
        {
            if (m_bApplyMatrix) { pOut[i] += adj; }
            else { pOut[i] = adj; }
            adj.x += 1.f; 
        }
    }
    return S_OK;
}

HRESULT CRollingShutterAddressGen::MapDstAddrToSrc(
    IN OUT CVec2f* pOut, int iSpan)
{
    if (m_bApplyMatrix)
    {
        if (m_bAffineMatrix)
        {
            for( int i = 0; i < iSpan; i++ )
            {
                CVec2f d = *pOut;
                pOut->x = m_xfrm(0,0) * d.x + m_xfrm(0,1) * d.y + m_xfrm(0,2);
                pOut->y = m_xfrm(1,0) * d.x + m_xfrm(1,1) * d.y + m_xfrm(1,2);
            }
        }
        else
        {
            for( int i = 0; i < iSpan; i++ )
            {
                CVec2f d = *pOut;
                pOut->x = m_xfrm(0,0) * d.x + m_xfrm(0,1) * d.y + m_xfrm(0,2);
                pOut->y = m_xfrm(1,0) * d.x + m_xfrm(1,1) * d.y + m_xfrm(1,2);
                float z = m_xfrm(2,0) * d.x + m_xfrm(2,1) * d.y + m_xfrm(2,2);
                if (z > 0.f)
                {
                    pOut->x /= z; pOut->y /= z;
                }
                else
                {
                    SetAddressOOB(*pOut);
                }
            }
        }
    }
    for( int i = 0; i < iSpan; i++ )
    {
        int iY = int(VtClamp(m_fScaleInv* pOut[i].y, 0.f, m_fScaleInv * m_iDstH - 1.f));    
        pOut[i] += m_pAdj[iY];
    }
    return S_OK;
}
