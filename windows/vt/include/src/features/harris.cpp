//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      specialized implementation of Harris feature detector
//
//------------------------------------------------------------------------------
#include "stdafx.h"

using namespace vt;

#include "vt_harrisdetect.h"
#include "vt_features_base.h"

// files containing inline functions for SSE and Neon processing
#if (defined(_M_IX86) || defined(_M_AMD64))
#include "harris_sse.h"
#elif defined(_M_ARM)
#include <arm_neon.h>
#include "harris_neon.h"
#endif

// scale unconditionally assigned to returned feature points
#define SCALE_CALIBRATION_HARRIS		1.35f

// max radius allowed by feature position refinement - detected points for
// which position refinement would be beyond this value are culled
#define HARRIS_REFINE_RADIUS            1.0f

#define HARRIS_DENOM_EPSILON			0.00001f

//------------------------------------------------------------------------------
//
// utility routines
//
//------------------------------------------------------------------------------
// compute one scanline of differences for Harris function; input is 8bpp unsigned
// and output is 16bpp signed (although always positive)
static inline void DoDrvLine(int w, 
    int16_t* res, size_t resStride,
    const uint8_t* srcP, const uint8_t* srcC, const uint8_t* srcN)
{
    srcC++; // advance so it points to x+1'th input
#if (defined(_M_IX86) || defined(_M_AMD64))
    HarrisDoDrvLineIntSSE(w,res,resStride,srcP,srcC,srcN);
#elif defined(_M_ARM) 
    HarrisDoDrvLineIntNeon(w,res,resStride,srcP,srcC,srcN);
#endif
    while (w > 0)
    {
        // subract .8 from .8 result is s.8
        int16_t dx = (*(srcC-2)) - (*(srcC++));
        int16_t dy = (*srcN++) - (*srcP++);
        // .8*.8 = .16; shift by one to get s.15 values (although they are always positive)
        int32_t xx = dx*dx; xx >>= 1;
        int32_t yy = dy*dy; yy >>= 1;
        int32_t xy = (dx>>1)*dy; // slightly less accurate, but matches SSE&Neon code
        *(res+(0*resStride)) = (int16_t)xx;
        *(res+(1*resStride)) = (int16_t)yy;
        *(res+(2*resStride)) = (int16_t)xy;
        res++;
        w--;
    }
}
//------------------------------------------------------------------------------
// operates on each Harris maxima to refine the feature position plus compute
// culling criteria for features based on the two eigenvalues
//
// copied from Sunflower then modified to do threshold culling and conditional
// engenvalue ratio culling
//
static void VtQuadraticFit2D(const float *pFm1, const float *pF0, const float *pFp1,
	CMtx2x2<float> &mHessRtn, CVec2<float> &vDerivRtn)
{
	float fRm1 = pFm1[-1] + pFm1[0] + pFm1[1];
	float fRp1 = pFp1[-1] + pFp1[0] + pFp1[1];
	float fCm1 = pFm1[-1] + pF0[-1] + pFp1[-1];
	float fCp1 = pFm1[1] + pF0[1] + pFp1[1];
	mHessRtn(0,0) = (fCm1 - 2*(pFm1[0] + pF0[0] + pFp1[0]) + fCp1)/(float)3.0;
	mHessRtn(1,1) = (fRm1 - 2*(pF0[-1] + pF0[0] + pF0[1]) + fRp1)/(float)3.0;
	mHessRtn(0,1) = (float)0.25 * (pFm1[-1] - pFm1[1] - pFp1[-1] + pFp1[1]);
	mHessRtn(1,0) = mHessRtn(0,1);
	vDerivRtn(0) = (-fCm1 + fCp1)/(float)6.0;
	vDerivRtn(1) = (-fRm1 + fRp1)/(float)6.0;
}

void RefineAndStoreResult(HRESULT& hr,
    vt::vector<HARRIS_FEATURE_POINT>& pts, int x, int y,
    const float *p0, const float *p1, const float *p2,
    const HARRIS_DETECTOR_PARAMS& prm)
{
    // improve position to sub-pixel accuracy
    CMtx2x2f mH;
    CVec2f vD;
    VtQuadraticFit2D(p0, p1, p2, mH, vD);
    float fDet = mH.Det();
    if (fDet<=0) { return; }

    float fAC = mH(0,0) + mH(1,1);
    float fD = sqrt(fabs(fAC*fAC - 4*fDet));
    float fSL1 = 0.5f * (fAC-fD);
    float fSL2 = 0.5f * (fAC+fD);
    if ((fSL1>=0) || (fSL2>=0)) { return; }

    float fL2 = fabs(fSL2); // smaller eigenvalue
    float fL1 = fabs(fSL1); // larger eigenvalue
    if (fL2 < prm.scoreThreshold) { return; }
    if (fL1 >= (prm.cornerEdgeThreshold * fL2)) { return; }

    CVec2f vdx;
    vdx.x = (-mH(1,1)*vD(0) + mH(0,1)*vD(1)) / fDet;
    vdx.y = (mH(0,1)*vD(0) - mH(0,0)*vD(1)) / fDet;
    float fdist = vdx.MagnitudeSq();
    if (fdist > (HARRIS_REFINE_RADIUS * HARRIS_REFINE_RADIUS)) { return; }

    HARRIS_FEATURE_POINT ip;
    ip.x = (float)x + vdx.x;
    ip.y = (float)y + vdx.y;
    ip.score = fL2; // smaller eigenvalue (rec'd by Shi&Thomasi)
  //  ip.score = *(p1+1); // just return the harris function value

    // add to interest point vector
    hr = pts.push_back(ip);
    return;
}

//------------------------------------------------------------------------------
//
// Harris detector operating on byte image data
//
// srcRect: rectangle over which features are detected; this may require that a
// larger area of the src image be present depending on the border: the Harris
// function is evaluated over a 3x3 area so expands srcRect-border by one; the
// filtering of the differences for the Harris function requires 3 pixels for
// each Harris function pixel; the computation of the differences prior to
// filtering is over a 3x3 area; thus the total area accessed outside
// srcRect-border is 5 (1 + 3 + 1)
//
//
//   R - result pixel (feature detect candidate)
//   H - Harris function pixel
//   f - pixels for filtering Harris function differences
//   d - pixels for computing differences prior to filtering
//   
//     | 0 1 2 3 4 5
//   --|------------
//   0 | - d d d d d  
//   1 | d f f f f f
//   2 | d f f f f f
//   3 | d f f f f f
//   4 | d f f f H H
//   5 | d f f f H R
//   
//------------------------------------------------------------------------------
HRESULT HarrisDetectBlock(vt::vector<HARRIS_FEATURE_POINT>& pts, 
    const CImg& srcImg, const CRect& srcRect, const CPoint& ptSrc,
    int iSrcWidth, int iSrcHeight,
    const HARRIS_DETECTOR_PARAMS& prm, 
    bool useHFBlockCache, int ftrPerBlk)
{
    VT_HR_BEGIN()

    // macro to access into source image
#define SRCPTR(_x, _y) srcImg.BytePtr((_x)-ptSrc.x,(_y)-ptSrc.y)

    // this code requires that the border be large enough to ensure that
    // the 3x3 Harris function computation area and the 7x7 filter area
    // never sample pixels outside the image extent - note that the border
    // also speeds up descriptor generation so is generally at least 6 anyway
    VT_HR_EXIT( (prm.border < 5) ? E_INVALIDARG : S_OK );

    // kernel width (from center)
    int kw = 3;

    // compute x,y start/end points for covered area incorporating prm.border;
    int xs = (srcRect.left < prm.border)?(prm.border):(srcRect.left);
    int xlast = iSrcWidth-prm.border;
    int xe =  VtMin((int)srcRect.right,VtMax(xlast,(int)srcRect.right-prm.border));
    int ys = (srcRect.top < prm.border)?(prm.border):(srcRect.top);
    int ylast = iSrcHeight-prm.border;
    int ye = VtMin((int)srcRect.bottom,VtMax(ylast,(int)srcRect.bottom-prm.border));

    // tile size for filtering result (harris function buffer) is rect expanded by 1
    // (since Harris function is sampled on a 3x3 grid for final result)
    int whf = (xe-xs+2);
    int hhf = (ye-ys+2);
    
    // tile width for filtering adds kernel width on both sides
    int wfilt = whf+2*kw;

    // create temporary buffers for derivatives and harris function
    CIntImg imgDrvH; VT_HR_EXIT( imgDrvH.Create(wfilt,3,1) );
    CLumaFloatImg imgH; 
    if (useHFBlockCache) { VT_HR_EXIT( imgH.Create(whf,hhf) ); } // full hf surface
    else                 { VT_HR_EXIT( imgH.Create(whf,3) ); }   // just 3 scanlines
    float* imgH0 = imgH.Ptr(1);
    float* imgH1 = imgH.Ptr(2);
    float* imgH2 = imgH.Ptr(0);

    // buffer for detector hits when using a full hf surface
    typedef struct { int x,y; float f; } detHit;
    vt::vector<detHit> vecDetHits;
    if (useHFBlockCache) { VT_HR_EXIT( vecDetHits.reserve(128) ); }

    // harris function detect threshold - will be initialized and then refined below
    float fDetectThreshold = 0.f;
    float fBlockAvg = 0.f;

    int ftrCnt = 0;

    // generate Harris function image; includes generation of the
    // squared differences, filtering, and computing the harris function;
    // difference and filtering operations are done with integer math;
    // runs a 7x7 filter with a 1 4 8 10 8 4 1 kernel
    //
    // uses a 7 scanline high window that steps down through the source
    // image; individual pointers are retained for each scanline in a rolling
    // buffer to enable reuse; the last row of the temp buffer is used for the
    // result of the vertical filtering and source for the horizontal filtering
    //
    // the filtered values are accumulated into 32 bit integers; the
    // normalization of the kernel weights is done after all of the contributions
    // have been summed and the value converted to floating point
    //
    // temp storage: ~tilewidth*7*3*2 bytes + width*3*4 bytes
    //
    // filtering data are non-interleaved because the computations for
    // each 'channel' are different and thus the SSE result is non-interleaved;
    // TODO: try interleaving the result during DoDrvLine and see if it
    //  makes the vertical filtering code enough faster to offset
    CLumaShortImg imgDrvV; VT_HR_EXIT( imgDrvV.Create(wfilt,21) );
    size_t srcStride = imgDrvV.StrideBytes()/sizeof(int16_t);
    CRollingBuffer<int16_t*> rbline;
    VT_HR_EXIT(rbline.resize(7));
    rbline.buffer(0) = (int16_t*)imgDrvV.Ptr( 0);
    rbline.buffer(1) = (int16_t*)imgDrvV.Ptr( 3);
    rbline.buffer(2) = (int16_t*)imgDrvV.Ptr( 6);
    rbline.buffer(3) = (int16_t*)imgDrvV.Ptr( 9);
    rbline.buffer(4) = (int16_t*)imgDrvV.Ptr(12);
    rbline.buffer(5) = (int16_t*)imgDrvV.Ptr(15);
    rbline.buffer(6) = (int16_t*)imgDrvV.Ptr(18);

    int y=(ys-1);
    int xsrc = (xs-1)-kw;
    DoDrvLine(wfilt,rbline.buffer(0),srcStride,SRCPTR(xsrc,y-3-1),SRCPTR(xsrc,y-3  ),SRCPTR(xsrc,y-3+1)); 
    DoDrvLine(wfilt,rbline.buffer(1),srcStride,SRCPTR(xsrc,y-2-1),SRCPTR(xsrc,y-2  ),SRCPTR(xsrc,y-2+1)); 
    DoDrvLine(wfilt,rbline.buffer(2),srcStride,SRCPTR(xsrc,y-1-1),SRCPTR(xsrc,y-1  ),SRCPTR(xsrc,y-1+1)); 
    DoDrvLine(wfilt,rbline.buffer(3),srcStride,SRCPTR(xsrc,y  -1),SRCPTR(xsrc,y    ),SRCPTR(xsrc,y  +1)); 
    DoDrvLine(wfilt,rbline.buffer(4),srcStride,SRCPTR(xsrc,y+1-1),SRCPTR(xsrc,y+1  ),SRCPTR(xsrc,y+1+1)); 
    DoDrvLine(wfilt,rbline.buffer(5),srcStride,SRCPTR(xsrc,y+2-1),SRCPTR(xsrc,y+2  ),SRCPTR(xsrc,y+2+1)); 

    for (; y<(ye+1); y++)
    {
        // load new scanline of derivatives and filter in vertical axis
        DoDrvLine(wfilt,rbline.buffer(6),srcStride,SRCPTR(xsrc,y+3-1),SRCPTR(xsrc,y+3  ),SRCPTR(xsrc,y+3+1));
        {
            int16_t* src0 = rbline.buffer(0); int16_t* src1 = rbline.buffer(1); int16_t* src2 = rbline.buffer(2);
            int16_t* src3 = rbline.buffer(3); int16_t* src4 = rbline.buffer(4); int16_t* src5 = rbline.buffer(5);
            int16_t* src6 = rbline.buffer(6); 
            int32_t* dst0 = imgDrvH.Ptr(0);
            int32_t* dst1 = imgDrvH.Ptr(1);
            int32_t* dst2 = imgDrvH.Ptr(2);
            int w = wfilt;
#if (defined(_M_IX86) || defined(_M_AMD64))
            HarrisFilterVert7SSE(w,int(srcStride),dst0,dst1,dst2,src0,src1,src2,src3,src4,src5,src6);
#elif defined(_M_ARM)
            HarrisFilterVert7Neon(w,int(srcStride),dst0,dst1,dst2,src0,src1,src2,src3,src4,src5,src6);
#endif
            while (w > 0)
            {
                int32_t acc0, acc1, acc2;
                acc0  = ((int32_t)(*(src0)))<<0; acc1  = ((int32_t)(*(src0+srcStride)))<<0; acc2  = ((int32_t)(*(src0+(2*srcStride))))<<0; // 1/64
                acc0 += ((int32_t)(*(src1)))<<2; acc1 += ((int32_t)(*(src1+srcStride)))<<2; acc2 += ((int32_t)(*(src1+(2*srcStride))))<<2; // 4/64
                acc0 += ((int32_t)(*(src2)))<<3; acc1 += ((int32_t)(*(src2+srcStride)))<<3; acc2 += ((int32_t)(*(src2+(2*srcStride))))<<3; // 8/64
                acc0 += ((int32_t)(*(src3)))<<3; acc1 += ((int32_t)(*(src3+srcStride)))<<3; acc2 += ((int32_t)(*(src3+(2*srcStride))))<<3; // 10/64 = 8/64+2/64
                acc0 += ((int32_t)(*(src3)))<<1; acc1 += ((int32_t)(*(src3+srcStride)))<<1; acc2 += ((int32_t)(*(src3+(2*srcStride))))<<1; //
                acc0 += ((int32_t)(*(src4)))<<3; acc1 += ((int32_t)(*(src4+srcStride)))<<3; acc2 += ((int32_t)(*(src4+(2*srcStride))))<<3;
                acc0 += ((int32_t)(*(src5)))<<2; acc1 += ((int32_t)(*(src5+srcStride)))<<2; acc2 += ((int32_t)(*(src5+(2*srcStride))))<<2;
                acc0 += ((int32_t)(*(src6)))<<0; acc1 += ((int32_t)(*(src6+srcStride)))<<0; acc2 += ((int32_t)(*(src6+(2*srcStride))))<<0;
                src0++; src1++; src2++; src3++; src4++; src5++; src6++;
                *dst0++ = acc0;
                *dst1++ = acc1;
                *dst2++ = acc2;
                w--;
            }
        }
        // filter in horizontal axis
        {
            int32_t* buf0 = imgDrvH.Ptr(0);
            int32_t* buf1 = imgDrvH.Ptr(1);
            int32_t* buf2 = imgDrvH.Ptr(2);
            float* dst = imgH2;

            // normalize for non-even power of two coefficient sum: 2(1+4+8)+10=36; 
            // (36./64.)**2 = .31640625 (squared because 36 summation applied twice); 
            // also scale 1<<6 twice (1<<12) to account for filtering shifts, 
            // and 2. since input is .15 instead of .16
            const float fscale = (1.f/.31640625f)*(1.f/(float)(1<<12))*2.f;
            const float fscaleSq = fscale*fscale;
            int w=whf;
#if 0
#if (defined(_M_IX86) || defined(_M_AMD64))
            HarrisFilterHoriz7SSE(w, dst, buf0, buf1, buf2, fscaleSq);
#elif defined(_M_ARM)
            HarrisFilterHoriz7Neon(w, dst, buf0, buf1, buf2, fscaleSq);
#endif
#endif
            while (w > 0)
            {
                int32_t xx,yy,xy;
                xx  = ((*buf0+0))<<0; xx += (*(buf0+1))<<2; xx += (*(buf0+2))<<3; xx += (*(buf0+3))<<3; xx += (*(buf0+3))<<1; xx += (*(buf0+4))<<3; xx += (*(buf0+5))<<2; xx += (*(buf0+6))<<0; buf0++; 
                yy  = ((*buf1+0))<<0; yy += (*(buf1+1))<<2; yy += (*(buf1+2))<<3; yy += (*(buf1+3))<<3; yy += (*(buf1+3))<<1; yy += (*(buf1+4))<<3; yy += (*(buf1+5))<<2; yy += (*(buf1+6))<<0; buf1++; 
                xy  = ((*buf2+0))<<0; xy += (*(buf2+1))<<2; xy += (*(buf2+2))<<3; xy += (*(buf2+3))<<3; xy += (*(buf2+3))<<1; xy += (*(buf2+4))<<3; xy += (*(buf2+5))<<2; xy += (*(buf2+6))<<0; buf2++; 

                float fxx = (float)xx;
                float fyy = (float)yy;
                float fxy = (float)xy;

                // compute Harris function and store
                *dst++ = ((fxx*fyy)-(fxy*fxy))/(fxx+fyy+HARRIS_DENOM_EPSILON);
                // compute numerator of Harris function and store
//                *dst++ = ((fxx*fyy)-(fxy*fxy))*fscaleSq;
                w--;
            }
        }
        // bump rolling buffer
        rbline.advance();

        // check if done three or more scanlines of harris function - if so then look
        // for maxima in harris function 3x3 neighborhoods for previous scanline
        if ((!useHFBlockCache) && (y == ys))
        {
            // run though the first row to get an initial value for the detect threshold
            float fScanlineAvg = 0.f;
            for(int xh=1; xh<(whf-1); xh++) { fScanlineAvg += imgH1[xh]; }
            fBlockAvg *= (1.f/(float)(VtMax(1,(whf-2))));
            fDetectThreshold = prm.detectThreshold*fBlockAvg;
        }
        if (y > ys)
        {
            float fScanlineAvg = 0.f;
            for(int xh=1; xh<(whf-1); xh++)
            {
                float f = imgH1[xh]; fScanlineAvg += f;
                if( (f > fDetectThreshold) &&
                    (f > imgH1[xh-1]) && (f > imgH1[xh+1]) && 
                    (f > imgH0[xh-1]) && (f > imgH0[xh]) && (f > imgH0[xh+1]) &&
                    (f > imgH2[xh-1]) && (f > imgH2[xh]) && (f > imgH2[xh+1]) &&
                    (*(SRCPTR(xh+xs-1, y-1)) > prm.lowSrcThreshold) )
                {
                    if (useHFBlockCache)
                    {
                        detHit dh; dh.x = xh; dh.y = y; dh.f = f;
                        vecDetHits.push_back(dh);
                    }
                    else
                    {
                        RefineAndStoreResult(hr, pts, xh+xs-1, y-1,
                            imgH0 + xh, imgH1 + xh, imgH2 + xh, prm);
                        ftrCnt++;
                    }
                }
            }
            fScanlineAvg *= (1.f/(float)(VtMax(1,(whf-2))));
            float yFrac = 1.f/(float)(y-ys);
            fBlockAvg = ((yFrac)*fScanlineAvg)+((1.f-yFrac)*fBlockAvg);
            fDetectThreshold = prm.detectThreshold*fBlockAvg;
        }
        // step harris function pointers
        if (useHFBlockCache)
        {
            imgH0 = imgH1; imgH1 = imgH2;
            imgH2 += imgH.StrideBytes()/sizeof(float);
        }
        else
        {
            float* fptmp = imgH0; imgH0 = imgH1; imgH1 = imgH2; imgH2 = fptmp;
        }
    }

    if (useHFBlockCache)
    {
        // find the 'n' detector hits with the highest harris function
        // value; does a partial sort by forming a list of the desired
        // number of values then comparing all of the hits to that list
        // and replacing the smallest when the compared hit is larger;
        // takes a guess at the number of features discarded during
        // refinement culling, so is inexact (but reasonably fast)

// keep this fraction of feature budget during block cache feature selection
// so ~1.0 * the feature budget will be returned after refinement culling;
// 1.5 guesses that 1/3 of the features are culled, which seems to be close
// to true in low-mid 100's of features but less true (with fewer culled)
// with higher counts
#define DETECTCOUNT_SCALE 1.5f 

        int nHits = int(vecDetHits.size());
        int nTopHits = (int)((float)ftrPerBlk*DETECTCOUNT_SCALE);
        vt::vector<detHit>& vecDHReturn = vecDetHits;
        if (nHits > nTopHits)
        {
            vt::vector<detHit> vecTopDetHits; VT_HR_EXIT( vecTopDetHits.reserve(nTopHits) );
            // fill 'top hits' with first nTopHits detector hits, and remember smallest
            int i = 0;
            float fSmallest = FLT_MAX;
            int idxSmallest = 0;
            for (; i<nTopHits; i++) 
            { 
                detHit& dh = vecDetHits[i];
                if (dh.f < fSmallest) { fSmallest = dh.f; idxSmallest = i; }
                vecTopDetHits.push_back(dh); 
            }
            // step through remainder of detector hits and replace smallest in top hits
            // if smaller, then recompute index and value of smallest of top hits
            for (; i<nHits; i++)
            {
                detHit& dh = vecDetHits[i];
                if (dh.f > fSmallest)
                {
                    vecTopDetHits[idxSmallest] = dh;
                    fSmallest = FLT_MAX;
                    for (int j = 0; j < (int)vecTopDetHits.size(); j++)
                    {
                        dh = vecTopDetHits[j];
                        if (dh.f < fSmallest) { fSmallest = dh.f; idxSmallest = j; }
                    }
                }
            }
            vecDHReturn = vecTopDetHits;
        }
        // refine and return top hits
        for (int i=0; i<(int)vecDHReturn.size(); i++)
        {
            detHit& dh = vecDHReturn[i];
            int xh = dh.x;
            int yh = dh.y;
            RefineAndStoreResult(hr, pts, xh+xs-1, yh-1,
                imgH.Ptr(xh,yh-ys-1), imgH.Ptr(xh,yh-ys), imgH.Ptr(xh,yh-ys+1), prm);
            ftrCnt++;
        }
    }

    VT_HR_END();
}
//------------------------------------------------------------------------------

// end