#include "stdafx.h"

#include "completion.h"

#ifdef DEBUG_OUTPUT
#include "vtfileio.h"
#endif

using namespace vt;

void VtScaleSpan(int* pD, const float* pS, float fScale, int iSpan)
{
    int iX = 0;

    if (g_SupportSSE2())
    {
        __m128 x0 = _mm_set1_ps(fScale);

        if (IsAligned16(pD) && IsAligned16(pS))
        {
            for ( ; iX < iSpan - 4; iX += 4, pD += 4, pS += 4 )
            {
                __m128 x4 = _mm_load_ps(pS); 
                x4 = _mm_mul_ps(x4, x0);            // 4 x scale
                __m128i x4i = _mm_cvtps_epi32(x4);  // 4 x float to i32
                _mm_store_si128((__m128i*) pD, x4i);
            }
        }
        else
        {
            for ( ; iX < iSpan - 4; iX += 4, pD += 4, pS += 4 )
            {
                __m128 x4 = _mm_loadu_ps(pS); 
                x4 = _mm_mul_ps(x4, x0);            // 4 x scale
                __m128i x4i = _mm_cvtps_epi32(x4);  // 4 x float to i32
                _mm_storeu_si128((__m128i*) pD, x4i);
            }
        }
    }

    for (; iX < iSpan; iX++, pD++, pS++)
        *pD = F2I( *pS * fScale );
}

//
// Multiply an image by a scalar
// 
void VtScaleImage(OUT CIntImg& imgDst, IN const CFloatImg& imgSrc,  
                  float fScale)
{
    int iSpan = imgSrc.Width() * imgSrc.Bands();

    for (int iY = 0; iY < imgSrc.Height(); iY++)
       VtScaleSpan(imgDst.Ptr(iY), imgSrc.Ptr(iY), fScale, iSpan);
}


template <typename TRGBA>
void patchAvg(float newCol[4], int w, int h, int dSrc,
    const float * fDistWeight, const int * pDist, const RGBAType<TRGBA> ** pPatch,
    int iDist, int iPatch, int srcOff)
{
    for (int y = 0; y <= h; y++)
    {
        for (int x = 0; x <= w; x++)
        {
            const RGBAType<TRGBA> * srcPix = pPatch[x];
            if (srcPix == NULL)
                continue;
            srcPix -= (dSrc + x);

            float dw = fDistWeight[pDist[x]];
            // float dw = powf(weightGamma, (dci - di) * 0.1f);

            newCol[3] += dw;
            newCol[2] += srcPix->r * dw;
            newCol[1] += srcPix->g * dw;
            newCol[0] += srcPix->b * dw;
        }

        pDist += iDist;

        pPatch += iPatch;

        dSrc += srcOff;
    }
}

inline void patchAvg(float newCol[4], int w, int h, int dSrc,
    const float * fDistWeight, const int * pDist, const RGBAPix ** pPatch,
    int iDist, int iPatch, int srcOff)
{
    if (!g_SupportSSE2())
        return patchAvg(newCol, w, h, dSrc, fDistWeight, pDist, pPatch,
                        iDist, iPatch, srcOff);

    __m128i xi0 = _mm_setzero_si128();
    __m128 x2 = _mm_setzero_ps();

    for (int y = 0; y <= h; y++)
    {
        for (int x = 0; x <= w; x++)
        {
            const RGBAPix * srcPix = pPatch[x];
            if (srcPix == NULL)
                continue;
            srcPix -= (dSrc + x);

            float dw = fDistWeight[pDist[x]];
            __m128 x0 = _mm_set1_ps(dw);        // weight

            __m128i xi1 = _mm_set1_epi32(*((int *) srcPix));    // load
            xi1 = _mm_unpacklo_epi8(xi1, xi0);  // convert 8bit to 16bit
            xi1 = _mm_insert_epi16(xi1, 1, 3);  // insert 1 for weight
            xi1 = _mm_unpacklo_epi16(xi1, xi0); // convert 16bit to 32bit
            __m128 x1 = _mm_cvtepi32_ps(xi1);   // convert 32bit to float
            x1 = _mm_mul_ps(x1, x0);            // mul by weight
            x2 = _mm_add_ps(x2, x1);            // add to total
        }

        pDist += iDist;

        pPatch += iPatch;

        dSrc += srcOff;
    }

	_mm_store_ps(newCol, x2);   // store floats
}


template <typename TRGBA, typename TUV>
HRESULT Completion::averageBlock(
    const CRect& rctDst, CCompositeImg<RGBAType<TRGBA>>& imgDst,
    const CRect& rctSrc, const CCompositeImg<RGBAType<TRGBA>>& imgSrc,
    const CRect& rctUv, const CTypedImg<TUV>& uv,
    const CRect& rctMrk, const CByteImg& marked)
{
    VT_HR_BEGIN()

    CImgInfo info = m_imgPyr->GetImgInfo(m_level);
    int w = info.width;
    int h = info.height;

    int dstOff = imgDst.StrideBytes() / sizeof(RGBAType<TRGBA>);
    int srcOff = imgSrc.StrideBytes() / sizeof(RGBAType<TRGBA>);
    int uvOff = uv.StrideBytes() / sizeof(CVec2<TUV>);
    int markOff = marked.StrideBytes();

    float bdDistFactor, invBdDistFactor;
    if (m_level >= m_bdDistLevel)
    {
        bdDistFactor = 1.0f / (2 << (m_level - m_bdDistLevel));
        invBdDistFactor = 0.5f / bdDistFactor;
    }
    else
    {
        bdDistFactor = (float) (2 << (m_bdDistLevel - m_level));
        invBdDistFactor = 1.0f / bdDistFactor;
    }

    CFloatImg fltDist;
    CIntImg intDist;
    CRect rctDist;

    int iPatch;
    const RGBAType<TRGBA> ** ppPatch = NULL;
    const RGBAType<TRGBA> ** ppPtr;
    int px0 = VtMax(s_prad, (int) rctDst.left - s_prad);
    int px1 = VtMin(w - s_prad, (int) rctDst.right + s_prad);
    int py0 = VtMax(s_prad, (int) rctDst.top - s_prad);
    int py1 = VtMin(h - s_prad, (int) rctDst.bottom + s_prad);

    if (m_mark == 0)
    {
        rctDist = CRect(0, 0, rctDst.Width() + 2 * s_prad, rctDst.Height() + 2 * s_prad);
        VT_HR_EXIT( CreateImageForTransform(fltDist, rctDist.Width(), rctDist.Height(), OBJ_FLOATIMG) );
        VT_HR_EXIT( CreateImageForTransform(intDist, rctDist.Width(), rctDist.Height(), OBJ_INTIMG) );

#ifdef DEBUG_MATCH
        if (m_level < m_bdDistLevel)
        {
            for (int y = 0; y < rctDist.Height(); y++)
                for (int x = 0; x < rctDist.Width(); x++)
                    VtSampleBilinear(m_bdDist,
                        (float) (rctDst.left - s_prad + x) * invBdDistFactor,
                        (float) (rctDst.top - s_prad + y) * invBdDistFactor,
                        (float *) NULL, fltDist.Ptr(x, y));

        }
        else
#endif DEBUG_MATCH
        VT_HR_EXIT( VtResizeImage(fltDist, rctDist, m_bdDist,
                                  invBdDistFactor, ((rctDst.left - s_prad)) * invBdDistFactor,
                                  invBdDistFactor, ((rctDst.top - s_prad)) * invBdDistFactor,
                                  m_level < m_bdDistLevel ? eSamplerKernelBilinear : eSamplerKernelNearest,
                                  IMAGE_EXTEND(Zero)) );

        CFloatRoundMode rm(CFloatRoundMode::RM_TRUNCATE);
        ::VtScaleImage(intDist, fltDist, bdDistFactor);

        iPatch = rctDist.Width();
        ppPatch = VT_NOTHROWNEW const RGBAType<TRGBA> *[iPatch * rctDist.Height()];

        CVec2<TUV>* uvPtr = (CVec2<TUV> *) uv.Ptr(py0 - rctUv.top);
        ppPtr = ppPatch;
        for (int py = py0; py < py1; py++)
        {
            for (int px = px0; px < px1; px++)
            {
                TUV nx = uvPtr[px - rctUv.left].x;
                TUV ny = uvPtr[px - rctUv.left].y;
                                
                const RGBAType<TRGBA> * srcPix =
                    nx == ElTraits<TUV>::MaxVal() ? NULL :
                    imgSrc.Ptr(nx - rctSrc.left, ny - rctSrc.top);
                ppPtr[px - px0] = srcPix;
            }

            uvPtr += uvOff;
            ppPtr += iPatch;
        }
    }
    else
    {
        rctDist = CRect(0, 0, 4 + 2 * s_prad, 4 + 2 * s_prad);
        VT_HR_EXIT( CreateImageForTransform(fltDist, rctDist.Width(), rctDist.Height(), OBJ_FLOATIMG) );
        VT_HR_EXIT( CreateImageForTransform(intDist, rctDist.Width(), rctDist.Height(), OBJ_INTIMG) );

        iPatch = 4 + 2 * s_prad;
        ppPatch = VT_NOTHROWNEW const RGBAType<TRGBA> *[iPatch * iPatch];
    }

    const int * pCent;
    int iDist = intDist.StrideBytes() / sizeof(int);

    const Byte * pMark = marked.Ptr();

    for (int my = rctMrk.top; my < rctMrk.bottom; my++)
    {
        for (int mx = rctMrk.left; mx < rctMrk.right; mx++)
        {
            if (m_mark == 0 || pMark[mx - rctMrk.left] == m_mark)
            {
                int x0 = mx * 4, x1 = VtMin(w, x0 + 4);
                int y0 = my * 4, y1 = VtMin(h, y0 + 4);

                if (m_mark != 0)
                {
#ifdef DEBUG_MATCH
                    if (m_level < m_bdDistLevel)
                    {
                        for (int y = 0; y < rctDist.Height(); y++)
                            for (int x = 0; x < rctDist.Width(); x++)
                                VtSampleBilinear(m_bdDist,
                                    (float) (x0 - s_prad + x) * invBdDistFactor,
                                    (float) (y0 - s_prad + y) * invBdDistFactor,
                                    (float *) NULL, fltDist.Ptr(x, y));

                    }
                    else
#endif DEBUG_MATCH
                    VT_HR_EXIT( VtResizeImage(fltDist, rctDist, m_bdDist,
                                              invBdDistFactor, (x0 - s_prad) * invBdDistFactor,
                                              invBdDistFactor, (y0 - s_prad) * invBdDistFactor,
                                              m_level < m_bdDistLevel ? eSamplerKernelBilinear : eSamplerKernelNearest,
                                              IMAGE_EXTEND(Zero)) );

                    CFloatRoundMode rm(CFloatRoundMode::RM_TRUNCATE);
                    ::VtScaleImage(intDist, fltDist, bdDistFactor);
                    pCent = intDist.Ptr(s_prad, s_prad);

                    px0 = VtMax(s_prad, x0 - s_prad);
                    px1 = VtMin(w - s_prad, x1 + s_prad);
                    py0 = VtMax(s_prad, y0 - s_prad);
                    py1 = VtMin(h - s_prad, y1 + s_prad);

                    CVec2<TUV>* uvPtr = (CVec2<TUV> *) uv.Ptr(py0 - rctUv.top);
                    ppPtr = ppPatch;
                    for (int py = py0; py < py1; py++)
                    {
                        for (int px = px0; px < px1; px++)
                        {
                            TUV nx = uvPtr[px - rctUv.left].x;
                            TUV ny = uvPtr[px - rctUv.left].y;
                                
                            const RGBAType<TRGBA> * srcPix =
                                nx == ElTraits<TUV>::MaxVal() ? NULL :
                                imgSrc.Ptr(nx - rctSrc.left, ny - rctSrc.top);
                            ppPtr[px - px0] = srcPix;
                        }

                        uvPtr += uvOff;
                        ppPtr += iPatch;
                    }
                }
                else
                {
                    pCent = intDist.Ptr(s_prad + x0 - rctDst.left, s_prad + y0 - rctDst.top);
                }

                ppPtr = ppPatch - py0 * iPatch - px0;

                RGBAType<TRGBA> * dstPix = imgDst.Ptr(x0 - rctDst.left, y0 - rctDst.top);

                for (int py = y0; py < y1; py++)
                {
                    for (int px = x0; px < x1; px++)
                    {
                        Byte alpha;
                        VtConv(&alpha, dstPix[px - x0].a);
                        if (alpha & s_knownBit)
                        {
                            continue;
                        }
                        
                        // fetch boundary distance at center pixel
                        float fDist = (float) pCent[px - x0];
                        int dci = (int) (fDist - 10 * (2 + s_prad) * (float) VT_SQRT2);
                        const float *fDistWeight = m_bdDistWeight.begin() - dci;

                        int rx0 = VtMax(s_prad, px - s_prad),
                            rx1 = VtMin(w - 1 - s_prad, px + s_prad);
                        int ry0 = VtMax(s_prad, py - s_prad),
                            ry1 = VtMin(h - 1 - s_prad, py + s_prad);

                        // fetch boundary distance at origin pixel
                        const int * pDist = pCent + rx0 - x0 - (py - ry0) * iDist;

                        const RGBAType<TRGBA> ** pPatch = ppPtr + rx0 + ry0 * iPatch;

                        int dSrc = (rx0 - px) + (ry0 - py) * srcOff;

                        VT_DECLSPEC_ALIGN(16) float newCol[4] = { 0 };
                        patchAvg(newCol, rx1 - rx0, ry1 - ry0, dSrc,
                            fDistWeight, pDist, pPatch, iDist, iPatch, srcOff);

#ifdef DEBUG_MATCH
                        newCol[0] /= newCol[3];
                        newCol[1] /= newCol[3];
                        newCol[2] /= newCol[3];
#else DEBUG_MATCH
                        newCol[3] = 1.f / newCol[3];
                        newCol[0] *= newCol[3];
                        newCol[1] *= newCol[3];
                        newCol[2] *= newCol[3];
#endif DEBUG_MATCH

                        dstPix[px - x0].r = (TRGBA) newCol[2];
                        dstPix[px - x0].g = (TRGBA) newCol[1];
                        dstPix[px - x0].b = (TRGBA) newCol[0];
                    }

                    dstPix += dstOff;

                    pCent += iDist;
                }
            }
        }

        pMark += markOff;
    }

    delete[] ppPatch;

    VT_HR_END()
}

#ifdef MULTI_CORE_AVG
class CAverageTransform : public IImageTransform
{
public:
    CAverageTransform(Completion* cpl) : m_cpl(cpl)
    { }

	virtual bool RequiresCloneForConcurrency()
	{ return false;	}

	virtual void    GetSrcPixFormat(IN OUT int* ptypeSrcs, 
									IN UINT /*uSrcCnt*/,
									IN int /*typeDst*/)
    {
        ptypeSrcs[0] = VT_IMG_FIXED( m_cpl->m_imgPyr->GetImgInfo().type );
        ptypeSrcs[1] = VT_IMG_FIXED( m_cpl->m_imgPyr->GetImgInfo().type );
        ptypeSrcs[2] = VT_IMG_FIXED( m_cpl->m_uvPyr.GetImgInfo().type );
        ptypeSrcs[3] = VT_IMG_FIXED( m_cpl->m_maskPyr.GetImgInfo().type );
    }

	virtual void    GetDstPixFormat(OUT int& typeDst,
									IN  const int* ptypeSrcs, 
									IN  UINT /*uSrcCnt*/)
	{ typeDst = ptypeSrcs[0]; }

	virtual HRESULT GetRequiredSrcRect(OUT TRANSFORM_SOURCE_DESC* pSrcReq,
									   OUT UINT& uSrcReqCount,
									   IN  UINT /*uSrcCnt*/,
									   IN  const CRect& rctLayerDst)
    {
        CImgInfo info = m_cpl->m_imgPyr->GetImgInfo(m_cpl->m_level);

        int xBlocks = (info.width + m_cpl->m_blockSize - 1) / m_cpl->m_blockSize;
        int block = (rctLayerDst.left / m_cpl->m_blockSize) + (rctLayerDst.top / m_cpl->m_blockSize) * xBlocks;

        CRect rctSrc = m_cpl->m_bboxes[2][block];
        if (rctSrc.IsRectEmpty())
            return E_UNEXPECTED;
        rctSrc.InflateRect(Completion::s_prad, Completion::s_prad);
        rctSrc &= info.Rect();

        CRect rctUv = rctLayerDst;
        rctUv.InflateRect(Completion::s_prad, Completion::s_prad);
        rctUv &= info.Rect();

        CRect rctMrk = CRect(rctLayerDst.left >> 2, rctLayerDst.top >> 2,
            (rctLayerDst.right + 3) >> 2, (rctLayerDst.bottom + 3) >> 2);
        rctMrk &= m_cpl->m_maskPyr.GetImgInfo(m_cpl->m_level).Rect();

        uSrcReqCount = 4;
        pSrcReq[0].uSrcIndex = 0;
        pSrcReq[0].bCanOverWrite = true;
        pSrcReq[0].rctSrc = rctLayerDst;
        pSrcReq[1].uSrcIndex = 1;
        pSrcReq[1].bCanOverWrite = false;
        pSrcReq[1].rctSrc = rctSrc;
        pSrcReq[2].uSrcIndex = 2;
        pSrcReq[2].bCanOverWrite = false;
        pSrcReq[2].rctSrc = rctUv;
        pSrcReq[3].uSrcIndex = 3;
        pSrcReq[3].bCanOverWrite = false;
        pSrcReq[3].rctSrc = rctMrk;

        return S_OK;
    }

	virtual HRESULT GetAffectedDstRect(OUT CRect& /*rctDst*/,
									   IN  const CRect& /*rctSrc*/,
									   IN  UINT /*uSrcIndex*/,
									   IN  UINT /*uSrcCnt*/)
    { return E_NOTIMPL; }

	virtual HRESULT GetResultingDstRect(OUT CRect& /*rctDst*/,
										IN  const CRect& /*rctSrc*/,
										IN  UINT /*uSrcIndex*/,
										IN  UINT /*uSrcCnt*/)
    { return E_NOTIMPL; }

    virtual HRESULT Transform(OUT CImg* pimgDstRegion, 
							  IN  const CRect& /*rctLayerDst*/,
							  IN  CImg *const *ppimgSrcRegions,
							  IN  const TRANSFORM_SOURCE_DESC* pSrcDesc,
							  IN  UINT /*uSrcCnt*/)
    {
        if (!ppimgSrcRegions[0]->IsSharingMemory(*pimgDstRegion))
            return E_UNEXPECTED;

        CImgInfo info = m_cpl->m_imgPyr->GetImgInfo(m_cpl->m_level);

        #define DO_AVG(EL, TRGBA, TUV) \
            case EL: \
                { \
                    return m_cpl->averageBlock( \
                        pSrcDesc[0].rctSrc, (CCompositeImg<RGBAType<TRGBA>> &) *ppimgSrcRegions[0], \
                        pSrcDesc[1].rctSrc, (const CCompositeImg<RGBAType<TRGBA>> &) *ppimgSrcRegions[1], \
                        pSrcDesc[2].rctSrc, (const CTypedImg<TUV> &) *ppimgSrcRegions[2], \
                        pSrcDesc[3].rctSrc, (CByteImg &) *ppimgSrcRegions[3]); \
                } \
                break;

        if (EL_FORMAT(ppimgSrcRegions[2]->GetType()) == EL_FORMAT_SHORT)
        {
            switch(EL_FORMAT(ppimgSrcRegions[1]->GetType()))
            {
                DO_AVG( EL_FORMAT_BYTE, Byte, UInt16 );
                DO_AVG( EL_FORMAT_SHORT, UInt16, UInt16 );
                DO_AVG( EL_FORMAT_FLOAT, float, UInt16 );
                default:
                    return E_NOTIMPL;
            }
        }
        else
        {
            switch(EL_FORMAT(ppimgSrcRegions[1]->GetType()))
            {
                DO_AVG( EL_FORMAT_BYTE, Byte, int );
                DO_AVG( EL_FORMAT_SHORT, UInt16, int );
                DO_AVG( EL_FORMAT_FLOAT, float, int );
                default:
                    return E_NOTIMPL;
            }
        }
    }

	HRESULT Clone(ITaskState **ppState)
	{
		return CloneTaskState<CAverageTransform>(ppState,
            VT_NOTHROWNEW CAverageTransform(m_cpl));
    }

private:
    Completion* m_cpl;
};
#endif MULTI_CORE_AVG

HRESULT Completion::average()
{
    VT_HR_BEGIN()

    CImgInfo info = m_imgPyr->GetImgInfo(m_level);

#ifdef MULTI_CORE_AVG
    CAverageTransform transform(this);

    CLevelChangeWriter wr(m_imgPyr, m_level);
    CLevelChangeReader rd(m_imgPyr, m_level);
    CLevelChangeReader src(m_imgPyr, m_level);
    CLevelChangeReader uv(&m_uvPyr, m_level);
    CLevelChangeReader mark(&m_maskPyr, m_level);

    CTransformGraphNaryNode node(&transform);

    CUnknownBlockMap bm(info.Rect(), m_blockSize, m_blockSize);
    VT_HR_EXIT( bm.SetKnown(m_bboxes[2]) );
    node.SetDest(NODE_DEST_PARAMS(&wr, &bm));

    VT_HR_EXIT( node.SetSourceCount(4) );
    VT_HR_EXIT( node.BindSourceToReader(0, &rd) );
    VT_HR_EXIT( node.BindSourceToReader(1, &src) );
    VT_HR_EXIT( node.BindSourceToReader(2, &uv) );
    VT_HR_EXIT( node.BindSourceToReader(3, &mark) );

    VT_HR_EXIT( PushTransformTaskAndWait(&node) );

#else MULTI_CORE_AVG
    CImg imgSrc, imgDst;

    CImg uv;
    CByteImg marked;

    #define DO_AVG(EL, TRGBA, TUV) \
        case EL: \
            VT_HR_EXIT( averageBlock( \
                rctDst, (CCompositeImg<RGBAType<TRGBA>> &) imgDst, \
                rctSrc, (CCompositeImg<RGBAType<TRGBA>> &) imgSrc, \
                rctUv, (CTypedImg<TUV> &) uv, rctMrk, marked) ); \
            break;

    int block = 0;
    for (CBlockIterator bi(info.Rect(), m_blockSize); !bi.Done(); bi.Advance())
    {
        // find source bounding box
        CRect rctSrc = m_bboxes[2][block++]; 
        if (rctSrc.IsRectEmpty())
            continue;

        CRect rctDst = bi.GetRect();
        VT_HR_EXIT( m_imgPyr->ReadRegion(rctDst, imgDst, m_level) );

        CRect rctMrk = CRect(rctDst.left >> 2, rctDst.top >> 2, (rctDst.right + 3) >> 2, (rctDst.bottom + 3) >> 2);
        rctMrk &= m_maskPyr.GetImgInfo(m_level).Rect();
        VT_HR_EXIT( m_maskPyr.ReadRegion(rctMrk, marked, m_level) );

        CRect rctUv = rctDst;
        rctUv.InflateRect(s_prad, s_prad);
        rctUv &= info.Rect();
        VT_HR_EXIT( m_uvPyr.ReadRegion(rctUv, uv, m_level) );
 
        // read source pixels
        rctSrc.InflateRect(s_prad, s_prad);
        rctSrc &= info.Rect();
        VT_HR_EXIT( m_imgPyr->ReadRegion(rctSrc, imgSrc, m_level) );

        if (EL_FORMAT(m_uvPyr.GetImgInfo().type) == EL_FORMAT_SHORT)
        {
            switch(EL_FORMAT(info.type))
            {
            DO_AVG( EL_FORMAT_BYTE, Byte, UInt16 );
            DO_AVG( EL_FORMAT_SHORT, UInt16, UInt16 );
            DO_AVG( EL_FORMAT_FLOAT, float, UInt16 );
            default:
                VT_HR_EXIT( E_NOTIMPL );
                break;
            }
        }
        else
        {
            switch(EL_FORMAT(info.type))
            {
            DO_AVG( EL_FORMAT_BYTE, Byte, int );
            DO_AVG( EL_FORMAT_SHORT, UInt16, int );
            DO_AVG( EL_FORMAT_FLOAT, float, int );
            default:
                VT_HR_EXIT( E_NOTIMPL );
                break;
            }
        }

        VT_HR_EXIT( m_imgPyr->WriteRegion(rctDst, imgDst, m_level) );
    }

#ifdef DEBUG_OUTPUT
    if (m_dbgSaveAverageData)
    {
        CByteImg vis;
        VT_HR_EXIT( m_maskPyr.ReadImg(vis, m_level) );
        for (int y = 0; y < vis.Height(); y++)
        {
            for (int x = 0; x < vis.Width(); x++)
            {
                if (m_mark == 0)
                {
                    vis(x, y) = ((x+y)%2)*255;
                }
                else
                {
                    if (vis(x, y) == m_mark + 1)
                    {
                        vis(x, y) = 128;
                    }
                    else if (vis(x, y) == m_mark +2 )
                    {
                        vis(x, y) = 255;
                    }
                    else if (vis(x, y) == m_mark)
                    {
                        vis(x, y) = 64;
                    }
                    else
                    {
                        vis(x, y) = 0;
                    }
                }
            }
        }

        wstring fileName = wformat(L"output\\d_avg_lvl%02d_gen%04d.png", m_level, m_dbgGen);
        VT_HR_EXIT( VtSaveImage(fileName, vis) );
    }
#endif DEBUG_OUTPUT

#endif MULTI_CORE_AVG

    // Clear average bounding boxes after average pass.
    for (int i = 0; i < (int) m_bboxes[2].size(); i++)
        m_bboxes[2][i] = CRect(LONG_MAX, LONG_MAX, 0, 0);

    VT_HR_END()
}
