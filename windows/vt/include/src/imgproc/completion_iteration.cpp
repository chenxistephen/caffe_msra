#include "stdafx.h"

#include "completion.h"

using namespace vt;

#ifdef DEBUG_OUTPUT
#include "vtfileio.h"

HRESULT Completion::saveUv(const wstring & fileName)
{
    VT_HR_BEGIN()

    CImg temp;
    VT_HR_EXIT( m_uvPyr.ReadImg(temp, m_level) );

    if (EL_FORMAT(temp.GetType()) == EL_FORMAT_SHORT)
    {
        CShortImg& imgShort = (CShortImg &) temp;
        for (int y = 0; y < temp.Height(); y++)
        {
            for (int x = 0; x < temp.Width(); x++)
            {
                if (imgShort(x, y, 0) == 0xffff)
                {
                    imgShort(x, y, 0) = imgShort(x, y, 1) = 0;
                }
            }
        }
    }
    else
    {
        CIntImg& imgInt = (CIntImg &) temp;
        for (int y = 0; y < temp.Height(); y++)
        {
            for (int x = 0; x < temp.Width(); x++)
            {
                if (imgInt(x, y, 0) == 0xffffffff)
                {
                    imgInt(x, y, 0) = imgInt(x, y, 1) = 0;
                }
            }
        }
    }

    VT_HR_EXIT( VtSaveImage(fileName, temp) );

    VT_HR_END()
}
#endif DEBUG_OUTPUT

template <typename TUV>
void UpdateBBox(CRect & rctSrc, const CTypedImg<TUV> & uv, const CRect & rctUv);

template <typename TELT> struct SADType;

template<> struct SADType<Byte> 
{
	typedef UINT64 T;
};
template<> struct SADType<UInt16> 
{
	typedef UINT64 T;
};
template<> struct SADType<float> 
{
	typedef double T;
};

// Load and SAD first 16 bytes of row.
#define SSE_SAD_16                          \
    xi1 = _mm_loadu_si128((__m128i *) res); \
    xi2 = _mm_loadu_si128((__m128i *) src); \
    xi1 = _mm_and_si128(xi1, xi16);         \
    xi2 = _mm_and_si128(xi2, xi16);         \
    xi2 = _mm_sad_epu8(xi2, xi1);           \
    xi0 = _mm_add_epi64(xi0, xi2);

// Load and SAD last 12 bytes of row.
#define SSE_SAD_12                                \
    xi1 = _mm_loadu_si128((__m128i *) (res + 3)); \
    xi2 = _mm_loadu_si128((__m128i *) (src + 3)); \
    xi1 = _mm_and_si128(xi1, xi12);               \
    xi2 = _mm_and_si128(xi2, xi12);               \
    xi2 = _mm_sad_epu8(xi2, xi1);                 \
    xi0 = _mm_add_epi64(xi0, xi2);

template <typename TRES, typename TRGB, typename TUV>
TRES patchDist(const RGBAType<TRGB> * res, int resoff,
               const RGBAType<TRGB> * src, int srcoff,
               const CVec2<TUV> * puv, int uvoff,
               const TUV sx, const TUV sy,
               TRES discontinuousPenalty)
{
    TRES cost = 0;

    if (g_SupportSSE2() && ElTraits<TRGB>::ElFormat() == EL_FORMAT_BYTE)
    {
        // alpha masks
        __m128i xi16 = _mm_set_epi32(0xffffff, 0xffffff, 0xffffff, 0xffffff);
        __m128i xi12 = _mm_set_epi32(0xffffff, 0xffffff, 0xffffff, 0x000000);
        // total
        __m128i xi0 = _mm_setzero_si128();
        __m128i xi1, xi2;

        SSE_SAD_16;
        SSE_SAD_12;
        res += resoff;
        src += srcoff;

        SSE_SAD_16;
        SSE_SAD_12;
        res += resoff;
        src += srcoff;

        SSE_SAD_16;
        SSE_SAD_12;
        res += resoff;
        src += srcoff;

        SSE_SAD_16;
        SSE_SAD_12;
        res += resoff;
        src += srcoff;

        SSE_SAD_16;
        SSE_SAD_12;
        res += resoff;
        src += srcoff;

        SSE_SAD_16;
        SSE_SAD_12;
        res += resoff;
        src += srcoff;

        SSE_SAD_16;
        SSE_SAD_12;

        int c0 = _mm_extract_epi16(xi0, 0);
        int c4 = _mm_extract_epi16(xi0, 4);
        cost = (TRES) (c0 + c4);
    }
    else
    {
        for (int i = 0; i < Completion::s_psize; i++)
        {
            cost += VtSADSpan<TRGB, TRES>(res, src, Completion::s_psize);
            res += resoff;
            src += srcoff;
        }
    }

    // penalize discontinous pixels
    if (puv[-1].x + 1 != sx || puv[-1].y != sy)
    {
        cost += discontinuousPenalty;
    }
    if (puv[1].x - 1 != sx || puv[1].y != sy)
    {
        cost += discontinuousPenalty;
    }
    if (puv[-uvoff].x != sx || puv[-uvoff].y + 1 != sy)
    {
        cost += discontinuousPenalty;
    }
    if (puv[+uvoff].x != sx || puv[+uvoff].y - 1 != sy)
    {
        cost += discontinuousPenalty;
    }

    return cost;
}

template <typename TRGBA, typename TUV>
HRESULT Completion::iterateBlock(int iter, int seed,
    const CRect& rctCur, const CRect& rctUv, CTypedImg<TUV>& uv,
    const CRect& rctDst, const CCompositeImg<RGBAType<TRGBA>>& imgDst,
    const CRect& rctSrc, const CCompositeImg<RGBAType<TRGBA>>& imgSrc
#ifdef DEBUG_OUTPUT
    , FILE * flog
#endif DEBUG_OUTPUT
    )
{
    VT_HR_BEGIN()

    SADType<TRGBA>::T discontinuousPenalty;
    SADType<TRGBA>::T windowEnergyThresh;
    if (ElTraits<TRGBA>::ElFormat() == EL_FORMAT_FLOAT)
    {
        discontinuousPenalty = (SADType<TRGBA>::T) m_discontinuousPenalty;
        windowEnergyThresh = (SADType<TRGBA>::T) m_windowEnergyThresh;
    }
    else
    {
        discontinuousPenalty = (SADType<TRGBA>::T)
            (m_discontinuousPenalty * (float) (ElTraits<TRGBA>::MaxVal() + 1));
        windowEnergyThresh = (SADType<TRGBA>::T)
            (m_windowEnergyThresh * (float) (ElTraits<TRGBA>::MaxVal() + 1));
    }

    int tileShift = m_tileLevel + Tiles::levelOffset - m_level;
    int tileWidth = m_restriction[0].Width();

    int srcoff = imgSrc.StrideBytes() / sizeof(RGBAType<TRGBA>);
    int dstoff = imgDst.StrideBytes() / sizeof(RGBAType<TRGBA>);
    int uvoff = uv.StrideBytes() / sizeof(CVec2<TUV>);

    short idxAdv;
    int uvAdv, dstAdv;
    bool downPass = (iter & 1) == 0;
    if (downPass)
    {
        idxAdv = +1;
        uvAdv = +uvoff;
        dstAdv = +dstoff;
    }
    else
    {
        idxAdv = -1;
        uvAdv = -uvoff;
        dstAdv = -dstoff;
    }

    // May need to read random trial images.
    CCompositeImg<RGBAType<TRGBA>> imgTrial;

    bool blockchanged = false;

    CImgInfo info = m_imgPyr->GetImgInfo(m_level);
    int w = info.width;
    int h = info.height;
    int xBlocks = (w + m_blockSize - 1) / m_blockSize;
    int block = (rctCur.left / m_blockSize) + (rctCur.top / m_blockSize) * xBlocks;

#ifdef FAST_RAND
    FastRand rnd;
#else FAST_RAND
    CRand rnd;
#endif FAST_RAND

#ifndef USE_TILE_BBOXES
    CRect rctBbx = m_bboxes[0][block] | m_bboxes[1][block];
    if (downPass)
        rctBbx.InflateRect(0, 0, 1, 1);
    else
        rctBbx.InflateRect(1, 1, 0, 0);
    rctBbx &= CRect(0, 0, w, h);
#endif USE_TILE_BBOXES

    CByteImg marked;
    CRect rctMrk = CRect((rctCur.left - s_prad) >> 2,
                         (rctCur.top  - s_prad) >> 2,
                         ((rctCur.right  + s_prad) >> 2) + 1,
                         ((rctCur.bottom + s_prad) >> 2) + 1);
    rctMrk &= m_maskPyr.GetImgInfo(m_level).Rect();
    VT_HR_EXIT( m_maskPyr.ReadRegion(rctMrk, marked, m_level) );
    int markoff = marked.StrideBytes();

    CRect rctOut = CRect(VtMax((int) rctCur.left, s_prad), VtMax((int) rctCur.top, s_prad),
                         VtMin((int) rctCur.right, w - s_prad), VtMin((int) rctCur.bottom, h - s_prad));
    rctOut &= info.Rect();

    int x0 = (downPass ? rctOut.left : rctOut.right - 1);
    int x1 = (downPass ? rctOut.right : rctOut.left - 1);
    int y0 = (downPass ? rctOut.top : rctOut.bottom - 1);
    int y1 = (downPass ? rctOut.bottom : rctOut.top - 1);

    CVec2<TUV> * puvy = (CVec2<TUV> *) uv.Ptr(y0 - rctUv.top);
    const RGBAType<TRGBA> * pdsty = imgDst.Ptr(y0 - s_prad - rctDst.top);

    for (int y = y0; y != y1; y += idxAdv, puvy += uvAdv, pdsty += dstAdv)
    {
        int tiley = tileWidth * (y >> tileShift);

        CVec2<TUV> * puv = puvy + x0 - rctUv.left;
        const RGBAType<TRGBA> * pdst = pdsty + x0 - s_prad - rctDst.left;

        for (int x = x0; x != x1; x += idxAdv, puv += idxAdv, pdst += idxAdv)
        {
            TUV nx = puv->x, ny = puv->y;
            if (nx == ElTraits<TUV>::MaxVal())
                continue;
            TUV nhx = puv[-idxAdv].x, nhy = puv[-idxAdv].y;
            TUV nvx = puv[-uvAdv].x, nvy = puv[-uvAdv].y;

            bool hcont = (nhx == ElTraits<TUV>::MaxVal() || (nhx + idxAdv == nx && nhy == ny));
            bool vcont = (nvx == ElTraits<TUV>::MaxVal() || (nvx == nx && nvy + idxAdv == ny));

            // only do work at patch edges
            if (hcont && vcont)
            {
                continue;
            }

            int tile = tiley + (x >> tileShift);
            CFloatImg& restriction = m_restriction[tile];

            SADType<TRGBA>::T curDist = INT_MAX;

            const RGBAType<TRGBA> * src = NULL;
            if (restriction(nx >> tileShift, ny >> tileShift) != FLT_MAX)
            {
                src = imgSrc.Ptr(nx - rctSrc.left - s_prad, ny - rctSrc.top - s_prad);

                curDist = patchDist<SADType<TRGBA>::T, TRGBA, TUV>
                    (pdst, dstoff, src, srcoff, puv, uvoff, nx, ny, discontinuousPenalty);
            }

#ifdef DEBUG_OUTPUT
            if (m_dbgSaveIterationData)
            {
                fprintf(flog, "examine %d, %d @ %d\n", x, y, curDist);
            }
#endif DEBUG_OUTPUT

            bool changed = false;

            // propagage horizontal
            if (nhx != ElTraits<TUV>::MaxVal())
            {
                const TUV px = nhx + idxAdv;
                const TUV py = nhy;

#ifdef DEBUG_OUTPUT
                if (m_dbgSaveIterationData)
                {
                    fprintf(flog, "horz %d, %d --> %d, %d\n", x, y, px, py);
                }
#endif DEBUG_OUTPUT

                if (px == nx && py == ny)
                {
                    goto skipHorizontalPropagation;
                }

                if (px < s_prad || px >= w - s_prad)
                {
                    goto skipHorizontalPropagation;
                }

                if (restriction(px >> tileShift, py >> tileShift) == FLT_MAX)
                {
                    goto skipHorizontalPropagation;
                }

                Byte alpha;
                int stride;
#ifndef USE_TILE_BBOXES
                if (!blockchanged || rctBbx.PtInRect(CPoint(px, py)))
#endif USE_TILE_BBOXES
                {
                    VtConv(&alpha, imgSrc(px - rctSrc.left, py - rctSrc.top, 3));
                    if ((alpha & s_validBit) == 0)
                        goto skipHorizontalPropagation;

                    src = imgSrc.Ptr(px - rctSrc.left - s_prad, py - rctSrc.top - s_prad);
                    stride = srcoff;
                }
#ifndef USE_TILE_BBOXES
                else
                {
                    VT_HR_EXIT( m_imgPyr->ReadRegion(CRect(px - s_prad, py - s_prad,
                                                           px + s_prad + 1, py + s_prad + 1),
                                                     imgTrial, m_level) );
                    VtConv(&alpha, imgTrial(s_prad, s_prad, 3));
                    if ((alpha & s_validBit) == 0)
                        goto skipHorizontalPropagation;

                    src = imgTrial.Ptr();
                    stride = imgTrial.StrideBytes() / sizeof(RGBAType<TRGBA>);
                }
#endif USE_TILE_BBOXES

                SADType<TRGBA>::T newDist =
                    patchDist<SADType<TRGBA>::T, TRGBA, TUV>
                        (pdst, dstoff, src, stride, puv, uvoff, px, py, discontinuousPenalty);

#ifdef DEBUG_OUTPUT
                if (m_dbgSaveIterationData)
                {
                    fprintf(flog, "    dist %d\n", newDist);
                }
#endif DEBUG_OUTPUT

                if (newDist < curDist)
                {
                    curDist = newDist;
                    puv->x = px;
                    puv->y = py;
                    nx = px;
                    ny = py;
                    blockchanged = changed = true;
#ifdef DEBUG_OUTPUT
                    if (m_dbgSaveIterationData)
                    {
                        fprintf(flog, "    accept (%d, %d) HORZ (%d, %d) @ %d\n", x, y, px, py, newDist);
                    }
#endif DEBUG_OUTPUT
                }
            }
skipHorizontalPropagation:

            // propagage vertical
            if (nvx != ElTraits<TUV>::MaxVal())
            {
                const TUV px = nvx;
                const TUV py = nvy + idxAdv;

#ifdef DEBUG_OUTPUT
                if (m_dbgSaveIterationData)
                {
                    fprintf(flog, "vert %d, %d --> %d, %d\n", x, y, px, py);
                }
#endif DEBUG_OUTPUT

                if (px == nx && py == ny)
                {
                    goto skipVerticalPropagation;
                }

                if (py < s_prad || py >= h - s_prad)
                {
                    goto skipVerticalPropagation;
                }

                if (restriction(px >> tileShift, py >> tileShift) == FLT_MAX)
                {
                    goto skipVerticalPropagation;
                }

                Byte alpha;
                int stride;
#ifndef USE_TILE_BBOXES
                if (!blockchanged || rctBbx.PtInRect(CPoint(px, py)))
#endif USE_TILE_BBOXES
                {
                    VtConv(&alpha, imgSrc(px - rctSrc.left, py - rctSrc.top, 3));
                    if ((alpha & s_validBit) == 0)
                        goto skipVerticalPropagation;

                    src = imgSrc.Ptr(px - rctSrc.left - s_prad, py - rctSrc.top - s_prad);
                    stride = srcoff;
                }
#ifndef USE_TILE_BBOXES
                else
                {
                    VT_HR_EXIT( m_imgPyr->ReadRegion(CRect(px - s_prad, py - s_prad,
                                                           px + s_prad + 1, py + s_prad + 1),
                                                     imgTrial, m_level) );
                    VtConv(&alpha, imgTrial(s_prad, s_prad, 3));
                    if ((alpha & s_validBit) == 0)
                        goto skipVerticalPropagation;

                    src = imgTrial.Ptr();
                    stride = imgTrial.StrideBytes() / sizeof(RGBAType<TRGBA>);
                }
#endif USE_TILE_BBOXES

                SADType<TRGBA>::T newDist =
                    patchDist<SADType<TRGBA>::T, TRGBA, TUV>
                        (pdst, dstoff, src, stride, puv, uvoff, px, py, discontinuousPenalty);

#ifdef DEBUG_OUTPUT
                if (m_dbgSaveIterationData)
                {
                    fprintf(flog, "    dist %d\n", newDist);
                }
#endif DEBUG_OUTPUT

                if (newDist < curDist)
                {
                    curDist = newDist;
                    puv->x = px;
                    puv->y = py;
                    nx = px;
                    ny = py;
                    blockchanged = changed = true;
#ifdef DEBUG_OUTPUT
                    if (m_dbgSaveIterationData)
                    {
                        fprintf(flog, "    accept (%d, %d) VERT (%d, %d) @ %d\n", x, y, px, py, newDist);
                    }
#endif DEBUG_OUTPUT
                }
            }
skipVerticalPropagation:{}

            // random search window based
            if (curDist >= windowEnergyThresh)
            {
                rnd.Seed(seed | (x << 16) | y);

                int radius = VtMax(w, h);
                for (; radius > 0; radius /= 2)
                {
                    if (curDist < windowEnergyThresh)
                    {
                        break;
                    }

                    const TUV px = (TUV) (nx - radius + rnd.IRand(2 * radius + 1));
                    const TUV py = (TUV) (ny - radius + rnd.IRand(2 * radius + 1));

#ifdef DEBUG_OUTPUT
                    if (m_dbgSaveIterationData)
                    {
                        fprintf(flog, "rnd %d, %d --> %d, %d\n", x, y, px, py);
                    }
#endif DEBUG_OUTPUT

                    if (px == nx && py == ny)
                    {
                        continue;
                    }

                    if (px < s_prad || px >= w - s_prad ||
                        py < s_prad || py >= h - s_prad)
                    {
                        continue;
                    }

                    if (restriction(px >> tileShift, py >> tileShift) == FLT_MAX)
                    {
                        continue;
                    }

                    Byte alpha;
                    int stride;
#ifndef USE_TILE_BBOXES
                    if (rctBbx.PtInRect(CPoint(px, py)))
#endif USE_TILE_BBOXES
                    {
                        VtConv(&alpha, imgSrc(px - rctSrc.left, py - rctSrc.top, 3));
                        if ((alpha & s_validBit) == 0)
                            continue;

                        src = imgSrc.Ptr(px - rctSrc.left - s_prad, py - rctSrc.top - s_prad);
                        stride = srcoff;
                    }
#ifndef USE_TILE_BBOXES
                    else
                    {
                        VT_HR_EXIT( m_imgPyr->ReadRegion(CRect(px - s_prad, py - s_prad,
                                                               px + s_prad + 1, py + s_prad + 1),
                                                         imgTrial, m_level) );
                        VtConv(&alpha, imgTrial(s_prad, s_prad, 3));
                        if ((alpha & s_validBit) == 0)
                            continue;

                        src = imgTrial.Ptr();
                        stride = imgTrial.StrideBytes() / sizeof(RGBAType<TRGBA>);
                    }
#endif USE_TILE_BBOXES

                    SADType<TRGBA>::T newDist =
                        patchDist<SADType<TRGBA>::T, TRGBA, TUV>
                            (pdst, dstoff, src, stride, puv, uvoff, px, py, discontinuousPenalty);

#ifdef DEBUG_OUTPUT
                    if (m_dbgSaveIterationData)
                    {
                        fprintf(flog, "    dist %d\n", newDist);
                    }
#endif DEBUG_OUTPUT

                    if (newDist < curDist)
                    {
                        curDist = newDist;
                        puv->x = px;
                        puv->y = py;
                        nx = px;
                        ny = py;
                        blockchanged = changed = true;
#ifdef DEBUG_OUTPUT
                        if (m_dbgSaveIterationData)
                        {
                            fprintf(flog, "    accept (%d, %d) RAND (%d, %d) @ %d\n", x, y, px, py, newDist);
                        }
#endif DEBUG_OUTPUT
                    }
                }
            }

            if (changed)
            {
                int mx0 = ((x - s_prad) >> 2) - rctMrk.left, mx1 = ((x + s_prad) >> 2) - rctMrk.left;
                int my0 = ((y - s_prad) >> 2) - rctMrk.top,  my1 = ((y + s_prad) >> 2) - rctMrk.top;

                Byte * pmark = marked.Ptr(my0);
                for (int my = my0; my <= my1; my++)
                {
                    for (int mx = mx0; mx <= mx1; mx++)
                    {
                        pmark[mx] = m_mark;
                    }
                    pmark += markoff;
                }
            }
        }
    }

    // Set next pass bounding box.
    int phase = downPass ? 0 : 1;
    if (blockchanged)
        UpdateBBox(m_bboxes[1 - phase][block], uv, rctCur - rctUv.TopLeft());
    else
        m_bboxes[1 - phase][block] |= m_bboxes[phase][block];

    // iteration bounding box is 1 larger along trailing edges

    // add right edge for down pass
    if (rctCur.right != w)
        UpdateBBox(m_bboxes[0][block + 1], uv, CRect(rctCur.right - 1, rctCur.top, rctCur.right, rctCur.bottom) - rctUv.TopLeft());
    // add bottom edge for down pass
    if (rctCur.bottom != h)
        UpdateBBox(m_bboxes[0][block + xBlocks], uv, CRect(rctCur.left, rctCur.bottom - 1, rctCur.right, rctCur.bottom) - rctUv.TopLeft());
    // add left edge for up pass
    if (rctCur.left != 0)
        UpdateBBox(m_bboxes[1][block - 1], uv, CRect(rctCur.left, rctCur.top, rctCur.left + 1, rctCur.bottom) - rctUv.TopLeft());
    // add top edge for up pass
    if (rctCur.top != 0)
        UpdateBBox(m_bboxes[1][block - xBlocks], uv, CRect(rctCur.left, rctCur.top, rctCur.right, rctCur.top + 1) - rctUv.TopLeft());

    if (iter == s_numIters - 1)
    {
        // Time for sparse average; find bounding box of marks.
        CRect bbox = CRect(LONG_MAX, LONG_MAX, 0, 0);

        Byte * pmark = marked.Ptr();
        for (int my = 0; my < rctMrk.Height(); my++)
        {
            for (int mx = 0; mx < rctMrk.Width(); mx++)
            {
                if (pmark[mx] == m_mark)
                {
                    bbox.left   = VtMin(bbox.left,   (long) mx);
                    bbox.top    = VtMin(bbox.top,    (long) my);
                    bbox.right  = VtMax(bbox.right,  (long) mx + 1);
                    bbox.bottom = VtMax(bbox.bottom, (long) my + 1);
                }
            }
            pmark += markoff;
        }

        if (!bbox.IsRectEmpty())
        {
            bbox += rctMrk.TopLeft();
            bbox.left   <<= 2;
            bbox.top    <<= 2;
            bbox.right  <<= 2;
            bbox.bottom <<= 2;
            bbox.InflateRect(s_prad, s_prad);
            bbox &= info.Rect();

            UpdateBBox(m_bboxes[2][block], uv, bbox - rctUv.TopLeft());

            // average bounding box is s_prad larger all around

            // add right edge
            if (rctCur.right != w)
                UpdateBBox(m_bboxes[2][block + 1], uv, (bbox &
                    CRect(rctCur.right - s_prad, rctUv.top, rctUv.right, rctUv.bottom)) - rctUv.TopLeft());
            // add bottom edge
            if (rctCur.bottom != h)
                UpdateBBox(m_bboxes[2][block + xBlocks], uv, (bbox &
                    CRect(rctUv.left, rctCur.bottom - s_prad, rctUv.right, rctUv.bottom)) - rctUv.TopLeft());
            // add bottom right corner
            if (rctCur.right != w && rctCur.bottom != h)
                UpdateBBox(m_bboxes[2][block + xBlocks + 1], uv, (bbox &
                    CRect(rctCur.right - s_prad, rctCur.bottom - s_prad, rctUv.right, rctUv.bottom)) - rctUv.TopLeft());
            // add bottom left corner
            if (rctCur.left != 0 && rctCur.bottom != h)
                UpdateBBox(m_bboxes[2][block + xBlocks - 1], uv, (bbox &
                    CRect(rctUv.left, rctCur.bottom - s_prad, rctCur.left + s_prad, rctUv.bottom)) - rctUv.TopLeft());
            // add left edge
            if (rctCur.left != 0)
                UpdateBBox(m_bboxes[2][block - 1], uv, (bbox &
                    CRect(rctUv.left, rctUv.top, rctCur.left + s_prad, rctUv.bottom)) - rctUv.TopLeft());
            // add top edge
            if (rctCur.top != 0)
                UpdateBBox(m_bboxes[2][block - xBlocks], uv, (bbox &
                    CRect(rctUv.left, rctUv.top, rctUv.right, rctCur.top + s_prad)) - rctUv.TopLeft());
            // add top left corner
            if (rctCur.left != 0 && rctCur.top != 0)
                UpdateBBox(m_bboxes[2][block - xBlocks - 1], uv, (bbox &
                    CRect(rctUv.left, rctUv.top, rctCur.left + s_prad, rctCur.top + s_prad)) - rctUv.TopLeft());
            // add top right corner
            if (rctCur.right != w && rctCur.top != 0)
                UpdateBBox(m_bboxes[2][block - xBlocks + 1], uv, (bbox &
                    CRect(rctCur.right - s_prad, rctUv.top, rctUv.right, rctCur.top + s_prad)) - rctUv.TopLeft());
        }
    }

    if (blockchanged)
    {
        VT_HR_EXIT( m_maskPyr.WriteRegion(rctMrk, marked, m_level) );

        m_changed = true;
    }

    VT_HR_END()
}

#ifdef MULTI_CORE_ITERATION
class CIterationTransform : public IImageTransform
{
public:
    CIterationTransform(Completion* cpl, int iter, int seed) :
      m_cpl(cpl), m_iter(iter), m_seed(seed)
    { }

	virtual bool RequiresCloneForConcurrency()
	{ return false;	}

	virtual void    GetSrcPixFormat(IN OUT int* ptypeSrcs, 
									IN UINT /*uSrcCnt*/,
									IN int /*typeDst*/)
    {
        ptypeSrcs[0] = VT_IMG_FIXED( m_cpl->m_uvPyr.GetImgInfo().type );
        ptypeSrcs[1] = VT_IMG_FIXED( m_cpl->m_imgPyr->GetImgInfo().type );
        ptypeSrcs[2] = VT_IMG_FIXED( m_cpl->m_imgPyr->GetImgInfo().type );
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

        CRect rctBbx = m_cpl->m_bboxes[0][block] | m_cpl->m_bboxes[1][block];
        if (rctBbx.IsRectEmpty())
            return E_UNEXPECTED;
        bool downPass = (m_iter & 1) == 0;
        if (downPass)
            rctBbx.InflateRect(0, 0, 1, 1);
        else
            rctBbx.InflateRect(1, 1, 0, 0);
        rctBbx &= info.Rect();

        CRect rctSrc = rctBbx;
        rctSrc.InflateRect(Completion::s_prad, Completion::s_prad);
        rctSrc &= info.Rect();

        CRect rctDst = rctLayerDst;
        rctDst.InflateRect(Completion::s_prad, Completion::s_prad);
        rctDst &= info.Rect();

        CRect rctUv;
        if (m_iter == Completion::s_numIters - 1)
        {
            // Time for sparse average; will need bounding box of marks.
            CRect rctMrk = CRect((rctLayerDst.left - Completion::s_prad) >> 2,
                                 (rctLayerDst.top  - Completion::s_prad) >> 2,
                                 ((rctLayerDst.right  + Completion::s_prad) >> 2) + 1,
                                 ((rctLayerDst.bottom + Completion::s_prad) >> 2) + 1);
            rctUv = CRect(rctMrk.left << 2, rctMrk.top << 2,
                          rctMrk.right << 2, rctMrk.bottom << 2);
            rctUv.InflateRect(Completion::s_prad, Completion::s_prad);
        }
        else
        {
            rctUv = rctLayerDst;
            rctUv.InflateRect(1, 1);
        }
        rctUv &= info.Rect();

        uSrcReqCount = 3;
        pSrcReq[0].uSrcIndex = 0;
        pSrcReq[0].bCanOverWrite = true;
        pSrcReq[0].rctSrc = rctUv;
        pSrcReq[1].uSrcIndex = 1;
        pSrcReq[1].bCanOverWrite = false;
        pSrcReq[1].rctSrc = rctDst;
        pSrcReq[2].uSrcIndex = 2;
        pSrcReq[2].bCanOverWrite = false;
        pSrcReq[2].rctSrc = rctSrc;

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
							  IN  const CRect& rctLayerDst,
							  IN  CImg *const *ppimgSrcRegions,
							  IN  const TRANSFORM_SOURCE_DESC* pSrcDesc,
							  IN  UINT /*uSrcCnt*/)
    {
        if (!ppimgSrcRegions[0]->IsSharingMemory(*pimgDstRegion))
            return E_UNEXPECTED;

        CImgInfo info = m_cpl->m_imgPyr->GetImgInfo(m_cpl->m_level);

        #define DO_IT(EL, TRGBA, TUV) \
            case EL: \
                { \
                    return m_cpl->iterateBlock(m_iter, m_seed, rctLayerDst, \
                        pSrcDesc[0].rctSrc, (CTypedImg<TUV> &) *ppimgSrcRegions[0], \
                        pSrcDesc[1].rctSrc, (const CCompositeImg<RGBAType<TRGBA>> &) *ppimgSrcRegions[1], \
                        pSrcDesc[2].rctSrc, (const CCompositeImg<RGBAType<TRGBA>> &) *ppimgSrcRegions[2]); \
                } \
                break;

        if (EL_FORMAT(ppimgSrcRegions[0]->GetType()) == EL_FORMAT_SHORT)
        {
            switch(EL_FORMAT(ppimgSrcRegions[1]->GetType()))
            {
                DO_IT( EL_FORMAT_BYTE, Byte, UInt16 );
                DO_IT( EL_FORMAT_SHORT, UInt16, UInt16 );
                DO_IT( EL_FORMAT_FLOAT, float, UInt16 );
                default:
                    return E_NOTIMPL;
            }
        }
        else
        {
            switch(EL_FORMAT(ppimgSrcRegions[1]->GetType()))
            {
                DO_IT( EL_FORMAT_BYTE, Byte, int );
                DO_IT( EL_FORMAT_SHORT, UInt16, int );
                DO_IT( EL_FORMAT_FLOAT, float, int );
                default:
                    return E_NOTIMPL;
            }
        }
    }

	HRESULT Clone(ITaskState **ppState)
	{
		return CloneTaskState<CIterationTransform>(ppState,
            VT_NOTHROWNEW CIterationTransform(m_cpl, m_iter, m_seed));
    }

private:
    Completion* m_cpl;
    int m_iter;
    int m_seed;
};

class CThreadAffinityWorkIdSequencer : public ITaskWorkIdSequencer
{
public:
	virtual int GetTotalWorkItems()
	{ return m_iTotalRows * m_iColsPerRow; }

	virtual HRESULT Advance(OUT bool& bDone, OUT int& iWorkId, int iThreadIndex)
	{
		bDone = false;
		if( iThreadIndex >= (int)m_vecColIds.size() )
		{
			bDone   = true;
			iWorkId = 0;
		}
		else
		{
			if( m_vecColIds[iThreadIndex] == 0 )
			{
				int iRowId;
                while (((iRowId = InterlockedIncrement(&m_iCurrentRowId)) & 1) == m_iOdd);
				if( iRowId > m_iTotalRows )
				{
					bDone = true;
				}
				m_vecRowIds[iThreadIndex] = (iRowId-1) * m_iColsPerRow;
			}
			iWorkId = m_vecRowIds[iThreadIndex] + (!m_bRev ? m_vecColIds[iThreadIndex] :
                (m_iColsPerRow - 1 - m_vecColIds[iThreadIndex]));

			if( m_vecColIds[iThreadIndex] >= m_iColsPerRow-1 )
			{
				m_vecColIds[iThreadIndex] = 0;
			}
			else
			{
				m_vecColIds[iThreadIndex]++;
			}
		}
		return S_OK;
	}

public:
	CThreadAffinityWorkIdSequencer() :
		m_iTotalRows(0), m_iColsPerRow(0), m_iOdd(0)
	{}

	HRESULT Initialize(int iCols, int iRows, int iThreads, bool bOdd, bool bRev)
	{
		VT_HR_BEGIN()

		VT_HR_EXIT( m_vecColIds.resize(iThreads) );
		VT_HR_EXIT( m_vecRowIds.resize(iThreads) );

		ZeroMemory(m_vecColIds.begin(), m_vecColIds.size()*sizeof(m_vecColIds[0]));
		
		m_iCurrentRowId = m_iOdd = bOdd ? 1 : 0;
		m_iTotalRows    = iRows;
		m_iColsPerRow   = iCols;
        m_bRev          = bRev;

		VT_HR_END()
	}

protected:
    bool m_bRev;
    int m_iOdd;
	int m_iTotalRows;
	int m_iColsPerRow;
	vt::vector<int> m_vecColIds;
	vt::vector<int> m_vecRowIds;
	volatile LONG   m_iCurrentRowId;
};

#else MULTI_CORE_ITERATION
//+-----------------------------------------------------------------------------
//
// Class: CSkipBlockIterator
// 
// Synposis: simple class to step block by block across a specified rectangular
//           region, odd or even rows, reverse or forward.
//
//------------------------------------------------------------------------------
class CSkipBlockIterator : public CBlockIterator
{
public:
    CSkipBlockIterator() : m_iSkip(0), m_iOdd(0), m_bRev(false),
        m_iBlkCntX(0), m_iBlkCntY(0)
    { m_bDone = false; }

	CSkipBlockIterator(const BLOCKITER_INIT& init, bool bSkip, bool bOdd, bool bRev)
    { Initialize(init, bSkip, bOdd, bRev); }

    CSkipBlockIterator(const vt::CRect& rctRegion, UInt32 uBlockSize, bool bSkip, bool bOdd, bool bRev)
    { Initialize(BLOCKITER_INIT(rctRegion, uBlockSize, uBlockSize), bSkip, bOdd, bRev); }

	int BlockCntX() const
	{ return m_iBlkCntX; }

	int BlockCntY() const
	{ return m_iBlkCntY; }

	bool Advance()
	{
		if( !m_bDone )
		{
            if( !m_bRev )
            {
			    m_ptCur.x += m_init.iBlockWidth;
			    if(	m_ptCur.x >= m_init.region.Width() )
			    {
				    m_ptCur.y += (m_init.iBlockHeight << m_iSkip);
				    if(	m_ptCur.y >= m_init.region.Height() )
				    {   
					    m_bDone = true; 
				    }
				    else
				    {
					    m_ptCur.x = 0;
				    }
			    }
            }
            else
            {
			    m_ptCur.x -= m_init.iBlockWidth;
			    if(	m_ptCur.x < 0 )
			    {
				    m_ptCur.y -= (m_init.iBlockHeight << m_iSkip);
				    if(	m_ptCur.y < 0 )
				    {   
					    m_bDone = true; 
				    }
				    else
				    {
					    m_ptCur.x = (m_iBlkCntX - 1) * m_init.iBlockWidth;
				    }
			    }
            }
		}
		return !m_bDone;
	}

protected:
	void Initialize(const BLOCKITER_INIT& init, bool bSkip, bool bOdd, bool bRev)
	{
        m_init  = init;
        m_iSkip = bSkip ? 1 : 0;
        m_iOdd  = bOdd  ? 1 : 0;
        m_bRev  = bRev;
        m_bDone = m_init.region.IsRectEmpty();

        m_iBlkCntX = (m_init.region.Width()  + m_init.iBlockWidth  - 1) / m_init.iBlockWidth;
        m_iBlkCntY = (m_init.region.Height() + m_init.iBlockHeight - 1) / m_init.iBlockHeight;

        m_ptCur = !m_bRev ? vt::CPoint(0,0) :
            vt::CPoint((m_iBlkCntX - 1) * m_init.iBlockWidth,
                       (m_iBlkCntY - 1) * m_init.iBlockHeight);

        if (m_iSkip)
        {
            if (!m_bRev)
            {
                m_ptCur.y += m_iOdd ? m_init.iBlockHeight : 0;
                m_bDone = m_bDone || m_ptCur.y >= m_init.region.Height();
            }
            else
            {
                m_ptCur.y -= (m_iBlkCntY & 1) == m_iOdd ? m_init.iBlockHeight : 0;
                m_bDone = m_bDone || m_ptCur.y < 0;
            }
        }
    }

    int  m_iSkip;
    int  m_iOdd;
	bool m_bRev;
	int  m_iBlkCntX;
	int  m_iBlkCntY;
};
#endif MULTI_CORE_ITERATION

HRESULT Completion::iteration(int iter)
{
    VT_HR_BEGIN()

    CImgInfo info = m_imgPyr->GetImgInfo(m_level);

    bool downPass = (iter & 1) == 0;

    int seed = m_rnd.IRand(INT_MAX);

    // Clear next pass bounding boxes before iteration.
    int phase = downPass ? 0 : 1;
    for (int i = 0; i < (int) m_bboxes[1 - phase].size(); i++)
        m_bboxes[1 - phase][i] = CRect(LONG_MAX, LONG_MAX, 0, 0);

#ifdef MULTI_CORE_ITERATION
    CIterationTransform transform(this, iter, seed);

    CLevelChangeWriter wr(&m_uvPyr, m_level);
    CLevelChangeReader rd(&m_uvPyr, m_level);
    CLevelChangeReader dst(m_imgPyr, m_level);
    CLevelChangeReader src(m_imgPyr, m_level);

    CTransformGraphNaryNode node(&transform);

    CUnknownBlockMap bm(info.Rect(), m_blockSize, m_blockSize, 1, 1);
    VT_HR_EXIT( bm.SetKnown(m_bboxes[phase]) );
    node.SetDest(NODE_DEST_PARAMS(&wr, &bm));

    VT_HR_EXIT( node.SetSourceCount(3) );
    VT_HR_EXIT( node.BindSourceToReader(0, &rd) );
    VT_HR_EXIT( node.BindSourceToReader(1, &dst) );
    VT_HR_EXIT( node.BindSourceToReader(2, &src) );

    VT_TRANSFORM_TASK_OPTIONS opts;

    // if user specified default then limit to numproc
	SYSTEM_INFO sysinfo;
	vt::CSystem::GetSystemInfo(&sysinfo);
	opts.maxthreads = sysinfo.dwNumberOfProcessors;

    CThreadAffinityWorkIdSequencer seq;
    opts.pSeq = &seq;

    CTaskProgress *prog = NULL;
    VT_HR_EXIT( seq.Initialize(bm.GetBlockCols(), bm.GetBlockRows(), opts.maxthreads, false, !downPass) );
    VT_HR_EXIT( PushTransformTaskAndWait(&node, prog, &opts) );
    if (m_blockSize < info.height)
    {
        VT_HR_EXIT( seq.Initialize(bm.GetBlockCols(), bm.GetBlockRows(), opts.maxthreads, true, !downPass) );
        VT_HR_EXIT( PushTransformTaskAndWait(&node, prog, &opts) );
    }

#else MULTI_CORE_ITERATION
    CImg uv;
    CImg imgDst, imgSrc;

#ifdef DEBUG_OUTPUT
    FILE * flog = 0;

    #define DO_IT(EL, TRGBA, TUV) \
        case EL: \
            VT_HR_EXIT( iterateBlock(iter, seed, rctCur, \
                                     rctUv, (CTypedImg<TUV> &) uv,\
                                     rctDst, (CCompositeImg<RGBAType<TRGBA>> &) imgDst, \
                                     rctSrc, (CCompositeImg<RGBAType<TRGBA>> &) imgSrc, \
                                     flog) ); \
            break;
#else DEBUG_OUTPUT
    #define DO_IT(EL, TRGBA, TUV) \
        case EL: \
            VT_HR_EXIT( iterateBlock(iter, seed, rctCur, \
                                     rctUv, (CTypedImg<TUV> &) uv,\
                                     rctDst, (CCompositeImg<RGBAType<TRGBA>> &) imgDst, \
                                     rctSrc, (CCompositeImg<RGBAType<TRGBA>> &) imgSrc) ); \
            break;
#endif DEBUG_OUTPUT

    int xBlocks = (info.width  + m_blockSize - 1) / m_blockSize;

#ifdef MCORE_MATCH
    for (int i = 0; i < 2; i++)
    {
    for (CSkipBlockIterator bi(info.Rect(), m_blockSize, true, i == 1, !downPass); !bi.Done(); bi.Advance())
#else MCORE_MATCH
#ifdef DEBUG_OUTPUT
    if (m_dbgSaveIterationData)
    {
        wstring fileName = wformat(L"output\\d_itlog_%04d_ref.txt", m_dbgGen);
        _wfopen_s(&flog, fileName.get_constbuffer(), L"wb");

        fileName = wformat(L"output\\d_uv_%04d_ref.vti", m_dbgGen++);
        saveUv(fileName);
    }
#endif DEBUG_OUTPUT

    for (CSkipBlockIterator bi(info.Rect(), m_blockSize, false, false, !downPass); !bi.Done(); bi.Advance())
#endif MCORE_MATCH
    {
        CRect rctCur = bi.GetRect();

#ifdef MCORE_MATCH
#ifdef DEBUG_OUTPUT
        if (m_dbgSaveIterationData && rctCur.left == (downPass ? 0 : (bi.BlockCntX() - 1) * m_blockSize))
        {
            if (flog)
                fclose(flog);

            wstring fileName = wformat(L"output\\d_itlog_%04d_ref.txt", m_dbgGen);
            _wfopen_s(&flog, fileName.get_constbuffer(), L"wb");

            fileName = wformat(L"output\\d_uv_%04d_ref.vti", m_dbgGen++);
            saveUv(fileName);
        }
#endif DEBUG_OUTPUT
#endif MCORE_MATCH

        int block = (rctCur.left / m_blockSize) + (rctCur.top / m_blockSize) * xBlocks;

        // find source bounding box
#ifdef USE_TILE_BBOXES
        CRect rctBbx = m_tileUnions[block];
        rctBbx = CRect(rctBbx.left   << (m_tileLevel + Tiles::levelOffset - m_level),
                       rctBbx.top    << (m_tileLevel + Tiles::levelOffset - m_level),
                       rctBbx.right  << (m_tileLevel + Tiles::levelOffset - m_level),
                       rctBbx.bottom << (m_tileLevel + Tiles::levelOffset - m_level));
#else USE_TILE_BBOXES
        CRect rctBbx = m_bboxes[0][block] | m_bboxes[1][block];
#endif USE_TILE_BBOXES
        if (rctBbx.IsRectEmpty())
            continue;

        if (downPass)
            rctBbx.InflateRect(0, 0, 1, 1);
        else
            rctBbx.InflateRect(1, 1, 0, 0);
        rctBbx &= info.Rect();

        CRect rctSrc = rctBbx;
        rctSrc.InflateRect(s_prad, s_prad);
        rctSrc &= info.Rect();

        // read source pixels
        VT_HR_EXIT( m_imgPyr->ReadRegion(rctSrc, imgSrc, m_level) );

        CRect rctUv;
        if (iter == s_numIters - 1)
        {
            // Time for sparse average; will need bounding box of marks.
            CRect rctMrk = CRect((rctCur.left - s_prad) >> 2,
                                 (rctCur.top  - s_prad) >> 2,
                                 ((rctCur.right  + s_prad) >> 2) + 1,
                                 ((rctCur.bottom + s_prad) >> 2) + 1);
            rctUv = CRect(rctMrk.left << 2, rctMrk.top << 2,
                          rctMrk.right << 2, rctMrk.bottom << 2);
            rctUv.InflateRect(s_prad, s_prad);
        }
        else
        {
            rctUv = rctCur;
            rctUv.InflateRect(1, 1);
        }
        rctUv &= info.Rect();
        VT_HR_EXIT( m_uvPyr.ReadRegion(rctUv, uv, m_level) );

        CRect rctDst = rctCur;
        rctDst.InflateRect(s_prad, s_prad);
        rctDst &= info.Rect();
        VT_HR_EXIT( m_imgPyr->ReadRegion(rctDst, imgDst, m_level) );

        if (EL_FORMAT(m_uvPyr.GetImgInfo().type) == EL_FORMAT_SHORT)
        {
            switch(EL_FORMAT(info.type))
            {
            DO_IT( EL_FORMAT_BYTE, Byte, UInt16 );
            DO_IT( EL_FORMAT_SHORT, UInt16, UInt16 );
            DO_IT( EL_FORMAT_FLOAT, float, UInt16 );
            default:
                VT_HR_EXIT( E_NOTIMPL );
                break;
            }
        }
        else
        {
            switch(EL_FORMAT(info.type))
            {
            DO_IT( EL_FORMAT_BYTE, Byte, int );
            DO_IT( EL_FORMAT_SHORT, UInt16, int );
            DO_IT( EL_FORMAT_FLOAT, float, int );
            default:
                VT_HR_EXIT( E_NOTIMPL );
                break;
            }
        }

        CImg uvCur;
        CRect rctUvCur = rctCur - rctUv.TopLeft();
        VT_HR_EXIT( uv.Share(uvCur, &rctUvCur) );
        VT_HR_EXIT( m_uvPyr.WriteRegion(rctCur, uvCur, m_level) );
    }
#ifdef MCORE_MATCH
    }
#endif MCORE_MATCH

#ifdef DEBUG_OUTPUT
    if (flog)
    {
        fclose(flog);
    }
#endif DEBUG_OUTPUT

#endif MULTI_CORE_ITERATION

    VT_HR_END()
}
