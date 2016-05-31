#include "stdafx.h"

#include "vtfileio.h"

#include "vt_completion.h"
#include "completion.h"
#include "completion_tiles.h"

using namespace vt;

#include "DistanceTransform.h"

HRESULT vt::VtCompleteImage(
    IImageReaderWriter * pPyramid,
    CTaskProgress* pProgress,
    FusionMethod /*fusionMethod*/,
    int numRepetitions, int numRepetitionsDec,
    float discontinuousPenalty, float windowEnergyThresh,
    const wchar_t * timingFile, CImg * pDebugUV)
{
    VT_HR_BEGIN()

    Completion cpl;

    if (pPyramid == NULL || numRepetitions <= 0 || numRepetitionsDec <= 0 ||
        discontinuousPenalty < 0.f || windowEnergyThresh < 0.f)
        VT_HR_EXIT( E_INVALIDARG );

    // set up progress reporting
    CPhasedTaskStatus report;
    report.SetOuterCallback(pProgress);
    float fRemainingProgress = 100.f;

    report.BeginPhase("Initializing", 10.f);
    VT_HR_EXIT( cpl.init(pPyramid, timingFile, discontinuousPenalty, windowEnergyThresh) );
    report.ReportProgress(100.f);
    fRemainingProgress -= 10.f;

    int startLevel = cpl.m_level;
    float fProgressL = 90.f / (float) (startLevel + 1);

    bool bDone = false;
    while (!bDone)
    {
        string phaseName;
		VT_HR_EXIT( phaseName.format_with_resize("Completing pyramid level %u", cpl.m_level) );
		report.BeginPhase(phaseName.get_constbuffer(), fProgressL);

        if (cpl.m_level == startLevel)
        {
            VT_HR_EXIT( cpl.pass(numRepetitions, true) );
            VT_HR_EXIT( cpl.pass(numRepetitions, false) );
        }
        else
        {
            VT_HR_EXIT( cpl.pass(numRepetitions, false) );
        }

        bDone = cpl.m_level == 0;
        if (!bDone)
        {
            if (EL_FORMAT(cpl.m_uvPyr.GetImgInfo().type) == EL_FORMAT_INT)
                VT_HR_EXIT( cpl.upsample<int>() );
            else
                VT_HR_EXIT( cpl.upsample<UInt16>() );
        }

        numRepetitions = VtMax(1, numRepetitions / numRepetitionsDec);

        report.ReportProgress(100.f);
    }

    VT_HR_EXIT( cpl.finish() );

    if (pDebugUV != NULL)
        VT_HR_EXIT( cpl.m_uvPyr.ReadImg(*pDebugUV) );

    VT_HR_END()
}

wstring vt::wformat(const wchar_t * str, ...)
{
    va_list argptr;
    va_start(argptr, str);
    int length = _vscwprintf(str, argptr) + 1;
    wstring res;
    res.resize(length);
    vswprintf_s(res.get_buffer(), length, str, argptr);
    return res;
};

Completion::Completion()
{
#ifdef DEBUG_OUTPUT
    m_dbgSaveIterationData = true;
    m_dbgSaveAverageData = true;
    m_dbgGen = 0;
#endif DEBUG_OUTPUT

    m_timePass = m_timePassIteration = m_timePassAverage = m_timeUpsample = m_timeUpsampleAverage = m_timeUpsampleInvalid = m_timeFusion = 0;

    m_timeCompletion = m_timer.GetTimeMilliSec() * 0.001f;
}

HRESULT Completion::init(IImageReaderWriter * pPyr, const wchar_t * pwszTimingFile,
    float discontinuousPenalty, float windowEnergyThresh)
{
    VT_HR_BEGIN()

    if (pwszTimingFile != NULL)
        VT_HR_EXIT( m_csvFileName.assign(pwszTimingFile) );

    m_weightGamma = 2.f;
    m_discontinuousPenalty = discontinuousPenalty;
    m_windowEnergyThresh = windowEnergyThresh;

    VT_DEBUG_LOG( "compiler switches:\n" );
#ifdef DEBUG_OUTPUT
    VT_DEBUG_LOG( "  DEBUG_OUTPUT:           enabled\n" );
#else
    VT_DEBUG_LOG( "  DEBUG_OUTPUT:           disabled\n" );
#endif

#ifdef DEBUG_MATCH
    VT_DEBUG_LOG( "  DEBUG_MATCH:            enabled\n" );
#else
    VT_DEBUG_LOG( "  DEBUG_MATCH:            disabled\n" );
#endif

#ifdef MCORE_MATCH
    VT_DEBUG_LOG( "  MCORE_MATCH:            enabled\n" );
#else
    VT_DEBUG_LOG( "  MCORE_MATCH:            disabled\n" );
#endif

#ifdef MULTI_CORE_INIT
    VT_DEBUG_LOG( "  MULTI_CORE_INIT:        enabled\n" );
#else
    VT_DEBUG_LOG( "  MULTI_CORE_INIT:        disabled\n" );
#endif

#ifdef MULTI_CORE_ITERATION
    VT_DEBUG_LOG( "  MULTI_CORE_ITERATION:   enabled\n" );
#else
    VT_DEBUG_LOG( "  MULTI_CORE_ITERATION:   disabled\n" );
#endif

#ifdef MULTI_CORE_AVG
    VT_DEBUG_LOG( "  MULTI_CORE_AVG:         enabled\n" );
#else
    VT_DEBUG_LOG( "  MULTI_CORE_AVG:         disabled\n" );
#endif

    VT_DEBUG_LOG( "\n" );

    VT_DEBUG_LOG( "init...\n" );

    m_timeInit = m_timer.GetTimeMilliSec() * 0.001f;

    if (pPyr->GetImgInfo().Bands() != 4)
    {
        VT_DEBUG_LOG( "only RGBA images supported." );
        VT_HR_EXIT( E_INVALIDARG );
    }
    
    m_imgPyr = pPyr;

    // find coarsest (start) level
    int startLevel = 0;
    CImgInfo info = pPyr->GetImgInfo(startLevel);
    while (info.width > 128 || info.height > 128)
    {
        startLevel++;
        info = pPyr->GetImgInfo(startLevel);
    }

    VT_DEBUG_LOG( "  start level %d\n", startLevel );

    // init boundary distance
    VT_HR_EXIT( initBoundaryDistance() );

    // process tiles
    m_timeInitTiles = m_timer.GetTimeMilliSec() * 0.001f;

    Tiles tiles;
    VT_HR_EXIT( tiles.process(pPyr, startLevel, m_tileLevel, m_restriction, m_tiles) );

    m_timeInitTiles = m_timer.GetTimeMilliSec() * 0.001f - m_timeInitTiles;

    IMG_CACHE_SOURCE_PROPERTIES props;
    props.bAutoLevelGenerate = false;
    props.iMaxLevel = startLevel;

    // init start level
    m_level = startLevel;
    info = m_imgPyr->GetImgInfo(2);
    VT_HR_EXIT( m_maskPyr.Create(info.width, info.height, 1, &props) );
    info = m_imgPyr->GetImgInfo();

    if (info.width  <= ElTraits<UInt16>::MaxVal() &&
        info.height <= ElTraits<UInt16>::MaxVal())
    {
        VT_HR_EXIT( m_uvPyr.Create(CImgInfo(info.width, info.height,
                                            VT_IMG_MAKE_TYPE(EL_FORMAT_SHORT, 2)),
                                   &props) );

	    vector<CVec2<UInt16>> validPixels, holePixels, uvPixels;
        VT_HR_EXIT( initLevel(&validPixels, &holePixels, &uvPixels) );

        VT_HR_EXIT( randomize(validPixels, holePixels, uvPixels) );
    }
    else
    {
        VT_HR_EXIT( m_uvPyr.Create(CImgInfo(info.width, info.height,
                                            VT_IMG_MAKE_TYPE(EL_FORMAT_INT, 2)),
                                   &props) );

	    vector<CVec2<int>> validPixels, holePixels, uvPixels;
        VT_HR_EXIT( initLevel(&validPixels, &holePixels, &uvPixels) );

        VT_HR_EXIT( randomize(validPixels, holePixels, uvPixels) );
    }

    m_timeInit = m_timer.GetTimeMilliSec() * 0.001f - m_timeInit;

    VT_HR_END()
}

HRESULT Completion::initBoundaryDistance()
{
    VT_HR_BEGIN()

    m_bdDistLevel = 0;
    CImgInfo info = m_imgPyr->GetImgInfo(m_bdDistLevel);
    while (info.width > 1024 || info.height > 1024)
    {
        m_bdDistLevel++;
        info = m_imgPyr->GetImgInfo(m_bdDistLevel);
    }
    
    VT_DEBUG_LOG( "  computing boundary distance at level %d\n", m_bdDistLevel );

    CByteImg bd;
    VT_HR_EXIT( m_imgPyr->ReadImg(bd, m_bdDistLevel) );

    const int w = bd.Width(), h = bd.Height();

    VT_HR_EXIT( m_bdDist.Create(w, h) );
    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            m_bdDist(x, y) = (bd(x, y, 3) == 255 ? 0 : FLT_MAX);
        }
    }

    VT_HR_EXIT( distanceTransform(m_bdDist) );

    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            m_bdDist(x, y) = 10.f * sqrtf(m_bdDist(x, y));
        }
    }

    // compute weights
    m_bdDistWeight.clear();
    for (int i = 0; i < s_bdDistWeightLength; i++)
    {
        float dist = i / 10.0f;
        float weight = powf(m_weightGamma, -dist);
        VT_HR_EXIT( m_bdDistWeight.push_back(weight) );
    }

    VT_HR_END()
}

HRESULT Completion::finish()
{
    VT_HR_BEGIN()

    // Set the alpha to opaque over the whole image.
    CImg img;
    Byte val[4];
    switch (EL_FORMAT(img.GetType()))
    {
    case EL_FORMAT_BYTE:  *((Byte *)   val) = ElTraits<Byte>::MaxVal();   break;
    case EL_FORMAT_SHORT: *((UInt16 *) val) = ElTraits<UInt16>::MaxVal(); break;
    case EL_FORMAT_FLOAT: *((float *)  val) = ElTraits<float>::MaxVal();  break;
    default: VT_HR_EXIT( E_NOTIMPL );
    }
     
    for (CBlockIterator bi(m_imgPyr->GetImgInfo().Rect(), m_blockSize); !bi.Done(); bi.Advance())
    {
        CRect rctCur = bi.GetRect();

        VT_HR_EXIT( m_imgPyr->ReadRegion(rctCur, img) );

        VT_HR_EXIT( img.Fill(val, NULL, 3) );

        VT_HR_EXIT( m_imgPyr->WriteRegion(rctCur, img) );
    }

    m_timeCompletion = m_timer.GetTimeMilliSec() * 0.001f - m_timeCompletion;

    VT_DEBUG_LOG( "\n" );

    VT_DEBUG_LOG( "TIME Completion:       %f\n", m_timeCompletion );
    VT_DEBUG_LOG( "TIME   Init:           %f\n", m_timeInit );
    VT_DEBUG_LOG( "TIME     Tiles:        %f\n", m_timeInitTiles );
    VT_DEBUG_LOG( "TIME     Other:        %f\n", m_timeInit - m_timeInitTiles );
    VT_DEBUG_LOG( "TIME   Pass:           %f\n", m_timePass );
    VT_DEBUG_LOG( "TIME     Iteration:    %f\n", m_timePassIteration );
    VT_DEBUG_LOG( "TIME     Average:      %f\n", m_timePassAverage );
    VT_DEBUG_LOG( "TIME     Other:        %f\n", m_timePass - (m_timePassIteration + m_timePassAverage) );
    VT_DEBUG_LOG( "TIME   Upsample:       %f\n", m_timeUpsample );
    VT_DEBUG_LOG( "TIME     Average:      %f\n", m_timeUpsampleAverage );
    VT_DEBUG_LOG( "TIME     Invalid:      %f\n", m_timeUpsampleInvalid );
    VT_DEBUG_LOG( "TIME     Other:        %f\n", m_timeUpsample - (m_timeUpsampleAverage + m_timeUpsampleInvalid) );
    VT_DEBUG_LOG( "TIME   Fusion:         %f\n", m_timeFusion );
    VT_DEBUG_LOG( "TIME   Other:          %f\n", m_timeCompletion - (m_timeInit + m_timePass + m_timeUpsample + m_timeFusion) );

    if (!m_csvFileName.empty())
    {
        FILE * fout = NULL;
        _wfopen_s(&fout, m_csvFileName.get_constbuffer(), L"ab");
        if (fout == NULL)
            VT_HR_EXIT( E_WRITEFAILED );

        fprintf(fout, "Completion,%f\n", m_timeCompletion);
        fprintf(fout, "  Init,%f\n", m_timeInit);
        fprintf(fout, "    Tiles,%f\n", m_timeInitTiles);
        fprintf(fout, "    Other,%f\n", m_timeInit - m_timeInitTiles);
        fprintf(fout, "  Pass,%f\n", m_timePass);
        fprintf(fout, "    Iteration,%f\n", m_timePassIteration);
        fprintf(fout, "    Average,%f\n", m_timePassAverage);
        fprintf(fout, "    Other,%f\n", m_timePass - (m_timePassIteration + m_timePassAverage));
        fprintf(fout, "  Upsample,%f\n", m_timeUpsample);
        fprintf(fout, "    Average,%f\n", m_timeUpsampleAverage);
        fprintf(fout, "    Invalid,%f\n", m_timeUpsampleInvalid);
        fprintf(fout, "    Other,%f\n", m_timeUpsample - (m_timeUpsampleAverage + m_timeUpsampleInvalid));
        fprintf(fout, "  Fusion,%f\n", m_timeFusion);
        fprintf(fout, "  Other,%f\n", m_timeCompletion - (m_timeInit + m_timePass + m_timeUpsample + m_timeFusion));

        fclose(fout);
    }

    VT_HR_END()
}

template <typename TUV>
HRESULT Completion::initLevel(vector<CVec2<TUV>>* validPixels, vector<CVec2<TUV>>* holePixels, vector<CVec2<TUV>>* uvPixels)
{
#ifdef MULTI_CORE_INIT
    Byte * row = NULL;
#endif MULTI_CORE_INIT

    VT_HR_BEGIN()

#ifdef DEBUG_OUTPUT
    float startTime, endTime;

    startTime = m_timer.GetTimeMilliSec() * 0.001f;
#endif DEBUG_OUTPUT

    m_rnd.Seed(0);

    CImgInfo info = m_imgPyr->GetImgInfo(m_level);

    int w = info.width;
    int h = info.height;

    m_blockSize = m_level <= m_tileLevel ? 256 : VtMax(w, h);

    if (holePixels != NULL)
        holePixels->clear();
    if (uvPixels != NULL)
        uvPixels->clear();
    if (validPixels != NULL)
        validPixels->clear();

#ifdef MULTI_CORE_INIT
    Byte parea = s_psize * s_psize;

    row = VT_NOTHROWNEW Byte[info.width + s_prad];
    VT_PTR_OOM_EXIT( row );
#endif MULTI_CORE_INIT

    CByteImg imgBlk;
    for (CBlockIterator bi(info.Rect(), m_blockSize); !bi.Done(); bi.Advance())
    {
        CRect rctDst = bi.GetRect();

        CRect rctSrc = rctDst;
        rctSrc.InflateRect(s_prad, s_prad);
        rctSrc &= info.Rect();
        VT_HR_EXIT( m_imgPyr->ReadRegion(rctSrc, imgBlk, m_level) );

#ifdef MULTI_CORE_INIT
        int stride = imgBlk.StrideBytes();

        for (int x = rctSrc.left; x < rctSrc.right; x++)
        {
            Byte sum = 0;
            Byte * colPtr = imgBlk.BytePtr(x - rctSrc.left, 0, 3);

            // pre-sum
            Byte * addPtr = colPtr;
            for (int y = 0; y < VtMin(s_prad, rctSrc.Height()); y++)
            {
                if ((rctDst.left == 0 || x >= rctDst.left + s_prad) && rctDst.top == 0)
                {
                    if (*addPtr == 255)
                        *addPtr = s_knownBit;
                    else
                        *addPtr = 0;
                }
                if (*addPtr & s_knownBit)
                    sum++;
                addPtr += stride;
            }

            // add + mask
            Byte * valPtr = colPtr;
            for (int y = 0; y <= VtMin(s_prad, rctSrc.Height() - 1 - s_prad); y++)
            {
                if ((rctDst.left == 0 || x >= rctDst.left + s_prad) &&
                    (rctDst.top == 0 || y == s_prad))
                {
                    if (*addPtr == 255)
                        *addPtr = s_knownBit;
                    else
                        *addPtr = 0;
                }
                if (*addPtr & s_knownBit)
                    sum++;

                *valPtr |= sum;

                addPtr += stride;
                valPtr += stride;
            }

            // add + del + mask
            Byte * delPtr = colPtr;
            for (int y = s_prad + 1; y < rctSrc.Height() - s_prad; y++)
            {
                if ((rctDst.left == 0 || x >= rctDst.left + s_prad))
                {
                    if (*addPtr == 255)
                        *addPtr = s_knownBit;
                    else
                        *addPtr = 0;
                }
                if (*addPtr & s_knownBit)
                    sum++;

                if ((*delPtr & s_knownBit) != 0)
                {
                    sum--;
                }

                *valPtr |= sum;

                addPtr += stride;
                delPtr += stride;
                valPtr += stride;
            }
        }

        for (int y = rctDst.top; y < rctDst.bottom; y++)
        {
            Byte * ptr = imgBlk.BytePtr(0, y - rctSrc.top, 3);
            Byte * copyPtr = &row[0];
            Byte * delPtr = copyPtr - s_prad - 1;
            Byte sum = 0;

            for (int x = 0; x < s_prad; x++)
            {
                sum += (*(ptr + (x << 2)) & 63);
            }

            for (int x = rctSrc.left; x < rctDst.right; x++)
            {
                Byte val = *ptr;
                *copyPtr = val;

                *ptr &= (s_knownBit | s_validBit);
                if (x >= rctDst.left)
                {
                    if ((val & s_knownBit) == 0)
                    {
                        if (holePixels != NULL)
                            VT_HR_EXIT( holePixels->push_back(CVec2<TUV>((TUV) x, (TUV) y)) );
                    }
                }

                int nx = x + s_prad, px = x - s_prad - 1;
                if (nx < w)
                {
                    sum += (*(ptr + (s_prad << 2)) & 63);
                }
                if (px >= rctSrc.left)
                {
                    sum -= (*delPtr & 63);
                }

                if (sum == parea)
                {
                    *ptr |= s_validBit;
                    if (validPixels != NULL)
                        VT_HR_EXIT( validPixels->push_back(CVec2<TUV>((TUV) x, (TUV) y)) );
                }
                else if (x >= rctSrc.left + s_prad && x < w - s_prad && y >= s_prad && y < h - s_prad)
                {
                    if (uvPixels != NULL)
                        VT_HR_EXIT( uvPixels->push_back(CVec2<TUV>((TUV) x, (TUV) y)) );
                }

                ptr += 4;
                copyPtr++;
                delPtr++;
            }
        }

        VT_HR_EXIT( m_imgPyr->WriteRegion(rctSrc, imgBlk, m_level) );

#else // MULTI_CORE_INIT
        for (int y = rctDst.top; y < rctDst.bottom; y++)
        {
            for (int x = rctDst.left; x < rctDst.right; x++)
            {
                // known or unknown?
                if (imgBlk(x - rctSrc.left, y - rctSrc.top, 3) < 255)
                {
                    if (holePixels != NULL)
                        VT_HR_EXIT( holePixels->push_back(CVec2<TUV>((TUV) x, (TUV) y)) );
                    imgBlk(x - rctSrc.left, y - rctSrc.top, 3) = 0;
                }
                else
                {
                    imgBlk(x - rctSrc.left, y - rctSrc.top, 3) = s_knownBit;
                }
            }
        }

        for (int y = VtMax((int) rctDst.top, s_prad); y < VtMin((int) rctDst.bottom, h - s_prad); y++)
        {
            for (int x = VtMax((int) rctDst.left, s_prad); x < VtMin((int) rctDst.right, w - s_prad); x++)
            {
                // valid source patch?
                bool hasHole = false;
                for (int ry = y - s_prad; ry <= y + s_prad && !hasHole; ry++)
                {
                    for (int rx = x - s_prad; rx <= x + s_prad && !hasHole; rx++)
                    {
                        if (ry < rctDst.top || (rx < rctDst.right && ry < rctDst.bottom))
                            hasHole = (imgBlk(rx - rctSrc.left, ry - rctSrc.top, 3) & s_knownBit) == 0;
                        else
                            hasHole = imgBlk(rx - rctSrc.left, ry - rctSrc.top, 3) < 255;
                    }
                }

                if (hasHole)
                {
                    if (uvPixels != NULL)
                        VT_HR_EXIT( uvPixels->push_back(CVec2<TUV>((TUV) x, (TUV) y)) );
                }
                else
                {
                    imgBlk(x - rctSrc.left, y - rctSrc.top, 3) |= s_validBit;
                    if (validPixels != NULL)
                        VT_HR_EXIT( validPixels->push_back(CVec2<TUV>((TUV) x, (TUV) y)) );
                }
            }
        }

        VT_HR_EXIT( m_imgPyr->WriteRegion(rctSrc, imgBlk, m_level) );
#endif MULTI_CORE_INIT
    }

#ifdef DEBUG_OUTPUT 
    endTime = m_timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "    INIT-LVL --- holepixels: %f\n", endTime - startTime );

    startTime = endTime;
#endif DEBUG_OUTPUT

#ifdef DEBUG_OUTPUT
    endTime = m_timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "    INIT-LVL --- uvpixels: %f\n", endTime - startTime );

    startTime = endTime;
#endif DEBUG_OUTPUT

#ifdef DEBUG_OUTPUT
    endTime = m_timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "    INIT-LVL --- mcore: %f\n", endTime - startTime );

    {
        wstring fileName = wformat(L"output\\d_initlevel_%02d.png", m_level);
        CByteImg mask;
        VT_HR_EXIT( m_maskPyr.ReadImg(mask, m_level) );
        VT_HR_EXIT( VtSaveImage(fileName, mask) );
    }
#endif DEBUG_OUTPUT

    VT_HR_EXIT( m_maskPyr.Fill(NULL, NULL, m_level) );
    m_mark = 0;

    VT_HR_EXIT_LABEL()

#ifdef MULTI_CORE_INIT
    delete[] row;
#endif MULTI_CORE_INIT

    return hr;
}

template <typename TUV>
HRESULT Completion::randomize(vector<CVec2<TUV>>& validPixels, vector<CVec2<TUV>>& holePixels, vector<CVec2<TUV>>& uvPixels)
{
    VT_HR_BEGIN()

    CByteImg img;
    VT_HR_EXIT( m_imgPyr->ReadImg(img, m_level) );

    int w = img.Width();
    int h = img.Height();

    CTypedImg<TUV> uv;
    VT_HR_EXIT( uv.Create(w, h, 2) );
    VT_HR_EXIT( uv.Fill(ElTraits<TUV>::MaxVal()) );

    int tileShift = m_tileLevel + Tiles::levelOffset - m_level;
    int tileWidth = m_restriction[0].Width();
    int tileHeight = m_restriction[0].Height();

    for (int i = 0; i < (int) uvPixels.size(); i++)
    {
        const CVec2<TUV> & t = uvPixels[i];
        CVec2<TUV> rs;
        int tile = tileWidth * (t.y >> tileShift) + (t.x >> tileShift);

        int maxTries = 100000, curTry = 0;

        for (curTry; curTry < maxTries; curTry++)
        {
            int tx = m_rnd.IRand(tileWidth);
            int ty = m_rnd.IRand(tileHeight);

            if (m_restriction[tile](tx, ty) < FLT_MAX)
            {
                int tx0 = tx << tileShift, ty0 = ty << tileShift;
                int tx1 = (tx+1) << tileShift, ty1 = (ty+1) << tileShift;
                tx1 = VtMin(tx1, w); ty1 = VtMin(ty1, h);

                rs.x = (TUV) (tx0 + m_rnd.IRand(tx1 - tx0));
                rs.y = (TUV) (ty0 + m_rnd.IRand(ty1 - ty0));

                if (img(rs.x, rs.y, 3) & s_validBit)
                {
                    break;
                }
            }
        }
        if (curTry == maxTries)
        {
            rs = (!validPixels.empty() ? validPixels[m_rnd.IRand((int) validPixels.size())] :
                CVec2<TUV>(ElTraits<TUV>::MaxVal(), ElTraits<TUV>::MaxVal()));
            VT_DEBUG_LOG( "failed to initalize pixel %d, %d; using %d, %d\n", t.x, t.y, rs.x, rs.y );
        }

        uv(t.x, t.y, 0) = rs.x;
        uv(t.x, t.y, 1) = rs.y;
    }

    VT_HR_EXIT( m_uvPyr.WriteImg(uv, m_level) );

    VT_HR_EXIT( m_bboxes[0].push_back(uv.Rect()) );
    VT_HR_EXIT( m_bboxes[1].push_back(CRect(LONG_MAX, LONG_MAX, 0, 0)) );
    VT_HR_EXIT( m_bboxes[2].push_back(CRect(LONG_MAX, LONG_MAX, 0, 0)) );
#ifdef USE_TILE_BBOXES
    VT_HR_EXIT( m_tileUnions.push_back(uv.Rect()) );
#endif USE_TILE_BBOXES

    CByteImg img2;
    VT_HR_EXIT( img.CopyTo(img2) );

    // init with wave diffusion
    CByteImg visited(w, h, 1);
    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            visited(x, y) = (img(x, y, 3) & s_knownBit) == 0 ? 0 : 255;
        }
    }

    CByteImg visited2;
    VT_HR_EXIT( visited.CopyTo(visited2) );

    bool foundPixels;
    do
    {
        foundPixels = false;

        for (int i = 0; i < (int) holePixels.size(); i++)
        {
            int x = holePixels[i].x, y = holePixels[i].y;
            if (visited(x, y))
            {
                continue;
            }

            int col[3] = {0, 0, 0};
            int nn = 0;

            if (x > 0 && visited(x-1, y))
            {
                for (int ch = 0; ch < 3; ch++)
                {
                    col[ch] += img(x-1, y, ch);
                }
                nn++;
            }

            if (x < w - 1 && visited(x + 1, y))
            {
                for (int ch = 0; ch < 3; ch++)
                {
                    col[ch] += img(x+1, y, ch);
                }
                nn++;
            }

            if (y > 0 && visited(x, y - 1))
            {
                for (int ch = 0; ch < 3; ch++)
                {
                    col[ch] += img(x, y-1, ch);
                }
                nn++;
            }

            if (y < h - 1 && visited(x, y + 1))
            {
                for (int ch = 0; ch < 3; ch++)
                {
                    col[ch] += img(x, y+1, ch);
                }
                nn++;
            }

            if (nn == 0)
            {
                continue;
            }

            foundPixels = true;
            visited2(x, y) = 255;

            for (int ch = 0; ch < 3; ch++)
            {
                img2(x, y, ch) = (Byte) (col[ch] / nn);
            }
        }

        if (!foundPixels)
        {
            break;
        }

        VT_HR_EXIT( img2.CopyTo(img) );
        VT_HR_EXIT( visited2.CopyTo(visited) );
    }
    while (foundPixels);

    // diffusion
    CFloatImg diff[2];
    VT_HR_EXIT( VtConvertImage(diff[0], img) );
    VT_HR_EXIT( diff[0].CopyTo(diff[1]) );

    for (int step = 0; step < 1000; step++)
    {
        for (int i = 0; i < (int) holePixels.size(); i++)
        {
            int x = holePixels[i].x, y = holePixels[i].y;

            float col[3] = {0, 0, 0};
            int nn = 0;

            if (x > 0)
            {
                nn++;
                for (int ch = 0; ch < 3; ch++)
                {
                    col[ch] += diff[step & 1](x-1, y, ch);
                }
            }
            if (x < w - 1)
            {
                nn++;
                for (int ch = 0; ch < 3; ch++)
                {
                    col[ch] += diff[step & 1](x+1, y, ch);
                }
            }
            if (y > 0)
            {
                nn++;
                for (int ch = 0; ch < 3; ch++)
                {
                    col[ch] += diff[step & 1](x, y-1, ch);
                }
            }
            if (y < h - 1)
            {
                nn++;
                for (int ch = 0; ch < 3; ch++)
                {
                    col[ch] += diff[step & 1](x, y+1, ch);
                }
            }

            for (int ch = 0; ch < 3; ch++)
            {
                diff[1 - (step & 1)](x, y, ch) = col[ch] / nn;
            }
        }
    }

    VT_HR_EXIT( VtConvertImage(img, diff[0]) );

    VT_HR_EXIT( m_imgPyr->WriteImg(img, m_level) );

    VT_HR_END()
}

HRESULT Completion::pass(const int numRepetitions, const bool lockImage)
{
    VT_HR_BEGIN()

    VT_DEBUG_LOG( "pass...   (lvl: %d, #: %d, lock image: %d)\n", m_level, numRepetitions, lockImage );

    float thisTimePass = m_timer.GetTimeMilliSec() * 0.001f;

#if defined(_DEBUG) || defined(_HRPRINT)
    float startTime = thisTimePass;
#endif

    int rep;
    for (rep = 1; rep <= numRepetitions; rep++)
    {
        if (!lockImage)
        {
            m_mark += 3;
            if (m_mark >= 253)
            {
                m_mark = 1;
                VT_HR_EXIT( m_maskPyr.Fill(NULL, NULL, m_level) );
            }
        }

        float thisTimePassIteration = m_timer.GetTimeMilliSec() * 0.001f;

        m_changed = false;

        for (int iter = 0; iter < s_numIters; iter++)
        {
            VT_HR_EXIT( iteration(iter) );
        }

        m_timePassIteration += m_timer.GetTimeMilliSec() * 0.001f - thisTimePassIteration;

        if (!m_changed)
            break;

        if (!lockImage)
        {
            float thisTimePassAverage = m_timer.GetTimeMilliSec() * 0.001f;
            VT_ASSERT( m_mark != 0);
            VT_HR_EXIT( average() );
            m_timePassAverage += m_timer.GetTimeMilliSec() * 0.001f - thisTimePassAverage;
        }
    }

    if (lockImage)
    {
        float thisTimePassAverage = m_timer.GetTimeMilliSec() * 0.001f;
        VT_ASSERT( m_mark == 0);
        VT_HR_EXIT( average() );
        m_timePassAverage += m_timer.GetTimeMilliSec() * 0.001f - thisTimePassAverage;
    }

#ifdef DEBUG_OUTPUT
    wstring fileName = wformat(L"output\\d_lvl_%02d_skip.vti", m_level);
    saveUv(fileName);
#endif DEBUG_OUTPUT

    m_timePass += m_timer.GetTimeMilliSec() * 0.001f - thisTimePass;

    VT_DEBUG_LOG( "  [done] level %d %d reps took %.2fs\n", m_level, rep - 1,
        m_timer.GetTimeMilliSec() * 0.001f - startTime );

#ifdef DEBUG_OUTPUT
    {
        wstring fileName = wformat(L"output\\d_lvl_%02d.png", m_level);
        CRGBImg img;
        VT_HR_EXIT( m_imgPyr->ReadImg(img, m_level) );
        VT_HR_EXIT( VtSaveImage(fileName, img) );
    }
#endif DEBUG_OUTPUT

    VT_HR_END()
}

template <typename TUV>
void UpdateBBox(CRect & bbox, const CTypedImg<TUV> & uv, const CRect & rctUv)
{
    CRect rct = rctUv & uv.Rect();
    if (rct.IsRectEmpty())
        return;

    CVec2<TUV> * uvPtr = (CVec2<TUV> *) uv.Ptr(rct.left, rct.top);
    int uvOff = uv.StrideBytes() / sizeof(CVec2<TUV>);

    int xl = 0;

    if (g_SupportSSE2() && ElTraits<TUV>::ElFormat() == EL_FORMAT_SHORT)
    {
        // Since 0xffff is invalid, we want it to be 0x7fff (32767) when doing min
        // and 0x8000 (-32768) when doing max
        __m128i xi0 = _mm_set1_epi16(-32768);   // convert to/from signed
        __m128i xi1 = _mm_set1_epi16(+32767);   // convert to/from signed + 1
        __m128i xi2 = _mm_set1_epi16(+32767);   // min
        __m128i xi3 = _mm_set1_epi16(-32768);   // max

        for (int y = 0; y < rct.Height(); y++)
        {
            for (xl = 0; xl < rct.Width() - 3; xl += 4)
            {
			    __m128i xi4 = _mm_loadu_si128 ((__m128i*) (uvPtr + xl));

                __m128i xi5 = _mm_sub_epi16(xi4, xi0);  // convert to signed
                xi2 = _mm_min_epi16(xi5, xi2);  // min

                __m128i xi6 = _mm_sub_epi16(xi4, xi1);  // convert to signed + 1
			    xi3 = _mm_max_epi16(xi6, xi3);  // max
            }

            uvPtr += uvOff;
        }

        xi2 = _mm_add_epi16(xi2, xi0);  // convert to unsigned
        int minx = 
            VtMin(_mm_extract_epi16(xi2, 0), VtMin(_mm_extract_epi16(xi2, 2),
            VtMin(_mm_extract_epi16(xi2, 4), _mm_extract_epi16(xi2, 6))));

        if (minx != ElTraits<TUV>::MaxVal())
        {
            int miny = 
                VtMin(_mm_extract_epi16(xi2, 1), VtMin(_mm_extract_epi16(xi2, 3),
                VtMin(_mm_extract_epi16(xi2, 5), _mm_extract_epi16(xi2, 7))));

            xi3 = _mm_add_epi16(xi3, xi0);  // convert to unsigned + 1
            int maxx =
                VtMax(_mm_extract_epi16(xi3, 0), VtMax(_mm_extract_epi16(xi3, 2),
                VtMax(_mm_extract_epi16(xi3, 4), _mm_extract_epi16(xi3, 6))));
            int maxy = 
                VtMax(_mm_extract_epi16(xi3, 1), VtMax(_mm_extract_epi16(xi3, 3),
                VtMax(_mm_extract_epi16(xi3, 5), _mm_extract_epi16(xi3, 7))));

            bbox.left   = VtMin(bbox.left,   (long) minx);
            bbox.top    = VtMin(bbox.top,    (long) miny);
            bbox.right  = VtMax(bbox.right,  (long) maxx);
            bbox.bottom = VtMax(bbox.bottom, (long) maxy);
        }
    }

    uvPtr = (CVec2<TUV> *) uv.Ptr(rct.left, rct.top);

    for (int y = 0; y < rct.Height(); y++)
    {
        for (int x = xl; x < rct.Width(); x++)
        {
            TUV nx = uvPtr[x].x;
            if (nx == ElTraits<TUV>::MaxVal())
                continue;
            TUV ny = uvPtr[x].y;

            bbox.left   = VtMin(bbox.left,   (long) nx);
            bbox.top    = VtMin(bbox.top,    (long) ny);
            bbox.right  = VtMax(bbox.right,  (long) nx + 1);
            bbox.bottom = VtMax(bbox.bottom, (long) ny + 1);
        }

        uvPtr += uvOff;
    }
}

template <typename TUV>
HRESULT Completion::upsample()
{
    VT_HR_BEGIN()

    if (m_level == 0)
    {
        return hr;
    }

    VT_DEBUG_LOG( "upsampling...   (level %d --> %d)\n", m_level, m_level-1 );

    float thisTimeUpsample = m_timer.GetTimeMilliSec() * 0.001f;

    // Save previous level's bounding boxes.
    CImgInfo info1 = m_uvPyr.GetImgInfo(m_level);
    int blockSize1 = m_blockSize;
    int xBlocks1 = (info1.width + blockSize1 - 1) / blockSize1;
    vector<CRect> bboxes1 = m_bboxes[0];

    m_level--;

#ifdef DEBUG_OUTPUT
    float startTime, endTime;

    startTime = m_timer.GetTimeMilliSec() * 0.001f;
#endif DEBUG_OUTPUT

    // Can store valid pixels for smaller levels.
	vector<CVec2<TUV>> validPixels;
    if (m_level > m_tileLevel)
        VT_HR_EXIT( initLevel(&validPixels) );
    else
        VT_HR_EXIT( initLevel<int>() );

#ifdef DEBUG_OUTPUT
    endTime = m_timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "UPS --- INIT: %f\n", endTime - startTime );

    startTime = endTime;
#endif DEBUG_OUTPUT

    // nn
    CRGBAImg img, src;
    CTypedImg<TUV> newUv, uv;

#ifdef DEBUG_OUTPUT
    endTime = m_timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "UPS --- FILL: %f\n", endTime - startTime );

    startTime = endTime;
#endif DEBUG_OUTPUT

    CImgInfo info = m_uvPyr.GetImgInfo(m_level);

    int w = info.width;
    int h = info.height;

    int xBlocks = (w + m_blockSize - 1) / m_blockSize;
    int yBlocks = (h + m_blockSize - 1) / m_blockSize;
    VT_HR_EXIT( m_bboxes[0].resize(xBlocks * yBlocks) );
    VT_HR_EXIT( m_bboxes[1].resize(xBlocks * yBlocks) );
    VT_HR_EXIT( m_bboxes[2].resize(xBlocks * yBlocks) );
    for (int i = 0; i < (int) m_bboxes[0].size(); i++)
        m_bboxes[0][i] = m_bboxes[1][i] = m_bboxes[2][i] = CRect(LONG_MAX, LONG_MAX, 0, 0);

#ifdef USE_TILE_BBOXES
    VT_HR_EXIT( m_tileUnions.resize(xBlocks * yBlocks) );
    for (int i = 0; i < (int) m_tileUnions.size(); i++)
        m_tileUnions[i] = CRect(LONG_MAX, LONG_MAX, 0, 0);
#endif USE_TILE_BBOXES

    int tileShift = m_tileLevel + Tiles::levelOffset - m_level;
    int tileWidth = m_restriction[0].Width();
    int tileHeight = m_restriction[0].Height();
    int iTile = -1;
    CRect rctTile;
    CRGBAImg imgTile;

    int block = 0;
    for (CBlockIterator bi(info.Rect(), m_blockSize); !bi.Done(); bi.Advance(), block++)
    {
        // read current image
        CRect rctNew = bi.GetRect();
        VT_HR_EXIT( m_imgPyr->ReadRegion(rctNew, img, m_level) );
        int imgoff = img.StrideBytes() / sizeof(RGBAPix);

        // create UV image for current level
        VT_HR_EXIT( CreateImageForTransform(newUv, rctNew.Width(), rctNew.Height(),
            VT_IMG_MAKE_TYPE(ElTraits<TUV>::ElFormat(), 2)) );
        VT_HR_EXIT( newUv.Fill(ElTraits<TUV>::MaxVal()) );
        int newoff = newUv.StrideBytes() / sizeof(CVec2<TUV>);

        // find source bounding box from previous level
        CRect rctOld = CRect(rctNew.left >> 1, rctNew.top >> 1, rctNew.right >> 1, rctNew.bottom >> 1);
        int block1 = (rctOld.left / blockSize1) + (rctOld.top / blockSize1) * xBlocks1;
        CRect rctSrc = bboxes1[block1];
        if (rctSrc.IsRectEmpty())
        {
            VT_HR_EXIT( m_uvPyr.WriteRegion(rctNew, newUv, m_level) );
            continue;
        }

        // make and read source bounding box for current level
        rctSrc.left <<= 1;
        if (rctNew.left == 0)
            rctSrc.left -= s_prad;
        rctSrc.top <<= 1;
        if (rctNew.top == 0)
            rctSrc.top -= s_prad;
        rctSrc.right <<= 1;
        if (rctNew.right > info.width - 2 * s_prad)
            rctSrc.right += s_prad;
        rctSrc.bottom <<= 1;
        if (rctNew.bottom > info.height - 2 * s_prad)
            rctSrc.bottom += s_prad;
        rctSrc &= info.Rect();
        VT_HR_EXIT( m_imgPyr->ReadRegion(rctSrc, src, m_level) );

        // read UV image from previous level
        if (rctOld.Width() <= s_prad)
            rctOld.left -= (s_prad + 1 - rctOld.Width());
        if (rctOld.Height() <= s_prad)
            rctOld.top -= (s_prad + 1 - rctOld.Height());
        rctOld &= m_uvPyr.GetImgInfo(m_level + 1).Rect();
        VT_HR_EXIT( m_uvPyr.ReadRegion(rctOld, uv, m_level + 1) );

        int x0 = VtMax((int) rctNew.left, s_prad);
        int x1 = VtMin((int) rctNew.right, w - s_prad);
        int y0 = VtMax((int) rctNew.top, s_prad);
        int y1 = VtMin((int) rctNew.bottom, h - s_prad);
        RGBAPix * pimg = img.Ptr(y0 - rctNew.top);
        CVec2<TUV> * pnew = (CVec2<TUV> *) newUv.Ptr(y0 - rctNew.top);

        for (int y = y0; y < y1; y++)
        {
            for (int x = x0; x < x1; x++)
            {
                if (pimg[x - rctNew.left].a & s_validBit)
                    continue;

                CVec2<TUV> t((TUV) x, (TUV) y);
                const int dx = clamp(t.x / 2, s_prad, info1.width  - s_prad - 1);
                const int dy = clamp(t.y / 2, s_prad, info1.height - s_prad - 1);
                CVec2<TUV> old = *((CVec2<TUV> *) uv.Ptr(dx - rctOld.left, dy - rctOld.top));
                int sx = old.x * 2 + (t.x - dx * 2);
                int sy = old.y * 2 + (t.y - dy * 2);

                if (sx >= 0 && sx < w && sy >= 0 && sy < h &&
                    (src(sx - rctSrc.left, sy - rctSrc.top, 3) & s_validBit) != 0)
                {
                    pnew[x - rctNew.left].x = (TUV) sx;
                    pnew[x - rctNew.left].y = (TUV) sy;
                }
                else
                {
                    CVec2<TUV> rs;
                    int tile = tileWidth * (t.y >> tileShift) + (t.x >> tileShift);
                    if (!imgTile.IsValid() || tile != iTile)
                    {
                        rctTile = CRect(
                            m_tiles[tile].left   << (m_tileLevel + Tiles::levelOffset - m_level),
                            m_tiles[tile].top    << (m_tileLevel + Tiles::levelOffset - m_level),
                            m_tiles[tile].right  << (m_tileLevel + Tiles::levelOffset - m_level),
                            m_tiles[tile].bottom << (m_tileLevel + Tiles::levelOffset - m_level));
                        VT_HR_EXIT( m_imgPyr->ReadRegion(rctTile, imgTile, m_level) );
                    }
                    iTile = tile;

                    int maxTries = 100000, curTry = 0;

                    for (curTry; curTry < maxTries; curTry++)
                    {
                        int tx = m_rnd.IRand(tileWidth);
                        int ty = m_rnd.IRand(tileHeight);

                        if (m_restriction[tile](tx, ty) < FLT_MAX)
                        {
                            int tx0 = tx << tileShift, ty0 = ty << tileShift;
                            int tx1 = (tx+1) << tileShift, ty1 = (ty+1) << tileShift;
                            tx1 = VtMin(tx1, w); ty1 = VtMin(ty1, h);

                            rs.x = (TUV) (tx0 + m_rnd.IRand(tx1 - tx0));
                            rs.y = (TUV) (ty0 + m_rnd.IRand(ty1 - ty0));

                            if (imgTile(rs.x - rctTile.left,
                                        rs.y - rctTile.top, 3) & s_validBit)
                            {
                                break;
                            }
                        }
                    }
                    if (curTry == maxTries)
                    {
                        rs = (!validPixels.empty() ? validPixels[m_rnd.IRand((int) validPixels.size())] :
                            CVec2<TUV>(ElTraits<TUV>::MaxVal(), ElTraits<TUV>::MaxVal()));
                        VT_DEBUG_LOG( "failed to upsample pixel %d, %d; using %d, %d\n", t.x, t.y, rs.x, rs.y );
                    }
                    else
                        VT_DEBUG_LOG( "invalid pixel %d, %d; using %d, %d\n", t.x, t.y, rs.x, rs.y );

                    pnew[x - rctNew.left].x = rs.x;
                    pnew[x - rctNew.left].y = rs.y;
                }
            }

            pimg += imgoff;
            pnew += newoff;
        }
 
        VT_HR_EXIT( m_uvPyr.WriteRegion(rctNew, newUv, m_level) );

        {
            UpdateBBox(m_bboxes[0][block], newUv, CRect(newUv.Rect()));

            m_bboxes[2][block] |= m_bboxes[0][block];

            // iteration bounding box is 1 larger along trailing edges

            // add right edge for down pass
            if (rctNew.right != w)
                UpdateBBox(m_bboxes[0][block + 1], newUv, CRect(rctNew.Width() - 1, 0, rctNew.Width(), rctNew.Height()));
            // add bottom edge for down pass
            if (rctNew.bottom != h)
                UpdateBBox(m_bboxes[0][block + xBlocks], newUv, CRect(0, rctNew.Height() - 1, rctNew.Width(), rctNew.Height()));

            // average bounding box is s_prad larger all around

            // add right edge
            if (rctNew.right != w)
                UpdateBBox(m_bboxes[2][block + 1], newUv, CRect(rctNew.Width() - s_prad, 0, rctNew.Width(), rctNew.Height()));
            // add bottom edge
            if (rctNew.bottom != h)
                UpdateBBox(m_bboxes[2][block + xBlocks], newUv, CRect(0, rctNew.Height() - s_prad, rctNew.Width(), rctNew.Height()));
            // add bottom right corner
            if (rctNew.right != w && rctNew.bottom != h)
                UpdateBBox(m_bboxes[2][block + xBlocks + 1], newUv, CRect(rctNew.Width() - s_prad, rctNew.Height() - s_prad, rctNew.Width(), rctNew.Height()));
            // add bottom left corner
            if (rctNew.left != 0 && rctNew.bottom != h)
                UpdateBBox(m_bboxes[2][block + xBlocks - 1], newUv, CRect(0, rctNew.Height() - s_prad, s_prad, rctNew.Height()));
            // add left edge
            if (rctNew.left != 0)
                UpdateBBox(m_bboxes[2][block - 1], newUv, CRect(0, 0, s_prad, rctNew.Height()));
            // add top edge
            if (rctNew.top != 0)
                UpdateBBox(m_bboxes[2][block - xBlocks], newUv, CRect(0, 0, rctNew.Width(), s_prad));
            // add top left corner
            if (rctNew.left != 0 && rctNew.top != 0)
                UpdateBBox(m_bboxes[2][block - xBlocks - 1], newUv, CRect(0, 0, s_prad, s_prad));
            // add top right corner
            if (rctNew.right != w && rctNew.top != 0)
                UpdateBBox(m_bboxes[2][block - xBlocks + 1], newUv, CRect(rctNew.Width() - s_prad, 0, rctNew.Width(), s_prad));

#ifdef USE_TILE_BBOXES
            for (int y = rctNew.top; y < rctNew.bottom; y += Tiles::blockSize << (m_tileLevel - m_level))
            {
                for (int x = rctNew.left; x < rctNew.right; x += Tiles::blockSize << (m_tileLevel - m_level))
                {
                    int tile = tileWidth * (y >> tileShift) + (x >> tileShift);
                    CRect bbox = m_tileBBoxes[tile];
                    bbox = CRect(bbox.left   << (m_tileLevel + Tiles::levelOffset - m_level),
                                 bbox.top    << (m_tileLevel + Tiles::levelOffset - m_level),
                                 bbox.right  << (m_tileLevel + Tiles::levelOffset - m_level),
                                 bbox.bottom << (m_tileLevel + Tiles::levelOffset - m_level));
                    m_tileUnions[block] |= m_tileBBoxes[tile];
                }
            }
#endif USE_TILE_BBOXES
        }
   }

#ifdef DEBUG_OUTPUT
    endTime = m_timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "UPS --- UV: %f\n", endTime - startTime );

    startTime = endTime;
#endif DEBUG_OUTPUT

    float thisTimeUpsampleAverage = m_timer.GetTimeMilliSec() * 0.001f;
    VT_ASSERT( m_mark == 0 );
    VT_HR_EXIT( average() );
    m_timeUpsampleAverage += m_timer.GetTimeMilliSec() * 0.001f - thisTimeUpsampleAverage;

#ifdef DEBUG_OUTPUT
    endTime = m_timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "UPS --- AVG: %f\n", endTime - startTime );

    startTime = endTime;
#endif DEBUG_OUTPUT

    thisTimeUpsample = m_timer.GetTimeMilliSec() * 0.001f - thisTimeUpsample;
    m_timeUpsample += thisTimeUpsample;

#ifdef DEBUG_OUTPUT
    endTime = m_timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "UPS --- PASS: %f\n", endTime - startTime );
#endif DEBUG_OUTPUT

    VT_DEBUG_LOG( "  [done] took %.2fs\n", thisTimeUpsample );

    VT_HR_END()
}
