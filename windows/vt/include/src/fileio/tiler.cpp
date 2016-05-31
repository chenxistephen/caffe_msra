//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Helper class for dealing with tiled images
//
//  History:
//      2006/10/16-mattu
//          Created
//
//------------------------------------------------------------------------
#include "stdafx.h"

#if !defined(VT_WINRT)

#include <shlwapi.h>

#include "wicio.h"
#include "tiler.h"

using namespace vt;

#define HDPHOTO_PAD    8

//+---------------------------------------------------------------------------
//
//  Utilities
//
//-----------------------------------------------------------------------------
// 
template <class T> HRESULT
vt::TileNameFormat(T& dst, eTileHierarchyFormat eFrmt, const WCHAR* base,
				   const WCHAR* ext, int l, int x, int y)
{
	// We use "/" rather than "\" for the path separator because it works both
	// for Windows APIs and within zip archives intended for use on UNIX systems.
    HRESULT hr = E_INVALIDARG;
	switch( eFrmt )
	{
	case eTileHierarchyFlat:
		hr = base != NULL ?
            dst.format(L"%s/%d_%d_%d.%s", base, l, x, y, ext) :
            dst.format(L"%d_%d_%d.%s", l, x, y, ext);
		break;
	case eTileHierarchyFolderPerCol:
		hr = base != NULL ?
            dst.format(L"%s/l_%d/c_%d/tile_%d.%s", base, l, x, y, ext) :
            dst.format(L"l_%d/c_%d/tile_%d.%s", l, x, y, ext);
		break;
	case eTileHierarchySeadragon:
	case eTileHierarchyPhotosynth:
		hr = base != NULL ?
            dst.format(L"%s/%d/%d_%d.%s", base, l, x, y, ext) :
            dst.format(L"%d/%d_%d.%s", l, x, y, ext);
		break;
    default:
        VT_ASSERT(0);
        break;
	}
    return hr;
}

HRESULT 
CreateLevelDirectory(const WCHAR* base, eTileHierarchyFormat eFrmt, int l)
{
    HRESULT hr = E_INVALIDARG;

    vt::wstring_b<MAX_PATH> dirname;

    switch ( eFrmt )
    {
    case eTileHierarchyFlat:
        VT_ASSERT( 0 );
        break;
    case eTileHierarchyFolderPerCol:
        hr = base != NULL ?
             dirname.format(L"%s\\l_%d", base, l) :
             dirname.format(L"l_%d", l);
        break;
    case eTileHierarchySeadragon:
        hr = base != NULL ?
             dirname.format(L"%s\\%d", base, l) :
             dirname.format(L"%d", l);
        break;
    }

    if (hr == S_OK && !CreateDirectoryW(dirname, NULL))
    {
        DWORD dwError = GetLastError();
        hr = dwError == ERROR_ALREADY_EXISTS ? S_OK :
            HRESULT_FROM_WIN32( dwError );
    }

    return hr;
}

HRESULT 
CreateColumnDirectory(const WCHAR* base, int l, int c)
{
    HRESULT hr = E_INVALIDARG;

    vt::wstring_b<MAX_PATH> dirname;
    hr = base != NULL ? dirname.format(L"%s\\l_%d\\c_%d", base, l, c) :
                        dirname.format(L"l_%d\\c_%d", l, c);
    if (hr == S_OK && !CreateDirectoryW(dirname, NULL))
    {
        DWORD dwError = GetLastError();
        hr = dwError == ERROR_ALREADY_EXISTS ? S_OK :
            HRESULT_FROM_WIN32( dwError );
    }

    return hr;
}

HRESULT
vt::LoadTiles(const TILE_SRC_FRMT& srcfrmt, int level, const CRect& rct, 
			  OUT CImg& imgBuff, bool bIgnoreMissing)
{
    HRESULT hr = NOERROR;

    int iEffectiveTileW = srcfrmt.iTileWidth  - 2*srcfrmt.iOverlap;
    int iEffectiveTileH = srcfrmt.iTileHeight - 2*srcfrmt.iOverlap;

    CImg imgIn;

    // load the first tile to determine its pixel format
    // TODO: handle missing tile case
    vt::wstring_b<MAX_PATH> fname;
    TileNameFormat(fname, srcfrmt.eFrmt, srcfrmt.pBase, srcfrmt.pExtension,
				   level, rct.left/iEffectiveTileW, rct.top/iEffectiveTileH);
    VT_HR_EXIT( VtLoadImage(fname, imgIn) );

    // create the appropriate imgBuff
	int iDefType = imgIn.GetType();
    VT_HR_EXIT( CreateImageForTransform(imgBuff, rct.Width(), rct.Height(), 
									    iDefType) );

    for ( int x = rct.left; x < rct.right; )
    {
        int iXTIn = x / iEffectiveTileW;

        for ( int y = rct.top; y < rct.bottom; )
        {
            int iYTIn = y / iEffectiveTileH;

            // read the input tile(s)
            TileNameFormat(fname, srcfrmt.eFrmt, srcfrmt.pBase,
                           srcfrmt.pExtension, level, iXTIn, iYTIn);

            if( bIgnoreMissing )
            {
                if( VtLoadImage(fname, imgIn) != S_OK )
                {
                    VT_HR_EXIT( imgIn.Create(srcfrmt.iTileWidth, 
											 srcfrmt.iTileHeight, iDefType, 
											 DEFAULT_ALIGN) );
                    imgIn.Clear();
                }
            }
            else
            {
                VT_HR_EXIT( VtLoadImage(fname, imgIn) );
            }

            // copy to appropriate sub-region
            CRect rctSubIn;
			int iOvlX = (srcfrmt.eFrmt == eTileHierarchyPhotosynth || iXTIn != 0)? 
				srcfrmt.iOverlap: 0;
			int iOvlY = (srcfrmt.eFrmt == eTileHierarchyPhotosynth || iYTIn != 0)? 
				srcfrmt.iOverlap: 0;
			rctSubIn.left   = x - iXTIn*iEffectiveTileW + iOvlX;
			rctSubIn.top    = y - iYTIn*iEffectiveTileH + iOvlY;
			rctSubIn.right  = iOvlX + min(iEffectiveTileW,
                                          rct.right - iXTIn*iEffectiveTileW);
			rctSubIn.bottom = iOvlY + min(iEffectiveTileH,
                                          rct.bottom - iYTIn*iEffectiveTileH);

            CRect rctSubBuff;
            rctSubBuff.left   = x - rct.left;
            rctSubBuff.top    = y - rct.top;
            rctSubBuff.right  = rctSubBuff.left + rctSubIn.Width();
            rctSubBuff.bottom = rctSubBuff.top  + rctSubIn.Height();

            CImg imgSubIn, imgSubBuff;
            VT_HR_EXIT( imgIn.Share(imgSubIn, &rctSubIn) );
            VT_HR_EXIT( imgBuff.Share(imgSubBuff, &rctSubBuff) );
            VT_HR_EXIT( VtConvertImage( imgSubBuff, imgSubIn ) );

            if( iYTIn * iEffectiveTileH != y )
            {
                y = iYTIn * iEffectiveTileH;
            }
			y += iEffectiveTileH;
        }
        if( iXTIn * iEffectiveTileW != x )
        {
            x = iXTIn * iEffectiveTileW;
        }
        x += iEffectiveTileW;
    }

Exit:
    return hr;
}

static HRESULT
SaveTileWithPadding(const CParams* pEncParams,
                    LPCWSTR tilename, const CRect& rctSrc, const CRect& rctPad,
                    IVTImageWriter* pWriter, IImageReader* pSrc,
                    IArchiveWriter* pArchive)
{
    HRESULT hr;

    CMemStream stream, memstream;

    CRect rctInTile = rctSrc - rctPad.TopLeft();

    // compress the padded tile
    VT_HR_EXIT( pWriter->OpenFile(&stream, L".wdp") );
    VT_HR_EXIT( pWriter->SetImage(pSrc, &rctPad, true, pEncParams) );
    VT_HR_EXIT( pWriter->CloseFile() );

    // extract the interior tile and save it
    // rewind the stream
    LARGE_INTEGER zero;
    zero.QuadPart = 0;
    VT_HR_EXIT( stream.Seek( zero, STREAM_SEEK_SET, NULL) );

    if (pArchive != NULL)
        VT_HR_EXIT( pWriter->OpenFile(&memstream, L".wdp") );
    else
        VT_HR_EXIT( pWriter->OpenFile(tilename) );
    VT_HR_EXIT( pWriter->SetImage(&stream, rctInTile, NULL) );
    VT_HR_EXIT( pWriter->CloseFile() );
    if (pArchive != NULL)
        VT_HR_EXIT( pArchive->AddFile(tilename, memstream) );

Exit:
    return hr;
}

static HRESULT
SaveTile(HTRANSFORM hTransform, const CParams* pEncParams,
         LPCWSTR tilename, const CRect& rctSrc, 
         IVTImageWriter* pWriter, IImageReader* pSrc,
         IArchiveWriter* pArchive)
{
    HRESULT hr;

    CMemStream memstream;

    CImg imgTile;
    if (hTransform != NULL)
    {
        VT_HR_EXIT( pSrc->ReadRegion(rctSrc, imgTile) );

        BMFORMAT bmFormat = imgTile.Bands() == 1 ? BM_GRAY :
            (imgTile.Bands() == 3 ? BM_RGBTRIPLETS : BM_xRGBQUADS);
        if (!TranslateBitmapBits(hTransform,
                                 imgTile.BytePtr(), bmFormat,
                                 rctSrc.Width(), rctSrc.Height(),
                                 (DWORD) imgTile.StrideBytes(),
                                 imgTile.BytePtr(), bmFormat,
                                 (DWORD) imgTile.StrideBytes(),
                                 NULL, 0L))
        {    
            VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );
        }
    }

    memstream.Clear();
    if (pArchive != NULL)
    {    
        VT_HR_EXIT( pWriter->OpenFile(&memstream, VtGetFileExt(tilename)) );
    }
    else
    {
        VT_HR_EXIT( pWriter->OpenFile(tilename) );
    }
    if (imgTile.IsValid())
        VT_HR_EXIT( pWriter->SetImage(imgTile, true, pEncParams) );
    else
        VT_HR_EXIT( pWriter->SetImage(pSrc, &rctSrc, true, pEncParams) ); 
    VT_HR_EXIT( pWriter->CloseFile() );
    if (pArchive != NULL)
    {
        VT_HR_EXIT( pArchive->AddFile(tilename, memstream) );
    }

Exit:
    return hr;
}

static HRESULT
SetSaveTilesParams(const TILE_DST_FRMT& dstfrmt, CParams& params,
                   HTRANSFORM& hTransform)
{
    HPROFILE hProfiles[2] = { NULL, NULL };

    VT_HR_BEGIN()

    if (dstfrmt.eProf != eTileProfileUnspecified &&
        dstfrmt.eProf != eTileProfileSRGB)
        VT_HR_EXIT( E_NOTIMPL );

    if (EL_FORMAT(dstfrmt.iImgType) == EL_FORMAT_BYTE &&
        dstfrmt.eProf == eTileProfileSRGB)
    {
        const CParamValue *pColorProfile = NULL;
        params.GetById(&pColorProfile, ImagePropertyICCProfile);
        if (pColorProfile != NULL)
        {
            BOOL bValid = FALSE;

            // Open the source color profile.
            PROFILE srcProfile = { PROFILE_MEMBUFFER, (PVOID) pColorProfile->GetDataPtr(), 
                                   (DWORD) pColorProfile->GetDataSize() };
            hProfiles[0] = OpenColorProfile(&srcProfile, PROFILE_READ, 
                                            FILE_SHARE_READ, OPEN_EXISTING);
            if ( hProfiles[0] == NULL ||
                 !IsColorProfileValid(hProfiles[0], &bValid) || !bValid )
            {
		        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );
            }

            // Open the destination color profile.
            WCHAR pwszPath[MAX_PATH];
            PROFILE dstProfile = { PROFILE_FILENAME, pwszPath, MAX_PATH };
            if (!GetStandardColorSpaceProfile(NULL, LCS_sRGB, (PWSTR) dstProfile.pProfileData, &dstProfile.cbDataSize) ||
                (hProfiles[1] = OpenColorProfile(&dstProfile, PROFILE_READ, FILE_SHARE_READ, OPEN_EXISTING)) == NULL ||
                !IsColorProfileValid(hProfiles[1], &bValid) || !bValid)
		        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

            // Create the ICM or WCS transform.
            DWORD dwIntents[2] = { INTENT_RELATIVE_COLORIMETRIC, INTENT_RELATIVE_COLORIMETRIC };
            if ((hTransform = CreateMultiProfileTransform(hProfiles, 2, dwIntents, 2,
                BEST_MODE | USE_RELATIVE_COLORIMETRIC, INDEX_DONT_CARE)) == NULL)
		        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );
        }

        params.DeleteById(ImagePropertyICCProfile);
    }

    if (_wcsicmp(dstfrmt.pExtension, L"png") != 0)
    {
		// set quality for everything BUT png
        CParamValue v;
        v.Set(float(dstfrmt.fQuality));
        VT_HR_EXIT( params.SetByName(L"ImageQuality", 0, v) );
    }
    if (_wcsicmp(dstfrmt.pExtension, L"wdp") == 0)
    {
        CParamValue v;
        v.Set(dstfrmt.fQuality >= 1.f);
        VT_HR_EXIT( params.SetByName(L"Lossless", 0, v) );
    }

    VT_HR_EXIT_LABEL()

    if (hProfiles[0] != NULL)
        CloseColorProfile(hProfiles[0]);
    if (hProfiles[1] != NULL)
        CloseColorProfile(hProfiles[1]);

    return hr;
}

HRESULT
vt::SaveTiles(const TILE_DST_FRMT& dstfrmt, const CParams* pEncParams, 
              UINT uTileLevel, const CPoint& ptSrc, const CRect& rctDst, 
              IVTImageWriter* pWriter, IImageReader* pSrc,
              IArchiveWriter* pArchive, CTaskProgress* pProgress)
{
    HTRANSFORM hTransform = NULL;

    VT_HR_BEGIN()

    VT_PTR_EXIT( pSrc );

    CParams params;
    if (pEncParams != NULL)
        params = *pEncParams;
    VT_HR_EXIT( SetSaveTilesParams(dstfrmt, params, hTransform) );

    CImgInfo info = pSrc->GetImgInfo();
    if ( ptSrc.x + rctDst.Width()  > info.width ||
         ptSrc.y + rctDst.Height() > info.height)
         VT_HR_EXIT( E_INVALIDARG );

    pArchive = dstfrmt.eFile == eTileArchiveZip ? pArchive : NULL;

    WCHAR chBase[MAX_PATH], *pBase = NULL;
    if (dstfrmt.pBase != NULL)
    {
        pBase = chBase;
        wcscpy_s(pBase, MAX_PATH, dstfrmt.pBase);
        PathStripPath(pBase);
    }

    if( dstfrmt.eFile == eTileArchiveFileSystem && 
        dstfrmt.eFrmt != eTileHierarchyFlat )
    {
        CreateLevelDirectory(dstfrmt.pBase, dstfrmt.eFrmt, uTileLevel);
    }

    int iStepX = dstfrmt.iTileWidth  - 2*dstfrmt.iOverlap;
    int iStepY = dstfrmt.iTileHeight - 2*dstfrmt.iOverlap;

    VT_ASSERT( (rctDst.left % iStepX) == 0 );
    VT_ASSERT( (rctDst.top  % iStepY) == 0 );

    // vars used for progress
    int iBlkCntX = (rctDst.Width()  + dstfrmt.iTileWidth - 1)/dstfrmt.iTileWidth;
    int iBlkCntY = (rctDst.Height() + dstfrmt.iTileHeight - 1)/dstfrmt.iTileHeight;
    float fProgStep = 1.f/(float(iBlkCntX)*float(iBlkCntY));
    float fProgress = 0;

    for ( CBlockIterator bi(BLOCKITER_INIT(rctDst, iStepX, iStepY)); !bi.Done(); 
          bi.Advance() )
    {
        CRect rctCmp = bi.GetCompRect();

        int xt = rctCmp.left / iStepX;
        int yt = rctCmp.top  / iStepY;

        if( dstfrmt.eFile == eTileArchiveFileSystem && 
            dstfrmt.eFrmt == eTileHierarchyFolderPerCol && yt == 0 )
        {
            CreateColumnDirectory(dstfrmt.pBase, uTileLevel, xt);
        }

        CRect rctSrc = bi.GetRect() + ptSrc;
        rctSrc.InflateRect(dstfrmt.iOverlap, dstfrmt.iOverlap);
		CRect infoRect = info.Rect();
        rctSrc.IntersectRect(&rctSrc, &infoRect);

        vt::wstring_b<MAX_PATH> tilename;
        TileNameFormat(tilename, dstfrmt.eFrmt,
                       dstfrmt.eFile == eTileArchiveFileSystem ? dstfrmt.pBase : pBase,
                       dstfrmt.pExtension, uTileLevel, xt, yt);

        // use padding when constructing wdp tiles unless lossless compression
        // is being done - then padding isn't necessary
        if( _wcsicmp( dstfrmt.pExtension, L"wdp" ) == 0 &&
	        dstfrmt.fQuality < 1.f )
        {
            CRect rctPad = rctSrc;
            rctPad.InflateRect(HDPHOTO_PAD, HDPHOTO_PAD);
			CRect ir = info.Rect();
	        rctPad.IntersectRect(&rctPad, &ir);
            VT_HR_EXIT( SaveTileWithPadding(&params, tilename, rctSrc, rctPad,
                                             pWriter, pSrc, pArchive) );
        }
        else
        {
            VT_HR_EXIT( SaveTile(hTransform, &params, tilename, rctSrc,
                                  pWriter, pSrc, pArchive) );
        }

        // report progress
        fProgress += fProgStep;
        VT_HR_EXIT( CheckTaskProgressCancel(pProgress, fProgress) );
    }

    VT_HR_EXIT_LABEL()

    if (hTransform != NULL)
        DeleteColorTransform(hTransform);

    SetProgress100Percent(pProgress);

    return hr;
}

HRESULT
vt::DetermineTileSetDim(const TILE_SRC_FRMT& srcfrmt, int level, 
                        OUT int& width, OUT int& height, OUT int& bands)
{
    // visit tiles in x direction
    int x = 0;
    width = 0;
    while( 1 )
    {
        vt::wstring_b<MAX_PATH> fname;
        TileNameFormat(fname, srcfrmt.eFrmt, srcfrmt.pBase, srcfrmt.pExtension, 
                       level, x, 0);

        CImg img;
        if( VtLoadImage(fname, img) != S_OK )
        {
            break;
        }
        width += img.Width() - 2*srcfrmt.iOverlap;
#if 0
        // TODO: when wdp supports filesrc then replace this
        CImageFileSrc fs;
        if( fs.OpenFile((const WCHAR*)fname) != S_OK )
        {
            break;
        }
        width += fs.Width();
#endif

        bands = img.Bands();
        x++;
    }
	// dont subtract off the first and last
	width += 2*srcfrmt.iOverlap;


    // visit tiles in y direction
    int y = 0;
    height = 0;
    while( 1 )
    {
        vt::wstring_b<MAX_PATH> fname;
        TileNameFormat(fname, srcfrmt.eFrmt, srcfrmt.pBase, srcfrmt.pExtension, 
                       level, 0, y);

        CImg img;
        if( VtLoadImage(fname, img) != S_OK )
        {
            break;
        }
        height += img.Height() - 2*srcfrmt.iOverlap;
#if 0
        // TODO: when wdp supports filesrc then replace this
        CImageFileSrc fs;
        if( fs.OpenFile(fname) != S_OK )
        {
            break;
        }
        height += fs.Height();
#endif
        
        y++;
    }	

    // dont subtract off the first and last one
    height += 2*srcfrmt.iOverlap;

    return S_OK;
}

//+-----------------------------------------------------------------------------
//
// Class: CSaveTilesTransform
// 
//------------------------------------------------------------------------------
HRESULT
CSaveTilesTransform::Begin()
{
    VT_HR_BEGIN()

    if (m_dstfrmt.eFile == eTileArchiveFileSystem && 
        m_dstfrmt.eFrmt != eTileHierarchyFlat)
    {
        VT_HR_EXIT( CreateLevelDirectory(m_dstfrmt.pBase, m_dstfrmt.eFrmt,
                                          m_uTileLevel) );
    }

    if (m_dstfrmt.eFile == eTileArchiveFileSystem && 
        m_dstfrmt.eFrmt == eTileHierarchyFolderPerCol)
    {
        int iStepX = m_dstfrmt.iTileWidth  - 2*m_dstfrmt.iOverlap;

		int xt = m_rctSrc.left - m_ptDstOrigin.x;
		VT_ASSERT( (xt % iStepX) == 0 );
		xt /= iStepX;

		for( int x = 0; x < m_rctSrc.Width(); x+=iStepX, xt++ )
        {   
            VT_HR_EXIT( CreateColumnDirectory(m_dstfrmt.pBase, m_uTileLevel, 
											   xt) );
        }
    }

    VT_HR_END()
}

HRESULT
CSaveTilesTransform::SetSaveTilesParams()
{
    return ::SetSaveTilesParams(m_dstfrmt, m_params, m_hTransform);
}

HRESULT
CSaveTilesTransform::Transform(OUT CImg*, 
                               IN  const CRect& rctLayerDst,
                               IN  const CImg& imgSrc,
                               IN  const vt::CPoint& ptSrc)
{
    VT_HR_BEGIN()

    int iStepX = m_dstfrmt.iTileWidth  - 2*m_dstfrmt.iOverlap;
    int iStepY = m_dstfrmt.iTileHeight - 2*m_dstfrmt.iOverlap;

    vt::CPoint ptTile = rctLayerDst.TopLeft() - m_ptDstOrigin; 
    VT_ASSERT( (ptTile.x % iStepX) == 0 );
    VT_ASSERT( (ptTile.y % iStepY) == 0 );
    ptTile.x /= iStepX;
    ptTile.y /= iStepY;

    WCHAR chBase[MAX_PATH], *pBase = NULL;
    if (m_dstfrmt.pBase != NULL)
    {
        pBase = chBase;
        wcscpy_s(pBase, MAX_PATH, m_dstfrmt.pBase);
        PathStripPath(pBase);
    }

    vt::wstring_b<MAX_PATH> tilename;
    TileNameFormat(tilename, m_dstfrmt.eFrmt,
                   m_dstfrmt.eFile == eTileArchiveFileSystem ? m_dstfrmt.pBase : pBase,
                   m_dstfrmt.pExtension, m_uTileLevel, ptTile.x, ptTile.y);

    // use padding when constructing wdp tiles unless lossless compression
    // is being done - then padding isn't necessary
    if (_wcsicmp(m_dstfrmt.pExtension, L"wdp") == 0 &&
        m_dstfrmt.fQuality < 1.f )
    {
        vt::CRect rctCoreSrc = rctLayerDst;
        rctCoreSrc.InflateRect(m_dstfrmt.iOverlap, m_dstfrmt.iOverlap);
        if (m_dstfrmt.eFrmt != eTileHierarchyPhotosynth)
        {
            // Crop rect if not Photosynth cube.
            rctCoreSrc.IntersectRect(&rctCoreSrc, &m_rctSrc);
        }

        // Move rect to current tile coords
        rctCoreSrc = rctCoreSrc - ptSrc;

        CImgReaderWriter<CImg> src;
		imgSrc.Share(src);
        VT_HR_EXIT( ::SaveTileWithPadding(&m_params, tilename, rctCoreSrc,
                                           imgSrc.Rect(), m_pWriter, &src, 
                                           m_pArchive) );
    }
    else
    {
        CImgReaderWriter<CImg> src; imgSrc.Share(src);
        VT_HR_EXIT( ::SaveTile(m_hTransform, &m_params, tilename, imgSrc.Rect(),
                                m_pWriter, &src, m_pArchive) );
    }

    VT_HR_END()
}

HRESULT
CSaveTilesTransform::Clone(ITaskState **ppState)

{
    VT_HR_BEGIN()

    CSaveTilesTransform* pTrans = 
        VT_NOTHROWNEW CSaveTilesTransform(m_dstfrmt, &m_params, m_uTileLevel, 
                                m_rctSrc, m_ptDstOrigin, m_pWriter, m_pArchive);

   VT_HR_EXIT( CloneTaskState(ppState, pTrans) );
   ANALYZE_ASSUME( pTrans != NULL );
   VT_HR_EXIT( pTrans->SetSaveTilesParams() );

   VT_HR_END()
}

#endif
