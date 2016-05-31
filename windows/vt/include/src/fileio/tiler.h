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
#pragma once

#if !defined(VT_WINRT)

#include "vtcommon.h"

namespace vt {

enum eTileHierarchyFormat
{
    eTileHierarchyFlat         =  0,
    eTileHierarchyFolderPerCol =  1,
    eTileHierarchySeadragon    =  2,
    eTileHierarchyPhotosynth   =  3,
    eTileHierarchyUnspecified  = -1
};

enum eTileArchiveFormat
{
    eTileArchiveFileSystem =  0,
    eTileArchiveZip        =  1
};

enum eTileDstProfile
{
    eTileProfileSRGB        =  0,
    eTileProfileProPhoto    =  1,
    eTileProfileAdobe       =  2,
    eTileProfileApple       =  3,
    eTileProfileUnspecified = -1
};

struct TILE_SRC_FRMT
{
    eTileHierarchyFormat  eFrmt;
    const WCHAR* pBase;
    const WCHAR* pExtension;
	int iTileWidth;
	int iTileHeight;
	int iOverlap;

    TILE_SRC_FRMT() : eFrmt(eTileHierarchyFolderPerCol),
                      pBase(NULL), pExtension(L"wdp"),
                      iTileWidth(256), iTileHeight(256),
                      iOverlap(0)
	{}
};

struct TILE_DST_FRMT
{
    eTileHierarchyFormat eFrmt;
    eTileArchiveFormat   eFile;
    eTileDstProfile      eProf;
    const WCHAR* pBase;
    const WCHAR* pExtension;
    float        fQuality;
    int          iImgType;
	int          iTileWidth;
	int          iTileHeight;
	int          iOverlap;

    TILE_DST_FRMT() : eFrmt(eTileHierarchyFolderPerCol),
                      pBase(NULL), pExtension(L"wdp"), 
		              fQuality(0.9f), iImgType(OBJ_RGBIMG), 
		              iTileWidth(256), iTileHeight(256),
		              iOverlap(0), eFile(eTileArchiveFileSystem),
                      eProf(eTileProfileUnspecified)
	{}
};

HRESULT
LoadTiles(const TILE_SRC_FRMT& srcfrmt, int level, const vt::CRect& rct, 
          OUT vt::CImg& imgBuff, bool bIgnoreMissing=false);

HRESULT
SaveTiles(const TILE_DST_FRMT& dstfrmt, const CParams* pEncParams, 
		  UINT uTileLevel, const CPoint& ptSrc, const CRect& rctDst, 
		  IVTImageWriter* pWriter, IImageReader* pSrc,
		  IArchiveWriter* pArchive=NULL,
		  CTaskProgress* pProgress=NULL);

// call format on wstring(dst) - output is the tile file name
template <class T> HRESULT
TileNameFormat(T& name, eTileHierarchyFormat eFrmt, const WCHAR* base,
               const WCHAR* ext, int l, int x, int y);

HRESULT
DetermineTileSetDim(const TILE_SRC_FRMT& srcfrmt, int level, 
                    OUT int& width, OUT int& height, OUT int& bands);

//+-----------------------------------------------------------------------------
//
// Class: CSaveTilesTransform
// 
// Synposis: IImageTransform implementation of SaveTiles()
// 
//------------------------------------------------------------------------------
#define HDPHOTO_PAD 8

class CSaveTilesTransform : public CImageTransformUnaryGeo<CSaveTilesTransform, true>
{
public:
	CSaveTilesTransform() : m_pWriter(NULL), m_hTransform(NULL)
	{}

	CSaveTilesTransform(const TILE_DST_FRMT& dstfrmt, const CParams* pEncParams, 
						UINT uTileLevel, const CRect& rctSrc, 
						const vt::CPoint& ptDstOffset,
						IVTImageWriter* pWriter, IArchiveWriter* pArchive=NULL) : 
		m_pWriter(NULL), m_hTransform(NULL)
	{
		m_dstfrmt    = dstfrmt;
		if (pEncParams != NULL)
		{
			m_params = *pEncParams;
		}
		m_uTileLevel  = uTileLevel;
		m_rctSrc      = rctSrc;
		m_ptDstOrigin = ptDstOffset;
		if ( pWriter != NULL )
		{
			pWriter->Clone(&m_pWriter);
		}
		m_pArchive   = dstfrmt.eFile == eTileArchiveZip? pArchive : NULL;
	}

	// need to call Begin() before pushing the transform - this will create
	// the necessary directory structure on disk
	HRESULT Begin();  

	~CSaveTilesTransform()
	{
		if (m_hTransform != NULL)
		{	
			DeleteColorTransform(m_hTransform);
		}
		delete m_pWriter;
	}

public:
	virtual vt::CRect GetRequiredSrcRect(const vt::CRect& rctDst)
	{
		vt::CRect rctPad = rctDst;

		// Add overlap.
		rctPad.InflateRect(m_dstfrmt.iOverlap, m_dstfrmt.iOverlap);
		if (_wcsicmp(m_dstfrmt.pExtension, L"wdp") == 0 &&
			m_dstfrmt.fQuality < 1.f )
		{
			// Add HD Photo padding.
			rctPad.InflateRect(HDPHOTO_PAD, HDPHOTO_PAD);
		}

		if (m_dstfrmt.eFrmt != eTileHierarchyPhotosynth)
		{
			// Crop rect if not Photosynth cube.
			rctPad.IntersectRect(&rctPad, &m_rctSrc);
		}

		return rctPad;
	}

	virtual vt::CRect GetAffectedDstRect(const vt::CRect& rctSrc)
	{ return rctSrc; }

	virtual vt::CRect GetResultingDstRect(const vt::CRect& rctSrc)
	{ return rctSrc; }

	virtual HRESULT Transform(OUT CImg* pimgDst, 
							  IN  const CRect& rctLayerDst,
							  IN  const CImg& imgSrc,
							  IN  const vt::CPoint& ptSrc);

    virtual HRESULT SetSaveTilesParams();

    virtual HRESULT Clone(ITaskState **ppState);

protected:
	TILE_DST_FRMT   m_dstfrmt;
	CParams         m_params;
	UINT            m_uTileLevel;
	vt::CRect       m_rctSrc;
	vt::CPoint      m_ptDstOrigin;
	IVTImageWriter* m_pWriter;
	IArchiveWriter* m_pArchive;
	HTRANSFORM      m_hTransform;
};

};

#endif
