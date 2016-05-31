//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//     Thumbnail reader
//
//  History:
//      2009/10/16-ericsto
//            Created
//
//------------------------------------------------------------------------

#pragma once

#if !defined(VT_WINRT)

#include "vtcommon.h"

namespace vt {

//+----------------------------------------------------------------------------
//
// Class CThumbnailReader
// 
// Synposis:
//     Thumbnail reading routines
// 
//-----------------------------------------------------------------------------

class CThumbnailReader
{
public:
    CThumbnailReader();
    ~CThumbnailReader();

    // ReadThumbnail extracts a thumbnail from the shell's thumbnail cache
	// or from the image's metadata.
    HRESULT ReadThumbnail(const WCHAR *pszFile, int desiredSize, CRGBAImg &cImg);

private:
	void Initialize();
	HRESULT GetThumbnailFromCache(const WCHAR *pszFile, int desiredSize, HBITMAP &hBitmap);
	HRESULT GetThumbnailFromImage(const WCHAR *pszFile, int desiredSize, HBITMAP &hBitmap);

	typedef HRESULT (STDAPICALLTYPE *SHCreateItemFromParsingNameFunc)(__in PCWSTR pszPath, __in_opt IBindCtx *pbc, __in REFIID riid, __deref_out void **ppv);

	bool m_bInitialized;
	SHCreateItemFromParsingNameFunc m_createItemFunc;
	IThumbnailCache *m_pThumbnailCache;
};

};

#endif
