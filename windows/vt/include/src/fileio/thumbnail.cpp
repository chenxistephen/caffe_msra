//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//     Thumbnail reading functions
//
//  History:
//      2009/10/16-ericsto
//            Created
//
//------------------------------------------------------------------------
#include "stdafx.h"

#if !defined(VT_WINRT)

#pragma warning( push )
#pragma warning( disable : 6387 )
#include <shlobj.h>
#include <shobjidl.h>
#pragma warning( pop ) 

#include "vt_global.h"
#include "vt_thumbnail.h"

using namespace vt;

CThumbnailReader::CThumbnailReader()
{
	m_bInitialized = false;
	m_createItemFunc = NULL;
	m_pThumbnailCache = NULL;
}

CThumbnailReader::~CThumbnailReader()
{
	if (m_pThumbnailCache)
	{
		m_pThumbnailCache->Release();
		m_pThumbnailCache = NULL;
	}
}

HRESULT
CThumbnailReader::ReadThumbnail(const WCHAR *pszFile, int desiredSize, CRGBAImg &cImg)
{
	HRESULT hr = S_OK;

	// Make sure the class is initialized.
	Initialize();

	// First, try to get a thumbnail from the shell's thumbnail cache.
	HBITMAP hBitmap = NULL;
	if (m_createItemFunc != NULL && m_pThumbnailCache != NULL)
	{
		hr = GetThumbnailFromCache(pszFile, desiredSize, hBitmap);
	}

	// Alternatively, try to extract a thumbnail from the image itself.
	if (hBitmap == NULL)
	{
		VT_HR_EXIT(GetThumbnailFromImage(pszFile, desiredSize, hBitmap));
	}

	// Allocate space for the pixels in the CImg.
	DIBSECTION dibSection;
	VT_HR_EXIT(GetObject(hBitmap, sizeof(DIBSECTION), &dibSection) == 0 ? E_FAIL : S_OK);
	VT_HR_EXIT(cImg.Create(dibSection.dsBm.bmWidth, dibSection.dsBm.bmHeight, align16Byte));

	// Convert the HBitmap to CImg pixels.
	if (dibSection.dsBm.bmBitsPixel != 32)
	{
		// The bitmap doesn't have 32 bits per pixel, so we have to convert.
		// Initialize a bitmap info header with the desired pixel format.
		// Use a negative height to invert the image (the bitmap produced
		// by IExtractImage::Extract has bottom-left origin, but we want
		// top-left origin).
		int width = dibSection.dsBm.bmWidth;
		int height = dibSection.dsBm.bmHeight;
		BITMAPINFOHEADER bitmapInfoHeader;
		bitmapInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
		bitmapInfoHeader.biWidth = width;
		bitmapInfoHeader.biHeight = -height;
		bitmapInfoHeader.biPlanes = 1;
		bitmapInfoHeader.biBitCount = 32;
		bitmapInfoHeader.biCompression = BI_RGB;
		bitmapInfoHeader.biSizeImage = 0;
		bitmapInfoHeader.biXPelsPerMeter = 0;
		bitmapInfoHeader.biYPelsPerMeter = 0;
		bitmapInfoHeader.biClrUsed = 0;
		bitmapInfoHeader.biClrImportant = 0;

		// Get the pixels from the bitmap.
		HDC hdc = GetDC(NULL);
		int result = GetDIBits(hdc, hBitmap, 0, (UINT)height, cImg.BytePtr(),
			(BITMAPINFO *)&bitmapInfoHeader, DIB_RGB_COLORS);
		ReleaseDC(NULL, hdc);
		VT_HR_EXIT(result == 0 ? E_FAIL : S_OK);

		// The alpha channel contains all zeros, so we have to fill it in.
		for (int j = 0; j < height; j++)
		{
			BYTE* ptr = cImg.BytePtr(j) + 3;
			for (int i = 0; i < width; i++)
			{
				*ptr = 255;
				ptr += 4;
			}
		}
	}
	else
	{
		VT_ASSERT(dibSection.dsBm.bmBitsPixel == 32);
		CRGBAImg tmpImg;
		// wrap 32bpp bitmap
		VT_HR_EXIT(tmpImg.Create((Byte*)dibSection.dsBm.bmBits, 
			dibSection.dsBm.bmWidth, dibSection.dsBm.bmHeight, 
			dibSection.dsBm.bmWidth*4));
		VT_HR_EXIT(tmpImg.CopyTo(cImg));
	}
	DeleteObject(hBitmap);

Exit:
	return hr;
}

void
CThumbnailReader::Initialize()
{
	// Only initialize this instance once.
	if (!m_bInitialized)
	{
		m_bInitialized = true;

		// Find the SHCreateItemFromParsingName function if possible (Vista or later).
		HMODULE module = GetModuleHandle(_T("Shell32.dll"));
		if (module != NULL)
		{
			m_createItemFunc = (SHCreateItemFromParsingNameFunc)GetProcAddress(module, "SHCreateItemFromParsingName");
		}

		// Get an implementation of IThumbnailCache if possible (Vista or later).
		CSystem::CoCreateInstance(CLSID_LocalThumbnailCache, NULL, CLSCTX_INPROC, IID_IThumbnailCache, (void**)&m_pThumbnailCache);
	}
}

HRESULT
CThumbnailReader::GetThumbnailFromCache(const WCHAR *pszFile, int desiredSize, HBITMAP &hBitmap)
{
	VT_ASSERT(m_bInitialized);
	VT_ASSERT(m_createItemFunc != NULL);
	VT_ASSERT(m_pThumbnailCache != NULL);
	ANALYZE_ASSUME(m_pThumbnailCache);

	HRESULT hr = S_OK;
	CComPtr<IShellItem> pShellItem;
	CComPtr<ISharedBitmap> pSharedBitmap;

	// Create a shell item from the filename.
	VT_HR_EXIT(m_createItemFunc(pszFile, NULL, IID_IShellItem, (void**)&pShellItem));

	// Try to get a shared bitmap from the thumbnail cache.  Extract it from the original image if necessary.
	WTS_CACHEFLAGS cacheFlags = (WTS_CACHEFLAGS)0;
	VT_HR_EXIT(m_pThumbnailCache->GetThumbnail(pShellItem, desiredSize, WTS_EXTRACT, &pSharedBitmap, &cacheFlags, NULL));

	// Detach an unshared bitmap.
	VT_HR_EXIT(pSharedBitmap->Detach(&hBitmap));

Exit:
	return hr;
}

HRESULT
CThumbnailReader::GetThumbnailFromImage(const WCHAR *pszFile, int desiredSize, HBITMAP &hBitmap)
{
	HRESULT hr = S_OK;
    CComPtr<IShellFolder> pShellFolder;
	CComPtr<IExtractImage> pExtractImage;
	CSize size(desiredSize, desiredSize);

	// Create a shell ID list from the filename.
	LPITEMIDLIST pidList = NULL;
	VT_HR_EXIT(SHILCreateFromPath(pszFile, &pidList, NULL));

	// Get the parent shell folder and relative ID list to the child item.
    LPCITEMIDLIST pidListChild;
    VT_HR_EXIT(SHBindToParent(pidList, IID_IShellFolder, (void**)&pShellFolder, &pidListChild));

	// Get an implementation of IExtractImage for the child item.
    VT_HR_EXIT(pShellFolder->GetUIObjectOf(NULL, 1, &pidListChild, IID_IExtractImage, NULL, (void**)&pExtractImage));

	// Locate and extract the thumbnail image.
	DWORD dwPriority = 0;
	DWORD dwFlags = IEIFLAG_SCREEN | IEIFLAG_ORIGSIZE | IEIFLAG_QUALITY;
	WCHAR buffer[_MAX_PATH];
	DWORD bitsPerPixel = 24;

	VT_HR_EXIT(pExtractImage->GetLocation(buffer, _MAX_PATH, &dwPriority, &size, bitsPerPixel, &dwFlags));
	VT_HR_EXIT(pExtractImage->Extract(&hBitmap));

Exit:
	ILFree(pidList);
	// According to the documentation for SHBindToParent, it's not neccessary to free pidListChild.
	return hr;
}

#endif
