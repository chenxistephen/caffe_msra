//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      General image format reader/writer
//
//  History:
//      2006/1/1-erudolph
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_io.h"
#include "vt_global.h"
#include "wicio.h"
#include "MFVideoSrc.h"

using namespace vt;

#if !defined(VT_WINRT)
void SetImageDebuggerImageNameToFilename(const vt::CImg& img)
{
	VT_HR_BEGIN();

	const vt::CParams* p = img.GetMetaData();

    if (p == NULL)
        return;

	wchar_t* fn;
	VT_HR_EXIT(p->GetByName(&fn, L"Filename"));
	
	string str;
	VT_HR_EXIT(VtWideCharToMultiByte(str, fn) > 0 ? S_OK : E_FAIL);

	char filename[MAX_PATH], extension[MAX_PATH];
	VT_HR_EXIT(_splitpath_s(str.get_buffer(), NULL, 0, NULL, 0, filename, 
		MAX_PATH, extension, MAX_PATH) == 0 ? S_OK : E_FAIL);

	string name(filename);
	VT_HR_EXIT(name.append(extension));
	
	VtImageDebuggerSetImageName(img, name.get_buffer());

	VT_HR_EXIT_LABEL();

	VT_ASSERT(!FAILED(hr));
}
#endif

const wchar_t* vt::VtIOErrorToString(HRESULT hr,  __out_ecount(numBufElem) wchar_t* buf, int numBufElem)
{
    wstring err;
    
	if (GetWICErrorString(hr, err) || GetMFErrorString(hr, err))
        wcscpy_s(buf, numBufElem, err.get_constbuffer());
    else
        VtErrorToString(hr, buf, numBufElem);

    return buf;
}

HRESULT vt::VtLoadImage(const wchar_t * pszFile, CImg &imgDst, 
                        bool bLoadMetadata)
{
    VT_HR_BEGIN();

	wchar_t filename[MAX_PATH], extension[MAX_PATH];
	VT_HR_EXIT(_wsplitpath_s(pszFile, NULL, 0, NULL, 0, filename, 
		MAX_PATH, extension, MAX_PATH) == 0 ? S_OK : E_FAIL);

    if (_wcsicmp(extension, L".vti") == 0)
    {
        VT_HR_EXIT(imgDst.Load(pszFile));
    }
    else
    {
        CWicReader fileSrc;
        VT_HR_EXIT( fileSrc.OpenFile( pszFile) );
        VT_HR_EXIT( fileSrc.GetImage( imgDst, NULL, bLoadMetadata ) );
    }

#if !defined(VT_WINRT)
	SetImageDebuggerImageNameToFilename(imgDst);
#endif

    VT_HR_END();
}

#if defined(VT_WINRT)
HRESULT vt::VtLoadImage(IStorageFile^ storageFile, CImg &imgDst, 
                        bool bLoadMetadata)
{
    VT_HR_BEGIN();

    CWicReader fileSrc;
    VT_HR_EXIT( fileSrc.OpenFile( storageFile ) );
    VT_HR_EXIT( fileSrc.GetImage( imgDst, NULL, bLoadMetadata ) );

    VT_HR_END();
}
#endif

HRESULT vt::VtLoadImage(const wchar_t * pszFile, IImageWriter* pWriter, 
                        const CPoint* ptDst, const CRect* pRect,
                        bool bLoadMetadata)
{
    VT_HR_BEGIN();

    CWicReader fileSrc;
    VT_HR_EXIT( fileSrc.OpenFile( pszFile) );
    VT_HR_EXIT( fileSrc.GetImage( pWriter, ptDst, pRect, bLoadMetadata ) );

	VT_HR_END();
}

HRESULT vt::VtSaveImage(const wchar_t * pszFile, const CImg &cImg, 
                        bool bSaveMetadata, float quality)
{
    VT_HR_BEGIN();

    if (quality < 0.f || quality > 1.f)
        return E_INVALIDARG;

    CParams p;
    CParamValue v;
    VT_HR_EXIT(v.Set(quality));
    VT_HR_EXIT(p.SetByName(L"ImageQuality", 0, v));

    VT_HR_EXIT(vt::VtSaveImage(pszFile, cImg, bSaveMetadata, &p, NULL));

    VT_HR_END();
}

#if defined(VT_WINRT)
HRESULT vt::VtSaveImage(IStorageFile^ storageFile, const CImg& imgSrc,
						bool bSaveMetadata, const CParams* pParams, CTaskProgress* pProgress)
{
	VT_HR_BEGIN();

    CWicWriter fileDst;
    VT_HR_EXIT( fileDst.OpenFile( storageFile) );
    VT_HR_EXIT( fileDst.SetImage( imgSrc, bSaveMetadata, pParams, pProgress) );
	VT_HR_EXIT( fileDst.CloseFile() );

    VT_HR_END();
}
#endif

HRESULT vt::VtSaveImage(const wchar_t * pszFile, const CImg &cImg, 
                        bool bSaveMetadata, const CParams* pParams,
                        CTaskProgress* pProgress)
{
    VT_HR_BEGIN();

	wchar_t filename[MAX_PATH], extension[MAX_PATH];
	VT_HR_EXIT(_wsplitpath_s(pszFile, NULL, 0, NULL, 0, filename, 
		MAX_PATH, extension, MAX_PATH) == 0 ? S_OK : E_FAIL);

    if (_wcsicmp(extension, L".vti") == 0)
    {
        VT_HR_EXIT(cImg.Save(pszFile));
    }
    else
    {
        CWicWriter fileDst;
        VT_HR_EXIT( fileDst.OpenFile( pszFile) );
        VT_HR_EXIT( fileDst.SetImage( cImg, bSaveMetadata, pParams, pProgress) );
    }

    VT_HR_END();
}

HRESULT vt::VtSaveImage(const wchar_t * pszFile, IImageReader* pReader, 
                        const CRect* pRect, 
                        bool bSaveMetadata,
                        const CParams* pParams, CTaskProgress* pProgress)
{
    VT_HR_BEGIN();

    CWicWriter fileDst;
    VT_HR_EXIT( fileDst.OpenFile( pszFile) );
    VT_HR_EXIT( fileDst.SetImage( pReader, pRect, bSaveMetadata, pParams, pProgress) );

    VT_HR_END();
}

HRESULT vt::VtSaveVideoImage(const wchar_t *pwcFilename,
			                 const CRGB32VideoImg &imgImage)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((imgImage.GetImg().IsValid()) ? S_OK : E_INVALIDARG);

    CRect rectCropRect = imgImage.GetVideoImgInfo().rectValidPixels;
    CRGBAImg imgCropped;
    VT_HR_EXIT(imgImage.GetImg().Share(imgCropped, &rectCropRect));
    VT_HR_EXIT(VtSaveImage(pwcFilename, imgCropped));

    VT_HR_END();
}

HRESULT vt::VtLoadVideoImage(const wchar_t *pwcFilename, CRGB32VideoImg &imgImage)
{
    VT_HR_BEGIN();

    CWicReader wicReader;
    VT_HR_EXIT(wicReader.OpenFile(pwcFilename));
	
    int iWidth  = wicReader.Width();
    int iHeight = wicReader.Height();
    int iVWidth  = iWidth + iWidth%2;
    int iVHeight = iHeight + iHeight%2;

    CVideoImgInfo vinfo;
    vinfo.rectValidPixels = CRect(0, 0, iWidth, iHeight);

    VT_HR_EXIT(imgImage.Create(iVWidth, iVHeight, vinfo));

    CImg imgValid;
    imgImage.GetImg().Share(imgValid, &vinfo.rectValidPixels);

    VT_HR_EXIT(wicReader.GetImage(imgValid));

    VT_HR_EXIT_LABEL();

    return hr;
}

HRESULT vt::VtSaveVideoImage(const wchar_t *pwcFilename, 
			                 const CNV12VideoImg &imgImage)
{
    VT_HR_BEGIN();
    CRGB32VideoImg imgTmp;
    VT_HR_EXIT(VtConvertVideoImage(imgTmp, imgImage));
    VT_HR_EXIT(VtSaveVideoImage(pwcFilename, imgTmp));
    VT_HR_END();
}

HRESULT vt::VtLoadVideoImage(const wchar_t *pwcFilename, CNV12VideoImg &imgImage)
{
    VT_HR_BEGIN();
    CRGB32VideoImg imgTmp;
    VT_HR_EXIT(VtLoadVideoImage(pwcFilename, imgTmp));
    VT_HR_EXIT(VtConvertVideoImage(imgImage, imgTmp));
    VT_HR_END();
}
