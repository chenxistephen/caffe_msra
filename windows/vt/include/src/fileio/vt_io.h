//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      General image format reader/writer
//
//  History:
//      2006/1/1-erudolph
//            Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vtcore.h"

namespace vt {

//+---------------------------------------------------------------------------
// 
// File loading/saving functions
//
//----------------------------------------------------------------------------

/// \ingroup error
/// <summary> Return i/o HRESULT as string </summary>
/// <param name="hr"> HRESULT from a VisionTools i/o operation </param>
/// <param name="buf"> Buffer to write the error text to </param>
/// <param name="numBufElem"> Buffer size, in characters </param>
/// <returns> Pointer to the beginning of the text buffer </returns>
const wchar_t* VtIOErrorToString(
    HRESULT hr, 
    __out_ecount(numBufElem) wchar_t* buf, 
    int numBufElem);

/// \ingroup loadsave
/// <summary> Load image from file </summary>
/// <param name="pszFile"> File name. </param>
/// <param name="imgDst"> Destination image. </param>
/// <param name="bLoadMetadata"> If true, metadata will be loaded. </param>
/// <returns> 
///		- S_OK on success
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///	The pixel format of `imgDst` is determined using \ref stdconvrules.
HRESULT VtLoadImage(const wchar_t * pszFile, CImg &imgDst, 
                    bool bLoadMetadata = true);

#if defined(MAKE_DOXYGEN_WORK) || defined(VT_WINRT)
/// \ingroup loadsave
/// <summary> Load image from file </summary>
/// <param name="storageFile"> Storage file. </param>
/// <param name="imgDst"> Destination image. </param>
/// <param name="bLoadMetadata"> If true, metadata will be loaded. </param>
/// <returns> 
///		- S_OK on success
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///	- Available only in WinRT configurations.
///	- The pixel format of `imgDst` is determined using \ref stdconvrules.
HRESULT VtLoadImage(Windows::Storage::IStorageFile^ storageFile, CImg &imgDst, 
                    bool bLoadMetadata = true);
#endif

HRESULT VtLoadImage(const wchar_t * pszFile, IImageWriter* pWriter, 
					const CPoint* ptDst = NULL, const CRect* pRectSrc = NULL,
                    bool bLoadMetadata = true);

/// \ingroup loadsave
/// <summary> Save image to file </summary>
/// <param name="pszFile"> File name. </param>
/// <param name="imgSrc"> Source image. </param>
/// <param name="bSaveMetadata"> If true, metadata will be saved. </param>
/// <param name="pParams"> (internal) </param>
/// <param name="pProgress"> (internal) </param>
/// <returns> 
///		- S_OK on success
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
HRESULT VtSaveImage(const wchar_t * pszFile, const CImg &imgSrc, 
                    bool bSaveMetadata = true,
                    const CParams* pParams = NULL,
                    CTaskProgress* pProgress = NULL);

#if defined(MAKE_DOXYGEN_WORK) || defined(VT_WINRT)
/// \ingroup loadsave
/// <summary> Save image to file </summary>
/// <param name="storageFile"> Storage file. </param>
/// <param name="imgSrc"> Source image. </param>
/// <param name="bSaveMetadata"> If true, metadata will be saved. </param>
/// <param name="pParams"> (internal) </param>
/// <param name="pProgress"> (internal) </param>
/// <returns> 
///		- S_OK on success
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///	Available only in WinRT configurations.

HRESULT VtSaveImage(Windows::Storage::IStorageFile^ storageFile, const CImg& imgSrc,
                    bool bSaveMetadata = true,
                    const CParams* pParams = NULL,
                    CTaskProgress* pProgress = NULL);
#endif

/// \ingroup loadsave
/// <summary> Save image to file, with quality setting </summary>
/// <param name="pszFile"> File name. </param>
/// <param name="imgSrc"> Source image. </param>
/// <param name="bSaveMetadata"> If true, metadata will be saved. </param>
/// <param name="quality"> Image quality, between 0 and 1. 
/// Ignored if output format doesn't support compression. </param>
/// <returns> 
///		- S_OK on success
///		- E_INVALIDARG if quality is out of range
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
HRESULT VtSaveImage(const wchar_t * pszFile, const CImg &imgSrc, 
                    bool bSaveMetadata, float quality);

HRESULT VtSaveImage(const wchar_t * pszFile, IImageReader* pReader, 
					const CRect* pRect = NULL, 
					bool bSaveMetadata = true,
                    const CParams* pParams = NULL,
                    CTaskProgress* pProgress = NULL );

/// \ingroup loadsave
/// <summary> Load video image from file into RGB32 format</summary>
/// <param name="pwcFilename"> File name. </param>
/// <param name="imgImage"> Destination image. </param>
/// <returns> 
///		- S_OK on success
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///		- imgImage will always have even width/height
///		- CVideoImgInfo::rectValidPixels contains the image size in the file
HRESULT VtLoadVideoImage(const wchar_t *pwcFilename, CRGB32VideoImg &imgImage);

/// \ingroup loadsave
/// <summary> Load video image from file into NV12 format</summary>
/// <param name="pwcFilename"> File name. </param>
/// <param name="imgImage"> Destination image. </param>
/// <returns> 
///		- S_OK on success
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///		- imgImage will always have even width/height
///		- CVideoImgInfo::rectValidPixels contains the image size in the file
HRESULT VtLoadVideoImage(const wchar_t *pwcFilename, CNV12VideoImg &imgImage);

/// \ingroup loadsave
/// <summary> Save RGB32 format video image to file</summary>
/// <param name="pwcFilename"> File name. </param>
/// <param name="imgImage"> Destination image. </param>
/// <returns> 
///		- S_OK on success
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///     - Only CVideoImgInfo::rectValidPixels are output
HRESULT VtSaveVideoImage(const wchar_t *pwcFilename, 
                         const CRGB32VideoImg &imgImage);

/// \ingroup loadsave
/// <summary> Save NV12 format video image to file</summary>
/// <param name="pwcFilename"> File name. </param>
/// <param name="imgImage"> Destination image. </param>
/// <returns> 
///		- S_OK on success
///		- Use VtIOErrorToString() to get extended error information.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
///     - Only CVideoImgInfo::rectValidPixels are output
HRESULT VtSaveVideoImage(const wchar_t *pwcFilename,
                         const CNV12VideoImg &imgImage);

/// \ingroup loadsave
/// <summary> Returns true if Media Foundation is available (and
/// VtCreateVideoSrc or  VtCreateVideoDst has any chance to work). </summary>
bool VtIsMediaFoundationAvailable();

/// \ingroup loadsave
/// <summary> Creates a video source and (optionally) opens a video file.
/// Returns E_NOINTERFACE if Media Foundation is unavailable. </summary>
/// <param name="ppVideoSrc"> Pointer to a pointer that receives the created object </param>
/// <param name="filename"> The video file to open. If NULL, no file is opened </param>
HRESULT VtCreateVideoSrc(IVideoSrc** ppVideoSrc, const WCHAR* filename = NULL);

#if defined(MAKE_DOXYGEN_WORK) || defined(VT_WINRT)
/// \ingroup loadsave
/// <summary> Creates a video source and (optionally) opens a video file.
/// Returns E_NOINTERFACE if Media Foundation is unavailable. </summary>
/// <param name="ppVideoSrc"> Pointer to a pointer that receives the created object </param>
/// <param name="storageFile"> The video file to open. If nullptr, no file is opened </param>
/// <DL><DT> Remarks: </DT></DL>
///	Available only in WinRT configurations.
HRESULT VtCreateVideoSrc(IVideoSrc** ppVideoSrc, Windows::Storage::IStorageFile^ storageFile);
#endif

/// \ingroup loadsave
/// <summary> Creates a video destination.
/// Returns E_NOINTERFACE if Media Foundation is unavailable. </summary>
/// <param name="ppVideoDst"> Pointer to a pointer that receives the created object </param>
HRESULT VtCreateVideoDst(IVideoDst** ppVideoDst);

};
