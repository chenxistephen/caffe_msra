//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Blob store implementation for Windows Store applications.
//
//  History:
//      2012/07/03-ericsto
//          Created
//
//------------------------------------------------------------------------
#include "stdafx.h"

#if defined(VT_WINRT)

#if WINAPI_FAMILY!=WINAPI_FAMILY_PHONE_APP // not for windows phone
#include <shlwapi.h>
#endif

#include "vt_utils.h"
#include "vt_blobstore.h"
#include "vt_winrtblobstore.h"

// macros from <intsafe.h>
#define LODWORD(_qw)    ((DWORD)(_qw))
#define HIDWORD(_qw)    ((DWORD)(((_qw) >> 32) & 0xffffffff))

using namespace vt;


typedef struct _Blob
{
    HANDLE hBlobFile;

    UINT uMaxTracts;
    UINT uBlobStoreIndex;
    FILE_BLOB_STORE BlobFile;

    _Blob() : hBlobFile(NULL), uBlobStoreIndex(UINT_MAX), uMaxTracts(0) {};
} Blob;

CWinRTBlobStore::CWinRTBlobStore() : m_tractSize(0)
{
    // Set temp directory as default backing file path.
    wstring_b<MAX_PATH> wstrTempPath;
	if (CSystem::GetTempPath(MAX_PATH, wstrTempPath.get_buffer()))
    {
        m_blobFilePath.uBlobStoreTracts = 0;
        m_blobFilePath.strBlobStorePath = wstrTempPath;
    }

    // From MSDN "File Buffering":
    // Therefore, in most situations, page-aligned memory will also be sector-aligned,
    // because the case where the sector size is larger than the page size is rare.
    SYSTEM_INFO si;
    vt::CSystem::GetSystemInfo(&si);
    m_tractSize = si.dwPageSize;
}

CWinRTBlobStore::~CWinRTBlobStore()
{
}

UINT
CWinRTBlobStore::GetTractSize()
{
    return m_tractSize;
}

HRESULT
CWinRTBlobStore::CreateBlob(LPVOID& pblob, LPCWSTR pwszName, UINT maxTracts)
{
    Blob* blob = NULL;

    VT_HR_BEGIN()

    if (m_tractSize == 0)
        VT_HR_EXIT( E_NOINIT );

    if (pwszName != NULL &&
        _wcsicmp(VtGetFileExt(pwszName), L".vti") != 0)
        VT_HR_EXIT( E_INVALIDARG );

    // See if there's room for a temporary blob.
    if (pwszName == NULL &&
        m_blobFilePath.uBlobStoreTracts > 0 &&
        m_blobFilePath.uBlobStoreTracts < maxTracts)
    {
		VT_HR_EXIT(HRESULT_FROM_WIN32(ERROR_DISK_FULL));
    }

    wstring_b<MAX_PATH> wstrBackingFile;
    if (pwszName == NULL)
    {
        // Get a temporary file name.
        if (!CSystem::GetTempFileName(m_blobFilePath.strBlobStorePath,
										L"IMG", 0, wstrBackingFile.get_buffer()))
        {
            VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));
        }
    }
    else
	{
		VT_HR_EXIT(wstrBackingFile.assign(pwszName));
	}

    // Create a temporary file.
    HANDLE hBlobFile = CSystem::CreateFile(wstrBackingFile, GENERIC_READ | GENERIC_WRITE,
                                    FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_ALWAYS,
                                    (pwszName == NULL ? FILE_FLAG_DELETE_ON_CLOSE : 0x0) |
                                    FILE_FLAG_NO_BUFFERING | FILE_FLAG_RANDOM_ACCESS,
                                    NULL);
    if (hBlobFile == INVALID_HANDLE_VALUE)
    {
        VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));
    }

	// Create the blob.
    blob = VT_NOTHROWNEW Blob;
    VT_PTR_OOM_EXIT(blob);

    blob->uMaxTracts = maxTracts;
    blob->hBlobFile = hBlobFile;

    if (pwszName == NULL)
        blob->uBlobStoreIndex = 0;
    else
        VT_HR_EXIT( blob->BlobFile.strBlobStorePath.assign(pwszName) );

    pblob = blob;

    VT_HR_EXIT_LABEL()

    if (hr != S_OK)
    {
        if (blob != NULL)
        {
            if (blob->hBlobFile != NULL)
                CloseHandle(blob->hBlobFile);
            delete blob;
        }
    }
        
    return hr;
}

UINT
CWinRTBlobStore::GetBlobSize(LPVOID pblob)
{
    Blob *blob = (Blob *) pblob;

    return blob == NULL ? 0 : blob->BlobFile.uBlobStoreTracts;
}

HRESULT
CWinRTBlobStore::ExtendBlob(LPVOID pblob, UINT numTracts)
{
    VT_HR_BEGIN()

    Blob *blob = (Blob *) pblob;

    if (blob == NULL || blob->hBlobFile == NULL)
        VT_HR_EXIT( E_NOINIT );

    if (numTracts == 0)
        VT_HR_EXIT( E_INVALIDARG );

    if (blob->uMaxTracts > 0 &&
        blob->uMaxTracts < blob->BlobFile.uBlobStoreTracts + numTracts)
        VT_HR_EXIT( E_ACCESSDENIED );

    if (blob->uBlobStoreIndex != UINT_MAX)
    {
        // Check if new tracts fit in blob path.
        if (m_blobFilePath.uBlobStoreTracts > 0 &&
            m_blobFilePath.uBlobStoreTracts < numTracts)
            VT_HR_EXIT( HRESULT_FROM_WIN32( ERROR_DISK_FULL ) );
    }

    LARGE_INTEGER liFileSize;
    liFileSize.QuadPart = (LONGLONG) m_tractSize *
        (blob->BlobFile.uBlobStoreTracts + numTracts);

    if (SetFilePointerEx(blob->hBlobFile, liFileSize, NULL, FILE_BEGIN) == INVALID_SET_FILE_POINTER ||
        !SetEndOfFile(blob->hBlobFile))
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

    if (blob->uBlobStoreIndex != UINT_MAX)
    {
        // Charge new tracts to blob path.
        if (m_blobFilePath.uBlobStoreTracts > 0)
            m_blobFilePath.uBlobStoreTracts -= numTracts;
    }

    blob->BlobFile.uBlobStoreTracts += numTracts;

    VT_HR_END()
}

HRESULT
CWinRTBlobStore::OpenBlob(LPVOID& pblob, UINT& numTracts, LPCWSTR pwszName,
                         bool bReadOnly)
{
    Blob* blob = NULL;

    VT_HR_BEGIN()

    if (pwszName == NULL ||
        _wcsicmp(VtGetFileExt(pwszName), L".vti") != 0)
        VT_HR_EXIT( E_INVALIDARG );

    blob = VT_NOTHROWNEW Blob;
    VT_PTR_OOM_EXIT(blob);

    // Open an existing file.
    HANDLE hBlobFile = CSystem::CreateFile(pwszName,
                                  GENERIC_READ | (bReadOnly ? 0x0 : GENERIC_WRITE),
                                  FILE_SHARE_READ | (bReadOnly ? 0x0 : FILE_SHARE_WRITE),
                                  NULL, OPEN_EXISTING,
                                  FILE_FLAG_NO_BUFFERING | FILE_FLAG_RANDOM_ACCESS,
                                  NULL);
    if (hBlobFile == INVALID_HANDLE_VALUE)
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

    LARGE_INTEGER liFileSize;
    if (!CSystem::GetFileSizeEx(hBlobFile, &liFileSize))
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

    VT_HR_EXIT( blob->BlobFile.strBlobStorePath.assign(pwszName) );

    blob->BlobFile.uBlobStoreTracts = (UINT) (liFileSize.QuadPart / m_tractSize);

    blob->hBlobFile = hBlobFile;

    pblob = blob;

    numTracts = blob->BlobFile.uBlobStoreTracts;

    VT_HR_EXIT_LABEL()

    if (hr != S_OK)
    {
        if (blob != NULL)
        {
            if (blob->hBlobFile != NULL)
                CloseHandle(blob->hBlobFile);
            delete blob;
        }
    }
        
    return hr;
}

HRESULT
CWinRTBlobStore::ReadBlob(LPVOID pblob, UINT tractNum, UINT numTracts,
                         LPBYTE tractBuf, IBlobCallback* cbfn)
{
    VT_HR_BEGIN()

    Blob *blob = (Blob *) pblob;

    if (blob == NULL || blob->hBlobFile == NULL)
        VT_HR_EXIT( E_NOINIT );

    if (blob->BlobFile.uBlobStoreTracts < tractNum + numTracts ||
        numTracts == 0 || tractBuf == NULL)
        VT_HR_EXIT( E_INVALIDARG );

    // Check tract buffer is multiple of tract size.
    UINT_PTR uTractMask = m_tractSize - 1;
    if ((UINT_PTR) tractBuf != (((UINT_PTR) tractBuf +
        uTractMask) & ~uTractMask))
        VT_HR_EXIT( E_POINTER );

	OVERLAPPED overlapped = { 0 };
    UINT64 uiOffset = (UINT64) tractNum * m_tractSize;
    overlapped.Offset = LODWORD(uiOffset);
    overlapped.OffsetHigh = HIDWORD(uiOffset);

    // Read synchronously.
    if (!ReadFile(blob->hBlobFile, tractBuf, numTracts * m_tractSize, NULL, &overlapped))
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() )  );

    // If a callback was supplied, call it now.
    if (cbfn != NULL)
    {
		cbfn->Callback(hr);
	}

    VT_HR_END()
}

HRESULT
CWinRTBlobStore::WriteBlob(LPVOID pblob, UINT tractNum, UINT numTracts,
                          const LPBYTE tractBuf, IBlobCallback* cbfn)
{
    VT_HR_BEGIN()

    Blob *blob = (Blob *) pblob;

    if (blob == NULL || blob->hBlobFile == NULL)
        VT_HR_EXIT( E_NOINIT );

    if (blob->BlobFile.uBlobStoreTracts < tractNum + numTracts ||
        numTracts == 0 || tractBuf == NULL)
        VT_HR_EXIT( E_INVALIDARG );

    // Check tract buffer is multiple of tract size.
    UINT_PTR uTractMask = m_tractSize - 1;
    if ((UINT_PTR) tractBuf != (((UINT_PTR) tractBuf +
        uTractMask) & ~uTractMask))
        VT_HR_EXIT( E_POINTER );

	OVERLAPPED overlapped = { 0 };
    UINT64 uiOffset = (UINT64) tractNum * m_tractSize;
    overlapped.Offset = LODWORD(uiOffset);
    overlapped.OffsetHigh = HIDWORD(uiOffset);

    // Write synchronously.
    if (!WriteFile(blob->hBlobFile, tractBuf, numTracts * m_tractSize, NULL, &overlapped))
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() )  );

    // If a callback was supplied, call it now.
    if (cbfn != NULL)
    {
		cbfn->Callback(hr);
	}

    VT_HR_END()
}

HRESULT
CWinRTBlobStore::CloseBlob(LPVOID pblob)
{
    VT_HR_BEGIN()

    Blob *blob = (Blob *) pblob;

    if (blob == NULL)
        VT_HR_EXIT( E_INVALIDARG );

    if (!CloseHandle(blob->hBlobFile))
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

    if (blob->uBlobStoreIndex != UINT_MAX)
    {
        // Credit blob tracts back to blob path.
        if (m_blobFilePath.uBlobStoreTracts > 0)
            m_blobFilePath.uBlobStoreTracts += blob->BlobFile.uBlobStoreTracts;
    }

    delete blob;

    VT_HR_END()
}

HRESULT
CWinRTBlobStore::DeleteBlob(LPCWSTR pwszName)
{
    VT_HR_BEGIN()

    if (pwszName == NULL ||
        _wcsicmp(VtGetFileExt(pwszName), L".vti") != 0)
        VT_HR_EXIT( E_INVALIDARG );

    if (!DeleteFile(pwszName))
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

    VT_HR_END()
}

#endif
