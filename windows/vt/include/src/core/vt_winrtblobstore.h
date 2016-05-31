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
#pragma once

#if defined(VT_WINRT)

#include "vtcommon.h"
#include "vt_blobstore.h"

namespace vt {

struct FILE_BLOB_STORE
{
	UINT        uBlobStoreTracts;
	vt::wstring strBlobStorePath;

    FILE_BLOB_STORE() : uBlobStoreTracts(0) {};
};

class CWinRTBlobStore : public vt::IBlobStore
{
public:
	CWinRTBlobStore();
    virtual ~CWinRTBlobStore();

	//+-------------------------------------------------------------------------
    //
    // From IBlobStore:
    //

    // Returns the tract size for this type of blob store.
    virtual UINT GetTractSize();

    // CreateBlob creates a blob for reading and writing.
    // Initial size is 0 tracts; maximum is maxTracts (0 == unlimited).
    // Blob is delete on close if not named.
    virtual HRESULT CreateBlob(OUT LPVOID& pblob, LPCWSTR pwszName = NULL,
                               UINT maxTracts = 0);

    // Returns the number of tracts in this blob.
    virtual UINT GetBlobSize(LPVOID pblob);

    // Extend an existing blob by a number of tracts.
    virtual HRESULT ExtendBlob(LPVOID pblob, UINT numTracts);

    // Open an existing named blob for reading and writing.
    virtual HRESULT OpenBlob(OUT LPVOID& pblob, OUT UINT& numTracts,
                             LPCWSTR pwszName, bool bReadOnly = false);

    // Read tracts from the blob.
    virtual HRESULT ReadBlob(LPVOID pblob, UINT tractNum, UINT numTracts,
                             LPBYTE tractBuf, IBlobCallback* cbfn = NULL);

    // Write tracts to the blob.
    virtual HRESULT WriteBlob(LPVOID pblob, UINT tractNum, UINT numTracts,
                              const LPBYTE tractBuf, IBlobCallback* cbfn = NULL);

    // Close the blob (and delete if not named).
    virtual HRESULT CloseBlob(LPVOID pblob);

    // Delete the named blob.
    virtual HRESULT DeleteBlob(LPCWSTR pwszName);

private:
    UINT m_tractSize;
    FILE_BLOB_STORE m_blobFilePath;
};

};

#endif
