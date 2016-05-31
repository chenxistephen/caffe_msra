//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Helper interface for archiving images.
//
//  History:
//      2009/5/20-v-hogood
//          Created
//
//------------------------------------------------------------------------
#pragma once

#include "vtcommon.h"
#include <ObjBase.h>

namespace vt {

//+---------------------------------------------------------------------------
// 
// Generic interfaces supported by various archive objects 
//
//----------------------------------------------------------------------------
class IArchiveWriter
{
public:
    virtual ~IArchiveWriter() {};

    // OpenArchive opens the archive for writing.
    virtual HRESULT OpenArchive(LPCWSTR pwszArchive, size_t cbPart) = 0;

    // Add a file to the archive.
    virtual HRESULT AddFile(LPCWSTR pwszFile, IStream& Stream) = 0;
    virtual HRESULT AddFile(LPCWSTR pwszFile, const LPBYTE pMemory, size_t cbSize) = 0;

	// Return the last part number.
	virtual int GetLastPartNumber() = 0;

    // Close the archive and return the last part number.
    virtual int CloseArchive() = 0;
};

class IArchiveReader
{
public:
    virtual ~IArchiveReader() {};

    // OpenArchive opens the archive for reading.
    virtual HRESULT OpenArchive(LPCWSTR pwszArchive, int iMaxPart) = 0;

    // Add a file to the archive.
    virtual HRESULT ReadFile(LPCWSTR pwszFile, IStream& Stream) = 0;

    // Close the archiver.
    virtual void CloseArchive() = 0;
};

};
