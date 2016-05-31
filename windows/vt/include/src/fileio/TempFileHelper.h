//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Helper class for creating and copying from a temporary file.
//
//  History:
//      2012/07/25-ericsto
//          Created
//
//------------------------------------------------------------------------

#pragma once

#if defined(VT_WINRT)

using namespace Windows::Storage;

namespace vt
{
	class CTempFileHelper
	{
	public:
		CTempFileHelper()
			: m_storageFile(nullptr)
			, m_tempStorageFile(nullptr)
			, m_fileAccessMode(FileAccessMode::Read)
		{
		}

		~CTempFileHelper()
		{
			Finalize();
		}

		HRESULT Initialize(IStorageFile^ storageFile, FileAccessMode fileAccessMode);

		IStorageFile^ GetTempStorageFile()
		{
			return m_tempStorageFile;
		}

		HRESULT Finalize();

	private:
		IStorageFile^ m_storageFile;
		IStorageFile^ m_tempStorageFile;
		FileAccessMode m_fileAccessMode;
	};
}

#endif