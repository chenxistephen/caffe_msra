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

#include "stdafx.h"

#if defined(VT_WINRT)

#include "TempFileHelper.h"

using namespace vt;
using namespace Platform;

HRESULT CTempFileHelper::Initialize(IStorageFile^ storageFile, FileAccessMode fileAccessMode)
{
	if (storageFile == nullptr)
	{
		return E_POINTER;
	}

	VT_HR_BEGIN();

	try
	{
		m_storageFile = storageFile;
		m_fileAccessMode = fileAccessMode;

		// Create the VisionTools temporary folder, if it doesn't exist yet.
		IStorageFolder^ tempFolder = CSystem::Await(ApplicationData::Current->TemporaryFolder->CreateFolderAsync(ref new String(L"VisionTools"),
			CreationCollisionOption::OpenIfExists));

		// Create a temporary file with the same extension as the original file.
		const WCHAR* format = L"%08x-%04x-%04x-%02x%02x-%02x%02x%02x%02x%02x%02x%s";
		GUID guid = { 0 };
		VT_HR_EXIT(CoCreateGuid(&guid));
		const WCHAR* inputFilename = m_storageFile->Path->Data();
		const WCHAR* extension = VtGetFileExt(inputFilename);
		WCHAR tempFilename[MAX_PATH];
		if (_snwprintf_s(tempFilename, MAX_PATH, MAX_PATH, format, guid.Data1, guid.Data2, guid.Data3, UINT(guid.Data4[0]), UINT(guid.Data4[1]),
			UINT(guid.Data4[2]), UINT(guid.Data4[3]), UINT(guid.Data4[4]), UINT(guid.Data4[5]), UINT(guid.Data4[6]), UINT(guid.Data4[7]), extension) < 0)
		{
			VT_HR_EXIT(E_FAIL);
		}
		String^ tempFilenameStr = ref new String(tempFilename);
		m_tempStorageFile = CSystem::Await(tempFolder->CreateFileAsync(tempFilenameStr, CreationCollisionOption::FailIfExists));

		// Copy the original file to the temporary file.
		CSystem::Await(m_storageFile->CopyAndReplaceAsync(m_tempStorageFile));
	}
	catch (...)
	{
		VT_HR_EXIT(E_FAIL);
	}

	VT_HR_EXIT_LABEL();

	if (FAILED(hr))
	{
		m_storageFile = nullptr;
		m_tempStorageFile = nullptr;
	}
	return hr;
}

HRESULT CTempFileHelper::Finalize()
{
	VT_HR_BEGIN();

	try
	{
		if (m_tempStorageFile != nullptr)
		{
			// If we're writing, copy the temporary file contents back to the original file.
			if (m_fileAccessMode == FileAccessMode::ReadWrite && m_storageFile != nullptr)
			{
				CSystem::Await(m_tempStorageFile->CopyAndReplaceAsync(m_storageFile));
			}

			// Delete the temporary file.
			CSystem::Await(m_tempStorageFile->DeleteAsync(StorageDeleteOption::PermanentDelete));
		}
	}
	catch (...)
	{
		VT_HR_EXIT(E_FAIL);
	}

	VT_HR_EXIT_LABEL();

	// Clear.
	m_storageFile = nullptr;
	m_tempStorageFile = nullptr;
	return hr;
}

#endif
