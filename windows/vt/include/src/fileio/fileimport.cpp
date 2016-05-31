//+-----------------------------------------------------------------------------
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//     File importer component.  Component is exposed as an IIndexedImageReader
//     It loads images on demand into the requested 'cache'.  It can optionally
//     prefetch images that will be used in the future.  
//
//  History:
//      2006/10/25-mattu
//            Created
//
//------------------------------------------------------------------------------

#include "stdafx.h"

#if !defined(VT_WINRT)

#include <shlwapi.h>

#include "wicio.h"
#include "fileimport.h"

using namespace vt;

HRESULT
CImageImport::GetMetaData(UINT uIndex, OUT CParams& params)
{
	HRESULT hr = Ok(uIndex);
	if( hr != S_OK )
	{
		return hr;
	}

	// If image has been loaded then get metadata out of the dest otherwise
	// load the dest metadata now.

	// NOTE: currently in the use of this in ICE the metadata is loaded
	//       into the m_pDst during the AllocateForImageImport step so it
	//       is OK to just defer to the m_pDst even if the pixels haven't
	//       been loaded. 
#if 0
	FILE_INFO& fi = m_vecFileInfo[uIndex];
	FILE_INFO::eState state = fi.state;
	MemoryBarrier();
	if( state != FILE_INFO::eLoaded )
	{ 
        // TODO: Need to fix metadata IO in the cache and other places.  
		//       It should be a shallow copy and return a refcounted meta-data
		//       class.  As things stand now, the metadata pointer returned
	    //       out of the cache isn't refcounted so it can go away before
		//       all references to it are done. 
		CWicReader rd;
		if( (hr = rd.OpenFile(fi.name)) == S_OK )
		{
			m_pDst->SetMetaData(uIndex, rd.GetMetaData());
		}
	}
#endif

	return m_pDst->GetMetaData(uIndex, params);
}

HRESULT
CImageImport::SetMetaData(UINT uIndex, IN const CParams* pParams)
{
	HRESULT hr = Ok(uIndex);
	if( hr != S_OK )
	{
		return hr;
	}

	// NOTE: assumption here is that the m_pDst is multi-thread safe (as the 
	//       image-cache is). It may be neccsary in the future to include a
	//       crit-section around this if this assumuption is not true
	return m_pDst->SetMetaData(uIndex, pParams);
}

void
CImageImport::Deallocate()
{
	// wait for anything outstanding
	if( m_hr == S_OK )
	{
		m_hr = E_ABORT;   
	}
	if( m_hPoolNotFullSubscribed )
	{
		SetEvent(m_hPoolNotFullSubscribed);
	}
	for( UINT i = 0; i < m_vecTaskInfoPool.size(); i++ )
	{
		// if there are outstanding refs to the tasks, wait for them to complete
		CTaskInfo* pInfo = m_vecTaskInfoPool.begin() + i;

		bool bWait = false;

		m_cs.Enter();
		if( pInfo->RefCount() > 0 )
		{
			pInfo->AddRef();
			bWait = true;
		}
		m_cs.Leave();

		if( bWait )
		{
			WaitForLoad(pInfo);
		}
		// TODO: this potentially still isn't quite right we'd really need to 
		//       wait until refs on the info are all 0.  but this is OK for now
		//       because all callers should shutdown client threads before calling
		//       deallocate
	}

	if( m_hPoolNotFullSubscribed )
	{
		CloseHandle(m_hPoolNotFullSubscribed);
		m_hPoolNotFullSubscribed = NULL;
	}
	m_vecFileInfo.clear();
	m_vecTaskInfoPool.clear();
	m_cs.Shutdown();
    m_filecs.Shutdown();
	m_pDst = NULL;
	m_hr = S_OK;  // clear the HR
}

HRESULT
CImageImport::Initialize(const vt::vector<IMPORT_IMAGE_FILE_INFO>& vecSources,   
						 IIndexedImageReaderWriter* pDst, 
						 UINT uConcurrentReadCount)
{
	HRESULT hr;

	VT_ASSERT( vecSources.size() == pDst->GetFrameCount() );

	// free stuff
	Deallocate();

	// allocate events
    VT_HR_EXIT( m_cs.Startup() );
    VT_HR_EXIT( m_filecs.Startup() );

	m_hPoolNotFullSubscribed = CSystem::CreateEvent(TRUE, TRUE);
	if( m_hPoolNotFullSubscribed == NULL )
	{
		VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );
	}

	VT_HR_EXIT( m_vecTaskInfoPool.resize(uConcurrentReadCount) );
	for( UINT i = 0; i < m_vecTaskInfoPool.size(); i++ )
	{
		VT_HR_EXIT( m_vecTaskInfoPool[i].Create() );
	}

	// allocate the file info structs
	VT_HR_EXIT( m_vecFileInfo.resize(vecSources.size()) );
	for( UINT i = 0; i < vecSources.size(); i++ )
	{
		IMPORT_IMAGE_FILE_INFO& dst = m_vecFileInfo[i].srcinfo;
		const IMPORT_IMAGE_FILE_INFO& src = vecSources[i];

		dst.eType = src.eType;
		dst.frameTime = src.frameTime;
		dst.quarterRotationCount = src.quarterRotationCount;
	    VT_HR_EXIT( dst.filename.assign(src.filename) );
	}

	m_pDst = pDst;

	// start the pipeline by pushing VtLoadImageTask tasks
	UINT uStartUpCnt = VtMin(uConcurrentReadCount, (UINT)m_vecFileInfo.size());
	for( UINT i = 0; i < uStartUpCnt; i++ )
	{
		// NOTE: normally fi.state needs to be updated in a thread-safe way
		//       but single threaded operation is assumed during initialization   
		FILE_INFO& fi = m_vecFileInfo[i];
		fi.state = FILE_INFO::eLoadStarted;
		fi.pTaskInfo = m_vecTaskInfoPool.begin() + i;
        fi.pTaskInfo->AddRef();
		VT_HR_EXIT( fi.pTaskInfo->StartTask(this, i) ); 
	}

Exit:
    if( hr != S_OK )
	{
		Deallocate();
		m_hr = hr;
	}
	return hr;
}

HRESULT
CImageImport::LoadImage(UINT index, bool bWait)
{
	FILE_INFO& fi = m_vecFileInfo[index];
	
	// first check outside of lock
	FILE_INFO::eState state = fi.state;
	MemoryBarrier();
	if( state == FILE_INFO::eLoaded || 
        (!bWait && state == FILE_INFO::eLoadStarted) )
	{
		return S_OK;
	}

	HRESULT hr = S_OK;

	bool bLoaded = false;
    do
    {
		bool bWaitForStatusPool = false;
		bool bStartTask         = false;
		bool bWaitForLoaded     = false;
	
		// if it hasn't been loaded lock the internal data structs
		// BEGIN LOCK  ====>
		m_cs.Enter();
			
		// if it has yet to be started
        if( fi.state == FILE_INFO::eInit )
		{
			// get an task out of the pool
			CTaskInfo* pTaskInfo = NULL;
			for( UINT i = 0; i < m_vecTaskInfoPool.size(); i++ )
			{
				CTaskInfo& el = m_vecTaskInfoPool[i];
				if( el.RefCount() == 0 )
				{
					pTaskInfo = &el;
					break;
				}
			}
	
			if( pTaskInfo == NULL )
			{
				if( bWait )
				{
					// if no event available need to wait on pool
					if( ResetEvent(m_hPoolNotFullSubscribed) )
					{
						bWaitForStatusPool = true;
					}
					else
					{
						hr = HRESULT_FROM_WIN32(GetLastError());
					}
				}
			}
			else
			{
				VT_ASSERT( fi.pTaskInfo == NULL );
				fi.state = FILE_INFO::eLoadStarted;
				fi.pTaskInfo = pTaskInfo;
				VT_ASSERT( fi.pTaskInfo->RefCount() == 0 );
				fi.pTaskInfo->AddRef();
				bStartTask = true;
			}
		}
        else if( fi.state == FILE_INFO::eLoadStarted )
		{
			// find the task info for this task so that we can wait on it
			// if necessary
			if( (bWaitForLoaded = bWait) == true )
			{
				fi.pTaskInfo->AddRef();
			}
		}
		else
		{
			// loaded came in since first check, so we are done
			bLoaded = true;
		}
	
		m_cs.Leave(); 
		// <==== END LOCK
	
		if( bWaitForStatusPool )
		{
			VT_ASSERT( bWait );
			hr = WaitForStatusPool();
		}

		if( bStartTask )
		{
			hr = fi.pTaskInfo->StartTask(this, index);
		}
 
		if( bWaitForLoaded )
		{
			VT_ASSERT( bWait );
			hr = WaitForLoad(fi.pTaskInfo);
			bLoaded = true;
		}

        SetError(hr);
	}
    while( bWait && !bLoaded && m_hr == S_OK );

	return m_hr;
}

HRESULT
CImageImport::UnloadImage(UINT index)
{
	HRESULT hr = S_OK;
	FILE_INFO& fi = m_vecFileInfo[index];
	bool bWait;
    do
    {
		// BEGIN LOCK  ====>
		m_cs.Enter();

		if( fi.state == FILE_INFO::eLoadStarted )
		{
			fi.pTaskInfo->AddRef();
			bWait = true;
		}
		else
		{
			fi.state = FILE_INFO::eInit;
			bWait = false;
		}

		m_cs.Leave(); 
		// <==== END LOCK

		if( bWait )
		{
			hr = WaitForLoad(fi.pTaskInfo);
		}
	} while( bWait && hr == S_OK );

	return hr;
}

HRESULT
CImageImport::SetTaskDone(UINT index)
{
	HRESULT hr = S_OK;

	VT_ASSERT( index < m_vecFileInfo.size() );

	FILE_INFO& fi = m_vecFileInfo[index];

	// switch state to eLoaded - it should currently be eLoadStarted
	VT_ASSERT( fi.state == FILE_INFO::eLoadStarted );

	// BEGIN LOCK  ====>
	m_cs.Enter();

	fi.state = FILE_INFO::eLoaded;

    // return info to the info pool
	hr = ReleaseTaskInfo(fi.pTaskInfo);

    fi.pTaskInfo = NULL;

	m_cs.Leave(); 
	// <==== END LOCK

	SetError(hr);
	return hr;
}

HRESULT
CImageImport::WaitForLoad(CTaskInfo* pTaskInfo)
{
	HRESULT hr = pTaskInfo->WaitForDone();

	// BEGIN LOCK  ====>
	m_cs.Enter();

	HRESULT hr2 = ReleaseTaskInfo(pTaskInfo);

	m_cs.Leave(); 
	// <==== END LOCK

	return (hr==S_OK)? hr2: hr; 
}

HRESULT
CImageImport::LoadImageTask(void* pArg, LONG, CTaskProgress* pProgress)
{
	CImageImport::CTaskInfo* pInfo = (CImageImport::CTaskInfo*)pArg;
	return pInfo->m_pImport->LoadImageTask(pInfo, pProgress);
}

HRESULT
CImageImport::LoadImageTask(CImageImport::CTaskInfo* pInfo,
							CTaskProgress* pProgress)
{
	const FILE_INFO& fi = m_vecFileInfo[pInfo->m_uIndex];
	if( fi.srcinfo.eType == eFileImportStill )
	{
		return LoadStillImageTask(fi.srcinfo.filename, pInfo->m_uIndex, pProgress);
	}
	else
	{
		return LoadVideoFrameTask(fi.srcinfo, pInfo, pProgress);
	}
}

HRESULT
CImageImport::LoadStillImageTask(const WCHAR* pFilename, UINT uIndex, 
								 CTaskProgress*)
{
	CTimer timer;
	// VT_DEBUG_LOG("Begin load %S\n", pFilename);

    HANDLE hFile = INVALID_HANDLE_VALUE;
    void* pMemory = NULL;

    VT_HR_BEGIN()

    if ((hFile = CreateFile(pFilename,
                            GENERIC_READ, FILE_SHARE_READ, NULL,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL |
                            FILE_FLAG_NO_BUFFERING, NULL)) ==
                            INVALID_HANDLE_VALUE)
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

    LARGE_INTEGER liFileSize;
    if (!GetFileSizeEx(hFile, &liFileSize))
        VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

	CIndexedImageWriterIterator wrcur(m_pDst, uIndex);

    if (liFileSize.HighPart == 0 &&
#ifdef _M_AMD64
        liFileSize.LowPart < 256*1024*1024)
#else
        liFileSize.LowPart < 128*1024*1024)
#endif
    {
        CWicReader rd;
        CMemStream ms;

        // if file less than 128 MB read into stream

        // find disk sector size
		wchar_t root[MAX_PATH];
		wcscpy_s(root, MAX_PATH, pFilename); 

		LPCWSTR lpRootPathName = PathStripToRoot(root)? root: NULL;
        DWORD dwSectorsPerCluster, dwBytesPerSector,
              dwNumberOfFreeClusters, dwTotalNumberOfClusters;
        if (!GetDiskFreeSpace(lpRootPathName, &dwSectorsPerCluster,
                              &dwBytesPerSector,
                              &dwNumberOfFreeClusters, &dwTotalNumberOfClusters))
		    VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

        // align to disk sector size for FILE_FLAG_NO_BUFFERING
        DWORD dwSectorMask = dwBytesPerSector - 1;
        DWORD dwSize = VtMin(16*1024*1024UL,
							 (liFileSize.LowPart + dwSectorMask) & ~dwSectorMask);
        pMemory = VirtualAlloc(NULL, dwSize, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
        VT_PTR_OOM_EXIT( pMemory );

        m_filecs.Enter();

        // read file 16 MB at a time into stream
        while (liFileSize.LowPart > 0 && hr == S_OK)
        {
            DWORD dwRead = VtMin(dwSize,
                (liFileSize.LowPart + dwSectorMask) & ~dwSectorMask);
            if (!ReadFile(hFile, pMemory, dwRead, &dwRead, NULL))
                hr = HRESULT_FROM_WIN32( GetLastError() );
            else if (dwRead != dwSize && dwRead != liFileSize.LowPart)
                hr = E_READFAILED;
            if (hr == S_OK)
                hr = ms.Write(pMemory, dwRead, NULL);
            liFileSize.LowPart -= dwRead;
        }

        m_filecs.Leave();

        VT_HR_EXIT( hr );

		LARGE_INTEGER liSeekPosition;
		liSeekPosition.LowPart = liSeekPosition.HighPart = 0;
		VT_HR_EXIT( ms.Seek(liSeekPosition, STREAM_SEEK_SET, NULL) );

	    VT_HR_EXIT( rd.OpenStream(&ms) );

	    VT_HR_EXIT( rd.GetImage(&wrcur) );
    }
    else
    {
        CWicReader rd;
        // otherwise read with Wic reader
        m_filecs.Enter();

	    if ((hr = rd.OpenFile(pFilename)) == S_OK)
		{
            hr = rd.GetImage(&wrcur);
		}

        m_filecs.Leave();

        VT_HR_EXIT( hr );
    }

    VT_HR_EXIT_LABEL()

    if (pMemory != NULL)
        VirtualFree(pMemory, 0, MEM_RELEASE);

    if (hFile != INVALID_HANDLE_VALUE)
        CloseHandle(hFile);
		
	// VT_DEBUG_LOG("End load %S - took %.2fsec\n", pFilename, timer.GetTimeMilliSec()/1000.);

	return hr;
}

HRESULT
CImageImport::LoadVideoFrameTask(const IMPORT_IMAGE_FILE_INFO& srcinfo,
								 CImageImport::CTaskInfo* pInfo, 
								 CTaskProgress*)
{
	VT_HR_BEGIN()

	VT_ASSERT( srcinfo.eType == eFileImportVideo );

	vt::CRGBImg img;

	if( pInfo->m_pVideoSrc==NULL || 
		_wcsicmp(pInfo->m_pVideoSrc->GetFileName(), srcinfo.filename)!=0 )
	{
		if( pInfo->m_pVideoSrc )
		{
		    pInfo->m_pVideoSrc->Release();
			pInfo->m_pVideoSrc = NULL;
		}
		
		VT_HR_EXIT( VtCreateVideoSrc(&(pInfo->m_pVideoSrc), srcinfo.filename) );
	}

	VT_HR_EXIT( pInfo->m_pVideoSrc->GetFrame(srcinfo.frameTime, img) );

	if ( (srcinfo.quarterRotationCount&3) == 0)
	{
		VT_HR_EXIT( m_pDst->WriteImg(pInfo->m_uIndex, img) );
	}
	else
	{
		vt::CRGBImg imgRot;
		VT_HR_EXIT( VtRotateImage(imgRot, img, 
								   srcinfo.quarterRotationCount) );
		VT_HR_EXIT( m_pDst->WriteImg(pInfo->m_uIndex, imgRot) );
	}

	VT_HR_END()
}

//+-----------------------------------------------------------------------------
//
// Class: CImageImport::CTaskInfo
// 
//-----------------------------------------------------------------------------
HRESULT
CImageImport::CTaskInfo::StartTask(CImageImport* pImport, UINT index)
{	
	HRESULT hr;
	
	VT_HR_EXIT( Begin() );
		
	// the m_iRefCount should have been bumped by the caller of StartTask
	// *** in  a thread-safe way ***
	VT_ASSERT( m_iRefCount != 0 );
	m_pImport = pImport;
	m_uIndex  = index;

	VT_HR_EXIT( PushTask(CImageImport::LoadImageTask, this, this) );

Exit:
	return hr;
} 

//+----------------------------------------------------------------------------
//
// Function: AllocateForImageImport
// 
//-----------------------------------------------------------------------------
HRESULT
vt::AllocateForImageImport(CImgListInCache* pList,
						   const vt::vector<IMPORT_IMAGE_FILE_INFO>& vecSources,
						   int iOverridePixelType,
						   bool bLoadMetadata,
						   IMG_CACHE_SOURCE_PROPERTIES* props)
{
	VT_HR_BEGIN()

	pList->Deallocate();

	for( UINT i = 0; i < vecSources.size(); i++ )
	{
		const IMPORT_IMAGE_FILE_INFO& ifi = vecSources[i];
		if( ifi.eType == eFileImportStill )
		{
			// TODO: the WIC reader should have an option to not load meta-data
			CWicReader rd;
			VT_HR_EXIT( rd.OpenFile(vecSources[i].filename) );
			CImgInfo info = rd.GetImgInfo();
			if( iOverridePixelType != OBJ_UNDEFINED )
			{
				info = CImgInfo(info.width, info.height, iOverridePixelType);
			}
			VT_HR_EXIT( pList->PushBack(info, props) );
			if( bLoadMetadata )
			{
				VT_HR_EXIT( pList->SetMetaData(i, rd.GetMetaData()) );
			}
			VT_HR_EXIT( rd.CloseFile() );
		}
		else
		{
			VT_ASSERT( ifi.eType == eFileImportVideo );

			int iPixType = (iOverridePixelType != OBJ_UNDEFINED)? 
				iOverridePixelType: OBJ_RGBIMG;
			CImgInfo info = ( ifi.quarterRotationCount & 1 )?
				CImgInfo(ifi.frameHeight,ifi.frameWidth,iPixType):
				CImgInfo(ifi.frameWidth,ifi.frameHeight,iPixType);

			VT_HR_EXIT( pList->PushBack(info, props) );
		}
	}

    VT_HR_EXIT_LABEL()

    if( hr != S_OK )
	{
		pList->Deallocate();
	}
	return hr;
}

#endif
