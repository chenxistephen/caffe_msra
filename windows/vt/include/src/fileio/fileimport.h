//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//     Image importer component.  Component is exposed as an IIndexedImageReader
//     It loads images on demand into the requested 'cache'.  It can optionally
//     prefetch images that will be used in the future.  
//
//  History:
//      2006/10/25-mattu
//            Created
//
//------------------------------------------------------------------------

#if !defined(VT_WINRT)

#include "vtcommon.h"

namespace vt {

//+----------------------------------------------------------------------------
//
// struct: IMPORT_IMAGE_FILE_INFO
// 
// Synposis: describes the images to be imported
// 
//-----------------------------------------------------------------------------
enum eFileImportType
{
	eFileImportStill = 0,
	eFileImportVideo = 1
};

struct IMPORT_IMAGE_FILE_INFO
{
	eFileImportType eType; // still image or video

	vt::wstring filename;  // filename can be a still image or a video

	// following 4 fields only apply to videos
	double     frameTime;     // frame time in seconds
	int        frameWidth;    // width and height of this particular frame
	int        frameHeight;
	int quarterRotationCount; // how to rotate the video

	IMPORT_IMAGE_FILE_INFO() :  
		eType(eFileImportStill), frameTime(0), frameWidth(0), frameHeight(0),
		quarterRotationCount(0)
	{} 
};

//+----------------------------------------------------------------------------
//
// Function: AllocateForImageImport
// 
// Synposis: given an CImgListInCache and a list of file names, this function
//           adds the dimension and pixel format info to the cache for each 
//           image
// 
//-----------------------------------------------------------------------------
HRESULT
AllocateForImageImport(CImgListInCache* pList,
					   const vt::vector<IMPORT_IMAGE_FILE_INFO>& vecSources,
					   int iOverridePixelType = OBJ_UNDEFINED,
					   bool bLoadMetadata     = true,
					   IMG_CACHE_SOURCE_PROPERTIES* props = NULL);

//+----------------------------------------------------------------------------
//
// Class: CImageImport
// 
// Synposis:
// 
//-----------------------------------------------------------------------------
class CImageImport: public IIndexedImageReaderWriter
{
public:
    CImageImport() : m_hPoolNotFullSubscribed(NULL), m_hr(S_OK)
	{}
 
    ~CImageImport()
	{ Deallocate(); }

	// implementation of IIndexedImageReader
public:
	virtual CLayerImgInfo GetImgInfo(UINT uIndex, UINT uLevel = 0)
	{ 
		return (Ok(uIndex)==S_OK)? m_pDst->GetImgInfo(uIndex, uLevel): 
			                       CLayerImgInfo();
	} 

	virtual HRESULT GetMetaData(UINT uIndex, OUT CParams& params);

	virtual HRESULT ReadRegion(UINT uIndex, IN const CRect& region, 
							   OUT CImg& dst, UINT uLevel = 0)
	{
		HRESULT hr = Ok(uIndex);
		if( hr == S_OK && (hr = LoadImage(uIndex)) == S_OK )
		{
			hr = m_pDst->ReadRegion(uIndex, region, dst, uLevel);
		}
		return hr;
	}

	virtual HRESULT Prefetch(UINT uIndex, IN const CRect&, 
							 UINT)
	{ 
		HRESULT hr = Ok(uIndex);
		return (hr != S_OK)? hr: LoadImage(uIndex, false);
    }

	virtual UINT GetFrameCount()
	{ return (UINT)m_vecFileInfo.size(); }

	// implementation of IIndexedImageWriter
public:
	virtual HRESULT  WriteRegion(UINT uIndex, IN const CRect& region,
								 IN const CImg& src, UINT uLevel=0)
	{
		HRESULT hr = Ok(uIndex);
		if( hr == S_OK && (hr = LoadImage(uIndex)) == S_OK )
		{
			hr = m_pDst->WriteRegion(uIndex, region, src, uLevel);
		}
		return hr;
	}

	virtual HRESULT  Fill(UINT uIndex, const void* pbValue, 
						  const RECT *prct = NULL, UINT uLevel = 0)
	{
		HRESULT hr = Ok(uIndex);
		if( hr == S_OK && (hr = LoadImage(uIndex)) == S_OK )
		{
			hr = m_pDst->Fill(uIndex, pbValue, prct, uLevel);
		}
		return hr;
	}

   	virtual HRESULT  SetMetaData(UINT uIndex, IN const CParams* pParams);

	// public methods specific to this class
public:
	void    Deallocate();

	HRESULT Initialize(const vt::vector<IMPORT_IMAGE_FILE_INFO>& vecSources,   
					   IIndexedImageReaderWriter* pDst, 
					   UINT uConcurrentReadCount);

	HRESULT LoadImage(UINT index)
	{ return LoadImage(index, true); }

	bool IsLoaded(UINT index)
	{ 
		return ( index >= m_vecFileInfo.size() )? 
			false: m_vecFileInfo[index].state == FILE_INFO::eLoaded;
    }

	HRESULT UnloadImage(UINT index);

	//+-------------------------------------------------------------------------
	//
	// Class: CTaskInfo
	// 
	// Synopsis: supplies the task progress handling and stores task specific
	//           information.
	// 
	//--------------------------------------------------------------------------
protected:
	class CTaskInfo: public CTaskStatusEvent
	{
	public:
		virtual HRESULT SetDone()
		{   
			HRESULT hr  = m_pImport->SetTaskDone(m_uIndex);
			HRESULT hr2 = CTaskStatusEvent::SetDone();
			return (hr == S_OK)? hr2: hr;
		}

		virtual void SetTaskError(HRESULT hr)
		{
			m_pImport->SetError(hr);
			CTaskStatusEvent::SetTaskError(hr);
		}

	public:
		CTaskInfo(): m_pImport(NULL), m_iRefCount(0), m_pVideoSrc(NULL)
		{}
		~CTaskInfo()
		{
			if( m_pVideoSrc )
			{
				m_pVideoSrc->Release();
			}
		}

		LONG RefCount()
		{ return m_iRefCount; }

		LONG AddRef()
		{ return ++m_iRefCount; }

		LONG Release()
		{ return --m_iRefCount; }

		HRESULT StartTask(CImageImport* pImport, UINT index);

	protected:  
		friend CImageImport;

		CImageImport* m_pImport;
		UINT          m_uIndex;
		LONG          m_iRefCount;

		// video specific members - have one video source per loading thread 
		IVideoSrc*    m_pVideoSrc;
	};

protected:
	struct FILE_INFO
	{  
		enum eState
		{
			eInit        = 0,
			eLoadStarted = 1,
			eLoaded      = 2
		};

		IMPORT_IMAGE_FILE_INFO srcinfo;
		volatile eState        state;          
		CTaskInfo*             pTaskInfo;
		FILE_INFO() : state(eInit), pTaskInfo(NULL)
		{}
	};

	// internal class methods
protected:
	HRESULT Ok(UINT uIndex)
	{ 
		return (m_hr != S_OK)? m_hr: 
			   (m_pDst==NULL)? E_NOINIT: 
			   (uIndex>=m_vecFileInfo.size())? E_INVALIDARG: S_OK;
	}

	HRESULT LoadImage(UINT index, bool bWait); 

	HRESULT WaitForStatusPool()
	{
		DWORD status = vt::CSystem::WaitForSingleObjectEx(m_hPoolNotFullSubscribed, INFINITE, FALSE);
		return (status==WAIT_FAILED)? 
			HRESULT_FROM_WIN32(GetLastError()): S_OK;
	}

	HRESULT WaitForLoad(CTaskInfo* pTaskInfo);

	static HRESULT LoadImageTask(void* pArg, LONG i, CTaskProgress* pProgress);

	HRESULT LoadImageTask(CImageImport::CTaskInfo* pInfo, 
						  CTaskProgress* pProgress);

	HRESULT LoadStillImageTask(const WCHAR* pFilename, UINT uIndex, 
							   CTaskProgress* pProgress);

	HRESULT LoadVideoFrameTask(const IMPORT_IMAGE_FILE_INFO& srcinfo, 
							   CImageImport::CTaskInfo* pInfo, 
							   CTaskProgress* pProgress);

	HRESULT SetTaskDone(UINT uIndex);

	void SetError(HRESULT hr)
	{ 
		if( m_hr == S_OK ) 
		{
			m_hr = hr;
		}
	}

	HRESULT ReleaseTaskInfo(CTaskInfo* pTaskInfo)
	{
		// NOTE: this must be called from within the critsection
		HRESULT hr = S_OK;
		if( pTaskInfo->Release() == 0 )
		{
			if( !SetEvent(m_hPoolNotFullSubscribed) )
			{
				hr = HRESULT_FROM_WIN32( GetLastError() );
			}
		}
		return hr;
	}


protected:
	volatile HRESULT           m_hr;
	vt::CCritSection           m_cs;
	HANDLE                     m_hPoolNotFullSubscribed;
	vt::vector<CTaskInfo>      m_vecTaskInfoPool;
	vt::vector<FILE_INFO>      m_vecFileInfo;
	IIndexedImageReaderWriter* m_pDst; 
	vt::CCritSection           m_filecs;
};

//+----------------------------------------------------------------------------
//
// Class: CImgInCacheImport
// 
// Synposis:
// 
//-----------------------------------------------------------------------------
class CImgInCacheImport: public CImageImport
{
public:
	~CImgInCacheImport()
	{ Deallocate(); }

	void Deallocate()
	{
		CImageImport::Deallocate();
		m_img.Deallocate();
	}

	HRESULT Initialize(const vt::vector<IMPORT_IMAGE_FILE_INFO>& vecSources,   
					   UINT uConcurrentReadCount,
					   int iOverridePixelType = OBJ_UNDEFINED,
					   bool bLoadMetadata     = true,
					   IMG_CACHE_SOURCE_PROPERTIES* props = NULL)
	{
		Deallocate();
		HRESULT hr = AllocateForImageImport(
			&m_img, vecSources, iOverridePixelType, bLoadMetadata, props);
		if( hr == S_OK )
		{
			hr = CImageImport::Initialize(vecSources, &m_img, uConcurrentReadCount);
		}
		return hr;
	}

	HRESULT UnloadImage(UINT index)
	{
		HRESULT hr = E_NOINIT;
		CImageCache* pCache = VtGetImageCacheInstance();
		if( pCache )
		{
			hr = pCache->ClearSource(m_img.GetCacheSourceId(index));
		}
		HRESULT hr2 = CImageImport::UnloadImage(index);
		return (hr==S_OK)? hr2: hr;
	}

protected:
	CImgListInCache m_img;
};

};

#endif
