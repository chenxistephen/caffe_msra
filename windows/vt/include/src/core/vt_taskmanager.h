//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Functions used to control the execution of 'tasks'. The most
//      important functionality is to farm tasks out to worker threads
// 
//  History:
//      2007/06/1-mattu
//          Created
//
//------------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_atltypes.h"

namespace vt {

class CTaskProgress; // forward declaration
class CTaskStatus;   // forward declaration
class ITaskState;    // forward declaration

//+-----------------------------------------------------------------------------
//
// Class: ITaskWorkIdSequencer
//
/// \ingroup transforms
/// <summary>The base interface for generating the iWorkId's passed into the
/// TASK_CALLBACK.
/// </summary>
//
//------------------------------------------------------------------------------
class ITaskWorkIdSequencer
{
public:
	virtual ~ITaskWorkIdSequencer() = 0
	{}

	virtual int     GetTotalWorkItems() = 0;

	virtual HRESULT Advance(OUT bool& bDone, OUT int& iWorkId, int iThreadIndex) = 0;
};

//+-----------------------------------------------------------------------------
//
// Class: CSimpleCounterWorkIdSequencer
//
/// \ingroup transforms
/// <summary>A simple implementation of ITaskWorkIdSequencer that just counts
/// up from 0 to total work items
/// </summary>
//
//------------------------------------------------------------------------------
class CSimpleCounterWorkIdSequencer: public ITaskWorkIdSequencer
{
public:
	virtual int GetTotalWorkItems()
	{ return m_iTotalWorkItems; }

	virtual HRESULT Advance(OUT bool& bDone, OUT int& iWorkId, int)
	{
		// -- below because we are 0-based
		iWorkId = InterlockedIncrement(&m_iCurrentWorkId);
		bDone   = (iWorkId > m_iTotalWorkItems); 
		iWorkId--;
		return S_OK;
	}

public:
	CSimpleCounterWorkIdSequencer() :
		m_iTotalWorkItems(0), m_iCurrentWorkId(0)
	{}

	CSimpleCounterWorkIdSequencer(int iTotalWorkItems) : 
		m_iTotalWorkItems(iTotalWorkItems), m_iCurrentWorkId(0)
	{}

	void Initialize(int iTotalWorkItems)
	{
		m_iTotalWorkItems = iTotalWorkItems;
		m_iCurrentWorkId  = 0;
	}

protected:
	LONG m_iTotalWorkItems;
	volatile LONG m_iCurrentWorkId;
};

//+-----------------------------------------------------------------------------
//
// Callback: TASK_CALLBACK, TASK_WITH_STATE_CALLBACK
// 
// Synposis: User provided callback that the task manager calls for each unit
//           of work.
// 
//           iWorkId: by default a number between 0 and count-1 that can be used
//                  to indicate which unit of work is to be performed, where
//                  'count' is the param provided in PushTask.  If an  
//                  ITaskWorkIdSequencer is passed to PushTask then this 
//                  component will generate the Ids.
//
//------------------------------------------------------------------------------
typedef HRESULT (*TASK_CALLBACK)(void* pArg, LONG iWorkId, CTaskProgress*);
typedef HRESULT (*TASK_WITH_STATE_CALLBACK)(ITaskState* pArg, LONG iWorkId, 
										    CTaskProgress*);

//+-----------------------------------------------------------------------------
//
// Function: PushTask
// 
// Synposis: functions used to submit work to the taskmanager
//
//------------------------------------------------------------------------------
struct VT_TASK_OPTIONS
{
	DWORD maxthreads;   // zero means don't throttle
	VT_TASK_OPTIONS() : maxthreads(0)
	{}
	VT_TASK_OPTIONS(DWORD t) : maxthreads(t)
	{} 
};

HRESULT	PushTask(TASK_CALLBACK callback, void* pArg, CTaskStatus* pStatus, 
				 LONG count = 1, const VT_TASK_OPTIONS* pOpts = NULL);

HRESULT	PushTaskAndWait(TASK_CALLBACK callback, void* pArg,	
						CTaskStatus* pStatus, LONG count = 1,
						const VT_TASK_OPTIONS* pOpts = NULL);

HRESULT	PushTask(TASK_WITH_STATE_CALLBACK callback, ITaskState* pArg, 
				 CTaskStatus* pStatus, LONG count = 1,
				 const VT_TASK_OPTIONS* pOpts = NULL);

HRESULT	PushTaskAndWait(TASK_WITH_STATE_CALLBACK callback, ITaskState* pArg,	
						CTaskStatus* pStatus, LONG count = 1,
						const VT_TASK_OPTIONS* pOpts = NULL);

HRESULT	PushTask(TASK_CALLBACK callback, void* pArg, CTaskStatus* pStatus, 
				 ITaskWorkIdSequencer* pSeq, const VT_TASK_OPTIONS* pOpts = NULL);

HRESULT	PushTaskAndWait(TASK_CALLBACK callback, void* pArg,	
						CTaskStatus* pStatus, ITaskWorkIdSequencer* pSeq,
						const VT_TASK_OPTIONS* pOpts = NULL);

HRESULT	PushTask(TASK_WITH_STATE_CALLBACK callback, ITaskState* pArg, 
				 CTaskStatus* pStatus, ITaskWorkIdSequencer* pSeq,
				 const VT_TASK_OPTIONS* pOpts = NULL);

HRESULT	PushTaskAndWait(TASK_WITH_STATE_CALLBACK callback, ITaskState* pArg,	
						CTaskStatus* pStatus, ITaskWorkIdSequencer* pSeq,
						const VT_TASK_OPTIONS* pOpts = NULL);

//+-----------------------------------------------------------------------------
//
// Function: 
// 
// Synposis: Startup and Shutdown functions for task manager and image cache
//           and file blob store
//
//------------------------------------------------------------------------------
HRESULT TaskManagerStartup();
void    TaskManagerShutdown();

HRESULT ImageCacheStartup();
void    ImageCacheShutdown();

//+-----------------------------------------------------------------------------
//
// Interface: ITaskState
// 
// Synopsis: base interface for concurrent task the need to be forked to
//           multiple threads (clone / thread) and joined merge is called
//           to collect the state after all threads are done
// 
//------------------------------------------------------------------------------
class ITaskState
{
public:
	virtual ~ITaskState() = 0
	{} 

	virtual bool RequiresCloneForConcurrency()
	{ return false;	}

	virtual HRESULT Clone(ITaskState** ppClone) = 0;

	virtual HRESULT Merge(const ITaskState*)
	{ return S_OK; }
};

// helper methods to implement ITaskState::Clone
template<class T>
HRESULT CloneTaskState(ITaskState **ppClone, T* pNew)
{ 
	HRESULT hr = E_POINTER;
	if( ppClone )
	{
		hr = pNew? S_OK: E_OUTOFMEMORY;
		*ppClone = pNew;
	}
	else
	{
		delete pNew;
	}
	return hr;
}

template <class T, class CloneFunc>
HRESULT CloneTaskState(ITaskState **ppState, 
					   const CloneFunc& f)
{
	if( ppState == NULL )
	{
		return E_POINTER;
	}

	*ppState = NULL;
	T* pNew  = VT_NOTHROWNEW T();
	if( pNew == NULL )
	{
		return E_OUTOFMEMORY;
	}

	HRESULT hr = f(pNew);
	if( hr != S_OK )
	{
		delete pNew;
		return hr;
	}
	*ppState = pNew;
	return S_OK;
}

//+-----------------------------------------------------------------------------
//
// Class: CTaskProgress
// 
// Synopsis: interface for a task to report its progress.
// 
//------------------------------------------------------------------------------
class CTaskProgress
{
public:
	virtual bool  GetCancel() = 0;
	virtual void  ReportProgress(float fPct) = 0;
};

//+-----------------------------------------------------------------------------
//
// Class: CPhasedTaskProgress
// 
// Synopsis: an extension of CTaskProgress for complex tasks that need to 
//           step through multiple phases
// 
//------------------------------------------------------------------------------
template<class Base>
class CPhasedTaskProgressImpl: public Base
{
public:
	CPhasedTaskProgressImpl() : m_pTaskProgOuter(NULL), 
        m_fBase(0), m_fProgressInPhase(0), m_fPhasePct(0)
	{}

    // Sets the callback function pointer and cookie
	void SetOuterCallback(CTaskProgress* pCallback)
	{
		m_pTaskProgOuter = pCallback;
	}

    // Call Start before any progress callbacks occur. Typically fBase is
    // 0, but it can be more if some previous operation has been under way and
    // progress should continue from there.
	void Start(float fBase = 0)
	{
		m_fBase = fBase;
		m_fPhasePct = 0;
		m_fProgressInPhase = 0;
	}

    // Start a new phase. The fPct parameter is the contribution this phase
    // will make to the overall 100%. The sum of the fPct parameters for 
	// all BeginPhase calls should add up to 100 - fBase.
	void BeginPhase(const char* phaseName, float fPct);

	virtual float GetProgress() const
	{ return m_fBase + m_fPhasePct*m_fProgressInPhase; }

	virtual float GetProgressInPhase() const
	{ return m_fProgressInPhase; }

    // callbacks to the outer task's CTaskProgress - accumulate the current
	// phase progress to report an overall progress
	virtual bool GetCancel()
	{
		return m_pTaskProgOuter? m_pTaskProgOuter->GetCancel(): false;
	}

	virtual void ReportProgress(float fPct);

protected:
	CTaskProgress* m_pTaskProgOuter;

    // Starting % for this phase
	float  m_fBase;

    // Full % contribution that this phase makes to the overall 100%. When
    // progress callbacks come in (which represent the percentage completion
    // 0 - 100 of the *current phase*), those values are scaled so that when
    // they reach 100 m_fPhasePct has been added to the overall 100%.
    float  m_fPhasePct;
    float  m_fProgressInPhase;
};

typedef CPhasedTaskProgressImpl<CTaskProgress> CPhasedTaskProgress;

//+-----------------------------------------------------------------------------
//
// function: CheckTaskCancel / CheckTaskProgressCancel
// 
// Synopsis: Simple implementation for a task to use to check for cancel request
//           Typical use is within a HR checking macro e.g. 
//           VT_HR_EXIT( CheckTaskCancel(pProg) );
// 
//------------------------------------------------------------------------------
inline HRESULT
CheckTaskCancel(CTaskProgress* pProg)
{
	return (pProg && pProg->GetCancel())? E_ABORT: S_OK;
}

inline HRESULT
CheckTaskProgressCancel(CTaskProgress* pProg, float fProgress)
{
	HRESULT hr = S_OK;
	if( pProg && (hr = (pProg->GetCancel()? E_ABORT: S_OK)) == S_OK )
	{
		pProg->ReportProgress(fProgress);
	}
	return hr;
}

inline void
SetProgress100Percent(CTaskProgress* pProg)
{
	if( pProg )
	{
		pProg->ReportProgress(100.f);
	}
}

//+-----------------------------------------------------------------------------
//
// Class: CTaskStatus
// 
// Synopsis: interface for the task manager to report status (progress, 
//           cancel, done, error codes) of a task.  Following are some simple
//           implementations used by code that issues tasks.
// 
//           CPhasedTaskStatus - for users who only want progress/cancel and
//                               don't need to query for done. 
//                               
//           CTaskStatusPoll - for users who simply want to poll for the 
//                             task completion/progress.
// 
//           CTaskStatusEvent - for users who want an event on task completion.
//------------------------------------------------------------------------------
class CTaskStatus: public CTaskProgress
{
public:
	CTaskStatus() : m_hrTask(S_OK)
	{}

	virtual HRESULT SetDone() = 0;

	virtual void SetTaskError(HRESULT hr)
	{ if(m_hrTask==S_OK) m_hrTask = hr; }

	virtual HRESULT GetTaskError()
	{ return m_hrTask; }

protected:
	HRESULT m_hrTask;
};

class CPhasedTaskStatus: public CPhasedTaskProgressImpl<CTaskStatus>
{
	virtual HRESULT SetDone()
	{ return S_OK; }
};

class CTaskStatusPoll: public CTaskStatus
{
public:
	CTaskStatusPoll() { Begin(); }

	// call before issuing a task
	virtual HRESULT Begin()
	{
		m_bDone   = false;
		m_bCancel = false;
		m_fProgressPct = 0;
		m_hrTask    = S_OK;
		return S_OK;
	}

	virtual void SetCancel()
	{ m_bCancel = true; }

	virtual float GetProgress()
	{ return m_fProgressPct; }

	virtual bool GetDone()
	{ return m_bDone; }

	// implementation of CTaskStatus
	virtual bool GetCancel()
	{ return m_bCancel; }

	virtual void ReportProgress(float fPct);

	virtual HRESULT WaitForDone();

	virtual HRESULT SetDone()
	{ m_bDone = true; return S_OK; }

protected:
	volatile bool m_bCancel;
	float   m_fProgressPct;
	bool    m_bDone;
};

class CTaskStatusEvent: public CTaskStatusPoll
{
public:
	CTaskStatusEvent() : m_hTaskIdle(NULL)
	{ }

	~CTaskStatusEvent()
	{
		if( m_hTaskIdle )
		{
			CloseHandle(m_hTaskIdle);
		}
	}

	// call to create the event 
	virtual HRESULT Create();

	// call before issuing a task
	virtual HRESULT Begin();

	virtual HRESULT WaitForDone();

	// callback to indicate done
	virtual HRESULT SetDone();

protected:
	HANDLE m_hTaskIdle;
};

//+-----------------------------------------------------------------------------
//
// Class: CPixelCountProgress
// 
// Synopsis: an specialization of CPhasedTaskStatus where the phase are regions
//           within a larger area.  The total size of the larger area is 
//           defined in Start() and each phase is described by its dimensions
//           in pixels. 
// 
//------------------------------------------------------------------------------
class CPixelCountProgress: public CPhasedTaskStatus
{
public:
	CPixelCountProgress() : m_fPctPerDstPixel(100.f)
	{}
	
	template <class T>
    void Start(IN T* pT, IN CTaskProgress* pOuterProgress)
	{
		SetOuterCallback(pOuterProgress);
		CPhasedTaskStatus::Start();
		m_fPctPerDstPixel = 100. / double(VtMax<UINT64>( TotalPixelCount(pT), 1 ));
	}

	// setup progress for a region of pixels
	void BeginRegion(UINT64 cnt, const char* regionName=NULL)
	{ BeginPhase(regionName, float(m_fPctPerDstPixel * double(cnt))); }

	void BeginRegion(int w, int h, const char* regionName=NULL)
	{ BeginRegion(UINT64(w)*UINT64(h), regionName); }

	void BeginRegion(const vt::CRect& r, const char* regionName=NULL)
	{ BeginRegion(r.Width(), r.Height(), regionName); }

protected:
	double m_fPctPerDstPixel;
};

//+-----------------------------------------------------------------------------
//
// Class: CCritSection
// 
// Synopsis: helper function to deal with common critical section pattern
// 
//------------------------------------------------------------------------------
class CCritSection
{
public:
	CCritSection() : m_bInit(false)
	{}

	~CCritSection()
	{ Shutdown(); }  

	HRESULT Startup(int spincount=4000);

	void Shutdown();

	void Enter();

	void Leave();

protected:
	bool             m_bInit;
	CRITICAL_SECTION m_cs;
};


//+-----------------------------------------------------------------------------
//
// Class: CVTResourceManagerInit
// 
// Synposis: An RAII class that handles startup and shutdown of the VT 
//           task manager and memory cache.  Typical use is to create a 
//           global instance of this, or to declare one on the stack.  This
//           calls shutdown when it goes out of scope.  Startup/Shutdown are
//           reference counted, so an application can safely use 
//           CVTResourceManagerInit more than once.
// 
//------------------------------------------------------------------------------
class CVTResourceManagerInit
{
public:
	CVTResourceManagerInit()
	{
		m_bTMStarted = false; 
		if( TaskManagerStartup() == S_OK )
		{
			m_bTMStarted = true;
		}

        m_bICStarted = false; 
		if( ImageCacheStartup() == S_OK )
		{
			m_bICStarted = true;
		}
	}
 
	~CVTResourceManagerInit()
	{
		if( m_bTMStarted )
		{
			TaskManagerShutdown();
		}

        if( m_bICStarted )
		{
			ImageCacheShutdown();
		}
	}

private:
	bool m_bTMStarted;
    bool m_bICStarted;
};

};
