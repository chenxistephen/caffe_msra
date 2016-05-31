//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Base class for singleton objects.
// 
//  History:
//      2007/06/1-mattu
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"

namespace vt {

template <class T>
class CVTSingleton: public T
{
public:
	CVTSingleton() : m_uRefCount(1) 
	{}
    ~CVTSingleton() 
    {}

	static HRESULT Startup();

	static void    Shutdown()
    {
        if( InterlockedDecrement(&m_singletonInstance->m_uRefCount) == 0 )
	    {
		    delete m_singletonInstance;
		    m_singletonInstance = NULL;
	    }
    }

    static T* GetInstance()
    { return m_singletonInstance; }

protected:
	LONG                    m_uRefCount;
    static CVTSingleton<T>* m_singletonInstance;
};

template <class T>
HRESULT CVTSingleton<T>::Startup()
{
	HRESULT hr = S_OK;

	//
	// This is a double-checked lock implementation of Singleton that 
	// is adapted from section 6 of: 
	// http://www.aristeia.com/Papers/DDJ_Jul_Aug_2004_revised.pdf
	// 
	// The CreateMutex solution for the lock is described in the 
	// "Windows: Named-Mutex Solution" of this article:
	// http://drdobbs.com/article/print?articleId=199203083&siteSectionName=    
	//

	CVTSingleton<T>* pS = m_singletonInstance;
	MemoryBarrier();

	if( pS == NULL )
	{
		//
		// Create a named mutex.  there can be only on instance of this in 
		// the system. but since this only locks once when a VT app is starting
		// up it shouldn't cause perf. issues due to blocking on it.
		//
		// All named mutexs are system wide - no way to create one of these
		// just in the scope of the current process.  The Local\ name is a 
		// small optimization that limits the scope of the mutex to just the 
		// current user's session.
		//
        HANDLE hMutex = CSystem::CreateMutexW(L"Local\\Microsoft:VisionTools:GlobalsMutex");
		if( hMutex == NULL )
		{
			return HRESULT_FROM_WIN32( GetLastError() );
		}

		// wait for access.  if the wait is successful we have an exclusive
		// lock on the one 'VisionTools:Globals' mutex
		DWORD dwStatus = vt::CSystem::WaitForSingleObjectEx(hMutex, INFINITE, FALSE);
		if( dwStatus == WAIT_OBJECT_0 )
		{
			// check for that this is still NULL, if it isn't another thread
			// created the singleton first
			if( m_singletonInstance == NULL )
			{
				pS = VT_NOTHROWNEW CVTSingleton<T>();
				if( pS == NULL )
				{
					hr = E_OUTOFMEMORY;
				}
				else
				{ 
					hr = pS->StartupInternal();
					MemoryBarrier();

					if( hr == S_OK )
                    {
						m_singletonInstance = pS;
                    }
                    else
					{
						delete pS;   
					}
				}
			}
            else
            {
			    // increment the refcount
			    m_singletonInstance->m_uRefCount++;
            }

			// we are done with the mutex
			ReleaseMutex(hMutex);
		}
		else
		{
			hr = (dwStatus==WAIT_FAILED)? 
				HRESULT_FROM_WIN32( GetLastError() ): E_UNEXPECTED;
		}
		CloseHandle(hMutex);
	}
    else
    {
	    InterlockedIncrement(&pS->m_uRefCount);
    }

	return hr;
}

};
