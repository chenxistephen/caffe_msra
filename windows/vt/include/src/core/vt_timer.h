//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Accurate timer class
//
//  History:
//      2004/11/8-mattu
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"

namespace vt {

// Timer class returns time as a 64 bit integer where each tick represents
// 100ns since the timer was created. Note that at this resolution, the timer
// can run for over 58,454 years without wrapping.

class CTimer
{
public:
    CTimer()
    {
        QueryPerformanceFrequency ( (LARGE_INTEGER *) &m_qwTicksPerSec );
        QueryPerformanceCounter   ( (LARGE_INTEGER *) &m_qwStartCount );
    }

    ~CTimer() {};

    const static int c_100nsTicksPerSecond = 10000000;  // 1 Tick per 100ns

	// this function normalizes the tick count from the perf counter to 
	// 100ns ticks
	UINT64 GetTime()
    {
        UINT64 qwCount;
		QueryPerformanceCounter((LARGE_INTEGER *) &qwCount);
		qwCount -= m_qwStartCount;
		qwCount *= c_100nsTicksPerSecond;
		qwCount = (qwCount+m_qwTicksPerSec/2) / m_qwTicksPerSec;
        return qwCount;
    }

	// return time in milli seconds
	float GetTimeMilliSec()
    {
        UINT64 qwCount;
		QueryPerformanceCounter((LARGE_INTEGER *) &qwCount);
		qwCount -= m_qwStartCount;
		double fMillSec = double(qwCount) * 1000. / double(m_qwTicksPerSec);  
		return (float)fMillSec;
	}

    // reset by setting start count to current count
    void Reset()
    {
        QueryPerformanceCounter((LARGE_INTEGER *) &m_qwStartCount );
    }

private:
    UINT64  m_qwTicksPerSec;
    UINT64  m_qwStartCount;
};

};
