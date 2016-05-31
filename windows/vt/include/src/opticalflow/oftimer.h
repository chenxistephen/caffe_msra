//+---------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation 2009.  All rights reserved.
//
//  Description:
//     Optical Flow: Timer Class
//
//  History:
//      2010/4/27-sbaker
//			Created
//
//----------------------------------------------------------------------------

#pragma once

// #define TIMER_ON

#include "vtcore.h"
using namespace vt;


class COFTimer : public CTimer
{
public:	
	COFTimer(wchar_t *pwcName) : CTimer()
	{
		m_i64Total = 0;
		m_iCount = 0;
		m_pwcName = pwcName;
	}

	void Start()
	{
		m_i64Start = GetTime();
	}

	void Stop()
	{
		m_i64Total += GetTime() - m_i64Start;
		m_iCount++;
	}

	void OutputTimes(int iCount)
	{
		double dTotalTime = double(m_i64Total)*0.0001;
		wprintf(L"%s %0.2lf\n", m_pwcName, dTotalTime/iCount);
	}

	static void OutputAllTimes(int iCount);

private:
	UINT64 m_i64Start;
	UINT64 m_i64Total;
	int m_iCount;
	wchar_t *m_pwcName;
};

#ifdef TIMER_ON
#ifdef TIMER_CPP
COFTimer cTimer_TOTAL(L"Total:  ");
COFTimer cTimer_RGB2Y(L"RGB2Y:  ");
COFTimer cTimer_PYRAMID(L"Pyramid:");
COFTimer cTimer_EXPAND(L"Expand: ");
COFTimer cTimer_W2(L"Warp:   ");
COFTimer cTimer_ERRORIM(L"Error:  ");
COFTimer cTimer_WEIGHT(L"Weights:");
COFTimer cTimer_DFLOW(L"DFLOW:  ");
COFTimer cTimer_LS(L"LS:     ");
COFTimer cTimer_LAP(L"LAP:    ");
COFTimer cTimer_SOR(L"SOR:    ");
COFTimer cTimer_DERIVS(L"Derivs: ");
void COFTimer::OutputAllTimes(int iCount)
{
	cTimer_TOTAL.OutputTimes(iCount);
	cTimer_RGB2Y.OutputTimes(iCount);
	cTimer_PYRAMID.OutputTimes(iCount);
	cTimer_EXPAND.OutputTimes(iCount);
	cTimer_DERIVS.OutputTimes(iCount);
	cTimer_W2.OutputTimes(iCount);
	cTimer_ERRORIM.OutputTimes(iCount);
	cTimer_WEIGHT.OutputTimes(iCount);
	cTimer_DFLOW.OutputTimes(iCount);
	cTimer_LS.OutputTimes(iCount);
	cTimer_LAP.OutputTimes(iCount);
	cTimer_SOR.OutputTimes(iCount);
}
#else
extern COFTimer cTimer_TOTAL;
extern COFTimer cTimer_RGB2Y;
extern COFTimer cTimer_PYRAMID;
extern COFTimer cTimer_EXPAND;
extern COFTimer cTimer_W2;
extern COFTimer cTimer_ERRORIM;
extern COFTimer cTimer_WEIGHT;
extern COFTimer cTimer_DFLOW;
extern COFTimer cTimer_LS;
extern COFTimer cTimer_LAP;
extern COFTimer cTimer_SOR;
extern COFTimer cTimer_DERIVS;
#define TIMER_START_TOTAL cTimer_TOTAL.Start();
#define TIMER_STOP_TOTAL cTimer_TOTAL.Stop();
#define TIMER_START_RGB2Y cTimer_RGB2Y.Start();
#define TIMER_STOP_RGB2Y cTimer_RGB2Y.Stop();
#define TIMER_START_PYRAMID cTimer_PYRAMID.Start();
#define TIMER_STOP_PYRAMID cTimer_PYRAMID.Stop();
#define TIMER_START_EXPAND cTimer_EXPAND.Start();
#define TIMER_STOP_EXPAND cTimer_EXPAND.Stop();
#define TIMER_START_W2 cTimer_W2.Start();
#define TIMER_STOP_W2 cTimer_W2.Stop();
#define TIMER_START_ERRORIM cTimer_ERRORIM.Start();
#define TIMER_STOP_ERRORIM cTimer_ERRORIM.Stop();
#define TIMER_START_WEIGHT cTimer_WEIGHT.Start();
#define TIMER_STOP_WEIGHT cTimer_WEIGHT.Stop();
#define TIMER_START_DFLOW cTimer_DFLOW.Start();
#define TIMER_STOP_DFLOW cTimer_DFLOW.Stop();
#define TIMER_START_LS cTimer_LS.Start();
#define TIMER_STOP_LS cTimer_LS.Stop();
#define TIMER_START_LAP cTimer_LAP.Start();
#define TIMER_STOP_LAP cTimer_LAP.Stop();
#define TIMER_START_SOR cTimer_SOR.Start();
#define TIMER_STOP_SOR cTimer_SOR.Stop();
#define TIMER_START_DERIVS cTimer_DERIVS.Start();
#define TIMER_STOP_DERIVS cTimer_DERIVS.Stop();
#endif
#else
#define TIMER_START_TOTAL
#define TIMER_STOP_TOTAL
#define TIMER_START_RGB2Y
#define TIMER_STOP_RGB2Y
#define TIMER_START_PYRAMID
#define TIMER_STOP_PYRAMID
#define TIMER_START_EXPAND
#define TIMER_STOP_EXPAND
#define TIMER_START_ERRORIM
#define TIMER_STOP_ERRORIM
#define TIMER_START_WEIGHT
#define TIMER_STOP_WEIGHT
#define TIMER_START_DFLOW
#define TIMER_STOP_DFLOW
#define TIMER_START_SOR
#define TIMER_STOP_SOR
#define TIMER_START_DERIVS
#define TIMER_STOP_DERIVS
#define TIMER_START_LS
#define TIMER_STOP_LS
#define TIMER_START_LAP
#define TIMER_STOP_LAP
#define TIMER_START_W2
#define TIMER_STOP_W2
#ifdef TIMER_CPP
void COFTimer::OutputAllTimes(int iCount) { UNREFERENCED_PARAMETER(iCount); }
#endif
#endif
