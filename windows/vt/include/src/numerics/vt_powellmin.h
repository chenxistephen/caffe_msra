//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Powell's minimization
//
//  History:
//      2006/11/20-swinder
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"

namespace vt {

HRESULT VtLineMinimize1D(float xmin, float xmid, float xmax, float &xrtn, float &fxrtn, 
                         HRESULT (*pFunc)(float x, float &fx, void *p), void *pUser);

// powell's multidimensional minimization without derivatives
HRESULT VtPowellSearch(CVecf &vStart, const CVecf &vDeltas, float &fRtn, 
                         HRESULT (*pFunc)(const CVecf &vParams, float &fValueRtn, void *pUserData),
                         void *pUser, int iMax = 20, float fFuncTolFrac = 0.0001f, float fXTolFrac = 0.0001f);

};
