//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Debugging support.
//
//  History:
//      2004/11/08-mattu
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vt_basetypes.h"

#if defined(VT_GCC)
#include <assert.h>
#define _ASSERT() VT_ASSERT()
#else
#include <crtdbg.h>
#endif

// warning format for Visual Studio
#define __STR2__(x) #x
#define __STR1__(x) __STR2__(x)
#define __LOC__ __FILE__ "("__STR1__(__LINE__)")"
#define __WARNING__ __LOC__" : warning: "
// #pragma message(__LOC__"Missing implementation for this target platform")

namespace vt {

#ifdef _DEBUG
#define VT_ASSERT(exp) _ASSERT(exp)
#define VT_VERIFY(exp) _ASSERT(exp)
#else
#define VT_ASSERT(exp) 
#define VT_VERIFY(exp) exp
#endif

#if defined(_DEBUG) || defined(VT_BOUNDS_CHECKING)
#define VT_ASSERT_BOUNDS(exp) _ASSERT(exp)
#else
#define VT_ASSERT_BOUNDS(exp) 
#endif

#ifdef __analysis_assume
#define ANALYZE_ASSUME(expr) __analysis_assume(expr)
#else
#define ANALYZE_ASSUME(expr)
#endif

void VtDebugLog(const char* msg, ...);
void VtDebugLog(const wchar_t* msg, ...);

#ifdef _DEBUG
#define VT_DEBUG_LOG_HR(hr) vt::VtDebugLog(__LOC__ " : HRESULT = 0x%08x\n", hr)
#define VT_DEBUG_LOG(fmt, ...) vt::VtDebugLog(fmt, __VA_ARGS__)
#else
#define VT_DEBUG_LOG_HR(hr)
#define VT_DEBUG_LOG(fmt, ...) 
#endif
};