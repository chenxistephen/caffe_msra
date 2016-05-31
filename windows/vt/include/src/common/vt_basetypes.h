//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      VisionTools typedefs
//
//  History:
//      2004/11/8-swinder
//			Created
//
//------------------------------------------------------------------------
#pragma once

#ifdef _WIN32_WINNT
#undef _WIN32_WINNT
#endif

#ifdef VT_WINRT
#define _WIN32_WINNT _WIN32_WINNT_WIN8
#else
#define _WIN32_WINNT _WIN32_WINNT_WIN7
#endif

#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
    #define NOMINMAX
#endif
#include <Windows.h>
#else
#include "vt_minwin.h"
#endif

#include <stddef.h>

// no stdint.h in Windows Win8/WP8 build
#if defined(MSRVT_WINDOWS_BUILD)

#define VT_STDCALL __stdcall

namespace vt {
typedef UINT8 Byte;
typedef UINT16 UInt16;
typedef INT32 Int32;
typedef UINT32 UInt32;
typedef INT64 Int64;
typedef UINT64 UInt64;
#define INT32_MIN (-0x7fffffff - 1)
#define INT32_MAX 0x7fffffff
#define UINT32_MAX 0xffffffff
}

typedef INT8   int8_t;  
typedef UINT8  uint8_t;  
typedef UINT16 uint16_t;
typedef INT16  int16_t;  
typedef UINT32 uint32_t;
typedef INT32  int32_t;  
typedef UINT64 uint64_t;
typedef INT64  int64_t;

#include <type_traits>
#else

#define VT_STDCALL

#include <stdint.h>

namespace vt {
typedef uint8_t Byte;
typedef uint16_t UInt16;
typedef int32_t Int32;
typedef uint32_t UInt32;
typedef int64_t Int64;
typedef uint64_t UInt64;

#ifndef INT32_MIN
#define INT32_MIN (-0x7fffffff - 1)
#endif

#ifndef INT32_MAX
#define INT32_MAX 0x7fffffff
#endif

#ifndef UINT32_MAX
#define UINT32_MAX 0xffffffff
#endif

}

#include <type_traits>

#endif

namespace vt {

/// <summary> Half precision (16 bit) float </summary>
struct HALF_FLOAT
{
	UInt16 v; ///< half-float bits
};

template <class T> T VtMin(T a, T b) { return (a<b) ? a : b; }
template <class T> T VtMax(T a, T b) { return (a>b) ? a : b; }

template <class T> T VtClamp(T a, T minv, T maxv) 
{ return VtMin(VtMax(a, minv), maxv); }

};

// dummy macros for SAL annotations on non-Windows compilers
#if defined(VT_GCC)
#ifndef __in_z
#define __in_z
#endif
#ifndef __out_z
#define __out_z
#endif
#ifndef __in_ecount
#define __in_ecount(__param__)
#endif
#ifndef __out_ecount
#define __out_ecount(__param__)
#endif
#endif

