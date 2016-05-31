#pragma once

#ifndef _WINDOWS_

#define UNREFERENCED_PARAMETER(P) (P)
	
typedef long HRESULT;

#define S_OK     ((HRESULT)0L)
#define S_FALSE  ((HRESULT)1L)
#define NOERROR  0

#define _HRESULT_TYPEDEF_(_sc) ((HRESULT)_sc)

#define SUCCEEDED(hr) (((HRESULT)(hr)) >= 0)
#define FAILED(hr)    (((HRESULT)(hr)) < 0)

#define E_NOTIMPL        ((HRESULT)0x80000001L)
#define E_OUTOFMEMORY    ((HRESULT)0x80000002L)
#define E_INVALIDARG     ((HRESULT)0x80000003L)
#define E_NOINTERFACE    ((HRESULT)0x00000004L)
#define E_POINTER        ((HRESULT)0x80000005L)
#define E_HANDLE         ((HRESULT)0x80000006L)
#define E_ABORT          ((HRESULT)0x80000007L)
#define E_FAIL           ((HRESULT)0x80000008L)
#define E_ACCESSDENIED   ((HRESULT)0x80000009L)
#define E_UNEXPECTED     ((HRESULT)0x8000FFFFL)

typedef wchar_t WCHAR;
typedef unsigned long ULONG;
typedef unsigned int UINT;
typedef unsigned long DWORD;

typedef struct tagRECT
{
    long left;
    long top;
    long right;
    long bottom;
} RECT;

typedef struct tagPOINT
{
    long x;
    long y;
} POINT;

typedef struct tagSIZE
{
    long cx;
    long cy;
} SIZE;

typedef struct _GUID {
    unsigned long  Data1;
    unsigned short Data2;
    unsigned short Data3;
    unsigned char  Data4[ 8 ];
} GUID;

//
// SAL
//

#ifndef IN
#define IN
#endif

#ifndef OUT
#define OUT
#endif

#ifndef OPTIONAL
#define OPTIONAL
#endif

#define _Printf_format_string_

/// 
/// int ptr
///

#ifdef VT_GCC
#include <stdint.h>
#define INT_PTR intptr_t
#endif
#else
#if defined(_WIN64)
    typedef __int64 INT_PTR;
#else
    typedef __w64 int INT_PTR;
#endif
#endif
