//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      stdafx for file IO project
//
//  History:
//      2004/11/8-swinder
//            Created
//
//------------------------------------------------------------------------

// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//
#pragma once

#define _FILEIO_STDAFX_H_

#include "vtfileio.h"

#include <mfidl.h>
#include <mferror.h>
#include <mfreadwrite.h>

//
// fileio needs some additional include files.  internal vt libraries redefine
// new to make sure it is the non throwing version - this can mess up system
// incldue files - so pop the macro state to make sure these other includes
// work.
// 
#pragma push_macro("OutputDebugStringA")
#undef OutputDebugStringA

#include <objbase.h>
#pragma warning(disable:4995)
#include <atlcomcli.h>    // for CComPtr
#pragma warning(default:4995)

#pragma pop_macro("OutputDebugStringA")

