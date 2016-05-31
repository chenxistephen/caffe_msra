//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//     Main include file for the file IO library
//
//  History:
//      2004/11/8-swinder
//            Created
//
//------------------------------------------------------------------------
#pragma once

#pragma comment(lib, "vfw32.lib")
#pragma comment(lib, "mscms.lib")
#pragma comment(lib, "windowscodecs.lib")

#include <mfapi.h>
#pragma comment(lib, "mfuuid.lib")

#if defined(VT_WINRT)
#undef WINVER
#define WINVER _WIN32_WINNT_WIN8
#endif

#include "vtcore.h"

#include <vfw.h>
#include <icm.h>

#if !defined(VT_WINRT)
#pragma warning( push )
#pragma warning( disable : 6387 )
#include <thumbcache.h>
#pragma warning( pop )
#endif

#pragma warning( push )
#pragma warning( disable : 4005 )
#include <wincodec.h>
#include <wincodecsdk.h>
#pragma warning( pop )

#ifndef _FILEIO_STDAFX_H_

#include "..\src\fileio\vt_global.h"
#include "..\src\fileio\vt_iointerfaces.h"
#include "..\src\fileio\vt_io.h"
#include "..\src\fileio\GifWriter.h"
#include "..\src\fileio\psd.h"
#include "..\src\fileio\wicio.h"
#include "..\src\fileio\aviapi.h"
#include "..\src\fileio\vt_thumbnail.h"
#include "..\src\fileio\fileimport.h"
#include "..\src\fileio\tiler.h"

#endif

