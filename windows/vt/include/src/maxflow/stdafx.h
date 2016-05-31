// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define WIN32_LEAN_AND_MEAN     // Exclude rarely-used stuff from Windows headers
#define NOSERVICE               // Exculude more stuff ... see Windows.h for specifics
#define NOMCX
#define NOIME
#define NOSOUND
#define NOCOMM
#define NORPC

#if !defined(VT_WINRT)
#ifndef _WIN32_WINNT
#define _WIN32_WINNT   0x0400   // windows nt 4.0 or later
#endif
#ifndef _WIN32_WINDOWS
#define _WIN32_WINDOWS 0x0410   // windows 98 or later
#endif
#endif

#include <windows.h>
#include <objbase.h>


// internal
//#include <strsafe.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

