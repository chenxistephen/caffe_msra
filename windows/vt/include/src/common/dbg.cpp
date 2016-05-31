//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Debugging support
//
//  History:
//      2004/11/08-mattu
//          Created
//		2011/06/09-wkienzle
//			Simplified
//
//------------------------------------------------------------------------

#include "vt_dbg.h"
#include "vt_string.h"

using namespace vt;

void vt::VtDebugLog(const char* msg, ...)
{
	va_list args;
	va_start(args, msg);
	string_b<512> buf;
#ifdef _MSC_VER
	vsprintf_s(buf.get_buffer(), 512, msg, args);
#else
	vsprintf(buf.get_buffer(), msg, args);
#endif	
	va_end(args);
#ifdef _MSC_VER
	OutputDebugStringA(buf.get_constbuffer());	
#else
#warning DebugLog not implemented yet
#endif
}

void vt::VtDebugLog(const wchar_t* msg, ...)
{
	va_list args;
	va_start(args, msg);
	wstring_b<512> buf;
#ifdef _MSC_VER
	vswprintf_s(buf.get_buffer(), 512, msg, args);
#else
	vswprintf(buf.get_buffer(), 512, msg, args);
#endif	
	va_end(args);
#ifdef _MSC_VER
	OutputDebugStringW(buf.get_constbuffer());
#else
#warning DebugLog not implemented yet
#endif
}
