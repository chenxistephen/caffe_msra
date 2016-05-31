//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      HRESULTS etc
//
//  History:
//      2004/11/08-swinder
//			Created
//
//------------------------------------------------------------------------

#pragma once

#include "vt_basetypes.h"
#include "vt_dbg.h"

namespace vt {

#define S_EOF			_HRESULT_TYPEDEF_(0x0fff0001L)
#define E_NOINIT		_HRESULT_TYPEDEF_(0x8fff0009L)
#define E_READFAILED	_HRESULT_TYPEDEF_(0x8fff0002L)
#define E_WRITEFAILED	_HRESULT_TYPEDEF_(0x8fff0003L)
#define E_BADFORMAT		_HRESULT_TYPEDEF_(0x8fff0004L)
#define E_TYPEMISMATCH	_HRESULT_TYPEDEF_(0x8fff0005L)
#define E_NOTFOUND		_HRESULT_TYPEDEF_(0x8fff0006L)
#define E_NOCOMPRESSOR	_HRESULT_TYPEDEF_(0x8fff0007L)
#define E_TOOCOMPLEX	_HRESULT_TYPEDEF_(0x8fff0008L)
#define E_INVALIDSRC	_HRESULT_TYPEDEF_(0x8fff0010L)
#define E_INVALIDDST	_HRESULT_TYPEDEF_(0x8fff0011L)

// Error macros
#pragma warning (disable: 4102) // no warning when function has no _EXIT instances
#define VT_HR_BEGIN() HRESULT hr = S_OK; {
#define VT_HR_EXIT_LABEL() } Exit:
#define VT_HR_END() } Exit: return hr;

#define VT_HR_RET(expression) \
if FAILED(hr = (expression)) \
{ VT_DEBUG_LOG_HR(hr); return hr;} \
else hr

#define VT_HR_EXIT(expression) \
if FAILED(hr = (expression)) \
{ VT_DEBUG_LOG_HR(hr); goto Exit;} \
else hr

#define VT_PTR_OOM_EXIT(expression) \
if (NULL == (expression)) \
{hr = E_OUTOFMEMORY; VT_DEBUG_LOG_HR(hr); goto Exit;} \
else hr

#define VT_PTR_EXIT(expression) \
if (NULL == (expression)) \
{hr = E_POINTER; VT_DEBUG_LOG_HR(hr); goto Exit;} \
else hr

/// \ingroup error
/// <summary> Return common HRESULT as string </summary>
/// <param name="hr"> HRESULT from a VisionTools operation </param>
/// <param name="buf"> Buffer to write the error text to </param>
/// <param name="numBufElem"> Buffer size, in characters </param>
/// <returns> Pointer to the beginning of the text buffer </returns>
const wchar_t* VtErrorToString(
    HRESULT hr, 
    __in_ecount(numBufElem) wchar_t* buf, 
    int numBufElem);

class CErrorBase
{
public:
	CErrorBase() : m_hr(NOERROR) {}

	HRESULT GetError() const { return m_hr; }
	bool    IsError()  const { return FAILED(m_hr); }
	HRESULT SetHR(HRESULT hr) { m_hr = hr; return hr; } // for setting S_EOF
	HRESULT SetError(HRESULT hr) { if(FAILED(hr)) m_hr = hr; return hr; } // only set error on failure
	HRESULT ClearError()         { m_hr = NOERROR; return m_hr; }

protected:
    HRESULT ClearError()         const { m_hr = NOERROR; return m_hr; }
	HRESULT SetError(HRESULT hr) const { m_hr = hr; return hr; }

private:
	mutable HRESULT m_hr;
};

};
