//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Class to write an animated GIF file using Windows Imaging
//      Components (WIC).
//
//  History:
//      2012/04/25-ericsto
//          Created
//
//------------------------------------------------------------------------
#pragma once

#if defined(VT_WINRT)
#include "TempFileHelper.h"
#endif

namespace vt {

class CGifWriter : public IVideoDst
{
public:
	CGifWriter();
	~CGifWriter();

	int GetRepeatCount()
	{
		return m_iRepeatCount;
	}

	HRESULT SetRepeatCount(int repeatCount)
	{
		if (repeatCount < 0 || repeatCount > 0xFFFF)
		{
			return E_INVALIDARG;
		}
		m_iRepeatCount = repeatCount;
		return S_OK;
	}

	int GetColorDifferenceTolerance()
	{
		return m_iColorDifferenceTolerance;
	}

	HRESULT SetColorDifferenceTolerance(int colorDifferenceTolerance)
	{
		if (colorDifferenceTolerance < 0)
		{
			return E_INVALIDARG;
		}
		m_iColorDifferenceTolerance = colorDifferenceTolerance;
		return S_OK;
	}

	float GetMinimumFrameDuration()
	{
		return m_fMinimumFrameDuration;
	}

	HRESULT SetMinimumFrameDuration(float minimumFrameDuration)
	{
		if (minimumFrameDuration < 0)
		{
			return E_INVALIDARG;
		}
		m_fMinimumFrameDuration = minimumFrameDuration;
		return S_OK;
	}

	const vector<RGBAPix>* GetPalette()
	{
		return m_paletteColors;
	}

	HRESULT SetPalette(const vector<RGBAPix>* palette)
	{
		if (palette != NULL && (palette->size() < 1 || palette->size() > 256))
		{
			return E_INVALIDARG;
		}
		m_paletteColors = palette;
		return S_OK;
	}

	virtual HRESULT Clone(IVideoDst** ppClone);

	virtual HRESULT OpenFile(__in_z const wchar_t* pwcFilename, int iValidPixelWidth, int iValidPixelHeight,
							 VideoFormat eVideoFormat, int iFramesPerSecond= 15, float fBitsPerPixel = 8.0f,
							 double dAspectRatio = 1.0, int iInterlaceMode = MFVideoInterlace_Progressive);
#if defined(MAKE_DOXYGEN_WORK) || defined(VT_WINRT)
	virtual HRESULT OpenFile(Windows::Storage::IStorageFile^ storageFile, int iValidPixelWidth, int iValidPixelHeight,
							 VideoFormat eVideoFormat, int iFramesPerSecond = 15, float fBitsPerPixel = 8.0f,
							 double dAspectRatio = 1.0, int iInterlaceMode = MFVideoInterlace_Progressive);
#endif
	virtual HRESULT Close();

	virtual HRESULT WriteFrame(CImg& img);
	virtual HRESULT WriteFrame(CRGB32VideoImg& img);
	virtual HRESULT WriteFrame(CNV12VideoImg& img);

	virtual ULONG AddRef();
	virtual ULONG Release();

private:
	LONG                   m_iRefCount;
	bool                   m_bComStarted;
	bool                   m_bOpened;
	IWICImagingFactory*    m_pFactory;
	IWICStream*            m_pStream;
	IWICBitmapEncoder*     m_pEncoder;
	IWICPalette*           m_pPalette;
	int                    m_iRepeatCount;
	int                    m_iColorDifferenceTolerance;
	float                  m_fMinimumFrameDuration; // in seconds
	float                  m_fFrameDuration; // in seconds
	int                    m_iWidth;
	int                    m_iHeight;
	const vector<RGBAPix>* m_paletteColors;
	CRGBAImg               m_previousFrame;
	CRGBAImg               m_currentFrame;
	CRGBAImg               m_frameDifference;
	CRGBAImg               m_pendingFrame;
	float                  m_fPendingFrameDuration; // in seconds
	int                    m_iPendingFrameLeftOffset;
	int                    m_iPendingFrameTopOffset;
	bool                   m_bPendingFrameHasTransparency;
	int                    m_iInputFrameCount;
	int                    m_iOutputFrameCount;
#if defined(VT_WINRT)
	CTempFileHelper        m_tempFileHelper;
#endif

	HRESULT CalculateFrameDifference();
	HRESULT CommitPendingFrame();
	HRESULT CopyPendingFrameToEncoderFrame(IWICBitmapFrameEncode* pEncoderFrame);
};

}
