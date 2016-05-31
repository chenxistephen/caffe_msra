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

#include "stdafx.h"
#include "GifWriter.h"

#include <wincodec.h>

#if defined(VT_WINRT)
using namespace Windows::Storage;
#endif

using namespace vt;

//===================================================================
// Helper functions
//===================================================================

template <class T> void SafeRelease(T*& pT)
{
	if (pT)
	{
		pT->Release();
		pT = NULL;
	}
}

//===================================================================
// CGifWriter class
//===================================================================

CGifWriter::CGifWriter()
	: m_iRefCount(1)
	, m_bOpened(false)
	, m_pFactory(NULL)
	, m_pStream(NULL)
	, m_pEncoder(NULL)
	, m_pPalette(NULL)
	, m_iRepeatCount(0)
	, m_iColorDifferenceTolerance(10)
	, m_fMinimumFrameDuration(0.06f) // frames shorter than 0.06 seconds are held for 0.1 seconds by older browsers
	, m_fFrameDuration(0)
	, m_iWidth(0)
	, m_iHeight(0)
	, m_paletteColors(NULL)
	, m_fPendingFrameDuration(0)
	, m_iPendingFrameLeftOffset(0)
	, m_iPendingFrameTopOffset(0)
	, m_bPendingFrameHasTransparency(false)
	, m_iInputFrameCount(0)
	, m_iOutputFrameCount(0)
{
	HRESULT hr = CoInitializeEx(NULL, COINIT_APARTMENTTHREADED);
	m_bComStarted = SUCCEEDED(hr);
}

CGifWriter::~CGifWriter()
{
	Close();
	SafeRelease(m_pFactory);
	if (m_bComStarted)
	{
		CoUninitialize();
	}
}

HRESULT CGifWriter::Clone(IVideoDst** /* ppClone */)
{
	// This class doesn't support cloning.
	return E_NOTIMPL;
}

HRESULT CGifWriter::OpenFile(__in_z const wchar_t* pwcFilename,
							 int iValidPixelWidth, int iValidPixelHeight,
							 VideoFormat eVideoFormat, int iFramesPerSecond, float fBitsPerPixel,
							 double /* dAspectRatio */, int /* iInterlaceMode */)
{
	// Check parameters.
	HRESULT hr = S_OK;
	if (m_bOpened)
	{
		VT_HR_EXIT(E_UNEXPECTED);
	}
	if (pwcFilename == NULL || iValidPixelWidth < 1 || iValidPixelHeight < 1 ||
		eVideoFormat != CIMG || iFramesPerSecond < 1 || fBitsPerPixel < 0.0f)
	{
		VT_HR_EXIT(E_INVALIDARG);
	}

	// Store the width, height, and frame duration.
	m_iWidth = iValidPixelWidth;
	m_iHeight = iValidPixelHeight;
	m_fFrameDuration = 1.0f / float(iFramesPerSecond);

	// Create a WIC imaging factory.
	if (m_pFactory == NULL)
	{
		VT_HR_EXIT(CSystem::CoCreateInstance(CLSID_WICImagingFactory, NULL, 
			CLSCTX_INPROC_SERVER, IID_IWICImagingFactory, (LPVOID*)&m_pFactory));
	}

	// Create a WIC stream.
	VT_HR_EXIT(m_pFactory->CreateStream(&m_pStream));
	VT_HR_EXIT(m_pStream->InitializeFromFilename(pwcFilename, GENERIC_WRITE));

	// Create a GIF encoder and tell it to use the stream.
	VT_HR_EXIT(m_pFactory->CreateEncoder(GUID_ContainerFormatGif, NULL, &m_pEncoder));
	VT_HR_EXIT(m_pEncoder->Initialize(m_pStream, WICBitmapEncoderNoCache));

	// If the caller has provided a palette, give it to the encoder.
	if (m_paletteColors != NULL)
	{
		// Convert to WIC colors.
		WICColor colors[256];
		UINT colorCount = VtMin((UINT)m_paletteColors->size(), 256u);
		for (int i = 0; i < (int)colorCount; i++)
		{
			const RGBAPix& color = (*m_paletteColors)[i];
			colors[i] = (color.a << 24) | (color.r << 16) | (color.g << 8) | color.b;
		}

		// Create and use the palette.
		VT_HR_EXIT(m_pFactory->CreatePalette(&m_pPalette));
		VT_HR_EXIT(m_pPalette->InitializeCustom(colors, colorCount));
		VT_HR_EXIT(m_pEncoder->SetPalette(m_pPalette));
	}

	// Write metadata.
	{
		// Create a metadata writer for the encoder.
		CComPtr<IWICMetadataQueryWriter> pQueryWriter;
		VT_HR_EXIT(m_pEncoder->GetMetadataQueryWriter(&pQueryWriter));

		// Write metadata indicating an animated GIF.
		UCHAR animatedGIFTag[12] = "NETSCAPE2.0";
		PROPVARIANT propVar;
		PropVariantInit(&propVar);
		propVar.vt = VT_UI1 | VT_VECTOR;
		propVar.caub.cElems = 11;
		propVar.caub.pElems = animatedGIFTag;
		VT_HR_EXIT(pQueryWriter->SetMetadataByName(L"/appext/application", &propVar));

		// Write metadata specifying the repeat count for the animation loop.
		if (m_iRepeatCount != 1)
		{
			UCHAR repeatCount[4];
			repeatCount[0] = 3;	// must be 3
			repeatCount[1] = 1; // must be 1
			repeatCount[2] = UCHAR(m_iRepeatCount & 0xFF); // low byte of repeat count (or zero for infinite)
			repeatCount[3] = UCHAR(m_iRepeatCount >> 8); // high byte of repeat count (or zero for infinite)
			propVar.vt = VT_UI1 | VT_VECTOR;
			propVar.caub.cElems = 4;
			propVar.caub.pElems = repeatCount;
			VT_HR_EXIT(pQueryWriter->SetMetadataByName(L"/appext/data", &propVar));
		}
	}

Exit:
	if (hr != S_OK)
	{
		Close();
	}
	else
	{
		m_bOpened = true;
	}

	return hr;
}

#if defined(VT_WINRT)
HRESULT CGifWriter::OpenFile(IStorageFile^ storageFile, int iValidPixelWidth, int iValidPixelHeight,
							 VideoFormat eVideoFormat, int iFramesPerSecond, float fBitsPerPixel,
							 double dAspectRatio, int iInterlaceMode)
{
	VT_HR_BEGIN();

	VT_HR_EXIT(m_tempFileHelper.Initialize(storageFile, FileAccessMode::ReadWrite));

	IStorageFile^ tempStorageFile = m_tempFileHelper.GetTempStorageFile();

	VT_HR_EXIT(OpenFile(tempStorageFile->Path->Data(), iValidPixelWidth, iValidPixelHeight,
		eVideoFormat, iFramesPerSecond, fBitsPerPixel, dAspectRatio, iInterlaceMode));

	VT_HR_END();
}
#endif

HRESULT CGifWriter::Close()
{
	VT_HR_BEGIN();

	// Finish writing.
	if (m_pEncoder)
	{
		// If there's a pending frame, commit it now.
		if (m_pendingFrame.IsValid())
		{
			VT_HR_EXIT(CommitPendingFrame());
		}

		// TODO: Only commit if frames were added.
		VT_HR_EXIT(m_pEncoder->Commit());
	}

	SafeRelease(m_pEncoder);
	SafeRelease(m_pStream);
	SafeRelease(m_pPalette);

#if defined(VT_WINRT)
	VT_HR_EXIT(m_tempFileHelper.Finalize());
#endif

	VT_HR_EXIT_LABEL();

	SafeRelease(m_pEncoder);
	SafeRelease(m_pStream);
	SafeRelease(m_pPalette);
	m_fFrameDuration = 0;
	m_iWidth = 0;
	m_iHeight = 0;
	m_paletteColors = NULL;
	m_previousFrame.Deallocate();
	m_currentFrame.Deallocate();
	m_frameDifference.Deallocate();
	m_pendingFrame.Deallocate();
	m_fPendingFrameDuration = 0;
	m_iPendingFrameLeftOffset = 0;
	m_iPendingFrameTopOffset = 0;
	m_bPendingFrameHasTransparency = false;
	m_iInputFrameCount = 0;
	m_iOutputFrameCount = 0;

	return hr;
}

HRESULT CGifWriter::WriteFrame(CImg& img)
{
	VT_HR_BEGIN();

	if (!m_bOpened)
	{
		VT_HR_EXIT(E_UNEXPECTED);
	}
	if (img.Width() != m_iWidth || img.Height() != m_iHeight)
	{
		VT_HR_EXIT(E_INVALIDARG);
	}

	// See if we should drop the frame in order to meet the minimum frame duration constraint.
	if (m_fFrameDuration < m_fMinimumFrameDuration &&
		m_iInputFrameCount++ * m_fFrameDuration / m_fMinimumFrameDuration < m_iOutputFrameCount)
	{
		return S_OK;
	}
	m_iOutputFrameCount++;

	// First, convert the current image to an RGBA byte image.
	VT_HR_EXIT(VtConvertImageToRGBA(m_currentFrame, img));
	
	// See if we have a previous frame.
	if (m_previousFrame.IsValid())
	{
		// Yes, there's a previous frame.
		// Create the pending frame from the difference between the current frame and the previous one.
		VT_HR_EXIT(CalculateFrameDifference());
	}
	else
	{
		// There's no previous frame, so just copy the current frame to the previous and pending frames.
		VT_HR_EXIT(m_currentFrame.CopyTo(m_previousFrame));
		VT_HR_EXIT(m_previousFrame.Share(m_pendingFrame));
		m_fPendingFrameDuration = VtMax(m_fFrameDuration, m_fMinimumFrameDuration);
		m_iPendingFrameLeftOffset = 0;
		m_iPendingFrameTopOffset = 0;
		m_bPendingFrameHasTransparency = false;
	}

	VT_HR_END();
}

HRESULT CGifWriter::WriteFrame(CRGB32VideoImg& /* img */)
{
	// This class doesn't support CRGB32VideoImg frames.
	return E_NOTIMPL;
}

HRESULT CGifWriter::WriteFrame(CNV12VideoImg& /* img */)
{
	// This class doesn't support CNV12VideoImg frames.
	return E_NOTIMPL;
}

ULONG CGifWriter::AddRef()
{
	return (ULONG)InterlockedIncrement(&m_iRefCount);
}

ULONG CGifWriter::Release()
{
	if (InterlockedDecrement(&m_iRefCount) == 0)
	{
		delete this;
		return 0;
	}
	return (ULONG)m_iRefCount;
}

HRESULT CGifWriter::CalculateFrameDifference()
{
	VT_HR_BEGIN();

	// Make sure the frame difference image is initialized.
	VT_HR_EXIT(m_frameDifference.Create(m_iWidth, m_iHeight));
	VT_HR_EXIT(m_frameDifference.Fill(RGBAPix(0, 0, 0, 0)));

	// Copy only the pixels that differ.
	int left = INT_MAX;
	int right = INT_MIN;
	int top = INT_MAX;
	int bottom = INT_MIN;
	for (int y = 0; y < m_iHeight; y++)
	{
		const RGBAPix* previousFramePtr = m_previousFrame.Ptr(y);
		const RGBAPix* currentFramePtr = m_currentFrame.Ptr(y);
		RGBAPix* frameDifferencePtr = m_frameDifference.Ptr(y);
		for (int x = 0; x < m_iWidth; x++)
		{
			// See if the previous frame and the current frame differ at this pixel.
			int dr = currentFramePtr->r - previousFramePtr->r;
			int dg = currentFramePtr->g - previousFramePtr->g;
			int db = currentFramePtr->b - previousFramePtr->b;
			int da = currentFramePtr->a - previousFramePtr->a;
			if (abs(dr) + abs(dg) + abs(db) + abs(da) > m_iColorDifferenceTolerance)
			{
				// Yes, they differ.  Store the current frame's pixel value and update the bounds.
				*frameDifferencePtr = *currentFramePtr;
				left = VtMin(left, x);
				right = VtMax(right, x + 1);
				top = VtMin(top, y);
				bottom = VtMax(bottom, y + 1);
			}
			previousFramePtr++;
			currentFramePtr++;
			frameDifferencePtr++;
		}
	}

	// See if any pixels differed.
	if (right > left && bottom > top)
	{
		// Yes, the current frame differs from the previous one.
		// Write out the pending frame.
		VT_HR_EXIT(CommitPendingFrame());

		// Create the next pending frame from the minimal region of differing pixels.
		m_frameDifference.CopyTo(m_pendingFrame, CRect(left, top, right, bottom));
		m_fPendingFrameDuration = VtMax(m_fFrameDuration, m_fMinimumFrameDuration);
		m_iPendingFrameLeftOffset = left;
		m_iPendingFrameTopOffset = top;
		m_bPendingFrameHasTransparency = true;

		// The current frame becomes the previous frame.
		m_currentFrame.CopyTo(m_previousFrame);
	}
	else
	{
		// The current frame is identical to the previous one.
		// Increase the duration of the pending frame.
		m_fPendingFrameDuration += VtMax(m_fFrameDuration, m_fMinimumFrameDuration);

		// Check for a frame duration that exceeds the capacity of the GIF format
		// (65535 hundredths of a second = 655 seconds = 10.9 minutes).
		if (F2I(100 * m_fPendingFrameDuration) > 0xFFFF)
		{
			// Go back to the previous duration and commit the pending frame.
			m_fPendingFrameDuration -= VtMax(m_fFrameDuration, m_fMinimumFrameDuration);
			VT_HR_EXIT(CommitPendingFrame());

			// Create a single-pixel transparent pending frame to use for the remaining time.
			VT_HR_EXIT(m_pendingFrame.Create(1, 1));
			VT_HR_EXIT(m_pendingFrame.Fill(RGBAPix(0, 0, 0, 0)));
			m_fPendingFrameDuration = VtMax(m_fFrameDuration, m_fMinimumFrameDuration);
			m_iPendingFrameLeftOffset = 0;
			m_iPendingFrameTopOffset = 0;
			m_bPendingFrameHasTransparency = true;
		}
	}

	VT_HR_END();
}

HRESULT CGifWriter::CommitPendingFrame()
{
	VT_HR_BEGIN();

	// Create and initialize an encoder frame.
	CComPtr<IWICBitmapFrameEncode> pEncoderFrame;
	VT_HR_EXIT(m_pEncoder->CreateNewFrame(&pEncoderFrame, NULL));
	VT_HR_EXIT(pEncoderFrame->Initialize(NULL));

	// Create a metadata writer for the encoder frame.
	CComPtr<IWICMetadataQueryWriter> pFrameQueryWriter;
	VT_HR_EXIT(pEncoderFrame->GetMetadataQueryWriter(&pFrameQueryWriter));

	// Write metadata specifying the frame delay (in hundredths of a second).
	float frameDuration = VtMax(m_fMinimumFrameDuration, m_fPendingFrameDuration);
	PROPVARIANT propVar;
	PropVariantInit(&propVar);
	propVar.vt = VT_UI2;
	propVar.uiVal = UINT16(F2I(100 * frameDuration));
	VT_HR_EXIT(pFrameQueryWriter->SetMetadataByName(L"/grctlext/Delay", &propVar));

	// Write metadata specifying the transparent color index, if needed.
	if (m_bPendingFrameHasTransparency)
	{
		propVar.vt = VT_BOOL;
		propVar.boolVal = true;
		VT_HR_EXIT(pFrameQueryWriter->SetMetadataByName(L"/grctlext/TransparencyFlag", &propVar));

		propVar.vt = VT_UI1;
		propVar.bVal = 0;
		VT_HR_EXIT(pFrameQueryWriter->SetMetadataByName(L"/grctlext/TransparentColorIndex", &propVar));
	}

	// Write metadata specifying the left offset.
	if (m_iPendingFrameLeftOffset > 0)
	{
		propVar.vt = VT_UI2;
		propVar.uiVal = UINT16(m_iPendingFrameLeftOffset);
		VT_HR_EXIT(pFrameQueryWriter->SetMetadataByName(L"/imgdesc/Left", &propVar));
	}

	// Write metadata specifying the top offset.
	if (m_iPendingFrameTopOffset > 0)
	{
		propVar.vt = VT_UI2;
		propVar.uiVal = UINT16(m_iPendingFrameTopOffset);
		VT_HR_EXIT(pFrameQueryWriter->SetMetadataByName(L"/imgdesc/Top", &propVar));
	}

	// Copy the pending frame's pixels to the encoder frame.
	VT_HR_EXIT(CopyPendingFrameToEncoderFrame(pEncoderFrame));

	// Commit the encoder frame.
	VT_HR_EXIT(pEncoderFrame->Commit());

	VT_HR_END();
}

HRESULT CGifWriter::CopyPendingFrameToEncoderFrame(IWICBitmapFrameEncode* pEncoderFrame)
{
	VT_HR_BEGIN();

	// Get information about the pending frame.
	const CImgInfo& imgInfo = m_pendingFrame.GetImgInfo();
	int width = imgInfo.width;
	int height = imgInfo.height;

	// Create a WIC bitmap.
	CComPtr<IWICBitmap> pBitmap;
	VT_HR_EXIT(m_pFactory->CreateBitmap(width, height, GUID_WICPixelFormat32bppPBGRA,
		WICBitmapCacheOnLoad, &pBitmap));

	// Lock the entire bitmap.
	WICRect lockRect = { 0, 0, width, height };
	CComPtr<IWICBitmapLock> pLock;
	VT_HR_EXIT(pBitmap->Lock(&lockRect, WICBitmapLockWrite, &pLock));

	// Obtain the pixel buffer from the lock.
	UINT cbStride = 0;
	UINT cbBufferSize = 0;
	BYTE* pData;
	VT_HR_EXIT(pLock->GetStride(&cbStride));
	VT_HR_EXIT(pLock->GetDataPointer(&cbBufferSize, &pData));

#ifdef _DEBUG
	// Check the stride and buffer size.
	UINT uByteCount = UINT(width * imgInfo.PixSize());
	VT_ASSERT(cbStride >= uByteCount);
	VT_ASSERT(cbBufferSize >= uByteCount * height);
#endif

	// Wrap the WIC pixel buffer in a VT image, then copy the pending frame's pixels.
	CImg imgDst;
	VT_HR_EXIT(imgDst.Create(pData, width, height, cbStride, imgInfo.type));
	VT_HR_EXIT(m_pendingFrame.CopyTo(imgDst));

	// Release the lock.
	pLock = NULL;

	// Decide which palette to use.
	CComPtr<IWICPalette> pPalette;
	if (m_pPalette != NULL)
	{
		// Use the caller's specified palette.
		pPalette = m_pPalette;
	}
	else
	{
		// Use WIC's color quantization algorithm to choose 256 colors from the bitmap.
		VT_HR_EXIT(m_pFactory->CreatePalette(&pPalette));
		VT_HR_EXIT(pPalette->InitializeFromBitmap(pBitmap, 256, m_bPendingFrameHasTransparency));

		// If the frame has transparency, make sure the first color in the palette
		// is the transparent color.  This may require swapping the first color
		// with another color in the palette (typically the last color).
		if (m_bPendingFrameHasTransparency)
		{
			WICColor colors[256];
			UINT colorCount = 0;
			VT_HR_EXIT(pPalette->GetColors(ARRAYSIZE(colors), colors, &colorCount));
			VT_ASSERT(colorCount <= 256);
			if (colorCount <= 256 && colors[0] != 0)
			{
				// Find the index of the transparent color.
				int transparentIndex = colorCount - 1;
				while (transparentIndex > 0 && colors[transparentIndex] != 0)
				{
					transparentIndex--;
				}

				// Swap the first color with the transparent one.
				colors[transparentIndex] = colors[0];
				colors[0] = 0;
				VT_HR_EXIT(pPalette->InitializeCustom(colors, colorCount));
			}
		}
	}

	// Create and initialize a format converter.
	CComPtr<IWICFormatConverter> pFormatConverter;
	VT_HR_EXIT(m_pFactory->CreateFormatConverter(&pFormatConverter));
	VT_HR_EXIT(pFormatConverter->Initialize(
		pBitmap,						   // input bitmap
		GUID_WICPixelFormat8bppIndexed,	   // destination pixel format
		WICBitmapDitherTypeErrorDiffusion, // dithering algorithm
		pPalette,                          // destination palette 
		50.0,                              // alpha threshold percent
		WICBitmapPaletteTypeCustom         // palette translation type
		));

	// Add the frame to the encoder.
	VT_HR_EXIT(pEncoderFrame->SetSize(UINT(width), UINT(height)));
	VT_HR_EXIT(pEncoderFrame->WriteSource(pFormatConverter, NULL));

	VT_HR_END();
}
