//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Debugging support
//
//  History:
//      2004/11/08-mattu
//			Created
//
//------------------------------------------------------------------------
#pragma once

#if !defined(VT_WINRT)

#include "vtcommon.h"
#include "vt_image.h"

namespace vt {

#define DEFAULT_PELS_METER	3937
#define VTWIDTHBYTES(bits)    ((((bits) + 31)>>5)<<2)

DECLARE_HANDLE(HDIB);
typedef LPBITMAPINFOHEADER	PDIB;

#define DIB_HEADER_MARKER   ((WORD) ('M' << 8) | 'B') 

#define DIBCMAP_GRAYSCALE	0
#define DIBCMAP_PSEUDO1		1
#define DIBCMAP_COLDHOT		2
#define DIBCMAP_COLORBARS   3
#define DIBCMAP_PSEUDO2     4 // blue - magenta - red - yellow - green - cyan


////////////////////////////////////////////////////////////////////////
// CDib - Dib wrapper class
////////////////////////////////////////////////////////////////////////

class CDib
{
public:
	CDib();
	~CDib();

	HRESULT Create32Bit(int iWidth, int iHeight);
	HRESULT Create24Bit(int iWidth, int iHeight);
	HRESULT Create8Bit(int iWidth, int iHeight);
	HRESULT Create(BITMAPINFOHEADER *pbmi);				// alloc and copy
	HRESULT CreateNoUpdate(BITMAPINFOHEADER *pbmi);		// just alloc the memory
	HRESULT Update(BITMAPINFOHEADER *pbmi);				// copy pix data over

	void Destroy();
	bool IsValid() const { return m_pDib!=NULL; }

	HDIB GetDibHandle() const { return m_hDib; }
	PDIB GetDibInfo() const { return m_pDib; }

	size_t GetDibHeaderSize()
		{ return (m_pDib==NULL) ? 0 : (size_t)m_iOffsetToBits; }
	int Width() const
		{ return (m_pDib==NULL) ? 0 : (int)m_pDib->biWidth; }
	int Height() const
		{ return (m_pDib==NULL) ? 0 : (int)m_pDib->biHeight; }
	size_t StrideBytes() const
		{ return (m_pDib==NULL) ? 0 : VTWIDTHBYTES(m_pDib->biWidth * m_pDib->biBitCount); }
	size_t ElSize() const
		{ return (m_pDib==NULL) ? 0 : (m_pDib->biBitCount + 7)>>3; }
	Byte *GetDibBits() const
		{ return (m_pDib==NULL) ? NULL : ((Byte *)m_pDib + m_iOffsetToBits); }
	RGBQUAD *GetCMapBits() const
		{ return (m_pDib==NULL) ? NULL : (RGBQUAD *)((Byte *)m_pDib + sizeof(BITMAPINFOHEADER)); }

	void SetCMap(int iType);
	void CopyCMap(CDib *pDibSrc);
	void SetCMapColour(int iIndex, RGBAPix &rgb);
	RGBAPix GetCMapColour(int iIndex) const;

	HRESULT SetDibBits(Byte *puchPix, int iSrcRgnX, int iSrcRgnY,
		int iSrcRgnW, int iSrcRgnH, int iSrcH, int iSrcStride, int iDstX, int iDstY);

	HRESULT GetPixels(CImg &cImg);
	HRESULT SetPixels(const CImg &cImg, int iDstX, int iDstY);

	void GetPixel(int iX, int iY, RGBAPix &cRGB) const;
	void SetPixel(int iX, int iY, RGBAPix &cRGB);

	HRESULT Load(FILE *fp);
	HRESULT Save(FILE *fp) const;

private:
	HDIB m_hDib;
	PDIB m_pDib;
	int m_iOffsetToBits;
};

};

#endif
