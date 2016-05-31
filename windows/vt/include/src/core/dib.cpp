//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Dib class for representing Bitmaps and converting between Windows
//       and VisionTools formats
//
//  History:
//      2004/11/08-swinder
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#if !defined(VT_WINRT)

#include "vt_dib.h"
#include "vt_convert.h"
#include "vt_utils.h"

using namespace vt;

////////////////////////////////////////////////////////////////////////
// CDib
////////////////////////////////////////////////////////////////////////

CDib::CDib() : m_hDib(NULL), m_pDib(NULL), m_iOffsetToBits(0)
{
}

CDib::~CDib()
{
    Destroy();
}

// create 32-bit Dib
HRESULT CDib::Create32Bit(int iWidth, int iHeight)
{
    if(iWidth<=0 || iHeight <=0)
        return E_INVALIDARG;

    Destroy();

    DWORD dwBitsSize = VTWIDTHBYTES(iWidth * 32) * (DWORD)iHeight;
    DWORD dwTotalSize = sizeof(BITMAPINFOHEADER) + dwBitsSize;

    m_hDib = (HDIB)::GlobalAlloc(GMEM_MOVEABLE | GMEM_ZEROINIT, dwTotalSize);
    if(m_hDib==NULL)
        return E_OUTOFMEMORY;
    
    m_pDib = (PDIB)::GlobalLock((HGLOBAL)m_hDib);
	if(m_pDib==NULL)
        return E_FAIL;

    m_pDib->biSize = sizeof(BITMAPINFOHEADER);
    m_pDib->biWidth = (LONG)iWidth;
    m_pDib->biHeight = (LONG)iHeight;
    m_pDib->biPlanes = 1;
    m_pDib->biBitCount = 32;
    m_pDib->biCompression = BI_RGB;
    m_pDib->biSizeImage = dwBitsSize;
    m_pDib->biXPelsPerMeter = DEFAULT_PELS_METER;
    m_pDib->biYPelsPerMeter = DEFAULT_PELS_METER;
    m_pDib->biClrUsed = 0;
    m_pDib->biClrImportant = 0;

    m_iOffsetToBits = sizeof(BITMAPINFOHEADER);

    return NOERROR;
}

// create 24-bit Dib
HRESULT CDib::Create24Bit(int iWidth, int iHeight)
{
    if(iWidth<=0 || iHeight <=0)
        return E_INVALIDARG;

    Destroy();

    DWORD dwBitsSize = VTWIDTHBYTES(iWidth * 24) * (DWORD)iHeight;
    DWORD dwTotalSize = sizeof(BITMAPINFOHEADER) + dwBitsSize;

    m_hDib = (HDIB)::GlobalAlloc(GMEM_MOVEABLE | GMEM_ZEROINIT, dwTotalSize);
    if(m_hDib==NULL)
        return E_OUTOFMEMORY;
    
    m_pDib = (PDIB)::GlobalLock((HGLOBAL)m_hDib);
    if(m_pDib==NULL)
        return E_FAIL;

    m_pDib->biSize = sizeof(BITMAPINFOHEADER);
    m_pDib->biWidth = (LONG)iWidth;
    m_pDib->biHeight = (LONG)iHeight;
    m_pDib->biPlanes = 1;
    m_pDib->biBitCount = 24;
    m_pDib->biCompression = BI_RGB;
    m_pDib->biSizeImage = dwBitsSize;
    m_pDib->biXPelsPerMeter = DEFAULT_PELS_METER;
    m_pDib->biYPelsPerMeter = DEFAULT_PELS_METER;
    m_pDib->biClrUsed = 0;
    m_pDib->biClrImportant = 0;

    m_iOffsetToBits = sizeof(BITMAPINFOHEADER);

    return NOERROR;
}

// create 8-bit colourmap Dib
HRESULT CDib::Create8Bit(int iWidth, int iHeight)
{
    if(iWidth<=0 || iHeight <=0)
        return E_INVALIDARG;

    Destroy();

    DWORD dwMapSize = 256 * sizeof(RGBQUAD);
    DWORD dwBitsSize = VTWIDTHBYTES(iWidth * 8) * (DWORD)iHeight;
    DWORD dwTotalSize = sizeof(BITMAPINFOHEADER) + dwMapSize + dwBitsSize;

    m_hDib = (HDIB)::GlobalAlloc(GMEM_MOVEABLE | GMEM_ZEROINIT, dwTotalSize);
    if(m_hDib==NULL)
        return E_OUTOFMEMORY;
    
    m_pDib = (PDIB)::GlobalLock((HGLOBAL)m_hDib);
	if(m_pDib==NULL)
		return E_FAIL;

    m_pDib->biSize = sizeof(BITMAPINFOHEADER);
    m_pDib->biWidth = (LONG)iWidth;
    m_pDib->biHeight = (LONG)iHeight;
    m_pDib->biPlanes = 1;
    m_pDib->biBitCount = 8;
    m_pDib->biCompression = BI_RGB;
    m_pDib->biSizeImage = dwBitsSize;
    m_pDib->biXPelsPerMeter = DEFAULT_PELS_METER;
    m_pDib->biYPelsPerMeter = DEFAULT_PELS_METER;
    m_pDib->biClrUsed = 0;
    m_pDib->biClrImportant = 0;

    m_iOffsetToBits = sizeof(BITMAPINFOHEADER) + dwMapSize;

    SetCMap(DIBCMAP_GRAYSCALE);
    return NOERROR;
}

HRESULT CDib::Create(BITMAPINFOHEADER *pbmi)
{
    Destroy();

    if(pbmi->biSize!=sizeof(BITMAPINFOHEADER))
        return E_INVALIDARG;

    if(pbmi->biBitCount==8)
    {
        if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        HRESULT hr = Create8Bit(pbmi->biWidth, abs(pbmi->biHeight));
        if(FAILED(hr))
            return hr;

        Byte *pClrMap = (Byte *)pbmi + pbmi->biSize;
        Byte *pBits;
        UINT uClrs = pbmi->biClrUsed;
        if(uClrs==0)
            uClrs = 256;
        pBits = pClrMap + sizeof(RGBQUAD) * uClrs;

        if(pbmi->biHeight>0)
            VtMemcpy(GetDibBits(), pBits, m_pDib->biSizeImage);
        else
        {
            int iY;
            Byte *p = GetDibBits() + (Height()-1) * (int)StrideBytes();
            for(iY=0; iY<Height(); iY++, p-=StrideBytes(), pBits+=StrideBytes())
                VtMemcpy(p, pBits, pbmi->biWidth);
        }
        VtMemcpy(GetCMapBits(), pClrMap, sizeof(RGBQUAD) * uClrs);
    }
    else if(pbmi->biBitCount==16)
    {
        bool bFiveBit = TRUE;
        if(pbmi->biCompression==BI_BITFIELDS) // may change this
        {
            DWORD *pClr = (DWORD *)((Byte *)pbmi + pbmi->biSize);
            if(pClr[0]==0xf800 && pClr[1]==0x07e0 && pClr[2]==0x001f)
                bFiveBit = FALSE;
            else if(pClr[0]!=0x7c00 || pClr[1]!=0x03e0 || pClr[2]!=0x001f)
                return E_NOTIMPL;
        }
        else if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        HRESULT hr = Create24Bit(pbmi->biWidth, abs(pbmi->biHeight));
        if(FAILED(hr))
            return hr;

        WORD *pBitsRow = (WORD *)((Byte *)pbmi + pbmi->biSize 
            + sizeof(RGBQUAD) * pbmi->biClrUsed);
        Byte *pDstBitsRow = GetDibBits();
        if(pbmi->biHeight<0)
            pDstBitsRow += StrideBytes() * (Height() - 1);

        int iX, iY;
        int iW = pbmi->biWidth;
        int iH = abs(pbmi->biHeight);
        int iDstInc = pbmi->biHeight<0 ? -(int)StrideBytes() : (int)StrideBytes();
        int iSrcInc = (VTWIDTHBYTES(iW * pbmi->biBitCount))>>1;

        if(bFiveBit)
            for(iY = 0; iY<iH; iY++, pDstBitsRow+=iDstInc, pBitsRow+=iSrcInc)
            {
                WORD *pBits = pBitsRow;
                Byte *pDstBits = pDstBitsRow;
                for(iX = 0; iX<iW; iX++, pBits++, pDstBits+=ElSize())
                {
                    pDstBits[0] = Byte((*pBits & 0x001f)<<3);  //b
                    pDstBits[1] = Byte((*pBits & 0x03e0)>>2);  //g
                    pDstBits[2] = Byte((*pBits & 0x7c00)>>7);  //r
                }
            }
        else
            for(iY = 0; iY<iH; iY++, pDstBitsRow+=iDstInc, pBitsRow+=iSrcInc)
            {
                WORD *pBits = pBitsRow;
                Byte *pDstBits = pDstBitsRow;
                for(iX = 0; iX<iW; iX++, pBits++, pDstBits+=ElSize())
                {
                    pDstBits[0] = Byte((*pBits & 0x001f)<<3);  //b
                    pDstBits[1] = Byte((*pBits & 0x07e0)>>3);  //g
                    pDstBits[2] = Byte((*pBits & 0xf800)>>8);  //r
                }
            }
    }
    else if(pbmi->biBitCount==24)
    {
        if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;
        HRESULT hr = Create24Bit(pbmi->biWidth, abs(pbmi->biHeight));
        if(FAILED(hr))
            return hr;

        Byte *pBits = (Byte *)pbmi + pbmi->biSize + sizeof(RGBQUAD) * pbmi->biClrUsed;
        if(pbmi->biHeight>0)
            VtMemcpy(GetDibBits(), pBits, m_pDib->biSizeImage);
        else
        {
            int iY;
            Byte *p = GetDibBits() + (Height()-1) * StrideBytes();
            for(iY=0; iY<Height(); iY++, p-=StrideBytes(), pBits+=StrideBytes())
                VtMemcpy(p, pBits, StrideBytes());
        }
    }
    else if(pbmi->biBitCount==32)
    {
        if(pbmi->biCompression==BI_BITFIELDS)
        {
            DWORD *pClr = (DWORD *)((Byte *)pbmi + pbmi->biSize);
            if(pClr[0]!=0x00ff0000 || pClr[1]!=0x0000ff00 || pClr[2]!=0x000000ff)
                return E_NOTIMPL;
        }
        else if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        HRESULT hr = Create32Bit(pbmi->biWidth, abs(pbmi->biHeight));
        if(FAILED(hr))
            return hr;

        Byte *pBits = (Byte *)pbmi + pbmi->biSize + sizeof(RGBQUAD) * pbmi->biClrUsed;
        if(pbmi->biHeight>0)
            VtMemcpy(GetDibBits(), pBits, m_pDib->biSizeImage);
        else
        {
            int iY;
            Byte *p = GetDibBits() + (Height()-1) * StrideBytes();
            for(iY=0; iY<Height(); iY++, p-=StrideBytes(), pBits+=StrideBytes())
                VtMemcpy(p, pBits, StrideBytes());
        }
    }
    else
        return E_NOTIMPL;

    if(pbmi->biXPelsPerMeter != 0)
    {
        m_pDib->biXPelsPerMeter = pbmi->biXPelsPerMeter;
        m_pDib->biYPelsPerMeter = pbmi->biYPelsPerMeter;
    }

    return NOERROR;
}

HRESULT CDib::CreateNoUpdate(BITMAPINFOHEADER *pbmi)
{
    Destroy();

    if(pbmi->biSize!=sizeof(BITMAPINFOHEADER))
        return E_INVALIDARG;

    if(pbmi->biBitCount==8)
    {
        if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        HRESULT hr = Create8Bit(pbmi->biWidth, abs(pbmi->biHeight));
        if(FAILED(hr))
            return hr;
    }
    else if(pbmi->biBitCount==16)
    {
        bool bFiveBit = TRUE;
        if(pbmi->biCompression==BI_BITFIELDS) // may change this
        {
            DWORD *pClr = (DWORD *)((Byte *)pbmi + pbmi->biSize);
            if(pClr[0]==0xf800 && pClr[1]==0x07e0 && pClr[2]==0x001f)
                bFiveBit = FALSE;
            else if(pClr[0]!=0x7c00 || pClr[1]!=0x03e0 || pClr[2]!=0x001f)
                return E_NOTIMPL;
        }
        else if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        HRESULT hr = Create24Bit(pbmi->biWidth, abs(pbmi->biHeight));
        if(FAILED(hr))
            return hr;
    }
    else if(pbmi->biBitCount==24)
    {
        if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;
        HRESULT hr = Create24Bit(pbmi->biWidth, abs(pbmi->biHeight));
        if(FAILED(hr))
            return hr;
    }
    else if(pbmi->biBitCount==32)
    {
        if(pbmi->biCompression==BI_BITFIELDS)
        {
            DWORD *pClr = (DWORD *)((Byte *)pbmi + pbmi->biSize);
            if(pClr[0]!=0x00ff0000 || pClr[1]!=0x0000ff00 || pClr[2]!=0x000000ff)
                return E_NOTIMPL;
        }
        else if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        HRESULT hr = Create32Bit(pbmi->biWidth, abs(pbmi->biHeight));
        if(FAILED(hr))
            return hr;
    }
    else
        return E_NOTIMPL;

    if(pbmi->biXPelsPerMeter != 0)
    {
        m_pDib->biXPelsPerMeter = pbmi->biXPelsPerMeter;
        m_pDib->biYPelsPerMeter = pbmi->biYPelsPerMeter;
    }

    return NOERROR;
}

HRESULT CDib::Update(BITMAPINFOHEADER *pbmi)
{
    if(pbmi->biSize!=sizeof(BITMAPINFOHEADER))
        return E_INVALIDARG;

    if(pbmi->biBitCount==8)
    {
        if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        if(ElSize()!=1 || Width()!=pbmi->biWidth || Height()!=abs(pbmi->biHeight))
            return E_FAIL;

        Byte *pClrMap = (Byte *)pbmi + pbmi->biSize;
        Byte *pBits;
        UINT uClrs = pbmi->biClrUsed;
        if(uClrs==0)
            uClrs = 256;
        pBits = pClrMap + sizeof(RGBQUAD) * uClrs;

        if(pbmi->biHeight>0)
            VtMemcpy(GetDibBits(), pBits, m_pDib->biSizeImage);
        else
        {
            int iY;
            Byte *p = GetDibBits() + (Height()-1) * (int)StrideBytes();
            for(iY=0; iY<Height(); iY++, p-=StrideBytes(), pBits+=StrideBytes())
                VtMemcpy(p, pBits, pbmi->biWidth);
        }
        VtMemcpy(GetCMapBits(), pClrMap, sizeof(RGBQUAD) * uClrs);
    }
    else if(pbmi->biBitCount==16)
    {
        bool bFiveBit = TRUE;
        if(pbmi->biCompression==BI_BITFIELDS) // may change this
        {
            DWORD *pClr = (DWORD *)((Byte *)pbmi + pbmi->biSize);
            if(pClr[0]==0xf800 && pClr[1]==0x07e0 && pClr[2]==0x001f)
                bFiveBit = FALSE;
            else if(pClr[0]!=0x7c00 || pClr[1]!=0x03e0 || pClr[2]!=0x001f)
                return E_NOTIMPL;
        }
        else if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        if(ElSize()!=3 || Width()!=pbmi->biWidth || Height()!=abs(pbmi->biHeight))
            return E_FAIL;

        WORD *pBitsRow = (WORD *)((Byte *)pbmi + pbmi->biSize 
            + sizeof(RGBQUAD) * pbmi->biClrUsed);
        Byte *pDstBitsRow = GetDibBits();
        if(pbmi->biHeight<0)
            pDstBitsRow += StrideBytes() * (Height() - 1);

        int iX, iY;
        int iW = pbmi->biWidth;
        int iH = abs(pbmi->biHeight);
        int iDstInc = pbmi->biHeight<0 ? -(int)StrideBytes() : (int)StrideBytes();
        int iSrcInc = (VTWIDTHBYTES(iW * pbmi->biBitCount))>>1;

        if(bFiveBit)
            for(iY = 0; iY<iH; iY++, pDstBitsRow+=iDstInc, pBitsRow+=iSrcInc)
            {
                WORD *pBits = pBitsRow;
                Byte *pDstBits = pDstBitsRow;
                for(iX = 0; iX<iW; iX++, pBits++, pDstBits+=ElSize())
                {
                    pDstBits[0] = Byte((*pBits & 0x001f)<<3);  //b
                    pDstBits[1] = Byte((*pBits & 0x03e0)>>2);  //g
                    pDstBits[2] = Byte((*pBits & 0x7c00)>>7);  //r
                }
            }
        else
            for(iY = 0; iY<iH; iY++, pDstBitsRow+=iDstInc, pBitsRow+=iSrcInc)
            {
                WORD *pBits = pBitsRow;
                Byte *pDstBits = pDstBitsRow;
                for(iX = 0; iX<iW; iX++, pBits++, pDstBits+=ElSize())
                {
                    pDstBits[0] = Byte((*pBits & 0x001f)<<3);  //b
                    pDstBits[1] = Byte((*pBits & 0x07e0)>>3);  //g
                    pDstBits[2] = Byte((*pBits & 0xf800)>>8);  //r
                }
            }
    }
    else if(pbmi->biBitCount==24)
    {
        if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        if(ElSize()!=3 || Width()!=pbmi->biWidth || Height()!=abs(pbmi->biHeight))
            return E_FAIL;

        Byte *pBits = (Byte *)pbmi + pbmi->biSize + sizeof(RGBQUAD) * pbmi->biClrUsed;
        if(pbmi->biHeight>0)
            VtMemcpy(GetDibBits(), pBits, m_pDib->biSizeImage);
        else
        {
            int iY;
            Byte *p = GetDibBits() + (Height()-1) * StrideBytes();
            for(iY=0; iY<Height(); iY++, p-=StrideBytes(), pBits+=StrideBytes())
                VtMemcpy(p, pBits, StrideBytes());
        }
    }
    else if(pbmi->biBitCount==32)
    {
        if(pbmi->biCompression==BI_BITFIELDS)
        {
            DWORD *pClr = (DWORD *)((Byte *)pbmi + pbmi->biSize);
            if(pClr[0]!=0x00ff0000 || pClr[1]!=0x0000ff00 || pClr[2]!=0x000000ff)
                return E_NOTIMPL;
        }
        else if(pbmi->biCompression!=BI_RGB)
            return E_NOTIMPL;

        if(ElSize()!=4 || Width()!=pbmi->biWidth || Height()!=abs(pbmi->biHeight))
            return E_FAIL;

        Byte *pBits = (Byte *)pbmi + pbmi->biSize + sizeof(RGBQUAD) * pbmi->biClrUsed;
        if(pbmi->biHeight>0)
            VtMemcpy(GetDibBits(), pBits, m_pDib->biSizeImage);
        else
        {
            int iY;
            Byte *p = GetDibBits() + (Height()-1) * StrideBytes();
            for(iY=0; iY<Height(); iY++, p-=StrideBytes(), pBits+=StrideBytes())
                VtMemcpy(p, pBits, StrideBytes());
        }
    }
    else
        return E_NOTIMPL;

    return NOERROR;
}

void CDib::Destroy()
{
    if(m_hDib!=NULL)
    {
        ::GlobalUnlock((HGLOBAL)m_hDib);
        ::GlobalFree((HGLOBAL)m_hDib);
        m_hDib = NULL;
    }

    m_pDib = NULL;
}

void CDib::SetCMap(int iType)
{
    if(m_hDib==NULL || m_pDib->biBitCount!=8)
        return;

    RGBQUAD *pCMap = GetCMapBits();
    int i;

    switch(iType)
    {
    case DIBCMAP_PSEUDO1:
        for(i = 0; i<256; i++)
        {
            VtClip(&pCMap[i].rgbRed, (i-85) * 3);
            VtClip(&pCMap[i].rgbGreen, (i-170) * 3);
            VtClip(&pCMap[i].rgbBlue, (i * 3));
        }

        break;

    case DIBCMAP_COLDHOT:
        for(i = 0; i<256; i++)
        {
            VtClip(&pCMap[i].rgbRed, (i-51) * 5);
            VtClip(&pCMap[i].rgbGreen, (i-153) * 5);
            if(i<100)
                VtClip(&pCMap[i].rgbBlue, (i * 5));
            else if(i>200)
                VtClip(&pCMap[i].rgbBlue, (i-204) * 5);
            else
                VtClip(&pCMap[i].rgbBlue, (153-i) * 5);
        }

        break;

    case DIBCMAP_COLORBARS:
        for(i = 0; i<256; i++)
        {
            RGBQUAD rgb;
            rgb.rgbRed = rgb.rgbGreen = rgb.rgbBlue = 0;
            if(i<=51)
            {
                rgb.rgbBlue = BYTE((i*255)/51);
            } else if(i<=102)
            {
                int q = (255*(i-51))/51;
                rgb.rgbRed = (BYTE)q;
                rgb.rgbBlue = BYTE(255-q);
            } else if(i<=153)
            {
                int q = (255*(i-102))/51;
                rgb.rgbRed = BYTE(255-q);
                rgb.rgbGreen = (BYTE)q;
            } else if(i<=204)
            {
                int q = (255*(i-153))/51;
                rgb.rgbRed = (BYTE)q;
                rgb.rgbGreen = 255;
            } else
            {
                int q = (255*(i-204))/51;
                rgb.rgbRed = 255;
                rgb.rgbGreen = 255;
                rgb.rgbBlue = (BYTE)q;
            }
            pCMap[i] = rgb;
        }
        break;

    case DIBCMAP_PSEUDO2:
        for(i=0; i<256; i++)
        {
            RGBQUAD rgb;
            rgb.rgbRed = rgb.rgbGreen = rgb.rgbBlue = 0;
            if(i<=43)
            {
                rgb.rgbBlue = 255;
                rgb.rgbRed = BYTE((i*255)/43);
            } else if(i<=86)
            {
                rgb.rgbRed = 255;
                rgb.rgbBlue = BYTE(255-((i-43)*255)/43);
            } else if(i<=128)
            {
                rgb.rgbRed = 255;
                rgb.rgbGreen = BYTE(((i-86)*255)/42);
            } else if(i<=171)
            {
                rgb.rgbGreen = 255;
                rgb.rgbRed = BYTE(255-((i-128)*255)/43);
            } else if(i<=214)
            {
                rgb.rgbGreen = 255;
                rgb.rgbBlue = BYTE(((i-171)*255)/43);
            } else
            {
                rgb.rgbBlue = 255;
                rgb.rgbGreen = BYTE(255-((i-214)*255)/42);
            }
            pCMap[i] = rgb;
        }

        break;

    case DIBCMAP_GRAYSCALE:
    default:
        for(i = 0; i<256; i++)
            pCMap[i].rgbRed = pCMap[i].rgbGreen = pCMap[i].rgbBlue = (Byte)i;

        break;
    }
}

void CDib::CopyCMap(CDib *pDib)
{
    if(pDib==NULL)
        return;
    if(pDib->ElSize()>1)
    {
        SetCMap(DIBCMAP_GRAYSCALE);
        return;
    }
    
    RGBQUAD *pCMapDst = GetCMapBits();
    RGBQUAD *pCMapSrc = pDib->GetCMapBits();
    int i;

    for(i = 0; i<256; i++)
        pCMapDst[i] = pCMapSrc[i];
}

void CDib::SetCMapColour(int iIndex, RGBAPix &rgb)
{
    if(m_hDib==NULL || m_pDib->biBitCount!=8 || iIndex<0 || iIndex>255)
        return;

    RGBQUAD *pCMap = GetCMapBits();
    pCMap[iIndex].rgbRed = rgb.r;
    pCMap[iIndex].rgbGreen = rgb.g;
    pCMap[iIndex].rgbBlue = rgb.b;
}

RGBAPix CDib::GetCMapColour(int iIndex) const
{
    RGBAPix rgb;

    rgb.r = rgb.g = rgb.b = 0;
    rgb.a = 255;

    if(m_hDib==NULL || m_pDib->biBitCount!=8)
        return rgb;
    if(iIndex>255)
        iIndex = 255;
    else if(iIndex<0)
        iIndex = 0;
    
    RGBQUAD *pCMap = GetCMapBits();
    rgb.r = pCMap[iIndex].rgbRed;
    rgb.g = pCMap[iIndex].rgbGreen;
    rgb.b = pCMap[iIndex].rgbBlue;

    return rgb;
}

HRESULT CDib::SetDibBits(Byte *puchPix, int iSrcRgnX, int iSrcRgnY,
    int iSrcRgnW, int iSrcRgnH, int iSrcH, int iSrcStride, int iDstX, int iDstY)
{
    if(m_hDib==NULL)
        return E_NOINIT;
    if(puchPix==NULL)
        return E_INVALIDARG;

    int iW = Width();
    int iH = Height();
    bool bInvertedSrc = FALSE;
    if(iSrcH<0)
    {
        bInvertedSrc = TRUE;
        iSrcH = -iSrcH;
    }
    if(iSrcH==0)
        iSrcH = iSrcRgnH;

    if(iSrcRgnH<=0 || iSrcRgnW<=0)
        return NOERROR;
    if(iDstX<0 || iDstX+iSrcRgnW>iW || iDstY<0 || iDstY+iSrcRgnH>iH)
        return NOERROR;

    size_t uPixSize = ElSize();
    if(iSrcStride<=0)
        iSrcStride = iSrcRgnW * (int)uPixSize;

    iDstY = iH - iDstY - 1; // invert Dst Y
    if(bInvertedSrc)
        iSrcRgnY = iSrcH - iSrcRgnY - 1; // invert Src Y

    int iDstStride = VTWIDTHBYTES(iW * m_pDib->biBitCount);
    size_t iRowCount = uPixSize * iSrcRgnW;
    int iY;

    Byte *puchSrc = puchPix + iSrcStride * iSrcRgnY + iSrcRgnX * (int)uPixSize;
    Byte *puchDst = (Byte *)m_pDib + m_iOffsetToBits 
        + iDstStride * iDstY + iDstX * (int)uPixSize;

    if(bInvertedSrc)
        iSrcStride = -iSrcStride;

    for(iY=0; iY<iSrcRgnH; iY++, puchSrc+=iSrcStride, puchDst-=iDstStride)
        VtMemcpy(puchDst, puchSrc, iRowCount);

    return NOERROR;
}

HRESULT CDib::GetPixels(CImg &cImg)
{
    if(m_hDib==NULL)
        return E_NOINIT;

    HRESULT hr = CreateImageForTransform(cImg, Width(), Height(), OBJ_RGBAIMG);
    if(FAILED(hr))
        return hr;

    int iMtxWidth  = cImg.Width();
    int iMtxHeight = cImg.Height();
    int iType      = IMG_FORMAT(cImg.GetType());

    if(EL_FORMAT(iType) != EL_FORMAT_BYTE && EL_FORMAT(iType) != EL_FORMAT_FLOAT)
        return E_NOTIMPL;

    if(cImg.Bands()!=1 && cImg.Bands()!=4)
        return E_INVALIDARG;

    int iWidth = Width();
    int iHeight = Height();
    size_t uSrcStride = StrideBytes();
    size_t uPixSize = ElSize();

    // get top source row (image is upside down)
    Byte *puchSrc = 
        (Byte *)m_pDib + m_iOffsetToBits + (iHeight - 1) * (int)uSrcStride;

    if(iMtxWidth<iWidth)
        iWidth = iMtxWidth;
    if(iMtxHeight<iHeight)
        iHeight = iMtxHeight;

    int iX, iY;

    for(iY=0; iY<iHeight; iY++, puchSrc-=uSrcStride)
    {
        Byte *puchSrcRow = puchSrc;
        
        if(EL_FORMAT(iType)== EL_FORMAT_BYTE && cImg.Bands() == 4)
        {
            RGBAPix *pRGB = ((CRGBAImg *)&cImg)->Ptr(iY);
                    
            if(m_pDib->biBitCount==32)
            {
                // rgba -> rgba direct conversion
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pRGB++)
                {
                    pRGB->b = puchSrcRow[0];
                    pRGB->g = puchSrcRow[1];
                    pRGB->r = puchSrcRow[2];
                    pRGB->a = puchSrcRow[3];
                }
            }
            else if(m_pDib->biBitCount==24)
            {
                // rgb -> rgba conversion 
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pRGB++)
                {
                    pRGB->b = puchSrcRow[0];
                    pRGB->g = puchSrcRow[1];
                    pRGB->r = puchSrcRow[2];
                    pRGB->a = 255;
                }
            }
            else
            {
                // index -> rgb conversion via colourmap
                RGBQUAD *pClr, *pCMap = GetCMapBits();
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pRGB++)
                {
                    pClr = pCMap + puchSrcRow[0];
                    pRGB->r = pClr->rgbRed;
                    pRGB->g = pClr->rgbGreen;
                    pRGB->b = pClr->rgbBlue;
                    pRGB->a = 255;
                }
            }
        }
        else if(EL_FORMAT(iType) == EL_FORMAT_FLOAT && cImg.Bands() == 4)
        {
            RGBAFloatPix *pRGB = ((CRGBAFloatImg *)&cImg)->Ptr(iY);

            if(m_pDib->biBitCount==32)
            {
                // rgba -> rgba direct conversion 
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pRGB++)
                {
                    pRGB->b = puchSrcRow[0];
                    pRGB->g = puchSrcRow[1];
                    pRGB->r = puchSrcRow[2];
                    pRGB->a = puchSrcRow[3];
                }
            }
            else if(m_pDib->biBitCount==24)
            {
                // rgb -> rgba conversion 
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pRGB++)
                {
                    pRGB->b = puchSrcRow[0];
                    pRGB->g = puchSrcRow[1];
                    pRGB->r = puchSrcRow[2];
                    pRGB->a = 255;
                }
            }
            else
            {
                // index -> rgb conversion via colourmap
                RGBQUAD *pClr, *pCMap = GetCMapBits();
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pRGB++)
                {
                    pClr = pCMap + puchSrcRow[0];
                    pRGB->r = pClr->rgbRed;
                    pRGB->g = pClr->rgbGreen;
                    pRGB->b = pClr->rgbBlue;
                    pRGB->a = 255;
                }
            }
        }
        else if(EL_FORMAT(iType) == EL_FORMAT_BYTE && cImg.Bands() == 1)
        {
            Byte *pDst = ((CByteImg *)&cImg)->Ptr(iY);

            if(m_pDib->biBitCount==24 || m_pDib->biBitCount==32)
            {
                // rgb -> byte conversion by calculating luminance
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pDst++)
                {
                    *pDst = (Byte)VtLumaFromRGB_CCIR601YPbPr(
                        puchSrcRow[2], puchSrcRow[1], puchSrcRow[0]);
                }
            }
            else
            {
                // index -> byte conversion via colourmap then calculating luminance
                RGBQUAD *pClr, *pCMap = GetCMapBits();
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pDst++)
                {
                    pClr = pCMap + puchSrcRow[0];
                    *pDst = (Byte)VtLumaFromRGB_CCIR601YPbPr(
                        pClr->rgbRed, pClr->rgbGreen, pClr->rgbBlue);
                }
            }
        }
        else // EL_FORMAT(iType) == EL_FORMAT_FLOAT && cImg.Bands() == 1
        {
            float *pDst = ((CFloatImg *)&cImg)->Ptr(iY);

            if(m_pDib->biBitCount==24 || m_pDib->biBitCount==32)
            {
                // rgb -> float conversion by calculating luminance
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pDst++)
                {
                    *pDst = VtLumaFromRGB_CCIR601YPbPr(
                        puchSrcRow[2], puchSrcRow[1], puchSrcRow[0]);
                }
            }
            else
            {
                // index -> byte conversion via colourmap then calculating luminance
                RGBQUAD *pClr, *pCMap = GetCMapBits();
                for(iX=0; iX<iWidth; iX++, puchSrcRow+=uPixSize, pDst++)
                {
                    pClr = pCMap + puchSrcRow[0];
                    *pDst = VtLumaFromRGB_CCIR601YPbPr(
                        pClr->rgbRed, pClr->rgbGreen, pClr->rgbBlue);
                }
            }
        }
    }

    return NOERROR;
}

HRESULT CDib::SetPixels(const CImg &cImg, int iDstX, int iDstY)
{
    if(m_hDib==NULL)
        return E_NOINIT;

    int iMtxWidth, iMtxHeight, iType;
    iMtxWidth = cImg.Width();
    iMtxHeight = cImg.Height();
    iType = cImg.GetType();

    if (EL_FORMAT(iType) != EL_FORMAT_BYTE || cImg.Bands() != 4)
    {
        return E_NOTIMPL;
    }

    int iWidth = Width();
    int iHeight = Height();
    size_t uDstStride = StrideBytes();
    size_t uPixSize = ElSize();

    if(iDstX+iMtxWidth<=0 || iDstX>=iWidth || iDstY+iMtxHeight<=0 || iDstY>=iHeight)
        return NOERROR;

    // clip regions
    int iLX = VtMax(iDstX, 0);
    int iLY = VtMax(iDstY, 0);
    int iRX = VtMin(iDstX + iMtxWidth, iWidth);
    int iRY = VtMin(iDstY + iMtxHeight, iHeight);
    iWidth = iRX - iLX;
    iHeight = iRY - iLY;
    int iSX = iLX - iDstX;
    int iSY = iLY - iDstY;

    iLY = (Height() - 1) - iLY; // invert Y
    
    Byte *puchDst = (Byte *)m_pDib + m_iOffsetToBits 
        + (int)uDstStride * iLY + iLX * (int)uPixSize;

    int iX, iY;
   
    for(iY=0; iY<iHeight; iY++, puchDst-=uDstStride)
    {
        Byte *puchDstRow = puchDst;

        if(cImg.Bands() == 4)
        {
            RGBAPix *pRGB = ((CRGBAImg *)&cImg)->Ptr(iSX, iY + iSY);

            if(m_pDib->biBitCount==32)
            {
                // rgb -> rgb direct conversion
                for(iX=0; iX<iWidth; iX++, puchDstRow+=uPixSize, pRGB++)
                {
                    puchDstRow[0] = pRGB->b;
                    puchDstRow[1] = pRGB->g;
                    puchDstRow[2] = pRGB->r;
                    puchDstRow[3] = pRGB->a;
                }
            }
            else if(m_pDib->biBitCount==24)
            {
                // rgb -> rgb direct conversion
                for(iX=0; iX<iWidth; iX++, puchDstRow+=uPixSize, pRGB++)
                {
                    puchDstRow[0] = pRGB->b;
                    puchDstRow[1] = pRGB->g;
                    puchDstRow[2] = pRGB->r;
                }
            }
            else
            {
                // rgb -> index convert to luma, assumes grayscale cmap
                for(iX=0; iX<iWidth; iX++, puchDstRow+=uPixSize, pRGB++)
                {
                    puchDstRow[0] = (Byte)VtLumaFromRGB_CCIR601YPbPr(pRGB);
                }
            }
        }
        else // Bands() == 1
        {
            Byte *pSrc = ((CByteImg *)&cImg)->Ptr(iSX, iY + iSY);
            if(m_pDib->biBitCount==32)
            {
                // byte -> rgb direct conversion
                for(iX=0; iX<iWidth; iX++, puchDstRow+=uPixSize, pSrc++)
                {
                    puchDstRow[0] = puchDstRow[1] = puchDstRow[2] = *pSrc;
                    puchDstRow[3] = 255;
                }
            }
            else if(m_pDib->biBitCount==24)
            {
                // byte -> rgb direct conversion
                for(iX=0; iX<iWidth; iX++, puchDstRow+=uPixSize, pSrc++)
                    puchDstRow[0] = puchDstRow[1] = puchDstRow[2] = *pSrc;
            }
            else
            {
                // byte -> index direct conversion, assumes grayscale cmap
                VtMemcpy(puchDst, pSrc, iWidth);
            }
        }
    }

    return NOERROR;
}

void CDib::GetPixel(int iX, int iY, RGBAPix &cRGB) const
{
    Byte *puchBits = GetDibBits() + iX * (int)ElSize()
        + ((Height() - 1) - iY) * (int)StrideBytes();

    if(m_pDib->biBitCount==8)
    {
        RGBQUAD *pClr = GetCMapBits() + puchBits[0];
        cRGB.r = pClr->rgbRed;
        cRGB.g = pClr->rgbGreen;
        cRGB.b = pClr->rgbBlue;
        cRGB.a = 255;
    }
    else if(m_pDib->biBitCount==24)
    {
        cRGB.r = puchBits[2];
        cRGB.g = puchBits[1];
        cRGB.b = puchBits[0];
        cRGB.a = 255;
    }
    else //32
    {
        cRGB.r = puchBits[2];
        cRGB.g = puchBits[1];
        cRGB.b = puchBits[0];
        cRGB.a = puchBits[3];
    }
}

void CDib::SetPixel(int iX, int iY, RGBAPix &cRGB)
{
    Byte *puchBits = GetDibBits() + iX * (int)ElSize()
        + ((Height() - 1) - iY) * (int)StrideBytes();

    if(m_pDib->biBitCount==8)
    {
        puchBits[0] = (Byte)VtLumaFromRGB_CCIR601YPbPr(&cRGB);
    }
    else if(m_pDib->biBitCount==24)
    {
        puchBits[2] = cRGB.r;
        puchBits[1] = cRGB.g;
        puchBits[0] = cRGB.b;
    }
    else //32
    {
        puchBits[2] = cRGB.r;
        puchBits[1] = cRGB.g;
        puchBits[0] = cRGB.b;
        puchBits[3] = cRGB.a;
    }
}

HRESULT CDib::Load(FILE *fp)
{
    if(fp==NULL)
        return E_INVALIDARG;

    BITMAPFILEHEADER bmfh;

    if(fread(&bmfh, sizeof(BITMAPFILEHEADER), 1, fp) != 1)
        return E_READFAILED;
    if(bmfh.bfType != DIB_HEADER_MARKER)
        return E_BADFORMAT;

    BITMAPINFOHEADER bmih;

    if(fread(&bmih, sizeof(BITMAPINFOHEADER), 1, fp) != 1)
        return E_READFAILED;

    bool bCoreHeader = FALSE;

    if(bmih.biSize != sizeof(BITMAPINFOHEADER))
    {
        // something else
        if(bmih.biSize == sizeof(BITMAPCOREHEADER))
        {
            bCoreHeader = TRUE;

            BITMAPCOREHEADER *pbmch = (BITMAPCOREHEADER *)&bmih;
            WORD wWidth = pbmch->bcWidth;
            WORD wHeight = pbmch->bcHeight;
            WORD wBitCount = pbmch->bcBitCount;
            WORD wPlanes = pbmch->bcPlanes;

            bmih.biSize = sizeof(BITMAPINFOHEADER);
            bmih.biWidth = (LONG) wWidth;
            bmih.biHeight = (LONG) wHeight;
            bmih.biPlanes = wPlanes;
            bmih.biBitCount = wBitCount;
            bmih.biCompression = BI_RGB;
            bmih.biSizeImage = 0;
            bmih.biXPelsPerMeter = 0;
            bmih.biYPelsPerMeter = 0;
            bmih.biClrUsed = 0;
            bmih.biClrImportant = 0;

            if(fseek(fp, (long)(long(sizeof(BITMAPCOREHEADER)) - 
				long(sizeof(BITMAPINFOHEADER))), SEEK_CUR))
                return E_READFAILED;
        }
        else
            return E_BADFORMAT;
    }

    if(bmih.biPlanes != 1)
        return E_BADFORMAT;
    if(bmih.biCompression != BI_RGB)
        return E_BADFORMAT;
    if(bmih.biBitCount == 24 || bmih.biBitCount == 32)
    {
        if(bmih.biClrUsed > 0)
        {
            // discard colour table
            if(!bCoreHeader)
            {
                if(fseek(fp, sizeof(RGBQUAD) * bmih.biClrUsed, SEEK_CUR))
                    return E_READFAILED;
            }
            else
            {
                if(fseek(fp, sizeof(RGBTRIPLE) * bmih.biClrUsed, SEEK_CUR))
                    return E_READFAILED;
            }
        }

        if(bmih.biBitCount == 24)
            Create24Bit(bmih.biWidth, bmih.biHeight);
        else
            Create32Bit(bmih.biWidth, bmih.biHeight);

        if(bmih.biXPelsPerMeter != 0)
        {
            m_pDib->biXPelsPerMeter = bmih.biXPelsPerMeter;
            m_pDib->biYPelsPerMeter = bmih.biYPelsPerMeter;
        }

        if(fread(GetDibBits(), m_pDib->biSizeImage, 1, fp) != 1)
            return E_READFAILED;

        return NOERROR;
    }
    else if(bmih.biBitCount == 8)
    {
        Create8Bit(bmih.biWidth, bmih.biHeight);

        if(bmih.biXPelsPerMeter != 0)
        {
            m_pDib->biXPelsPerMeter = bmih.biXPelsPerMeter;
            m_pDib->biYPelsPerMeter = bmih.biYPelsPerMeter;
        }

        int iClrs = 256;
        if(bmih.biClrUsed > 0)
            iClrs = bmih.biClrUsed;

        if(iClrs>256)
            return E_BADFORMAT;

        if(!bCoreHeader)
        {
            if(fread(GetCMapBits(), sizeof(RGBQUAD), iClrs, fp) != (size_t)iClrs)
                return E_READFAILED;
        }
        else
        {
            RGBTRIPLE *pRGB = VT_NOTHROWNEW RGBTRIPLE [256];
            if(pRGB==NULL)
                return E_OUTOFMEMORY;
            if(fread(pRGB, sizeof(RGBTRIPLE), iClrs, fp) != (size_t)iClrs)
            {
                delete [] pRGB;
                return E_READFAILED;
            }

            int i;
            RGBQUAD *pRGBQ = GetCMapBits();
            for(i = 0; i<iClrs; i++)
            {
                pRGBQ->rgbRed = pRGB->rgbtRed;
                pRGBQ->rgbGreen = pRGB->rgbtGreen;
                pRGBQ->rgbBlue = pRGB->rgbtBlue;
            }
            delete [] pRGB;
        }

        if(fread(GetDibBits(), m_pDib->biSizeImage, 1, fp) != 1)
            HRESULT E_READFAILED;

        return NOERROR;
    }
    else
        return E_BADFORMAT;
}

HRESULT CDib::Save(FILE *fp) const
{
    if(fp==NULL)
        return E_INVALIDARG;

    if(m_hDib==NULL)
        return E_NOINIT;
    
    BITMAPFILEHEADER bmfh;

    bmfh.bfType = DIB_HEADER_MARKER;
    bmfh.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)
        + m_pDib->biSizeImage;
    bmfh.bfReserved1 = 0;
    bmfh.bfReserved2 = 0;
    bmfh.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);

    if(m_pDib->biBitCount == 8)
    {
        bmfh.bfSize += 256 * sizeof(RGBQUAD);
        bmfh.bfOffBits += 256 * sizeof(RGBQUAD);
    }

    if(fwrite(&bmfh, sizeof(BITMAPFILEHEADER), 1, fp) != 1)
        return E_WRITEFAILED;

    if(fwrite(m_pDib, sizeof(BITMAPINFOHEADER), 1, fp) != 1)
        return E_WRITEFAILED;

    if(m_pDib->biBitCount == 8)
        if(fwrite(GetCMapBits(), sizeof(RGBQUAD), 256, fp) != 256)
            return E_WRITEFAILED;

    if(fwrite(GetDibBits(), m_pDib->biSizeImage, 1, fp) != 1)
        return E_WRITEFAILED;

    return NOERROR;
}

////////////////////////////////////////////////////////////////////////
// Library Routines
////////////////////////////////////////////////////////////////////////

/*
HRESULT ImgLoadDib(LPCTSTR pchName, CxMatrix **ppMtxRtn, int iType)
{
    if(pchName==NULL || ppMtxRtn==NULL)
        return ERES_PARAMS;
    if(iType!=MTX_RGB && iType!=MTX_BYTE && iType!=MTX_FLOAT)
        return ERES_UNSUP;
    
    FILE *fp = _tfopen(pchName, _T("rb"));
    if(fp==NULL)
        return ERES_NOFILE;

    CDib dib;
    HRESULT ir = dib.Load(fp);

    fclose(fp);

    if(ir.IsError())
        return ir;

    CxMatrix *pMtx = ImgCreateMatrix(dib.Width(), dib.Height(), iType);

    if(pMtx==NULL)
        return ERES_MEMORY;

    ir = dib.GetPixels(pMtx);
    if(ir.IsError())
    {
        delete pMtx;
        return ir;
    }

    *ppMtxRtn = pMtx;
    return ERES_OK;
}

HRESULT ImgSaveDib(LPCTSTR pchName, CxMatrix *pMtx, int iCMapType)
{
    if(pchName==NULL || pMtx==NULL)
        return ERES_PARAMS;

    int iType = pMtx->GetType();
    if(iType != MTX_RGB && iType != MTX_BYTE)
        return ERES_UNSUP;

    FILE *fp = _tfopen(pchName, _T("wb"));
    if(fp==NULL)
        return ERES_NOFILE;

    int iWidth, iHeight;
    pMtx->GetSize(&iWidth, &iHeight);

    CDib dib;
    HRESULT ir;

    if(iType == MTX_RGB)
        ir = dib.Create24Bit(iWidth, iHeight);
    else
    {
        ir = dib.Create8Bit(iWidth, iHeight);
        if(ir.IsOK())
            dib.SetCMap(iCMapType);
    }

    if(ir.IsOK())
        ir = dib.SetPixels(pMtx, 0, 0, iWidth, iHeight, 0, 0);
    if(ir.IsOK())
        ir = dib.Save(fp);

    fclose(fp);

    return ir;
}
*/

#endif
