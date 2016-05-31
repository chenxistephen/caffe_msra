//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Classes for handling video frames
//
//  History:
//      2011/07/05-sbaker
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_vimage.h"
#include "vt_convert.h"

using namespace vt;


//+-----------------------------------------------------------------------
//
//  CRGB32VideoImg
//
//------------------------------------------------------------------------

HRESULT CRGB32VideoImg::Create(int iWidth, int iHeight)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((iWidth%2 != 0 || iHeight%2 != 0) ? E_INVALIDARG : S_OK);

    VT_HR_EXIT(m_imgRGB.Create(iWidth, iHeight));
    m_imgVidInfo.rectValidPixels = CRect(0, 0, iWidth, iHeight);
    m_imgVidInfo.eValidPixelsRectType = CVideoImgInfo::MinDisplay;
    m_imgVidInfo.dAspectRatio = 1.0;
    m_imgVidInfo.iInterlaceMode = 2;                         // Progressive

    VT_HR_EXIT_LABEL()
    if (FAILED(hr))
    {
        Deallocate();
    }
    return hr;
}

HRESULT CRGB32VideoImg::Create(Byte *pbBuffer, int iWidth, int iHeight, int iStrideBytes)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((iWidth%2 != 0 || iHeight%2 != 0) ? E_INVALIDARG : S_OK);

    VT_HR_EXIT(m_imgRGB.Create(pbBuffer, iWidth, iHeight, iStrideBytes));
    m_imgVidInfo.rectValidPixels = CRect(0, 0, iWidth, iHeight);
    m_imgVidInfo.eValidPixelsRectType = CVideoImgInfo::MinDisplay;
    m_imgVidInfo.dAspectRatio = 1.0;
    m_imgVidInfo.iInterlaceMode = 2;                         // Progressive

    VT_HR_EXIT_LABEL()
    if (FAILED(hr))
    {
        Deallocate();
    }
    return hr;
}

HRESULT CRGB32VideoImg::Create(int iWidth, int iHeight, const CVideoImgInfo& vi)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((iWidth%2 != 0 || iHeight%2 != 0) ? E_INVALIDARG : S_OK);

    VT_HR_EXIT(m_imgRGB.Create(iWidth, iHeight));
    m_imgVidInfo = vi;

    VT_HR_EXIT_LABEL()
    if (FAILED(hr))
    {
        Deallocate();
    }
    return hr;
}

HRESULT CRGB32VideoImg::Create(Byte *pbBuffer, int iWidth, int iHeight, int iStrideBytes, const CVideoImgInfo& vi)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((iWidth%2 != 0 || iHeight%2 != 0) ? E_INVALIDARG : S_OK);

    VT_HR_EXIT(m_imgRGB.Create(pbBuffer, iWidth, iHeight, iStrideBytes));
    m_imgVidInfo = vi;

    VT_HR_EXIT_LABEL()
    if (FAILED(hr))
    {
        Deallocate();
    }
    return hr;
}


//+-----------------------------------------------------------------------
//
//  CNV12VImg
//
//------------------------------------------------------------------------
	
// BUG BUG XXX: why is this code repeated 4X - need to share more of it, 
//
HRESULT CNV12VideoImg::Create(int iWidth, int iHeight)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((iWidth%2 != 0 || iHeight%2 != 0) ? E_INVALIDARG : S_OK);

    // allocate the data via m_imgData - make sure to not use the stride
    // from m_imgData since it may not be the same as the width
    VT_HR_EXIT( m_imgData.Create(iWidth, iHeight+(iHeight/2)) );

    // wrap the imgY and imgUV around the imgData allocated memory
    int iUVWidth = iWidth/2;
    int iUVHeight = iHeight/2;
	vt::CRect rctY(0,0,iWidth,iHeight);
	m_imgData.Share(m_imgY, &rctY);
	// BUG BUG XXX:
	//     m_imgUV needs to use Share() instead of just grabbing the pointer
	//     this is currently not possible due to the conflicting image types
	VT_HR_EXIT(m_imgUV.Create(m_imgData.BytePtr()+(iHeight*m_imgData.StrideBytes()),
		                      iUVWidth, iUVHeight, iWidth));

    m_imgVidInfo.rectValidPixels = CRect(0, 0, iWidth, iHeight);
    m_imgVidInfo.eValidPixelsRectType = CVideoImgInfo::MinDisplay;
    m_imgVidInfo.dAspectRatio = 1.0;
    m_imgVidInfo.iInterlaceMode = 2;                         // Progressive

    VT_HR_EXIT_LABEL()
    if (FAILED(hr))
    {
        Deallocate();
    }
    return hr;
}

HRESULT CNV12VideoImg::Create(Byte *pbBuffer, int iWidth, int iHeight, int iStrideBytes)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((iWidth%2 != 0 || iHeight%2 != 0) ? E_INVALIDARG : S_OK);

    // wrap the m_imgData around the buffer
    VT_HR_EXIT( m_imgData.Create(pbBuffer, iWidth, iHeight+(iHeight/2), iStrideBytes) );

    int iUVWidth = iWidth/2;
    int iUVHeight = iHeight/2;
	vt::CRect rctY(0,0,iWidth,iHeight);
	m_imgData.Share(m_imgY, &rctY);
	// BUG BUG XXX:
	//     m_imgUV needs to use Share() instead of just grabbing the pointer
	//     this is currently not possible due to the conflicting image types
    VT_HR_EXIT(m_imgUV.Create(pbBuffer+iHeight*iStrideBytes, iUVWidth, iUVHeight, iStrideBytes));
    m_imgVidInfo.rectValidPixels = CRect(0, 0, iWidth, iHeight);
    m_imgVidInfo.eValidPixelsRectType = CVideoImgInfo::MinDisplay;
    m_imgVidInfo.dAspectRatio = 1.0;
    m_imgVidInfo.iInterlaceMode = 2;                         // Progressive

    VT_HR_EXIT_LABEL()
    if (FAILED(hr))
    {
        Deallocate();
    }
    return hr;
}

HRESULT CNV12VideoImg::Create(int iWidth, int iHeight, const CVideoImgInfo& vi)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((iWidth%2 != 0 || iHeight%2 != 0) ? E_INVALIDARG : S_OK);

    // allocate the data with m_imgData - make sure to not use the stride
    // from m_imgData since it may not be the same as the width
    VT_HR_EXIT( m_imgData.Create(iWidth, iHeight+(iHeight/2)) );

    // wrap the imgY and imgUV around the imgData allocated memory
    int iUVWidth = iWidth/2;
    int iUVHeight = iHeight/2;
	vt::CRect rctY(0,0,iWidth,iHeight);
	m_imgData.Share(m_imgY, &rctY);
	// BUG BUG XXX:
	//     m_imgUV needs to use Share() instead of just grabbing the pointer
	//     this is currently not possible due to the conflicting image types
	VT_HR_EXIT(m_imgUV.Create(m_imgData.BytePtr()+(iHeight*m_imgData.StrideBytes()), 
		                      iUVWidth, iUVHeight, iWidth));
    m_imgVidInfo = vi;

    VT_HR_EXIT_LABEL()
    if (FAILED(hr))
    {
        Deallocate();
    }
    return hr;
}

HRESULT CNV12VideoImg::Create(Byte *pbBuffer, int iWidth, int iHeight, int iStrideBytes, const CVideoImgInfo& vi)
{
    VT_HR_BEGIN();

    VT_HR_EXIT((iWidth%2 != 0 || iHeight%2 != 0) ? E_INVALIDARG : S_OK);

    // wrap the m_imgData around the buffer
    VT_HR_EXIT( m_imgData.Create(pbBuffer, iWidth, iHeight+(iHeight/2), iStrideBytes) );

    int iUVWidth = iWidth/2;
    int iUVHeight = iHeight/2;
	vt::CRect rctY(0,0,iWidth,iHeight);
	m_imgData.Share(m_imgY, &rctY);
	// BUG BUG XXX:
	//     m_imgUV needs to use Share() instead of just grabbing the pointer
	//     this is currently not possible due to the conflicting image types
    VT_HR_EXIT(m_imgUV.Create(pbBuffer+iHeight*iStrideBytes, iUVWidth, iUVHeight, iStrideBytes));
    m_imgVidInfo = vi;

    VT_HR_EXIT_LABEL()
    if (FAILED(hr))
    {
        Deallocate();
    }
    return hr;
}
//+-----------------------------------------------------------------------
//
//  VtConvertImage: NV12 <-> RGB32
//
//------------------------------------------------------------------------

HRESULT vt::VtConvertVideoImage(CNV12VideoImg &imgDst, const CRGB32VideoImg &imgSrc)
{
    VT_HR_BEGIN();

    int iWidth  = imgSrc.GetImg().Width();
    int iHeight = imgSrc.GetImg().Height();
    VT_HR_EXIT((iWidth<=0 || iHeight<=0) ? E_INVALIDARG : S_OK);
    VT_HR_EXIT(imgDst.Create(iWidth, iHeight,imgSrc.GetVideoImgInfo()));

    for(int iY=0; iY<iHeight; iY+=2)
    {
        const RGBAPix* pRGBPixels1    = imgSrc.GetImg().Ptr(iY);
        const RGBAPix* pRGBPixels1End = pRGBPixels1 + iWidth;
        const RGBAPix* pRGBPixels2    = imgSrc.GetImg().Ptr(iY+1);
        Byte* pucYPixels1 = imgDst.GetYImg().Ptr(iY);
        Byte* pucYPixels2 = imgDst.GetYImg().Ptr(iY+1);
        UVPix* pUVPixels  = imgDst.GetUVImg().Ptr(iY/2);
        while(pRGBPixels1 < pRGBPixels1End)
        {
            // Read in RGB
            Byte iB11 = pRGBPixels1->b;
            Byte iG11 = pRGBPixels1->g;
            Byte iR11 = pRGBPixels1->r;
            pRGBPixels1++;
            Byte iB12 = pRGBPixels1->b;
            Byte iG12 = pRGBPixels1->g;
            Byte iR12 = pRGBPixels1->r;
            pRGBPixels1++;
            Byte iB21 = pRGBPixels2->b;
            Byte iG21 = pRGBPixels2->g;
            Byte iR21 = pRGBPixels2->r;
            pRGBPixels2++;
            Byte iB22 = pRGBPixels2->b;
            Byte iG22 = pRGBPixels2->g;
            Byte iR22 = pRGBPixels2->r;
            pRGBPixels2++;

            // Y Channel -> Each Pixel Separate Estimate
            *pucYPixels1++ = VtLumaFromRGB_CCIR601YCbCr(iR11, iG11, iB11);
            *pucYPixels1++ = VtLumaFromRGB_CCIR601YCbCr(iR12, iG12, iB12);
            *pucYPixels2++ = VtLumaFromRGB_CCIR601YCbCr(iR21, iG21, iB21);
            *pucYPixels2++ = VtLumaFromRGB_CCIR601YCbCr(iR22, iG22, iB22);

            // UV Channels -> Average to subsample
            Byte iAveB = Byte((int(iB11)+int(iB12)+int(iB21)+int(iB22)+2) >> 2);
            Byte iAveG = Byte((int(iG11)+int(iG12)+int(iG21)+int(iG22)+2) >> 2);
            Byte iAveR = Byte((int(iR11)+int(iR12)+int(iR21)+int(iR22)+2) >> 2);
            pUVPixels->u = VtCbFromRGB_CCIR601YCbCr(iAveR, iAveG, iAveB);
            pUVPixels->v = VtCrFromRGB_CCIR601YCbCr(iAveR, iAveG, iAveB);
            pUVPixels++;
        }
    }

    VT_HR_END();
}


// Core Conversion Functions
unsigned char RFromCE(int iC, int iE)
{
    int iR = (298*iC + 409*iE + 128) >> 8;
    if (iR <= 0) 
    {    
        return 0;
    }
    else if (iR >= 255) 
    {
        return 255;
    }
    else
    {
        return (unsigned char) iR;
    }
}
unsigned char GFromCDE(int iC, int iD, int iE)
{
    int iG = (298*iC - 100*iD - 208*iE + 128) >> 8;
    if (iG <= 0) 
    {    
        return 0;
    }
    else if (iG >= 255) 
    {
        return 255;
    }
    else
    {
        return (unsigned char) iG;
    }
}
unsigned char BFromCD(int iC, int iD)
{
    int iB = (298*iC + 516*iD + 128) >> 8;
    if (iB <= 0) 
    {    
        return 0;
    }
    else if (iB >= 255) 
    {
        return 255;
    }
    else
    {
        return (unsigned char) iB;
    }
}

HRESULT vt::VtConvertVideoImage(CRGB32VideoImg &imgDst, const  CNV12VideoImg &imgSrc)
{
    VT_HR_BEGIN();

    int iWidth  = imgSrc.GetYImg().Width();
    int iHeight = imgSrc.GetYImg().Height();
    VT_HR_EXIT((iWidth<=0 || iHeight<=0) ? E_INVALIDARG : S_OK);
    VT_HR_EXIT(imgDst.Create(iWidth, iHeight, imgSrc.GetVideoImgInfo()));

    for(int iY=0; iY<iHeight; iY+=2)
    {
        unsigned char* pucRGBPixels1 = (unsigned char*) imgDst.GetImg().Ptr(iY);
        unsigned char* pucRGBPixels1End = pucRGBPixels1 + 4*iWidth;
        unsigned char* pucRGBPixels2 = (unsigned char*) imgDst.GetImg().Ptr(iY+1);
        unsigned char* pucYPixels1 = (unsigned char*) imgSrc.GetYImg().Ptr(iY);
        unsigned char* pucYPixels2 = (unsigned char*) imgSrc.GetYImg().Ptr(iY+1);
        unsigned char* pucUVPixels = (unsigned char*) imgSrc.GetUVImg().Ptr(iY/2);
        while(pucRGBPixels1 < pucRGBPixels1End)
        {
            // Read in YUV data
            int iY11 = int(*pucYPixels1++);
            int iY12 = int(*pucYPixels1++);
            int iY21 = int(*pucYPixels2++);
            int iY22 = int(*pucYPixels2++);
            int iU = int(*pucUVPixels++);
            int iV = int(*pucUVPixels++);

            // Shift
            int iC11 = iY11 - 16;
            int iC12 = iY12 - 16;
            int iC21 = iY21 - 16;
            int iC22 = iY22 - 16;
            int iD = iU - 128;
            int iE = iV - 128;

            // Compute RGB
            *pucRGBPixels1++ = BFromCD(iC11, iD);
            *pucRGBPixels1++ = GFromCDE(iC11, iD, iE);
            *pucRGBPixels1++ = RFromCE(iC11, iE);
            *pucRGBPixels1++ = 255;
            *pucRGBPixels1++ = BFromCD(iC12, iD);
            *pucRGBPixels1++ = GFromCDE(iC12, iD, iE);
            *pucRGBPixels1++ = RFromCE(iC12, iE);
            *pucRGBPixels1++ = 255;
            *pucRGBPixels2++ = BFromCD(iC21, iD);
            *pucRGBPixels2++ = GFromCDE(iC21, iD, iE);
            *pucRGBPixels2++ = RFromCE(iC21, iE);
            *pucRGBPixels2++ = 255;
            *pucRGBPixels2++ = BFromCD(iC22, iD);
            *pucRGBPixels2++ = GFromCDE(iC22, iD, iE);
            *pucRGBPixels2++ = RFromCE(iC22, iE);
            *pucRGBPixels2++ = 255;
        }    
    }

    VT_HR_END();
}

//+----------------------------------------------------------------------------
//
// Fast NV12->RGBA Conversion
//
//-----------------------------------------------------------------------------

#if defined(_M_ARM)
extern "C" void NeonYUVToRGBALine(int width, const Byte* pYRow, const Byte* pUVRow, Byte* pRgbRow);
extern "C" void NeonYUVToRGBALineFast(int width, const Byte* pYRow, const Byte* pUVRow, Byte* pRgbRow);
#endif
static void ConvertYUVtoARGB(
    int width,
    int height,
    const Byte* ypixels, int ystride,
    const Byte* uvpixels, int uvstride, 
    Byte* outPixels, int outStride)
{
    // TODO: two scanlines at a time, compute single base rgb with u,v components
    // then incorporate y; yp = (y|16)-16 
    for (int ya = 0; ya < height; ya++)
    {
        const Byte* pYRow = ypixels + (ya * ystride);
        const Byte* pUVRow = uvpixels + ((ya>>1) * uvstride);
        Byte* pRgbRow = (Byte*)(outPixels + (ya * outStride));

        int x = 0;
#if 1
        // faster but slightly less accurate YUV->RGB conversion; 
        //
        // the primary optimization is using 16 bit accumulation instead of 32
        // bit (which speeds up the Neon/SSE versions), which is enabled by
        // not using the bottom two bits of the color transform factors to
        // enable the color transform accumulation to fit in signed 16 bits;
        // tests show this version varies by <= 2.001f/255.f
#if defined(_M_ARM)
        int width32 = width & ~(0x1f);
        NeonYUVToRGBALineFast(width32,pYRow,pUVRow,pRgbRow);
        x += width32;
        if (x < width)
        {
            pYRow += width32;
            pUVRow += width32;
            pRgbRow += (4*width32);
        }
#endif
        // combine all constant values plus rounding bias to constants to
        // initialize accumulators; need to do rounded divide by 4 to
        // accumulate to fit signed value in 16 bit
        // TODO: tweak these constants to get a better match with reference
        const short rb = ((((-298*16)           +(-409*128))+0x2)>>2) + 0x20;
        const short gb = ((((-298*16)+( 100*128)+( 208*128))+0x2)>>2) + 0x20;
        const short bb = ((((-298*16)+(-516*128)           )+0x2)>>2) + 0x20;
        // same divide by 4 for channel data scale coefficients
        const short ys  = ((298)+0x2)>>2; // Y scale for all three result channels
        const short vsr = ((409)+0x2)>>2; // V scale for red
        const short usg = ((100)+0x2)>>2; // U scale for green (to be subtracted)
        const short vsg = ((208)+0x2)>>2; // V scale for green (to be subtracted)
        const short usb = ((516)+0x2)>>2; // U scale for blue
#if (defined(_M_IX86) || defined(_M_AMD64))
        while (x <= (width-8))
        {
            __m128i t0,t1,t2,t3;

            // load 8 Y values and expand to 16bpp
            t0 = _mm_loadl_epi64((const __m128i*)pYRow); pYRow += 8;
            t0 = _mm_unpacklo_epi8(t0,_mm_setzero_si128());
            //t0 =  _mm_cvtepu8_epi16(t0); // SSE4.1

            // initialize accumulators
            __m128i racc = _mm_set1_epi16(rb);
            __m128i gacc = _mm_set1_epi16(gb);
            __m128i bacc = _mm_set1_epi16(bb);

            // scale Y by ys and accumulate
            t0 = _mm_mullo_epi16(t0, _mm_set1_epi16(ys));
            racc = _mm_add_epi16(racc, t0);
            gacc = _mm_add_epi16(gacc, t0);
            bacc = _mm_add_epi16(bacc, t0);

            // load 4 UV values; de-interleave, replicate, and convert to 16bpp
            t1 = _mm_loadl_epi64((const __m128i*)pUVRow); pUVRow += 8;
            t2 = _mm_setr_epi8(0,-1,0,-1, 2,-1,2,-1, 4,-1,4,-1, 6,-1,6,-1);
            t0 = _mm_shuffle_epi8( t1, t2 ); // replicated U to 8x16bpp
            t2 = _mm_setr_epi8(1,-1,1,-1, 3,-1,3,-1, 5,-1,5,-1, 7,-1,7,-1);
            t1 = _mm_shuffle_epi8( t1, t2 ); // replicated V to 8x16bpp

            // accumulate U contribution
            t2 = _mm_mullo_epi16(t0, _mm_set1_epi16(usg));
            gacc = _mm_sub_epi16(gacc, t2);
            t2 = _mm_mullo_epi16(t0, _mm_set1_epi16(usb));
            bacc = _mm_add_epi16(bacc, t2);

            // accumulate V contribution
            t2 = _mm_mullo_epi16(t1, _mm_set1_epi16(vsr));
            racc = _mm_add_epi16(racc, t2);
            t2 = _mm_mullo_epi16(t1, _mm_set1_epi16(vsg));
            gacc = _mm_sub_epi16(gacc, t2);

            // align to 8bpp
            racc = _mm_srai_epi16(racc, 6);
            gacc = _mm_srai_epi16(gacc, 6);
            bacc = _mm_srai_epi16(bacc, 6);

            // clamp to 0..255 and pack
            t0 = _mm_packus_epi16(bacc,racc); // b0..b7 | r0..r7
            // TODO: can use any available constant for alpha that is >= 255
            t1 = _mm_packus_epi16(gacc,_mm_set1_epi16(0xff)); // g0..g7 | ff..ff

            // shuffle so pixel 0..3 is grouped in bottom half and 4..7 in top
            t2 = _mm_shuffle_epi32(t0,0xd8); // b0..b3|r0..r3|b4..b7|r4..r7
            t3 = _mm_shuffle_epi32(t1,0xd8); // g0..g3|ff..ff|g4..g7|ff..ff

            // unpack to get all of rgba 0..3 or 4..7 in one register, then
            // shuffle to interleave
            t1 = _mm_setr_epi8(0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15);
            t0 = _mm_unpacklo_epi32(t2,t3); // b0..b3|g0..g3|r0..r3|ff..ff
            t0 = _mm_shuffle_epi8(t0,t1);
            _mm_storeu_si128((__m128i*)pRgbRow, t0); pRgbRow += 16;
            t0 = _mm_unpackhi_epi32(t2,t3); // b4..b7|g4..g7|r4..r7|ff..ff
            t0 = _mm_shuffle_epi8(t0,t1);
            _mm_storeu_si128((__m128i*)pRgbRow, t0); pRgbRow += 16;

            x += 8;
        }
#endif
        for (; x < width; x++)
        {
            // load yuv data
            short y = (short)(*pYRow);
            short u = (short)(*(pUVRow+0));
            short v = (short)(*(pUVRow+1));

            // initialize accumulators
            short r = rb;
            short g = gb;
            short b = bb;

            // accumulate Y contribution
            short yp = y*ys;
            r += yp;
            g += yp;
            b += yp;

            // accumulate U contribution
            g -= u*usg;
            b += u*usb;

            // accumulate V contribution
            r += v*vsr;
            g -= v*vsg;

            // normalize for 8 bit output (rounding bias already included)
            r >>= 6;
            g >>= 6;
            b >>= 6;
#else // slower but more accurate
#if defined(_M_ARM)
        int width32 = width & ~(0x1f);
        NeonYUVToRGBALine(width32,pYRow,pUVRow,pRgbRow);
        x += width32;
        if (x < width)
        {
            pYRow += width32;
            pUVRow += width32;
            pRgbRow += (4*width32);
        }
#endif
        for (; x < width; x++)
        {
            Byte y = *pYRow;
            Byte u = *(pUVRow+0);
            Byte v = *(pUVRow+1);

            int c = y - 16;
            int d = u - 128;
            int e = v - 128;
            
            int r = (298 * c           + 409 * e + 128) >> 8;
            int g = (298 * c - 100 * d - 208 * e + 128) >> 8;
            int b = (298 * c + 516 * d           + 128) >> 8;
#endif // faster/slower
#define CLAMP_U8(__val) (((__val)<0)?(0):((__val>0xff)?(0xff):(__val)))
            pRgbRow[0] = (Byte)CLAMP_U8(b);
            pRgbRow[1] = (Byte)CLAMP_U8(g);
            pRgbRow[2] = (Byte)CLAMP_U8(r);
            pRgbRow[3] = 255; // Alpha                
            pYRow++;
            if (x & 0x1) { pUVRow += 2; }
            pRgbRow += 4;
        }
    }
}

HRESULT vt::VtConvertImageNV12ToRGBA(CRGBAByteImg& imgDst, 
                                 const CLumaByteImg& imgSrcY,
                                 const CUVByteImg& imgSrcUV)
{
    VT_HR_BEGIN()

    VT_HR_EXIT( (!imgSrcY.IsValid()) ? E_INVALIDSRC : S_OK );
    VT_HR_EXIT( (!imgSrcUV.IsValid()) ? E_INVALIDSRC : S_OK );
    VT_HR_EXIT( (imgSrcY.IsSharingMemory(imgDst)) ? E_INVALIDDST : S_OK );
    VT_HR_EXIT( (imgSrcUV.IsSharingMemory(imgDst)) ? E_INVALIDDST : S_OK );
    // create the destination image if it needs to be
    VT_HR_EXIT( CreateImageForTransform(imgDst, imgSrcY.Width(), imgSrcY.Height(),
                                        imgDst.GetType()) );


    int w = imgSrcY.Width();
    int h = imgSrcY.Height();
    ConvertYUVtoARGB(w,h,
        imgSrcY.BytePtr(),imgSrcY.StrideBytes(),
        imgSrcUV.BytePtr(),imgSrcUV.StrideBytes(),
        imgDst.BytePtr(), imgDst.StrideBytes());

    VT_HR_END()
}

#ifndef VT_NO_XFORMS

void CConvertImageNV12toRGBATransform::GetSrcPixFormat(IN OUT int* pfrmtSrcs, 
                                 IN UINT uSrcCnt,
                                 IN int /*frmtDst*/)
{
    if (!m_NoRW)
    { 
        VT_ASSERT(uSrcCnt == 2); uSrcCnt;
        pfrmtSrcs[0] = VT_IMG_FIXED(VT_IMG_MAKE_COMP_TYPE(PIX_FORMAT_LUMA, EL_FORMAT_BYTE, 1)); 
        pfrmtSrcs[1] = VT_IMG_FIXED(VT_IMG_MAKE_COMP_TYPE(PIX_FORMAT_UV, EL_FORMAT_BYTE, 2)); 
    }
}
void CConvertImageNV12toRGBATransform::GetDstPixFormat(OUT int& frmtDst,
                                IN  const int* /*pfrmtSrcs*/, 
                                IN  UINT uSrcCnt)
{
    if (!m_NoRW)
    { 
        VT_ASSERT(uSrcCnt == 2); uSrcCnt;
        frmtDst = VT_IMG_FIXED(VT_IMG_MAKE_COMP_TYPE(PIX_FORMAT_RGBA, EL_FORMAT_BYTE, 4));
    }
}

HRESULT CConvertImageNV12toRGBATransform::GetRequiredSrcRect(OUT TRANSFORM_SOURCE_DESC* pSrcReq,
                                    OUT UINT& uSrcReqCount,
                                    IN  UINT uSrcCnt,
                                    IN  const CRect& rctLayerDst
                                    )
{ 
    if (!m_NoRW)
    {
        VT_ASSERT(uSrcCnt == 2); uSrcCnt;
        uSrcReqCount = 2;

        pSrcReq[0].bCanOverWrite = true;
		pSrcReq[0].rctSrc        = rctLayerDst;
		pSrcReq[0].uSrcIndex     = 0;
    	pSrcReq[1].bCanOverWrite = true;
		pSrcReq[1].rctSrc        = 
            vt::CRect(rctLayerDst.left/2, rctLayerDst.top/2, 
            rctLayerDst.right/2, rctLayerDst.bottom/2);
		pSrcReq[1].uSrcIndex     = 1;
    }
    else
    {
        uSrcReqCount = 0;
    }
    return S_OK;
}
HRESULT CConvertImageNV12toRGBATransform::GetResultingDstRect(OUT CRect& rctDst,
	                                IN  const CRect& rctSrc,
	                                IN  UINT uSrcIndex,
									IN  UINT /*uSrcCnt*/)
{
    if (uSrcIndex == 0) { rctDst = rctSrc; }
    else 
    {
        rctDst = vt::CRect(rctSrc.left*2,rctSrc.top*2,
            rctSrc.right*2,rctSrc.bottom*2);
    }
	return S_OK;
}
		
HRESULT CConvertImageNV12toRGBATransform::GetAffectedDstRect(OUT CRect& rctDst,
                                    IN  const CRect& rctSrc,
                                    IN  UINT uSrcIndex,
                                    IN  UINT uSrcCnt)
{ 
    GetResultingDstRect(rctDst,rctSrc,uSrcIndex,uSrcCnt);
    return S_OK;
}


HRESULT CConvertImageNV12toRGBATransform::Transform(OUT CImg* pimgDstRegion, 
                              IN  const CRect& rctLayerDst,
                              IN  CImg *const *ppimgSrcRegions,
                              IN  const TRANSFORM_SOURCE_DESC* /*pSrcDesc*/,
                              IN  UINT  /*uSrcCnt*/)
{ 
    VT_HR_BEGIN()

    if (!m_NoRW)
    {
        VT_HR_EXIT( VtConvertImageNV12ToRGBA((CRGBAByteImg&)*pimgDstRegion,
            (CLumaByteImg&)*ppimgSrcRegions[0],
            (CUVByteImg&)*ppimgSrcRegions[1]) );
    }
    else
    {
        // transform framework has no source or dest reader/writer in this case,
        // so need to share out the dest rect only from the destination for the
        // processing routine
        CRGBAByteImg dstImgR; VT_HR_EXIT( m_dstImg.Share(dstImgR, rctLayerDst) );
        CLumaByteImg srcImgYR; VT_HR_EXIT( m_srcYImg.Share(srcImgYR, rctLayerDst) );
        CUVByteImg srcImgUVR; VT_HR_EXIT( m_srcUVImg.Share(srcImgUVR, 
            CRect(rctLayerDst.left/2,rctLayerDst.top/2,
                  rctLayerDst.right/2,rctLayerDst.bottom/2)) );

        VT_HR_EXIT( VtConvertImageNV12ToRGBA(dstImgR,srcImgYR,srcImgUVR) );
    }
    VT_HR_END()
}

#endif

