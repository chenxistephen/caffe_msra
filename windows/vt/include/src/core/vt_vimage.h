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
#pragma once

#include "vtcommon.h"
#include "vt_image.h"
#include "vt_atltypes.h"

#ifndef VT_NO_XFORMS
#include "vt_transform.h"
#endif

namespace vt {

//+-----------------------------------------------------------------------
//
//  CVideoImgInfo
//
//  See http://msdn.microsoft.com/en-us/library/aa376629(v=VS.85).aspx
//  for details of media type attributes
//------------------------------------------------------------------------

/// <summary> Holds information about a video image
/// </summary>
class CVideoImgInfo
{
public:
    /// <summary> The type of the valid pixels rectangle from Media Foundation </summary>
    enum ValidPixelsRectType
    {
        /// <summary> See http://msdn.microsoft.com/en-us/library/ms705634(v=VS.85).aspx </summary>
        PanScan,                      
        /// <summary> See http://msdn.microsoft.com/en-us/library/ms700173(v=VS.85).aspx </summary>
        MinDisplay, 
        /// <summary> See http://msdn.microsoft.com/en-us/library/ms701632(v=VS.85).aspx </summary>
        Geometric           
    };

    /// <summary> Constructor </summary>
    CVideoImgInfo() : eValidPixelsRectType(MinDisplay), dAspectRatio(1.0),
            iInterlaceMode(2)
    {
        rectValidPixels = CRect(0,0,0,0);
    }

    /// <summary> A rectangle specifying which pixels are valid </summary>
    RECT  rectValidPixels;
    /// <summary> The type of the valid pixels rectangle from Media Foundation </summary>
    ValidPixelsRectType eValidPixelsRectType;        
    /// <summary> The aspect ratio. See http://msdn.microsoft.com/en-us/library/ms704767(v=VS.85).aspx </summary>
    double dAspectRatio; 
    /// <summary> The Media Foundation interlace model. See http://msdn.microsoft.com/en-us/library/ms694269(v=VS.85).aspx </summary>
    int    iInterlaceMode;             
};


//+-----------------------------------------------------------------------
//
//  CRGB32VideoImg
//
//------------------------------------------------------------------------

/// <summary> Holds an RGB32 video image
/// </summary>
class CRGB32VideoImg
{
public:
    /// <summary> Constructor </summary>
    CRGB32VideoImg()  {}
    /// <summary> Destructor </summary>
    ~CRGB32VideoImg() { CheckInvariant(); }

    /// <summary> Create an RGB32 format video image</summary>
    /// <param name="iWidth"> Width </param>
    /// <param name="iHeight"> Height </param>
    /// <returns> 
    ///     - S_OK on success
    ///     - Use VtIOErrorToString() to get extended error information.
    /// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    ///     - Width and height must both be even
    HRESULT Create(int iWidth, int iHeight);
    /// <summary> Create an RGB32 format video image, wrap existing pixel buffer.  </summary>
    /// <param name="pbBuffer"> Pointer to Pixel Buffer </param>
    /// <param name="iWidth"> Width </param>
    /// <param name="iHeight"> Height </param>
    /// <param name="iStrideBytes"> Number of bytes per stride</param>
    /// <returns> 
    ///     - S_OK on success
    ///     - E_INVALIDARG if requested properties are invalid or if pbBuffer is NULL.
    /// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    ///     - No memory is allocated
    ///     - Width and height must both be even
    HRESULT Create(Byte *pbBuffer, int iWidth, int iHeight, int iStrideBytes);
    /// <summary> Create an RGB32 format video image</summary>
    /// <param name="iWidth"> Width </param>
    /// <param name="iHeight"> Height </param>
    /// <param name="vi"> A CVideoImgInfo object</param>
    /// <returns> 
    ///     - S_OK on success
    ///     - Use VtIOErrorToString() to get extended error information.
    /// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    ///     - Width and height must both be even
    HRESULT Create(int iWidth, int iHeight, const CVideoImgInfo& vi);
    /// <summary> Create an RGB32 format video image, wrap existing pixel buffer.  </summary>
    /// <param name="pbBuffer"> Pointer to Pixel Buffer </param>
    /// <param name="iWidth"> Width </param>
    /// <param name="iHeight"> Height </param>
    /// <param name="iStrideBytes"> Number of bytes per stride </param>
    /// <param name="vi"> A CVideoImgInfo object</param>
    /// <returns> 
    ///     - S_OK on success
    ///     - E_INVALIDARG if requested properties are invalid or if pbBuffer is NULL.
    /// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    ///     - No memory is allocated
    ///     - Width and height must both be even
    HRESULT Create(Byte *pbBuffer, int iWidth, int iHeight, int iStrideBytes, const CVideoImgInfo& vi);

    /// <summary> Deallocates memory </summary>
    void Deallocate()
    {
        CheckInvariant();
        m_imgVidInfo = CVideoImgInfo();
        m_imgRGB.Deallocate();
    }

    /// <summary> Read Access to the CVideoImgInfo object </summary>
    const CVideoImgInfo& GetVideoImgInfo() const
    { 
        CheckInvariant(); 
        return m_imgVidInfo;
    }

    /// <summary> Access to the CRGBAImg that holds the pixels </summary>
    CRGBAImg&  GetImg()
    { 
        CheckInvariant(); 
        return m_imgRGB;
    }
    /// <summary> Access to the CRGBAImg that holds the pixels </summary>
    const CRGBAImg&  GetImg() const
    {
        CheckInvariant(); 
        return m_imgRGB;
    }

protected:
    /// <summary> Checks the CRGBAImg is consistent with the CVideoImgInfo </summary>
    void CheckInvariant() const
    {
#ifdef _DEBUG
        VT_ASSERT( vt::CRect(m_imgRGB.Rect()).RectInRect(
                &m_imgVidInfo.rectValidPixels) );
#endif
    }

protected:
    /// <summary> The CVideoImgInfo object </summary>
    CVideoImgInfo m_imgVidInfo;
    /// <summary> A CRGBAImg that holds the pixels </summary>
    CRGBAImg      m_imgRGB;
};


//+-----------------------------------------------------------------------
//
//  CYV12VideoImg
//
//------------------------------------------------------------------------

/// <summary> Holds an NV12 video image
/// </summary>
class CNV12VideoImg
{
public:
    /// <summary> Constructor </summary>
    CNV12VideoImg()  {}
    /// <summary> Destructor </summary>
    ~CNV12VideoImg() { CheckInvariant(); } 

    /// <summary> Create an NV12 format video image</summary>
    /// <param name="iWidth"> Width </param>
    /// <param name="iHeight"> Height </param>
    /// <returns> 
    ///     - S_OK on success
    ///     - Use VtIOErrorToString() to get extended error information.
    /// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    ///     - Width and height must both be even
    HRESULT Create(int iWidth, int iHeight);
     /// <summary> Create an NV12 format video image, wrap existing pixel buffer.  </summary>
    /// <param name="pbBuffer"> Pointer to Pixel Buffer </param>
    /// <param name="iWidth"> Width </param>
    /// <param name="iHeight"> Height </param>
    /// <param name="iStrideBytes"> Number of bytes per stride</param>
    /// <returns> 
    ///     - S_OK on success
    ///     - E_INVALIDARG if requested properties are invalid or if pbBuffer is NULL.
    /// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    ///     - No memory is allocated
    ///     - Width and height must both be even
    HRESULT Create(Byte *pbBuffer, int iWidth, int iHeight, int iStrideBytes);
    /// <summary> Create an NV12 format video image</summary>
    /// <param name="iWidth"> Width </param>
    /// <param name="iHeight"> Height </param>
    /// <param name="vi"> A CVideoImgInfo object</param>
    /// <returns> 
    ///     - S_OK on success
    ///     - Use VtIOErrorToString() to get extended error information.
    /// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    ///     - Width and height must both be even
    HRESULT Create(int iWidth, int iHeight, const CVideoImgInfo& vi);
    /// <summary> Create an NV12 format video image, wrap existing pixel buffer.  </summary>
    /// <param name="pbBuffer"> Pointer to Pixel Buffer </param>
    /// <param name="iWidth"> Width </param>
    /// <param name="iHeight"> Height </param>
    /// <param name="iStrideBytes"> Number of bytes per stride</param>
    /// <param name="vi"> A CVideoImgInfo object</param>
    /// <returns> 
    ///     - S_OK on success
    ///     - E_INVALIDARG if requested properties are invalid or if pbBuffer is NULL.
    /// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    ///     - No memory is allocated
    ///     - Width and height must both be even
    HRESULT Create(Byte *pbBuffer, int iWidth, int iHeight, int iStrideBytes, const CVideoImgInfo& vi);

    /// <summary> Deallocates memory </summary>
    void Deallocate()
    {
        CheckInvariant();
        m_imgVidInfo = CVideoImgInfo();
        m_imgY.Deallocate();
        m_imgUV.Deallocate();
    }
 
   /// <summary> Read Access to the CVideoImgInfo object </summary>
    const CVideoImgInfo& GetVideoImgInfo() const
    { 
        CheckInvariant(); 
        return m_imgVidInfo;
    }
 
    /// <summary> Access to the CLumaByteImg that holds the Y pixels </summary>
    CLumaByteImg&  GetYImg()
    {
        CheckInvariant(); 
        return m_imgY;
    }
    /// <summary> Access to the CLumaByteImg that holds the Y pixels </summary>
    const CLumaByteImg&  GetYImg() const
    { 
        CheckInvariant(); 
        return m_imgY;
    }
    /// <summary> Access to the CUVByteImg that holds the UV pixels </summary>
    CUVByteImg&  GetUVImg()
    {
        CheckInvariant(); 
        return m_imgUV;
    }
    /// <summary> Access to the CUVByteImg that holds the UV pixels </summary>
    const CUVByteImg&  GetUVImg() const
    { 
        CheckInvariant(); 
        return m_imgUV;
    }

    /// <summary> Pixel memory access </summary>
    const Byte *BytePtr() const { return m_imgData.BytePtr(); }
    /// <summary> Pixel memory access </summary>
    Byte* BytePtr() { return const_cast<Byte*>( static_cast<const CNV12VideoImg&>(*this).BytePtr()); }

protected:
    /// <summary> Checks the Y and UV images are consistent with the CVideoImgInfo </summary>
    void CheckInvariant() const
    {
#ifdef _DEBUG
        VT_ASSERT( vt::CRect(m_imgY.Rect()).RectInRect(&m_imgVidInfo.rectValidPixels) );
        CRect validPixels = m_imgVidInfo.rectValidPixels;
        CRect uvValidPixels = CRect(validPixels.left/2, validPixels.top/2, validPixels.Width()/2, validPixels.Height()/2);
        VT_ASSERT( vt::CRect(m_imgUV.Rect()).RectInRect(&uvValidPixels) );
#endif
    }

protected:
    /// <summary> The CVideoImgInfo object </summary>
    CVideoImgInfo m_imgVidInfo;
    /// <summary> A CLumaByteImg that wraps the Y pixels </summary>
    CLumaByteImg m_imgY;
    /// <summary> A CUVByteImg that wraps the UV pixels </summary>
    CUVByteImg   m_imgUV;
    /// <summary> A CLumaImg that holds the NV12 data block </summary>
    CLumaImg     m_imgData; 
};


//+-----------------------------------------------------------------------
//
//  Conversion Routines
//  Relatively Slow (No SSE)
//  Use color conversion DMO where possible for efficient color processing 
//
//------------------------------------------------------------------------

/// \ingroup conversion
/// <summary> Converts an RGB32 format video image into a NV12 one</summary>
/// <param name="imgSrc"> The source RGB32 video image </param>
/// <param name="imgDst"> The destination NV12 video image </param>
/// <returns> 
///        - S_OK on success
///        - Use VtIOErrorToString() to get extended error information.
/// </returns>
HRESULT VtConvertVideoImage(CNV12VideoImg &imgDst,  const CRGB32VideoImg &imgSrc);

/// \ingroup conversion
/// <summary> Converts an NV12 format video image into a RGB32 one</summary>
/// <param name="imgSrc"> The source NV12 video image </param>
/// <param name="imgDst"> The destination RGB32 video image </param>
/// <returns> 
///     - S_OK on success
///     - Use VtIOErrorToString() to get extended error information.
/// </returns>
HRESULT VtConvertVideoImage(CRGB32VideoImg &imgDst, const CNV12VideoImg &imgSrc);

//+-----------------------------------------------------------------------------
//
// Fast NV12->RGBA Conversion
// 
//------------------------------------------------------------------------------

/// \ingroup conversion
/// <summary> Converts an NV12 image with separate Y and UV images to an RGBA image.</summary>
/// <param name="imgSrcY"> The source NV12 Y (Luminance) image </param>
/// <param name="imgSrcUV"> The source NV12 UV (Chrominance) image </param>
/// <param name="imgDst"> The destination RGB32 image </param>
/// <returns> 
///     - S_OK on success
///     - Use VtIOErrorToString() to get extended error information.
/// </returns>
HRESULT VtConvertImageNV12ToRGBA(CRGBAByteImg& imgDst, 
                                 const CLumaByteImg& imgSrcY,
                                 const CUVByteImg& imgSrcUV);

#ifndef VT_NO_XFORMS

/// \ingroup filtertransforms
/// <summary>Implementation of \link IImageTransform image transform \endlink that
/// converts an NV12 image to RGBA
/// </summary> 

//+-----------------------------------------------------------------------------
//
// Class: CConvertNV12ImagetoRGBATransform
// 
// Synposis: IImageTransform to convert an NV12 image to RGBA
// 
//------------------------------------------------------------------------------
class CConvertImageNV12toRGBATransform: public IImageTransform
{
	// IImageTransform implementation
public:
    virtual bool RequiresCloneForConcurrency()
    { return false; }

    virtual void GetSrcPixFormat(IN OUT int* pfrmtSrcs, 
                                 IN UINT /*uSrcCnt*/,
                                 IN int /*frmtDst*/);
    virtual void    GetDstPixFormat(OUT int& frmtDst,
                                    IN  const int* /*pfrmtSrcs*/, 
                                    IN  UINT  /*uSrcCnt*/);
    virtual HRESULT GetRequiredSrcRect(OUT TRANSFORM_SOURCE_DESC* pSrcReq,
                                      OUT UINT& uSrcReqCount,
                                      IN  UINT /*uSrcCnt*/,
                                      IN  const CRect& rctLayerDst
                                      );
	virtual HRESULT GetResultingDstRect(OUT CRect& rctDst,
	                                    IN  const CRect& rctSrc,
	                                    IN  UINT uSrcIndex,
										IN  UINT /*uSrcCnt*/);
    virtual HRESULT GetAffectedDstRect(OUT CRect& rctDst,
                                      IN  const CRect& rctSrc,
                                      IN  UINT uSrcIndex,
                                      IN  UINT uSrcCnt);
    virtual HRESULT Transform(OUT CImg* pimgDstRegion, 
                              IN  const CRect& rctLayerDst,
                              IN  CImg *const *ppimgSrcRegions,
                              IN  const TRANSFORM_SOURCE_DESC* pSrcDesc,
                              IN  UINT /*uSrcCnt*/ );

	virtual HRESULT Clone(ITaskState ** /*ppState*/ )
	{
		return E_FAIL;
    }

public:
    CConvertImageNV12toRGBATransform()
        :m_NoRW(false)
    {}

	/// <summary> Initialize transform</summary> 
	/// <returns> 
	///		- always returns S_OK
	/// </returns>
    HRESULT Initialize(void) { m_NoRW = false; return S_OK; };

	/// <summary> Initialize transform for use with CTransformGraphNoSrcNode</summary> 
	/// <returns> 
	///		- always returns S_OK
	/// </returns>
    HRESULT Initialize(CRGBAByteImg& dst, CLumaByteImg& srcY, CUVByteImg& srcUV)
    { m_NoRW = true; dst.Share(m_dstImg); srcY.Share(m_srcYImg); srcUV.Share(m_srcUVImg); return S_OK; }

protected:
    bool m_NoRW; // true if using no-copy
    CLumaByteImg m_srcYImg;
    CUVByteImg m_srcUVImg;
    CRGBAByteImg m_dstImg;
};

#endif


};
