//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Class to read in video frames using Media Foundation (MF).
//
//  History:
//      2010/06/25-t-parasa
//          Created
//      2011/07/17-sbaker
//          Added support for Video Images C*VImg
//          Added CMFVideoDst
//
//------------------------------------------------------------------------
#pragma once

#include "vtcommon.h"

#if defined(VT_WINRT)
#include "TempFileHelper.h"
#endif

// enable use of the DXVA2 hardware video decode; only for WinRT for now
// since dxgi device manager is not available on Win7 so need to use d3d9
// device manager for desktop
#if defined(VT_WINRT)
#define ENABLE_HARDWARE_DECODER
#if defined(_M_ARM)
// DXVA2 decode with NV12->RGBA conversion is very slow on Surface - much faster
// to do NV12 decode and then invoke Transform based NV12->RGBA conversion
#define MFVIDEOSRC_DECODE_TO_NV12
//  DXVA resize works well on Surface, but exhibits very poor quality filtering
//  on many x86 devices, both for NV12 and RGB32, so is enabled for ARM only
#define ENABLE_HARDWARE_RESIZE
#endif
#endif

#if ( defined(ENABLE_HARDWARE_DECODER) && (_MSC_VER >= 1700) )
#include "d3d11.h"
#endif

namespace vt {

//=============================================================================
// CCopyTransform - used for faster CopyTo on Arm; 
// TODO: move this elsewhere for more general use
//=============================================================================

class CCopyTransform: public IImageTransform
{
public:
    CCopyTransform()
        :m_bypassCache(false), m_noCopy(false)
    {}

    virtual bool RequiresCloneForConcurrency()
    { return false; }

    virtual void GetSrcPixFormat(IN OUT int* pfrmtSrcs, 
                                 IN UINT /*uSrcCnt*/,
                                 IN int frmtDst)
    { if (!m_noCopy) { pfrmtSrcs[0] = frmtDst; } }
    virtual void    GetDstPixFormat(OUT int& frmtDst,
                                    IN  const int* pfrmtSrcs, 
                                    IN  UINT  /*uSrcCnt*/)
    { if (!m_noCopy) { frmtDst = pfrmtSrcs[0]; } }
    virtual HRESULT GetRequiredSrcRect(OUT TRANSFORM_SOURCE_DESC* pSrcReq,
                                      OUT UINT& uSrcReqCount,
                                      IN  UINT /*uSrcCnt*/,
                                      IN  const CRect& rctLayerDst
                                      )
    {
        if (m_noCopy)
        {
            uSrcReqCount = 0;
        }
        else
        {
            vt::CRect rctRqdSrc = rctLayerDst;
            uSrcReqCount = 1;
            pSrcReq[0].bCanOverWrite = true;
            pSrcReq[0].rctSrc = rctRqdSrc;
            pSrcReq[0].uSrcIndex = 0;
        }
        return S_OK;
    }
    virtual HRESULT GetResultingDstRect(OUT CRect& rctDst,
                                        IN  const CRect& rctSrc,
                                        IN  UINT /*uSrcIndex*/,
                                        IN  UINT /*uSrcCnt*/)
    {
        rctDst = rctSrc;
        return S_OK;
    }
    virtual HRESULT GetAffectedDstRect(OUT CRect& rctDst,
                                      IN  const CRect& rctSrc,
                                      IN  UINT uSrcIndex,
                                      IN  UINT uSrcCnt)
    { 
        GetResultingDstRect(rctDst,rctSrc,uSrcIndex,uSrcCnt);
        return S_OK;
    }

    virtual HRESULT Transform(OUT CImg* pimgDstRegion, 
                              IN  const CRect& rctLayerDst,
                              IN  CImg *const * ppimgSrcRegions,
                              IN  const TRANSFORM_SOURCE_DESC* /*pSrcDesc*/,
                              IN  UINT /*uSrcCnt*/ )
    { 
        if (m_noCopy)
        {
            int wb = rctLayerDst.Width()*m_imgd.Bands()*m_imgd.ElSize();
            for (int y = rctLayerDst.top; y < rctLayerDst.bottom; y++)
            {
                void* pSrc = (void*)m_imgs.BytePtr(rctLayerDst.left,y);
                void* pDst = (void*)m_imgd.BytePtr(rctLayerDst.left,y);
                VtMemcpy(pDst,pSrc,wb,m_bypassCache);
            }
        }
        else
        {
            int wb = rctLayerDst.Width()*pimgDstRegion->Bands()*pimgDstRegion->ElSize();
            for (int y = 0; y < rctLayerDst.Height(); y++)
            {
                void* pSrc = (void*)ppimgSrcRegions[0]->BytePtr(0,y);
                void* pDst = (void*)pimgDstRegion->BytePtr(0,y);
                VtMemcpy(pDst,pSrc,wb,m_bypassCache);
            }
        }
        return S_OK; 
    }

    virtual HRESULT Clone(ITaskState ** /*ppState*/ )
    {
        return E_FAIL;
    }

    HRESULT Initialize(bool bypassCache = false)
    { 
        m_noCopy = false;
        m_bypassCache = bypassCache; 
        return S_OK;
    }
    HRESULT Initialize(CImg& dst, CImg& src, bool bypassCache = false)
    { 
        VT_HR_BEGIN();
        m_noCopy = true;
        VT_HR_EXIT( (dst.Width() == src.Width())?(S_OK):(E_INVALIDARG) );
        VT_HR_EXIT( (dst.Height() == src.Height())?(S_OK):(E_INVALIDARG) );
        VT_HR_EXIT( (dst.GetType() == src.GetType())?(S_OK):(E_INVALIDARG) );
        VT_HR_EXIT( dst.Share(m_imgd) ); 
        VT_HR_EXIT( src.Share(m_imgs) ); 
        m_bypassCache = bypassCache; 
        VT_HR_END();
    }

private:
    CImg m_imgs;
    CImg m_imgd;
    bool m_bypassCache;
    bool m_noCopy;
};

class CMFVideoSrc : public IVideoSrc
{
private:

    // A data structure to hold misc data about the video that might 
    //  need resetting once Close() is called. The resetting is
    //	implicitly done as we create a new instance of this object.
    struct VIDEO_FRAME_DATA
    {
        DOUBLE		duration;			// in seconds
        UINT32		width;				// in native (non-square) pixels
        UINT32		height;				// in native (non-square) pixels
        CRect		cropRect;			// in native (non-square) pixels
        int         cropRectType;
        DOUBLE		frameRate;			// in frames per second
        LONG		stride;				// in bytes per scan line
        DOUBLE		pixelAspectRatio;	// width/height of a native pixel
        LONGLONG	seekTolerance;		// hundreds of nanoseconds per frame
        UINT32      bitrate; 
        UINT32      interlaceMode;

        VIDEO_FRAME_DATA() :
            duration(0),
            width(0),
            height(0),
            frameRate(0),
            stride(0),
            pixelAspectRatio(0),
            seekTolerance(DEFAULT_SEEK_TOLERANCE)
        { }
        
        static const LONGLONG DEFAULT_SEEK_TOLERANCE = 1000000;
    };
    CVideoImgInfo& CMFVideoSrc::VideoImgInfo(int w = -1, int h = -1);
    CVideoImgInfo       m_VideoImgInfo;

    bool				m_bMFInitialized;
    IMFSourceReader*	m_pReader;
    VIDEO_FRAME_DATA	m_videoFrameData;
    vt::wstring         m_strCurrentFilename;
    bool				m_bCanSeek;
    LONG				m_iRefCount;
    VideoFormat         m_eVideoFormat;
    GUID                m_gFrameFormat;
    bool                m_bComStarted;
    bool                m_bSoftwareResize;
    int                 m_resizeWidth;
    int                 m_resizeHeight;
    int                 m_origWidth;
    int                 m_origHeight;
    int                 m_rotationMultipleOf90;

#if ( defined(ENABLE_HARDWARE_DECODER) && (_MSC_VER >= 1700) )
    CComPtr<ID3D11Device>          m_pD3DDev;
    CComPtr<ID3D11DeviceContext>   m_pD3DDevCon;
    CComPtr<IMFDXGIDeviceManager>  m_pdxgiDevMan;
#endif

public:

    CMFVideoSrc();
    ~CMFVideoSrc();

    // Methods to open and close a video file.
    virtual HRESULT OpenFile(IN const WCHAR* pwszFileName, VideoFormat eVideoFormat = CIMG);
#if defined(MAKE_DOXYGEN_WORK) || defined(VT_WINRT)
    virtual HRESULT OpenFile(Windows::Storage::IStorageFile^ storageFile, VideoFormat eVideoFormat = CIMG);
#endif
    virtual HRESULT Close();

    // Return a clone of this video source that references the same opened
    // file (if a file is open)
    virtual HRESULT Clone(IVideoSrc** ppClone);

    // Methods to get information about the video file.
    virtual VideoFormat GetVideoFormat() { return m_eVideoFormat; }
    virtual const WCHAR* GetFileName();
    virtual	HRESULT CanSeek(OUT bool& canSeek);
    virtual HRESULT GetPixelAspectRatio(OUT double& pixelAspectRatio);
    virtual HRESULT GetDuration(OUT double& durationInSeconds);
    virtual HRESULT GetFrameCount(OUT vt::Int32& frameCount);
    virtual HRESULT GetFrameRate(OUT double& framesPerSecond);
    virtual HRESULT GetBitrate(OUT float& bitsPerPixel);
    virtual HRESULT GetInterlaceMode(OUT int& interlaceMode);
    virtual HRESULT GetFrameSize(OUT vt::Int32& width, OUT vt::Int32& height);
    virtual HRESULT GetRawFrameSize(OUT vt::Int32& width, OUT vt::Int32& height);

    // Methods to modify video source
    virtual HRESULT SetFrameSize(IN vt::Int32 width, IN vt::Int32 height, bool bSoftwareResize = false);
    virtual HRESULT SetFrameRotation(IN vt::Int32 multipleOf90);

    // Methods to seek and get frames.
    virtual HRESULT Seek(IN double desiredFrameTimeInSeconds);
    virtual HRESULT GetNextFrame(IN OUT vt::CImg& img, 
                                 OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetNextFrame(IN OUT vt::CRGB32VideoImg& img, 
                                 OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetNextFrame(IN OUT vt::CNV12VideoImg& img, 
                                 OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetNextFrame(IN OUT vt::IImageWriter* pDst, 
                                 OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetNextFrame(IN OUT vt::IImageWriter* pDstY, 
                                 IN OUT vt::IImageWriter* pDstUV,
                                 OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CRGB32VideoImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CNV12VideoImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN vt::Int32 frameIndex, IN OUT vt::IImageWriter* pDst, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN vt::Int32 frameIndex, 
                             IN OUT vt::IImageWriter* pDstY, 
                             IN OUT vt::IImageWriter* pDstUV, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::CImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::CRGB32VideoImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::CNV12VideoImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, 
                             IN OUT vt::IImageWriter* pDst, 
                             OUT double* pActualFrameTimeInSeconds = NULL);
    virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, 
                             IN OUT vt::IImageWriter* pDstY, 
                             IN OUT vt::IImageWriter* pDstUV, 
                             OUT double* pActualFrameTimeInSeconds = NULL);

    // Reference counting.
    virtual ULONG AddRef();
    virtual ULONG Release();

private:
    template <typename T> HRESULT GetNextFrameInternal(T** img, 
        double* pActualFrameTimeInSeconds, VideoFormat fmt);
    template <typename T> HRESULT GetFrameInternal(
        double desiredFrameTimeInSeconds, T** img,
        double* pActualFrameTimeInSeconds, VideoFormat fmt);
    template <typename T> HRESULT GetFrameInternal(
        Int32 frameIndex, T** img,
        double* pActualFrameTimeInSeconds, VideoFormat fmt);
    HRESULT OpenByteStream(IMFByteStream* pByteStream, VideoFormat eVideoFormat);
    HRESULT ReadPresentationAttributes();
    HRESULT SelectVideoStream();
    HRESULT ReadMediaTypeAttributes();

    bool m_bTMStarted; // true if task manager has been started
    struct VIDEO_PROCESS_TRANSFORM_GRAPH
    {
        CTransformGraphNode                 *pTopNode;
        CTransformGraphNode                 *pTopNodeUV;
        CTransformGraphNaryNode             nodeCvt; // Nary node for NV12toRGBA
        CTransformGraphUnaryNode            node[2]; // Unary nodes for copy, resize, rotate
        CTransformGraphUnaryNode            nodeuv[2];
        CCopyTransform                      xformCpy;
        CConvertImageNV12toRGBATransform    xformCvt;
        CRotateTransform                    xformRot[2];
        CWarpTransform                      xformRsz[2];
        IMAGE_EXTEND                        ex;
        VIDEO_PROCESS_TRANSFORM_GRAPH():pTopNode(nullptr),pTopNodeUV(nullptr),ex(Extend) {}
    } m_vpgraph;
    HRESULT SetProcessVideoFrameTransform(IImageWriter* pDst, IImageWriter* pDstUV,
        IImageReader* pSrc, IImageReader* pSrcUV, bool bFlip);
    HRESULT ProcessVideoFrame(IImageWriter* pDst, IImageWriter* pDstUV, BYTE* pSrc, int signedStride);

    HRESULT ConvertSampleToImage(vt::CImg** ppimg, BYTE* pBitmapData, LONG actualStride);
    HRESULT ConvertSampleToImage(vt::CRGB32VideoImg** ppimg, BYTE* pBitmapData, LONG actualStride);
    HRESULT ConvertSampleToImage(vt::CNV12VideoImg** ppimg, BYTE* pBitmapData, LONG actualStride);
    HRESULT ConvertSampleToImage(vt::IImageWriter** ppDst, BYTE* pBitmapData, LONG actualStride);
    HRESULT InitializeMediaFoundation();
    HRESULT ShutdownMediaFoundation();
};


class CMFVideoDst : public IVideoDst
{
public:
    CMFVideoDst();
    ~CMFVideoDst();

    virtual HRESULT Clone(IVideoDst** ppClone);

    virtual HRESULT OpenFile(__in_z const wchar_t *pwcFilename, int iValidPixelWidth, int iValidPixelHeight,
                             VideoFormat eVideoFormat, int iFramesPerSecond = 25, float fBitsPerPixel = 4.0f,
                             double dAspectRatio = 1.0, int iInterlaceMode = MFVideoInterlace_Progressive);
#if defined(MAKE_DOXYGEN_WORK) || defined(VT_WINRT)
    virtual HRESULT OpenFile(Windows::Storage::IStorageFile^ storageFile, int iValidPixelWidth, int iValidPixelHeight,
                             VideoFormat eVideoFormat, int iFramesPerSecond = 25, float fBitsPerPixel = 4.0f,
                             double dAspectRatio = 1.0, int iInterlaceMode = MFVideoInterlace_Progressive);
#endif
    virtual HRESULT Close();

    virtual HRESULT WriteFrame(CImg &img);
    virtual HRESULT WriteFrame(CRGB32VideoImg &img);
    virtual HRESULT WriteFrame(CNV12VideoImg &img);

    // Reference counting.
    virtual ULONG AddRef();
    virtual ULONG Release();

protected:
    HRESULT OpenByteStreamOrFile(__in_opt IMFByteStream *pByteStream, __in_z const wchar_t *pwcFilename,
                                 int iValidPixelWidth, int iValidPixelHeight,
                                 VideoFormat eVideoFormat, int iFramesPerSecond, float fBitsPerPixel,
                                 double dAspectRatio, int iInterlaceMode);
    HRESULT InitializeMediaFoundation();
    HRESULT ShutdownMediaFoundation();

private:
    // Data Members
    LONG m_iRefCount;
    bool m_bComStarted;
    bool m_bMFInitialized;
    IMFSinkWriter* m_pSinkWriter;
    UINT32 m_iWidth;
    UINT32 m_iHeight;
    UINT32 m_iFramesPerSecond;
    DWORD m_streamIndex;
    LONGLONG m_rtStart;
    UINT64 m_rtDuration;
    UINT32 m_iBitRate;
    VideoFormat m_eVideoFormat;
    GUID m_videoOutputFormat;
#if defined(VT_WINRT)
    CTempFileHelper m_tempFileHelper;
#endif
};

bool GetMFErrorString(HRESULT hr, wstring& errmsg);

};

