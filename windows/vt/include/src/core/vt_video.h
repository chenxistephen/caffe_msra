//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Interface for reading frames from a video file.
//
//  History:
//      2010/08/24-t-parasa
//          Created
//      2011/07/17-sbaker
//          Added IVideoDst Interface
//      2011/08/16-sbaker
//          Added documentation
//
//------------------------------------------------------------------------
#pragma once

#if !(defined(WINAPI_FAMILY) && (WINAPI_FAMILY == WINAPI_FAMILY_PHONE_APP))
#include <Mferror.h>
#endif
#include <mfobjects.h>
#include "vtcommon.h"

namespace vt {

#define FIO_VIDEO_E_NOT_INITIALIZED MF_E_NOT_INITIALIZED
#define FIO_VIDEO_E_END_OF_STREAM MF_E_END_OF_STREAM
#define FIO_VIDEO_E_INVALID_VIDEO_FORMAT MF_E_INVALID_FORMAT
#define FIO_VIDEO_E_BYTESTREAM_NOT_SEEKABLE MF_E_BYTESTREAM_NOT_SEEKABLE

#ifndef VideoFormat
	/// <summary> A enumeration of all the currently supported video formats. </summary>	
	typedef enum {
		/// <summary> VFNone is an invalid format, used when initializing the 
		/// format in reader/writer objects. </summary>	
		VFNone,
		/// <summary> CIMG means that the video will be input/output to/from
		/// Vision Tools CImg objects. </summary>	
		CIMG,
		/// <summary> RGB32 corresponds to the Media Foundation RGB32 video
		/// subtype. Input/output is through CRGB32VideoImg objects. </summary>	
		RGB32,
		/// <summary> NV12 corresponds to the Media Foundation NV12 video
		/// subtype. Input/output is through CNV12VideoImg objects. </summary>	
		NV12
	} VideoFormat;
#endif

/// <summary> Interface to open a video file and read images. Most common
/// video formats are supported. </summary>	
class IVideoSrc
{
public:

	/// <summary> Opens a video file for reading. Sets up to be able to read
	/// the pixels into a type of object specified by the eVideoFormat. </summary>	
    /// <param name="pwszFileName"> Video filename </param>
    /// <param name="eVideoFormat"> The format that the pixels will be requested in.
	/// CIMG = Want to recieve pixels in a CImg. RGB32 = Want to received pixels in
	/// a CRGB32VideoImg object. NV12 = Want to receive pixels in a CNV12VideoImg
	/// object.</param>
	/// <DL><DT> Remarks: </DT></DL>
	/// Generally speaking, when processing video, it is best to first try to open
	/// the video file for NV12 and then if that fails, try RGB32
	virtual HRESULT OpenFile(IN const WCHAR* pwszFileName, 
                             VideoFormat eVideoFormat = CIMG) = 0;

#if defined(MAKE_DOXYGEN_WORK) || defined(VT_WINRT)
	/// <summary> Opens a video stream for reading. Sets up to be able to read
	/// the pixels into a type of object specified by the eVideoFormat. </summary>	
    /// <param name="storageFile"> Video file </param>
    /// <param name="eVideoFormat"> The format that the pixels will be requested in.
	/// CIMG = Want to recieve pixels in a CImg. RGB32 = Want to received pixels in
	/// a CRGB32VideoImg object. NV12 = Want to receive pixels in a CNV12VideoImg
	/// object.</param>
	/// <DL><DT> Remarks: </DT></DL>
	///	- Available only in WinRT configurations.
	/// - Generally speaking, when processing video, it is best to first try to open
	/// the video file for NV12 and then if that fails, try RGB32
	virtual HRESULT OpenFile(IN Windows::Storage::IStorageFile^ storageFile,
							 VideoFormat eVideoFormat = CIMG) = 0;
#endif

	/// <summary> Closes the video file or stream. </summary>	
	virtual HRESULT Close() = 0;

	/// <summary> Return a clone of this video source that references the same opened
    /// file (if a file is open). </summary>	
    virtual HRESULT Clone(IVideoSrc** ppClone) = 0;

	/// <summary> Returns the video format that the file has been opened for. </summary>	
    virtual VideoFormat GetVideoFormat() = 0;
	/// <summary> Returns the filename of the opened video file. </summary>	
    virtual const WCHAR* GetFileName() = 0;
	/// <summary> Returns whether seeking is possible. </summary>	
	virtual	HRESULT CanSeek(OUT bool& canSeek) = 0;
	/// <summary> Returns the pixel aspect ratio (1 means square pixels). </summary>	
	virtual HRESULT GetPixelAspectRatio(OUT double& pixelAspectRatio) = 0;
	/// <summary> Returns the duration of the video in seconds. </summary>	
	virtual HRESULT GetDuration(OUT double& durationInSeconds) = 0;
	/// <summary> Returns an approximate estimate of the number of frames in the 
	/// video. </summary>	
	/// <DL><DT> Remarks: </DT></DL>
	/// Determining the exact number of frames from the header of a video is in
	/// general not possible to do exactly. The estimate return by GetFrameCount()
	/// is probably an under-estimate. In general, the caller should always check
	/// whether GetNextFrame() (etc) succeeded. And also continue to try to process
	/// frames until GetNextFrame() fails. Relying on GetFrameCount() when 
	/// processing a video will often end up truncating the video (most likely
	/// by just 1 frame.
	virtual HRESULT GetFrameCount(OUT vt::Int32& frameCount) = 0;
	/// <summary> Returns the frame rate of the video in frames per second. </summary>	
	virtual HRESULT GetFrameRate(OUT double& framesPerSecond) = 0;
	/// <summary> Returns the approximate bitrate of the video in bits per pixel. </summary>	
	virtual HRESULT GetBitrate(OUT float& bitsPerPixel) = 0;
	/// <summary> Returns the media foundation interlace mode. 
	/// See http://msdn.microsoft.com/en-us/library/ms694868(v=VS.85).aspx. </summary>	
	virtual HRESULT GetInterlaceMode(OUT int& interlaceMode) = 0;
	/// <summary> Returns the frame size of the video frames that can be expected 
	/// after compensating for any cropping aperture and the aspect ratio. 
	/// GetFrameSize() is the size of the image returned in the CImg mode.</summary>	
	/// <DL><DT> Remarks: </DT></DL>
	/// In general, the frame size of a video may not be constant throughout the video.
	/// GetFrameSize() just returns the current frame size. 
	/// Even if the frame size is constant, the frame size specified in the header may
	/// well not equal the frame size of the actual frames. So, GetFrameSize() may
	/// return a different result after the first call to GetNextFrame() (compared
	/// to the value returned after opening the file. The best way to determine the
	/// frame size is to call GetNextFrame(). GetFrameSize() can be wrong.
	virtual HRESULT GetFrameSize(OUT vt::Int32& width, OUT vt::Int32& height) = 0;
	/// <summary> Returns the frame size of the raw video frames. 
	/// GetRawFrameSize() is the size of the image returned in the NV12 and RGB32 modes.
	/// </summary>	
	/// <DL><DT> Remarks: </DT></DL>
	/// In general, the frame size of a video may not be constant throughout the video.
	/// GetRawFrameSize() just returns the current frame size. 
	/// Even if the frame size is constant, the frame size specified in the header may
	/// well not equal the frame size of the actual frames. So, GetFrameSize() may
	/// return a different result after the first call to GetNextFrame() (compared
	/// to the value returned after opening the file. The best way to determine the
	/// frame size is to call GetNextFrame(). GetRawFrameSize() can be wrong when
	/// first called after opening the video.
	virtual HRESULT GetRawFrameSize(OUT vt::Int32& width, OUT vt::Int32& height) = 0;

	/// <summary> Sets a new frame size to be returned from subsequent Get*Frame
    /// calls.  If bSoftwareResize is true, then the resizing will alwoays be done
    /// via VtResizeImage or the resize Transform.  If bSoftware is false, the
    /// resize will be done with the DXVA video processing pipeline when possible.
    /// Width and Height must be even when requesting the DXVA pipeline resize.
	/// </summary>	
	/// <DL><DT> Remarks: </DT></DL>
    /// To set the frame size back to that of the source stream, call SetFrameSize
    /// passing -1 for the width.
	virtual HRESULT SetFrameSize(IN vt::Int32 width, IN vt::Int32 height, 
        IN bool bSoftwareResize = false) = 0;

	/// <summary> Sets a new frame rotation to be applied to the result for 
    /// subsequent Get*Frame calls.
	/// </summary>	
    /// <param name="multipleOf90">The image will be rotated clockwise by 
    /// multipleOf90 * 90 degrees</param>
	/// <DL><DT> Remarks: </DT></DL>
  	/// <DL><DT> Remarks: </DT></DL>
	///		- Rotates around (0,0) and translates the upper left corner of the
	///	resulting image back to (0,0)
    virtual HRESULT SetFrameRotation(IN vt::Int32 multipleOf90) = 0;

	/// <summary> Seek to a specific location in the video </summary>
	/// <param name="desiredFrameTimeInSeconds"> Position in the video in seconds </param>
	virtual HRESULT Seek(IN double desiredFrameTimeInSeconds) = 0; 
	/// <summary> Get the next frame and writes into a CImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a CIMG eVideoFormat.
	virtual HRESULT GetNextFrame(IN OUT vt::CImg& img, 
                                 OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Get the next frame and writes into a CRGB32VideoImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a RGB32 eVideoFormat.
	virtual HRESULT GetNextFrame(IN OUT vt::CRGB32VideoImg& img, 
                                 OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Get the next frame and writes into a CNV12VideoImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a NV12 eVideoFormat.
	virtual HRESULT GetNextFrame(IN OUT vt::CNV12VideoImg& img, 
                                 OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Get the next frame and writes into a CRGBAImg based IImageWriter. </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	virtual HRESULT GetNextFrame(IN OUT vt::IImageWriter* pDst, 
                                 OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Get the next frame and writes the Y result to a CLumaByteImg based
    /// IImageWriter, and the UV result to a CUVByteImg based IIImageWriter.</summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a NV12 eVideoFormat.
	virtual HRESULT GetNextFrame(IN OUT vt::IImageWriter* pDstY, 
                                 IN OUT vt::IImageWriter* pDstUV,
                                 OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given frame number based on estimating using
	/// framerate, get the next frame and write into a CImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a CIMG eVideoFormat.
	/// <param name="frameIndex"> The frame number to try to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	virtual HRESULT GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given frame number based on estimating using
	/// framerate, get the next frame and write into a CRGB32VideoImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a RGB32 eVideoFormat.
	/// <param name="frameIndex"> The frame number to try to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	virtual HRESULT GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CRGB32VideoImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given frame number based on estimating using
	/// framerate, get the next frame and write into a CNV12VideoImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a NV12 eVideoFormat.
	/// <param name="frameIndex"> The frame number to try to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	virtual HRESULT GetFrame(IN vt::Int32 frameIndex, IN OUT vt::CNV12VideoImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given frame number based on estimating using
	/// framerate, get the next frame and write to a CRGBAImg based IImageWriter object.
    /// </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// <param name="frameIndex"> The frame number to try to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	/// OpenFile() must have been called with a RGB32 or NV12 eVideoFormat.  If
    /// NV12 format is used, the video decoder will generate NV12 and a fast
    /// software conversion to RGB will be invoked for the returned frame data.
	virtual HRESULT GetFrame(IN vt::Int32 frameIndex, IN OUT vt::IImageWriter* pDstY,
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given frame number based on estimating using
	/// framerate, get the next frame and writes the Y result to a CLumaByteImg based
    /// IImageWriter, and the UV result to a CUVByteImg based IIImageWriter.</summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a NV12 eVideoFormat.
	/// <param name="frameIndex"> The frame number to try to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	virtual HRESULT GetFrame(IN vt::Int32 frameIndex, 
                             IN OUT vt::IImageWriter* pDstY, IN OUT vt::IImageWriter* pDstUV,
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given time, get the next frame and write into a 
	/// CImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a CIMG eVideoFormat.
	/// <param name="desiredFrameTimeInSeconds"> The time to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::CImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given time, get the next frame and write into a 
	/// CRGB32VideoImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a RGB32 eVideoFormat.
	/// <param name="desiredFrameTimeInSeconds"> The time to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::CRGB32VideoImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given time, get the next frame and write into a 
	/// CNV12VideoImg object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a NV12 eVideoFormat.
	/// <param name="desiredFrameTimeInSeconds"> The time to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::CNV12VideoImg& img, 
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given time, get the next frame and write into a 
	/// CRGBAImg based IImageWriter object </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// <param name="desiredFrameTimeInSeconds"> The time to seek to </param>
	/// <param name="img"> The image to write the result into </param>
	/// <param name="pActualFrameTimeInSeconds"> The actual time of the frame transfered </param>
	/// OpenFile() must have been called with a RGB32 or NV12 eVideoFormat.  If
    /// NV12 format is used, the video decoder will generate NV12 and a fast
    /// software conversion to RGB will be invoked for the returned frame data.
	virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, IN OUT vt::IImageWriter* pDst, 
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;
	/// <summary> Attempt to seek to a given time, get the next frame and writes 
	/// the Y result to a CLumaByteImg based IImageWriter, and the UV result to a
    /// CUVByteImg based IIImageWriter. </summary>
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a NV12 eVideoFormat.
	virtual HRESULT GetFrame(IN double desiredFrameTimeInSeconds, 
                             IN OUT vt::IImageWriter* pDstY, IN OUT vt::IImageWriter* pDstUV,
                             OUT double* pActualFrameTimeInSeconds = NULL) = 0;

	/// <summary> Add 1 to the reference counter </summary>
	virtual ULONG AddRef() = 0;
	/// <summary> Subtract 1 from the reference counter. 
	/// If the reference count is now 0, delete the object </summary>
	virtual ULONG Release() = 0;
};

/// <summary> Interface to open a video file and write images. The current 
/// implementaiton of this interface supports writing .mp4/.mov (H.264 encoded),
/// .wmv (VC1 encoded), and .3gp (H.264 encoded) videos. </summary>	
class IVideoDst
{
public:
	/// <summary> Opens a video file for writing. Sets up to be able to write
	/// the pixels from a type of object specified by the eVideoFormat. </summary>	
    /// <param name="pwcFilename"> Video filename. The extension determines the
	/// type of video format and encoding. The H.264 encoder is used for files
	/// with an extension of ".mp4", ".mov", or ".3gp".  The VC1 encoder is used
	/// for files with an extension of ".wmv". </param>
    /// <param name="iValidPixelWidth"> The width of the valid region of pixels.
	/// Media Foundation only supported even widths and heights. So iValidPixelWidth
	/// is rounded up to the nearest even number. When using CRGB32VideoImgs or
	/// CNV12VideoImg, the width should be the width of the rectValidPixels 
	/// rectangle in the CVideoImgInfo. </param>
    /// <param name="iValidPixelHeight"> The height of the valid region of pixels.
	/// Media Foundation only supported even widths and heights. So iValidPixelHeight
	/// is rounded up to the nearest even number. When using CRGB32VideoImgs or
	/// CNV12VideoImg, the height should be the height of the rectValidPixels 
	/// rectangle in the CVideoImgInfo. </param>
    /// <param name="eVideoFormat"> The format that the pixels will be written.
	/// CIMG = Want to write pixels in a CImg. RGB32 = Want to write pixels in
	/// a CRGB32VideoImg object. NV12 = Want to write pixels in a CNV12VideoImg
	/// object.
	/// See http://msdn.microsoft.com/en-us/library/ms694868(v=VS.85).aspx. </param>
    /// <param name="iFramesPerSecond"> The desired frame rate in frames per second. </param>
    /// <param name="fBitsPerPixel"> The desired bitrate in bits per pixel. </param>
    /// <param name="dAspectRatio"> The aspect ratio to write to Media Foundation. 
	/// See http://msdn.microsoft.com/en-us/library/ms704767(v=VS.85).aspx. </param>
    /// <param name="iInterlaceMode"> The interlace mode to write to Media Foundation. </param>
	virtual HRESULT OpenFile(
        __in_z const wchar_t *pwcFilename, 
        int iValidPixelWidth, int iValidPixelHeight, VideoFormat eVideoFormat,  
		int iFramesPerSecond = 25, float fBitsPerPixel = 4.0f,
		double dAspectRatio = 1.0, int iInterlaceMode = MFVideoInterlace_Progressive) = 0;

#if defined(MAKE_DOXYGEN_WORK) || defined(VT_WINRT)
	/// <summary> Opens a video file for writing. Sets up to be able to write
	/// the pixels from a type of object specified by the eVideoFormat. </summary>	
    /// <param name="storageFile"> Video file. The file's extension determines the
	/// type of video format and encoding. The H.264 encoder is used for files
	/// with an extension of ".mp4", ".mov", or ".3gp".  The VC1 encoder is used
	/// for files with an extension of ".wmv". </param>
    /// <param name="iValidPixelWidth"> The width of the valid region of pixels.
	/// Media Foundation only supported even widths and heights. So iValidPixelWidth
	/// is rounded up to the nearest even number. When using CRGB32VideoImgs or
	/// CNV12VideoImg, the width should be the width of the rectValidPixels 
	/// rectangle in the CVideoImgInfo. </param>
    /// <param name="iValidPixelHeight"> The height of the valid region of pixels.
	/// Media Foundation only supported even widths and heights. So iValidPixelHeight
	/// is rounded up to the nearest even number. When using CRGB32VideoImgs or
	/// CNV12VideoImg, the height should be the height of the rectValidPixels 
	/// rectangle in the CVideoImgInfo. </param>
    /// <param name="eVideoFormat"> The format that the pixels will be written.
	/// CIMG = Want to write pixels in a CImg. RGB32 = Want to write pixels in
	/// a CRGB32VideoImg object. NV12 = Want to write pixels in a CNV12VideoImg
	/// object.
	/// See http://msdn.microsoft.com/en-us/library/ms694868(v=VS.85).aspx. </param>
    /// <param name="iFramesPerSecond"> The desired frame rate in frames per second </param>
    /// <param name="fBitsPerPixel"> The desired bitrate in bits per pixel. </param>
    /// <param name="dAspectRatio"> The aspect ratio to write to Media Foundation. 
	/// See http://msdn.microsoft.com/en-us/library/ms704767(v=VS.85).aspx. </param>
    /// <param name="iInterlaceMode"> The interlace mode to write to Media Foundation. </param>
	/// <dl class="section remarks"><dt>Remarks</dt><dd>
	///	Available only in WinRT configurations.
	/// </dd></dl>
	virtual HRESULT OpenFile(
		Windows::Storage::IStorageFile^ storageFile, int iValidPixelWidth, int iValidPixelHeight,
		VideoFormat eVideoFormat, int iFramesPerSecond = 25, float fBitsPerPixel = 4.0f,
		double dAspectRatio = 1.0, int iInterlaceMode = MFVideoInterlace_Progressive) = 0;
#endif

	/// <summary> Writes all data and closes the video file or stream. </summary>	
	virtual HRESULT Close() = 0;

	/// <summary> Return a clone of this video dest that references the same opened
    /// file (if a file is open). </summary>	
	virtual HRESULT Clone(IVideoDst** ppClone) = 0;

	/// <summary> Writes a CImg to the video. </summary>	
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a CImg eVideoFormat. The width/height 
	/// of the CImg must match the values of iValidPixelWidth/iValidPixelHeight
	/// when OpenFile() was called
	virtual HRESULT WriteFrame(CImg &img) = 0;
	/// <summary> Writes a CRGB32VideoImg to the video. </summary>	
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a RGB32 eVideoFormat. WriteFrame()
	/// attempts to write the pixels in the rectValidPixels rectangle. If the
	/// width/height of that rectangle is too big/small, WriteFrame() attempts to
	/// place the top-left of the written pixel rectangle as close as possible
	/// to the top-left of the rectValidPixels rectangle.
	/// Currently the XVP color converter appears to perform the wrong color
	/// conversion when writing RGB32 pixels to an H.264 encoded .mp4/.mov file.
	/// As a work around, when WriteFrame() is called with a CRGB32VideoImg, first
	/// the video image is converted to a CNV12VideoImg and then WriteFrame() is
	/// called on that CNV12VideoImg frame.
	virtual HRESULT WriteFrame(CRGB32VideoImg &img) = 0;
	/// <summary> Writes a CNV12VideoImg to the video. </summary>	
 	/// <DL><DT> Remarks: </DT></DL>
	/// OpenFile() must have been called with a NV12 eVideoFormat. WriteFrame()
	/// attempts to write the pixels in the rectValidPixels rectangle. If the
	/// width/height of that rectangle is too big/small, WriteFrame() attempts to
	/// place the top-left of the written pixel rectangle as close as possible
	/// to the top-left of the rectValidPixels rectangle.
	virtual HRESULT WriteFrame(CNV12VideoImg & img) = 0;

	/// <summary> Add 1 to the reference counter </summary>
	virtual ULONG AddRef() = 0;
	/// <summary> Subtract 1 from the reference counter. 
	/// If the reference count is now 0, delete the object </summary>
	virtual ULONG Release() = 0;
};

};
