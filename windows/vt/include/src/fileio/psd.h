//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      PSD writer
//
//  History:
//      2005/1/13-swinder
//            Created
//
//------------------------------------------------------------------------

#pragma once

#if !defined(VT_WINRT)

#include "vt_iointerfaces.h"

#include "vtcommon.h"

namespace vt {

#define PSD_SIG(a, b, c, d)         ((a) | ((b)<<8) | ((c)<<16) | ((d)<<24))
#define PSD_FILE_SIGNATURE          PSD_SIG('8','B','P','S')
#define PSD_RESOURCE_SIGNATURE      PSD_SIG('8','B','I','M')

class CPascalString {
public:
    CPascalString() { memset(m_rgch, 0, 256); }
	
    HRESULT Create(const WCHAR * pwsz);
    UInt32 SizeBytes() { if(Length()==0) return 2; else return Length() + 1; }
    UInt32 SizeBytesPad2() { return (SizeBytes() + 1) & (~1); }
    UInt32 SizeBytesPad4() { return (SizeBytes() + 3) & (~3); }
    Byte Length() { return m_rgch[0]; }
    HRESULT Write(HANDLE h);
    HRESULT WritePad2(HANDLE h);
    HRESULT WritePad4(HANDLE h);
    char *Ptr() { return m_rgch + 1; }

private:
    char m_rgch[256];
};

// ensure that data is not padded to UInt32 boundaries
#pragma pack(1)

typedef struct {
    UInt32 uiSignature;
    UInt16 uiVersion;
    Byte rgbReserved[6];
    UInt16 uiChannels;
    UInt32 uiHeight;
    UInt32 uiWidth;
    UInt16 uiDepth;
    UInt16 uiColorMode;
    UInt32 uiColorTableLength;
    UInt32 uiResourcesLength;
} PSDFileHeader;

typedef struct {
    UInt32 uiLayerMaskLength;
    UInt32 uiLayerLength;
    UInt16 uiLayerCount;
} PSDFileLayerMask;

typedef struct {
    UInt64 uiLayerMaskLength;
    UInt64 uiLayerLength;
    UInt16 uiLayerCount;
} PSBFileLayerMask;

typedef struct {
    UInt32 uiTop;
    UInt32 uiLeft;
    UInt32 uiBottom;
    UInt32 uiRight;
    UInt16 uiChannels;
} PSDFileLayerPartA;

typedef struct {
    UInt16 uiID;
    UInt32 uiChannelDataLength;
} PSDFileChannelInfo;

typedef struct {
    UInt16 uiID;
    UInt64 uiChannelDataLength;
} PSBFileChannelInfo;

typedef struct {
    UInt32 uiSignature;
    Byte rgbBlendMode[4];
    Byte bOpacity;
    Byte bClipping;
    Byte bFlags;
    Byte bFiller;
    UInt32 uiExtraDataLength;
    UInt32 uiLayerMaskLength;
    UInt32 uiBlendingRangesLength; 
} PSDFileLayerPartB;

typedef struct {
    UInt32 uiFormat;
    UInt32 uiWidth;
    UInt32 uiHeight;
    UInt32 uiWidthBytes;
    UInt32 uiSize;
    UInt32 uiCompressedSize;
    UInt16 uiBitsPerPixel;
    UInt16 uiPlanes;
} PSDThumbnail;

typedef struct {
    UInt32  uiVersion;
    Byte    bHasRealMergedData;
    UInt32  cchWriterName;
    WCHAR   wszWriterName[15];
    UInt32  cchReaderName;
    WCHAR   wszReaderName[19];
    UInt32  uiFileVersion;
} PSDVersion;

typedef struct {
    UInt32 uiHRes;
    UInt16 uiHResUnit;
    UInt16 uiWidthUnit;
    UInt32 uiVRes;
    UInt16 uiVResUnit;
    UInt16 uiHeightUnit;
} PSDResolution;

// restore default packing

#pragma pack()

class CPSDRow {
public:
    // no copy constructor so this is fragile and requires care with push_back
    CPSDRow() { m_uiSize = 0; m_pbBuffer = NULL; }
    ~CPSDRow() { if(m_pbBuffer!=NULL) delete m_pbBuffer; }
    HRESULT Alloc(UInt32 uiSize);
    void NoAlloc(UInt32 uiSize) { m_uiSize = uiSize; }
    UInt32 Size() { return m_uiSize; }
    Byte *Ptr() { return m_pbBuffer; }

private:
    UInt32 m_uiSize;
    Byte *m_pbBuffer;
};

typedef struct {
    UInt16 uiID;
    CPascalString sName;
    CPSDRow cData;
} PSDResource;

typedef struct {
    vector<UInt32> vChannelLength;
    int iRows;
    vector<CPSDRow> vRows;
} PSDRowSet;

typedef struct {
    int iWidth, iHeight;
    int iX, iY;
    wstring wsName;
    UInt16 uiChannels; // 3 for RGB, 4 for RGBA
    Byte bOpacity; // 0 - 255
    bool bVisible;
    bool bAlpha; // first channel is an alpha channel
    int iValidRows;
    int iSrcType;
    Byte rgbBlendMode[4];
    vector<LARGE_INTEGER> vFileLoc;
    vector<PSDRowSet> vRowSet;
} PSDLayer;

// psd file writer
//
// call init, enable disk cacheing if your images are large.
// max image size is 32000x32000 (per psd spec).
// images can be 16 bit or 8 bit per pixel. if they are four band, bUseAlpha controls whether the alpha plane
// is written or ignored. the file consists of zero or more layer images and exactly one composed image.
// the images can be written as one, or a few rows at a time by using the "ByRows" versions of the interface.
// the last call must be WriteFile. after this, the class cannot be used. the alphablend parameters control
// the transparency of the layer as a whole in the range 0-255. 16 bit files are not compressed.
//

// layer blend modes
#define PSD_BLEND_NORMAL        "norm"
#define PSD_BLEND_DARKEN        "dark"
#define PSD_BLEND_LIGHTEN        "lite"
#define PSD_BLEND_HUE            "hue "
#define PSD_BLEND_SATURATION    "sat "
#define PSD_BLEND_COLOR            "colr"
#define PSD_BLEND_LUMINOSITY    "lum "
#define PSD_BLEND_MULTIPLY        "mul "
#define PSD_BLEND_SCREEN        "scrn"
#define PSD_BLEND_DISSOLVE        "diss"
#define PSD_BLEND_OVERLAY        "over"
#define PSD_BLEND_HARDLIGHT        "hLit"
#define PSD_BLEND_SOFTLIGHT        "sLit"
#define PSD_BLEND_DIFFERENCE    "diff"
#define PSD_BLEND_EXCLUSION        "smud"
#define PSD_BLEND_COLORDODGE    "div "
#define PSD_BLEND_COLORBURN        "idiv"


class CPSDWriter : public IVTImageWriter
{
public:
    CPSDWriter();
    ~CPSDWriter();

    // Clone the writer
    HRESULT Clone( IVTImageWriter **ppWriter );

    // OpenFile opens the file or stream for writing
    HRESULT OpenFile(const WCHAR * pwcName);
    HRESULT OpenFile( IStream*, const WCHAR *)
        { return E_NOTIMPL; }

    // Save the image
    HRESULT SetImage(IStream*, CRect&, const CParams*)
	{ return E_NOTIMPL; }

    HRESULT SetImage(const CImg &img,
                     bool bSaveMetadata = true,
                     const CParams* pParams = NULL,
                     CTaskProgress* pProgress = NULL)
        { CImgReaderWriter<CImg> src; img.Share(src);
          return SetImage(&src, NULL, bSaveMetadata, pParams, pProgress); }

    HRESULT SetImage(IImageReader* pComp,
                     const CRect* pRect = NULL,
                     bool bSaveMetadata = true,
                     const CParams* pParams = NULL,
                      CTaskProgress* pProgress = NULL)
        { return SetImageWithLayers(pComp, pRect, bSaveMetadata, pParams, pProgress); }

    HRESULT SetImageWithLayers(IImageReader* pComp = NULL,
                               const CRect* pRect = NULL,
                               bool bSaveMetadata = true,
                               const CParams* pParams = NULL,
                               CTaskProgress* pProgress = NULL,
                               const CRGBAImg* pThumb = NULL,
                               IIndexedImageReader* pLayers = NULL,
                               vt::wstring* pwstrNames = NULL,
                               IIndexedImageReader* pMasks = NULL);

	// finish writing and close
    HRESULT CloseFile( );

protected:
    HRESULT Init(int iX, int iY, int iW, int iH, bool b16Bit = false, bool bCompress = true);

    HRESULT AddVersion(bool bComp);
    HRESULT AddResolution(const CParams &params);
    HRESULT AddProfile(const CParams &params);
    HRESULT AddThumbnail(const CRGBAImg &img);
    HRESULT AddMetaData(const CParams *pParams);
    
    // layer image must have 4 bands if alpha is enabled
    // if present, the layer mask consists of a 1 band image
    HRESULT AddLayerImageByRowsInit(int *piLayerIDRtn, const WCHAR * pwcName, 
                                    CRect rect, int iType, int iBands, 
                                    bool bUseAlpha, bool bVisible, bool bUseLayerMask = false,
                                    int iAlphaBlend = 255, 
                                    const char *pchBlendMode = NULL);
    HRESULT AddLayerImageByRows(int iLayerID, Byte *pData, int iStride, 
                                int iRows = 0, Byte *pLayerMaskData = NULL, int iLayerMaskStride = 0);
    HRESULT AddComposedImageByRowsInit(int iType, int iBands);
    HRESULT AddComposedImageByRows(Byte *pData, int iStride, int iRows = 0);

    HRESULT WriteHeader();
    HRESULT WriteLayers();

private:
    UInt16 Swap16(UInt16 usVal);
    UInt32 Swap32(UInt32 uiVal);
    UInt64 Swap64(UInt64 uiVal);
    HRESULT ProcessLayer(int iLayer, int iRowSet, int iChan,
                         const CImg &img, const CImg &imgLayerMask);
    HRESULT CompressBand(const Byte *pbStart, int iElSize,
                         int iPixSize, int iWidth, UInt32 *puiSize);
    void CopyBand(const Byte *pbSrc, Byte *pbDst, int iElSize, 
                  int iPixSize, int iWidth);
    HRESULT WriteChannelRowData(int iLayer, int iChan);
    HRESULT WriteBytes(HANDLE h, const void *pData, DWORD dwSize);

    vector<PSDResource> m_resources;
    vector<PSDLayer> m_layers; // layer 0 is the composed image
    bool m_b16Bit;
    bool m_bInit;
    bool m_bCompress;
    bool m_bPSB;
    HANDLE m_hOutput;
    Byte *m_pbBuffer;

    wstring m_wstrName;
};

};

#endif
