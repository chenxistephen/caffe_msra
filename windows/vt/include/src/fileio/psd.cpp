//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      PSD writer
//
//  History:
//      2005/1/13-swinder
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#if !defined(VT_WINRT)

#include "vt_global.h"
#include "vt_io.h"
#include "psd.h"
#include "wicio.h"

using namespace vt;

#define PSD_MAX_LAYERS      32767
#define COMPRESS_BUFFER_SIZE    655360

HRESULT CPSDRow::Alloc(UInt32 uiSize)
{
    HRESULT hr = NOERROR;

    m_pbBuffer = VT_NOTHROWNEW Byte [uiSize];
    if(m_pbBuffer==NULL)
        VT_HR_EXIT(E_OUTOFMEMORY);

    m_uiSize = uiSize;

Exit:
    return hr;
}

CPSDWriter::CPSDWriter()
: m_bInit(false), m_hOutput(INVALID_HANDLE_VALUE), m_pbBuffer(NULL)
{
}

CPSDWriter::~CPSDWriter()
{
    CloseFile();
}

HRESULT CPSDWriter::CloseFile()
{
    m_bInit = false;
    if(m_hOutput!=INVALID_HANDLE_VALUE)
        CloseHandle(m_hOutput);
    m_hOutput = INVALID_HANDLE_VALUE;
    if(m_pbBuffer)
        delete m_pbBuffer;
    m_pbBuffer = NULL;
    m_layers.clear();
	return S_OK;
}

HRESULT
CPSDWriter::Clone(IVTImageWriter **ppWriter)
{
    if (ppWriter == NULL)
        return E_INVALIDARG;

    if ((*ppWriter = VT_NOTHROWNEW CPSDWriter) == NULL)
        return E_OUTOFMEMORY;

    return S_OK;
}

HRESULT CPSDWriter::OpenFile(const WCHAR * pwcName)
{
    HRESULT hr = NOERROR;

    VT_PTR_EXIT(pwcName);

    m_bPSB = _wcsicmp(VtGetFileExt(pwcName), L".psb") == 0;

    VT_HR_EXIT( m_wstrName.assign(pwcName) );

Exit:
    return hr;
}

#define SCANLINES 16

HRESULT CPSDWriter::SetImageWithLayers(IImageReader *pComp,
                                       const CRect* pRect,
                                       bool bSaveMetadata,
                                       const CParams*,
                                       CTaskProgress *pProgress,
                                       const CRGBAImg* pThumb,
                                       IIndexedImageReader *pLayers,
                                       wstring *pwstrNames,
                                       IIndexedImageReader *pMasks)
{
    HRESULT hr = NOERROR;

    CParams params;

    LARGE_INTEGER liHead = { 0 };

	// compute some progress params
	CPhasedTaskStatus prog;
	prog.SetOuterCallback(pProgress);

	// compute the total number of bytes in image
	// this is used to compute progress and check for PSB
	UInt64 uiByteCount = 0;

    // Initialize composed image (if any).
    if (pComp != NULL)
    {
        CImgInfo info = pComp->GetImgInfo();

		int iBands  = info.Bands() == 4 ? 4 : 3;
        int iElType = EL_FORMAT(info.type) == EL_FORMAT_BYTE ?
            EL_FORMAT_BYTE : EL_FORMAT_SHORT;

        CRect rctComp = info.Rect();

        VT_HR_EXIT( pComp->GetMetaData(params) );

        VT_HR_EXIT( Init(rctComp.left, rctComp.top,
                         rctComp.Width(), rctComp.Height(),
                         iElType == EL_FORMAT_SHORT,
                         true) );

        VT_HR_EXIT( AddComposedImageByRowsInit(iElType, iBands) );

        uiByteCount += (UInt64) m_layers[0].iWidth * m_layers[0].iHeight *
            m_layers[0].uiChannels * (m_b16Bit ? 2 : 1);
    }

    // Initialize layered images (if any).
    if (pLayers != NULL)
    {
        if (pMasks != NULL && pMasks->GetFrameCount() != pLayers->GetFrameCount())
            VT_HR_EXIT( E_INVALIDARG );

        for (UINT iIndex = 0; iIndex < pLayers->GetFrameCount(); iIndex++)
        {
            CLayerImgInfo info = pLayers->GetImgInfo(iIndex);

            int iBands = info.Bands() == 4 ? 4 : 3;
            int iElType = EL_FORMAT(info.type) == EL_FORMAT_BYTE ?
                EL_FORMAT_BYTE : EL_FORMAT_SHORT;

            CParams layers;
            VT_HR_EXIT( pLayers->GetMetaData(iIndex, layers) );
            VT_HR_EXIT( params.Merge(&layers) );

            if (!m_bInit)
            {
                VT_ASSERT( pComp == NULL );

				vt::CRect rctComp = pRect? *pRect: 
					vt::CRect(0,0,info.compositeWidth, info.compositeHeight);

                // If no composed image make bogus one of composite size.
                VT_HR_EXIT( Init(0, 0, rctComp.Width(), rctComp.Height(),
                                  iElType == EL_FORMAT_SHORT, true) );

                VT_HR_EXIT( AddComposedImageByRowsInit(iElType, 4) );

                // Bogus black composed image compresses well with RLE so don't count it.
                // uiByteCount += (UInt64) m_layers[0].iWidth * m_layers[0].iHeight *
                //     m_layers[0].uiChannels * (m_b16Bit ? 2 : 1);
            }

			vt::CRect rctLayer = info.LayerRect();
			if (pRect != NULL)
				rctLayer -= pRect->TopLeft();

            int iLayerID;
            VT_HR_EXIT( AddLayerImageByRowsInit(
				&iLayerID, pwstrNames == NULL? NULL : pwstrNames[iIndex],
				rctLayer, iElType, iBands, iBands == 4, true, pMasks != NULL, 
				255) );

            uiByteCount += (UInt64) m_layers[iLayerID].iWidth * m_layers[iLayerID].iHeight *
                m_layers[iLayerID].uiChannels * (m_b16Bit ? 2 : 1);
        }
    }

    // TBD:  Switch to PSB if size too big for PSD.
    // Ideally, compute exact PSD size and start whole process over with PSB if >= 2 GB (2147483648).
    // Instead, compute uncompressed image size and swtich to PSB if >= 2000000000
    // (conservative to account for headers, metadata, bad RLE compression, etc.)
    if (uiByteCount >= 2000000000)
    {
        m_bPSB = true;
        LPWSTR pwszExt = (LPWSTR) VtGetFileExt(m_wstrName);
        if (pwszExt != NULL && _wcsicmp(pwszExt, L".psd") == 0)
            pwszExt[3] = iswlower(pwszExt[3]) ? L'b' : L'B';
    }

    UInt64 uiLayerLength = 0;
    if (pLayers != NULL)
    {
        for (UINT iIndex = 0; iIndex < pLayers->GetFrameCount(); iIndex++)
        {
            int iLayerID = iIndex + 1;

            uiLayerLength += sizeof(PSDFileLayerPartA);
            uiLayerLength += (m_bPSB ? sizeof(PSBFileChannelInfo) : sizeof(PSDFileChannelInfo)) *
                m_layers[iLayerID].uiChannels;
            uiLayerLength += sizeof(PSDFileLayerPartB);
            if (pwstrNames != NULL)
            {
                CPascalString sLayerName;
                VT_HR_EXIT( sLayerName.Create((const WCHAR *) m_layers[iLayerID].wsName) );
                uiLayerLength += sLayerName.SizeBytesPad4();
            }
        }

        uiLayerLength += (m_bPSB ? sizeof(PSBFileLayerMask) : sizeof(PSDFileLayerMask));
    }

    // Initialize metadata.
    if (bSaveMetadata)
        VT_HR_EXIT( AddMetaData(&params) );

    // Initialize thumbnail (if any).
    if (pThumb != NULL)
        VT_HR_EXIT( AddThumbnail(*pThumb) );

    // Initialize version.
    VT_HR_EXIT( AddVersion(pComp != NULL) );

    // Initialize resolution.
    VT_HR_EXIT( AddResolution(params) );

    // Initialize color profile (if any).
    VT_HR_EXIT( AddProfile(params) );

    // Write header and resources (thumbnail, metadata) to file.
    VT_HR_EXIT( WriteHeader() );

    // Skip over layer headers.
    liHead.QuadPart = uiLayerLength;
    if (!SetFilePointerEx(m_hOutput, liHead, &liHead, FILE_CURRENT))
        VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));

	// setup progress
	float fProgressScale = 100.f / (uiByteCount == 0 ? 1.f : (float) uiByteCount);

    // Write layered images (if any).
    if (pLayers != NULL)
    {
        for (UINT iIndex = 0; iIndex < pLayers->GetFrameCount(); iIndex++)
        {
            CImg imgLayer, imgMask;

            int iLayerID = iIndex + 1;

            vt::string_b<64> phaseName;
            if (!m_layers[iLayerID].wsName.empty())
                phaseName.format("Saving layer %S", m_layers[iLayerID].wsName.get_constbuffer());
            else
                phaseName.format("Saving layer %d", iLayerID);
            prog.BeginPhase(phaseName, fProgressScale *
                (float) m_layers[iLayerID].iWidth * (float) m_layers[iLayerID].iHeight *
                (float) m_layers[iLayerID].uiChannels * (m_b16Bit ? 2.f : 1.f));

            for (int iChan = 0; iChan < m_layers[iLayerID].uiChannels; iChan++)
            {
                UInt16 uiCompMode = Swap16(m_bCompress ? 1 : 0);
                VT_HR_EXIT( WriteBytes(m_hOutput, &uiCompMode, sizeof(uiCompMode)) );

                // copy 16 scan-lines at a time
                int iRowSet = 0;
                for (CBlockIterator bi(BLOCKITER_INIT(CRect(0, 0, m_layers[iLayerID].iWidth, m_layers[iLayerID].iHeight),
                                                      m_layers[iLayerID].iWidth, SCANLINES));
                    !bi.Done(); bi.Advance(), iRowSet++)
                {
                    const CRect rctSubLayer = bi.GetRect();

                    VT_HR_EXIT( CreateImageForTransform(imgLayer, rctSubLayer.Width(), rctSubLayer.Height(),
                                                        m_layers[iLayerID].iSrcType) );

                    VT_HR_EXIT( pLayers->ReadRegion(iIndex, rctSubLayer, imgLayer) );

                    if (pMasks != NULL)
                    {
                        CImgInfo mask = pMasks->GetImgInfo(iIndex);

                        VT_HR_EXIT( CreateImageForTransform(imgMask, imgLayer.Width(), imgLayer.Height(),
                                                            mask.type) );

                        VT_HR_EXIT( pMasks->ReadRegion(iIndex, rctSubLayer, imgMask) );
                    }

                    if (iChan == 0)
                        VT_HR_EXIT( AddLayerImageByRows(iLayerID, imgLayer.BytePtr(), imgLayer.StrideBytes(),
                                                         imgLayer.Height(),
                                                         imgMask.BytePtr(), imgMask.StrideBytes()) );
                    else
                        VT_HR_EXIT( ProcessLayer(iLayerID, iRowSet, iChan, imgLayer, imgMask) );
                }

                if (!FlushFileBuffers(m_hOutput))
                    VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

                prog.ReportProgress(100.f * (float) (iChan + 1) / (float) m_layers[iLayerID].uiChannels);
            }
        }

        if (m_bCompress)
        {
            for(int iLayer=1; iLayer<(int)m_layers.size(); iLayer++)
            {
                for(int iChan = 0; iChan<m_layers[iLayer].uiChannels; iChan++)
                {
                    // Go back to channel row table before first row set.
                    LARGE_INTEGER liZero = m_layers[iLayer].vFileLoc[iChan];
                    liZero.QuadPart -=
                        (m_bPSB ? sizeof(UInt32) : sizeof(UInt16)) *
                            m_layers[iLayer].iHeight;
 
                    if (!SetFilePointerEx(m_hOutput, liZero, NULL, FILE_BEGIN))
                        VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));

                    VT_HR_EXIT( WriteChannelRowData(iLayer, iChan) );
                }
            }
        }
    }

    // Go back to layer headers.
    liHead.QuadPart -= uiLayerLength;
    if (!SetFilePointerEx(m_hOutput, liHead, NULL, FILE_BEGIN))
        VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));

    // Write layer headers to file.
    VT_HR_EXIT( WriteLayers() );

    // Write composed image.
    if (!m_layers.empty())
    {
        CImg imgComp;

        CRect rctComp = CRect(m_layers[0].iX, m_layers[0].iY, m_layers[0].iWidth, m_layers[0].iHeight);

        if (pComp != NULL)
            prog.BeginPhase("Saving composed image", fProgressScale *
                (float) m_layers[0].iWidth * (float) m_layers[0].iHeight *
                (float) m_layers[0].uiChannels * (m_b16Bit ? 2.f : 1.f));

        UInt16 uiCompMode = Swap16(m_bCompress ? 1 : 0);
        VT_HR_EXIT( WriteBytes(m_hOutput, &uiCompMode, sizeof(uiCompMode)) );

        for (int iChan = 0; iChan < m_layers[0].uiChannels; iChan++)
        {
            // copy 16 scan-lines at a time
            int iRowSet = 0;
            for (CBlockIterator bi(BLOCKITER_INIT(rctComp, rctComp.Width(),
                                                  pComp != NULL ? SCANLINES : rctComp.Height()));
                !bi.Done(); bi.Advance(), iRowSet++)
            {
                const CRect rctSubComp = bi.GetCompRect();

                if (pComp != NULL)
                {
                    VT_HR_EXIT( CreateImageForTransform(imgComp, rctSubComp.Width(), rctSubComp.Height(),
                                                         m_layers[0].iSrcType) );

                    VT_HR_EXIT( pComp->ReadRegion(rctSubComp, imgComp) );
                }

                if (iChan == 0)
                    VT_HR_EXIT( AddComposedImageByRows(imgComp.BytePtr(), imgComp.StrideBytes(),
                                                        imgComp.Height()) );
                else
                    VT_HR_EXIT( ProcessLayer(0, iRowSet, iChan, imgComp, CImg()) );
            }

            if (!FlushFileBuffers(m_hOutput))
                VT_HR_EXIT( HRESULT_FROM_WIN32( GetLastError() ) );

            if (pComp != NULL)
                prog.ReportProgress(100.f * (float) (iChan + 1) / (float) m_layers[0].uiChannels);
        }

        if (m_bCompress)
        {
            // Go back to channel row table before first row set.
            // determine current file pointer
            LARGE_INTEGER liZero = m_layers[0].vFileLoc[0];
            liZero.QuadPart -=
                (m_bPSB ? sizeof(UInt32) : sizeof(UInt16)) *
                    m_layers[0].uiChannels * m_layers[0].iHeight;

            if (!SetFilePointerEx(m_hOutput, liZero, NULL, FILE_BEGIN))
                VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));

            // composed image data
            for(int iChan = 0; iChan<m_layers[0].uiChannels; iChan++)
                VT_HR_EXIT( WriteChannelRowData(0, iChan) );
        }
    }

Exit:
    return hr;
}

HRESULT CPSDWriter::Init(int iX, int iY, int iW, int iH, bool b16Bit, bool bCompress)
{
    PSDLayer psd;

    HRESULT hr = NOERROR;

    // Supported range is 1 to 30,000. (**PSB** max of 300,000.)
    if (iW < 1 || iW > 300000 || iH < 1 || iH > 300000)
        VT_HR_EXIT(E_INVALIDARG);

    if (m_bInit)
        VT_HR_EXIT( E_FAIL );

    m_bCompress = false;
    if(bCompress)
        m_bCompress = true;
    
    {
        // for row compression and uncompressed re-order buffer for write ops
        m_pbBuffer = VT_NOTHROWNEW Byte [COMPRESS_BUFFER_SIZE];
        if(m_pbBuffer==NULL)
            VT_HR_EXIT(E_OUTOFMEMORY);
    }

    psd.iX = iX;
    psd.iY = iY;
    psd.iWidth = iW;
    psd.iHeight = iH;
    psd.iValidRows = 0;
    psd.uiChannels = 3;
    psd.bAlpha = false;
   
    VT_HR_EXIT( m_layers.push_back(psd) );

    m_b16Bit = b16Bit;
    m_bInit = true;

Exit:
    return hr;
}

HRESULT CPSDWriter::AddVersion(bool bComp)
{
    HRESULT hr = NOERROR;
    PSDResource res;
    PSDVersion version;

    version.uiVersion = Swap32(m_bPSB ? 2 : 1); // version 1 (**PSB** version is 2.)
    version.bHasRealMergedData = (bComp ? 1 : 0);
    version.cchWriterName = Swap32(15);
    for (int i = 0; i < 15; i++)
        version.wszWriterName[i] = Swap16(L"Adobe Photoshop"[i]);
    version.cchReaderName = Swap32(19);
    for (int i = 0; i < 19; i++)
        version.wszReaderName[i] = Swap16(L"Adobe Photoshop CS3"[i]);
    version.uiFileVersion = Swap32(1);

    VT_HR_EXIT( m_resources.push_back(res) );

    m_resources.back().uiID = 0x421;

    VT_HR_EXIT( m_resources.back().cData.Alloc(sizeof(version)) );

    Byte *pbData = m_resources.back().cData.Ptr();
    memcpy(pbData, &version, sizeof(version));

Exit:
    return hr;
}

HRESULT CPSDWriter::AddResolution(const CParams &params)
{
    HRESULT hr = NOERROR;
    PSDResource res;
    PSDResolution resolution;

    const CParamValue *value = NULL;

    // Get the X resolution (if any).
    UInt32 uiXRes = 0;
    // sometimes rationals are stored as U8s
    if (params.GetById(&value, ImagePropertyXResolution) == S_OK &&
        value != NULL && value->GetDataPtr() != NULL &&
        (value->GetType() == ParamType_Rational ||
         value->GetType() == ParamType_U8) &&
        value->GetDataSize() >= sizeof(UInt64))
        uiXRes = (UInt32)
            ((*((RATIONAL *) value->GetDataPtr())).AsFloat() * 65536.f);

    // Get the Y resolution (if any).
    UInt32 uiYRes = 0;
    // sometimes rationals are stored as U8s
    if (params.GetById(&value, ImagePropertyYResolution) == S_OK &&
        value != NULL && value->GetDataPtr() != NULL &&
        (value->GetType() == ParamType_Rational ||
         value->GetType() == ParamType_U8) &&
        value->GetDataSize() >= sizeof(UInt64))
        uiYRes = (UInt32)
            ((*((RATIONAL *) value->GetDataPtr())).AsFloat() * 65536.f);

    // Get the resolution unit (if any).
    UInt16 uiUnit = 0;  // 2 = in, 3 = cm
    if (params.GetById(&value, ImagePropertyResolutionUnit) == S_OK &&
        value != NULL && value->GetType() == ParamType_UShort &&
        value->GetDataPtr() != NULL &&
        value->GetDataSize() >= sizeof(UInt16))
        uiUnit = *((UInt16 *) value->GetDataPtr());

    if (uiXRes != 0 && uiYRes != 0 && uiUnit != 0)
    {
        resolution.uiHRes = Swap32(uiXRes);
        resolution.uiHResUnit = Swap16(uiUnit - 1);
        resolution.uiWidthUnit = Swap16(uiUnit - 1);
        resolution.uiVRes = Swap32(uiYRes);
        resolution.uiVResUnit = Swap16(uiUnit - 1);
        resolution.uiHeightUnit = Swap16(uiUnit - 1);

        VT_HR_EXIT( m_resources.push_back(res) );

        m_resources.back().uiID = 0x3ed;

        VT_HR_EXIT( m_resources.back().cData.Alloc(sizeof(resolution)) );

        Byte *pbData = m_resources.back().cData.Ptr();
        memcpy(pbData, &resolution, sizeof(resolution));
    }

Exit:
    return hr;
}

HRESULT CPSDWriter::AddProfile(const CParams &params)
{
    HRESULT hr = NOERROR;
    PSDResource res;

    const CParamValue *value = NULL;

    // Get the ICC profile (if any).
    if (params.GetById(&value, ImagePropertyICCProfile) == S_OK &&
        value != NULL && value->GetDataPtr() != NULL &&
        value->GetDataSize() > 0)
    {
        VT_HR_EXIT( m_resources.push_back(res) );

        m_resources.back().uiID = 0x40f;

        VT_HR_EXIT( m_resources.back().cData.Alloc((UInt32) value->GetDataSize()) );

        Byte *pbData = m_resources.back().cData.Ptr();
        memcpy(pbData, value->GetDataPtr(), value->GetDataSize());
    }

Exit:
    return hr;
}

HRESULT CPSDWriter::AddThumbnail(const CRGBAImg &img)
{
    PSDResource res;
    CMemStream cMemStream;
    PSDThumbnail thumb;
    CWicWriter wr;

	VT_HR_BEGIN()

    VT_HR_EXIT( cMemStream.Create() );

    if(!img.IsValid())
        VT_HR_EXIT(E_INVALIDARG);

    VT_HR_EXIT( wr.OpenFile(&cMemStream, L".jpg") );
    VT_HR_EXIT( wr.SetImage(img, false) );
	wr.CloseFile();  // need to do this here to ensure that stream is flushed

    UInt32 uiCompressedSize = (UInt32)cMemStream.Size();
    if(uiCompressedSize==0)
        VT_HR_EXIT( E_FAIL );

    UInt32 uiWidthBytes = ((img.Width() * 24 + 31)/32) * 4;
    UInt32 uiTotalSize = uiWidthBytes * img.Height();

    thumb.uiFormat = Swap32(1); //  kJpegRGB
    thumb.uiWidth = Swap32(img.Width());
    thumb.uiHeight = Swap32(img.Height());
    thumb.uiWidthBytes = Swap32(uiWidthBytes);
    thumb.uiSize = Swap32(uiTotalSize);
    thumb.uiCompressedSize = Swap32(uiCompressedSize);
    thumb.uiBitsPerPixel = Swap16(24);
    thumb.uiPlanes = Swap16(1);

    VT_HR_EXIT( m_resources.push_back(res) );

    m_resources.back().uiID = 0x40c;

    VT_HR_EXIT( m_resources.back().cData.Alloc(sizeof(thumb) + uiCompressedSize) );

    Byte *pbData = m_resources.back().cData.Ptr();
    memcpy(pbData, &thumb, sizeof(thumb));
    pbData += sizeof(thumb);
    
    LARGE_INTEGER iZero;
    iZero.QuadPart = 0;
    VT_HR_EXIT( cMemStream.Seek(iZero, STREAM_SEEK_SET, NULL) );
    VT_HR_EXIT( cMemStream.Read(pbData, uiCompressedSize, NULL) );

    VT_HR_END()
}

HRESULT CPSDWriter::AddLayerImageByRowsInit(int *piLayerIDRtn, const WCHAR * pwcName,
                                            CRect rct, int iType, int iBands, bool bUseAlpha,
                                bool bVisible, bool bUseLayerMask, int iAlphaBlend, const char *pchBlendMode)
{
    PSDLayer psd;

    HRESULT hr = NOERROR;
    if(!m_bInit)
        VT_HR_EXIT(E_NOINIT);
    if(m_layers.size() > PSD_MAX_LAYERS)
        VT_HR_EXIT(E_FAIL);
    if(piLayerIDRtn==NULL)
        VT_HR_EXIT(E_POINTER);

    //if(m_layers.size()>1 && m_layers.back().iHeight != m_layers.back().iValidRows)
    //    VT_HR_EXIT(E_FAIL); // didnt finish the previous layer

    iType = EL_FORMAT(iType);

    if(pchBlendMode==NULL || strlen(pchBlendMode)<4)
        pchBlendMode = PSD_BLEND_NORMAL;

    psd.iX = rct.left;
    psd.iY = rct.top;
    psd.iWidth = rct.Width();
    psd.iHeight = rct.Height();
    psd.iValidRows = 0;
    psd.uiChannels = 3;
    psd.bAlpha = false;
    psd.bOpacity = (Byte)iAlphaBlend;
    psd.bVisible = bVisible;
    psd.iSrcType = VT_IMG_MAKE_TYPE(iType, iBands);
    memcpy(psd.rgbBlendMode, pchBlendMode, 4);

    // only allow 3 or 4 banded byte or short images
    if(iBands!=3 && iBands!=4)
        VT_HR_EXIT(E_INVALIDARG);
    if(iBands==4 && bUseAlpha)
    {
        if(bUseLayerMask)
            psd.uiChannels = 5;
        else
            psd.uiChannels = 4;

        psd.bAlpha = true;
    }
    else
        psd.bOpacity = 255;

    if(iType==EL_FORMAT_BYTE)
    {
        if(m_b16Bit)
            VT_HR_EXIT(E_INVALIDARG);
    }
    else if(iType==EL_FORMAT_SHORT)
    {
        if(!m_b16Bit)
            VT_HR_EXIT(E_INVALIDARG);
    }
    else
        VT_HR_EXIT(E_INVALIDARG);

    VT_HR_EXIT( m_layers.push_back(psd) );

    int iLayer = (int)(m_layers.size() - 1);

    if(pwcName!=NULL)
        VT_HR_EXIT( m_layers[iLayer].wsName.assign(pwcName) );

    *piLayerIDRtn = iLayer;

    // Switch to PSB if width or height too big for PSD.
    if (m_layers[iLayer].iWidth > 30000 || m_layers[iLayer].iHeight > 30000)
    {
        m_bPSB = true;
        LPWSTR pwszExt = (LPWSTR) VtGetFileExt(m_wstrName);
        if (pwszExt != NULL && _wcsicmp(pwszExt, L".psd") == 0)
            pwszExt[3] = iswlower(pwszExt[3]) ? L'b' : L'B';
    }

Exit:
    return hr;
}

HRESULT CPSDWriter::AddLayerImageByRows(int iLayer, Byte *pData, int iStride, int iRows,
                                        Byte *pLayerMaskData, int iLayerMaskStride)
{
    CImg img, imgLayerMask;
    PSDRowSet rowset;
    HRESULT hr = NOERROR;

    if(!m_bInit)
        VT_HR_EXIT(E_NOINIT);

    if(iLayer<1 || iLayer > (int)m_layers.size() - 1)
        VT_HR_EXIT(E_INVALIDARG);

    if(iRows<1)
        iRows = m_layers[iLayer].iHeight - m_layers[iLayer].iValidRows; // write remainder

    if(iRows<1 || m_layers[iLayer].iValidRows + iRows > m_layers[iLayer].iHeight)
        VT_HR_EXIT(E_INVALIDARG);

    VT_HR_EXIT( m_layers[iLayer].vFileLoc.resize(m_layers[iLayer].uiChannels) );

    // create new row set
    VT_HR_EXIT( m_layers[iLayer].vRowSet.push_back(rowset) );

    // create a set of row buffers for all the channels
    VT_HR_EXIT( m_layers[iLayer].vRowSet.back().vRows.resize(iRows * m_layers[iLayer].uiChannels) );

    m_layers[iLayer].vRowSet.back().iRows = iRows;

    int iSrcType = m_layers[iLayer].iSrcType;
    // create image to wrap the data
    VT_HR_EXIT( img.Create(pData, m_layers[iLayer].iWidth, iRows, iStride, iSrcType) );

    if(m_layers[iLayer].uiChannels > 4)
        VT_HR_EXIT( imgLayerMask.Create(pLayerMaskData, m_layers[iLayer].iWidth, iRows,
            iLayerMaskStride, VT_IMG_MAKE_TYPE(EL_FORMAT(iSrcType), 1)) );

    VT_HR_EXIT( ProcessLayer(iLayer, (int)(m_layers[iLayer].vRowSet.size() - 1), 0,
                              img, imgLayerMask) );

    m_layers[iLayer].iValidRows += img.Height();

Exit:
    return hr;
}

HRESULT CPSDWriter::AddComposedImageByRowsInit(int iType, int iBands)
{
    HRESULT hr = NOERROR;
    if(!m_bInit)
        VT_HR_EXIT(E_NOINIT);
    if(m_layers[0].iValidRows > 0)
        VT_HR_EXIT(E_FAIL);

    // only allow 3 or 4 banded byte or short images
    if(iBands!=3 && iBands!=4)
        VT_HR_EXIT(E_INVALIDARG);
    if(iBands==4)
    {
        m_layers[0].uiChannels = 4;
        m_layers[0].bAlpha = true;
    }

    iType = EL_FORMAT(iType);

    m_layers[0].iSrcType = VT_IMG_MAKE_TYPE(iType, iBands);

    if(iType==EL_FORMAT_BYTE)
    {
        if(m_b16Bit)
            VT_HR_EXIT(E_INVALIDARG);
    }
    else if(iType==EL_FORMAT_SHORT)
    {
        if(!m_b16Bit)
            VT_HR_EXIT(E_INVALIDARG);
    }
    else
        VT_HR_EXIT(E_INVALIDARG);

    // Switch to PSB if width or height too big for PSD.
    if (m_layers[0].iWidth > 30000 || m_layers[0].iHeight > 30000)
    {
        m_bPSB = true;
        LPWSTR pwszExt = (LPWSTR) VtGetFileExt(m_wstrName);
        if (pwszExt != NULL && _wcsicmp(pwszExt, L".psd") == 0)
            pwszExt[3] = iswlower(pwszExt[3]) ? L'b' : L'B';
    }

Exit:
    return hr;
}

HRESULT CPSDWriter::AddComposedImageByRows(Byte *pData, int iStride, int iRows)
{
    CImg img, imgLayerMask;
    PSDRowSet rowset;
    HRESULT hr = NOERROR;

    if(!m_bInit)
        VT_HR_EXIT(E_NOINIT);

    if(iRows<1)
        iRows = m_layers[0].iHeight - m_layers[0].iValidRows; // write remainder

    if(iRows<1 || m_layers[0].iValidRows + iRows > m_layers[0].iHeight)
        VT_HR_EXIT(E_INVALIDARG);

    VT_HR_EXIT( m_layers[0].vFileLoc.resize(1) );

    // create new row set
    VT_HR_EXIT( m_layers[0].vRowSet.push_back(rowset) );

    // create a set of row buffers for all the channels
    VT_HR_EXIT( m_layers[0].vRowSet.back().vRows.resize(iRows * m_layers[0].uiChannels) );

    m_layers[0].vRowSet.back().iRows = iRows;

    // create image to wrap the data
    if (pData != NULL)
        VT_HR_EXIT( img.Create(pData, m_layers[0].iWidth, iRows,
                               iStride, m_layers[0].iSrcType) );

    VT_HR_EXIT( ProcessLayer(0, (int)(m_layers[0].vRowSet.size() - 1), 0,
                              img, imgLayerMask) );

    m_layers[0].iValidRows += img.Height();

Exit:
    return hr;
}

HRESULT CPSDWriter::ProcessLayer(int iLayer, int iRowSet, int c, const CImg &img, const CImg &imgLayerMask)
{
    HRESULT hr = NOERROR;
    PSDLayer *pL = &m_layers[iLayer];
    PSDRowSet *pRS = &pL->vRowSet[iRowSet];

    {
        // determine current file pointer
        LARGE_INTEGER liZero;
        liZero.QuadPart = 0;

        if (iRowSet == 0 && (c == 0 || iLayer > 0))
        {
            // add size of row lengths table if we are compressing
            if(m_bCompress)
            {
                liZero.QuadPart += (m_bPSB ? sizeof(UInt32) : sizeof(UInt16)) *
                    pL->iHeight * (iLayer == 0 ? pL->uiChannels : 1);
            }

            if (!SetFilePointerEx(m_hOutput, liZero, &pL->vFileLoc[c], FILE_CURRENT))
                VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));
        }

        // write image planes to disk in RGB order for 3 channels, ARGB order for 4, and ARGB-Layermask order for 5
        int iHeight = (img.IsValid() ? img.Height() : pL->iHeight);
        int j=c * iHeight;
        {
            int index = c;
            if (iLayer > 0)
                index = 3-c;    // 3, 2, 1, 0, -1
            else if (c < 3)
                index = 2-c;    // 2, 1, 0, 3

            UInt32 uiTmp = 0;
            VT_HR_EXIT( pRS->vChannelLength.push_back(uiTmp) );

            int i;
            const Byte *pbStart;
            int iInc;
            int iPixSize;
            int iElSize;

            if(index<0)
            {
                pbStart = imgLayerMask.BytePtr();
                iInc = imgLayerMask.StrideBytes();
                iPixSize = imgLayerMask.PixSize();
                iElSize = imgLayerMask.ElSize();
            }
            else if (img.IsValid())
            {
                pbStart  = img.BytePtr() + index * (int)img.ElSize();
                iInc     = img.StrideBytes();
                iPixSize = img.PixSize();
                iElSize  = img.ElSize();
            }
            else
            {
				iElSize  = VT_IMG_ELSIZE(pL->iSrcType);
                iPixSize = iElSize;
                iInc = 0;
                pbStart = (Byte *) calloc(iPixSize, pL->iWidth);
            }
            VT_PTR_EXIT( pbStart );

            UInt32 uiSize = 0;
            for(i=0; i<iHeight; i++, j++, pbStart += iInc)
            {
                if(m_bCompress && (img.IsValid() || uiSize == 0))
                {
                    // 8 bit or 16 bit compressed
                    // compress band to temporary buffer
                    VT_HR_EXIT( CompressBand(pbStart, iElSize, iPixSize, pL->iWidth, &uiSize) );
                }
                else if (img.IsValid() || uiSize == 0)
                {
                    // 8 bit or 16 bit uncompressed
                    uiSize = (UInt32)(pL->iWidth * iElSize);
                    // copy directly uncompressed to buffer
                    CopyBand(pbStart, m_pbBuffer, iElSize, iPixSize, pL->iWidth);
                }

                // save size in structure but don't alloc a row buffer
                pRS->vRows[j].NoAlloc(uiSize);

                // write data to output file
                VT_HR_EXIT( WriteBytes(m_hOutput, m_pbBuffer, uiSize) );

                pRS->vChannelLength.back() += uiSize;
            }

            if (index >= 0 && !img.IsValid())
                free(const_cast<Byte*>(pbStart));
        }
    }

Exit:
    return hr;
}

HRESULT CPSDWriter::AddMetaData(const CParams *pParams)
{
    VT_HR_BEGIN()

    CMemStream stream;

    PSDResource res;
    res.uiID = 0x422;

    LARGE_INTEGER liZero = { 0 };

    if(!m_bInit)
        VT_HR_EXIT(E_NOINIT);

    // copy the metadata via WIC
    CWicMetadataWriter wr;
    VT_HR_EXIT( wr.OpenStream(&stream) );

    VT_HR_EXIT( wr.SetMetadata(pParams) );

    wr.CloseStream();

    VT_HR_EXIT( stream.Seek(liZero, STREAM_SEEK_SET, NULL) );

    STATSTG statstg;
    ZeroMemory(&statstg, sizeof(statstg));
    VT_HR_EXIT( stream.Stat(&statstg, STATFLAG_NONAME) );

    if (statstg.cbSize.HighPart != 0)
        VT_HR_EXIT( E_WRITEFAILED );

    VT_HR_EXIT( m_resources.push_back(res) );

    VT_HR_EXIT( m_resources.back().cData.Alloc(statstg.cbSize.LowPart) );

    VT_HR_EXIT( stream.Read(m_resources.back().cData.Ptr(),
                            statstg.cbSize.LowPart,
                            &statstg.cbSize.LowPart) );

    VT_HR_END()
}

HRESULT CPSDWriter::WriteHeader()
{
    HRESULT hr = NOERROR;

    PSDFileHeader hdr;

    if(!m_bInit)
        VT_HR_EXIT(E_NOINIT);
    
    memset(&hdr, 0, sizeof(hdr));

    // calculate length of resources
    int i;
    hdr.uiResourcesLength = 0;
    for(i=0; i<(int)m_resources.size(); i++)
    {
        UInt32 uiSize = 4 + 2 + m_resources[i].sName.SizeBytesPad2() + 4 + m_resources[i].cData.Size();
        if(uiSize & 1)
            uiSize++; // for pad byte
        hdr.uiResourcesLength += uiSize;
    }

    hdr.uiSignature = PSD_FILE_SIGNATURE;
    hdr.uiVersion = Swap16(m_bPSB ? 2 : 1); // version 1 (**PSB** version is 2.)
    hdr.uiChannels = Swap16(m_layers[0].uiChannels);
    hdr.uiHeight = Swap32(m_layers[0].iHeight);
    hdr.uiWidth = Swap32(m_layers[0].iWidth);
    hdr.uiDepth = Swap16(m_b16Bit ? 16 : 8);
    hdr.uiColorMode = Swap16(3); // color format is RGB
    hdr.uiColorTableLength = 0;
    hdr.uiResourcesLength = Swap32(hdr.uiResourcesLength);

    m_hOutput = CreateFileW(m_wstrName,
                            GENERIC_WRITE, // open for read/write 
                            0,                            // do not share 
                            NULL,                         // default security 
                            CREATE_ALWAYS,                // overwrite existing file
                            FILE_ATTRIBUTE_NORMAL,        // normal attributes
                            NULL);                        // no template
    if(m_hOutput==INVALID_HANDLE_VALUE)
        VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));

    VT_HR_EXIT( WriteBytes(m_hOutput, &hdr, sizeof(hdr)) );

    // resources go here
    for(i=0; i<(int)m_resources.size(); i++)
    {
        UInt32 uiSig = PSD_RESOURCE_SIGNATURE;
        VT_HR_EXIT( WriteBytes(m_hOutput, &uiSig, sizeof(uiSig)) );
        UInt16 uiID = Swap16(m_resources[i].uiID);
        VT_HR_EXIT( WriteBytes(m_hOutput, &uiID, sizeof(uiID)) );
        VT_HR_EXIT( m_resources[i].sName.WritePad2(m_hOutput) );
        UInt32 uiSize = m_resources[i].cData.Size();
        UInt32 uiSizeRev = Swap32(uiSize);
        VT_HR_EXIT( WriteBytes(m_hOutput, &uiSizeRev, sizeof(uiSizeRev)) );
        VT_HR_EXIT( WriteBytes(m_hOutput, m_resources[i].cData.Ptr(), uiSize) );
        if(uiSize & 1)
        {
            Byte bPad = 0;
            VT_HR_EXIT( WriteBytes(m_hOutput, &bPad, sizeof(bPad)) );
        }
    }

Exit:
    return hr;
}

HRESULT CPSDWriter::WriteLayers()
{
    HRESULT hr = NOERROR;

    PSBFileLayerMask lmh;
    PSDFileLayerPartA lma;
    PSDFileLayerPartB lmb;
    CPascalString sLayerName;
    int iLayer, iChan;

    for(iLayer=1; iLayer<(int)m_layers.size(); iLayer++)
        if(m_layers[iLayer].iValidRows < m_layers[iLayer].iHeight)
            VT_HR_EXIT(E_FAIL);

    if(m_layers.size() > 1)
    {
        vector<PSBFileChannelInfo> vlci;
        vector<PSDFileChannelInfo> vlci1;

        UInt64 uiDataLength = 0;

        // some layers present
        bool bAddPadByte = false;
        lmh.uiLayerLength = 2; // size of layer count field

        for(iLayer=1; iLayer<(int)m_layers.size(); iLayer++)
        {
            lmh.uiLayerLength += sizeof(PSDFileLayerPartA);
            lmh.uiLayerLength += (m_bPSB ? sizeof(PSBFileChannelInfo) : sizeof(PSDFileChannelInfo)) *
                m_layers[iLayer].uiChannels;
            lmh.uiLayerLength += sizeof(PSDFileLayerPartB);
            VT_HR_EXIT( sLayerName.Create((const WCHAR *)m_layers[iLayer].wsName) );
            lmh.uiLayerLength += sLayerName.SizeBytesPad4();

            // channel image data section length calculation
            for(iChan=0; iChan<m_layers[iLayer].uiChannels; iChan++)
            { 
                PSBFileChannelInfo lci;
                switch(m_layers[iLayer].uiChannels)
                {
                case 3:
                    lci.uiID = (UInt16)iChan;
                    break;
                case 4:
                    lci.uiID = (UInt16)(iChan-1);
                    break;
                case 5:
                default:
                    lci.uiID = (UInt16)(iChan==4 ? -2 : iChan-1);
                    break;
                }
                lci.uiChannelDataLength = 2; // compression type
                
                int iRowSet;
                // add length of all the image rows for this channel
                for(iRowSet = 0; iRowSet<(int)m_layers[iLayer].vRowSet.size(); iRowSet++)
                   lci.uiChannelDataLength += m_layers[iLayer].vRowSet[iRowSet].vChannelLength[iChan];

                // add size of row lengths table if we are compressing
                if(m_bCompress)
                    lci.uiChannelDataLength += (m_bPSB ? sizeof(UInt32) : sizeof(UInt16)) *
                        m_layers[iLayer].iHeight;

                lmh.uiLayerLength += lci.uiChannelDataLength;   
                uiDataLength += lci.uiChannelDataLength;

                lci.uiID = Swap16(lci.uiID);
                if (m_bPSB)
                {
                    lci.uiChannelDataLength = Swap64(lci.uiChannelDataLength);
                    VT_HR_EXIT( vlci.push_back(lci) );
                }
                else
                {
                    PSDFileChannelInfo lci1;
                    lci1.uiID = lci.uiID;
                    lci1.uiChannelDataLength = Swap32((UInt32) lci.uiChannelDataLength);
                    VT_HR_EXIT( vlci1.push_back(lci1) );
                }
            }
        }

        if(lmh.uiLayerLength & 1)
        {
            bAddPadByte = true;
            lmh.uiLayerLength++;
        }

        // add 8 for length of layers section plus length of global mask
        lmh.uiLayerMaskLength = lmh.uiLayerLength + (m_bPSB ? 12 : 8);
        lmh.uiLayerCount = (UInt16)m_layers.size() - 1;

        // swap byte order
        lmh.uiLayerCount = Swap16(lmh.uiLayerCount);
        if (m_bPSB)
        {
            lmh.uiLayerLength = Swap64(lmh.uiLayerLength);
            lmh.uiLayerMaskLength = Swap64(lmh.uiLayerMaskLength);
            VT_HR_EXIT( WriteBytes(m_hOutput, &lmh, sizeof(PSBFileLayerMask)) );
        }
        else
        {
            PSDFileLayerMask lmh1;
            lmh1.uiLayerCount = lmh.uiLayerCount;
            lmh1.uiLayerLength = Swap32((UInt32) lmh.uiLayerLength);
            lmh1.uiLayerMaskLength = Swap32((UInt32) lmh.uiLayerMaskLength);
            VT_HR_EXIT( WriteBytes(m_hOutput, &lmh1, sizeof(PSDFileLayerMask)) );
        }

        // layer and mask information
        int iLCI = 0; // for indexing vLCI
        for(iLayer=1; iLayer<(int)m_layers.size(); iLayer++)
        {
            // write layer headers
            lma.uiLeft = Swap32(m_layers[iLayer].iX);
            lma.uiTop = Swap32(m_layers[iLayer].iY);
            lma.uiRight = Swap32(m_layers[iLayer].iX + m_layers[iLayer].iWidth);
            lma.uiBottom = Swap32(m_layers[iLayer].iY + m_layers[iLayer].iHeight);
            lma.uiChannels = Swap16(m_layers[iLayer].uiChannels);

            VT_HR_EXIT( sLayerName.Create(m_layers[iLayer].wsName) );
           
            lmb.uiSignature = PSD_RESOURCE_SIGNATURE;
            memcpy(lmb.rgbBlendMode, m_layers[iLayer].rgbBlendMode, 4);
            lmb.bOpacity = m_layers[iLayer].bOpacity;
            lmb.bClipping = 0;
            lmb.bFlags = m_layers[iLayer].bVisible ? 0 : 2;
            lmb.bFiller = 0;
            lmb.uiExtraDataLength = Swap32(8 + sLayerName.SizeBytesPad4());
            lmb.uiLayerMaskLength = 0;
            lmb.uiBlendingRangesLength = 0;

            VT_HR_EXIT( WriteBytes(m_hOutput, &lma, sizeof(PSDFileLayerPartA)) );

            if (m_bPSB)
                VT_HR_EXIT( WriteBytes(m_hOutput, &(vlci[iLCI]), sizeof(PSBFileChannelInfo) * m_layers[iLayer].uiChannels) );
            else
                VT_HR_EXIT( WriteBytes(m_hOutput, &(vlci1[iLCI]), sizeof(PSDFileChannelInfo) * m_layers[iLayer].uiChannels) );
            iLCI += m_layers[iLayer].uiChannels; // move to next LCI

            VT_HR_EXIT( WriteBytes(m_hOutput, &lmb, sizeof(PSDFileLayerPartB)) );

            VT_HR_EXIT( sLayerName.WritePad4(m_hOutput) );
        }

        // Skip layer row and image data which was written already.
        LARGE_INTEGER liCurr;
        liCurr.QuadPart = uiDataLength;
        if (!SetFilePointerEx(m_hOutput, liCurr, NULL, FILE_CURRENT))
            VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));

        if(bAddPadByte)
        {
            Byte bZero = 0;
            VT_HR_EXIT( WriteBytes(m_hOutput, &bZero, sizeof(bZero)) );
        }

        // length of global layer mask info = 0
        UInt32 uiTmp = 0;
        VT_HR_EXIT( WriteBytes(m_hOutput, &uiTmp, sizeof(uiTmp)) );
    }
    else
    {
        // layer and mask length = 0
        UInt64 uiTmp = 0;
        VT_HR_EXIT( WriteBytes(m_hOutput, &uiTmp, m_bPSB ? sizeof(UInt64) : sizeof(UInt32)) );
    }

Exit:
    return hr;
}

HRESULT CPSDWriter::WriteChannelRowData(int iLayer, int iChan)
{
    VT_HR_BEGIN()

    int iRowSet;
    PSDLayer *pL = &m_layers[iLayer];

    if(m_bCompress)
    {
        // write table of compressed row lengths
        vt::vector<UInt32> vecRowSizeTable;
		VT_HR_EXIT( vecRowSizeTable.resize(pL->iHeight) );

        int iRow, iY = 0;
        for(iRowSet = 0; iRowSet < (int)pL->vRowSet.size(); iRowSet++)
        {
            int iRowOffset = iChan * pL->vRowSet[iRowSet].iRows; // offset to start of correct channel
            for(iRow = 0; iRow < pL->vRowSet[iRowSet].iRows; iRow++, iY++)
            {
                if (m_bPSB)
                    vecRowSizeTable[iY] = 
                        Swap32((UInt32)(pL->vRowSet[iRowSet].vRows[iRowOffset + iRow].Size()));
                else
                    ((UInt16 *)vecRowSizeTable.begin())[iY] = 
                        Swap16((UInt16)(pL->vRowSet[iRowSet].vRows[iRowOffset + iRow].Size()));
            }
        }

		hr = WriteBytes(m_hOutput, vecRowSizeTable.begin(), 
                        (m_bPSB ? sizeof(UInt32) : sizeof(UInt16)) * pL->iHeight);
    }

    VT_HR_END()
}

HRESULT CPSDWriter::WriteBytes(HANDLE h, const void *pData, DWORD dwSize)
{
    HRESULT hr = NOERROR;
    DWORD dwWritten;
    if(!WriteFile(h, pData, dwSize, &dwWritten, NULL))
        VT_HR_EXIT(HRESULT_FROM_WIN32(GetLastError()));
    if(dwWritten!=dwSize)
        VT_HR_EXIT(E_WRITEFAILED);

Exit:
    return hr;
}

UInt16 CPSDWriter::Swap16(UInt16 ui)
{
    return ((ui & 0xff)<<8) | ((ui & 0xff00)>>8);
}

UInt32 CPSDWriter::Swap32(UInt32 ui)
{
    return ((ui & 0xff)<<24) | ((ui & 0xff00)<<8) | ((ui & 0xff0000)>>8) | ((ui & 0xff000000)>>24);
}

UInt64 CPSDWriter::Swap64(UInt64 ui)
{
    return ((ui & 0xff)<<56) | ((ui & 0xff00)<<40) | ((ui & 0xff0000)<<24) | ((ui & 0xff000000)<<8) |
           ((ui>>56) & 0xff) | ((ui>>40) & 0xff00) | ((ui>>24) & 0xff0000) | ((ui>>8) & 0xff000000);
}

enum
{
    BASE,
    LITERAL,
    RUN,
    LITERAL_RUN
} state;

HRESULT CPSDWriter::CompressBand(const Byte *pbStart, int iElSize, int iPixSize, int iWidth, UInt32 *puiSize)
{
    HRESULT hr = NOERROR;

    int iUnpackedLength = iWidth;
    char *pPackedData = (char *)m_pbBuffer;
    const Byte *pUnpacked = pbStart + iElSize - 1; // reverse order of short
    char *pPacked = pPackedData;
    *puiSize = 0;
    char *lastliteral;
    char *pEnd = pPackedData + COMPRESS_BUFFER_SIZE;

    state = BASE;
    lastliteral = 0;
    while( iUnpackedLength > 0 )
    {
        // Find the longest string of identical bytes.
        int b = *pUnpacked;
        int n = 1;
        if(iElSize==1)
        {
            pUnpacked += iPixSize;
            iUnpackedLength--;
            for( ; ( iUnpackedLength > 0) && ( b == *pUnpacked );
                iUnpackedLength--, pUnpacked+=iPixSize)
            {
                n++;
            }
        }
        else if(iElSize==2)
        {
            if (((UINT_PTR) pUnpacked & 0x1) == 1)
                pUnpacked--;
            else
            {
                pUnpacked += iPixSize + 1;
                iUnpackedLength--;
            }
            for( ; ( iUnpackedLength > 0) && ( b == *pUnpacked ); )
            {
                n++;

                if (((UINT_PTR) pUnpacked & 0x1) == 1)
                    pUnpacked--;
                else
                {
                    pUnpacked += iPixSize + 1;
                    iUnpackedLength--;
                }
            }
        }
    again:
        if( pPacked + 2 >= pEnd )
        {
            VT_HR_EXIT(E_FAIL);

            /* Note - normally if this happens some more elaborate
               handling can occur with writing to the front of
               the buffer - see TIFF implementation.
               Since this is being used to compress highly compressable
               data, we will not worry.
            */
        }

        switch( state )
        {
        case BASE: // initial state, set run/literal
            if( n > 1 )
            {
                state = RUN;
                if( n > 128 )
                {
                    *pPacked++ = -127;
                    *pPacked++ = (char)b;
                    n -= 128;
                    goto again;
                }
                *pPacked++ = char(-(n-1));
                *pPacked++ = (char)b;
            }
            else
            {
                lastliteral = pPacked;
                *pPacked++ = 0;
                *pPacked++ = (char)b;
                state = LITERAL;
            }
            break;
        case LITERAL: // last object was literal string
            if( n > 1 )
            {
                state = LITERAL_RUN;
                if( n > 128 )
                {
                    *pPacked++ = -127;
                    *pPacked++ = (char)b;
                    n -= 128;
                    goto again;
                }
                *pPacked++ = char(-(n-1)); // encode run
                *pPacked++ = (char)b;
            }
            else
            {           // extend literal
                if (++(*lastliteral) == 127)
                {
                    state = BASE;
                }
                *pPacked++ = (char)b;
            }
            break;
        case RUN: // last object was run
            if( n > 1 )
            {
                if( n > 128 )
                {
                    *pPacked++ = -127;
                    *pPacked++ = (char)b;
                    n -= 128;
                    goto again;
                }
                *pPacked++ = char(-(n-1));
                *pPacked++ = (char)b;
            }
            else
            {
                lastliteral = pPacked;
                *pPacked++ = 0;
                *pPacked++ = (char)b;
                state = LITERAL;
            }
            break;
        case LITERAL_RUN: // literal followed by a run
            /*
             * Check to see if previous run should
             * be converted to a literal, in which
             * case we convert literal-run-literal
             * to a single literal.
             */
            if( n == 1 && pPacked[-2] == (char)-1 &&
                *lastliteral < 126 )
            {
                state = (((*lastliteral) += 2) == 127 ?
                    BASE : LITERAL);
                pPacked[-2] = pPacked[-1];  // replicate
            }
            else
            {
                state = RUN;
            }
            goto again;
        }
    }
    *puiSize += (UInt32)(pPacked - pPackedData);

Exit:
    return hr;
}

void CPSDWriter::CopyBand(const Byte *pbSrc, Byte *pbDst, int iElSize, int iPixSize, int iWidth)
{
    int i;
    if(iElSize==1)
    {
        for(i=0; i<iWidth; i++, pbSrc += iPixSize)
            *pbDst++ = pbSrc[0];
    } 
    else if(iElSize==2)
    {
        for(i=0; i<iWidth; i++, pbSrc += iPixSize)
        {
            *pbDst++ = pbSrc[1]; // reverse order of short
            *pbDst++ = pbSrc[0];
        }
    }
}


HRESULT CPascalString::Create(const WCHAR * pwsz)
{
    HRESULT hr = NOERROR;

    memset(m_rgch, 0, 256);

    vt::string strMB;
    VtWideCharToMultiByte(strMB, pwsz);

    size_t uiLen = strlen(strMB);
    if(uiLen>255)
        uiLen = 255;
    m_rgch[0] = (Byte)uiLen;

    memcpy(m_rgch + 1, (const char*)strMB, Length());

    return hr;
}

HRESULT CPascalString::Write(HANDLE h)
{
    HRESULT hr = NOERROR;
    DWORD dwWritten;
    if(!WriteFile(h, m_rgch, SizeBytes(), &dwWritten, NULL))
        VT_HR_EXIT( HRESULT_FROM_WIN32(GetLastError()) );
    if(dwWritten!=SizeBytes())
        VT_HR_EXIT(E_WRITEFAILED);
Exit:
    return hr;
}

HRESULT CPascalString::WritePad2(HANDLE h)
{
    HRESULT hr = NOERROR;
    DWORD dwWritten;
    if(!WriteFile(h, m_rgch, SizeBytesPad2(), &dwWritten, NULL))
        VT_HR_EXIT( HRESULT_FROM_WIN32(GetLastError()) );
    if(dwWritten!=SizeBytesPad2())
        VT_HR_EXIT(E_WRITEFAILED);
Exit:
    return hr;
}

HRESULT CPascalString::WritePad4(HANDLE h)
{
    HRESULT hr = NOERROR;
    DWORD dwWritten;
    if(!WriteFile(h, m_rgch, SizeBytesPad4(), &dwWritten, NULL))
        VT_HR_EXIT( HRESULT_FROM_WIN32(GetLastError()) );
    if(dwWritten!=SizeBytesPad4())
        VT_HR_EXIT(E_WRITEFAILED);
Exit:
    return hr;
}

#endif
