//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      AVI file reader/writer class
//
//  History:
//      2004/11/8-swinder
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#if !defined(VT_WINRT)

#include "vt_global.h"
#include "aviapi.h"

using namespace vt;

////////////////////////////////////////////////////////////////////////
// CAviFileSrc
////////////////////////////////////////////////////////////////////////

CAviFileSrc::CAviFileSrc()
    : m_pAviFile(NULL), m_pAviStream(NULL), m_pGetFrame(NULL), m_iLoc(0)
{
	AVIFileInit();
}

CAviFileSrc::~CAviFileSrc()
{
    CloseFile();
    AVIFileExit();
}

HRESULT CAviFileSrc::Seek(Int64 iPos)
{
    if(m_pAviStream==NULL)
        return SetError(E_NOINIT);
    if(iPos<0)
        return SetError(E_INVALIDARG);
    if((DWORD)iPos >= m_asi.dwLength)
        return S_EOF;
    m_iLoc = iPos;
    return NOERROR;
}

HRESULT CAviFileSrc::OpenFile(const WCHAR * pchFile)
{
    if(pchFile==NULL)
        return SetError(E_POINTER);

    ClearError();

    HRESULT hr = AVIFileOpenW(&m_pAviFile, pchFile, OF_READ, NULL);
    if(FAILED(hr))
        return SetError(hr);

    hr = AVIFileInfoW(m_pAviFile, &m_afi, sizeof(m_afi));
    if(FAILED(hr))
    {
        CloseFile();
        return SetError(hr);
    }

    int i;
    bool bFound = false;
    for(i=0; !bFound && (DWORD)i < m_afi.dwStreams; i++)
    {
        hr = AVIFileGetStream(m_pAviFile, &m_pAviStream, 0L, i);
        if(FAILED(hr))
        {
            CloseFile();
            return SetError(hr);
        }
        hr = AVIStreamInfoW(m_pAviStream, &m_asi, sizeof(m_asi));
        if(FAILED(hr))
        {
            CloseFile();
            return SetError(hr);
        }
        if(m_asi.fccType == streamtypeVIDEO)
            bFound = true;
        else
            AVIStreamRelease(m_pAviStream);
    }
    if(!bFound)
    {
        CloseFile();
        return SetError(AVIERR_NODATA);
    }

    m_pGetFrame = AVIStreamGetFrameOpen(m_pAviStream, NULL);
    if(m_pGetFrame==NULL)
        return SetError(AVIERR_NOCOMPRESSOR);

    m_dib.Destroy(); // ensure no mem allocd

    return NOERROR;
}

void CAviFileSrc::CloseFile()
{
    if(m_pGetFrame)
    {
        AVIStreamGetFrameClose(m_pGetFrame);
        m_pGetFrame = NULL;
    }
    if(m_pAviStream)
    {
        AVIStreamRelease(m_pAviStream);
        m_pAviStream = NULL;
    }
    if(m_pAviFile)
    {
        AVIFileRelease(m_pAviFile);
        m_pAviFile = NULL;
    }
}

void CAviFileSrc::GetSize(int *piWidthRtn, int *piHeightRtn)
{
    if(piWidthRtn!=NULL)
        *piWidthRtn = m_asi.rcFrame.right - m_asi.rcFrame.left;
    if(piHeightRtn!=NULL)
        *piHeightRtn = m_asi.rcFrame.bottom - m_asi.rcFrame.top;
}

HRESULT CAviFileSrc::GetFrame(Int64 iFrame, CImg &img)
{
    if(iFrame<0 || (DWORD)iFrame >= m_asi.dwLength)
        return SetError(E_INVALIDARG);
    if(m_pAviStream==NULL || m_pGetFrame==NULL)
        return SetError(E_NOINIT);

    // dont care if img doesnt match size of frame
//  int iW, iH, iSW, iSH;
//  iW = img.Width();
//  iH = img.Height();
//  GetSize(&iSW, &iSH);
//  if(iW!=iSW || iH!=iSH)
//      return E_INVALIDARG;

    BITMAPINFOHEADER *pbmi = (BITMAPINFOHEADER *)AVIStreamGetFrame(m_pGetFrame, (LONG)iFrame);
    if(pbmi==NULL)      
        return SetError(AVIERR_ERROR);

    HRESULT hr;
    
    if(m_dib.IsValid())
        hr = m_dib.Update(pbmi);
    else
        hr = m_dib.Create(pbmi); // first time create

    if(SUCCEEDED(hr))
        hr = m_dib.GetPixels(img);

    if(FAILED(hr))
        return SetError(hr);

    return hr;
}


////////////////////////////////////////////////////////////////////////
// CAviFileDst
////////////////////////////////////////////////////////////////////////

CAviFileDst::CAviFileDst()
    : m_pAviFile(NULL), m_pAviStream(NULL), m_iLoc(0), m_iEnumComp(0), 
      m_bCompressStart(false), m_fFrameRate(30.0), m_iCMap(DIBCMAP_GRAYSCALE),
      m_bUseAlphaPlane(false)
{
    AVIFileInit();
    memset(&m_cvars, 0, sizeof(COMPVARS));
    m_cvars.cbSize = sizeof(COMPVARS);
}

CAviFileDst::~CAviFileDst()
{
    CloseFile();
    if(m_cvars.dwFlags & ICMF_COMPVARS_VALID)
    {
        if(m_cvars.hic)
            ICClose(m_cvars.hic);
        m_cvars.hic = 0;
        ICCompressorFree(&m_cvars);
    }
    AVIFileExit();
}

HRESULT CAviFileDst::StartEnumCompressors(CImg &img)
{
    m_iEnumComp = 0;

    // initialise bmih
    memset(&m_bmih, 0, sizeof(BITMAPINFOHEADER));
    m_bmih.biSize = sizeof(BITMAPINFOHEADER);
    m_bmih.biPlanes = 1;
    m_bmih.biCompression = BI_RGB;

    if (EL_FORMAT(img.GetType()) == EL_FORMAT_BYTE && img.Bands() == 1)
        m_bmih.biBitCount = 8;
    else if (EL_FORMAT(img.GetType()) == EL_FORMAT_BYTE && (img.Bands() == 4 || img.Bands() == 3))
        m_bmih.biBitCount = 24;
    else
        return SetError(E_BADFORMAT);

    if(img.IsValid())
    {
        m_bmih.biWidth = img.Width();
        m_bmih.biHeight = img.Height();
    }

    return NOERROR;
}

bool CAviFileDst::NextEnumCompressors(ICINFO *pInfoRtn)
{
    if(pInfoRtn==NULL)
        return false;

    while(ICInfo(ICTYPE_VIDEO, m_iEnumComp++, pInfoRtn))
    {
        HIC hic = ICOpen(pInfoRtn->fccType, pInfoRtn->fccHandler, ICMODE_QUERY); 
        if (hic) 
        {
            // check that it handles the input format
            if(ICCompressQuery(hic, &m_bmih, NULL) != ICERR_OK)
            {
                ICClose(hic);
                continue;
            }
     
            // Find out the compressor name and extra info. 
            LRESULT lr = ICGetInfo(hic, pInfoRtn, sizeof(ICINFO));
            ICClose(hic);

            if(lr==0)
                continue;

            return true;
        }
    }

    return false;
}

HRESULT CAviFileDst::SetCompressor(FOURCC fcc, int iKey, int iQuality, int iRate, bool bUseDialog, HWND hwndApp)
{
    // cannot set compressor if file is open
    if(m_pAviFile)
        return SetError(E_FAIL);

    // free compressor if already set
    if(m_cvars.dwFlags & ICMF_COMPVARS_VALID)
    {
        if(m_cvars.hic)
            ICClose(m_cvars.hic);
        m_cvars.hic = 0;
        ICCompressorFree(&m_cvars);
    }
    memset(&m_cvars, 0, sizeof(COMPVARS));
    m_cvars.cbSize = sizeof(COMPVARS);

    // if fcc == 0 then avi is uncompressed
    m_cvars.fccType = ICTYPE_VIDEO;
    m_cvars.fccHandler = fcc;

    if(fcc!=0)
    {
        m_cvars.hic = ICOpen(ICTYPE_VIDEO, fcc, ICMODE_COMPRESS);
        if(m_cvars.hic==0)
            return SetError(E_NOCOMPRESSOR);

        if(iKey>=0)
            m_cvars.lKey = iKey;
        else
            m_cvars.lKey = ICGetDefaultKeyFrameRate(m_cvars.hic); 
        
        if(iQuality>=0)
            m_cvars.lQ = iQuality;
        else
            m_cvars.lQ = ICGetDefaultQuality(m_cvars.hic);
        
        m_cvars.lDataRate = iRate;

        if(bUseDialog && ICQueryConfigure(m_cvars.hic))
        {
            // not quite sure  if the stuff selected here gets overwritten by COMPVARS settings later
            LRESULT iStatus = ICConfigure(m_cvars.hic, hwndApp);
            if(iStatus!=1)
            {
                // error or cancel button hit
                ICClose(m_cvars.hic);
                m_cvars.hic = 0;
                return SetError(E_NOINIT);
            }
        }

        m_cvars.dwFlags = ICMF_COMPVARS_VALID;
    }

    return NOERROR;
}

HRESULT CAviFileDst::ChooseCompressorUsingDialog(CImg &img, HWND hwndApp)
{
    // cannot set compressor if file is open
    if(m_pAviFile)
        return SetError(E_FAIL);

    HRESULT hr = StartEnumCompressors(img);
    if(FAILED(hr))
        return hr;

    // free compressor if already set
    if(m_cvars.dwFlags & ICMF_COMPVARS_VALID)
    {
        if(m_cvars.hic)
            ICClose(m_cvars.hic);
        m_cvars.hic = 0;
        ICCompressorFree(&m_cvars);
    }
    memset(&m_cvars, 0, sizeof(COMPVARS));
    m_cvars.cbSize = sizeof(COMPVARS);

    m_cvars.fccType = ICTYPE_VIDEO;
    m_cvars.fccHandler = 0;

    if(ICCompressorChoose(hwndApp, ICMF_CHOOSE_ALLCOMPRESSORS | ICMF_CHOOSE_DATARATE | ICMF_CHOOSE_KEYFRAME,
        &m_bmih, NULL, &m_cvars, "Compression Settings")==FALSE)
        return SetError(E_NOINIT);

    return NOERROR;
}

HRESULT CAviFileDst::OpenFile(const WCHAR * pchFile)
{
    if(pchFile==NULL)
        return SetError(E_POINTER);

    // avi write needs com started
    VtStartCOM();

    CloseFile();

    HRESULT hr;
    hr = AVIFileOpenW(&m_pAviFile, pchFile, 
                      OF_CREATE | OF_WRITE | OF_SHARE_DENY_WRITE, NULL);
    if(FAILED(hr))
        return SetError(hr);

    m_iLoc = 0;

    return NOERROR;
}

void CAviFileDst::CloseFile()
{
    if(m_bCompressStart && m_cvars.hic!=0)
    {
        // finish compression
        ICSeqCompressFrameEnd(&m_cvars);
    }

    m_bCompressStart = false;

    if(m_pAviStream)
    {
        AVIStreamRelease(m_pAviStream);
        m_pAviStream = NULL;
    }
    if(m_pAviFile)
    {
        AVIFileRelease(m_pAviFile);
        m_pAviFile = NULL;
    }
}

const CImg& CAviFileDst::operator+=(const CImg &cImg)
{
    HRESULT hr;

    const CImg *pImg = &cImg;

    if (EL_FORMAT(pImg->GetType()) != EL_FORMAT_BYTE || (pImg->Bands() != 1 && pImg->Bands()!=4))
	{
        SetError(E_NOTIMPL);
		goto Exit;
	}

    if(m_iLoc==0)
    {
        // first frame - construct the DIB
        if(pImg->Bands()==1 && EL_FORMAT(pImg->GetType())==EL_FORMAT_BYTE)
        {
            hr = m_cdib.Create8Bit(pImg->Width(), pImg->Height());
            if(SUCCEEDED(hr))
                if(m_iCMap!=DIBCMAP_GRAYSCALE)
                    m_cdib.SetCMap(m_iCMap);
        }
        else if(m_bUseAlphaPlane)
            hr = m_cdib.Create32Bit(pImg->Width(), pImg->Height());
        else
            hr = m_cdib.Create24Bit(pImg->Width(), pImg->Height());
        if(FAILED(hr))
		{
			SetError(hr);
			goto Exit;
		}
    }

    // copy image to DIB
    hr = m_cdib.SetPixels(*pImg, 0, 0);
    if(FAILED(hr))
	{
        SetError(hr);
		goto Exit;
	}
    
    PDIB pbmi = m_cdib.GetDibInfo();

    if(m_iLoc==0)
    {
        // First time to put object, so open new stream
        AVISTREAMINFOW asi;

        // fill in the header for the video stream
        memset(&asi, 0, sizeof(AVISTREAMINFO));
        asi.fccType = streamtypeVIDEO;
        asi.fccHandler = m_cvars.hic ? m_cvars.fccHandler : 0;
        asi.dwScale = SCALE_FACTOR;
        asi.dwRate = (DWORD)(m_fFrameRate * SCALE_FACTOR);
        asi.dwSuggestedBufferSize = pbmi->biSizeImage;
        SetRect(&asi.rcFrame, 0, 0, pbmi->biWidth, pbmi->biHeight);

        hr = AVIFileCreateStreamW(m_pAviFile, &m_pAviStream, &asi);
        if(FAILED(hr))
		{
			SetError(hr);
			return cImg;
		}		

        if(m_cvars.hic!=0)
        {
            // initialise compression
            DWORD dwFormatSize = ICCompressGetFormatSize(m_cvars.hic, pbmi);
            LPVOID lpFormat = VT_NOTHROWNEW Byte [dwFormatSize];
            if (NULL == lpFormat)
			{
				hr = E_OUTOFMEMORY;
				goto StreamError;
			}

            if(ICCompressGetFormat(m_cvars.hic, pbmi, lpFormat)!=ICERR_OK)
            {
                delete [] lpFormat;
                hr = E_NOCOMPRESSOR;
                goto StreamError;
            }
            
            hr = AVIStreamSetFormat(m_pAviStream, 0, lpFormat, dwFormatSize);

            delete [] lpFormat;

            if(FAILED(hr))
                goto StreamError;
            
            if(!ICSeqCompressFrameStart(&m_cvars, (BITMAPINFO *)pbmi))
            {
                hr = E_NOCOMPRESSOR;
                goto StreamError;
            }

            m_bCompressStart = true;
        }
        else
        {
            // no compression
            hr = AVIStreamSetFormat(m_pAviStream, 0, pbmi, 
                                    (LONG)m_cdib.GetDibHeaderSize());
            if(FAILED(hr))
                goto StreamError;
        }
    }

    if(m_cvars.hic!=0)
    {
        // write one compressed frame
        BOOL bKey = FALSE;
        LONG lSize = pbmi->biSizeImage;
        LPVOID lpBits = ICSeqCompressFrame(&m_cvars, 0, m_cdib.GetDibBits(),
                                           &bKey, &lSize);
        if(lpBits==NULL)
            hr = E_FAIL;
        else
            hr = AVIStreamWrite(m_pAviStream, (LONG)m_iLoc, 1L, lpBits, lSize,
                (bKey ? AVIIF_KEYFRAME : 0), NULL, NULL);
    }
    else
    {
        // write one uncompressed frame
        hr = AVIStreamWrite(m_pAviStream, (LONG)m_iLoc, 1L, m_cdib.GetDibBits(), 
                            pbmi->biSizeImage, AVIIF_KEYFRAME, NULL, NULL);
    }
    m_iLoc++;
    
    if(FAILED(hr))
    {
StreamError:
        if(m_pAviStream)
            AVIStreamRelease(m_pAviStream);
        m_pAviStream = NULL;
		SetError(hr);
    }

Exit:
	return cImg;
}

#endif
