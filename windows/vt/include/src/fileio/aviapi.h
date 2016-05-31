//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      AVI file reader/writer class
//
//  History:
//      2004/11/8-swinder
//            Created
//
//------------------------------------------------------------------------

#pragma once

#if !defined(VT_WINRT)

#include "vtcommon.h"

namespace vt {

// frame rate scale factor
#define SCALE_FACTOR 10000

// reading avi file code example:
//  CAviFileSrc myavisrc;
//  hr = myavisrc.OpenFile(L"myfile.avi");
//  int w, h;
//  myavisrc.GetSize(&w, &h);
//  int framecount = myavisrc.GetLength();
//  CRGBAImg myimg(w,h);
//  hr = myavisrc.GetFrame(framenum, myimg);
//  myavisrc.CloseFile();
//
// writing avi file code example:
//  CRGBAImg myimg(w,h);
//   ...
//  CAviFileDst myavidst;
//  hr = myavidst.OpenFile(L"myfile.avi");
//  myavidst.SetFrameRate(framespersecond);
//  [optionally hr = ChooseCompressorUsingDialog(myimg, NULL or hwndapp)]
//  myavidst += myimg;
//  myavidst.CloseFile();

////////////////////////////////////////////////////////////////////////
// CAviFileSrc - Avi file access class
////////////////////////////////////////////////////////////////////////

class CAviFileSrc : public CErrorBase
{
public:
    CAviFileSrc();
    ~CAviFileSrc();

    HRESULT OpenFile(const WCHAR * pszFile);
    void CloseFile();

    virtual Int64 GetLength() { return (Int64)m_asi.dwLength; }
    virtual HRESULT Seek(Int64 iPos);

    // other members
    double GetFrameRate() { return m_asi.dwRate/(double)m_asi.dwScale; }
    AVIFILEINFOW   *GetFileInfo()   { return &m_afi; }
    AVISTREAMINFOW *GetStreamInfo() { return &m_asi; }
    void GetSize(int *piWidthRtn, int *piHeightRtn);
    HRESULT GetFrame(Int64 iFrame, CImg &img);
    Int64 GetTimeOfFrame(Int64 iFrame) 
        { return 10000*(Int64)AVIStreamSampleToTime(m_pAviStream, (int)iFrame); }

private:

    PAVIFILE       m_pAviFile;
    PAVISTREAM     m_pAviStream;
    PGETFRAME      m_pGetFrame;
    AVIFILEINFOW   m_afi;
    AVISTREAMINFOW m_asi;
    CDib  m_dib;
    Int64 m_iLoc;
};

class CAviFileDst : public CErrorBase
{
public:
    CAviFileDst();
    ~CAviFileDst();

    HRESULT OpenFile(const WCHAR * pszFile);
    void CloseFile();

    // other members
    const CImg &operator+=(const CImg &cImg);
    double GetFrameRate() { return m_fFrameRate; }
    void SetFrameRate(double fRate) { m_fFrameRate = fRate; }
    Int64 GetLength() { return m_iLoc; }

    void SetColorMap(int iType) { m_iCMap = iType; } // one of DIBCMAP_GRAYSCALE, DIBCMAP_PSEUDO1, DIBCMAP_COLDHOT

    // call if you want uncompressed 32-bit ARGB format
    void EnableAlphaPlane(bool bUse) { m_bUseAlphaPlane = bUse; }

    HRESULT ChooseCompressorUsingDialog(CImg &img, HWND hwndApp);
    HRESULT StartEnumCompressors(CImg &img);
    bool NextEnumCompressors(ICINFO *pInfoRtn);
    // key is in milliseconds and can be -1 for default value
    // quality is from 0-10000 and can be -1 for default value
    // rate is in kilobits per second and must be specified
    HRESULT SetCompressor(FOURCC fcc, int iKey, int iQuality, int iRate, 
                          bool bUseDialog = false, HWND hwndApp = NULL);

private:
    bool m_bUseAlphaPlane;
    int  m_iCMap;
    CDib m_cdib;
    bool m_bCompressStart;
    COMPVARS m_cvars;
    int m_iEnumComp;
    PAVIFILE m_pAviFile;
    PAVISTREAM m_pAviStream;
    double m_fFrameRate;
    Int64 m_iLoc;
    BITMAPINFOHEADER m_bmih;
};

};

#endif
