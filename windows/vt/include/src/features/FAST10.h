//
// Copyright (c) Microsoft Corporation. All rights reserved.
//
#pragma once

#define Assert(cond)

class FAST10Corners
{
public:
    int count;
    int *offset;
    int *score;

protected:
    HRESULT hr; // return from vt::vector allocations
    static const size_t reserveSize = 4096;
    size_t allocSize;
    vt::vector<int> offsetVec;
    vt::vector<int> scoreVec;
    vt::vector<int> offsetMaxVec;
    int countMax;
    int *offsetMax;

public:
    HRESULT GetHr(void) { return hr; }
    
    int GetCount(void) { return count; }
    void GetCoordinate(int i, int &x, int &y, int &scoreRet, int stride)
    {
        if (i>=count) { x=-1; i=-1; return; }
        int offs = offset[i];
        x = offs % stride;
        y = offs / stride;
        scoreRet = score[i];
    }
    
    int GetMaxCount(void) { return countMax; }
    void GetMaxCoordinate(int i, int &x, int &y, int &scoreRet, int stride)
    {
        if (i>=countMax) { x=-1; i=-1; return; }
        int offs = offsetMax[i];
        x = offs % stride;
        y = offs / stride;
        scoreRet = score[i];
    }

    inline void AddMax(int i)
    {
        // code structure dictates that there can't be more Max corners than corners, so no bounds check needed here
        offsetMax[countMax] = i;
        countMax++;
    }
    inline void CheckAdd8(void) { if ((count+8)>=(int)allocSize) { IncreaseSize(); } }
    inline void CheckAdd16(void) { if ((count+16)>=(int)allocSize) { IncreaseSize(); } }
    inline void Add(int i) 
    {  
        if ((count+1)>=(int)allocSize) { IncreaseSize(); }
        offset[count] = i;
        count++;
    }
    inline void AddFast(int i) 
    {  
        offset[count] = i;
        count++;
    }
    void Reset(void) 
    { 
        count = 0;
        countMax = 0;
    }
    void ResetMax(void) { countMax = 0; } // for testing only
    void IncreaseSize(void)
    {
        allocSize += reserveSize;

        HRESULT hr0 = offsetVec.resize(allocSize);
        HRESULT hr1 = scoreVec.resize(allocSize);
        HRESULT hr2 = offsetMaxVec.resize(allocSize);
        if (hr0!=S_OK) { hr = hr0; }
        if (hr1!=S_OK) { hr = hr1; }
        if (hr2!=S_OK) { hr = hr2; }

        offset = &offsetVec[0];
        score = &scoreVec[0];
        offsetMax = &offsetMaxVec[0];
    }
    FAST10Corners(size_t initSize = reserveSize) 
    { 
        allocSize = initSize;

        hr = S_OK;
        HRESULT hr0 = offsetVec.resize(allocSize);
        HRESULT hr1 = scoreVec.resize(allocSize);
        HRESULT hr2 = offsetMaxVec.resize(allocSize);
        if (hr0!=S_OK) { hr = hr0; }
        if (hr1!=S_OK) { hr = hr1; }
        if (hr2!=S_OK) { hr = hr2; }

        offset = &offsetVec[0];
        score = &scoreVec[0];
        offsetMax = &offsetMaxVec[0];
        Reset();
    }
    ~FAST10Corners()
    {
    }
};

extern void FASTCornerDetect10C(FAST10Corners &corners, const unsigned char *src, int w, int h, int stride, int threshold);
#ifdef __ARM_NEON__
extern void FASTCornerDetect10Neon(FAST10Corners &corners, const unsigned char *src, int w, int h, int stride, int threshold);
#endif
#if (defined(_M_IX86) || defined(_M_AMD64))
extern void FASTCornerDetect10SSE(FAST10Corners &corners, const unsigned char *src, int w, int h, int stride, int threshold);
#endif
extern void FASTNonMaxSuppression(FAST10Corners &corners, const unsigned char *src, int stride, int threshold);

HRESULT FAST10Detect(vt::vector<HARRIS_FEATURE_POINT>& pts, 
    const vt::CByteImg& src, const vt::CRect& srcRect,
    int border = 5, float threshold = 20.f);


// end