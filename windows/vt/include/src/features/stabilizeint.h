//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      stabilization internal header
//
//------------------------------------------------------------------------------
#pragma once
#include "brief.h"
#include "tracker.h"

//------------------------------------------------------------------------------
// internal class for stabilizer frame data
//------------------------------------------------------------------------------
class CStabilizerFrameData::CStabilizerFrameDataInternal
{
public:
    CStabilizerFrameDataInternal()
        :m_blurLevel(0.f)
        ,m_imgWidth(0),m_imgHeight(0)
    {};
    ~CStabilizerFrameDataInternal() {};

    int m_imgWidth;
    int m_imgHeight;
    float m_blurLevel;
    vt::vector<HARRIS_FEATURE_POINT> curFP;
    vt::vector<BriefDesc<BRIEF_DESC_SIZE>> curDesc;
};

