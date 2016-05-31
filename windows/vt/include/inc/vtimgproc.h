//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Main include file for imgproc routines
// 
//  History:
//      2009/05/25 - mattu
//			Created
//
//------------------------------------------------------------------------

// vtimgproc.h

#pragma once

#include "vtcore.h"

#include "../src/imgproc/vt_edgedetect_common.h"
#include "../src/imgproc/vt_edgedetect.h"
#ifndef VT_GCC
#include "..\src\imgproc\vt_denoise.h"
#include "..\src\imgproc\vt_linedetector.h"
#include "..\src\imgproc\vt_fitEllipse.h"
#include "..\src\imgproc\vt_vanishingpoint_extractor.h"
#include "..\src\imgproc\vt_edgedetect_util.h"
#include "..\src\imgproc\vt_lucasKanade.h"
#include "..\src\imgproc\vt_CAremove.h"
#include "..\src\imgproc\vt_VignetteRemove.h"
#include "..\src\imgproc\vt_completion.h"
#endif
