//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for drawing onto images
//
//  History:
//      2005/7 - swinder
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_image.h"

namespace vt {
	
// draws a simple 1 pixel wide line (clipped to image boundaries)
template <class T>
void VtDrawLine(CTypedImg<T> &img, float x1, float y1, float x2, float y2,
                const T *pColor, const RECT *prctclip = NULL);


template <class T>
void VtDrawCircle(CTypedImg<T> &img, float x, float y, float r,
                      const T *pColor, const RECT *prctclip = NULL);

// draw ellipse based on inverse covariance matrix
// scale=1 draws the unit variance ellipse. make it bigger by specifying larger scale
template <class T>
void VtDrawEllipse(CTypedImg<T> &img, float x, float y, CMtx2x2f mInvCov, float scale, const T *pColor, const RECT *prct = NULL);

};
