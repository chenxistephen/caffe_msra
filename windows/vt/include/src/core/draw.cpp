//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      drawing functions
//
//  History:
//      2005/7-swinder
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_image.h"
#include "vt_atltypes.h"
#include "vt_draw.h"

using namespace vt;

#define CLIP_TOP        1
#define CLIP_BOTTOM     2
#define CLIP_LEFT       4
#define CLIP_RIGHT      8

int GetOutCode(float x, float y, float xmin, float ymin, float xmax, float ymax)
{
    int code = 0;
    if(y>ymax)
        code = CLIP_TOP;
    else if(y<ymin)
        code = CLIP_BOTTOM;
    if(x>xmax)
        code |= CLIP_RIGHT;
    else if(x<xmin)
        code |= CLIP_LEFT;
    return code;
}

// draws a simple 1 pixel wide line (clipped to image boundaries)
template <class T>
void vt::VtDrawLine(CTypedImg<T> &img, float x1, float y1, float x2, float y2, const T *pColor, const RECT *prct)
{
    int w = img.Width();
    int h = img.Height();
    int b = img.Bands();

    float ymin, ymax, xmin, xmax;
    if(prct)
    {
        CRect rctclip = img.ClipRect(prct);
        ymin = (float)rctclip.top;
        ymax = (float)(rctclip.bottom - 1);
        xmin = (float)rctclip.top;
        xmax = (float)(rctclip.right - 1);
    }
    else
    {
        ymin = 0;
        ymax = (float)(h - 1);
        xmin = 0;
        xmax = (float)(w - 1);
    }

    if(GetOutCode(x1, y1, xmin, ymin, xmax, ymax) || GetOutCode(x2, y2, xmin, ymin, xmax, ymax))
    {
        // clip
        bool accept = false;
        bool done = false;
        int outcode1 = GetOutCode(x1, y1, xmin, ymin, xmax, ymax);
        int outcode2 = GetOutCode(x2, y2, xmin, ymin, xmax, ymax);
        do
        {
            if(outcode1==0 && outcode2==0)
            {
                accept = true;
                done = true;
            }
            else if(outcode1 & outcode2)
                done = true;
            else
            {
                int outcode_out;
                float x = 0.f, y = 0.f;
                if(outcode1)
                    outcode_out = outcode1;
                else
                    outcode_out = outcode2;
                if(outcode_out & CLIP_TOP)
                {
                    x = x1 + (x2 - x1) * (ymax - y1)/(y2 - y1);
                    y = ymax;
                }
                else if(outcode_out & CLIP_BOTTOM)
                {
                    x = x1 + (x2 - x1) * (ymin - y1)/(y2 - y1);
                    y = ymin;
                }
                else if(outcode_out & CLIP_RIGHT)
                {
                    x = xmax;
                    y = y1 + (y2 - y1) * (xmax - x1)/(x2 - x1);
                }
                else if(outcode_out & CLIP_LEFT)
                {
                    x = xmin;
                    y = y1 + (y2 - y1) * (xmin - x1)/(x2 - x1);
                }
                if(outcode_out==outcode1)
                {
                    x1 = x;
                    y1 = y;
                    outcode1 = GetOutCode(x1, y1, xmin, ymin, xmax, ymax);
                }
                else
                {
                    x2 = x;
                    y2 = y;
                    outcode2 = GetOutCode(x2, y2, xmin, ymin, xmax, ymax);
                }
            }
        } while(!done);
        if(!accept)
            return;
    }

    float dx = x2 - x1;
    float dy = y2 - y1;
    if(fabs(dx) >= fabs(dy))
    {
        if(dx==0)
        {
            T *p = img.Ptr((int)floor(x1+0.5f), (int)floor(y1+0.5f));
            memcpy(p, pColor, b * sizeof(T));
            return;
        }

        if(dx<0)
        {
            float tmp = x1;
            x1 = x2;
            x2 = tmp;
            y1 = y2;
            dx = -dx;
            dy = -dy;
        }

        // step dx
        float x;
        float y = y1;
        float dydx = dy/dx;
        for(x = x1; x<=x2; x+=1.0f, y += dydx)
        {
            T *p = img.Ptr((int)floor(x+0.5f), (int)floor(y+0.5f));
            memcpy(p, pColor, b * sizeof(T));
        }
    }
    else
    {
        if(dy<0)
        {
            float tmp = y1;
            y1 = y2;
            y2 = tmp;
            x1 = x2;
            dx = -dx;
            dy = -dy;
        }

        // step dy
        float x = x1;
        float y;
        float dxdy = dx/dy;
        for(y = y1; y<=y2; y+=1.0f, x += dxdy)
        {
            T *p = img.Ptr((int)floor(x+0.5f), (int)floor(y+0.5f));
            memcpy(p, pColor, b * sizeof(T));
        }
    }
}

template void vt::VtDrawLine(CTypedImg<Byte> &img, float x1, float y1, float x2, float y2,
                             const Byte *pColor, const RECT *prct);
template void vt::VtDrawLine(CTypedImg<float> &img, float x1, float y1, float x2, float y2,
                             const float *pColor, const RECT *prct);

template <class T>
static void VtPlot(CTypedImg<T> &img, float fx, float fy, const T *pColor, const RECT *prct)
{
    int x = (int)floor(fx + 0.5f);
    int y = (int)floor(fy + 0.5f);
    if(x>=prct->left && x<prct->right && y>=prct->top && y<prct->bottom)
    {
        memcpy(img.Ptr(x,y), pColor, img.PixSize());
    }
}

template <class T>
static void CirclePoints(CTypedImg<T> &img, float x, float y, float dx, float dy, const T *pColor, const RECT *prct)
{
    float xpx = x + dx;
    float xmx = x - dx;
    float xpy = x + dy;
    float xmy = x - dy;
    float ypy = y + dy;
    float ymy = y - dy;
    float ypx = y + dx;
    float ymx = y - dx;

    VtPlot(img, xpx, ypy, pColor, prct);
    VtPlot(img, xpy, ypx, pColor, prct);
    VtPlot(img, xpy, ymx, pColor, prct);
    VtPlot(img, xpx, ymy, pColor, prct);
    VtPlot(img, xmx, ymy, pColor, prct);
    VtPlot(img, xmy, ymx, pColor, prct);
    VtPlot(img, xmy, ypx, pColor, prct);
    VtPlot(img, xmx, ypy, pColor, prct);
}

template <class T>
void vt::VtDrawCircle(CTypedImg<T> &img, float x, float y, float r, const T *pColor, const RECT *prct)
{
    CRect rctclip = img.Rect();
    if(prct)
        rctclip = img.ClipRect(prct);

    float dx = 0;
    float dy = r;
    float rs = r*r;
    while(dy>=dx)
    {
        CirclePoints(img, x, y, dx, dy, pColor, &rctclip);
        dx+=1;
        float fd = rs - dx*dx;
        if(fd<0)
        {
            dx = 0;
            break;
        }
        dy = sqrt(fd);
    }
    if(dx!=0)
        CirclePoints(img, x, y, dx, dy, pColor, &rctclip);
}

template void vt::VtDrawCircle(CTypedImg<float> &img, float x, float y, float r,
                               const float *pColor, const RECT *prct);
template void vt::VtDrawCircle(CTypedImg<Byte> &img, float x, float y, float r,
                               const Byte *pColor, const RECT *prct);

template <class T>
void vt::VtDrawEllipse(CTypedImg<T> &img, float x, float y, CMtx2x2f mInvCov, float r, const T *pColor, const RECT *prct)
{
    // rather inefficient ellipse drawing - needs improvement
    CRect rctclip = img.Rect();
    if(prct)
        rctclip = img.ClipRect(prct);

    // make covariance matrix
    float det = mInvCov(0,0) * mInvCov(1,1) - mInvCov(0,1) * mInvCov(0,1);
    if(det==0)
        return;
    float xx = mInvCov(1,1) / det;
    float yy = mInvCov(0,0) / det;
    float xy = -mInvCov(0,1) / det;

    // 2x2 cholesky decompose U'U=C
    if(xx<=0)
        return;
    float u0 = sqrt(xx);
    float u1 = xy/u0;
    float tmp = yy - u1*u1;
    float u3 = sqrt(VtMax(tmp, 0.f));
    u0 *= r;
    u1 *= r;
    u3 *= r;

    float fm = u0;
    fm = VtMax(fm, fabs(u1));
    fm = VtMax(fm, u3);
    fm = VtMax(fm, 1.f);
    int count = 10 * (int)fm;
    int i;
    for(i=0; i<=count; i++)
    {
        float arg = ((float)VT_PI * 2 * i)/count;
        float xp = sin(arg);
        float yp = cos(arg);
        VtPlot(img, x + u0 * xp, y + u1 * xp + u3 * yp, pColor, &rctclip);
    }
}

template void vt::VtDrawEllipse(CTypedImg<float> &img, float x, float y, CMtx2x2f mInvCov, float r, 
                               const float *pColor, const RECT *prct);
template void vt::VtDrawEllipse(CTypedImg<Byte> &img, float x, float y, CMtx2x2f mInvCov, float r, 
                               const Byte *pColor, const RECT *prct);

