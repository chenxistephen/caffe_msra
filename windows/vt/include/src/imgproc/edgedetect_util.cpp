//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Utility routines for visualizing and dumping out edge/curve/line/ellipse detection results
//
//  History:
//      2012/02/05-szeliski
//          Split off from main EdgeDetectTest.cpp file, for easier re-use
//
//------------------------------------------------------------------------

#include "stdafx.h"
#include <errno.h>
#include <memory>
#include <stdio.h>
#include "vt_edgedetect_util.h"

using namespace vt;


// Forward declarations for local heper functions
void StackImages(int blockSize, vt::CRGBImg& I, vt::vector<std::shared_ptr<vt::CRGBImg> >& images);
RGBPix blend(const RGBPix& c, RGBPix& d, float fd);
RGBPix pseudoColor(int i);
LineSegment ExtendLine(CVec2d& p0, CVec2d& p1, double s);
HRESULT BresenhamLineTraverse(vector<int>& px, vector<int>& py, int p1x, int p1y, int p2x, int p2y);



HRESULT vt::VtDumpEdgels(CRGBImg& dst, const CRGBImg& src, float attenuation, const vector<EdgeSegment> &el, bool useZeroCrossings)
{
    VT_HR_BEGIN()

        VT_HR_EXIT( VtScaleImage(dst, src, attenuation) );
    int w = src.Width();
    int h = src.Height();
    for (int k=0; k<(int)el.size(); k++) 
    {
        RGBPix tag(255,255,0);	// yellow
        const EdgeSegment &edgel = el[k];
        int x, y;
        if(useZeroCrossings)
        {
            x = (int)edgel.x;
            y = (int)edgel.y;
        }
        else
        {
            x = (int)(edgel.x+0.5f);
            y = (int)(edgel.y+0.5f);
        }

        if (x>=0 && x<w && y>=0 && y<h) 
            dst(x,y) = tag;
    }

    VT_HR_END()
}



HRESULT vt::VtDumpCurves(CRGBImg& dst, const CRGBImg& src, float attenuation, const vector< vector<int> > &curves, const vector<EdgeSegment> &el, bool useZeroCrossings)
{
    VT_HR_BEGIN()

        VT_HR_EXIT( VtScaleImage(dst, src, attenuation) );
    int w = src.Width();
    int h = src.Height();

    for (int k=0; k<(int)curves.size(); k++) 
    {
        RGBPix tag = pseudoColor(k);
        const vector<int>& line = curves[k];
        for (int j=0; j<(int)line.size(); j++) 
        {
            const EdgeSegment &edgel = el[line[j]];
            int x, y;
            if(useZeroCrossings)
            {
                x = (int)edgel.x;
                y = (int)edgel.y;
            }
            else
            {
                x = (int)(edgel.x+0.5f);
                y = (int)(edgel.y+0.5f);
            }
            if (x>=0 && x<w && y>=0 && y<h) 
                dst(x,y) = tag;
        }
    }
    VT_HR_END()
}



HRESULT vt::VtDumpLines(CRGBImg& dst, const CRGBImg& src, float attenuation, const vector<LineSegment>& lines)
{
    VT_HR_BEGIN()
        VT_HR_EXIT( VtScaleImage(dst, src, attenuation) );
    int w = src.Width();
    int h = src.Height();
    RGBPix  white, bkg;
    white.r = white.g = white.b = 255;
    bkg.r   = bkg.g   = bkg.b   = 0; 
    for (int k=0;k<(int)lines.size();k++)
    {
        const LineSegment& lseg = lines[k];
        int p1x, p1y, p2x, p2y;
        p1x = int(lseg.start[0]); p1y = int(lseg.start[1]);
        p2x = int(lseg.end[0]); p2y = int(lseg.end[1]);

        vector<int> px, py;
        BresenhamLineTraverse( px, py, p1x, p1y, p2x, p2y);
        for (int kk=0;kk<(int)px.size();kk++)
            if (px[kk]>=0 && px[kk]<=w-1 && py[kk]>=0 && py[kk]<=h-1)
                dst(px[kk],py[kk]) = white;
    }
    VT_HR_END()
}


HRESULT vt::VtDumpEllipses(CRGBImg& dst, const CRGBImg& src, float attenuation, 
    const vector<EdgeSegment> &edgel_list,  const vector<EllipseSegment> &ellipses, bool useZeroCrossings)
{
    VT_HR_BEGIN()
        VT_HR_EXIT( VtScaleImage(dst, src, attenuation) );
    int w = src.Width();
    int h = src.Height();
    for (int k=0, n = (int)ellipses.size(); k<n; k++) 
    {
        EllipseSegment ell = ellipses[k];
        if (ell.n_points > 0)
        {
            // Generate synthetic points
            vector<EdgeSegment> edgels;
            vector<int> l;
            ell.GeneratePoints(0, edgels, l);

            // Draw the fitted ellipse
            RGBPix tag = pseudoColor(k);
            for (int j=0; j<(int)l.size(); j++) 
            {
                EdgeSegment &edgel = edgels[l[j]];
                int x, y;
                if(useZeroCrossings)
                {
                    x = (int)edgel.x;
                    y = (int)edgel.y;
                }
                else
                {
                    x = (int)(edgel.x+0.5f);
                    y = (int)(edgel.y+0.5f);
                }

                if (x>=0 && x<w && y>=0 && y<h) 
                    dst(x,y) = tag;
            }

            // Brighten the inlier points
            RGBPix white(255, 255, 255);
            float bfac = 0.3f;   // brightening factor
            vector<int> &filtered = ell.points;
            for (int j=0; j<(int)filtered.size(); j++) 
            {
                const EdgeSegment &edgel = edgel_list[filtered[j]];
                int x, y;
                if(useZeroCrossings)
                {
                    x = (int)edgel.x;
                    y = (int)edgel.y;
                }
                else
                {
                    x = (int)(edgel.x+0.5f);
                    y = (int)(edgel.y+0.5f);
                }

                if (x>=0 && x<w && y>=0 && y<h) 
                    dst(x,y) = blend(src(x,y), white, bfac);
            }
        }
    }
    VT_HR_END()
}



HRESULT vt::VtDumpBitangents(CRGBImg& dst, LineSegment bt[2])
{
    VT_HR_BEGIN()
    // Draw the lines, extended by a constant factor
    vector<LineSegment> lines;
    const double extend = 4;    
    VT_HR_EXIT( lines.push_back(ExtendLine(bt[0].start, bt[0].end  , extend)) );
    VT_HR_EXIT( lines.push_back(ExtendLine(bt[1].start, bt[1].end  , extend)) );
    VT_HR_EXIT( lines.push_back(ExtendLine(bt[0].start, bt[1].start, extend)) );
    VT_HR_EXIT( lines.push_back(ExtendLine(bt[0].end,   bt[1].end,   extend)) );
    VT_HR_EXIT( VtDumpLines(dst, dst, 1.0f, lines) );
    VT_HR_END()
}

HRESULT vt::VtDumpVanishingPoints(CRGBImg& dst, const CRGBImg& src, int& numValid, float attenuation, const vector<VanishingPoint>& vanishingPoints, const vector<LineSegment>& lineSegments)
{
    VT_HR_BEGIN()

    int w = src.Width();
    int h = src.Height();
    RGBPix  white, bkg;
    white.r = white.g = white.b = 255;
    bkg.r   = bkg.g   = bkg.b   = 0;

    vt::vector<std::shared_ptr<vt::CRGBImg> > vpOutputImages;

    int nVPs = (int) vanishingPoints.size();
    numValid = 0;
    for (int j = 0; j < nVPs; j++)
    {
        float minLength = 0.f;
        if(vanishingPoints[j].score > 0)
        {
            std::shared_ptr<CRGBImg> tmpPtr = std::make_shared<CRGBImg>();
            VT_HR_EXIT( VtScaleImage(*tmpPtr, src, attenuation) );
            RGBPix pixcolor;
            if (j==0)
            {
                pixcolor.r = 255; pixcolor.g = pixcolor.b = 0;
            }
            else if (j==1)
            {
                pixcolor.g = 255; pixcolor.r = pixcolor.b = 0;
            }
            else if (j==2)
            {
                pixcolor.b = 255; pixcolor.r = pixcolor.g = 0;
            }
            else if (j==3)
            {
                pixcolor.b = pixcolor.r = 255; pixcolor.g = 0;
            }
            else if (j==4)
            {
                pixcolor.g = pixcolor.r = 255;  pixcolor.b = 0;
            }

            for (int k = 0; k < (int)vanishingPoints[j].idx.size(); k++)
            {
                int vpi = vanishingPoints[j].idx[k];
                const LineSegment& lseg = lineSegments[vpi];
                int p1x, p1y, p2x, p2y;
                p1x = int(lseg.start[0]); p1y = int(lseg.start[1]);
                p2x = int(lseg.end[0]); p2y = int(lseg.end[1]);
                vector<int> px, py;
                BresenhamLineTraverse( px, py, p1x, p1y, p2x, p2y);
                for (int kk=0;kk<(int)px.size();kk++)
                    if (px[kk]>=0 && px[kk]<=w-1 && py[kk]>=0 && py[kk]<=h-1)
                        (*tmpPtr)(px[kk],py[kk]) = pixcolor;
                if(lineSegments[vanishingPoints[j].idx[k]].length > minLength)
                    minLength = lineSegments[vanishingPoints[j].idx[k]].length;
            }

            VT_HR_EXIT( vpOutputImages.push_back(tmpPtr) );
            numValid++;
        }
    }
    // Produces a stack of images each for a single vanishing point direction
    StackImages(3, dst, vpOutputImages); 

    VT_HR_END()
}



LineSegment ExtendLine(CVec2d& p0, CVec2d& p1, double s)
{
    // Lengthen the line by a factor of s
    LineSegment l;
    CVec2d v = p1 - p0;
    l.start = p0 - s * v;
    l.end   = p1 + s * v;
    l.length = (float) ((1 + s) * hypot(v[0], v[1]));
    return l;
}


RGBPix pseudoColor(int i)
{
    // De-interleave index into 3 channel bit patterns, reverse, and subtract from white
    RGBPix c(255,255, 255);
    unsigned char *b = (unsigned char *) &c;
    for (int j = 0; j < 3; j++)
    {
        int v = 128;
        for (int k = i >> j; k > 0; k = k >> 3, v = v >> 1)
        {
            b[j] -= static_cast<unsigned char>( (k & 1) * v );
        }
    }
    return c;
}


RGBPix blend(const RGBPix& c, RGBPix& d, float fd)
{
    float fc = 1.0f - fd;
    unsigned char r = (unsigned char) (fc*c.r + fd*d.r + 0.5f);
    unsigned char g = (unsigned char) (fc*c.g + fd*d.g + 0.5f);
    unsigned char b = (unsigned char) (fc*c.b + fd*d.b + 0.5f);
    return RGBPix(r, g, b);
}




void StackImages(int blockSize, vt::CRGBImg& I, vt::vector<std::shared_ptr<vt::CRGBImg> >& images)
{
    int n = (int) images.size();
    int nRows = int(ceil(double(n)/double(blockSize)));
    int nCols = blockSize;     if (n < nCols)      nCols = n;
    int wMax = 0;
    int hMax = 0;
    for (int i=0;i<n;i++) 
    {
        wMax = VtMax(wMax, images[i]->Width());
        hMax = VtMax(hMax, images[i]->Height());
    }

    int W = wMax * nCols;
    int H = hMax * nRows;

    RGBPix blk(0,0,0);
    I.Create(W,H); I.Fill(blk);
    int xOff, yOff;
    for (int i=0;i<n;i++) 
    {
        int by = i / blockSize;
        int bx = i % blockSize;
        xOff = bx * wMax;
        yOff = by * hMax;
        I.Paste(xOff,yOff, *(images[i]));
    }
    return;
}


HRESULT
    BresenhamLineTraverse(vector<int>& px, vector<int>& py,
    int p1x, int p1y, int p2x, int p2y)
{
    VT_HR_BEGIN()
        int d, x, y, ax, ay, sx, sy, dx, dy;
    bool keepon = true;

    dx = p2x - p1x; dy = p2y - p1y;
    ax = abs(dx)*2; ay = abs(dy)*2;
    sx = (dx>0) ? 1 : (-1);
    sy = (dy>0) ? 1 : (-1);

    x = p1x; y = p1y;
    if (ax>ay) {
        // x dominant
        d = ay - ax/2;
        while (keepon) {
            // do something at (x,y)
            VT_HR_EXIT( px.push_back(x) );
            VT_HR_EXIT( py.push_back(y) );
            if (x==p2x) { keepon = false; continue; }
            if (d>=0) { y += sy; d -= ax; }
            x += sx; d += ay;
        }
    }
    else {
        // y dominant
        d = ax - ay/2;
        while (keepon) {
            // do something at (x,y)
            VT_HR_EXIT( px.push_back(x) );
            VT_HR_EXIT( py.push_back(y) );
            if (y==p2y) { keepon = false; continue; }
            if (d>=0) { x += sx; d -= ay; }
            y += sy; d += ax;
        }
    }
    VT_HR_END()
}

void vt::VtDumpEdgels(const TCHAR* filename, const vector<EdgeSegment> &edgelList, bool Matlab, bool extraColumns)
{
    FILE *stream = NULL;
    _wfopen_s(&stream, filename, L"w");
    if(NULL != stream)
    {
        if (Matlab)
        {
        }
        else // tab separated text for Excel
        {
            fwprintf(stream, L"x\ty\tn_x\tn_y\tlength\tstrength\tphi%s\n",
                extraColumns ? L"\ttheta\tcurvature" : L"");
        }
        for (int k=0, n = (int) edgelList.size(); k < n; k++) 
        {
            const EdgeSegment &e = edgelList[k];
            fwprintf(stream, L"%s%.2f\t%.2f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f", 
                (Matlab) ? L"[ " : L"",
                e.x, e.y, e.n_x, e.n_y, e.length, e.strength, e.phi);
            if (extraColumns)
                fwprintf(stream, L"\t%.2f\t%.4f", 
                e.theta, e.curvature);

            fwprintf(stream, L"%s\n", (Matlab) ? L" ]" : L"");
        }
        fclose(stream);
    }
    return;
}


void vt::VtDumpCurves(const TCHAR* filename, const vector<EdgeSegment> &edgelList, const vector< vector<int> > &curves, bool Matlab)
{
    FILE *stream = NULL;
    _wfopen_s(&stream, filename, L"w");
    if(NULL != stream)
    {
        if (Matlab)
        {
        }
        else // tab separated text for Excel
        {
            fwprintf(stream, L"x\ty\tn_x\tn_y\tlength\tstrength\tphi\n");
        }
        for (int k=0, n = (int) curves.size(); k < n; k++) 
        {
            const vector<int> &line = curves[k];
            for (int j=0, m = (int)line.size(); j < m; j++) 
            {
                const EdgeSegment &e = edgelList[line[j]];
                fwprintf(stream, L"%s%.2f\t%.2f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f%s\n", 
                    (Matlab) ? L"[ " : L"",
                    e.x, e.y, e.n_x, e.n_y, e.length, e.strength, e.phi, e.curvature,
                    (Matlab) ? L" ]" : L"");
            }
            fwprintf(stream, (Matlab) ?  L"\n" : L"\n");
        }
        fclose(stream);
    }
    return;
}


void vt::VtDumpEllipses(const TCHAR* filename, const vector<EllipseSegment> &ellipses, bool Matlab)
{
    FILE *stream = NULL;
    _wfopen_s(&stream, filename, L"w");
    if(NULL != stream)
    {
        if (Matlab)
        {
        }
        else // tab separated text for Excel
        {
            fwprintf(stream, L"cx\tcy\tax\tay\ttheta\tmin_t\tmax_t\tdiff_t\tstrength\tdensity\tscore\tn_points\n");
        }
        for (int k=0, n = (int) ellipses.size(); k < n; k++) 
        {
            const EllipseSegment &e = ellipses[k];
            if (e.n_points > 0)
            {
                fwprintf(stream, L"%s%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d%s\n",
                    (Matlab) ? L"[ " : L"",
                    e.cx, e.cy, e.ax, e.ay, e.theta, e.min_t, e.max_t, e.max_t-e.min_t,
                    e.strength, e.density, e.score, e.n_points,
                    (Matlab) ? L" ]" : L"");
            }
        }
        fclose(stream);
    }
    return;
}


void vt::VtDumpBitangents(const TCHAR* filename, const LineSegment bitangents[2], const  CVec3f vanishingPoints[2], bool Matlab)
{
    FILE *stream = NULL;
    _wfopen_s(&stream, filename, L"w");
    if(NULL != stream)
    {
        if (Matlab)
        {
        }
        else // tab separated text for Excel
        {
            fwprintf(stream, L"px/vx\tpy/vy\tvz\n");
        }
        for (int k=0; k < 2; k++) 
        {
            fwprintf(stream, L"%.2f\t%.2f\n", bitangents[k].start[0], bitangents[k].start[1]);
            fwprintf(stream, L"%.2f\t%.2f\n", bitangents[k].end  [0], bitangents[k].end  [1]);
        }
        for (int l=0; l < 2; l++) 
        {
            CVec3f u = vanishingPoints[l].Unit();
            const double z_min = 1.0e-6;
            if (fabs(u[3]) > z_min)
                u = u / u[2];
            fwprintf(stream, L"%.6f\t%.6f\t%.6f\n", u[0], u[1], u[2]);
        }
        fclose(stream);
    }
    return;
}
