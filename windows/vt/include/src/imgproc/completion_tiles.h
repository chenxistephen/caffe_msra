#pragma once

#include "vtcommon.h"

namespace vt {

// #define USE_TILE_BBOXES

template <class TA, class TB, class TC> TA clamp(const TA & x, const TB & minv, const TC & maxv) { if (x < (TA)minv) return (TA)minv; else if (x > (TA)maxv) return (TA)maxv; return x; }
template <class T> T sqr(const T & x) { return x*x; }

class Tiles
{
public:
    HRESULT process(IImageReader * pPyr, int startLevel,
        int & tileLevel, vector<CFloatImg> & restriction, vector<CRect> & tiles);

    const static int levelOffset = 4;
    const static int blockSize = 1 << levelOffset;
    const static int maxLevel = 5;
    const static int maxRestrictExtent = 10;

private:
    const static int INTERIOR_TILE = 1;
    const static int HOLE_TILE = 2;
    const static int BOUNDARY_TILE = 3;

    HRESULT computeHistograms();
    HRESULT computeAffinities();
    float histoDistEmdScan(float * h0, float * h1);
    HRESULT thresholdDistances(const int tile);
    HRESULT geodesicDistance(const int tile, const int x0, const int y0);

    const static int histoBins;
    const static float transportCost, validThresh, epsilon;
    const static int minValidCells;
    const static float distanceScale;

    CFloatImg m_histo, m_affH, m_affV, m_segments;
    CByteImg m_img, m_type;
};
}
