#include "stdafx.h"

#include "completion_tiles.h"

using namespace vt;

// constants
const int   Tiles::histoBins = 16;
const float Tiles::transportCost = 0.01f;
const int   Tiles::minValidCells = 32;
const float Tiles::validThresh = 0.05f;
const float Tiles::distanceScale = 1.f;
const float Tiles::epsilon = 0.001f;

typedef CVec2<int> CVec2i;
typedef CVec3<int> CVec3i;

static int qcompVec3i(void * /*context*/, const void *elem1, const void *elem2)
{
    int i = ((CVec3i *) elem1)->x - ((CVec3i *) elem2)->x;
    if (i == 0)
        i = ((CVec3i *) elem1)->y - ((CVec3i *) elem2)->y;
    if (i == 0)
        i = ((CVec3i *) elem1)->z - ((CVec3i *) elem2)->z;
    return i;
}

HRESULT Tiles::process(IImageReader * pPyr, int startLevel,
    int & tileLevel, vector<CFloatImg> & restriction, vector<CRect> & tiles)
{
    VT_HR_BEGIN()

    VT_DEBUG_LOG( "  computing restricted search spaces...\n" );

#if defined(_DEBUG) || defined(_HRPRINT) || defined(DEBUG_OUTPUT)
    CTimer timer;
#endif
#if defined(_DEBUG) || defined(_HRPRINT)
    float startTime = timer.GetTimeMilliSec() * 0.001f;
#endif

	tileLevel = VtClamp(startLevel - levelOffset, 0, maxLevel);

    VT_DEBUG_LOG( "    processing tiles at level %d\n", tileLevel );

    VT_HR_EXIT( pPyr->ReadImg(m_img, tileLevel) );

    int w = (m_img.Width()  + blockSize - 1) / blockSize;
    int h = (m_img.Height() + blockSize - 1) / blockSize;

    // compute tile types
    VT_HR_EXIT( m_type.Create(w, h, 1) );
    VT_HR_EXIT( m_type.Clear() );

    for (int y = 0; y < m_img.Height(); y++)
    {
        for (int x = 0; x < m_img.Width(); x++)
        {
            int cx = x/blockSize, cy = y/blockSize;
            byte * ptr = &m_type(cx, cy);

            if (!m_img(x, y, 3))
            {
                *ptr |= HOLE_TILE;
            }
            else
            {
                *ptr |= INTERIOR_TILE;
            }
        }
    }

    // compute histograms
#ifdef DEBUG_OUTPUT
    float dbgStartTime, dbgEndTime;

    dbgStartTime = timer.GetTimeMilliSec() * 0.001f;
#endif DEBUG_OUTPUT

    VT_HR_EXIT( computeHistograms() );

#ifdef DEBUG_OUTPUT
    dbgEndTime = timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "RESTR --- histo: %f\n", dbgEndTime-dbgStartTime );

    dbgStartTime = dbgEndTime;
#endif DEBUG_OUTPUT

    VT_HR_EXIT( computeAffinities() );

#ifdef DEBUG_OUTPUT
    dbgEndTime = timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "RESTR --- aff: %f\n", dbgEndTime-dbgStartTime );

    dbgStartTime = dbgEndTime;
#endif DEBUG_OUTPUT

    // find boundary & hole tiles
    int numBoundaryTiles = 0, numTiles = 0;

    CIntImg bdIdx;
    VT_HR_EXIT( bdIdx.Create(w, h, 1) );
    VT_HR_EXIT( bdIdx.Fill(-1) );

    vector<CVec2i> boundaryTiles;

    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            numTiles++;
            if (m_type(x, y) == BOUNDARY_TILE)
            {
                bdIdx(x, y) = numBoundaryTiles++;
                VT_HR_EXIT( boundaryTiles.push_back(CVec2i(x, y)) );
            }
        }
    }

    // compute segments
    VT_HR_EXIT( m_segments.Create(w, h, numBoundaryTiles, alignAny) );
    VT_HR_EXIT( m_segments.Fill(FLT_MAX) );

    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            if (m_type(x, y) != BOUNDARY_TILE)
            {
                continue;
            }

            int tile = bdIdx(x, y);

            VT_HR_EXIT( geodesicDistance(tile, x, y) );
        }
    }

#ifdef DEBUG_OUTPUT
    dbgEndTime = timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "RESTR --- segm: %f\n", dbgEndTime-dbgStartTime );

    dbgStartTime = dbgEndTime;
#endif DEBUG_OUTPUT

    //
    // compute restricted search space
    //
    VT_HR_EXIT( restriction.resize(numTiles) );
    for (int i = 0; i < numTiles; i++)
        VT_HR_EXIT( restriction[i].Create(w, h, 1, alignAny) );
    VT_HR_EXIT( tiles.resize(numTiles) );
    for (int i = 0; i < numTiles; i++)
        tiles[i] = CRect(LONG_MAX, LONG_MAX, 0, 0);

    for (int y = 0; y < h; y++)
    {
        vector<CVec3i> neighbors;

        for (int x = 0; x < w; x++)
        {
            int tile = w * y + x;

            const float distThresh = 1.5f;
            const int minUnionSize = 7;

            int minDist = INT_MAX, minX = -1, minY = -1;
            neighbors.clear();
            for (int i = 0; i < (int) boundaryTiles.size(); i++)
            {
                const CVec2i & bd = boundaryTiles[i];

                int dist = sqr(x-bd.x) + sqr(y-bd.y);
				if (dist < minDist)
				{
					minDist = dist;
					minX = bd.x;
					minY = bd.y;
				}

                VT_HR_EXIT( neighbors.push_back(CVec3i(dist, bd.x, bd.y)) );
            }

            float actualDistThresh = ceilf(minDist * sqr(distThresh));

            qsort_s(&neighbors[0], neighbors.size(), sizeof(CVec3i), qcompVec3i, NULL);

            vector<float *> segPtr;
            for (int i = 0; i < (int) neighbors.size(); i++)
            {
                int dist = neighbors[i].x;
                if (i >= minUnionSize && dist > actualDistThresh)
                {
                    break;                
                }

                int nx = neighbors[i].y;
                int ny = neighbors[i].z;
                
                int nidx = bdIdx(nx, ny);

                VT_HR_EXIT( segPtr.push_back(m_segments.Ptr(0, 0, nidx)) );
            }

            float * resPtr = restriction[tile].Ptr();
            int segStride = m_segments.Bands();

			for (int ry = 0; ry < h; ry++)
			{
				for (int rx = 0; rx < w; rx++)
				{
                    if (rx < minX-maxRestrictExtent || rx > minX+maxRestrictExtent || ry < minY-maxRestrictExtent || ry > minY+maxRestrictExtent)
                    {
                        for (int j = 0; j < (int)segPtr.size(); j++)
                        {
                            segPtr[j] += segStride;
                        }

                        *resPtr++ = FLT_MAX;
                    }
                    else
                    {
						float val = *segPtr[0];
						segPtr[0] += segStride;

						for (int j = 1; j < (int)segPtr.size(); j++)
						{
							if (*segPtr[j] < val)
							{
								val = *segPtr[j];
							}
							segPtr[j] += segStride;
						}

						*resPtr++ = val;

                        if (val != FLT_MAX)
						{
							tiles[tile].left   = VtMin(tiles[tile].left,   (long) rx);
							tiles[tile].top    = VtMin(tiles[tile].top,    (long) ry);
							tiles[tile].right  = VtMax(tiles[tile].right,  (long) rx + 1);
							tiles[tile].bottom = VtMax(tiles[tile].bottom, (long) ry + 1);
						}
					}
				}
            }
        }
    }

#ifdef DEBUG_OUTPUT
    dbgEndTime = timer.GetTimeMilliSec() * 0.001f;

    VT_DEBUG_LOG( "RESTR --- restr: %f\n", dbgEndTime-dbgStartTime );
#endif DEBUG_OUTPUT

    VT_DEBUG_LOG( "  [done] took %.2fs\n", timer.GetTimeMilliSec() * 0.001f - startTime );

    VT_HR_END()
}

HRESULT Tiles::computeHistograms()
{
    VT_HR_BEGIN()

    const int numChannels = histoBins * 3;

    int w = m_type.Width();
    int h = m_type.Height();

    VT_HR_EXIT( m_histo.Create(w, h, numChannels) );
    VT_HR_EXIT( m_histo.Clear() );

    CFloatImg wsum;
    VT_HR_EXIT( wsum.Create(w, h, 1) );
    VT_HR_EXIT( wsum.Clear() );

    for (int y = 0; y < m_img.Height(); y++)
    {
        for (int x = 0; x < m_img.Width(); x++)
        {
            if (!m_img(x, y, 3))
            {
                continue;
            }

            for (int ch = 0; ch < 3; ch++)
            {
                int bin = ch * histoBins + clamp(m_img(x, y, ch) * histoBins / 255, 0, histoBins - 1);
                m_histo(x / blockSize, y / blockSize, bin)++;
            }

            wsum(x / blockSize, y / blockSize)++;
        }
    }

    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            if (wsum(x, y) > 0)
            {
                for (int bin = 0; bin < numChannels; bin++)
                {
                    m_histo(x, y, bin) /= wsum(x, y);
                }
            }
        }
    }

    VT_HR_END()
}

float Tiles::histoDistEmdScan(float * h0, float * h1)
{
    float dist = 0;

    const float epsilon = 1e-4f;

    for (int ch = 0; ch < 3; ch++)
    {
        int cpos = 0;
        float cval = h0[ch*histoBins+0];

        for (int qpos = 0; qpos < histoBins; qpos++)
        {
            float togo = h1[ch*histoBins+qpos];
            if (togo < epsilon)
            {
                continue;
            }

            while (togo > epsilon)
            {
                if (cval >= togo-epsilon)
                {
                    dist += abs(cpos-qpos)*transportCost*togo;
                    cval -= togo;
                    break;
                }
                else if (cval > 0)
                {
                    dist += abs(cpos-qpos)*transportCost*cval;
                    togo -= cval;
                    cpos++;
//                    assert(cpos < histoBins);
                    if (cpos >= histoBins)
                    {
                        for (int i = 0; i < histoBins; i++)
                        {
                            VT_DEBUG_LOG( "bin %02d --> %f --> %f\n", i, h0[ch*histoBins+i], h1[ch*histoBins+i] );
                        }
                        VT_ASSERT( false );
                    }
                    cval = h0[ch*histoBins+cpos];
                }
                else
                {
                    cpos++;
//                    assert(cpos < histoBins);
                    if (cpos >= histoBins)
                    {
                        for (int i = 0; i < histoBins; i++)
                        {
                            VT_DEBUG_LOG( "bin %02d --> %f --> %f\n", i, h0[ch*histoBins+i], h1[ch*histoBins+i] );
                        }
                        VT_ASSERT( false );
                    }
                    cval = h0[ch*histoBins+cpos];
                }
            }
        }
    }

    return dist;
}

HRESULT Tiles::computeAffinities()
{
    VT_HR_BEGIN()

    int w = m_type.Width();
    int h = m_type.Height();

    VT_HR_EXIT( m_affH.Create(w - 1, h, 1) );
    VT_HR_EXIT( m_affH.Clear() );
    VT_HR_EXIT( m_affV.Create(w, h - 1, 1) );
    VT_HR_EXIT( m_affV.Clear() );

    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            if ((m_type(x, y) & INTERIOR_TILE) == 0)
            {
                continue;
            }

            if (x < w - 1 && (m_type(x+1, y) & INTERIOR_TILE))
            {
                m_affH(x, y) = histoDistEmdScan(&m_histo(x, y), &m_histo(x+1, y));
                VT_ASSERT( _finite(m_affH(x, y)) && m_affH(x, y) >= 0 );
            }

            if (y < h - 1 && (m_type(x, y+1) & INTERIOR_TILE))
            {
                m_affV(x, y) = histoDistEmdScan(&m_histo(x, y), &m_histo(x, y+1));
                VT_ASSERT( _finite(m_affV(x, y)) && m_affV(x, y) >= 0 );
            }
        }
    }

    VT_HR_END()
}

struct FII
{
    FII(float dist0, int x0, int y0) : dist(dist0), x(x0), y(y0) {};
    float dist;
    int x, y;
};

static int qcompFII(void * /*context*/, const void *elem1, const void *elem2)
{
    float f = ((FII *) elem1)->dist - ((FII *) elem2)->dist;
    int i = f == 0 ? 0 : (f > 0 ? 1 : -1);
    if (i == 0)
        i = ((FII *) elem1)->x - ((FII *) elem2)->x;
    if (i == 0)
        i = ((FII *) elem1)->y - ((FII *) elem2)->y;
    return i;
}

HRESULT Tiles::thresholdDistances(const int tile)
{
    VT_HR_BEGIN()

    const int dw = m_histo.Width(), dh = m_histo.Height();

    vector<FII> cells;

    for (int y = 0; y < dh; y++)
    {
        for (int x = 0; x < dw; x++)
        {
            VT_HR_EXIT( cells.push_back(FII(m_segments(x, y, tile), x, y)) );
        }
    }

    qsort_s(&cells[0], cells.size(), sizeof(FII), qcompFII, NULL);

    for (int i = 0; i < (int) cells.size(); i++)
    {
        const FII & c = cells[i];

        const float dist = c.dist;
        const int x = c.x, y = c.y;

        if (i >= minValidCells && dist > validThresh)
        {
            m_segments(x, y, tile) = FLT_MAX;
        }
    }

    VT_HR_END()
}

HRESULT Tiles::geodesicDistance(const int tile, const int x0, const int y0)
{
    VT_HR_BEGIN()

    if ((m_type(x0, y0) & INTERIOR_TILE) == 0)
    {
        return hr;
    }

    m_segments(x0, y0, tile) = 0;

    struct E
    {
        E(float dist0, CVec2i& p0) : dist(dist0), p(p0) {};
        float dist;
        CVec2i p;
    };

    struct CompareElement
    {
        bool operator()(const E & lhs, const E & rhs)
        {
            return lhs.dist > rhs.dist;
        }
    };

    vector<E> Q;
    CVec2i v0(x0, y0);
    VT_HR_EXIT( Q.push_back(E(0.0f, v0)) );

    while (!Q.empty())
    {
        // find smallest dist in queue
        float myDist = FLT_MAX;
        int myI = 0;
        for (int i = (int) Q.size() - 1; i >= 0; i--)
        {
            if (Q[i].dist < myDist)
            {
                myDist = Q[i].dist;
                myI = i;
            }
        }
        CVec2i p = Q[myI].p;

        // swap smallest with back and pop
        Q[myI] = Q.back();
        Q.pop_back();

        if (myDist > m_segments(p.x, p.y, tile))
        {
            continue;
        }

        if (p.x > 0 && (m_type(p.x-1, p.y) & INTERIOR_TILE))
        {
            float newDist = myDist + m_affH(p.x-1, p.y) * distanceScale + epsilon;
            if (newDist < m_segments(p.x-1, p.y, tile))
            {
                m_segments(p.x-1, p.y, tile) = newDist;
                CVec2i v(p.x-1, p.y);
                VT_HR_EXIT( Q.push_back(E(newDist, v)) );
            }
        }

        if (p.y > 0 && (m_type(p.x, p.y-1) & INTERIOR_TILE))
        {
            float newDist = myDist + m_affV(p.x, p.y-1) * distanceScale + epsilon;
            if (newDist < m_segments(p.x, p.y-1, tile))
            {
                m_segments(p.x, p.y-1, tile) = newDist;
                CVec2i v(p.x, p.y-1);
                VT_HR_EXIT( Q.push_back(E(newDist, v)) );
            }
        }

        if (p.x < m_type.Width() - 1 && (m_type(p.x+1, p.y) & INTERIOR_TILE))
        {
            float newDist = myDist + m_affH(p.x, p.y) * distanceScale + epsilon;
            if (newDist < m_segments(p.x+1, p.y, tile))
            {
                m_segments(p.x+1, p.y, tile) = newDist;
                CVec2i v(p.x+1, p.y);
                VT_HR_EXIT( Q.push_back(E(newDist, v)) );
            }
        }

        if (p.y < m_type.Height() - 1 && (m_type(p.x, p.y+1) & INTERIOR_TILE))
        {
            float newDist = myDist + m_affV(p.x, p.y) * distanceScale + epsilon;
            if (newDist < m_segments(p.x, p.y+1, tile))
            {
                m_segments(p.x, p.y+1, tile) = newDist;
                CVec2i v(p.x, p.y+1);
                VT_HR_EXIT( Q.push_back(E(newDist, v)) );
            }
        }
    }

    VT_HR_EXIT( thresholdDistances(tile) );

    VT_HR_END()
}
