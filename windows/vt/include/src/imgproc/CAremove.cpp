
#include "stdafx.h"

#include "vtimgproc.h"

using namespace vt;

#define N_LEVELS      3
#define MIN_STRENGTH  (0.14f/2.f)

//
// flags used in the edge image computation to filter spurious edges
// 
#define EDGE_HORI		1
#define EDGE_VERT		2

#if 0 // enable some diagnostics

#include "vtfileio.h"

void DebugDumpDoG(const CFloatImg& DoG)
{
	CByteImg imgDump(DoG.Width(), DoG.Height());
	for( int y = 0; y < DoG.Height(); y++ )
	{
		const float* pS = DoG.Ptr(y);
		Byte* pD = imgDump.Ptr(y);
		VtScaleOffsetSpan(pS, pD, 4.f, 128.f, DoG.Width());
	}
	VtGammaMap(imgDump, imgDump, eGammaSRGBLinearToNonLinear);
	VtSaveImage(L"c:\\temp\\dogimg.png", imgDump);
}

void DebugDumpEdge(const CByteImg& edge_img, WCHAR* base, int i)
{
	CByteImg imgDump(edge_img.Width(), edge_img.Height());
	edge_img.CopyTo(imgDump);
	VtGammaMap(imgDump, imgDump, eGammaSRGBLinearToNonLinear);

	vt::wstring_b<MAX_PATH> fname;
	fname.format(L"c:\\temp\\%s_%d.png",base,i);
	VtSaveImage(fname, imgDump);
}

#endif

///////////////////////////////////////////////////
// helper functions

inline float LengthNorm( int w, int h )
{
	int length = (w>h) ? w : h;
	
	return (length/2.0f);
}

template <class T>
T _interpcnv(float v)
{
	return F2I(v);
}

template <>
float _interpcnv(float v)
{
	return (float)v;
}

template <class T>
inline T _interpolate( const CTypedImg<T>& img, vt::CPoint c, 
					   const vt::CVec2f& frac )
{
    int w = img.Width(), h = img.Height();

	// deal with out-of-bounds conditions
	if (c.x < 0 )
	{ 
		if ( c.x > -1 ) 
			c.x = 0;
		else
			return 0;
	}
	else if( c.x > w-1 )
	{
		return 0;
	}
	if (c.y < 0)
	{ 
		if ( c.y > -1 ) 
			c.y = 0;
		else
			return 0;
	}
	else if( c.y > h-1 )
	{
		return 0;
	}

	// run interpolation
	const T *ptr00, *ptr01, *ptr10, *ptr11;
	ptr00 = img.Ptr(c.y) + c.x;
	ptr10 = (c.y >= h-1)? ptr00: img.Ptr(c.y+1)+c.x;
	if( c.x >= w-1 )
	{
		ptr01 = ptr00;
		ptr11 = ptr10;
	}
	else
	{
		ptr01 = ptr00+1;
		ptr11 = ptr10+1;
	}

	float c0, c1, c2, c3, d0, d1, d2;
    c0 = (float)ptr00[0], c1 = (float)ptr01[0]; 
	c2 = (float)ptr10[0], c3 = (float)ptr11[0];
  	d0 = c1 - c0, d1 = c2 - c0, d2 = c3 - c2 - d0;
    c0 += d0 * frac.x;
	d1 += d2 * frac.x; 
	c0 += d1 * frac.y; 

	return _interpcnv<T>(c0);
}

void 
ComputeRowWarpMap(vt::CPoint* pC, vt::CVec2f* pF, int w, int y,
				  const CVec2f& cntr, const CVecf& poly)
{
	VT_ASSERT( poly.Size() == 5 );
	float yd  = float(y) - cntr.y;
	float yd2 = yd*yd; 

#define WARP_MAP_SPAN 64

	float arfXAddr[WARP_MAP_SPAN];

	for( int s = 0; s < w; s+=WARP_MAP_SPAN )
	{
		int x = s;
        float* pScaleFactor = arfXAddr;
		int end = VtMin(s+WARP_MAP_SPAN, w);

#if (defined(_M_IX86) || defined(_M_AMD64))
		// compute the scale factor at each pixel - use SSE if supported
        BOOL bSSE2 = g_SupportSSE2();
		if( bSSE2 )
		{
			// setup a 4 way poly computation
			float p0[4] = {poly[0], poly[0], poly[0], poly[0]};
			float p1[4] = {poly[1], poly[1], poly[1], poly[1]};
			float p2[4] = {poly[2], poly[2], poly[2], poly[2]};
			float p3[4] = {poly[3], poly[3], poly[3], poly[3]};
			float p4[4] = {poly[4], poly[4], poly[4], poly[4]};
			__m128 x0 = _mm_set1_ps(4.f);
			__m128 x1 = _mm_set_ps(x-cntr.x+3, x-cntr.x+2, x-cntr.x+1, x-cntr.x+0);
			__m128 x2 = _mm_set1_ps(yd2);

			while( x+4 <= end )
			{
				__m128 x3 = x1;
				x3 = _mm_mul_ps(x3, x1);  // (x-cntr.x)^2
				x3 = _mm_add_ps(x3, x2);  // (x-cntr.x)^2 + (y-cntr.y)^2
				x3 = _mm_sqrt_ps(x3);     // sqrt => radius
										  
				__m128 x4 = x3;                // radius
				__m128 x5 =	_mm_loadu_ps(p1);  // poly[1]
				x5 = _mm_mul_ps(x5, x4);       // poly[1]*radius
	
				x4 = _mm_mul_ps(x4, x3);       // radius^2
				__m128 x6 = _mm_loadu_ps(p2);  // poly[2]
				x6 = _mm_mul_ps(x6, x4);       // poly[2]*radius^2
				x5 = _mm_add_ps(x5, x6);       // accumulate   
	
				x4 = _mm_mul_ps(x4, x3);       // radius^3
				x6 = _mm_loadu_ps(p3);         // poly[3]
				x6 = _mm_mul_ps(x6, x4);       // poly[3]*radius^3
				x5 = _mm_add_ps(x5, x6);       // accumulate   
	
				x4 = _mm_mul_ps(x4, x3);       // radius^4
				x6 = _mm_loadu_ps(p4);         // poly[4]
				x6 = _mm_mul_ps(x6, x4);       // poly[4]*radius^4
				x5 = _mm_add_ps(x5, x6);       // accumulate   
	
				x6 = _mm_loadu_ps(p0);         // poly[0]
				x5 = _mm_add_ps(x5, x6);       // accumulate   

				_mm_storeu_ps(pScaleFactor, x5);       // store the scale factor

				x1 = _mm_add_ps(x1, x0);       // increment addresses
	
				x+=4;
				pScaleFactor += 4;
			}
		}
#endif

		// compute the scale factor at each pixel - if SSE is supported then
		// this loop simply finishes of any values for spans that aren't a
		// multiple of 4, for non SSE case this does all of the work
		for( ; x < end; x++, pScaleFactor++ )
		{
			float xd = float(x)-cntr.x;
			float radius = sqrtf(xd*xd+yd2);
	
			float rpow   = radius;
			float factor = (float)poly[0];
			for (int i=1; i<poly.Size(); i++)
			{
				factor += rpow*(float)poly[i];
				rpow *= radius;
			}
			*pScaleFactor = factor;
		}

		// compute the adjusted addresses - use SSE if supported
        pScaleFactor = arfXAddr;
		x = s;
#if (defined(_M_IX86) || defined(_M_AMD64))
		if( bSSE2 )
		{
	        float arfYAddr[WARP_MAP_SPAN];

            // 1) compute the float address in the source
			__m128 x0 = _mm_set1_ps(4.f);
			__m128 x1 = _mm_set_ps(x-cntr.x+3, x-cntr.x+2, x-cntr.x+1, x-cntr.x);
			__m128 x2 = _mm_set1_ps(yd);
			__m128 x3 = _mm_set1_ps(cntr.x);
			__m128 x4 = _mm_set1_ps(cntr.y);

			float* pXAddr = arfXAddr;
			float* pYAddr = arfYAddr;
			while( x+4 <= end )
			{
				__m128 x5 =	_mm_loadu_ps(pScaleFactor);  // load radial scale factor

				__m128 x6 = x1;                  // (x-cntr.x)*factor + cntr.x
				x6 = _mm_mul_ps(x6, x5);
				x6 = _mm_add_ps(x6, x3);   

				__m128 x7 = x2;                  // (y-cntr.y)*factor + cntr.y
				x7 = _mm_mul_ps(x7, x5);
				x7 = _mm_add_ps(x7, x4);   

				_mm_storeu_ps(pXAddr, x6);       // store the X address
				_mm_storeu_ps(pYAddr, x7);       // store the Y address

				x1 = _mm_add_ps(x1, x0);       // increment addresses

				x+=4;
				pScaleFactor += 4;
				pXAddr += 4;
				pYAddr += 4;
			}

			// 2) compute the int address and scale factor
			pXAddr = arfXAddr;
			pYAddr = arfYAddr;
			x = s;
			while( x+4 <= end )
			{
			    x0 = _mm_loadu_ps(pXAddr);            // load the x address
				__m128i x0i = _mm_cvttps_epi32(x0);   // int(x)
				x1 = _mm_cvtepi32_ps(x0i);
				x0 = _mm_sub_ps(x0, x1);              // x - float(int(x))

			    x2 = _mm_loadu_ps(pYAddr);            // load the y address
				__m128i x2i = _mm_cvttps_epi32(x2);   // int(y)
				x3 = _mm_cvtepi32_ps(x2i);
				x2 = _mm_sub_ps(x2, x3);              // y - float(int(y))

                __m128i x3i = _mm_unpackhi_epi32(x0i, x2i);  // interleave address
                x0i = _mm_unpacklo_epi32(x0i, x2i);

				_mm_storeu_si128((__m128i*)pC,     x0i); // store address
				_mm_storeu_si128((__m128i*)(pC+2), x3i);

                x3 = _mm_unpackhi_ps(x0, x2);        // interleave fractions
                x0 = _mm_unpacklo_ps(x0, x2);

				_mm_storeu_ps((float*)pF,     x0);   // store fractions
				_mm_storeu_ps((float*)(pF+2), x3);

				x+=4;
				pXAddr += 4;
				pYAddr += 4;
				pC += 4;
				pF += 4;
			}
		}
#endif

		// finish computing the adjusted addresses - if SSE supporte then this just
		// does any remaining due to non-multiple of 4 span, otherwise this
		// does all of the work
		for( ; x < end; x++, pC++, pF++, pScaleFactor++ )
		{
			float factor = *pScaleFactor;
			float xsrc = factor*(float(x)-cntr.x) + cntr.x;
			float ysrc = factor*(yd) + cntr.y;
			pC->x = (int)xsrc;
			pC->y = (int)ysrc;
			pF->x = xsrc-float(pC->x);
			pF->y = ysrc-float(pC->y);
		}
	}
}

// Compute radial warp only in designated areas
template <class T>
HRESULT RadialWarp( CTypedImg<T> &warped, const CTypedImg<T> &img, 
					const CVecd& mag_poly_norm )
{
	HRESULT hr = NOERROR;
	int i, r, c, w = img.Width(), h = img.Height();
	float length_norm = LengthNorm( w, h );
	CVec2f cntr(float(w)/2.f, float(h)/2.f);

	vt::vector<vt::CPoint> vecSrcCoords;
	vt::vector<vt::CVec2f> vecSrcCoef;
    CVecf mag_poly;

	VT_HR_EXIT( warped.Create(w,h) );
	VT_HR_EXIT( vecSrcCoords.resize(w) );
	VT_HR_EXIT( vecSrcCoef.resize(w) );
    VT_HR_EXIT( mag_poly.Create(mag_poly_norm.Size()) );

	// scale the poly coeffs for this particular image size
    mag_poly[0] = float(mag_poly_norm[0]);
	double lenpow = length_norm;
    for( i = 1; i < mag_poly.Size(); i++ )
	{
		mag_poly[i] = float(mag_poly_norm[i] / lenpow);
		lenpow *= length_norm;
	}

	//int n_mag_poly = mag_poly_norm.Size();
	for (r=0; r<h; r++)
	{
        vt::CPoint* pC = vecSrcCoords.begin();
        vt::CVec2f* pF = vecSrcCoef.begin();
		ComputeRowWarpMap(pC, pF, w, r, cntr, mag_poly);

		T *q = warped.Ptr(r);
		for (c=0; c<w; c++, q++, pC++, pF++)
		{
			*q = _interpolate( img, *pC, *pF );
		}
	}

Exit:
	return hr;
}

//#define MIN_DISTANCE_FROM_CENTER				10
#define MAX_COS_ANGLE							0.707

void
FilterEdgeList(IN OUT vector<EdgeSegment>& el)
{
#ifdef MIN_DISTANCE_FROM_CENTER
	int wr = 0;
	for (int i=0; i<(int)el.size(); i++)
	{
		const EdgeSegment &e = el[i];

		// remove edges too close to radial direction
		float c_x = w/2.0f, c_y = h/2.0f; // assume image center is CA center
		float dx = e.x - c_x, dy = e.y - c_y;
		float mag = sqrt( dx*dx + dy*dy );

		if (mag<MIN_DISTANCE_FROM_CENTER)
			continue;

		float cos_angle = (dx*e.n_x + dx*e.n_y)/mag;
		if (cos_angle<0) cos_angle = -cos_angle;
		if (cos_angle>=MAX_COS_ANGLE) 
			continue;

		if( i != wr )
		{
			VT_ASSERT( wr < i );
			el[wr] = e;
		}
		wr++;
	}
	el.resize(wr);
#else 
    UNREFERENCED_PARAMETER(el);
#endif
}

//+-----------------------------------------------------------------------------
// 
// Function: KeepAdjacentEdgels
// 
// Synposis: retain edges where R, G, and B are very close to each other
//           and remove the rest
//
//------------------------------------------------------------------------------
HRESULT
KeepAdjacentEdgels( IN OUT vector<EdgeSegment> (&el)[3], int w, int h, 
					bool bRemoveSpeckle = true )
{
    HRESULT hr;

	int i, j, r2, c2;
	CByteImg count;
	VT_HR_EXIT( count.Create(w,h) );
    count.Clear();

	// extract a mask for each pixel that indicates which edgels are present 
	for( i = 0; i < 3; i++ )
	{
		const vector<EdgeSegment>& elchan = el[i];
		int chanshift = 2 * i;
		for( j = 0; j < (int)elchan.size(); j++ )
		{
			// because many edgels have coords that fall extremely near the 0.5
			// it is better to do these operations on 'dual' coords - those
			// halfway between orig image pixels
			const EdgeSegment &e = elchan[j];
			int xdual = F2I(e.x-0.5f), ydual = F2I(e.y-0.5f);
			count.Pix(xdual, ydual) |= 
				((fabs(e.n_y)>fabs(e.n_x)) ? EDGE_VERT : EDGE_HORI) << chanshift;
		}
	}

	if ( bRemoveSpeckle )
	{
		// delete edgels that have no nearby neighbors of similar orientation
		// in this color channel
		for ( i = 0; i < 3; i++ )
		{
			int wr = 0;
			vector<EdgeSegment>& elchan = el[i];
			int chanshift = 2 * i;
			int chanmask  = (EDGE_VERT|EDGE_HORI) << chanshift;
			for (j=0; j < (int)elchan.size(); j++)
			{
				const EdgeSegment &e = elchan[j];
				int xdual = F2I(e.x-0.5f);
				int ydual = F2I(e.y-0.5f);
				BYTE *p0 = count.Ptr(ydual-1);
				BYTE *p1 = count.Ptr(ydual);
				BYTE *p2 = count.Ptr(ydual+1);
   
				// VT_ASSERT( p1[xdual] & chanmask ); // harmless to remove
				Byte adjflag = p0[xdual-1] | p0[xdual] | p0[xdual+1] |
							   p1[xdual-1] |             p1[xdual+1] |
							   p2[xdual-1] | p2[xdual] | p2[xdual+1];
				if ( (p1[xdual] & chanmask & adjflag ) == 0 )
				{
					p1[xdual] &= ~chanmask;
					continue;
				}

				if ( wr != j )
				{
					VT_ASSERT( wr < j );
					elchan[wr] = e;
				}
				wr++;
			}
			elchan.resize(wr);
		}
	}

	// delete edgels that don't have edgels from other color channels nearby
	for( i = 0; i < 3; i++ )
	{
		vector<EdgeSegment>& elchan = el[i];

		int wr = 0;
		for( j = 0; j < (int)elchan.size(); j++ )
		{
			const EdgeSegment &e = elchan[j];
			int xdual = F2I(e.x-0.5f), ydual = F2I(e.y-0.5f);

			BYTE status = 0;
			for (r2=ydual-3; r2<=ydual+3; r2++)
			{
				const Byte* pC = count.Ptr(r2);
				for (c2=xdual-3; c2<=xdual+3; c2++)
					status |= pC[c2];
			}

			// if all channels aren't present nearby, then this edgel gets
			// deleted
			if ( (status & (3<<0)) == 0 || 
				 (status & (3<<2)) == 0 || 
				 (status & (3<<4)) == 0 )
			{
				continue;
			}

			if( wr != j )
			{
				VT_ASSERT( wr < j );
				elchan[wr] = e;
			}
			wr++;
		}
		elchan.resize(wr);
	}

Exit:
    return hr;
}

class CFootprintSubTable
{
public:
	CFootprintSubTable(int iCenterFlr, int iHalfWidth, float* pCurSubMap) :
		m_iCenterFlr(iCenterFlr), m_iHalfWidth(iHalfWidth),
		m_pCurSubMap(pCurSubMap)
	{}
	 
	float Map(int c)
	{
		VT_ASSERT( c >= m_iCenterFlr-m_iHalfWidth && 
				   c <= m_iCenterFlr+m_iHalfWidth );
        return m_pCurSubMap[c - m_iCenterFlr + m_iHalfWidth];
	}

private:
	int    m_iCenterFlr;
	int    m_iHalfWidth;
	float* m_pCurSubMap;
};

class CFootprintMasterTable
{
public:
	HRESULT Initialize(int iHalfWidth, int iSubPixels = 256)
	{
		HRESULT hr;

		m_iSubPixMapSize = (2*iHalfWidth+1);
		m_iSubPixels     = iSubPixels;
		m_iHalfWidth     = iHalfWidth;

		VT_HR_EXIT( m_map.resize((iSubPixels+1)*m_iSubPixMapSize) );

		float fNorm = 1.f / float(iHalfWidth);
		for( int s = 0, i = 0; s <= iSubPixels; s++ )
		{
			float fSub = float(s)/float(iSubPixels);
			for( int x = -iHalfWidth; x <= iHalfWidth; x++, i++ )
			{
				// the edge weight in refine edges used an additional sqrt
				// so apply ^0.25 to avoid that calc.
				m_map[i] = 
					pow((iHalfWidth - fabs(fSub - float(x))) * fNorm, 0.25f);
			}
            VT_ASSERT( i <= (int)m_map.size() );
		}

    Exit:
        return hr;
	}

	CFootprintSubTable GetSubTableForCenter(float c)
	{
		float fFlr   = floor(c);
		int iSubPix  = F2I((c-fFlr)*float(m_iSubPixels));
		VT_ASSERT( iSubPix <= m_iSubPixels );

		return CFootprintSubTable(F2I(fFlr), m_iHalfWidth,
                                  m_map.begin() + iSubPix*m_iSubPixMapSize);
	}

private:
	vt::vector<float> m_map;
	int               m_iSubPixMapSize;
	int               m_iSubPixels;
	int               m_iHalfWidth;
};

HRESULT
ComputeEdgeFootprints(OUT CByteImg &edge_img, const vector<EdgeSegment>& el)
{
    HRESULT hr = NOERROR;

    int w = edge_img.Width(), h = edge_img.Height(); 

    // for splatting
	// note: count can be a Byte image unless foot_half^2 gets to be larger
	//       than 256, because there can be no more than one edgel per pixel
	int foot_half = 3, foot_half2 = foot_half + 1;
	CFootprintMasterTable fh;
	CFloatImg w_sum;
	CByteImg  count;

	// splat for better localization
	VT_HR_EXIT( w_sum.Create(w,h) );
	VT_HR_EXIT( count.Create(w,h) );
	VT_HR_EXIT( w_sum.Clear() );
	VT_HR_EXIT( count.Clear() );
	VT_HR_EXIT( fh.Initialize(foot_half2) );

	for (int i=0; i<(int)el.size(); i++)
	{
		const EdgeSegment &e = el[i];

		// have some footprint to facilitate tracking
		CFootprintSubTable footx = fh.GetSubTableForCenter(e.x);
		CFootprintSubTable footy = fh.GetSubTableForCenter(e.y);

		int x = F2I(e.x);
		int y = F2I(e.y);
		for (int r=y-foot_half; r<=y+foot_half; r++)
		{
			float dy2 = footy.Map(r);

			int c = x-foot_half;
			float* pw = w_sum.Ptr(c,r);
			Byte*  pc = count.Ptr(c,r);

			for (; c<=x+foot_half; c++, pw++, pc++)
			{
				float dx2 = footx.Map(c);
  
				VT_ASSERT( dx2 >=0 && dy2 >= 0 );
				VT_ASSERT( count.Pix(c,r) < 255 );

				*pw += dx2*dy2;
				(*pc)++;
                                        
				// test code to see how well the lookup table is working
				/*{
				    float dxt = (foot_half2 - fabs(c - e.x))/foot_half2;
					float dyt = (foot_half2 - fabs(r - e.y))/foot_half2;
					float dxyt = sqrt(dxt*dyt);
					VT_ASSERT( fabs(dxyt-dx2*dy2)/float(dxyt) < 0.005f );
				}*/
			}
		}
	}

	for (int r=1; r<h-1; r++)
	{
		BYTE *q = edge_img.Ptr(r);
		const float *p = w_sum.Ptr(r);
		const Byte  *n = count.Ptr(r);

		for (int c=1; c<w-1; c++)
		{
			q[c] = n[c]? (Byte) F2I(255.f * p[c]/float(n[c])) : 0;
		}
	}

Exit:
	return hr;
}

///////////////////////////////////////////////////

CChromaticAberrationCorrectFast::CChromaticAberrationCorrectFast(void)
{
	m_width = m_height = 0;
	m_reference = REFERENCE_DEF;
	m_bmaskcomputed = false;
	m_baligned = false;

	m_spread = 5;
	m_margin = 2*m_spread; // area for evaluation to avoid image boundary effects
}

CChromaticAberrationCorrectFast::~CChromaticAberrationCorrectFast(void)
{
}

HRESULT 
CChromaticAberrationCorrectFast::SetImage( CFloatImg *fchannels, int reference )
{
	HRESULT hr = NOERROR;

	m_fchannels = fchannels;
	m_width = fchannels[0].Width();
	m_height = fchannels[0].Height();
	m_c[1] = (float)m_height/2.0f;
	m_c[0] = (float)m_width/2.0f;

	m_reference = reference;

	VT_DEBUG_LOG( "VT_LOG_CA: SetImage: w = %d, h = %d, spread = %d\n", 
		m_width, m_height, m_spread);

	for (int i=0; i<3; i++)
	{
		VT_HR_EXIT( m_fchannels_aligned[i].Create( m_width, m_height ) );
	}

Exit:
	return hr;
}

HRESULT
CChromaticAberrationCorrectFast::ComputeEdgeImg( )
{
	HRESULT hr = NOERROR;

    int w = m_fchannels[0].Width(), h = m_fchannels[0].Height();

	EdgeDetectParams edParams;
	edParams.numLevelsForDoG = N_LEVELS;
	edParams.gradientMagnitudeThreshold = MIN_STRENGTH/2;

	// extract the edge list for each channel
	vector<EdgeSegment> el[3];
	for( int i = 0; i < 3; i++ )
	{
		VT_DEBUG_LOG( "VT_LOG_CA: Debug: edge channel %d\n", i ); 

		//CFloatImg DoG;
		//VT_HR_EXIT( GetDifferenceOfGaussianImage( DoG, m_fchannels[i], N_LEVELS ) );
		//VT_DEBUG_LOG( "VT_LOG_CA: Debug: DoG[%d] (%d x %d)\n", i,
		//		 DoG.Width(), DoG.Height());

		VT_HR_EXIT( VtCreateEdgeSegmentListUsingDoG( el[i], m_fchannels[i], edParams) );

		VT_DEBUG_LOG( "VT_LOG_CA: #Edge points detected[%d] = %d\n", i, (int)el[i].size() );

		FilterEdgeList(el[i]);
	}

	// keep edges that are near edges in other color channels
	VT_HR_EXIT( KeepAdjacentEdgels( el, w, h ) );

	// build the footprints for each edgel
	for( int i = 0; i < 3; i++ )
	{
		VT_HR_EXIT( m_edges[i].Create(w, h) );
		VT_HR_EXIT( ComputeEdgeFootprints(m_edges[i], el[i]) );
	}

	//DebugDumpEdge(m_edges[0], L"edge_filt",0);
	//DebugDumpEdge(m_edges[1], L"edge_filt",1);
	//DebugDumpEdge(m_edges[2], L"edge_filt",2);

Exit:
	return hr;
}

#define POLY_5

// Find radial function that maps dst radius to src radius
HRESULT 
CChromaticAberrationCorrectFast::ComputeRadialWarp( CVecd &poly, 
													CByteImg &src_edges, 
													CByteImg &dst_edges )
{
	HRESULT hr = NOERROR;
	/*
	a = poly[1], b = poly[2], c = poly[3]:

	a \sum r_dst^4 + b \sum r_dst^5 + c \sum r_dst^6 = \sum r_dst^2 (r_src - r_dst)
	a \sum r_dst^5 + b \sum r_dst^6 + c \sum r_dst^7 = \sum r_dst^3 (r_src - r_dst)
	a \sum r_dst^6 + b \sum r_dst^7 + c \sum r_dst^8 = \sum r_dst^4 (r_src - r_dst)
	or
	M (a b c)^T = m
	*/

	LucasKanadeTrack LKT;
	int r, c, w = src_edges.Width(), h = src_edges.Height();
	double slope[2], d[2];
	double length_norm = LengthNorm( m_width, m_height );
#ifdef POLY_5
	CVec4d ans;
	CMtx4x4d M(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	CVec4d m; m[0] = m[1] = m[2] = m[3] = 0;
#else // POLY_5
	CVec3d ans;
	CMtx3x3d M(0,0,0,0,0,0,0,0,0);
	CVec3d m; m[0] = m[1] = m[2] = 0;
#endif // POLY_5

#ifdef POLY_5
	VT_HR_EXIT( poly.Create(5) );
#else // POLY_5
	VT_HR_EXIT( poly.Create(4) );
#endif // POLY_5
	LKT.Setup( src_edges, dst_edges ); // ok, since Setup() only updates variables

	for (r=m_margin; r<h-m_margin; r+=5)
	{
		BYTE *psrc[3] = { src_edges.Ptr(r-1), src_edges.Ptr(r), src_edges.Ptr(r+1) };

		double r2 = r - m_c[1];
		for (c=m_margin; c<w-m_margin; c++)
		{
			// track only for maximal points
			if (psrc[1][c]>=226 &&
				((psrc[1][c]>=psrc[1][c-1] && psrc[1][c]>=psrc[1][c+1]) ||
				 (psrc[1][c]>=psrc[0][c  ] && psrc[1][c]>=psrc[2][c  ])))
			{
				double c2 = c - m_c[0];
				double mag = sqrt(c2*c2 + r2*r2);
				double d_radius = 0;
				double radius_src = mag/length_norm;
				double radius_dst = radius_src;

				if (mag>10.0f)
				{
					slope[0] = c2/mag; slope[1] = r2/mag;
					LKT.Track( d, c, r, (double)c, (double)r, slope );
					d_radius = (slope[0]*d[0] + slope[1]*d[1])/length_norm;
					radius_dst += d_radius;
				}

				double radius_2 = radius_dst*radius_dst;
				double radius_3 = radius_2*radius_dst;
				double radius_4 = radius_2*radius_2;
				double radius_5 = radius_2*radius_3;
				double radius_6 = radius_3*radius_3;
				double radius_7 = radius_3*radius_4;
				double radius_8 = radius_4*radius_4;

				M[0][0] += radius_4;
				M[0][1] += radius_5;
				M[0][2] += radius_6;
				M[1][2] += radius_7;
				M[2][2] += radius_8;
				m[0] -= radius_2*d_radius;
				m[1] -= radius_3*d_radius;
				m[2] -= radius_4*d_radius;
#ifdef POLY_5
				//double radius_9 = radius_5*radius_4;
				double radius_10 = radius_5*radius_5;

				M[0][3] += radius_7;
				M[2][3] += radius_8;
				M[3][3] += radius_10;
				m[3] -= radius_5*d_radius;
#endif // POLY_5
			}
		}
	}
	M[1][0] = M[0][1];
	M[1][1] = M[2][0] = M[0][2];
	M[2][1] = M[1][2];
#ifdef POLY_5
	M[3][0] = M[0][3];
	M[3][1] = M[1][3] = M[2][2];
	M[3][2] = M[2][3];
	ans = M.Inv() * m;
#else // POLY_5
	ans = M.Inv() * m;
#endif // POLY_5

	poly[0] = 1.0;
	poly[1] = ans[0];
	poly[2] = ans[1];
	poly[3] = ans[2];
#ifdef POLY_5
	poly[4] = ans[3];
	VT_DEBUG_LOG( "VT_LOG_CA: ans: (1, %f, %f, %f, %f)\n", ans[0], ans[1], ans[2], ans[3]);
#endif // POLY_5

Exit:
	return hr;
}

HRESULT CChromaticAberrationCorrectFast::AlignColorChannels()
{
	HRESULT hr = NOERROR;

	if (m_baligned) return hr;

	VT_DEBUG_LOG( "VT_LOG_CA: AlignColorChannels\n" );

	if ( !m_bmaskcomputed )
	{
		VT_HR_EXIT( ComputeEdgeImg( ) );
		m_bmaskcomputed = true;
	}

	for( int i = 1; i < 3; i++ )
	{
		CVecd mag_poly_norm; // polynomial function maps dst radius to src radius
		int ii = (m_reference+i)%3;
		VT_HR_EXIT( ComputeRadialWarp(mag_poly_norm, m_edges[ii], m_edges[m_reference]) );
		VT_HR_EXIT( RadialWarp(m_fchannels_aligned[ii], m_fchannels[ii], mag_poly_norm) );
	}
	VT_HR_EXIT( m_fchannels[m_reference].CopyTo(m_fchannels_aligned[m_reference]) );

	int n_refine = 0;
	if (n_refine>=1)
	{
		VT_DEBUG_LOG( "VT_LOG_CA: Now refining edges\n" );
		for (int i=0; i<n_refine; i++)
		{
			VT_HR_EXIT( RefineEdges() );
		}
	}
	VT_DEBUG_LOG( "VT_LOG_CA: AlignColorChannels done\n");

	m_baligned = true;

Exit:
	return hr;
}

class CRefineStatsIntegral
{
public:
	HRESULT Create(int iMaxW, int iHeight);
	void    StartRegion(const CFloatImg& imgSrc, const CFloatImg& imgDst,
						const CRect& region);
	void    GetSubRegionSums(int c0, int c1,
							 float& s, float& ssq, float& d, float& ds) const;

private:
	CFloatImg m_imgS;
	CFloatImg m_imgSsq;
	CFloatImg m_imgD;
	CFloatImg m_imgDS;
};

HRESULT
CRefineStatsIntegral::Create(int iMaxW, int iHeight)
{
	HRESULT hr;
	VT_HR_EXIT( m_imgS.Create(iMaxW, iHeight) );
	VT_HR_EXIT( m_imgSsq.Create(iMaxW, iHeight) );
	VT_HR_EXIT( m_imgD.Create(iMaxW, iHeight) );
	VT_HR_EXIT( m_imgDS.Create(iMaxW, iHeight) );
Exit:
	return hr;
}

void
CRefineStatsIntegral::StartRegion(const CFloatImg& imgSrc, const CFloatImg& imgDst,
								  const CRect& region)
{
	VT_ASSERT( region.Width() <= m_imgS.Width() &&
			   region.Height() == m_imgS.Height() );

	for( int y = 0; y < region.Height(); y++ )
	{
		const float* pS = imgSrc.Ptr(region.left, y+region.top);
		const float* pD = imgDst.Ptr(region.left, y+region.top);

		float* pSumS   = m_imgS.Ptr(y);
		float* pSumSsq = m_imgSsq.Ptr(y);
		float* pSumD   = m_imgD.Ptr(y);
		float* pSumDS  = m_imgDS.Ptr(y);

		float fSR = 0, fSsqR = 0, fDR = 0, fDSR  = 0;
		if( y == 0 )
		{
			for( int x = 0; x < region.Width(); x++ )
			{
				fSR   += pS[x];
				fSsqR += pS[x]*pS[x];
				fDR   += pD[x]; 
				fDSR  += pD[x]*pS[x];
				*pSumS++   = fSR;
				*pSumSsq++ = fSsqR;
				*pSumD++   = fDR;
				*pSumDS++  = fDSR;
			}
		}
		else
		{
			float* pSumS1   = m_imgS.Ptr(y-1);
			float* pSumSsq1 = m_imgSsq.Ptr(y-1);
			float* pSumD1   = m_imgD.Ptr(y-1);
			float* pSumDS1  = m_imgDS.Ptr(y-1);

			for( int x = 0; x < region.Width(); x++ )
			{
				fSR   += pS[x];
				fSsqR += pS[x]*pS[x];
				fDR   += pD[x];
				fDSR  += pD[x]*pS[x];
				*pSumS++   = fSR   + pSumS1[x];
				*pSumSsq++ = fSsqR + pSumSsq1[x];
				*pSumD++   = fDR   + pSumD1[x];
				*pSumDS++  = fDSR  + pSumDS1[x];
			}
		}
	}
}

void
CRefineStatsIntegral::GetSubRegionSums(int c0, int c1, 
									   float& s, float& ssq, float& d, float& ds) const
{
    VT_ASSERT( c0 < m_imgS.Width() && c1 < m_imgS.Width() );

	const float* pS   = m_imgS.Ptr(m_imgS.Height()-1);
	const float* pSsq = m_imgSsq.Ptr(m_imgSsq.Height()-1);
	const float* pD   = m_imgD.Ptr(m_imgD.Height()-1);
	const float* pDS  = m_imgDS.Ptr(m_imgDS.Height()-1);

	s   = pS[c1];
	ssq = pSsq[c1];
	d   = pD[c1];
	ds  = pDS[c1];
	if( c0 )
	{
		c0--;
		s   -= pS[c0];
		ssq -= pSsq[c0];
		d   -= pD[c0];
		ds  -= pDS[c0];
	}
}


// Rather than trying to account for spatially varying blur and other
// nonlinearities, at the edges, replicate the reference channel
// signal, but with the appropriate scale and shift
HRESULT 
CChromaticAberrationCorrectFast::RefineEdgesSingleChannel( const CFloatImg& imgSrc,
														   CFloatImg& imgDst ) const
{
	HRESULT hr = NOERROR;
	const CByteImg &edge_wt = m_edges[m_reference];
	const CFloatImg &ref = m_fchannels[m_reference];

	int r, c;
	int footprintdim = (2*m_spread + 1);
	float size = (float)(footprintdim*footprintdim);

	CRefineStatsIntegral si;
	int iMaxIntgrl = 16*footprintdim;
	VT_HR_EXIT( si.Create(iMaxIntgrl+m_margin+1, footprintdim) );

	// set up lookup table for weights associated with difference with reference value
	int min_intensity = 2, max_diff = 50;
	float fmin_intensity_norm = min_intensity/255.0f, fmax_diff_norm = max_diff/255.0f;
	float fdiff_LUT[256];
	for (r=0; r<=max_diff; r++)
	{
		fdiff_LUT[r] = sqrt(1.0f - (float)(r)/(float)(max_diff+1));
	}
	for (r=max_diff+1; r<256; r++)
	{
		fdiff_LUT[r] = 0;
	}

	for (r=0; r<m_height; r++)
	{
		const float *q = imgSrc.Ptr(r);
        float *o = imgDst.Ptr(r);

		// default is to just copy
		VtMemcpy(o, q, m_width*sizeof(float));

		if( r < m_margin || r >= m_height-m_margin )
		{
			continue;
		}

		const BYTE *w = edge_wt.Ptr(r);
		const float *t = ref.Ptr(r);
		for (c=m_margin; c<m_width-m_margin; )
		{
            // don't trust weak edges or very small ref values
			// or source is much larger than reference
			if (w[c]<50 || t[c]<fmin_intensity_norm || 
				(q[c-1]-t[c]>fmax_diff_norm && q[c]-t[c]>fmax_diff_norm && q[c+1]-t[c]>fmax_diff_norm))
			{
				c++;
				continue;
			}

			int cend = c + 1;
			for( ; cend<m_width-m_margin && cend < c+iMaxIntgrl; cend++ )
			{
				if (w[cend]<50 || t[cend]<fmin_intensity_norm || 
					(q[cend-1]-t[cend]>fmax_diff_norm && q[cend]-t[cend]>fmax_diff_norm && q[cend+1]-t[cend]>fmax_diff_norm))
				{
					break;
				}
			}

			// compute the integral images for a span
			si.StartRegion(ref, imgSrc, CRect(c-m_spread, r-m_spread, 
											 cend-1+m_spread+1, r+m_spread+1) );

            int f = 0;
			for( int cc = c; cc < cend; cc++, f++ )
			{			
				int idiff = (int)(255.0f*(q[cc] - t[cc]));
				float fdiff_wt = (idiff>=0) ? fdiff_LUT[idiff] : 1;
				CMtx2x2f M;
				CVec2f m;
				si.GetSubRegionSums(f, f+footprintdim-1, 
									M[0][0], M[1][0], m[0], m[1]);
				M[0][1] = size;
				M[1][1] = M[0][0];

				CVec2f ans = M.Inv() * m;
				float wt = fdiff_wt*w[cc]/255.f;
				o[cc] = q[cc] + wt*(ans[0]*t[cc] + ans[1] - q[cc]);
			}

			// move to the end of this region
			c = cend;
		}
	}

Exit:
	return hr;
}

HRESULT CChromaticAberrationCorrectFast::RefineEdges()
{
	HRESULT hr = NOERROR;
	int i;

	for (i=1; i<3; i++)
	{
		int index = (m_reference+i)%3;

		VT_HR_EXIT( RefineEdgesSingleChannel( m_fchannels_aligned[index],
											  m_fchannels[index] ) );
		VT_HR_EXIT( m_fchannels[index].CopyTo( m_fchannels_aligned[index] ) );
	}

Exit:
	return hr;
}
