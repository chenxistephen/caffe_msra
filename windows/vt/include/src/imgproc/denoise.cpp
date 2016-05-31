
#include "stdafx.h"

#include "vt_denoise.h"

using namespace vt;

#define N_ITER_DEFAULT					2
#define ALPHA_DEFAULT					1.0

// constants
const float CBayesDenoiseAD::kDfltEdgeSigma = 5.0f;

// anisotropic diffusion
//

CBayesDenoiseAD::CBayesDenoiseAD()
{
	m_niter = N_ITER_DEFAULT;
	m_alpha = (float)ALPHA_DEFAULT;
	m_edge_sigma = 0.0f;
}

CBayesDenoiseAD::~CBayesDenoiseAD()
{
}

HRESULT CBayesDenoiseAD::Set_edge_sigma( float edge_sigma )
{
	HRESULT hr = NOERROR;
	if (edge_sigma<=0) return hr;

	if( m_edge_sigma != edge_sigma )
	{
		// most entries in this look-up table will end up being 0, so
		// compute where the zeros start and only do exp() up to that point
		UINT uPrevZeroIndex;
		if( m_edge_mult_array.size() == 0 )
		{
			uPrevZeroIndex = ElTraits<UInt16>::MaxVal();
			VT_HR_EXIT( m_edge_mult_array.resize( ElTraits<UInt16>::MaxVal()+1 ) );
		}
		else
		{
			// FLT_MIN = exp(-(i/edge)*(i/edge)) =>
			// i = edge*sqrt(-ln(FLT_MIN))
			//     sqrt(-ln(FLT_MIN)) = 9.3454
			uPrevZeroIndex = (int)ceil(m_edge_sigma*9.3454f);
			uPrevZeroIndex = VtMin(uPrevZeroIndex, (UINT)ElTraits<UInt16>::MaxVal());
		}

		m_edge_sigma = edge_sigma;
		UINT uZeroIndex = (int)ceil(m_edge_sigma*9.3454f);
		uZeroIndex = VtMin(uZeroIndex, (UINT)ElTraits<UInt16>::MaxVal());
		if( uZeroIndex < uPrevZeroIndex )
		{
			ZeroMemory(m_edge_mult_array.begin()+uZeroIndex+1, 
					   (uPrevZeroIndex-uZeroIndex)*sizeof(float) );
		}

		// incrementally compute exp(-(i/m_edge_sigma)^2);
        float res = 1.f;
        float v   = 1.f / m_edge_sigma;
        float b = exp(-v*v);
		float c = b*b;
		for (UINT i=0; i<=uZeroIndex; i++)
		{
			m_edge_mult_array[i] = res;
            res *= b;
            b   *= c;
		}
	}

Exit:
	return hr;
}

void CBayesDenoiseAD::Set_niter( int niter )
{
	if (niter<=0) return;
	m_niter = niter;
}

void CBayesDenoiseAD::Set_alpha( float alpha )
{
	if (alpha<0.0f || alpha>1.0f) return;
	m_alpha = alpha;
}

HRESULT CBayesDenoiseAD::ComputeBayesDenoised(const CShortImg& src, 
											  CShortImg& dst)
{
	HRESULT hr = NOERROR;

	// note the src and dst can be the same - because the loop below has a 
	// rolling buffer that caches 3 lines ahead

	if( src.Width()  != dst.Width()   || 
		src.Height() != dst.Height()  ||
		src.Bands()  != 1 || dst.Bands() != 1 )
	{
		return E_INVALIDARG;
	}

	if( m_edge_mult_array.size() == 0 )
	{
		VT_HR_RET( Set_edge_sigma( kDfltEdgeSigma ) );
	}

	int iter, r, c, w = src.Width(), h = src.Height();

    CFloatImg buf, dx, dy, res;
	VT_HR_EXIT( buf.Create(w,3) );
	VT_HR_EXIT( res.Create(w,1) );
	VT_HR_EXIT( dx.Create(w-1,3) );
	VT_HR_EXIT( dy.Create(w-1,3) );

	// copy first and last rows
	if( src.Ptr() != dst.Ptr() )
	{
		VtConvertSpan(dst.Ptr(0),   src.Ptr(0), w); 
		VtConvertSpan(dst.Ptr(h-1), src.Ptr(h-1), w);
	}

    const CShortImg* pSrc = &src;
	for (iter=0; iter<m_niter; iter++)
	{
		VT_DEBUG_LOG( "IP_LOG_DENOISE: CBayesDenoiseAD: Iter = %d\n", iter );

		int sidx = 0;
		int rr = 0;
		for (r=1; r<h-1; r++)
		{   
			for( ; rr < r+2 && rr < h-1; rr++, sidx=(sidx+1)%3)
			{
				// compute the local gradients
				float *gx = dx.Ptr(sidx);
				float *gy = dy.Ptr(sidx);
				const UINT16 *p[2] = { pSrc->Ptr(rr), pSrc->Ptr(rr+1) };
				for (c=0; c<w-1; c++)
				{
					gx[c] = m_edge_mult_array[abs(p[0][c] - p[0][c+1])];
					gy[c] = m_edge_mult_array[abs(p[0][c] - p[1][c])];
				}

				// convert source line(s) to float
				VtConvertSpan(buf.Ptr(sidx), pSrc->Ptr(rr), w);
			}

			// source pointers
			int si1 = (sidx+1)%3;
			int si2 = (sidx+2)%3;
			const float *gx  = dx.Ptr(si1)+1;
        	const float *gy  = dy.Ptr(si1)+1;
			const float *gyp = dy.Ptr(sidx)+1;
			const float *pf[3] = { buf.Ptr(sidx)+1, buf.Ptr(si1)+1, buf.Ptr(si2)+1 };

			// dest pointer
			float *q = res.Ptr() + 1;   

			// just do a copy of pixel 0
            res.Ptr()[0] = buf.Ptr(si1)[0];

			// run the denoise filter
			c = 1;
#if (defined(_M_IX86) || defined(_M_AMD64))
			if( g_SupportSSE2() )
			{
				__m128 x7 = _mm_set1_ps(m_alpha);

				while( (c+4) < (w-1) )
				{
					// read previous row, mul by prev grads
					__m128 x0 = _mm_loadu_ps(pf[0]);
					__m128 x1 = _mm_loadu_ps(gyp);
					x0 = _mm_mul_ps(x0, x1);

					// read current row, mul by x-1, m_offset, x+1 grads
					__m128 x2 = _mm_loadu_ps(pf[1]-1);
					__m128 x3 = _mm_loadu_ps(gx-1);
					x2 = _mm_mul_ps(x2, x3);
					x0 = _mm_add_ps(x0, x2);
					x1 = _mm_add_ps(x1, x3);

					x2 = _mm_loadu_ps(pf[1]);
					x2 = _mm_mul_ps(x2, x7);
					x0 = _mm_add_ps(x0, x2);
					x1 = _mm_add_ps(x1, x7);

					x2 = _mm_loadu_ps(pf[1]+1);
					x3 = _mm_loadu_ps(gx);
					x2 = _mm_mul_ps(x2, x3);
					x0 = _mm_add_ps(x0, x2);
					x1 = _mm_add_ps(x1, x3);

					// read next row, mul by next grads
					x2 = _mm_loadu_ps(pf[2]);
					x3 = _mm_loadu_ps(gy);
					x2 = _mm_mul_ps(x2, x3);
					x0 = _mm_add_ps(x0, x2);
					x1 = _mm_add_ps(x1, x3);

					// normalize the result and store
					x0 = _mm_div_ps(x0, x1);
					_mm_storeu_ps(q, x0);

                    // increment stuff
					c+=4;
                    gx+=4; gy+=4; gyp+=4; 
                    pf[0]+=4; pf[1]+=4; pf[2]+=4;
                    q+=4;
				}
			}
#endif

			for (; c<w-1; c++, gx++, gy++, gyp++, pf[0]++, pf[1]++, pf[2]++, q++ )
			{
				float factor = (gyp[0]+gx[-1]+m_alpha+gx[0]+gy[0]);
				*q = ( gyp[0]*pf[0][0]   + 
					   gx[-1]*pf[1][-1]  + 
					   m_alpha*pf[1][0]  + 
                       gx[0]*pf[1][1]    + 
					   gy[0]*pf[2][0] )  / factor;
			}

			// just do a copy of the last pixel
            res.Ptr()[w-1] = buf.Ptr(si1)[w-1];

			// it is safe to copy back to dst because we 
			// are buffering the previous values a line ahead
			VtConvertSpan(dst.Ptr(r), res.Ptr(), w);
		}

		// for subsequent iterations dst is the source
		pSrc = &dst; 
	}

Exit:
	return hr;
}

