
#include "stdafx.h"

#include "vt_lucasKanade.h"

using namespace vt;

//-----------------------------------------------------------------------------

LucasKanadeTrack::LucasKanadeTrack(void)
{
	m_pimg0 = m_pimg1 = 0;
	m_width = m_height = 0;
	m_patch_size = 0;
	m_niters = 0;
}

LucasKanadeTrack::~LucasKanadeTrack(void)
{
}

#define LARGE_ERR		1.0e8

static double pixel_gradient( double x, double y, int iu, int iv, 
	        				  CByteImg &img, double *grad )
{
	int width = img.Width(), height = img.Height();
	double c0, c1, c2, c3;
#ifdef PIX_GRAD_3x3
	// handle boundaries
	if (iu<1) iu = 1; else if (iu>=width-1) iu = width - 2;
	if (iv<1) iv = 1; else if (iv>=height-1) iv = height - 2;
	c0 = (float)img.Pix(iu-1,iv-1);
	c1 = (float)img.Pix(iu+1,iv-1);
	c2 = (float)img.Pix(iu-1,iv+1);
	c3 = (float)img.Pix(iu+1,iv+1);
#else // PIX_GRAD_3x3
	// handle boundaries
	if (iu<0) iu = 0; else if (iu>=width-1) iu = width - 2;
	if (iv<0) iv = 0; else if (iv>=height-1) iv = height - 2;
	c0 = (double)img.Pix(iu,iv);
	c1 = (double)img.Pix(iu+1,iv);
	c2 = (double)img.Pix(iu,iv+1); 
	c3 = (double)img.Pix(iu+1,iv+1); 
#endif // PIX_GRAD_3x3
	double d0 = c1 - c0,
		   d1 = c2 - c0,
		   d2 = c3 - c2 - d0;

	c0 += d0 * x;
	d1 += d2 * x;
	d0 += d2 * y;
	c0 += d1 * y;

	grad[0] = 0.5*d0;
	grad[1] = 0.5*d1;

	return c0;
}

double LucasKanadeTrack::LK_step(CByteImg &img0, int ix0, int iy0,
                                 CByteImg &img1, double x1, double y1,
                                 double d[2], double *min_max_eigen)
{
	int width = img0.Width(), height = img0.Height();
	int ix1 = (int) x1, iy1 = (int) y1;
    double fx1 = x1 - ix1, fy1 = y1 - iy1;
    double Z[2][2], e[2], grad[2];
    double error = 0.0f, det, invdet, discr;
    int k, l, s = m_patch_size / 2;

    /* Clear out the Hessian and gradients */
    Z[0][0] = Z[0][1] = Z[1][0] = Z[1][1] = 0.0f;
    e[0] = e[1] = 0.0f;

    /* Accumulate the Hessian and gradients */
    for (k = -s; k < m_patch_size-s; k++) {
        int iy2 = iy1 + k;

        for (l = -s; l < m_patch_size-s; l++) {
            int ix2 = ix1 + l;
			int uu = ix0+l, vv = iy0+k;
			if (uu<0 || uu>=width || vv<0 || vv>=height) continue;
            double v0 = img0.Pix(uu,vv);
            double v1;
            double err;

            if (ix2>=0 && iy2>=0 &&
                (ix2+1)<width && (iy2+1)<height) {
                v1 = pixel_gradient(fx1, fy1, ix2, iy2,
                                    img1, grad);
                err = v0 - v1;

                error += err*err;
                e[0] += err * grad[0];
                e[1] += err * grad[1];
                Z[0][0] += grad[0] * grad[0];
                Z[0][1] += grad[0] * grad[1]; /* == Z[1][0] by symmetry */
                Z[1][1] += grad[1] * grad[1];
            }
            else {
                return( LARGE_ERR );
            }
        }
    }

    /* Solve for the correction d */
    det = Z[0][0]*Z[1][1] - Z[0][1]*Z[0][1];
    invdet = (det > 0.0) ? (1.0 / det) : 0.0;
    d[0] = 0.5 * (Z[1][1] * e[0] - Z[0][1] * e[1]) * invdet;
    d[1] = 0.5 * (Z[0][0] * e[1] - Z[0][1] * e[0]) * invdet;

	if (min_max_eigen) {
		discr = sqrt((Z[0][0] - Z[1][1])*(Z[0][0] - Z[1][1]) +
						4.0 * Z[0][1] * Z[0][1]);
		min_max_eigen[0] = 0.5 * (Z[0][0] + Z[1][1] - discr); /* min eigenvalue */
		min_max_eigen[1] = 0.5 * (Z[0][0] + Z[1][1] + discr); /* max eigenvalue */
		min_max_eigen[2] = (Z[0][0] > Z[1][1]) ?  /* eigenvector direction */
				(180.0 / VT_PI * atan2(Z[0][0] - min_max_eigen[0], Z[0][1])) :
				(180.0 / VT_PI * atan2(Z[0][1], Z[1][1] - min_max_eigen[0]));
	}

    return error;
}

void LucasKanadeTrack::Setup( CByteImg &img0, CByteImg &img1,
							  int patch_size, int niters )
{
	m_pimg0 = &img0;
	m_pimg1 = &img1;
	m_width = img0.Width(); 
	m_height = img0.Height(); 
	m_patch_size = patch_size;
	m_niters = niters;

#ifdef USE_HALF_RES
	m_img0_half.Create(m_width/2,m_height/2);
	m_img1_half.Create(m_width/2,m_height/2);
	VtZoom( m_img0_half, *m_pimg0, 2.0f, 0.0f, 2.0f, 0.0f);
	VtZoom( m_img1_half, *m_pimg1, 2.0f, 0.0f, 2.0f, 0.0f);
#endif // USE_HALF_RES
}

inline void EnforceRadialGradient( double d[2], double slope[2] )
{
	double s = d[0]*slope[0] + d[1]*slope[1];

	d[0] = s*slope[0]; d[1] = s*slope[1];
}

double LucasKanadeTrack::Track( double d_best[2], int ix0, int iy0, double x1, double y1,
							    double *slope )
{
#ifdef USE_HALF_RES
	ix0 /= 2; iy0 /= 2; x1 /= 2.0f; y1 /= 2.0f;
	CByteImg &img0 = m_img0_half, &img1 = m_img1_half;
#else // USE_HALF_RES
	CByteImg &img0 = *m_pimg0, &img1 = *m_pimg1;
#endif // USE_HALF_RES

	double d[2], RMS_best = LK_step( img0, ix0, iy0, img1, x1, y1, d );
	if (slope)
		EnforceRadialGradient( d, slope );
	d_best[0] = d[0]; d_best[1] = d[1];

	for (int i=1; i<m_niters; i++)
	{
		x1 += d[0]; y1 += d[1];
		double RMS = LK_step( img0, ix0, iy0, img1, x1, y1, d );
		if (slope)
			EnforceRadialGradient( d, slope );

		if (RMS_best>RMS)
		{
			RMS_best = RMS;
			d_best[0] += d[0]; d_best[1] += d[1];
		}
		else
			break;
	}

#ifdef USE_HALF_RES
	d_best[0] *= 2.0f;
	d_best[1] *= 2.0f;
#endif // USE_HALF_RES

	return RMS_best;
}
