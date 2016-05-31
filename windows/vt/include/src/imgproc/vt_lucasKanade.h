#pragma once

#include "vtcommon.h"

namespace vt {

#define PATCH_SIZE_DEF					7
#define N_ITERS_DEF						5
// #define USE_HALF_RES

class LucasKanadeTrack
{
public:
    LucasKanadeTrack();
    ~LucasKanadeTrack();

	void Setup( vt::CByteImg &img0, vt::CByteImg &img1, 
		int patch_size=PATCH_SIZE_DEF, int niters=N_ITERS_DEF );
	// given (ix0,iy0) in img0, (x1,y1) in img0
	// compute shift d
	// returns the RMS error
	double Track( double d[2], int ix0, int iy0, double x1, double y1,
		double *slope=0 );

private:
	double LK_step(vt::CByteImg &img0, int ix0, int iy0,
                   vt::CByteImg &img1, double x1, double y1,
                   double d[2], double *min_max_eigen=0);

	vt::CByteImg *m_pimg0, *m_pimg1;
#ifdef USE_HALF_RES
	vt::CByteImg m_img0_half, m_img1_half;
#endif // USE_HALF_RES
	int m_width, m_height;
	int m_patch_size;
	int m_niters;
};

};
