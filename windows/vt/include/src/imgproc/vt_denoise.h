#pragma once

#include "vtcommon.h"

namespace vt {

//
// denoise by anisotropic diffusion
//

class CBayesDenoiseAD
{
public:
    static const float kDfltEdgeSigma;

public:
    CBayesDenoiseAD(void);
    ~CBayesDenoiseAD(void);

	HRESULT ComputeBayesDenoised(const CShortImg& src, CShortImg& dst);
	HRESULT Set_edge_sigma( float edge_sigma );
	float   Get_edge_sigma() { return m_edge_sigma; }
	void    Set_alpha( float alpha );
	float   Get_alpha() { return m_alpha; }
	void    Set_niter( int niter );
	int     Get_niter() { return m_niter; }

private:
	int   m_niter;
	float m_alpha;
	float m_edge_sigma;
	vector<float> m_edge_mult_array;
};

};
