#pragma once

#include "vtcommon.h"

namespace vt {

#define REFERENCE_DEF					1 // green?

class CChromaticAberrationCorrectFast
{
public:
    CChromaticAberrationCorrectFast(void);
    ~CChromaticAberrationCorrectFast(void);

	HRESULT SetImage( CFloatImg *fimgs, int reference=REFERENCE_DEF );
	HRESULT ComputeEdgeImg( );
	HRESULT AlignColorChannels();
	HRESULT RefineEdgesSingleChannel(const CFloatImg& imgSrc, 
									 CFloatImg& imgDst ) const;
	HRESULT RefineEdges();
	bool IsAligned() { return m_baligned; }
	CFloatImg *GetChannelsAligned() { return m_fchannels_aligned; }

protected:
	int m_width, m_height;
	double m_c[2]; // image center
	int m_reference; // reference channel
	int m_margin; // area for evaluation to avoid image boundary effects
	int m_spread; // neighborhood size for edge (for searching parameters)
	bool m_bmaskcomputed;
	bool m_baligned;

	CFloatImg *m_fchannels; // input images
	CByteImg  m_edges[3];   // edge strength image (separate channels)
	CFloatImg m_fchannels_aligned[3]; // output

	HRESULT ComputeRadialWarp( CVecd &poly, CByteImg &src_edges, CByteImg &dst_edges );
};

};
