#pragma once

#include "vtcore.h"
#include "vt_channels.h"
#include "vt_dpyramid.h"

namespace vt
{
	struct CChannelsPyramidParams
	{
		CChannelsPyramidParams()
			: approxSubOctaveGradients(true),
			approxSubOctaveGradientLambda(0.2f),
			numPostSmoothIterations(1)
		{
		}

		bool approxSubOctaveGradients;
		float approxSubOctaveGradientLambda;
		unsigned int numPostSmoothIterations;
	};
		
	template <typename StaticParams>
	class CChannelsPyramidEngine
	{
	public:
		CChannelsPyramidEngine();

	private:
		CChannelsPyramidEngine(const CChannelsPyramidEngine&);
		CChannelsPyramidEngine& operator=(const CChannelsPyramidEngine&);

	public:
		HRESULT Initialize(unsigned int inputWidth, unsigned int inputHeight, 
			const CChannelsEngine<StaticParams>& engine, 
			const CChannelsPyramidParams& pyramidParams);
		
		bool IsInitialized() const;
				
		HRESULT ComputeChannelsPyramid(CDetectionPyramid& channelsPyramid, 
			const CByteImg& input) const;

	private:
		unsigned int GetNextLargerInputLevel(unsigned int inputWidth, 
			unsigned int inputHeight, bool& isPerfectMatch) const;

		HRESULT ApplyPostSmoothing(CDetectionPyramid& channelsPyramid) const;

	private:
		unsigned int baseInputWidth;
		unsigned int baseInputHeight;

		CChannelsPyramidParams params;
	
		bool isInitialized;
		
		const CChannelsEngine<StaticParams>* channelsEngine;
		mutable CBytePyramid inputPyramid;		
		mutable CByteImg inputBuffer;
		mutable CLumaByteImg oneBandInputBuffer1;
		mutable CLumaByteImg oneBandInputBuffer2;
	};
}
