#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class AcrossChannelNormLayer : public Layer
	{
	public:
		AcrossChannelNormLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~AcrossChannelNormLayer();
		virtual void Forward();
		virtual void InitActs();

	private:
		float m_fAlpha;
		float m_fBeta;
		int m_iWindowSize;
		CPUData m_dataScale;
	};
}
