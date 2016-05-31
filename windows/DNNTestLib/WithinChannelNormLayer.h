#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class WithinChannelNormLayer : public Layer
	{
	public:
		WithinChannelNormLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~WithinChannelNormLayer();
		virtual void Forward();
		virtual void InitActs();
		virtual int GetInputPadX();
		virtual int GetInputPadY();

	private:
		float m_alpha;
		float m_beta;
		int m_local_x;
		int m_local_y;
	};
}