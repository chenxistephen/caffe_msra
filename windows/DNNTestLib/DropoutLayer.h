#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class DropoutLayer : public Layer
	{
	public:
		DropoutLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~DropoutLayer();
		virtual void Forward();
		virtual void InitActs();

	private:
		float m_dropout_ratio;
	};
}