#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class ConcatLayer : public Layer
	{
	public:
		ConcatLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~ConcatLayer();
		virtual void Forward();
		virtual void InitActs();
	};
}