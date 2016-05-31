#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class FlattenLayer : public Layer
	{
	public:
		FlattenLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~FlattenLayer();
		virtual void Forward();
		virtual void InitActs();
	};
}


