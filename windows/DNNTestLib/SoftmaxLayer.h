#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class SoftmaxLayer : public Layer
	{
	public:
		SoftmaxLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~SoftmaxLayer();
		virtual void InitActs();
		virtual void Forward();
		virtual void Dump(std::ostream &stream);
	};
}



