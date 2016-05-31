#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	template <class NeuronType>
	class NeuronLayer : public Layer
	{
	public:
		NeuronLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual void Forward();
		virtual int GetInputPadX();
		virtual int GetInputPadY();
		virtual void InitActs();
		virtual ~NeuronLayer();
	};
}

