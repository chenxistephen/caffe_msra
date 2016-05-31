#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class FCLayer : public Layer
	{
	public:
		FCLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~FCLayer();
		virtual void Dump(std::ostream &stream);
		virtual void Forward();
		virtual void InitActs();

	protected:
		// properties
		CPUData *m_weight, *m_biases;
	};

}


