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
		// Added by Stephen Chen
		virtual void CopyWeight(float *buffer, int buffer_size) const;
		virtual void CopyBias(float *buffer, int buffer_size) const;
		virtual CPUData* GetWeight() const;
		virtual CPUData* GetBias() const;
	protected:
		// properties
		CPUData *m_weight, *m_biases;
	};

}


