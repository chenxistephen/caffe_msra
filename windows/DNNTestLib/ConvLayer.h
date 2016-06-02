#pragma once

#include "Layer.h"
#include "Blas.h"

namespace DNNTestLib
{

	class ConvLayer : public Layer
	{
	public:
		ConvLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~ConvLayer();
		virtual void Dump(std::ostream &stream);
		virtual void Forward();
		virtual int GetInputPadX();
		virtual int GetInputPadY();
		virtual void InitActs();
		// Added by Stephen Chen
		virtual void CopyWeight(float *buffer, int buffer_size) const;
		virtual void CopyBias(float *buffer, int buffer_size) const;
		virtual CPUData* GetWeight() const;
		virtual CPUData* GetBias() const;
	public://private:
		// properties
		int m_padding_x, m_padding_y, m_stride_x, m_stride_y;
		CPUData *m_weight, *m_biases;
#if USE_BLAS
		static CPUData m_im2col;
		static CPUData m_output_buffer;
#endif
	};

}

