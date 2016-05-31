#pragma once

#include "Layer.h"
#include "Blas.h"

namespace DNNTestLib
{

	class SkipConvLayer : public Layer
	{
	public:
		SkipConvLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~SkipConvLayer();
		virtual void Dump(std::ostream &stream);
		virtual void Forward();
		virtual int GetInputPadX();
		virtual int GetInputPadY();
		virtual void InitActs();

	private:
		// properties
		int m_padding_x, m_padding_y, m_stride_x, m_stride_y, m_skip_x, m_skip_y;
		CPUData *m_weight, *m_biases;
#if USE_BLAS
		static CPUData m_im2col;
		static CPUData m_output_buffer;
#endif
	};

}
