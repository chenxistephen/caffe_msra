#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class LocalShareLayer : public Layer
	{
	public:
		LocalShareLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~LocalShareLayer();
		virtual void Dump(std::ostream &stream);
		virtual void Forward();
		virtual int GetInputPadX();
		virtual int GetInputPadY();
		virtual void InitActs();

	private:
		// properties
		int m_padding_x, m_padding_y, m_stride_x, m_stride_y;
		int m_share_x, m_share_y, m_filter_num_x, m_filter_num_y;
		int m_output_channels, m_input_channels, m_filter_width, m_filter_height;
		std::vector<CPUData*> m_weight, m_biases;
	};
}

