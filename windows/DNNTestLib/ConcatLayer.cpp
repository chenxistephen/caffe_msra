#include "ConcatLayer.h"
#include <assert.h>

using namespace DNNTestLib;

ConcatLayer::ConcatLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
}

ConcatLayer::~ConcatLayer()
{
}

void ConcatLayer::InitActs()
{
	assert(GetInputNumber() > 0);

	int width = m_inputs[0]->GetActs()->GetWidth();
	int height = m_inputs[0]->GetActs()->GetHeight();
	int num = m_inputs[0]->GetActs()->GetNum();
	int channels = 0;
	for (size_t i = 0; i < m_inputs.size(); i++)
	{
		CPUData *input = m_inputs[i]->GetActs();
		assert(input->GetDataStorageType() == CPUData::CWHN);
		assert(input->GetWidth() == width);
		assert(input->GetHeight() == height);
		assert(input->GetNum() == num);
		channels += input->GetChannels();
	}

	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(num, channels, width, height, m_min_pad_x, m_min_pad_y);
	else
		m_layerActs->Resize(num, channels, width, height, m_min_pad_x, m_min_pad_y);
}

void ConcatLayer::Forward()
{
	for (int n = 0; n < m_inputs[0]->GetActs()->GetNum(); ++n)
	{
		int start_channel = 0;
		for (size_t i = 0; i < m_inputs.size(); i++)
		{
			const CPUData *input = m_inputs[i]->GetActs();
			CPUData *output = m_layerActs;

			const float *src = input->GetDataPointer(n, 0, 0, 0);
			const int input_channels = input->GetChannels();
			const int input_width = input->GetWidth();
			const int input_height = input->GetHeight();
			const int input_det_stride = input->GetStrideH() - input_width * input_channels;
			float *dst = output->GetDataPointer(n, start_channel, 0, 0);
			const int output_channels = output->GetChannels();
			const int output_det_stride = output->GetStrideH() - input_width * output_channels;
			for (int y = 0; y < input_height; y++)
			{
				for (int x = 0; x < input_width; x++)
				{
					for (int i = 0; i < input_channels; i++)
					{
						dst[i] = src[i];
					}
					src += input_channels;
					dst += output_channels;
				}
				src += input_det_stride;
				dst += output_det_stride;
			}
			start_channel += input_channels;
		}
	}
}

