#include "FlattenLayer.h"
#include <assert.h>

using namespace DNNTestLib;

FlattenLayer::FlattenLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
}

FlattenLayer::~FlattenLayer()
{
}

void FlattenLayer::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);
	assert(m_min_pad_x == 0 && m_min_pad_y == 0);

	CPUData *input = m_inputs[0]->GetActs();
	int num = input->GetNum();
	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(num, input->GetDataSize() / num, 1, 1, 0, 0);
	else
		m_layerActs->Resize(num, input->GetDataSize() / num, 1, 1, 0, 0);
}

void FlattenLayer::Forward()
{
	const CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;
	const int input_channels = input->GetChannels();
	const int input_width = input->GetWidth();
	const int input_height = input->GetHeight();
	const int input_stride = input->GetStrideH();
	const int output_stride = input_width * input_height;

	for (int n = 0; n < input->GetNum(); ++n)
	{
		const float *src = input->GetDataPointer(n, 0, 0, 0);
		float *dst = output->GetDataPointer(n, 0, 0, 0);
		for (int y = 0; y < input_height; y++)
		{
			for (int x = 0; x < input_width; x++)
			{
				for (int i = 0; i < input_channels; i++)
				{
					dst[(i*input_height + y)*input_width + x] = src[i + x*input_channels];
				}
			}
			src += input_stride;
		}
	}

}