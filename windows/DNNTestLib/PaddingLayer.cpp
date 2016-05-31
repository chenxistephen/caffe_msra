#include "PaddingLayer.h"
#include "NetUtils.h"

using namespace DNNTestLib;

PaddingLayer::PaddingLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_padding_x = ReadData<int>(stream);
	m_padding_y = ReadData<int>(stream);
}

PaddingLayer::~PaddingLayer()
{
}

int PaddingLayer::GetInputPadX()
{
	return m_padding_x;
}

int PaddingLayer::GetInputPadY()
{
	return m_padding_y;
}

void PaddingLayer::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);
	CPUData *input = m_inputs[0]->GetActs();
	int num = input->GetNum();
	if (m_layerActs == NULL)
	{
		m_layerActs = CPUData::Create(num, input->GetChannels(),
			input->GetWidth() + m_padding_x * 2,
			input->GetHeight() + m_padding_y * 2,
			m_min_pad_x, m_min_pad_y);
	}
	else
	{
		m_layerActs->Resize(num, input->GetChannels(),
			input->GetWidth() + m_padding_x * 2,
			input->GetHeight() + m_padding_y * 2,
			m_min_pad_x, m_min_pad_y);
	}
}


void PaddingLayer::Forward()
{
	CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;

	const int num = input->GetNum();
	const int width_channels = input->GetChannels() * input->GetWidth();
	const int height = input->GetHeight();
	const int input_stride = input->GetStrideH();
	const int output_stride = output->GetStrideH();
	for (int n = 0; n < num; n++)
	{
		const float *data = input->GetDataPointer(n, 0, 0, 0);
		float *target = output->GetDataPointer(n, 0, m_padding_x, m_padding_y);
		for (int h = 0; h < height; h++)
		{
			memcpy(target, data, width_channels*sizeof(float));
			data += input_stride;
			target += output_stride;
		}
	}
	
}