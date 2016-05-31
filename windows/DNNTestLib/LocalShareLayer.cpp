#include "LocalShareLayer.h"
#include "NetUtils.h"
#include "LayerActs.h"

using namespace DNNTestLib;

LocalShareLayer::LocalShareLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_padding_x = ReadData<int>(stream);
	m_padding_y = ReadData<int>(stream);
	m_stride_x = ReadData<int>(stream);
	m_stride_y = ReadData<int>(stream);
	m_share_x = ReadData<int>(stream);
	m_share_y = ReadData<int>(stream);
	m_filter_num_x = ReadData<int>(stream);
	m_filter_num_y = ReadData<int>(stream);
	m_output_channels = ReadData<int>(stream);
	m_input_channels = ReadData<int>(stream);
	m_filter_width = ReadData<int>(stream);
	m_filter_height = ReadData<int>(stream);
	int filter_num = m_filter_num_x * m_filter_num_y;
	m_weight.resize(filter_num);
	m_biases.resize(filter_num);
	for (int i = 0; i < filter_num; i++)
		m_weight[i] = CPUData::CreateFromDump(stream);
	for (int i = 0; i < filter_num; i++)
		m_biases[i] = CPUData::CreateFromDump(stream);

	// check parameters
	for (int i = 0; i < filter_num; i++)
	{
#if USE_SSE
		assert(m_weight[i]->GetNum() % 4 == 0);
#endif
		assert(m_weight[i]->GetNum() == m_output_channels);
		assert(m_weight[i]->GetChannels() == m_input_channels);
		assert(m_weight[i]->GetWidth() == m_filter_width);
		assert(m_weight[i]->GetHeight() == m_filter_height);
		assert(m_biases[i]->GetNum() == m_weight[i]->GetNum());
		assert(m_biases[i]->GetChannels() == 1);
		assert(m_biases[i]->GetWidth() == 1);
		assert(m_biases[i]->GetHeight() == 1);
	}
}

void LocalShareLayer::Dump(std::ostream &stream)
{
	Layer::Dump(stream);
	WriteData(stream, m_padding_x);
	WriteData(stream, m_padding_y);
	WriteData(stream, m_stride_x);
	WriteData(stream, m_stride_y);
	WriteData(stream, m_share_x);
	WriteData(stream, m_share_y);
	WriteData(stream, m_filter_num_x);
	WriteData(stream, m_filter_num_y);
	WriteData(stream, m_output_channels);
	WriteData(stream, m_input_channels);
	WriteData(stream, m_filter_width);
	WriteData(stream, m_filter_height);
	int filter_num = m_filter_num_x * m_filter_num_y;
	for (int i = 0; i < filter_num; i++)
		m_weight[i]->Dump(stream);
	for (int i = 0; i < filter_num; i++)
		m_biases[i]->Dump(stream);
}

LocalShareLayer::~LocalShareLayer()
{
	int filter_num = m_filter_num_x * m_filter_num_y;
	for (int i = 0; i < filter_num; i++)
		SaveRelease(m_weight[i]);
	for (int i = 0; i < filter_num; i++)
		SaveRelease(m_biases[i]);
}

int LocalShareLayer::GetInputPadX()
{
	return m_padding_x;
}

int LocalShareLayer::GetInputPadY()
{
	return m_padding_y;
}

void LocalShareLayer::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);

	const CPUData *input = m_inputs[0]->GetActs();
	int num = input->GetNum();
	int filter_num = m_filter_num_x * m_filter_num_y;
	assert(input->GetChannels() == m_input_channels);
	int output_width = (input->GetWidth() + 2 * m_padding_x - m_filter_width) / m_stride_x + 1;
	int output_height = (input->GetHeight() + 2 * m_padding_y - m_filter_height) / m_stride_y + 1;
	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(num, m_output_channels, output_width, output_height, m_min_pad_x, m_min_pad_y);
	else
		m_layerActs->Resize(num, m_output_channels, output_width, output_height, m_min_pad_x, m_min_pad_y);
}

void LocalShareLayer::Forward()
{
	CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;

	int output_width = output->GetWidth();
	int output_height = output->GetHeight();

	const int input_channels = input->GetChannels();
	const int input_stride = input->GetStrideH();
	const int output_stride = output->GetStrideH();

	for (int n = 0; n < input->GetNum(); ++n)
	{
		const float *data = input->GetDataPointer(n, 0, -m_padding_x, -m_padding_y);
		float *target = output->GetDataPointer(n, 0, 0, 0);
		for (int y = 0; y < output_height; y++)
		{
			int filter_index_y = y / m_share_y;
			for (int x = 0; x < output_width; x++)
			{
				int filter_index_x = x / m_share_x;
				int filter_index = filter_index_y * m_filter_num_x + filter_index_x;
				const float *weight = m_weight[filter_index]->GetDataPointer();
				const float *biases = m_biases[filter_index]->GetDataPointer();
				ConvDotProduct(
					data + x * m_stride_x * input_channels,
					weight,
					biases,
					target + x * m_output_channels,
					input_channels,
					input_stride,
					m_output_channels,
					m_filter_width,
					m_filter_height);
			}
			data += m_stride_y * input_stride;
			target += output_stride;
		}
	}
}
