#include "DataLayer.h"
#include "NetUtils.h"

using namespace DNNTestLib;

DataLayer::DataLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_init_num = ReadData<int>(stream);
	m_init_channels = ReadData<int>(stream);
	m_init_width = ReadData<int>(stream);
	m_init_height = ReadData<int>(stream);
}


DataLayer::~DataLayer()
{
}

void DataLayer::Dump(std::ostream &stream)
{
	Layer::Dump(stream);

	WriteData(stream, m_init_num);
	WriteData(stream, m_init_channels);
	WriteData(stream, m_init_width);
	WriteData(stream, m_init_height);
}

void DataLayer::InitActs()
{
	assert(GetInputNumber() == 0);

	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(m_init_num, m_init_channels, m_init_width, m_init_height, m_min_pad_x, m_min_pad_y);
}

void DataLayer::FeedImage(const unsigned char *img,
	int width, int height, int channels, int stride, void *param)
{
	m_layerActs->Resize(1, channels, width, height, m_min_pad_x, m_min_pad_y);

	const int width_channels = width * channels;
	float *output = m_layerActs->GetDataPointer();
	const int output_stride = m_layerActs->GetStrideH();
	for (int h = 0; h < height; h++)
	{
		for (int i = 0; i < width_channels; i++)
			output[i] = (float)img[i];
		img += stride;
		output += output_stride;
	}
}

void DataLayer::FeedCPUData(const CPUData *input)
{
	m_layerActs->Resize(input->GetNum(), input->GetChannels(), input->GetWidth(), input->GetHeight(), m_min_pad_x, m_min_pad_y);
	m_layerActs->CopyDataFrom(input);
}
