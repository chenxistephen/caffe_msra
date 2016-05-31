#include "DataLayerMinusMean.h"
#include "NetUtils.h"

using namespace DNNTestLib;

DataLayerMinusMean::DataLayerMinusMean(std::istream &stream, std::vector<std::string> &input_links)
	: DataLayer(stream, input_links)
{
	m_mean = CPUData::CreateFromDump(stream);

	assert(m_init_num == 1);
	assert(m_mean->GetNum() == 1);
	assert(m_mean->GetChannels() == m_init_channels);
	assert(m_mean->GetHeight() == m_init_height);
	assert(m_mean->GetWidth() == m_init_width);
}

DataLayerMinusMean::~DataLayerMinusMean()
{
	SaveRelease(m_mean);
}

void DataLayerMinusMean::Dump(std::ostream &stream)
{
	DataLayer::Dump(stream);

	m_mean->Dump(stream);
}

void DataLayerMinusMean::FeedImage(const unsigned char *img, int width, int height, int channels, int stride, void *param)
{
	assert(channels == m_layerActs->GetChannels());
	assert(width == m_layerActs->GetWidth());
	assert(height = m_layerActs->GetHeight());

	const int width_channels = width * channels;
	float *output = m_layerActs->GetDataPointer();
	const int output_stride = m_layerActs->GetStrideH();
	const float *mean = m_mean->GetDataPointer();
	const int mean_stride = m_mean->GetStrideH();
	for (int h = 0; h < height; h++)
	{
		for (int i = 0; i < width_channels; i++)
			output[i] = (float)img[i] - mean[i];
		img += stride;
		output += output_stride;
		mean += mean_stride;
	}
}


void DataLayerMinusMean::FeedCPUData(const CPUData *input)
{
	std::cout << "ERROR: DataLayerMinusMean do not support CPUData input!" << std::endl;
	throw;
}
