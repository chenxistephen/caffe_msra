#include "DataLayerMinusConst.h"
#include "NetUtils.h"

using namespace DNNTestLib;

DataLayerMinusConst::DataLayerMinusConst(std::istream &stream, std::vector<std::string> &input_links)
	: DataLayer(stream, input_links)
{
	m_mean = ReadData<float>(stream);
}

DataLayerMinusConst::~DataLayerMinusConst()
{
}

void DataLayerMinusConst::Dump(std::ostream &stream)
{
	DataLayer::Dump(stream);

	WriteData(stream, m_mean);
}

void DataLayerMinusConst::FeedImage(const unsigned char *img,
	int width, int height, int channels, int stride, void *param)
{
	m_layerActs->Resize(1, channels, width, height, m_min_pad_x, m_min_pad_y);

	const int width_channels = width * channels;
	float *output = m_layerActs->GetDataPointer();
	const int output_stride = m_layerActs->GetStrideH();
	for (int h = 0; h < height; h++)
	{
		for (int i = 0; i < width_channels; i++)
		{
			output[i] = (float)img[i] - m_mean;
			//std::cout << (float)output[i] << " ";
		}
			
		img += stride;
		output += output_stride;
	}
}

void DataLayerMinusConst::FeedCPUData(const CPUData *input)
{
	std::cout << "ERROR: DataLayerMinusConst do not support CPUData input!" << std::endl;
	throw;
}
