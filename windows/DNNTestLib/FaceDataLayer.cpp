#include "FaceDataLayer.h"
#include "NetUtils.h"
#include "Geometry.h"

using namespace DNNTestLib;

FaceDataLayer::FaceDataLayer(std::istream &stream, std::vector<std::string> &input_links)
	: DataLayerMinusMean(stream, input_links)
{
	m_landmark_id = ReadData<int>(stream);
	m_resize_width = ReadData<int>(stream);
	m_resize_height = ReadData<int>(stream);

	m_buffer = new unsigned char[m_init_channels*m_init_width*m_init_height];
	ThrowException(MemoryAllocationError, "Out of memory", m_buffer != NULL);
}


FaceDataLayer::~FaceDataLayer()
{
	SaveRelease(m_buffer);
}


void FaceDataLayer::Dump(std::ostream &stream)
{
	DataLayerMinusMean::Dump(stream);

	WriteData(stream, m_landmark_id);
	WriteData(stream, m_resize_width);
	WriteData(stream, m_resize_height);
}


void FaceDataLayer::FeedImage(const unsigned char *img,
	int width, int height, int channels, int stride, void *param)
{
	assert(param != NULL);
	const float *facial_points = (const float *)param;
	float center_x = facial_points[2 * m_landmark_id];
	float center_y = facial_points[2 * m_landmark_id + 1];
	float scale_x = (float)(width - 1) / (float)(m_resize_width - 1);
	float scale_y = (float)(height - 1) / (float)(m_resize_height - 1);

	// scale and crop image
	unsigned char *output = m_buffer;
	const int crop_channels = m_init_channels;
	const int crop_width = m_init_width;
	const int crop_height = m_init_height;
	for (int y = 0; y < crop_height; y++)
	{
		float y_ori = ((float)y - (float)(crop_height - 1) / 2.f) * scale_y + center_y;
		for (int x = 0; x < crop_width; x++)
		{
			float x_ori = ((float)x - (float)(crop_width - 1) / 2.f) * scale_x + center_x;
			Bilinear(img, width, height, channels, stride, x_ori, y_ori, output);
			output += crop_channels;
		}
	}
	DataLayerMinusMean::FeedImage(m_buffer, crop_width, crop_height, crop_channels, crop_width*crop_channels, NULL);
}
