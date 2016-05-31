#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class DataLayer : public Layer
	{
	public:
		DataLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~DataLayer();
		virtual void InitActs();
		virtual void Dump(std::ostream &stream);
		virtual void FeedImage(const unsigned char *img,
			int width, int height, int channels, int stride, void *param);
		virtual void FeedCPUData(const CPUData *input);

	protected:
		int m_init_num;
		int m_init_channels;
		int m_init_width;
		int m_init_height;
	};
}

