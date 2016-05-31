#pragma once

#include "DataLayer.h"

namespace DNNTestLib
{
	class DataLayerMinusConst : public DataLayer
	{
	public:
		DataLayerMinusConst(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~DataLayerMinusConst();
		virtual void Dump(std::ostream &stream);
		virtual void FeedImage(const unsigned char *img,
			int width, int height, int channels, int stride, void *param);
		virtual void FeedCPUData(const CPUData *input);

	protected:
		// properties
		float m_mean;
	};
}
