#pragma once

#include "DataLayerMinusMean.h"

namespace DNNTestLib
{
	class FaceDataLayer : public DataLayerMinusMean
	{
	public:
		FaceDataLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~FaceDataLayer();
		virtual void Dump(std::ostream &stream);
		virtual void FeedImage(const unsigned char *img,
			int width, int height, int channels, int stride, void *param);

	protected:
		// properties
		int m_landmark_id;
		int m_resize_width;
		int m_resize_height;
		unsigned char *m_buffer;
	};
}

