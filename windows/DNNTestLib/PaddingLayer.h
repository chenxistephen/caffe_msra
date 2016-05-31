#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class PaddingLayer : public Layer
	{
	public:
		PaddingLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~PaddingLayer();
		virtual void Forward();
		virtual void InitActs();
		virtual int GetInputPadX();
		virtual int GetInputPadY();

	private:
		int m_padding_x;
		int m_padding_y;
	};
}