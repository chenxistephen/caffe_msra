#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	template <class PoolType>
	class PoolLayerWithPadding : public Layer
	{
	public:
		PoolLayerWithPadding(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~PoolLayerWithPadding();
		virtual void Forward();
		virtual void InitActs();
		virtual void Dump(std::ostream &stream);

	private:
		// attributes
		int m_stride_x, m_stride_y, m_pool_x, m_pool_y, m_pad_x, m_pad_y;
	};

}