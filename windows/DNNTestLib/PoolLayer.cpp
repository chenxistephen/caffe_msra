#include "PoolLayer.h"
#include "NetUtils.h"
#include "LayerActs.h"

using namespace DNNTestLib;


template class PoolLayer < MaxPoolActs >;
template class PoolLayer < AvgPoolActs >;


template <class PoolType>
PoolLayer<PoolType>::PoolLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_stride_x = ReadData<int>(stream);
	m_stride_y = ReadData<int>(stream);
	m_pool_x = ReadData<int>(stream);
	m_pool_y = ReadData<int>(stream);
}

template <class PoolType>
void PoolLayer<PoolType>::Dump(std::ostream &stream)
{
	Layer::Dump(stream);

	WriteData(stream, m_stride_x);
	WriteData(stream, m_stride_y);
	WriteData(stream, m_pool_x);
	WriteData(stream, m_pool_y);
}

template <class PoolType>
PoolLayer<PoolType>::~PoolLayer()
{
}

template <class PoolType>
void PoolLayer<PoolType>::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);

	const CPUData *input = m_inputs[0]->GetActs();
	int num = input->GetNum();
#if USE_SSE
	assert(input->GetChannels() % 4 == 0);
#endif
	int output_width = (int)ceil((float)(input->GetWidth() - m_pool_x) / (float)m_stride_x) + 1;
	int output_height = (int)ceil((float)(input->GetHeight() - m_pool_x) / (float)m_stride_y) + 1;
	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(num, input->GetChannels(), output_width, output_height, m_min_pad_x, m_min_pad_y);
	else
		m_layerActs->Resize(num, input->GetChannels(), output_width, output_height, m_min_pad_x, m_min_pad_y);
}

template <class PoolType>
void PoolLayer<PoolType>::Forward()
{
	CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;

	const int output_height = output->GetHeight();
	const int output_width = output->GetWidth();
	const int output_stride = output->GetStrideH();
	const int input_channels = input->GetChannels();
	const int input_stride = input->GetStrideH();
	const int input_width = input->GetWidth();
	const int input_height = input->GetHeight();
	for (int n = 0; n < input->GetNum(); ++n)
	{
		const float *data = input->GetDataPointer(n, 0, 0, 0);
		float *target = output->GetDataPointer(n, 0, 0, 0);
		for (int y = 0; y < output_height; y++)
		{
			int remain_y = input_height - y * m_stride_y;
			for (int x = 0; x < output_width; x++)
			{
				int remain_x = input_width - x * m_stride_x;
#if USE_SSE
				PoolType::ActsSSE(
					reinterpret_cast<const __m128 *>(data + x * m_stride_x * input_channels),
					reinterpret_cast<__m128 *>(target + x * input_channels),
					input_channels,
					input_stride,
					min(m_pool_x, remain_x),
					min(m_pool_y, remain_y)
					);
#else
				PoolType::Acts(
					data + x * m_stride_x * input_channels,
					target + x * input_channels,
					input_channels,
					input_stride,
					min(m_pool_x, remain_x),
					min(m_pool_y, remain_y)
					);
#endif
			}
			data += m_stride_y * input_stride;
			target += output_stride;
		}
	}

}
