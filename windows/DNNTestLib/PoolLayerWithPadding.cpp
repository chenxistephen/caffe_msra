#include "PoolLayerWithPadding.h"
#include "NetUtils.h"
#include "LayerActs.h"

using namespace DNNTestLib;


template class PoolLayerWithPadding < MaxPoolActs >;
template class PoolLayerWithPadding < AvgPoolActs >;


template <class PoolType>
PoolLayerWithPadding<PoolType>::PoolLayerWithPadding(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_stride_x = ReadData<int>(stream);
	m_stride_y = ReadData<int>(stream);
	m_pool_x = ReadData<int>(stream);
	m_pool_y = ReadData<int>(stream);
	m_pad_x = ReadData<int>(stream);
	m_pad_y = ReadData<int>(stream);

}

template <class PoolType>
void PoolLayerWithPadding<PoolType>::Dump(std::ostream &stream)
{
	Layer::Dump(stream);

	WriteData(stream, m_stride_x);
	WriteData(stream, m_stride_y);
	WriteData(stream, m_pool_x);
	WriteData(stream, m_pool_y);
	WriteData(stream, m_pad_x);
	WriteData(stream, m_pad_y);
}

template <class PoolType>
PoolLayerWithPadding<PoolType>::~PoolLayerWithPadding()
{
}

template <class PoolType>
void PoolLayerWithPadding<PoolType>::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);

	const CPUData *input = m_inputs[0]->GetActs();
	int num = input->GetNum();
#if USE_SSE
	assert(input->GetChannels() % 4 == 0);
#endif
	int output_width = (int)ceil((float)(input->GetWidth() + m_pad_x * 2 - m_pool_x) / (float)m_stride_x) + 1;
	int output_height = (int)ceil((float)(input->GetHeight() + m_pad_y * 2 - m_pool_x) / (float)m_stride_y) + 1;
	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(num, input->GetChannels(), output_width, output_height, m_min_pad_x, m_min_pad_y);
	else
		m_layerActs->Resize(num, input->GetChannels(), output_width, output_height, m_min_pad_x, m_min_pad_y);
}

template <class PoolType>
void PoolLayerWithPadding<PoolType>::Forward()
{
	CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;

	// debug
	if (0)
	{
		std::vector<float> temp(110);
		const float *p = input->GetDataPointer();
		for (int i = 0; i < 110; ++i)
		{
			temp[i] = p[i];
			std::cout << temp[i] << " ";
		}


	}

	const int output_height = output->GetHeight();
	const int output_width = output->GetWidth();
	const int output_stride = output->GetStrideH();
	const int input_channels = input->GetChannels();
	const int input_stride = input->GetStrideH();
	const int input_width = input->GetWidth();
	const int input_height = input->GetHeight();
	for (int n = 0; n < input->GetNum(); ++n)
	{
		float *target = output->GetDataPointer(n, 0, 0, 0);
		for (int y = 0; y < output_height; y++)
		{
			int start_y = max(0, min(input_height, y * m_stride_y - m_pad_y));
			int end_y = max(0, min(input_height, y * m_stride_y + m_pool_y - m_pad_y));

			if (end_y - start_y <= 0)
				continue;

			const float *data = input->GetDataPointer(n, 0, 0, start_y);
			for (int x = 0; x < output_width; x++)
			{
				int start_x = max(0, min(input_width, x * m_stride_x - m_pad_x));
				int end_x = max(0, min(input_width, x * m_stride_x + m_pool_x - m_pad_x));

				if (end_x - start_x <= 0)
					continue;
#if USE_SSE
				PoolType::ActsSSE(
					reinterpret_cast<const __m128 *>(data + start_x * input_channels),
					reinterpret_cast<__m128 *>(target + x * input_channels),
					input_channels,
					input_stride,
					end_x - start_x,
					end_y - start_y
					);
#else
				PoolType::Acts(
					data + start_x * input_channels,
					target + x * input_channels,
					input_channels,
					input_stride,
					end_x - start_x,
					end_y - start_y
					);
#endif
			}
			target += output_stride;
		}
	}

}
