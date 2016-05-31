#include "WithinChannelNormLayer.h"
#include "NetUtils.h"
#include "LayerActs.h"

using namespace DNNTestLib;

WithinChannelNormLayer::WithinChannelNormLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_alpha = ReadData<float>(stream);
	m_beta = ReadData<float>(stream);
	m_local_x = ReadData<int>(stream);
	m_local_y = ReadData<int>(stream);
}

WithinChannelNormLayer::~WithinChannelNormLayer()
{
}

int WithinChannelNormLayer::GetInputPadX()
{
	return (m_local_x - 1) / 2;
}

int WithinChannelNormLayer::GetInputPadY()
{
	return (m_local_y - 1) / 2;
}

void WithinChannelNormLayer::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);

	CPUData *input = m_inputs[0]->GetActs();
	int num = input->GetNum();
	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(num, input->GetChannels(), input->GetWidth(), input->GetHeight(), m_min_pad_x, m_min_pad_y);
	else
		m_layerActs->Resize(num, input->GetChannels(), input->GetWidth(), input->GetHeight(), m_min_pad_x, m_min_pad_y);
}

void WithinChannelNormLayer::Forward()
{
	CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;
	const int channels = input->GetChannels();
	const int width = input->GetWidth();
	const int height = input->GetHeight();
	const int input_stride = input->GetStrideH();
	const int output_stride = output->GetStrideH();

	const int pad_x = (m_local_x - 1) / 2;
	const int pad_y = (m_local_y - 1) / 2;

	for (int n = 0; n < input->GetNum(); ++n)
	{
		const float *data = input->GetDataPointer(n, 0, -pad_x, -pad_y);
		const float *center = input->GetDataPointer(n, 0, 0, 0);
		float *target = output->GetDataPointer(n, 0, 0, 0);
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
#if USE_SSE
				WithinChannelNormActs(
					reinterpret_cast<const __m128 *>(data + x * channels),
					reinterpret_cast<__m128 *>(target + x * channels),
					reinterpret_cast<const __m128 *>(center + x * channels),
					m_local_x,
					m_local_y,
					channels,
					input_stride,
					m_alpha,
					m_beta
					);
#else
				WithinChannelNormActs(
					data + x * channels,
					target + x * channels,
					center + x * channels,
					m_local_x,
					m_local_y,
					channels,
					input_stride,
					m_alpha,
					m_beta
					);
#endif
			}
			data += input_stride;
			target += output_stride;
			center += input_stride;
		}
	}
}
