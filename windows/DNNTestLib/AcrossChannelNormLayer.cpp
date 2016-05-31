#include "AcrossChannelNormLayer.h"
#include "NetUtils.h"
#include "LayerActs.h"

using namespace DNNTestLib;

AcrossChannelNormLayer::AcrossChannelNormLayer(std::istream &stream, std::vector<std::string> &input_links)
: Layer(stream, input_links)
{
	m_fAlpha = ReadData<float>(stream);
	m_fBeta = ReadData<float>(stream);
	m_iWindowSize = ReadData<int>(stream);
	assert(m_iWindowSize % 2 == 1);
}

AcrossChannelNormLayer::~AcrossChannelNormLayer()
{
}

void AcrossChannelNormLayer::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);

	CPUData *input = m_inputs[0]->GetActs();
	int num = input->GetNum();
	assert(input->GetChannels() >= m_iWindowSize);
	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(num, input->GetChannels(), input->GetWidth(), input->GetHeight(), m_min_pad_x, m_min_pad_y);
	else
		m_layerActs->Resize(num, input->GetChannels(), input->GetWidth(), input->GetHeight(), m_min_pad_x, m_min_pad_y);
	m_dataScale.Resize(1, input->GetChannels(), 1, 1, 0, 0);
}

void AcrossChannelNormLayer::Forward()
{
	CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;
	const int channels = input->GetChannels();
	const int width = input->GetWidth();
	const int height = input->GetHeight();
	const int input_stride = input->GetStrideH();
	const int output_stride = output->GetStrideH();

	float *scale_vec = m_dataScale.GetDataPointer();
	const float alpha = m_fAlpha / (float)m_iWindowSize;

	for (int n = 0; n < input->GetNum(); ++n)
	{
		const float *data = input->GetDataPointer(n, 0, 0, 0);
		float *target = output->GetDataPointer(n, 0, 0, 0);
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				AcrossChannelNormActs(data + x * channels, target + x * channels,
					scale_vec, channels, m_iWindowSize, alpha, m_fBeta);
			}
			data += input_stride;
			target += output_stride;
		}
	}
}

