#include "SoftmaxLayer.h"
#include <assert.h>

using namespace DNNTestLib;

SoftmaxLayer::SoftmaxLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
}

void SoftmaxLayer::Dump(std::ostream &stream)
{
	Layer::Dump(stream);
}

SoftmaxLayer::~SoftmaxLayer()
{
	m_layerActs = NULL;
}

void SoftmaxLayer::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);

	CPUData *input = m_inputs[0]->GetActs();
	//assert(input->GetWidth() == 1);
	//assert(input->GetHeight() == 1);
	assert(input->GetPaddingX() == 0);
	assert(input->GetPaddingY() == 0);
	m_layerActs = input;
}


void SoftmaxLayer::Forward()
{
	const int channels = m_layerActs->GetChannels();
	for (int n = 0; n < m_layerActs->GetNum(); ++n)
	{
		for (int h = 0; h < m_layerActs->GetHeight(); ++h)
		{
			for (int w = 0; w < m_layerActs->GetWidth(); ++w)
			{
				float exp_sm = 0;

				float* input_data = m_layerActs->GetDataPointer(n, 0, w, h);

				// find max value
				float max_v = -FLT_MAX;
				for (int i = 0; i < channels; i++)
				{
					if (max_v < input_data[i])
						max_v = input_data[i];
				}

				for (int i = 0; i < channels; i++)
				{
					float t = exp(input_data[i] - max_v);
					exp_sm += t;
					input_data[i] = t;
				}
				for (int i = 0; i < channels; i++)
				{
					input_data[i] /= exp_sm;
				}
			}
		}
	}
}