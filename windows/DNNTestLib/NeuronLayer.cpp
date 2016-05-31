#include "NeuronLayer.h"
#include "LayerActs.h"

using namespace DNNTestLib;

template class NeuronLayer < ReluNeuronActs >;
template class NeuronLayer < SigmoidNeuronActs >;

template <class NeuronType>
NeuronLayer<NeuronType>::NeuronLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
}

template <class NeuronType>
NeuronLayer<NeuronType>::~NeuronLayer()
{
	m_layerActs = NULL;
}

template <class NeuronType>
int NeuronLayer<NeuronType>::GetInputPadX()
{
	return m_min_pad_x;
}

template <class NeuronType>
int NeuronLayer<NeuronType>::GetInputPadY()
{
	return m_min_pad_y;
}

template <class NeuronType>
void NeuronLayer<NeuronType>::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);
	m_layerActs = m_inputs[0]->GetActs();
#if USE_SSE
	assert((m_layerActs->GetChannels() * m_layerActs->GetWidth()) % 4 == 0);
#endif
}

template <class NeuronType>
void NeuronLayer<NeuronType>::Forward()
{
	assert(m_inputs[0]->GetActs() == GetActs()); // "neuron" contains inplace operations

	const int num = m_layerActs->GetNum();
	const int height = m_layerActs->GetHeight();
#if USE_SSE
	const int line_num = m_layerActs->GetChannels() * m_layerActs->GetWidth() >> 2;
	const int line_stride = m_layerActs->GetStrideH() >> 2;
#else
	const int line_num = m_layerActs->GetChannels() * m_layerActs->GetWidth();
	const int line_stride = m_layerActs->GetStrideH();
#endif
	for (int n = 0; n < num; n++)
	{
#if USE_SSE
		__m128 *input_data = reinterpret_cast<__m128 *>(m_layerActs->GetDataPointer(n, 0, 0, 0));
#else
		float *input_data = m_layerActs->GetDataPointer(n, 0, 0, 0);
#endif
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < line_num; x++)
			{
#if USE_SSE
				input_data[x] = NeuronType::ActsSSE(input_data[x]);
#else
				input_data[x] = NeuronType::Acts(input_data[x]);
#endif
			}
			input_data += line_stride;
		}
	}
}